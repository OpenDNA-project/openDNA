//! Annotation tracks from GFF3 / GTF / BED files, indexed for fast "what is in this
//! window" queries (millions of lines, e.g. a human annotation, stay interactive).
//!
//! Index: per sequence, items sorted by start, plus `max_end[i]` = the largest end among
//! items 0..=i. To find items overlapping [start, end): binary-search the first item that
//! starts at or after `end`, then walk left while `max_end` says something further left
//! can still reach `start`. Simple, no tree, and fast for typical annotations.

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

/// One annotation, 0-based end-exclusive coordinates. On a circular sequence an item may cross the origin: its `end` is then
/// past the sequence's length (the GFF3 convention), see `Track::keep`.
#[derive(Clone, Debug, PartialEq)]
pub struct Item {
    pub start: u64,
    pub end: u64,
    pub strand: i8, // +1, -1, 0 = unknown
    pub name: String,
    pub kind: String, // GFF type (gene, CDS...) or "bed"
    /// The colour as the file writes it (GFF `color=`/`colour=`, BED itemRgb), not interpreted here; None when there is none.
    pub color: Option<Box<str>>,
    /// GFF `Note`, `product` or `description`; None when there is none.
    pub note: Option<Box<str>>,
    /// A BED line's columns past the ninth (itemRgb), tab-separated, the block columns of a BED12 left out: what a file adds
    /// to each item (TSS height, TU length...); None when there are none.
    pub extra: Option<Box<str>>,
}

pub struct Track {
    items: HashMap<String, Vec<Item>>, // sequence name -> items sorted by start
    max_end: HashMap<String, Vec<u64>>,
    /// Items whose end comes before their start (they cross the origin of a circular sequence), per sequence, kept apart until
    /// `keep` knows the sequence's length.
    wrapping: HashMap<String, Vec<Item>>,
    /// What the file held that was not read (lines that are not valid), as sentences; empty when all was read.
    pub notes: Vec<String>,
    /// Sequence lengths the file states (GFF `##sequence-region`).
    lengths: HashMap<String, u64>,
}


impl Track {
    /// Read a .bed file, or a GFF3/GTF file (anything else). Every line is kept, whole-sequence entries (GFF "region" naming the
    /// chromosome) too: `types` flags them and `keep` chooses what stays. A line that is not valid is skipped and
    /// counted in `notes` (the first one is named); a file with no valid line at all is an error.
    pub fn open(path: &Path) -> io::Result<Self> {
        let lower = path.to_string_lossy().to_lowercase();
        if is_bigbed(path) {   // (by its contents, whatever its name): its entries are BED lines
            return Self::read(path, "bigBed", true, false, bigbed_lines(path)?.into_iter().map(Ok));
        }
        let bed = lower.ends_with(".bed");
        let vcf = lower.ends_with(".vcf") || lower.ends_with(".vcf.gz");
        // `Box<dyn BufRead>` = "any line reader": plain file, or bgzip-compressed VCF
        // (the .vcf.gz that bcftools writes).
        let reader: Box<dyn BufRead> = if lower.ends_with(".gz") {
            Box::new(BufReader::new(noodles::bgzf::io::Reader::new(File::open(path)?)))
        } else {
            Box::new(BufReader::new(File::open(path)?))
        };
        let format = if bed { "BED" } else if vcf { "VCF" } else { "GFF/GTF" };
        Self::read(path, format, bed, vcf, reader.lines())
    }

    /// The items of `lines` (named `format` in the notes; `bed`, `vcf`: how to read them, else GFF/GTF).
    fn read(path: &Path, format: &str, bed: bool, vcf: bool, lines: impl Iterator<Item = io::Result<String>>) -> io::Result<Self> {
        let mut items: HashMap<String, Vec<Item>> = HashMap::new();
        let mut wrapping: HashMap<String, Vec<Item>> = HashMap::new();
        let mut lengths = HashMap::new();
        let (mut skipped, mut first_bad) = (0usize, None);
        for (n, line) in lines.enumerate() {
            let line = line?;
            if line.starts_with("##FASTA") {
                break; // GFF3 may end with the sequences themselves
            }
            if let Some(rest) = line.strip_prefix("##sequence-region") {
                let f: Vec<&str> = rest.split_whitespace().collect();
                if let (Some(seq), Some(Ok(end))) = (f.first(), f.get(2).map(|e| e.parse::<u64>())) {
                    lengths.insert(seq.to_string(), end);
                }
                continue;
            }
            if line.is_empty() || line.starts_with('#') || line.starts_with("track") || line.starts_with("browser") {
                continue;
            }
            let parsed = if bed { parse_bed(&line) } else if vcf { parse_vcf(&line) } else { parse_gff(&line) };
            match parsed {
                Ok(Some((seq, item))) if !bed && !vcf && item.end <= item.start => wrapping.entry(seq).or_default().push(item),
                Ok(Some((seq, item))) => items.entry(seq).or_default().push(item),
                Ok(None) => {} // a line we deliberately skip
                Err(why) => {
                    skipped += 1;
                    first_bad.get_or_insert((n + 1, why));
                }
            }
        }
        let notes = skipped_notes(path, format, skipped, first_bad, items.is_empty() && wrapping.is_empty())?;
        let mut track = Track::from_items(items);
        track.notes = notes;
        track.wrapping = wrapping;
        track.lengths = lengths;
        Ok(track)
    }

    /// The types of the items (GFF column 3; "bed", "variant") in order of name, each with its number of items and whether every
    /// one of them covers its whole sequence (a `region` line that only says which sequence this is: not worth showing). An item
    /// that crosses the origin counts for its type.
    pub fn types(&self) -> Vec<(String, usize, bool)> {
        let mut out: HashMap<&str, (usize, bool)> = HashMap::new();
        for (seq, list) in self.items.iter().chain(self.wrapping.iter()) {
            // The sequence's length: as the file states it, else as far as its items reach (then only a type that names a
            // sequence counts: a gene from the first base to the last item's end is still a gene).
            let stated = self.lengths.get(seq).copied();
            let length = stated.or_else(|| self.max_end.get(seq).and_then(|m| m.last().copied())).unwrap_or(0);
            for item in list {
                let names_it = stated.is_some() || matches!(item.kind.as_str(), "region" | "chromosome" | "source" | "contig");
                let t = out.entry(&item.kind).or_insert((0, true));
                t.0 += 1;
                t.1 &= names_it && item.start == 0 && item.end >= length;
            }
        }
        let mut v: Vec<(String, usize, bool)> = out.into_iter().map(|(k, (n, whole))| (k.to_string(), n, whole)).collect();
        v.sort();
        v
    }

    /// How many items have a colour, and a note (to offer them at import).
    pub fn field_counts(&self) -> (usize, usize) {
        let all = self.items.values().chain(self.wrapping.values()).flatten();
        all.fold((0, 0), |(c, n), i| (c + i.color.is_some() as usize, n + i.note.is_some() as usize))
    }

    /// The colours as the file writes them, each with its number of items (few distinct ones: they are read once each).
    pub fn colors(&self) -> Vec<(String, usize)> {
        let mut out: HashMap<&str, usize> = HashMap::new();
        for c in self.items.values().chain(self.wrapping.values()).flatten().filter_map(|i| i.color.as_deref()) {
            *out.entry(c).or_default() += 1;
        }
        let mut v: Vec<(String, usize)> = out.into_iter().map(|(c, n)| (c.to_string(), n)).collect();
        v.sort();
        v
    }

    /// Keep only the items of `types`. The items of `seq` that cross the origin: with `circular_length` (the sequence is a
    /// circle of that many bases), one item each, its end past the length; without (a linear sequence), left out. Returns how many
    /// were left out that way (of the kept types). The other sequences' origin-crossing items are dropped (they cannot be drawn).
    pub fn keep(&mut self, types: &[String], seq: &str, circular_length: Option<u64>) -> usize {
        let wanted: HashSet<&str> = types.iter().map(String::as_str).collect();
        for list in self.items.values_mut() {
            list.retain(|i| wanted.contains(i.kind.as_str()));
        }
        let mut dropped = 0;
        let wrapping = std::mem::take(&mut self.wrapping).remove(seq).unwrap_or_default();
        for item in wrapping.into_iter().filter(|i| wanted.contains(i.kind.as_str())) {
            match circular_length {
                Some(n) => self.items.entry(seq.to_string()).or_default().push(Item { end: item.end + n, ..item }),
                None => dropped += 1,
            }
        }
        self.items.retain(|_, list| !list.is_empty());
        self.index();
        dropped
    }

    pub fn from_items(items: HashMap<String, Vec<Item>>) -> Self {
        let mut track = Track { items, max_end: HashMap::new(), wrapping: HashMap::new(), notes: Vec::new(), lengths: HashMap::new() };
        track.index();
        track
    }

    /// Sort each sequence's items by start and build `max_end` (see the top of the file).
    fn index(&mut self) {
        self.max_end.clear();
        for (seq, list) in self.items.iter_mut() {
            list.sort_by_key(|i| (i.start, i.end));
            // `scan` = running fold: the running maximum of `end`.
            let ends: Vec<u64> = list.iter().scan(0, |m, i| { *m = (*m).max(i.end); Some(*m) }).collect();
            self.max_end.insert(seq.clone(), ends);
        }
    }

    /// Sequence names with their number of items (those that cross the origin included).
    pub fn sequences(&self) -> Vec<(String, usize)> {
        let mut counts: HashMap<&String, usize> = HashMap::new();
        for (k, l) in self.items.iter().chain(self.wrapping.iter()) {
            *counts.entry(k).or_default() += l.len();
        }
        let mut v: Vec<_> = counts.into_iter().map(|(k, n)| (k.clone(), n)).collect();
        v.sort();
        v
    }

    /// Items overlapping [start, end) on `seq`, in start order, at most `limit` (the bool
    /// is true when there were more).
    pub fn query(&self, seq: &str, start: u64, end: u64, limit: usize) -> (Vec<&Item>, bool) {
        let (Some(list), Some(max_end)) = (self.items.get(seq), self.max_end.get(seq)) else {
            return (Vec::new(), false);
        };
        // First index whose item starts at or after `end`: everything before may overlap.
        let stop = list.partition_point(|i| i.start < end);
        // `max_end[i]` is the largest end among items 0..=i, so it never decreases: the first index where it
        // passes `start` is the first item that can overlap (a binary search, not a walk).
        let first = max_end.partition_point(|&m| m <= start);
        let mut out = Vec::new();
        // Walk forward, so that when there are more than `limit` items it is the LAST ones that are left out.
        for item in &list[first.min(stop)..stop] {
            if item.end > start {
                if out.len() == limit {
                    return (out, true);
                }
                out.push(item);
            }
        }
        (out, false)
    }

    /// Number of items overlapping each of `bins` equal bins of [start, end) (for a
    /// density plot when there are too many items to draw).
    pub fn density(&self, seq: &str, start: u64, end: u64, bins: usize) -> Vec<u32> {
        let mut counts = vec![0u32; bins];
        if end <= start || bins == 0 {
            return counts;
        }
        let width = (end - start) as f64 / bins as f64;
        for item in self.query(seq, start, end, usize::MAX).0 {
            let a = ((item.start.max(start) - start) as f64 / width) as usize;
            let b = (((item.end.min(end) - start) as f64 / width).ceil() as usize).clamp(a + 1, bins);
            for c in &mut counts[a.min(bins - 1)..b] {
                *c += 1;
            }
        }
        counts
    }

    /// Items whose name contains `needle` (case-insensitive), at most `limit`; with `only`, only those on that sequence (in
    /// order of start), so the limit never counts matches elsewhere.
    pub fn search(&self, needle: &str, limit: usize, only: Option<&str>) -> Vec<(&str, &Item)> {
        let needle = needle.to_lowercase();
        let mut out = Vec::new();
        for (seq, list) in self.items.iter().filter(|(s, _)| only.map_or(true, |o| s.as_str() == o)) {
            for item in list {
                if item.name.to_lowercase().contains(&needle) {
                    out.push((seq.as_str(), item));
                    if out.len() == limit {
                        return out;
                    }
                }
            }
        }
        out
    }
}

/// A track line's result: Ok(Some) an item, Ok(None) a line skipped on purpose, Err why the line is not valid.
type Parsed = Result<Option<(String, Item)>, &'static str>;

const BIGBED_MAGIC: u32 = 0x8789_F2EB;

/// A file's first 4 bytes as a little-endian number (what tells a bigWig or a bigBed whatever its name), None when shorter.
pub fn magic(path: &Path) -> Option<u32> {
    let mut bytes = [0u8; 4];
    File::open(path).and_then(|mut f| std::io::Read::read_exact(&mut f, &mut bytes)).ok().map(|_| u32::from_le_bytes(bytes))
}

pub fn is_bigbed(path: &Path) -> bool {
    magic(path) == Some(BIGBED_MAGIC)
}

/// The notes of a text file read line by line: "N lines were skipped as not valid FORMAT (first: line L, why).", or an error
/// when no line at all was valid (`nothing`).
pub fn skipped_notes(path: &Path, format: &str, skipped: usize, first_bad: Option<(usize, &str)>, nothing: bool) -> io::Result<Vec<String>> {
    let Some((line, why)) = first_bad else { return Ok(Vec::new()) };
    if nothing {
        return Err(io::Error::new(io::ErrorKind::InvalidData, format!("{}: line {line} is not valid {format} ({why})", path.display())));
    }
    let count = if skipped == 1 { "1 line was".to_string() } else { format!("{skipped} lines were") };
    Ok(vec![format!("{count} skipped as not valid {format} (first: line {line}, {why}).")])
}

/// Every entry of a bigBed file as the BED line it stands for ("sequence start end" and the rest of its columns).
fn bigbed_lines(path: &Path) -> io::Result<Vec<String>> {
    let bad = |e: &dyn std::fmt::Display| io::Error::new(io::ErrorKind::InvalidData, format!("{}: {e}", path.display()));
    let mut reader = bigtools::BigBedRead::open_file(path).map_err(|e| bad(&e))?;
    let chroms: Vec<(String, u32)> = reader.chroms().iter().map(|c| (c.name.clone(), c.length)).collect();
    let mut lines = Vec::new();
    for (name, length) in chroms {
        for entry in reader.get_interval(&name, 0, length).map_err(|e| bad(&e))? {
            let entry = entry.map_err(|e| bad(&e))?;
            lines.push(if entry.rest.is_empty() { format!("{name}\t{}\t{}", entry.start, entry.end) }
                       else { format!("{name}\t{}\t{}\t{}", entry.start, entry.end, entry.rest) });
        }
    }
    Ok(lines)
}

fn number(s: &str, what: &'static str) -> Result<u64, &'static str> {
    s.parse().map_err(|_| what)
}

/// GFF3 / GTF line -> (sequence, item).
fn parse_gff(line: &str) -> Parsed {
    let f: Vec<&str> = line.split('\t').collect();
    if f.len() < 9 {
        return Err("fewer than 9 columns");
    }
    let start = number(f[3], "start is not a number")?;
    let end = number(f[4], "end is not a number")?;
    if start == 0 || end == 0 {
        return Err("start or end is 0 (GFF counts from 1)");
    }
    // An end before the start: the feature crosses the origin of a circular sequence (Track::keep decides).
    let strand = match f[6] { "+" => 1, "-" => -1, _ => 0 };
    // GFF3: key=value;...   GTF: key "value"; ...
    let mut attrs: HashMap<String, String> = HashMap::new();
    for part in f[8].split(';') {
        let part = part.trim();
        let Some((k, v)) = part.split_once('=').or_else(|| part.split_once(' ')) else { continue };
        attrs.insert(k.trim().to_string(), decode(v.trim().trim_matches('"')));
    }
    let first = |keys: &[&str]| keys.iter().find_map(|k| attrs.get(*k).filter(|v| !v.is_empty())).cloned(); // an empty `Name=` falls through
    let name = first(&["Name", "gene_name", "gene", "locus_tag", "gene_id", "ID", "transcript_id"]).unwrap_or_else(|| f[2].to_string());
    let color = first(&["color", "colour", "Color", "Colour"]).map(String::into_boxed_str);
    let note = first(&["Note", "product", "description"]).map(String::into_boxed_str);
    Ok(Some((f[0].to_string(), Item { start: start - 1, end, strand, name, kind: f[2].to_string(), color, note, extra: None })))
}

/// BED line (0-based start, end-exclusive) -> (sequence, item).
fn parse_bed(line: &str) -> Parsed {
    let f: Vec<&str> = line.split(|c| c == '\t' || c == ' ').filter(|s| !s.is_empty()).collect();
    if f.len() < 3 {
        return Err("fewer than 3 columns");
    }
    let start = number(f[1], "start is not a number")?;
    let end = number(f[2], "end is not a number")?;
    if end < start {
        return Err("end before start");
    }
    let strand = match f.get(5) { Some(&"+") => 1, Some(&"-") => -1, _ => 0 };
    let name = f.get(3).map(|s| s.to_string()).unwrap_or_default();
    // itemRgb (column 9): "r,g,b"; "0" means none.
    let color = f.get(8).filter(|c| **c != "0" && **c != ".").map(|c| c.to_string().into_boxed_str());
    // Past itemRgb: a BED12's blocks (a count, then lists of sizes and starts), then whatever columns the file adds.
    let blocks = f.len() >= 12 && f[9].parse::<u32>().is_ok() && f[10].split(',').all(|v| v.is_empty() || v.parse::<u64>().is_ok());
    let rest = f.get(if blocks { 12 } else { 9 }..).unwrap_or_default();
    let extra = (!rest.is_empty()).then(|| rest.join("\t").into_boxed_str());
    Ok(Some((f[0].to_string(), Item { start, end, strand, name, kind: "bed".into(), color, note: None, extra })))
}

/// VCF line (CHROM POS ID REF ALT ...) -> (sequence, item spanning the REF bases), named
/// "ID REF>ALT" (long alleles shortened), kind "variant".
fn parse_vcf(line: &str) -> Parsed {
    let f: Vec<&str> = line.split('\t').collect();
    if f.len() < 5 {
        return Err("fewer than 5 columns");
    }
    let pos = number(f[1], "position is not a number")?;
    if pos == 0 {
        return Err("position is 0 (VCF counts from 1)");
    }
    let short = |s: &str| if s.len() > 12 { format!("{}…", s.chars().take(10).collect::<String>()) } else { s.to_string() };
    let change = format!("{}>{}", short(f[3]), short(f[4]));
    let name = if f[2] == "." { change } else { format!("{} {change}", f[2]) };
    let item = Item { start: pos - 1, end: pos - 1 + f[3].len().max(1) as u64, strand: 0, name, kind: "variant".into(), color: None, note: None, extra: None };
    Ok(Some((f[0].to_string(), item)))
}

/// GFF3 escapes special characters as %XX (e.g. "%2C" = ",").
fn decode(s: &str) -> String {
    if !s.contains('%') {
        return s.to_string();
    }
    let b = s.as_bytes();
    let mut out = Vec::with_capacity(b.len());
    let mut i = 0;
    while i < b.len() {
        if b[i] == b'%' && i + 2 < b.len() {
            // Bytes, not &s[..]: slicing a str mid-character would panic.
            let hex = std::str::from_utf8(&b[i + 1..i + 3]).ok();
            if let Some(v) = hex.and_then(|h| u8::from_str_radix(h, 16).ok()) {
                out.push(v);
                i += 3;
                continue;
            }
        }
        out.push(b[i]);
        i += 1;
    }
    String::from_utf8_lossy(&out).into_owned()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gff_bed_and_queries() {
        let dir = std::env::temp_dir();
        let gff = dir.join("ge_test.gff3");
        std::fs::write(&gff, "##gff-version 3\n\
chr1\tsrc\tregion\t1\t1000\t.\t+\t.\tID=chr1\n\
chr1\tsrc\tgene\t10\t100\t.\t+\t.\tID=g1;Name=lacZ\n\
chr1\tsrc\tgene\t50\t900\t.\t-\t.\tID=g2;Name=big%2Cgene\n\
chr1\tsrc\tCDS\t200\t300\t.\t+\t0\tID=c1;gene=araC\n\
chr2\tsrc\tgene\t5\t6\t.\t.\t.\tgene_id \"x1\"; gene_name \"tiny\"\n\
##FASTA\n>chr1\nACGT\n").unwrap();
        let mut t = Track::open(&gff).unwrap();
        // --- every line is read; the region line, which only covers its sequence, is the one type not worth showing ---
        assert_eq!(t.sequences(), vec![("chr1".to_string(), 4), ("chr2".to_string(), 1)]);
        assert_eq!(t.types(), vec![("CDS".into(), 1, false), ("gene".into(), 3, false), ("region".into(), 1, true)]);
        assert_eq!(t.keep(&["gene".into(), "CDS".into()], "chr1", None), 0);
        assert_eq!(t.sequences(), vec![("chr1".to_string(), 3), ("chr2".to_string(), 1)]);
        let names = |v: Vec<&Item>| v.iter().map(|i| i.name.clone()).collect::<Vec<_>>();
        // 150..160 is inside the big reverse gene only (found through max_end).
        assert_eq!(names(t.query("chr1", 150, 160, 10).0), vec!["big,gene"]);
        assert_eq!(names(t.query("chr1", 0, 1000, 10).0), vec!["lacZ", "big,gene", "araC"]);
        // The limit keeps the first items in start order and says the list was cut.
        let (first, cut) = t.query("chr1", 0, 1000, 2);
        assert_eq!((names(first), cut), (vec!["lacZ".to_string(), "big,gene".to_string()], true));
        assert!(!t.query("chr1", 0, 1000, 3).1, "3 items fit a limit of 3: not cut");
        assert_eq!(t.query("chr2", 4, 5, 10).0[0].name, "tiny"); // GTF attributes
        assert_eq!(t.query("chr1", 9, 10, 10).0[0].start, 9); // 1-based -> 0-based
        assert_eq!(t.density("chr1", 0, 1000, 10)[0], 2); // lacZ + big gene in bin 0
        // --- every bin (100 bp each): lacZ 9..100 -> bin 0; big 49..900 -> bins 0..8; araC 199..300 -> bins 1, 2
        // (the region line was left out) ---
        assert_eq!(t.density("chr1", 0, 1000, 10), vec![2, 2, 2, 1, 1, 1, 1, 1, 1, 0]);
        assert_eq!(t.search("ARA", 10, None)[0].1.name, "araC");

        let bed = dir.join("ge_test.bed");
        std::fs::write(&bed, "track name=x\nchr1\t0\t10\tfirst\t0\t-\nchr1 20 30\n").unwrap();
        let b = Track::open(&bed).unwrap();
        let hits = b.query("chr1", 5, 25, 10).0;
        assert_eq!((hits[0].strand, hits.len()), (-1, 2));

        std::fs::write(&bed, "chr1\tnot-a-number\t10\n").unwrap();
        assert!(Track::open(&bed).is_err());

        let vcf = dir.join("ge_test.vcf");
        std::fs::write(&vcf, "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\n\
chr1\t100\t.\tA\tG\t50\nchr1\t200\trs7\tACGT\tA\t.\n").unwrap();
        let v = Track::open(&vcf).unwrap();
        let hits = v.query("chr1", 0, 1000, 10).0;
        assert_eq!((hits[0].start, hits[0].end, hits[0].name.as_str()), (99, 100, "A>G"));
        assert_eq!((hits[1].end, hits[1].name.as_str(), hits[1].kind.as_str()), (203, "rs7 ACGT>A", "variant"));
    }

    #[test]
    fn partial_regions_colours_notes_and_features_across_the_origin() {
        let gff = std::env::temp_dir().join("ge_test_wrap.gff3");
        std::fs::write(&gff, "##gff-version 3
##sequence-region c 1 1000
c\tsrc\tregion\t100\t500\t.\t.\t.\tID=a;Name=one;color=#257215
c\tsrc\tregion\t950\t50\t.\t.\t.\tID=w;Name=wraps;colour=185,21,21;Note=across%20the%20origin
c\tsrc\tgene\t600\t700\t.\t+\t.\tID=g;product=lacZ
").unwrap();
        let mut t = Track::open(&gff).unwrap();
        // --- a region that does not cover the sequence is worth showing; the wrapping line is counted with the others ---
        assert_eq!(t.types(), vec![("gene".into(), 1, false), ("region".into(), 2, false)]);
        assert_eq!(t.field_counts(), (2, 2));
        assert_eq!(t.colors(), vec![("#257215".to_string(), 1), ("185,21,21".to_string(), 1)]);
        assert_eq!(t.sequences(), vec![("c".to_string(), 3)]);
        // --- on a circle of 1000 bases, the wrapping feature is ONE item, 949..1050, found from both sides of the origin ---
        assert_eq!(t.keep(&["region".into(), "gene".into()], "c", Some(1000)), 0);
        let (items, _) = t.query("c", 990, 1000, 10);
        assert_eq!((items[0].start, items[0].end, items[0].note.as_deref()), (949, 1050, Some("across the origin")));
        assert_eq!(t.query("c", 1000, 1010, 10).0.len(), 1, "the part after the origin is at 1000.. in the item's coordinates");
        assert_eq!(t.query("c", 600, 601, 10).0[0].note.as_deref(), Some("lacZ"));
        // --- on a linear sequence it is left out and counted; a type not kept goes with its items ---
        let mut linear = Track::open(&gff).unwrap();
        assert_eq!(linear.keep(&["region".into()], "c", None), 1);
        assert_eq!(linear.sequences(), vec![("c".to_string(), 1)]);
        // --- with no stated length, a gene from the first base to the end is a gene, not "the whole sequence" ---
        let one = std::env::temp_dir().join("ge_test_one.gff3");
        std::fs::write(&one, "c\tsrc\tgene\t1\t5000\t.\t+\t.\tName=lacZ\n").unwrap();
        assert_eq!(Track::open(&one).unwrap().types(), vec![("gene".into(), 1, false)]);
        // --- a BED colour (itemRgb), "0" meaning none ---
        let bed = std::env::temp_dir().join("ge_test_rgb.bed");
        std::fs::write(&bed, "c\t0\t10\tx\t0\t+\t0\t10\t255,0,0
c\t20\t30\ty\t0\t+\t20\t30\t0
").unwrap();
        assert_eq!(Track::open(&bed).unwrap().colors(), vec![("255,0,0".to_string(), 1)]);
        // --- the columns past itemRgb are kept (a BED12's blocks left out) ---
        std::fs::write(&bed, "c\t0\t10\tx\t0\t+\t0\t10\t0\theight=81\tRSPM=12\n\
c\t20\t30\ty\t0\t+\t20\t30\t0\t2\t3,4,\t0,6,\tTU=1\nc\t40\t50\tz\n").unwrap();
        let extras: Vec<Option<String>> = Track::open(&bed).unwrap().query("c", 0, 100, 10).0.iter()
            .map(|i| i.extra.as_deref().map(str::to_string)).collect();
        assert_eq!(extras, vec![Some("height=81\tRSPM=12".into()), Some("TU=1".into()), None]);
    }
}
