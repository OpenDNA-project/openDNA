//! Annotation tracks from GFF3 / GTF / BED files, indexed for fast "what is in this
//! window" queries (millions of lines, e.g. a human annotation, stay interactive).
//!
//! Index: per sequence, items sorted by start, plus `max_end[i]` = the largest end among
//! items 0..=i. To find items overlapping [start, end): binary-search the first item that
//! starts at or after `end`, then walk left while `max_end` says something further left
//! can still reach `start`. Simple, no tree, and fast for typical annotations.

use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

/// One annotation, 0-based end-exclusive coordinates.
#[derive(Clone, Debug, PartialEq)]
pub struct Item {
    pub start: u64,
    pub end: u64,
    pub strand: i8, // +1, -1, 0 = unknown
    pub name: String,
    pub kind: String, // GFF type (gene, CDS...) or "bed"
}

pub struct Track {
    items: HashMap<String, Vec<Item>>, // sequence name -> items sorted by start
    max_end: HashMap<String, Vec<u64>>,
}

impl Track {
    /// Read a .bed file, or a GFF3/GTF file (anything else). Whole-sequence entries
    /// (GFF "region"/"chromosome"/"source") are skipped, like the GenBank `source`.
    pub fn open(path: &Path) -> io::Result<Self> {
        let lower = path.to_string_lossy().to_lowercase();
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
        let mut items: HashMap<String, Vec<Item>> = HashMap::new();
        for (n, line) in reader.lines().enumerate() {
            let line = line?;
            if line.starts_with("##FASTA") {
                break; // GFF3 may end with the sequences themselves
            }
            if line.is_empty() || line.starts_with('#') || line.starts_with("track") || line.starts_with("browser") {
                continue;
            }
            let parsed = if bed { parse_bed(&line) } else if vcf { parse_vcf(&line) } else { parse_gff(&line) };
            match parsed {
                Some(Some((seq, item))) => items.entry(seq).or_default().push(item),
                Some(None) => {} // a line we deliberately skip
                None => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("{}: line {} is not valid {format}", path.display(), n + 1),
                    ))
                }
            }
        }
        Ok(Track::from_items(items))
    }

    pub fn from_items(mut items: HashMap<String, Vec<Item>>) -> Self {
        let mut max_end = HashMap::new();
        for (seq, list) in items.iter_mut() {
            list.sort_by_key(|i| (i.start, i.end));
            // `scan` = running fold: the running maximum of `end`.
            let ends: Vec<u64> = list.iter().scan(0, |m, i| { *m = (*m).max(i.end); Some(*m) }).collect();
            max_end.insert(seq.clone(), ends);
        }
        Track { items, max_end }
    }

    /// Sequence names with their number of items.
    pub fn sequences(&self) -> Vec<(String, usize)> {
        let mut v: Vec<_> = self.items.iter().map(|(k, l)| (k.clone(), l.len())).collect();
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

    /// Items whose name contains `needle` (case-insensitive), at most `limit`.
    pub fn search(&self, needle: &str, limit: usize) -> Vec<(&str, &Item)> {
        let needle = needle.to_lowercase();
        let mut out = Vec::new();
        for (seq, list) in &self.items {
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

/// GFF3 / GTF line -> (sequence, item). `None` = malformed, `Some(None)` = skipped type.
fn parse_gff(line: &str) -> Option<Option<(String, Item)>> {
    let f: Vec<&str> = line.split('\t').collect();
    if f.len() < 9 {
        return None;
    }
    if matches!(f[2], "region" | "chromosome" | "source" | "contig") {
        return Some(None);
    }
    let start: u64 = f[3].parse().ok()?;
    let end: u64 = f[4].parse().ok()?;
    if start == 0 || end < start {
        return None;
    }
    let strand = match f[6] { "+" => 1, "-" => -1, _ => 0 };
    // GFF3: key=value;...   GTF: key "value"; ...
    let mut attrs: HashMap<String, String> = HashMap::new();
    for part in f[8].split(';') {
        let part = part.trim();
        let (k, v) = match part.split_once('=') {
            Some(kv) => kv,
            None => match part.split_once(' ') {
                Some(kv) => kv,
                None => continue,
            },
        };
        attrs.insert(k.trim().to_string(), decode(v.trim().trim_matches('"')));
    }
    let name = ["Name", "gene_name", "gene", "locus_tag", "gene_id", "ID", "transcript_id"]
        .iter()
        .find_map(|k| attrs.get(*k).filter(|v| !v.is_empty())) // an empty `Name=` falls through to the next attribute
        .cloned()
        .unwrap_or_else(|| f[2].to_string());
    Some(Some((f[0].to_string(), Item { start: start - 1, end, strand, name, kind: f[2].to_string() })))
}

/// BED line (0-based start, end-exclusive) -> (sequence, item).
fn parse_bed(line: &str) -> Option<Option<(String, Item)>> {
    let f: Vec<&str> = line.split(|c| c == '\t' || c == ' ').filter(|s| !s.is_empty()).collect();
    if f.len() < 3 {
        return None;
    }
    let start: u64 = f[1].parse().ok()?;
    let end: u64 = f[2].parse().ok()?;
    if end < start {
        return None;
    }
    let strand = match f.get(5) { Some(&"+") => 1, Some(&"-") => -1, _ => 0 };
    let name = f.get(3).map(|s| s.to_string()).unwrap_or_default();
    Some(Some((f[0].to_string(), Item { start, end, strand, name, kind: "bed".into() })))
}

/// VCF line (CHROM POS ID REF ALT ...) -> (sequence, item spanning the REF bases), named
/// "ID REF>ALT" (long alleles shortened), kind "variant".
fn parse_vcf(line: &str) -> Option<Option<(String, Item)>> {
    let f: Vec<&str> = line.split('\t').collect();
    if f.len() < 5 {
        return None;
    }
    let pos: u64 = f[1].parse().ok()?;
    if pos == 0 {
        return None;
    }
    let short = |s: &str| if s.len() > 12 { format!("{}…", s.chars().take(10).collect::<String>()) } else { s.to_string() };
    let change = format!("{}>{}", short(f[3]), short(f[4]));
    let name = if f[2] == "." { change } else { format!("{} {change}", f[2]) };
    let item = Item { start: pos - 1, end: pos - 1 + f[3].len().max(1) as u64, strand: 0, name, kind: "variant".into() };
    Some(Some((f[0].to_string(), item)))
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
        let t = Track::open(&gff).unwrap();
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
        assert_eq!(t.search("ARA", 10)[0].1.name, "araC");

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
}
