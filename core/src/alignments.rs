//! BAM alignment files (reads mapped to a reference, e.g. the output of `bwa mem | samtools
//! sort`), read through their `.bai` / `.csi` index so only the reads in view are decoded.
//!
//! Two things are offered to the genome view:
//! - `reads(window)`: each read's span, strand, mapping quality and where it differs from
//!   the reference (mismatched bases, insertions, deletions), for drawing reads IGV-style.
//! - `coverage(window, bins)`: mean read depth per bin. For wide windows it comes from
//!   per-block totals computed once per sequence by streaming all its reads (`depth_blocks`);
//!   narrow windows are counted exactly from the reads in view.
//!
//! Reads that are unmapped, duplicates or failed QC are ignored (IGV's defaults).

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

use noodles::bam;
use noodles::bgzf;
use noodles::core::{Position, Region};
use noodles::sam;
use noodles::sam::alignment::record::cigar::op::Kind;

/// Reference bases per coverage block (wide-zoom coverage is averaged over these).
pub const BLOCK: u64 = 128;

/// One read, in 0-based end-exclusive reference coordinates.
#[derive(Clone, Debug, PartialEq)]
pub struct Read {
    pub start: u64,
    pub end: u64,
    pub forward: bool,
    pub mapq: u8,
    pub name: String,
    pub mismatches: Vec<(u64, u8)>,  // (reference position, read base) where they differ
    pub insertions: Vec<(u64, u32)>, // (reference position the bases go before, length)
    pub deletions: Vec<(u64, u64)>,  // deleted (or spliced-out) reference span
    // Pairing, for IGV-style colouring (insert size, orientation, where the mate is).
    pub flags: u16,                  // SAM flags (paired 0x1, mate unmapped 0x8, …)
    pub mate_start: Option<u64>,     // 0-based, None if no mate position
    pub mate_elsewhere: bool,        // mate on another reference sequence
    pub insert: i64,                 // SAM TLEN: signed template length, 0 if unknown
}

/// Per-position counts of a pileup: A, C, G, T, N (other), deletion, insertion (before
/// this position). Depth = the first six (reads covering the base, deleted ones included).
pub type Counts = [u32; 7];
pub const DEL: usize = 5;
pub const INS: usize = 6;
/// Clips shorter than this are left out of "clip" variants (aligners clip a few bad bases at read ends).
pub const MIN_CLIP: usize = 10;
/// Clipped reads this close (bp) are counted as one junction.
pub const CLIP_GROUP: u64 = 10;
/// Bases of a clip kept to say where the rest of the read came from, and how many reads of one place are kept.
const TAIL: usize = 24;
const MAX_SAMPLES: usize = 32;

/// What a pileup of a window gathers besides the counts.
struct Piled {
    counts: Vec<Counts>,
    inserted: HashMap<u64, HashMap<Vec<u8>, u32>>,
    clipped: HashMap<u64, u32>,
    covered: Vec<bool>,
    /// (where the aligned part starts or ends, clip on the right of it?) -> the TAIL bases of the clip next to the aligned part.
    samples: HashMap<(u64, bool), Vec<Vec<u8>>>,
    /// Clip positions where a hard clip (bases not in the file) was counted near an end of a circular reference.
    unchecked: HashSet<u64>,
}

/// True when the clipped bases `clip` continue the circle `reference` from where the aligned part ended (`right`) or started: a read
/// that runs over the origin, not a break. Reads have errors, mostly missing or extra bases, so the TAIL bases next to the aligned part
/// are compared by edit distance (a quarter of them may differ), not base by base.
fn continues_circle(reference: &[u8], r: u64, right: bool, clip: &[u8]) -> bool {
    let (n, m) = (reference.len(), clip.len().min(TAIL));
    if n == 0 || m == 0 {
        return false;
    }
    // both read away from the aligned part's edge: the clip, and the circle from where the aligned part stopped (or started)
    let read: Vec<u8> = if right { clip[..m].to_vec() } else { clip[clip.len() - m..].iter().rev().cloned().collect() };
    let width = m + 4;
    let circle: Vec<u8> = (0..width)
        .map(|i| if right { reference[(r as usize + i) % n] } else { reference[(r as usize + n * (1 + width / n) - 1 - i) % n] })
        .collect();
    let mut prev: Vec<usize> = (0..=width).collect();
    for i in 1..=m {
        let mut cur = vec![i; width + 1];
        for j in 1..=width {
            let cost = (!read[i - 1].eq_ignore_ascii_case(&circle[j - 1])) as usize;
            cur[j] = (prev[j - 1] + cost).min(prev[j] + 1).min(cur[j - 1] + 1);
        }
        prev = cur;
    }
    prev.iter().min().is_some_and(|&d| d <= (m / 4).max(2))
}

/// The base at each of the TAIL places where at least 60 percent of at least 3 `samples` agree; None when they do not.
fn consensus(samples: &[Vec<u8>]) -> Option<Vec<u8>> {
    let good: Vec<&Vec<u8>> = samples.iter().filter(|s| s.len() == TAIL).collect();
    if good.len() < 3 {
        return None;
    }
    (0..TAIL)
        .map(|i| {
            let mut tally = [0usize; 4];
            for s in &good {
                tally[base_slot(s[i]).min(3)] += (base_slot(s[i]) < 4) as usize;
            }
            let (best, &n) = tally.iter().enumerate().max_by_key(|(_, n)| **n).unwrap();
            (n * 10 >= good.len() * 6).then_some(b"ACGT"[best])
        })
        .collect()
}

/// Where on the circle `reference` the bases `t` lie (start of the best place, at most 2 bases different), None when no place or
/// two places are equally good.
fn find_on_circle(reference: &[u8], t: &[u8]) -> Option<usize> {
    let (n, k) = (reference.len(), t.len());
    let (mut best, mut at, mut tied) = (3usize, 0usize, false);
    for start in 0..n {
        let mut bad = 0;
        for (i, b) in t.iter().enumerate() {
            if !reference[(start + i) % n].eq_ignore_ascii_case(b) {
                bad += 1;
                if bad > 2 {
                    break;
                }
            }
        }
        if bad <= 2 && k > 0 {
            if bad < best {
                (best, at, tied) = (bad, start, false);
            } else if bad == best {
                tied = true;
            }
        }
    }
    (best <= 2 && !tied).then_some(at)
}

/// A position where enough reads disagree with the reference.
#[derive(Clone, Debug, PartialEq)]
pub struct Variant {
    pub pos: u64,      // 0-based; for a deletion its first deleted base, for an insertion the base it goes before
    pub kind: &'static str, // "snp", "del", "ins"
    pub reference: String,  // reference base(s) ("" for an insertion)
    pub alt: String,        // read base, "" for a deletion, the inserted bases for an insertion
    pub count: u32,         // reads showing it (for a deletion: at its first base)
    pub depth: u32,
}

fn base_slot(b: u8) -> usize {
    match b {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => 4,
    }
}

/// `bgzf::io::Reader<File>` = the gzip-block decompressor over the file; the BAM reader
/// sits on top of it.
type Reader = bam::io::Reader<bgzf::io::Reader<File>>;

pub struct Bam {
    path: PathBuf,
    header: sam::Header,
    index: bam::Index,
    // A Mutex lets methods take `&self` (shared access from Python) while still being able
    // to move the file cursor; queries take turns on it.
    reader: Mutex<Reader>,
    // Per-sequence block totals, filled lazily. `Arc` = shared pointer, so the coverage
    // code can keep using a block list without holding the lock.
    depth: Mutex<HashMap<String, Arc<Vec<u64>>>>,
}

impl Bam {
    /// Open `path` and its index (`path.bai` or `path.csi`). A coordinate-sorted BAM without
    /// an index gets one written next to it (`path.bai`), like `samtools index`.
    pub fn open(path: &Path) -> io::Result<Self> {
        let index = match bam::fs::read_associated_index(path) {
            Ok(index) => index,
            Err(e) if e.kind() == io::ErrorKind::NotFound => {
                let index = bam::fs::index(path).map_err(|e| {
                    io::Error::new(e.kind(), format!("no .bai index next to the file and it cannot be indexed ({e}); sort it with `samtools sort` first"))
                })?;
                if let bam::Index::Bai(bai) = &index {
                    // Not fatal if the folder is read-only: the index is kept in memory.
                    let _ = bam::bai::fs::write(path.with_added_extension("bai"), bai);
                }
                index
            }
            Err(e) => return Err(e),
        };
        let (reader, header) = open_reader(path)?;
        Ok(Bam { path: path.to_path_buf(), header, index, reader: Mutex::new(reader), depth: Mutex::new(HashMap::new()) })
    }

    /// Reference sequences named in the BAM header, with their lengths.
    pub fn sequences(&self) -> Vec<(String, u64)> {
        self.header
            .reference_sequences()
            .iter()
            .map(|(name, map)| (name.to_string(), map.length().get() as u64))
            .collect()
    }

    fn length(&self, seq: &str) -> Option<u64> {
        self.header.reference_sequences().get(seq.as_bytes()).map(|m| m.length().get() as u64)
    }

    /// Reads overlapping [start, end) of `seq`, in start order, at most `limit` (the bool is
    /// true when there were more). `reference` = the reference bases of [start, end), used to
    /// find mismatches (outside the window nothing is drawn, so nothing is compared).
    pub fn reads(&self, seq: &str, start: u64, end: u64, reference: &[u8], limit: usize) -> io::Result<(Vec<Read>, bool)> {
        let mut out = Vec::new();
        let mut more = false;
        self.for_each_record(&mut self.reader.lock().unwrap(), seq, start, end, |record, rstart| {
            if out.len() == limit {
                more = true;
                return Ok(false); // stop reading
            }
            let read = decode(record, rstart, start, reference)?;
            if read.end > start && read.start < end {
                out.push(read);
            }
            Ok(true)
        })?;
        Ok((out, more))
    }

    /// Mean depth in each of `bins` equal bins of [start, end). Wide bins use the block
    /// totals (computed on first use, can take seconds on a big BAM: call `has_depth` first
    /// to know); narrow windows are counted exactly from the reads.
    pub fn coverage(&self, seq: &str, start: u64, end: u64, bins: usize) -> io::Result<Vec<f32>> {
        if end <= start || bins == 0 {
            return Ok(vec![0.0; bins]);
        }
        let width = (end - start) as f64 / bins as f64;
        if width >= 4.0 * BLOCK as f64 {
            let blocks = self.depth_blocks(seq)?;
            return Ok((0..bins)
                .map(|b| {
                    let a = (start as f64 + b as f64 * width) as u64 / BLOCK;
                    let z = ((start as f64 + (b + 1) as f64 * width) as u64 / BLOCK).max(a + 1);
                    let z = (z as usize).min(blocks.len());
                    let a = (a as usize).min(z);
                    let total: u64 = blocks[a..z].iter().sum();
                    total as f32 / ((z - a).max(1) as u64 * BLOCK) as f32
                })
                .collect());
        }
        // Exact: per-base depth over the window, then averaged per bin.
        let mut depth = vec![0u32; (end - start) as usize];
        self.for_each_record(&mut self.reader.lock().unwrap(), seq, start, end, |record, rstart| {
            for_each_aligned(record, rstart, |pos, _| {
                if pos >= start && pos < end {
                    depth[(pos - start) as usize] += 1;
                }
            })?;
            Ok(true)
        })?;
        Ok((0..bins)
            .map(|b| {
                let a = (b as f64 * width) as usize;
                let z = (((b + 1) as f64 * width) as usize).clamp(a + 1, depth.len());
                let a = a.min(z - 1);
                depth[a..z].iter().map(|&d| d as u64).sum::<u64>() as f32 / (z - a) as f32
            })
            .collect())
    }

    /// Base counts per position of [start, end) (see `Counts`). Bases below `min_baseq`
    /// and reads below `min_mapq` are not counted.
    pub fn pileup(&self, seq: &str, start: u64, end: u64, min_baseq: u8, min_mapq: u8) -> io::Result<Vec<Counts>> {
        let piled = self.pile(&mut self.reader.lock().unwrap(), seq, start, end, min_baseq, min_mapq, false, None)?;
        Ok(piled.counts)
    }

    /// Positions of `seq` where at least `min_fraction` of at least `min_depth` reads
    /// disagree with `reference` (the whole reference sequence): substitutions, deletions
    /// (consecutive deleted bases merged into one) and insertions (the most frequent
    /// inserted bases), plus two signs of larger changes the reads cannot spell out: "clip", where
    /// reads are clipped by at least MIN_CLIP bases (a junction: the edge of a large insertion, an
    /// inversion or a rearrangement), and "gap", a stretch no read covers (a deletion longer than the reads
    /// can span, or an uncovered end). Streams the reads 1 Mb at a time with its own file handle.
    ///
    /// `circular`: the reference is a circle. A read that runs over the origin (its clipped bases continue the
    /// circle) is not a break; a break anywhere, the ends included, is reported. Reads that jump forward over the origin
    /// (their clip is found again on the circle, fewer than half a circle on) are one "del" across the origin: its `alt` is
    /// "across the origin" and its `reference` the deleted bases, from `pos` to the end and on from the start. A gap at both
    /// ends is one gap through the origin in the same way. Nothing is dropped for being near an end.
    #[allow(clippy::too_many_arguments)]
    pub fn variants(&self, seq: &str, reference: &[u8], min_depth: u32, min_fraction: f32, min_baseq: u8, min_mapq: u8, circular: bool) -> io::Result<Vec<Variant>> {
        let len = self.length(seq).ok_or_else(|| unknown(seq))?.min(reference.len() as u64);
        let reference = &reference[..len as usize];
        let circle = circular.then_some(reference);
        let (mut reader, _) = open_reader(&self.path)?;
        let mut out: Vec<Variant> = Vec::new();
        // Clipped reads within CLIP_GROUP bp are one junction (both sides of a breakpoint, a target-site duplication):
        // (position with the most clips, its clips, all clips, deepest depth, last position seen).
        let mut group: Option<(u64, u32, u32, u32, u64)> = None;
        // Kept apart from `out` until the end: a junction pushed there would split a gap still being extended.
        let mut junctions: Vec<Variant> = Vec::new();
        let mut samples: HashMap<(u64, bool), Vec<Vec<u8>>> = HashMap::new();
        let mut unchecked: HashSet<u64> = HashSet::new();
        // Reads crossing a junction split their evidence between its two sides: half the variant share is enough.
        let flush = |g: Option<(u64, u32, u32, u32, u64)>, out: &mut Vec<Variant>| {
            if let Some((pos, _, n, depth, _)) = g {
                if n >= 3 && n as f32 >= 0.5 * min_fraction * depth.max(n) as f32 {
                    out.push(Variant { pos, kind: "clip", reference: String::new(), alt: String::new(), count: n, depth: depth.max(n) });
                }
            }
        };
        const WINDOW: u64 = 1_000_000;
        let mut start = 0;
        while start < len {
            let end = (start + WINDOW).min(len);
            let Piled { counts, inserted, clipped, covered, samples: found, unchecked: hard } = self.pile(&mut reader, seq, start, end, min_baseq, min_mapq, true, circle)?;
            for (key, list) in found {
                let kept = samples.entry(key).or_default();
                kept.extend(list.into_iter().take(MAX_SAMPLES.saturating_sub(kept.len())));
            }
            unchecked.extend(hard);
            for (i, c) in counts.iter().enumerate() {
                let pos = start + i as u64;
                let depth: u32 = c[..6].iter().sum();
                let want = reference[pos as usize].to_ascii_uppercase();
                if !covered[i] {   // no read at all (not merely none passing the quality filters)
                    // shortcut: a gap keeps every reference base it spans (a whole uncovered contig is one long string)
                    match out.last_mut() {
                        Some(v) if v.kind == "gap" && v.pos + v.reference.len() as u64 == pos => v.reference.push(want as char),
                        _ => out.push(Variant { pos, kind: "gap", reference: (want as char).to_string(), alt: String::new(), count: 0, depth: 0 }),
                    }
                }
                let clips = clipped.get(&pos).copied().unwrap_or(0);
                if clips > 0 {
                    group = match group {
                        Some((best, most, n, d, last)) if pos - last <= CLIP_GROUP => {
                            let (best, most) = if clips > most { (pos, clips) } else { (best, most) };
                            Some((best, most, n + clips, d.max(depth), pos))
                        }
                        g => {
                            flush(g, &mut junctions);
                            Some((pos, clips, clips, depth, pos))
                        }
                    };
                }
                if depth < min_depth.max(1) {
                    continue;
                }
                let enough = |n: u32| n as f32 >= min_fraction * depth as f32 && n > 0;
                // Substitution: the most frequent base that is not the reference one.
                let (slot, &n) = c[..4].iter().enumerate().filter(|(s, _)| *s != base_slot(want)).max_by_key(|(_, n)| **n).unwrap();
                if base_slot(want) < 4 && enough(n) {
                    out.push(Variant { pos, kind: "snp", reference: (want as char).to_string(), alt: (b"ACGT"[slot] as char).to_string(), count: n, depth });
                }
                if enough(c[DEL]) {
                    // Extend the previous deletion if it ends right here.
                    match out.last_mut() {
                        Some(v) if v.kind == "del" && v.pos + v.reference.len() as u64 == pos => v.reference.push(want as char),
                        _ => out.push(Variant { pos, kind: "del", reference: (want as char).to_string(), alt: String::new(), count: c[DEL], depth }),
                    }
                }
                if enough(c[INS]) {
                    let best = inserted.get(&pos).and_then(|m| m.iter().max_by_key(|(_, n)| **n)).map(|(s, _)| s.clone());
                    out.push(Variant { pos, kind: "ins", reference: String::new(), alt: String::from_utf8_lossy(&best.unwrap_or_default()).into_owned(), count: c[INS], depth });
                }
            }
            start = end;
        }
        flush(group, &mut junctions);
        if let Some(rf) = circle {
            let n = rf.len() as u64;
            // A hard clip has no bases in the file: near an end of a circle it cannot be told from a read running over the origin.
            // It is reported, with that said, never dropped.
            for j in junctions.iter_mut() {
                if unchecked.iter().any(|&u| u.abs_diff(j.pos) <= CLIP_GROUP) {
                    j.alt = "hard-clipped reads near an end of a circular sequence: they may be reads running over the origin (their bases are not in the file, so this could not be checked)".into();
                }
            }
            // Reads that leave the reference at p1 and go on at q < p1 jumped forward over the origin: a deletion across it,
            // when that is the shorter reading (the other is a duplication of q..p1).
            let near = |r: u64, p: u64| r.min(n - 1).abs_diff(p) <= CLIP_GROUP;
            let mut consumed = vec![false; junctions.len()];
            let mut crossing: Vec<Variant> = Vec::new();
            for i in 0..junctions.len() {
                if consumed[i] {
                    continue;
                }
                let p = junctions[i].pos;
                let gather = |right: bool| -> Vec<Vec<u8>> {
                    samples.iter().filter(|((r, side), _)| *side == right && near(*r, p)).flat_map(|(_, v)| v.iter().cloned()).collect()
                };
                let jump = match consensus(&gather(true)) {
                    Some(t) => find_on_circle(rf, &t).map(|q| (p, q as u64)),
                    None => consensus(&gather(false)).and_then(|h| find_on_circle(rf, &h)).map(|s| ((s as u64 + TAIL as u64) % n, p)),
                };
                if let Some((p1, q)) = jump {
                    if q < p1 && ((n - p1) + q) * 2 < n {
                        let deleted: String = rf[p1 as usize..].iter().chain(&rf[..q as usize]).map(|b| b.to_ascii_uppercase() as char).collect();
                        let (count, depth) = (junctions[i].count, junctions[i].depth);
                        crossing.push(Variant { pos: p1, kind: "del", reference: deleted, alt: "across the origin".into(), count, depth });
                        consumed[i] = true;
                        for k in 0..junctions.len() {
                            if !consumed[k] && (junctions[k].pos.abs_diff(q) <= CLIP_GROUP || junctions[k].pos.abs_diff(p1) <= CLIP_GROUP) {
                                consumed[k] = true;
                            }
                        }
                    }
                }
            }
            junctions = junctions.into_iter().zip(consumed).filter(|(_, gone)| !gone).map(|(j, _)| j).collect();
            let deleted: Vec<(u64, u64)> = crossing.iter().map(|d| (d.pos, d.pos + d.reference.len() as u64)).collect();
            out.extend(crossing);
            // An uncovered stretch at the start and at the end of a circle is one stretch through the origin.
            let first = out.iter().position(|v| v.kind == "gap" && v.pos == 0);
            let last = out.iter().rposition(|v| v.kind == "gap" && v.pos + v.reference.len() as u64 == n);
            if let (Some(f), Some(l)) = (first, last) {
                if f != l {
                    let head = out.remove(f);
                    let l = if l > f { l - 1 } else { l };
                    out[l].reference.push_str(&head.reference);
                    out[l].alt = "across the origin".into();
                }
            }
            // The uncovered stretch of a deletion through the origin is that deletion, not a second finding.
            out.retain(|v| {
                let (gs, ge) = (v.pos, v.pos + v.reference.len() as u64);
                v.kind != "gap" || !deleted.iter().any(|&(ds, de)| [0, n].iter().any(|&shift| gs + shift + CLIP_GROUP >= ds && ge + shift <= de + CLIP_GROUP))
            });
        }
        out.extend(junctions);
        out.sort_by_key(|v| v.pos); // stable: the order of variants at one position is kept
        Ok(out)
    }

    /// Pileup of [start, end) with the given reader; with `with_ins`, also the inserted
    /// bases seen before each position ({position: {bases: reads}}) and how many reads are
    /// clipped (soft or hard, at least MIN_CLIP bases) at each position ({position: reads}), and
    /// which positions any read spans, whatever its base or mapping quality (for gaps).
    /// `circle`: the reference when it is circular: a clip near an end whose bases continue the circle is a read over the
    /// origin, and is not counted.
    #[allow(clippy::too_many_arguments)]
    fn pile(&self, reader: &mut Reader, seq: &str, start: u64, end: u64, min_baseq: u8, min_mapq: u8, with_ins: bool, circle: Option<&[u8]>) -> io::Result<Piled> {
        let len = self.length(seq).ok_or_else(|| unknown(seq))?;
        let end = end.min(len);
        let mut counts = vec![[0u32; 7]; end.saturating_sub(start) as usize];
        let mut inserted: HashMap<u64, HashMap<Vec<u8>, u32>> = HashMap::new();
        let mut clipped: HashMap<u64, u32> = HashMap::new();
        let mut samples: HashMap<(u64, bool), Vec<Vec<u8>>> = HashMap::new();
        let mut unchecked: HashSet<u64> = HashSet::new();
        let mut covered = vec![false; if with_ins { counts.len() } else { 0 }];
        let inside = |p: u64| p >= start && p < end;
        self.for_each_record(reader, seq, start, end, |record, rstart| {
            if with_ins {
                let mut r = rstart;
                for op in record.cigar().iter() {
                    let op = op?;
                    if matches!(op.kind(), Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Deletion) {
                        for pos in r.max(start)..(r + op.len() as u64).min(end) {
                            covered[(pos - start) as usize] = true;
                        }
                    }
                    if matches!(op.kind(), Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Deletion | Kind::Skip) {
                        r += op.len() as u64;
                    }
                }
            }
            if record.mapping_quality().map_or(255, |m| m.get()) < min_mapq {
                return Ok(true);
            }
            let bases = record.sequence();
            let quals = record.quality_scores();
            let quals = quals.as_bytes();
            let (mut r, mut q) = (rstart, 0usize);
            let mut clip_counted_at = None;
            let mut aligned = false;   // (a clip before any aligned base is on the left of the read's alignment)
            for op in record.cigar().iter() {
                let op = op?;
                let n = op.len();
                match op.kind() {
                    Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                        aligned = true;
                        for k in 0..n {
                            let pos = r + k as u64;
                            // 255 = "no qualities stored": counted.
                            let qual = quals.get(q + k).copied().unwrap_or(255);
                            if inside(pos) && (qual == 255 || qual >= min_baseq) {
                                counts[(pos - start) as usize][base_slot(bases.get(q + k).unwrap_or(b'N'))] += 1;
                            }
                        }
                        r += n as u64;
                        q += n;
                    }
                    Kind::Insertion => {
                        aligned = true;
                        if inside(r) {
                            counts[(r - start) as usize][INS] += 1;
                            if with_ins {
                                let ins: Vec<u8> = (q..q + n).map(|i| bases.get(i).unwrap_or(b'N')).collect();
                                *inserted.entry(r).or_default().entry(ins).or_default() += 1;
                            }
                        }
                        q += n;
                    }
                    Kind::Deletion => {
                        aligned = true;
                        for k in 0..n as u64 {
                            if inside(r + k) {
                                counts[(r + k - start) as usize][DEL] += 1;
                            }
                        }
                        r += n as u64;
                    }
                    Kind::Skip => {
                        aligned = true;
                        r += n as u64;
                    }
                    Kind::SoftClip | Kind::HardClip => {
                        // a clip at the very end of the reference (r == len) is a break there too: counted at the last base
                        let at_end = aligned && r == len && end == len;
                        if with_ins && n >= MIN_CLIP && (inside(r) || at_end) && clip_counted_at != Some(r) {
                            let clip: Vec<u8> = if op.kind() == Kind::SoftClip { (q..q + n).map(|i| bases.get(i).unwrap_or(b'N')).collect() } else { Vec::new() };
                            let near_end = if aligned { r + CLIP_GROUP >= len } else { r <= CLIP_GROUP };
                            let mut counted = true;
                            if let (Some(rf), true) = (circle, near_end) {
                                if op.kind() == Kind::HardClip {
                                    unchecked.insert(r.min(len - 1));
                                } else if continues_circle(rf, r, aligned, &clip) {
                                    counted = false;   // a read running over the origin
                                }
                            }
                            if counted {
                                *clipped.entry(r.min(len - 1)).or_default() += 1;   // r: where the aligned part starts (or ended)
                                clip_counted_at = Some(r);
                                if clip.len() >= TAIL {
                                    let kept = samples.entry((r, aligned)).or_default();
                                    if kept.len() < MAX_SAMPLES {
                                        kept.push(if aligned { clip[..TAIL].to_vec() } else { clip[clip.len() - TAIL..].to_vec() });
                                    }
                                }
                            }
                        }
                        if op.kind() == Kind::SoftClip {
                            q += n;
                        }
                    }
                    Kind::Pad => {}
                }
            }
            Ok(true)
        })?;
        Ok(Piled { counts, inserted, clipped, covered, samples, unchecked })
    }

    /// True once the block totals of `seq` are computed (wide coverage is then instant).
    pub fn has_depth(&self, seq: &str) -> bool {
        self.depth.lock().unwrap().contains_key(seq)
    }

    /// Aligned bases per BLOCK-sized block of `seq`, computed once by streaming every read
    /// of that sequence. Uses its own file handle, so drawing reads is not blocked meanwhile.
    pub fn depth_blocks(&self, seq: &str) -> io::Result<Arc<Vec<u64>>> {
        if let Some(d) = self.depth.lock().unwrap().get(seq) {
            return Ok(d.clone()); // cloning an Arc copies the pointer, not the list
        }
        let len = self.length(seq).ok_or_else(|| unknown(seq))?;
        let mut blocks = vec![0u64; len.div_ceil(BLOCK) as usize];
        let (mut reader, _) = open_reader(&self.path)?;
        self.for_each_record(&mut reader, seq, 0, len, |record, rstart| {
            for_each_aligned(record, rstart, |pos, _| {
                if let Some(b) = blocks.get_mut((pos / BLOCK) as usize) {
                    *b += 1;
                }
            })?;
            Ok(true)
        })?;
        let blocks = Arc::new(blocks);
        self.depth.lock().unwrap().insert(seq.to_string(), blocks.clone());
        Ok(blocks)
    }

    /// Call `f(record, 0-based start)` for each kept read overlapping [start, end) of `seq`,
    /// until `f` returns false.
    fn for_each_record(
        &self,
        reader: &mut Reader,
        seq: &str,
        start: u64,
        end: u64,
        mut f: impl FnMut(&bam::Record, u64) -> io::Result<bool>,
    ) -> io::Result<()> {
        let len = self.length(seq).ok_or_else(|| unknown(seq))?;
        let end = end.min(len);
        if end <= start {
            return Ok(());
        }
        // noodles positions are 1-based and inclusive: [start, end) -> start+1..=end.
        let interval = Position::try_from(start as usize + 1).unwrap()..=Position::try_from(end as usize).unwrap();
        let region = Region::new(seq, interval);
        for result in reader.query(&self.header, &self.index, &region)?.records() {
            let record = result?;
            let flags = record.flags();
            if flags.is_unmapped() || flags.is_duplicate() || flags.is_qc_fail() {
                continue;
            }
            let Some(pos) = record.alignment_start().transpose()? else { continue };
            if !f(&record, pos.get() as u64 - 1)? {
                break;
            }
        }
        Ok(())
    }
}

impl Bam {
    /// Write this whole BAM file as plain-text SAM at `out` (File > Export > Alignment): the
    /// header, then every record unchanged, in the BAM's order.
    ///
    /// Rust notes (see docs/RUST_PRIMER.md):
    /// - `io::Result<()>` = "nothing on success, or an I/O error"; each `?` returns the error
    ///   to the caller at once (Python then sees an OSError).
    /// - `use ... as _` brings a trait's methods into scope without its name: noodles' own
    ///   `Write` trait (for alignment records) would otherwise clash with `std::io::Write`.
    /// - Opening a fresh reader (not the shared `self.reader`) keeps the view's drawing
    ///   free to read the same file while this runs in a background thread.
    pub fn write_sam(&self, out: &Path) -> io::Result<()> {
        use noodles::sam::alignment::io::Write as _;
        let (mut reader, header) = open_reader(&self.path)?;
        // BufWriter gathers small writes into big ones (a SAM line per read would be slow).
        let mut writer = sam::io::Writer::new(io::BufWriter::new(File::create(out)?));
        writer.write_header(&header)?;
        // `records()` is an iterator of io::Result<Record>: one per read, read lazily.
        for result in reader.records() {
            let record = result?;
            writer.write_alignment_record(&header, &record)?;
        }
        // Flush what BufWriter still holds; an error here (disk full) must reach the user.
        std::io::Write::flush(writer.get_mut())?;
        Ok(())
    }
}

fn open_reader(path: &Path) -> io::Result<(Reader, sam::Header)> {
    let mut reader = File::open(path).map(bam::io::Reader::new)?;
    let header = reader.read_header()?;
    Ok((reader, header))
}

fn unknown(seq: &str) -> io::Error {
    io::Error::new(io::ErrorKind::NotFound, format!("sequence '{seq}' is not in the BAM file"))
}

/// Call `f(reference position, read base)` for every read base aligned to a reference base
/// (CIGAR M, = and X), walking the CIGAR from the read's start.
fn for_each_aligned(record: &bam::Record, rstart: u64, mut f: impl FnMut(u64, u8)) -> io::Result<()> {
    let seq = record.sequence();
    let (mut r, mut q) = (rstart, 0usize); // reference and read cursors
    for op in record.cigar().iter() {
        let op = op?;
        let n = op.len();
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                for k in 0..n {
                    // A read stored without bases ("*") still covers the reference.
                    f(r + k as u64, seq.get(q + k).unwrap_or(b'N'));
                }
                r += n as u64;
                q += n;
            }
            Kind::Insertion | Kind::SoftClip => q += n,
            Kind::Deletion | Kind::Skip => r += n as u64,
            Kind::HardClip | Kind::Pad => {}
        }
    }
    Ok(())
}

/// Record -> Read, comparing aligned bases with `reference` (the bases of [win_start, ...)).
fn decode(record: &bam::Record, rstart: u64, win_start: u64, reference: &[u8]) -> io::Result<Read> {
    let mut read = Read {
        start: rstart,
        end: rstart,
        forward: !record.flags().is_reverse_complemented(),
        mapq: record.mapping_quality().map_or(255, |m| m.get()),
        name: record.name().map(|n| n.to_string()).unwrap_or_default(),
        mismatches: Vec::new(),
        insertions: Vec::new(),
        deletions: Vec::new(),
        flags: record.flags().bits(),
        mate_start: record.mate_alignment_start().transpose()?.map(|p| p.get() as u64 - 1),
        mate_elsewhere: match (record.reference_sequence_id().transpose()?, record.mate_reference_sequence_id().transpose()?) {
            (Some(a), Some(b)) => a != b,
            _ => false,
        },
        insert: record.template_length() as i64,
    };
    let mut r = rstart;
    for op in record.cigar().iter() {
        let op = op?;
        let n = op.len() as u64;
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => r += n,
            Kind::Insertion => read.insertions.push((r, n as u32)),
            Kind::Deletion | Kind::Skip => {
                read.deletions.push((r, r + n));
                r += n;
            }
            _ => {}
        }
    }
    read.end = r;
    for_each_aligned(record, rstart, |pos, base| {
        let Some(i) = pos.checked_sub(win_start) else { return };
        if let Some(&want) = reference.get(i as usize) {
            let want = want.to_ascii_uppercase();
            if base != want && base != b'N' && want != b'N' {
                read.mismatches.push((pos, base));
            }
        }
    })?;
    Ok(read)
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles::sam::alignment::io::Write;
    use noodles::sam::alignment::record::cigar::Op;
    use noodles::sam::alignment::record::{Flags, MappingQuality};
    use noodles::sam::alignment::record_buf::{Cigar, Sequence};
    use noodles::sam::alignment::RecordBuf;
    use noodles::sam::header::record::value::map::header::{sort_order::COORDINATE, tag::SORT_ORDER};
    use noodles::sam::header::record::value::map::{self, ReferenceSequence};
    use noodles::sam::header::record::value::Map;
    use std::num::NonZero;

    /// Write a small sorted BAM on a 1000-bp reference "ref" (no index: `open` makes one).
    fn write_bam(path: &Path, reads: &[(usize, &str, Vec<Op>, bool)]) {
        let header = sam::Header::builder()
            .set_header(Map::<map::Header>::builder().insert(SORT_ORDER, COORDINATE).build().unwrap())
            .set_reference_sequences(
                [("ref".into(), Map::<ReferenceSequence>::new(NonZero::new(1000).unwrap()))].into_iter().collect(),
            )
            .build();
        let mut writer = bam::io::Writer::new(File::create(path).unwrap());
        writer.write_header(&header).unwrap();
        for (i, (start, bases, cigar, reverse)) in reads.iter().enumerate() {
            let flags = if *reverse { Flags::REVERSE_COMPLEMENTED } else { Flags::empty() };
            let record = RecordBuf::builder()
                .set_name(format!("r{i}"))
                .set_flags(flags)
                .set_reference_sequence_id(0)
                .set_alignment_start(Position::try_from(start + 1).unwrap())
                .set_mapping_quality(MappingQuality::new(60).unwrap())
                .set_cigar(Cigar::from(cigar.clone()))
                .set_sequence(Sequence::from(bases.as_bytes().to_vec()))
                .build();
            writer.write_alignment_record(&header, &record).unwrap();
        }
        writer.try_finish().unwrap();
    }

    #[test]
    fn reads_mismatches_and_coverage() {
        let path = std::env::temp_dir().join("ge_test_reads.bam");
        let _ = std::fs::remove_file(path.with_added_extension("bai"));
        let op = |k, n| Op::new(k, n);
        // Reference: 1000 x 'A'.
        write_bam(&path, &[
            (10, "AAAAAAAAAA", vec![op(Kind::Match, 10)], false),           // 10..20, perfect
            (15, "AAGAAAAAAA", vec![op(Kind::Match, 10)], true),            // 15..25, G at 17
            (30, "AAAAATTAAAAA", vec![op(Kind::Match, 5), op(Kind::Insertion, 2), op(Kind::Match, 5)], false),
            (50, "AAAAAAAAAA", vec![op(Kind::Match, 5), op(Kind::Deletion, 3), op(Kind::Match, 5)], false),
            // three reads clipped by 12 bases at 700: a junction
            (700, "CCCCCCCCCCCCAAAAAAAA", vec![op(Kind::SoftClip, 12), op(Kind::Match, 8)], false),
            (700, "CCCCCCCCCCCCAAAAAAAA", vec![op(Kind::SoftClip, 12), op(Kind::Match, 8)], true),
            (700, "CCCCCCCCCCCCAAAAAAAA", vec![op(Kind::SoftClip, 12), op(Kind::Match, 8)], false),
            // a fourth whose end is hard- and soft-clipped: one read, one clip
            (700, "CCCCCCCCCCCCAAAAAAAA", vec![op(Kind::HardClip, 5), op(Kind::SoftClip, 12), op(Kind::Match, 8)], false),
            (900, "CCAAAAAAAA", vec![op(Kind::SoftClip, 2), op(Kind::Match, 8)], false),
            // three reads clipped 5 bases from the end (running over the origin of a circular plasmid): no junction
            (985, "AAAAAAAAAACCCCCCCCCCCC", vec![op(Kind::Match, 10), op(Kind::SoftClip, 12)], false),
            (985, "AAAAAAAAAACCCCCCCCCCCC", vec![op(Kind::Match, 10), op(Kind::SoftClip, 12)], false),
            (985, "AAAAAAAAAACCCCCCCCCCCC", vec![op(Kind::Match, 10), op(Kind::SoftClip, 12)], false),
        ]);
        let bam = Bam::open(&path).unwrap();
        assert!(path.with_added_extension("bai").exists(), "index written next to the BAM");
        // Export as SAM: a text header naming the reference, then one line per read.
        let sam_path = std::env::temp_dir().join("ge_test_reads.sam");
        bam.write_sam(&sam_path).unwrap();
        let sam_text = std::fs::read_to_string(&sam_path).unwrap();
        assert!(sam_text.contains("@SQ	SN:ref	LN:1000"));
        assert_eq!(sam_text.lines().filter(|l| !l.starts_with('@')).count(), 12);
        assert!(sam_text.contains("r2	0	ref	31	60	5M2I5M"));
        assert_eq!(bam.sequences(), vec![("ref".to_string(), 1000)]);

        let reference = vec![b'a'; 100]; // lowercase on purpose: compared case-insensitively
        let (reads, more) = bam.reads("ref", 0, 100, &reference, 100).unwrap();
        assert!(!more);
        assert_eq!(reads.len(), 4);
        assert_eq!((reads[0].start, reads[0].end, reads[0].forward, reads[0].mapq), (10, 20, true, 60));
        assert_eq!((reads[1].forward, reads[1].mismatches.clone()), (false, vec![(17, b'G')]));
        assert_eq!((reads[2].end, reads[2].insertions.clone()), (40, vec![(35, 2)]));
        assert!(reads[2].mismatches.is_empty(), "inserted bases are not mismatches");
        assert_eq!((reads[3].end, reads[3].deletions.clone()), (63, vec![(55, 58)]));
        assert_eq!(bam.reads("ref", 0, 100, &reference, 2).unwrap().1, true);
        // Soft-clipped bases (CC) are not aligned, so no mismatch and the read starts at 900.
        let far = bam.reads("ref", 890, 1000, &vec![b'A'; 110], 10).unwrap().0;
        assert_eq!((far[0].start, far[0].end, far[0].mismatches.len()), (900, 908, 0));

        // Exact coverage: bins of 5 bp over 10..30: 5, 10 (two reads overlap 15..20), 5, 0.
        assert_eq!(bam.coverage("ref", 10, 30, 4).unwrap(), vec![1.0, 2.0, 1.0, 0.0]);
        // --- per-base coverage across the indel reads: r2 (5M2I5M) covers 30..40, r3 (5M3D5M) 50..55 and 58..63 ---
        let want: Vec<f32> = (28..65).map(|p| if (30..40).contains(&p) || (50..55).contains(&p) || (58..63).contains(&p) { 1.0 } else { 0.0 }).collect();
        assert_eq!(bam.coverage("ref", 28, 65, 37).unwrap(), want);
        // --- the mismatch scan skips the deleted bases too: C only under the deletion (55..58), so a read
        // whose second half is not shifted past it would show C mismatches ---
        let mut gap_ref = vec![b'A'; 100];
        gap_ref[55..58].copy_from_slice(b"CCC");
        let (r, _) = bam.reads("ref", 0, 100, &gap_ref, 100).unwrap();
        assert_eq!((r[3].start, r[3].end), (50, 63));
        assert!(r[3].mismatches.is_empty(), "{:?}", r[3].mismatches);
        // Block coverage (bins wide enough): 10+10+10+10 aligned bases in block 0, 4x8 + 8 + 3x10 in the second bin
        // (the clipped reads at 700, the read at 900 and the three at 985).
        assert!(!bam.has_depth("ref"));
        let wide = bam.coverage("ref", 0, 1024, 2).unwrap();
        assert!(bam.has_depth("ref"));
        assert_eq!(wide, vec![40.0 / 512.0, 70.0 / 512.0]);
        assert!(bam.reads("chrX", 0, 10, b"", 10).is_err());

        // Pairing info: these test reads are unpaired.
        assert_eq!((reads[1].flags & 0x10, reads[1].mate_start, reads[1].insert), (0x10, None, 0));

        // Pileup: at 17 one read says A, the reverse one says G.
        let pile = bam.pileup("ref", 15, 20, 0, 0).unwrap();
        assert_eq!(pile[2], [1, 0, 1, 0, 0, 0, 0]);
        assert_eq!(bam.pileup("ref", 50, 60, 0, 0).unwrap()[5][DEL], 1); // 55 is deleted
        assert_eq!(bam.pileup("ref", 15, 20, 0, 61).unwrap()[2], [0; 7]); // MAPQ 60 < 61

        // Variants against an all-A reference: the G (half of 2 reads), the TT insertion and
        // the 3-base deletion merged into one; nothing from the soft-clipped CC (under MIN_CLIP).
        let found = bam.variants("ref", &vec![b'A'; 1000], 1, 0.5, 0, 0, false).unwrap();
        let small = |kind: &str| ["snp", "ins", "del"].contains(&kind);
        let brief: Vec<_> = found.iter().filter(|v| small(v.kind)).map(|v| (v.pos, v.kind, v.reference.as_str(), v.alt.as_str(), v.count, v.depth)).collect();
        assert_eq!(brief, vec![(17, "snp", "A", "G", 1, 2), (35, "ins", "", "TT", 1, 1), (55, "del", "AAA", "", 1, 1)]);
        // The 12-base clips at 700 are a junction, and so are the three at 995 (this reference is linear: nothing is dropped for being near
        // an end); stretches no read covers are gaps, the ends included.
        let big: Vec<_> = found.iter().filter(|v| !small(v.kind)).map(|v| (v.pos, v.kind, v.reference.len(), v.count)).collect();
        assert_eq!(big, vec![(0, "gap", 10, 0), (25, "gap", 5, 0), (40, "gap", 10, 0), (63, "gap", 637, 0), (700, "clip", 0, 4), (708, "gap", 192, 0), (908, "gap", 77, 0), (995, "gap", 5, 0), (995, "clip", 0, 3)]);
        // The same reads on a circle whose first 12 bases and last 5 are the three reads clip: they run over the origin (no break there), and the
        // uncovered start and end are one gap through the origin; the clips at 700 are still a break.
        let mut circle = vec![b'A'; 1000];
        circle[..12].copy_from_slice(b"CCCCCCCCCCCC");
        circle[995..].copy_from_slice(b"CCCCC");
        let round: Vec<_> = bam.variants("ref", &circle, 1, 0.5, 0, 0, true).unwrap().iter().filter(|v| !small(v.kind)).map(|v| (v.pos, v.kind, v.reference.len(), v.alt.clone(), v.count)).collect();
        assert_eq!(round[0], (25, "gap", 5, String::new(), 0));
        assert!(round.contains(&(700, "clip", 0, String::new(), 4)) && !round.iter().any(|v| v.1 == "clip" && v.0 > 900), "{round:?}");
        assert!(round.contains(&(995, "gap", 15, "across the origin".to_string(), 0)), "{round:?}");
        assert!(!round.iter().any(|v| v.1 == "gap" && v.0 == 0));
        // Gaps are where no read is, whatever the quality filters: none extra when every read is below min_mapq.
        let filtered = bam.variants("ref", &vec![b'A'; 1000], 1, 0.5, 0, 61, false).unwrap();
        let gaps = |vs: &[Variant]| vs.iter().filter(|v| v.kind == "gap").map(|v| (v.pos, v.reference.len())).collect::<Vec<_>>();
        assert_eq!(gaps(&filtered), gaps(&found));
        // depth < 3: no small variants (the clip has its 3 reads, gaps need none)
        assert!(!bam.variants("ref", &vec![b'A'; 1000], 3, 0.5, 0, 0, false).unwrap().iter().any(|v| small(v.kind)));
    }
}
