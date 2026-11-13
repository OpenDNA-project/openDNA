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

use std::collections::HashMap;
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
        let (counts, _) = self.pile(&mut self.reader.lock().unwrap(), seq, start, end, min_baseq, min_mapq, false)?;
        Ok(counts)
    }

    /// Positions of `seq` where at least `min_fraction` of at least `min_depth` reads
    /// disagree with `reference` (the whole reference sequence): substitutions, deletions
    /// (consecutive deleted bases merged into one) and insertions (the most frequent
    /// inserted bases). Streams the reads 1 Mb at a time with its own file handle.
    pub fn variants(&self, seq: &str, reference: &[u8], min_depth: u32, min_fraction: f32, min_baseq: u8, min_mapq: u8) -> io::Result<Vec<Variant>> {
        let len = self.length(seq).ok_or_else(|| unknown(seq))?.min(reference.len() as u64);
        let (mut reader, _) = open_reader(&self.path)?;
        let mut out: Vec<Variant> = Vec::new();
        const WINDOW: u64 = 1_000_000;
        let mut start = 0;
        while start < len {
            let end = (start + WINDOW).min(len);
            let (counts, inserted) = self.pile(&mut reader, seq, start, end, min_baseq, min_mapq, true)?;
            for (i, c) in counts.iter().enumerate() {
                let pos = start + i as u64;
                let depth: u32 = c[..6].iter().sum();
                if depth < min_depth.max(1) {
                    continue;
                }
                let enough = |n: u32| n as f32 >= min_fraction * depth as f32 && n > 0;
                let want = reference[pos as usize].to_ascii_uppercase();
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
        Ok(out)
    }

    /// Pileup of [start, end) with the given reader; with `with_ins`, also the inserted
    /// bases seen before each position ({position: {bases: reads}}).
    #[allow(clippy::type_complexity)]
    fn pile(&self, reader: &mut Reader, seq: &str, start: u64, end: u64, min_baseq: u8, min_mapq: u8, with_ins: bool) -> io::Result<(Vec<Counts>, HashMap<u64, HashMap<Vec<u8>, u32>>)> {
        let len = self.length(seq).ok_or_else(|| unknown(seq))?;
        let end = end.min(len);
        let mut counts = vec![[0u32; 7]; end.saturating_sub(start) as usize];
        let mut inserted: HashMap<u64, HashMap<Vec<u8>, u32>> = HashMap::new();
        let inside = |p: u64| p >= start && p < end;
        self.for_each_record(reader, seq, start, end, |record, rstart| {
            if record.mapping_quality().map_or(255, |m| m.get()) < min_mapq {
                return Ok(true);
            }
            let bases = record.sequence();
            let quals = record.quality_scores();
            let quals = quals.as_bytes();
            let (mut r, mut q) = (rstart, 0usize);
            for op in record.cigar().iter() {
                let op = op?;
                let n = op.len();
                match op.kind() {
                    Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
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
                        for k in 0..n as u64 {
                            if inside(r + k) {
                                counts[(r + k - start) as usize][DEL] += 1;
                            }
                        }
                        r += n as u64;
                    }
                    Kind::Skip => r += n as u64,
                    Kind::SoftClip => q += n,
                    Kind::HardClip | Kind::Pad => {}
                }
            }
            Ok(true)
        })?;
        Ok((counts, inserted))
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
            (900, "CCAAAAAAAA", vec![op(Kind::SoftClip, 2), op(Kind::Match, 8)], false),
        ]);
        let bam = Bam::open(&path).unwrap();
        assert!(path.with_added_extension("bai").exists(), "index written next to the BAM");
        // Export as SAM: a text header naming the reference, then one line per read.
        let sam_path = std::env::temp_dir().join("ge_test_reads.sam");
        bam.write_sam(&sam_path).unwrap();
        let sam_text = std::fs::read_to_string(&sam_path).unwrap();
        assert!(sam_text.contains("@SQ	SN:ref	LN:1000"));
        assert_eq!(sam_text.lines().filter(|l| !l.starts_with('@')).count(), 5);
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
        // Block coverage (bins wide enough): 10+10+10+10+8 aligned bases, all in block 0
        // except the last read (block 7).
        assert!(!bam.has_depth("ref"));
        let wide = bam.coverage("ref", 0, 1024, 2).unwrap();
        assert!(bam.has_depth("ref"));
        assert_eq!(wide, vec![40.0 / 512.0, 8.0 / 512.0]);
        assert!(bam.reads("chrX", 0, 10, b"", 10).is_err());

        // Pairing info: these test reads are unpaired.
        assert_eq!((reads[1].flags & 0x10, reads[1].mate_start, reads[1].insert), (0x10, None, 0));

        // Pileup: at 17 one read says A, the reverse one says G.
        let pile = bam.pileup("ref", 15, 20, 0, 0).unwrap();
        assert_eq!(pile[2], [1, 0, 1, 0, 0, 0, 0]);
        assert_eq!(bam.pileup("ref", 50, 60, 0, 0).unwrap()[5][DEL], 1); // 55 is deleted
        assert_eq!(bam.pileup("ref", 15, 20, 0, 61).unwrap()[2], [0; 7]); // MAPQ 60 < 61

        // Variants against an all-A reference: the G (half of 2 reads), the TT insertion and
        // the 3-base deletion merged into one; nothing from the soft-clipped CC.
        let found = bam.variants("ref", &vec![b'A'; 1000], 1, 0.5, 0, 0).unwrap();
        let brief: Vec<_> = found.iter().map(|v| (v.pos, v.kind, v.reference.as_str(), v.alt.as_str(), v.count, v.depth)).collect();
        assert_eq!(brief, vec![(17, "snp", "A", "G", 1, 2), (35, "ins", "", "TT", 1, 1), (55, "del", "AAA", "", 1, 1)]);
        assert!(bam.variants("ref", &vec![b'A'; 1000], 3, 0.5, 0, 0).unwrap().is_empty()); // depth < 3
    }
}
