//! Signal tracks: a value along the sequence (coverage, RNA-seq, insertion density...) from a bigWig, bedGraph or wig file,
//! read as bins for the view: for each bin the mean of the values over the bases that have one, and the highest value; NaN
//! where the file has no value at all.
//!
//! A bigWig is read from the file at each query (never loaded whole: a track hub can bring 200 of them), through its zoom
//! levels when a bin spans many bases. bedGraph and wig are text: read once into intervals, with zoom levels built here the
//! same way (each level sums blocks 8 times wider than the one below), so a whole-genome view adds up blocks, not bases.

use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;
use std::sync::Mutex;

use bigtools::utils::file::reopen::ReopenableFile;
use bigtools::BigWigRead;

const BIGWIG_MAGIC: u32 = 0x888F_FC26;
const BASE_LEVEL: u64 = 64; // bp of the finest zoom level of a text file

/// The values added up in one bin (or one zoom block): their sum over the bases covered, those bases, and the extremes.
#[derive(Clone, Copy, Debug)]
pub struct Acc {
    pub sum: f64,
    pub covered: f64,
    pub max: f64,
    pub min: f64,
}

impl Default for Acc {
    fn default() -> Self {
        Acc { sum: 0.0, covered: 0.0, max: f64::NEG_INFINITY, min: f64::INFINITY }
    }
}

impl Acc {
    /// (mean over the bases covered, highest), NaN for both where nothing was added.
    pub fn result(&self) -> (f64, f64) {
        if self.covered > 0.0 {
            (self.sum / self.covered, self.max)
        } else {
            (f64::NAN, f64::NAN)
        }
    }

    fn merge(&mut self, o: &Acc) {
        self.sum += o.sum;
        self.covered += o.covered;
        self.max = self.max.max(o.max);
        self.min = self.min.min(o.min);
    }
}

/// Bins of [start, end): add a record spanning [s, e) whose bases sum to `acc` (the share of it inside each bin goes there).
struct Bins {
    start: f64,
    width: f64,
    out: Vec<Acc>,
}

impl Bins {
    fn new(start: u64, end: u64, n: usize) -> Self {
        let n = n.max(1);
        Bins { start: start as f64, width: (end.saturating_sub(start)).max(1) as f64 / n as f64, out: vec![Acc::default(); n] }
    }

    fn add(&mut self, s: u64, e: u64, acc: &Acc) {
        let (s, e) = (s as f64, e as f64);
        let end = self.start + self.width * self.out.len() as f64;
        if e <= self.start || s >= end || e <= s || acc.covered <= 0.0 {
            return;
        }
        let first = ((s.max(self.start) - self.start) / self.width).floor() as usize;
        let last = (((e.min(end) - self.start) / self.width).ceil() as usize).min(self.out.len());
        for b in first..last {
            let (b0, b1) = (self.start + b as f64 * self.width, self.start + (b + 1) as f64 * self.width);
            let part = (e.min(b1) - s.max(b0)) / (e - s);
            if part <= 0.0 {
                continue;
            }
            let bin = &mut self.out[b];
            bin.sum += acc.sum * part;
            bin.covered += acc.covered * part;
            bin.max = bin.max.max(acc.max);
            bin.min = bin.min.min(acc.min);
        }
    }
}

fn value(s: u64, e: u64, v: f64) -> Acc {
    let n = e.saturating_sub(s) as f64;
    Acc { sum: v * n, covered: n, max: v, min: v }
}

/// One sequence of a text file: its intervals (sorted by start), the furthest end among the first i + 1 of them (to find the
/// first one that reaches a position), and zoom levels (block size, blocks from 0).
struct Mem {
    intervals: Vec<(u64, u64, f32)>,
    reach: Vec<u64>,
    levels: Vec<(u64, Vec<Acc>)>,
}

impl Mem {
    fn new(mut intervals: Vec<(u64, u64, f32)>) -> Self {
        intervals.sort_by_key(|iv| (iv.0, iv.1));
        let reach: Vec<u64> = intervals.iter().scan(0, |far, iv| { *far = (*far).max(iv.1); Some(*far) }).collect();
        let length = reach.last().copied().unwrap_or(0);
        let mut levels = Vec::new();
        let mut size = BASE_LEVEL;
        let mut blocks = vec![Acc::default(); (length / size + 1) as usize];
        for &(s, e, v) in &intervals {
            let mut b = s / size;
            while b * size < e {
                let (b0, b1) = (b * size, (b + 1) * size);
                let part = value(s.max(b0), e.min(b1), v as f64);
                blocks[b as usize].merge(&part);
                b += 1;
            }
        }
        loop {
            let done = blocks.len() <= 1;
            levels.push((size, blocks));
            if done {
                break;
            }
            let below = &levels.last().unwrap().1;
            let next: Vec<Acc> = below.chunks(8).map(|c| c.iter().fold(Acc::default(), |mut a, b| { a.merge(b); a })).collect();
            size *= 8;
            blocks = next;
        }
        Mem { intervals, reach, levels }
    }

    fn bins(&self, start: u64, end: u64, n: usize) -> Vec<Acc> {
        let mut bins = Bins::new(start, end, n);
        let half = bins.width / 2.0;
        // The coarsest level whose blocks are at most half a bin wide; none: the intervals themselves.
        match self.levels.iter().rev().find(|(size, _)| *size as f64 <= half) {
            Some((size, blocks)) => {
                let (first, last) = ((start / size) as usize, ((end + size - 1) / size) as usize);
                for (i, acc) in blocks.iter().enumerate().take(last.min(blocks.len())).skip(first) {
                    bins.add(i as u64 * size, (i as u64 + 1) * size, acc);
                }
            }
            None => {
                let first = self.reach.partition_point(|&far| far <= start);   // the first interval that reaches past start
                for &(s, e, v) in self.intervals[first..].iter().take_while(|iv| iv.0 < end) {
                    bins.add(s, e, &value(s, e, v as f64));
                }
            }
        }
        bins.out
    }
}

enum Source {
    Big(Mutex<BigWigRead<ReopenableFile>>),
    Mem(HashMap<String, Mem>),
}

pub struct Signal {
    source: Source,
    sequences: Vec<(String, u64)>,
    /// What the file held that was not read (lines that are not valid), as sentences; empty when all was read.
    pub notes: Vec<String>,
}

fn invalid(msg: String) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg)
}

/// Is this file a bigWig (by its first 4 bytes, whatever its name)?
pub fn is_bigwig(path: &Path) -> bool {
    crate::tracks::magic(path) == Some(BIGWIG_MAGIC)
}

impl Signal {
    /// A bigWig (by its contents), else a bedGraph or wig file (fixedStep / variableStep sections, else 4 columns: sequence,
    /// start, end, value). Lines that are not valid are skipped and counted in `notes`; no valid line at all is an error.
    pub fn open(path: &Path) -> io::Result<Self> {
        if is_bigwig(path) {
            let reader = BigWigRead::open_file(path).map_err(|e| invalid(format!("{}: {e}", path.display())))?;
            let sequences = reader.chroms().iter().map(|c| (c.name.clone(), c.length as u64)).collect();
            return Ok(Signal { source: Source::Big(Mutex::new(reader)), sequences, notes: Vec::new() });
        }
        let mut per: HashMap<String, Vec<(u64, u64, f32)>> = HashMap::new();
        let mut order: Vec<String> = Vec::new();
        // a wig section: (sequence, next start for fixedStep, step, span); None: bedGraph lines
        let mut wig: Option<(String, Option<u64>, u64, u64)> = None;
        let (mut skipped, mut first_bad) = (0usize, None);
        for (n, line) in BufReader::new(File::open(path)?).lines().enumerate() {
            let line = line?;
            let t = line.trim();
            if t.is_empty() || t.starts_with('#') || t.starts_with("track") || t.starts_with("browser") {
                continue;
            }
            let parsed: Result<Option<(String, u64, u64, f32)>, &str> = if t.starts_with("fixedStep") || t.starts_with("variableStep") {
                let field = |k: &str| t.split_whitespace().find_map(|w| w.strip_prefix(k).map(str::to_string));
                let num = |k: &str| field(k).map(|v| v.parse::<u64>().map_err(|_| "a number in the step line is not valid"));
                match (field("chrom="), num("start="), num("step="), num("span=")) {
                    (None, ..) => Err("no chrom= in the step line"),
                    (Some(_), Some(Err(why)), ..) | (Some(_), _, Some(Err(why)), _) | (Some(_), .., Some(Err(why))) => Err(why),
                    (Some(chrom), start, step, span) => {
                        let fixed = t.starts_with("fixedStep");
                        let step = step.and_then(Result::ok).unwrap_or(1).max(1);
                        let span = span.and_then(Result::ok).unwrap_or(if fixed { step } else { 1 }).max(1);
                        let start = start.and_then(Result::ok);
                        if fixed && start.is_none() {
                            Err("fixedStep without start=")
                        } else {
                            wig = Some((chrom, start.map(|s| s.saturating_sub(1)), step, span)); // wig counts from 1
                            Ok(None)
                        }
                    }
                }
            } else if let Some((chrom, next, step, span)) = wig.as_mut() {
                let f: Vec<&str> = t.split_whitespace().collect();
                match next {
                    Some(at) => match f.first().and_then(|v| v.parse::<f32>().ok()) {
                        Some(v) => {
                            let s = *at;
                            *at += *step;
                            Ok(Some((chrom.clone(), s, s + *span, v)))
                        }
                        None => Err("value is not a number"),
                    },
                    None => match (f.first().and_then(|p| p.parse::<u64>().ok()), f.get(1).and_then(|v| v.parse::<f32>().ok())) {
                        (Some(p), Some(v)) if p >= 1 => Ok(Some((chrom.clone(), p - 1, p - 1 + *span, v))),
                        _ => Err("not a position and a value"),
                    },
                }
            } else {
                let f: Vec<&str> = t.split_whitespace().collect();
                match (f.len() >= 4, f.get(1).and_then(|v| v.parse::<u64>().ok()), f.get(2).and_then(|v| v.parse::<u64>().ok()),
                       f.get(3).and_then(|v| v.parse::<f32>().ok())) {
                    (false, ..) => Err("fewer than 4 columns"),
                    (true, Some(s), Some(e), Some(v)) if e >= s => Ok(Some((f[0].to_string(), s, e, v))),
                    (true, Some(_), Some(_), Some(_)) => Err("end before start"),
                    _ => Err("start, end or value is not a number"),
                }
            };
            match parsed {
                Ok(Some((seq, s, e, v))) => {
                    if !per.contains_key(&seq) {
                        order.push(seq.clone());
                    }
                    per.entry(seq).or_default().push((s, e, v));
                }
                Ok(None) => {}
                Err(why) => {
                    skipped += 1;
                    first_bad.get_or_insert((n + 1, why));
                }
            }
        }
        let notes = crate::tracks::skipped_notes(path, "bedGraph or wig", skipped, first_bad, per.is_empty())?;
        let sequences = order.iter().map(|s| (s.clone(), per[s].iter().map(|iv| iv.1).max().unwrap_or(0))).collect();
        let mems = per.into_iter().map(|(seq, list)| (seq, Mem::new(list))).collect();
        Ok(Signal { source: Source::Mem(mems), sequences, notes })
    }

    /// (name, length) of each sequence: a bigWig states it, a text file reaches as far as its last value.
    pub fn sequences(&self) -> &[(String, u64)] {
        &self.sequences
    }

    /// Read from a bigWig (else from a bedGraph or wig text file).
    pub fn is_bigwig(&self) -> bool {
        matches!(self.source, Source::Big(_))
    }

    /// `n` bins over [start, end) of sequence `seq`: (mean over the bases with a value, highest value) each, NaN where none.
    pub fn bins(&self, seq: &str, start: u64, end: u64, n: usize) -> io::Result<Vec<(f64, f64)>> {
        let accs = match &self.source {
            Source::Mem(mems) => match mems.get(seq) {
                Some(m) => m.bins(start, end, n),
                None => vec![Acc::default(); n.max(1)],
            },
            Source::Big(reader) => {
                let Some(&(_, length)) = self.sequences.iter().find(|(name, _)| name == seq) else {
                    return Ok(vec![(f64::NAN, f64::NAN); n.max(1)]);
                };
                let mut r = reader.lock().unwrap_or_else(|e| e.into_inner());
                let (s, e) = (start.min(length) as u32, end.min(length) as u32);
                let mut bins = Bins::new(start, end, n);
                let half = bins.width / 2.0;
                let zoom = r.info().zoom_headers.iter().map(|z| z.reduction_level).filter(|&l| l as f64 <= half).max();
                let err = |e: &dyn std::fmt::Display| invalid(format!("{seq}:{start}-{end}: {e}"));
                if s < e {
                    match zoom {
                        Some(level) => {
                            for rec in r.get_zoom_interval(seq, s, e, level).map_err(|x| err(&x))? {
                                let rec = rec.map_err(|x| err(&x))?;
                                let m = &rec.summary;
                                bins.add(rec.start as u64, rec.end as u64,
                                         &Acc { sum: m.sum, covered: m.bases_covered as f64, max: m.max_val, min: m.min_val });
                            }
                        }
                        None => {
                            for v in r.get_interval(seq, s, e).map_err(|x| err(&x))? {
                                let v = v.map_err(|x| err(&x))?;
                                bins.add(v.start as u64, v.end as u64, &value(v.start as u64, v.end as u64, v.value as f64));
                            }
                        }
                    }
                }
                bins.out
            }
        };
        Ok(accs.iter().map(Acc::result).collect())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn write(dir: &Path, name: &str, text: &str) -> std::path::PathBuf {
        let p = dir.join(name);
        File::create(&p).unwrap().write_all(text.as_bytes()).unwrap();
        p
    }

    #[test]
    fn bedgraph_and_wig_bins_by_hand() {
        let dir = std::env::temp_dir().join(format!("opendna_signals_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        // c: [0,100) = 1, [100,200) = 3, nothing in [200,500), [500,600) = 2
        let bg = write(&dir, "a.bedGraph", "track type=bedGraph\nc\t0\t100\t1\nc\t100\t200\t3\nc\t500\t600\t2\nc\tx\t1\t1\n");
        let s = Signal::open(&bg).unwrap();
        assert_eq!(s.sequences(), &[("c".to_string(), 600)]);
        assert_eq!(s.notes, vec!["1 line was skipped as not valid bedGraph or wig (first: line 5, start, end or value is not a number).".to_string()]);
        // 6 bins of 100: means 1, 3, none x3, 2; the max of each the same
        let b = s.bins("c", 0, 600, 6).unwrap();
        assert_eq!((b[0], b[1], b[5]), ((1.0, 1.0), (3.0, 3.0), (2.0, 2.0)));
        assert!(b[2].0.is_nan() && b[4].1.is_nan());
        // one bin over [50, 150): 50 bases at 1 and 50 at 3 -> mean 2, max 3; one over [150, 550): 50 at 3, 50 at 2 -> 2.5
        assert_eq!(s.bins("c", 50, 150, 1).unwrap(), vec![(2.0, 3.0)]);
        assert_eq!(s.bins("c", 150, 550, 1).unwrap(), vec![(2.5, 3.0)]);
        // the whole of it in 1 bin goes through the zoom levels (64 bp blocks): (100 + 300 + 200) / 300 bases = 2
        assert_eq!(s.bins("c", 0, 600, 1).unwrap(), vec![(2.0, 3.0)]);
        assert!(s.bins("other", 0, 10, 2).unwrap()[0].0.is_nan());
        // wig: fixedStep from 1 (0 here) every 10 with span 5, then variableStep at 101 (100) span 10
        let w = write(&dir, "a.wig", "fixedStep chrom=c start=1 step=10 span=5\n4\n8\nvariableStep chrom=c span=10\n101 6\n");
        let s = Signal::open(&w).unwrap();
        let b = s.bins("c", 0, 20, 4).unwrap();
        assert!(b[0] == (4.0, 4.0) && b[1].0.is_nan() && b[2] == (8.0, 8.0) && b[3].0.is_nan());
        assert_eq!(s.bins("c", 100, 110, 1).unwrap(), vec![(6.0, 6.0)]);
        assert!(Signal::open(&write(&dir, "bad.bedGraph", "c\t1\n")).is_err());
        assert!(!is_bigwig(&bg));
    }

    #[test]
    fn bigwig_and_bigbed_of_the_demo_hub() {
        // demo/hub/pBAD30_eilR/plus.bw: [0,2000) = 10, [2000,4000) = 30, [4000,5936) = 20 (written by pybigtools)
        let hub = Path::new(env!("CARGO_MANIFEST_DIR")).join("../demo/hub/pBAD30_eilR");
        let s = Signal::open(&hub.join("plus.bw")).unwrap();
        assert!(s.is_bigwig() && s.notes.is_empty());
        assert_eq!(s.sequences(), &[("pBAD30_eilR".to_string(), 5936)]);
        assert_eq!(s.bins("pBAD30_eilR", 1000, 3000, 2).unwrap(), vec![(10.0, 10.0), (30.0, 30.0)]);
        // the whole of it in one bin: (2000 x 10 + 2000 x 30 + 1936 x 20) / 5936 = 20 exactly, highest 30
        assert_eq!(s.bins("pBAD30_eilR", 0, 5936, 1).unwrap(), vec![(20.0, 30.0)]);
        // [1990, 2010) in 2 bins: 10 bases at 10 then 10 at 30; a sequence the file does not have: no value
        assert_eq!(s.bins("pBAD30_eilR", 1990, 2010, 2).unwrap(), vec![(10.0, 10.0), (30.0, 30.0)]);
        assert!(s.bins("other", 0, 10, 1).unwrap()[0].0.is_nan());
        // sites.bb: three BED 9 entries, read as the BED lines they stand for (name, strand, itemRgb colour)
        let t = crate::tracks::Track::open(&hub.join("sites.bb")).unwrap();
        let (items, more) = t.query("pBAD30_eilR", 0, 5936, 10);
        let got: Vec<_> = items.iter().map(|i| (i.start, i.end, i.strand, i.name.as_str(), i.color.as_deref())).collect();
        assert!(!more);
        assert_eq!(got, vec![(100, 400, 1, "site_a", Some("255,0,0")), (2500, 2600, -1, "site_b", Some("0,0,255")),
                             (5000, 5100, 1, "site_c", Some("0,128,0"))]);
    }
}
