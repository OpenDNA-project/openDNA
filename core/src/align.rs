//! Quick alignment ("BLAST against the open file"): where does a query sequence best match
//! the template, on either strand, with mismatches and gaps?
//!
//! Seed and extend, per strand:
//!   1. Index the query's k-mers (k = 11; 9 or less for queries under 200 bases).
//!   2. Scan the template once; every shared k-mer is a "seed", sitting on a diagonal
//!      (diagonal = template position - query position; an ungapped match stays on one).
//!   3. Seeds on nearby diagonals are grouped into clusters. The clusters with the most
//!      seeds get a banded Smith-Waterman (local alignment with affine gaps) over the
//!      query and the few diagonals around the cluster, which finds mismatches, insertions
//!      and deletions.
//!   3b. A cluster can hold SEVERAL separate matching pieces of the query (for example the two
//!      ends of a query whose middle is a long stretch that matches nothing: both ends sit on
//!      the same diagonal). After the best local alignment of a cluster is found, the query
//!      to its left and to its right are aligned again, and so on, so that every piece scoring
//!      above the reporting threshold is reported on its own (`pieces_of_cluster`).
//!   4. Hits are sorted by score and the best ones are kept, dropping a hit only when it is a
//!      repeat of a better one (same template place AND same part of the query).
//! Circular templates are scanned with their first bases appended, so hits can cross the origin.
//!
//! This makes a query that is a patchwork (a structural variant: pieces in another order,
//! an inverted piece, a deleted or inserted stretch) come out as one row per piece, each with
//! its own strand, query range and template range.
//!
//! shortcut: a gap longer than the band (~32 bases) splits a hit in two instead of being
//! bridged; very repetitive queries skip their most repeated k-mers.

use std::collections::HashMap;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::{Duration, Instant};

use crate::primers::{clean, revcomp};

// Scores: +2 match, -3 mismatch, a gap of n bases costs 5 + 2 * (n - 1).
const MATCH: i32 = 2;
const MISMATCH: i32 = -3;
const GAP_OPEN: i32 = -5;
const GAP_EXTEND: i32 = -2;
/// Diagonals searched on each side of a seed cluster.
const BAND: i64 = 16;
/// Separate pieces reported from one cluster of diagonals (keeps a repetitive query fast).
const MAX_PIECES_PER_CLUSTER: usize = 64;
/// A hit is never allowed to cross a stretch where its running score falls this many points below
/// its best earlier value (roughly 14 mismatches in a row, or 40 bases of unrelated DNA: unrelated DNA costs about
/// 0.7 points per base even when gaps are used to hide the differences). Without this rule a
/// 200 base unmatched insert between two 100 base matches would be "bridged" into one low-identity
/// hit (400 points of matches minus ~150 for the middle still beats one piece alone), and the
/// pieces could not be told apart. A short unmatched stretch, a SNP or a small indel costs far
/// less than this and is still bridged.
const MAX_DRAWDOWN: i32 = 40;
/// Window (in alignment columns) in which fewer than half the columns may match (see
/// `smith_waterman`): such a stretch is cut out of a hit.
const LOW_WINDOW: usize = 40;
/// An alignment must start and end with this many exact matches in a row (see the trimming in
/// `smith_waterman`).
const TRIM_RUN: usize = 12;
/// Exact matches needed beyond a lone mismatch for the trimming to keep it (chance gives 5 in a row once in 1000).
const KEEP_RUN: usize = 5;
/// Safety stop for degenerate input (e.g. poly-A query on a poly-A template).
const MAX_SEEDS: usize = 5_000_000;
/// A query k-mer present more often than this in the query is not used as a seed.
const MAX_KMER_COPIES: usize = 100;
/// "Minus infinity" that cannot overflow when a gap penalty is added.
const NEG: i32 = i32::MIN / 2;

/// One local alignment of the query on the template. Positions are 0-based on the top
/// strand; `end` (and the positions inside) may pass the sequence length for a hit that
/// crosses the origin of a circular template.
#[derive(Clone, Debug, PartialEq)]
pub struct Hit {
    pub start: usize,
    pub end: usize, // exclusive
    pub forward: bool, // false: the reverse complement of the query matches the top strand
    pub score: i32,
    pub identity: f32, // matching bases / alignment columns, 0..1
    pub qstart: usize, // part of the query aligned, in the query as given
    pub qend: usize,
    /// (template position, query base) like the mismatches of a BAM read.
    pub mismatches: Vec<(usize, u8)>,
    /// (template position the extra bases sit before, number of bases)
    pub insertions: Vec<(usize, usize)>,
    /// (first, last + 1) template positions missing from the query
    pub deletions: Vec<(usize, usize)>,
}

/// 2-bit code of a base; None for N and other letters (they break k-mers).
fn code(b: u8) -> Option<u32> {
    match b {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// Seeds as (diagonal, template position) for one strand of the query, the number of query k-mers
/// left out as too repeated, and the template position where the seed limit stopped the scan (if it did).
fn seeds(text: &[u8], query: &[u8], k: usize) -> (Vec<(i64, usize)>, usize, Option<usize>) {
    // k-mer (2 bits per base) -> positions in the query.
    let mut index: HashMap<u32, Vec<usize>> = HashMap::new();
    let mask: u32 = if k >= 16 { u32::MAX } else { (1u32 << (2 * k)) - 1 };
    let (mut word, mut run) = (0u32, 0usize); // `run` = valid bases in a row so far
    for (i, &b) in query.iter().enumerate() {
        match code(b) {
            Some(c) => { word = ((word << 2) | c) & mask; run += 1; }
            None => run = 0,
        }
        if run >= k {
            index.entry(word).or_default().push(i + 1 - k);
        }
    }
    // Very repeated k-mers would flood the seed list: leave them out (and say so).
    let distinct = index.len();
    index.retain(|_, v| v.len() <= MAX_KMER_COPIES);
    let repeated = distinct - index.len();
    // Bit filter: one bit per possible k-mer, so almost every template position is
    // rejected by a single memory read instead of a hash lookup.
    let mut filter = vec![0u64; ((mask as usize) >> 6) + 1];
    for &w in index.keys() {
        filter[(w >> 6) as usize] |= 1 << (w & 63);
    }
    let mut out = Vec::new();
    let (mut word, mut run) = (0u32, 0usize);
    for (i, &b) in text.iter().enumerate() {
        match code(b) {
            Some(c) => { word = ((word << 2) | c) & mask; run += 1; }
            None => run = 0,
        }
        if run >= k && filter[(word >> 6) as usize] & (1 << (word & 63)) != 0 {
            let t = i + 1 - k;
            for &q in &index[&word] {
                out.push((t as i64 - q as i64, t));
            }
            if out.len() > MAX_SEEDS {
                return (out, repeated, Some(t));
            }
        }
    }
    (out, repeated, None)
}

/// `work(i)` for every i in 0..count, spread over `threads` threads (0 = every CPU core); the results come back
/// in index order, so the answer is the same as on one thread. Tasks are started in index order; one not yet
/// started when `deadline` passes is not run (None).
fn on_cores<T: Send>(threads: usize, count: usize, deadline: Option<Instant>, work: impl Fn(usize) -> T + Sync) -> Vec<Option<T>> {
    let next = AtomicUsize::new(0);
    let all = std::thread::available_parallelism().map_or(1, |n| n.get());
    let threads = if threads == 0 { all } else { threads }.min(count.max(1));
    let done: Vec<Vec<(usize, T)>> = std::thread::scope(|s| {
        let workers: Vec<_> = (0..threads)
            .map(|_| {
                s.spawn(|| {
                    let mut mine = Vec::new();
                    loop {
                        let i = next.fetch_add(1, Ordering::Relaxed);
                        if i >= count || deadline.is_some_and(|d| Instant::now() >= d) {
                            return mine;
                        }
                        mine.push((i, work(i)));
                    }
                })
            })
            .collect();
        workers.into_iter().map(|w| w.join().expect("an alignment thread panicked")).collect()
    });
    let mut out: Vec<Option<T>> = (0..count).map(|_| None).collect();
    for (i, t) in done.into_iter().flatten() {
        out[i] = Some(t);
    }
    out
}

/// Groups seeds on nearby diagonals: (lowest diagonal, highest diagonal, number of seeds), most seeds first.
fn clusters(mut seeds: Vec<(i64, usize)>, min_seeds: usize) -> Vec<(i64, i64, usize)> {
    seeds.sort_unstable();
    let mut out: Vec<(i64, i64, usize)> = Vec::new();
    let mut cur: Option<(i64, i64, usize)> = None;
    for (d, _) in seeds {
        cur = match cur {
            // Same cluster while diagonals stay close (indels shift them a little).
            Some((lo, hi, n)) if d - hi <= BAND && d - lo <= 3 * BAND => Some((lo, d, n + 1)),
            Some(c) => { out.push(c); Some((d, d, 1)) }
            None => Some((d, d, 1)),
        };
    }
    out.extend(cur);
    out.retain(|c| c.2 >= min_seeds);
    out.sort_unstable_by(|a, b| b.2.cmp(&a.2));
    out
}

/// What one column of an alignment is, kept while tracing back.
enum Op {
    Match,
    Mismatch(usize, u8), // template position, query base
    Insert(usize),       // template position it sits before
    Delete(usize),       // template position
}

/// Banded Smith-Waterman of the whole `query` against `text` along diagonals lo..=hi
/// (diagonal = template position - query position). Returns the best local alignment.
///
/// Returns the hit and, when the alignment crosses a stretch that scores badly (see
/// `MAX_DRAWDOWN`), `Some(q)`: a position of the query (counted in `query`) inside that stretch,
/// where the caller should cut the query in two and align the halves separately.
fn smith_waterman(text: &[u8], query: &[u8], lo: i64, hi: i64) -> Option<(Hit, Option<usize>)> {
    let m = query.len();
    let width = (hi - lo + 1) as usize;
    // Cell (i, j): i query bases and t = i + lo + j template bases used so far.
    // h = best score ending here; e = ending in a deletion; f = ending in an insertion.
    // `trace` remembers how each cell was reached (bits 0-1: where h came from, bit 2:
    // the deletion extends a previous one, bit 3: the insertion does).
    let mut trace = vec![0u8; (m + 1) * width];
    let mut h_prev = vec![0i32; width];
    let mut f_prev = vec![NEG; width];
    let mut h_cur = vec![0i32; width];
    let mut f_cur = vec![NEG; width];
    let (mut best, mut best_i, mut best_j) = (0i32, 0usize, 0usize);
    for i in 1..=m {
        let mut e = NEG; // deletion score of the previous cell in this row
        for j in 0..width {
            let t = i as i64 + lo + j as i64;
            let cell = i * width + j;
            if t < 0 || t > text.len() as i64 {
                h_cur[j] = 0; f_cur[j] = NEG; e = NEG;
                continue;
            }
            let mut tb = 0u8;
            // Deletion: uses template base t-1, from the cell to the left (j - 1, same row).
            let mut e_here = NEG;
            if j > 0 {
                let (open, ext) = (h_cur[j - 1] + GAP_OPEN, e + GAP_EXTEND);
                e_here = open.max(ext);
                if ext > open { tb |= 4; }
            }
            // Insertion: uses query base i-1, from the cell above-right (i - 1, j + 1).
            let mut f_here = NEG;
            if j + 1 < width {
                let (open, ext) = (h_prev[j + 1] + GAP_OPEN, f_prev[j + 1] + GAP_EXTEND);
                f_here = open.max(ext);
                if ext > open { tb |= 8; }
            }
            // Diagonal: query base i-1 against template base t-1.
            let mut h = 0;
            if t >= 1 {
                let (q, b) = (query[i - 1], text[t as usize - 1]);
                let s = h_prev[j] + if q == b && code(q).is_some() { MATCH } else { MISMATCH };
                if s > h { h = s; tb |= 1; }
            }
            if e_here > h { h = e_here; tb = (tb & !3) | 2; }
            if f_here > h { h = f_here; tb = (tb & !3) | 3; }
            h_cur[j] = h; f_cur[j] = f_here; e = e_here;
            trace[cell] = tb;
            if h > best { best = h; best_i = i; best_j = j; }
        }
        std::mem::swap(&mut h_prev, &mut h_cur);
        std::mem::swap(&mut f_prev, &mut f_cur);
    }
    if best <= 0 {
        return None;
    }
    // Trace back from the best cell to where the score was 0, collecting columns end-first.
    let (mut i, mut j) = (best_i, best_j);
    let end = best_i as i64 + lo + best_j as i64; // template bases used = exclusive end
    let mut state = 0; // 0 = in h, 1 = in a deletion, 2 = in an insertion
    let mut ops: Vec<Op> = Vec::new();
    loop {
        let t = (i as i64 + lo + j as i64) as usize;
        let tb = trace[i * width + j];
        match state {
            0 => match tb & 3 {
                1 => {
                    let (q, b) = (query[i - 1], text[t - 1]);
                    ops.push(if q == b && code(q).is_some() { Op::Match } else { Op::Mismatch(t - 1, q) });
                    i -= 1;
                }
                2 => state = 1,
                3 => state = 2,
                _ => break,
            },
            1 => {
                ops.push(Op::Delete(t - 1));
                state = if tb & 4 != 0 { 1 } else { 0 };
                j -= 1;
            }
            _ => {
                ops.push(Op::Insert(t));
                state = if tb & 8 != 0 { 2 } else { 0 };
                i -= 1;
                j += 1;
            }
        }
        if i == 0 && state == 0 {
            break;
        }
    }
    let start = (i as i64 + lo + j as i64) as usize;
    ops.reverse();
    // Trim ragged ends. The best-scoring alignment may reach a few bases into an unmatched
    // neighbour through chance matches and a small gap (e.g. 15 bases of random DNA next to a
    // real piece), which makes the reported coordinates of a piece fuzzy. So the alignment must
    // begin and end with TRIM_RUN exact matches in a row; whatever lies outside the first and
    // the last such run is cut off. (Alignments with no such run, or fewer than 40 columns,
    // are left alone.) The tuple below shadows the earlier names with ones we may change:
    // `let` can re-declare a name, here to make it `mut` (changeable).
    let (mut start, mut end, mut qstart, mut qend) = (start, end as usize, i, best_i);
    if ops.len() >= 40 {
        // `is_run(k)`: do the TRIM_RUN columns from k on all match? (A closure: a small
        // inline function that can read the surrounding variable `ops`.)
        let is_run = |k: usize| ops[k..k + TRIM_RUN].iter().all(|o| matches!(o, Op::Match));
        let first = (0..=ops.len() - TRIM_RUN).find(|&k| is_run(k));
        let last = (0..=ops.len() - TRIM_RUN).rev().find(|&k| is_run(k)).map(|k| k + TRIM_RUN);
        if let (Some(mut f), Some(mut l)) = (first, last) {
            // A lone mismatch near an end is a real difference (a SNP 8 bases from the end of a read), not a
            // ragged end: keep it when KEEP_RUN exact matches lie beyond it, and the matches before those.
            let exact = |o: &Op| matches!(o, Op::Match);
            while f > KEEP_RUN && matches!(ops[f - 1], Op::Mismatch(..)) && ops[f - 1 - KEEP_RUN..f - 1].iter().all(exact) {
                f -= 1 + KEEP_RUN;
                while f > 0 && exact(&ops[f - 1]) { f -= 1; }
            }
            while l + KEEP_RUN < ops.len() && matches!(ops[l], Op::Mismatch(..)) && ops[l + 1..l + 1 + KEEP_RUN].iter().all(exact) {
                l += 1 + KEEP_RUN;
                while l < ops.len() && exact(&ops[l]) { l += 1; }
            }
            // A ragged end is a few columns. A cut of over a quarter of them is a short query with several
            // mismatches near its ends (a 45-mer with 3), which would be reported as a stub or not at all: keep it whole.
            if (l - f) * 4 < ops.len() * 3 {
                (f, l) = (0, ops.len());
            }
            // Every cut column used query and/or template bases; give them back.
            for op in &ops[..f] {
                match op {
                    Op::Match | Op::Mismatch(..) => { start += 1; qstart += 1; }
                    Op::Insert(_) => qstart += 1,
                    Op::Delete(_) => start += 1,
                }
            }
            for op in &ops[l..] {
                match op {
                    Op::Match | Op::Mismatch(..) => { end -= 1; qend -= 1; }
                    Op::Insert(_) => qend -= 1,
                    Op::Delete(_) => end -= 1,
                }
            }
            ops.truncate(l); // keep columns 0..l
            ops.drain(..f); // remove columns 0..f
        }
    }
    // Cutting can stop one base short when a gap sat just before the last exact bases (a gap
    // inside a run of equal bases can slide either way, e.g. "GGT" against "GGTT"). So grow the
    // alignment again, without gaps, over query bases that equal the template bases next to it.
    let col = |q: usize, t: usize| if query[q] == text[t] && code(query[q]).is_some() { Op::Match } else { Op::Mismatch(t, query[q]) };
    let one_mismatch = |ops: &[Op]| ops.iter().filter(|o| matches!(o, Op::Mismatch(..))).count() == 1;
    while qend < m && end < text.len() && matches!(col(qend, end), Op::Match) {
        ops.push(Op::Match); end += 1; qend += 1;
    }
    // A lone SNP in the query's last few bases (too few exact bases beyond it for the trimming to keep it):
    // take the query to its end and report the SNP, rather than end the hit short of it.
    if qend < m && m - qend <= KEEP_RUN + 1 && end + (m - qend) <= text.len() {
        let tail: Vec<Op> = (0..m - qend).map(|k| col(qend + k, end + k)).collect();
        if one_mismatch(&tail) {
            (end, qend) = (end + tail.len(), m);
            ops.extend(tail);
        }
    }
    let mut head = Vec::new(); // columns added before the start, nearest first
    while qstart > 0 && start > 0 && matches!(col(qstart - 1, start - 1), Op::Match) {
        head.push(Op::Match); start -= 1; qstart -= 1;
    }
    if qstart > 0 && qstart <= KEEP_RUN + 1 && start >= qstart {
        let lead: Vec<Op> = (0..qstart).map(|k| col(k, start - qstart + k)).collect();
        if one_mismatch(&lead) {
            start -= qstart;
            qstart = 0;
            head.extend(lead.into_iter().rev());
        }
    }
    head.reverse();
    ops.splice(0..0, head);
    let (mut mismatches, mut insertions, mut deletions) = (Vec::new(), Vec::new(), Vec::new());
    let (mut matches, mut columns) = (0usize, 0usize);
    // Running score along the alignment, its best value so far, the query position reached, and
    // what the previous column was (0 = match/mismatch, 1 = deletion, 2 = insertion: a gap column
    // right after a gap column of the same kind is a cheap "extend", otherwise a dear "open").
    let (mut running, mut peak, mut qpos, mut prev) = (0i32, 0i32, qstart, 0u8);
    let mut split: Option<usize> = None;
    let mut after: Vec<usize> = Vec::with_capacity(ops.len()); // query position after each column
    for op in &ops {
        columns += 1;
        match *op {
            Op::Match => { matches += 1; running += MATCH; qpos += 1; prev = 0; }
            Op::Mismatch(t, q) => { mismatches.push((t, q)); running += MISMATCH; qpos += 1; prev = 0; }
            // Consecutive inserted bases / deleted positions merge into one event.
            Op::Insert(t) => {
                match insertions.last_mut() {
                    Some((pos, n)) if *pos == t => *n += 1,
                    _ => insertions.push((t, 1)),
                }
                running += if prev == 2 { GAP_EXTEND } else { GAP_OPEN };
                qpos += 1;
                prev = 2;
            }
            Op::Delete(t) => {
                match deletions.last_mut() {
                    Some((_, e)) if *e == t => *e += 1,
                    _ => deletions.push((t, t + 1)),
                }
                running += if prev == 1 { GAP_EXTEND } else { GAP_OPEN };
                prev = 1;
            }
        }
        after.push(qpos);
        peak = peak.max(running);
        // The first place where the score has fallen too far below its best: remember it.
        if split.is_none() && peak - running > MAX_DRAWDOWN {
            split = Some(qpos);
        }
    }
    // A second way to spot an unmatched stretch: any run of LOW_WINDOW columns of which fewer than
    // half match is not homology, however cheaply gaps let the alignment cross it. The cut goes
    // in the middle of such a window, which then lies in the unmatched stretch.
    if split.is_none() && ops.len() >= LOW_WINDOW {
        // `is_match[k]` = 1 when column k is a match; `ones` is a running count over the window.
        let is_match: Vec<usize> = ops.iter().map(|o| matches!(o, Op::Match) as usize).collect();
        let mut ones: usize = is_match[..LOW_WINDOW].iter().sum();
        for k in 0..=ops.len() - LOW_WINDOW {
            if ones * 2 < LOW_WINDOW {
                split = Some(after[k + LOW_WINDOW / 2]);
                break;
            }
            if k + LOW_WINDOW < ops.len() {
                ones += is_match[k + LOW_WINDOW];
                ones -= is_match[k];
            }
        }
    }
    let hit = Hit {
        start,
        end,
        forward: true,
        score: running, // the score of the (possibly trimmed) alignment
        identity: matches as f32 / columns.max(1) as f32,
        qstart,
        qend,
        mismatches,
        insertions,
        deletions,
    };
    Some((hit, split))
}

/// Every separate piece of the query that matches inside one cluster of diagonals.
///
/// `lo..=hi` are the diagonals of the band. The idea is "find the best, then look left and
/// right of it": align the whole query, keep the best local alignment, then do the same
/// for the part of the query before it and the part after it, recursively. A piece is only
/// kept when its score reaches `min_score`; when the best piece of a stretch is already below
/// that, nothing in the stretch can be better, so that branch stops.
///
/// Two pieces on the same diagonal separated by a long run of mismatches (a 200 base insert
/// in the query, say) stay two pieces: bridging 200 mismatching bases costs far more than the
/// pieces earn, so the best alignment is one of them and the other is found in the next round.
/// The bool is true when MAX_PIECES_PER_CLUSTER stopped the search.
fn pieces_of_cluster(text: &[u8], query: &[u8], lo: i64, hi: i64, min_score: i32) -> (Vec<Hit>, bool) {
    let mut found: Vec<Hit> = Vec::new();
    // A to-do list of query ranges (from, to) still to be searched. `vec![(0, query.len())]`
    // builds a Vec holding that one tuple; `pop()` takes the last item (or None when empty).
    let mut todo = vec![(0usize, query.len())];
    while let Some((from, to)) = todo.pop() {
        if found.len() >= MAX_PIECES_PER_CLUSTER {
            return (found, true);
        }
        if to - from < 4 {
            continue; // nothing worth aligning in so little query
        }
        // Aligning the slice query[from..to] instead of the whole query shifts every diagonal by
        // `from` (template position - query position uses the position inside the slice), so the
        // band is shifted too. The positions inside the returned hit are template positions and
        // stay valid as they are; only the query positions need `from` added back.
        let shift = from as i64;
        let Some((mut hit, split)) = smith_waterman(text, &query[from..to], lo + shift, hi + shift) else { continue };
        if hit.score < min_score {
            continue;
        }
        // The alignment crossed a long stretch that matches nothing (see MAX_DRAWDOWN): it is not
        // one piece. Cut the query inside that stretch and align each half on its own; the halves
        // are shorter than `from..to`, so this always ends.
        if let Some(cut) = split.map(|q| from + q).filter(|&c| c > from && c < to) {
            todo.push((from, cut));
            todo.push((cut, to));
            continue;
        }
        hit.qstart += from;
        hit.qend += from;
        todo.push((from, hit.qstart)); // the query left of this piece
        todo.push((hit.qend, to)); // and right of it
        found.push(hit);
    }
    (found, false)
}

/// True when `a` and `b` are the same piece found twice: they overlap on the template (a circular
/// template: also one lap before or after) AND on the query. Overlapping on the template only is
/// kept: it is a duplicated piece of the query (both copies align to the same place). Overlapping
/// on the query only is kept too: one stretch of the query matching two places of the template
/// (a repeat).
fn redundant(a: &Hit, b: &Hit, n: usize) -> bool {
    let on_query = a.qstart < b.qend && b.qstart < a.qend;
    on_query && [0, n as i64, -(n as i64)].iter().any(|&s| (a.start as i64) < b.end as i64 + s && (b.start as i64 + s) < a.end as i64)
}

/// The best `max_hits` local alignments of `query` on `template` (either strand), best first.
/// A query that is a patchwork gives one hit per matching piece. `forward == false` hits are the
/// reverse complement of the query; their `qstart..qend` is still counted in the query as given.
/// The notes say what a limit left out (never silently); the caller shows them with the hits.
/// `time_limit` (seconds, 0 = none): candidate places not yet aligned when it runs out are left out, the
/// ones with fewest seeds. The places are aligned on `threads` threads (0 = every CPU core).
pub fn align(template: &[u8], query: &[u8], circular: bool, max_hits: usize, time_limit: f64, threads: usize) -> (Vec<Hit>, Vec<String>) {
    let deadline = (time_limit > 0.0).then(|| Instant::now() + Duration::from_secs_f64(time_limit));
    let (template, query) = (clean(template), clean(query));
    let (n, m) = (template.len(), query.len());
    let mut notes = Vec::new();
    if n == 0 || m < 4 || max_hits == 0 {
        return (Vec::new(), notes);
    }
    if circular && m > n {
        notes.push(format!("The query ({m} bp) is longer than the circular template ({n} bp): it is matched at most once round the circle."));
    }
    // Circular: append the start of the template so hits can run through the origin.
    let mut text = template;
    if circular {
        text.extend_from_within(..n.min(m + 2 * BAND as usize));
    }
    // Short queries get short k-mers: a 60-mer with 3 mismatches often has no two shared 11-mers.
    let k = if m >= 200 { 11 } else { (m / 4).clamp(4, 9) };
    let min_score = (m as i32 * 6 / 5).min(40); // random matches score ~20-30 on long queries
    // Every candidate place of both strands: (strand, lowest diagonal, highest diagonal, seeds).
    let mut places: Vec<(bool, i64, i64, usize)> = Vec::new();
    let strands = [query.clone(), revcomp(&query)];
    for forward in [true, false] {
        let q = &strands[!forward as usize];
        let strand = if forward { "forward" } else { "reverse" };
        let (seeds, repeated, stopped) = seeds(&text, q, k);
        if forward && repeated > 0 {
            notes.push(format!("{repeated} query {k}-mers occur more than {MAX_KMER_COPIES} times in the query and were not used as seeds: a match made only of that repeat can be missed."));
        }
        if let Some(t) = stopped {
            notes.push(format!("The seed limit ({MAX_SEEDS}) was reached on the {strand} strand at template base {}: the template after it was not searched on that strand (the query is too repetitive for this template).", t % n + 1));
        }
        places.extend(clusters(seeds, 2).into_iter().map(|(lo, hi, s)| (forward, lo, hi, s)));
    }
    // Most seeds first (both strands mixed), so a time limit leaves out the weakest places. A stable sort keeps
    // each strand's own order among equals: the answer does not depend on the number of cores.
    places.sort_by(|a, b| b.3.cmp(&a.3));
    let done = on_cores(threads, places.len(), deadline, |i| {
        let (forward, lo, hi, _) = places[i];
        pieces_of_cluster(&text, &strands[!forward as usize], lo - BAND, hi + BAND, min_score)
    });
    let skipped = done.iter().filter(|d| d.is_none()).count();
    if skipped > 0 {
        notes.push(format!("The time limit ({time_limit} s) was reached: {skipped} of {} candidate places (those with fewest seeds) were not aligned, so hits may be missing. Raise or clear the time limit to align them all.", places.len()));
    }
    let mut hits: Vec<Hit> = Vec::new();
    let mut capped = 0;
    for (&(forward, ..), (pieces, full)) in places.iter().zip(done).filter_map(|(p, d)| d.map(|d| (p, d))) {
        capped += full as usize;
        for mut hit in pieces {
            if hit.start >= n {
                continue; // a repeat of a hit at the start of the template
            }
            hit.forward = forward;
            if !forward {
                // Report the query part in the query's own orientation.
                (hit.qstart, hit.qend) = (m - hit.qend, m - hit.qstart);
            }
            hits.push(hit);
        }
    }
    if capped > 0 {
        notes.push(format!("{capped} places hold more than {MAX_PIECES_PER_CLUSTER} matching pieces of the query: only the first {MAX_PIECES_PER_CLUSTER} of each were listed."));
    }
    hits.sort_by(|a, b| b.score.cmp(&a.score));
    let mut kept: Vec<Hit> = Vec::new();
    for h in hits {
        if kept.len() < max_hits && !kept.iter().any(|k| redundant(k, &h, if circular { n } else { 0 })) {
            kept.push(h);
        }
    }
    (kept, notes)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The hits alone (the notes are checked in `every_limit_reached_is_a_note`).
    fn align(t: &[u8], q: &[u8], circular: bool, max_hits: usize) -> Vec<Hit> {
        super::align(t, q, circular, max_hits, 0.0, 0).0
    }

    /// Reproducible pseudo-random DNA (no `rand` crate needed).
    fn random_dna(len: usize, mut seed: u64) -> Vec<u8> {
        (0..len)
            .map(|_| {
                seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                b"ACGT"[(seed >> 33) as usize % 4]
            })
            .collect()
    }

    #[test]
    fn exact_hit_forward_and_reverse() {
        let t = random_dna(2000, 1);
        let q = t[500..800].to_vec();
        let hits = align(&t, &q, false, 10);
        assert_eq!((hits[0].start, hits[0].end, hits[0].forward), (500, 800, true));
        assert_eq!(hits[0].identity, 1.0);
        assert_eq!((hits[0].qstart, hits[0].qend), (0, 300));
        let rc = revcomp(&q);
        let hits = align(&t, &rc, false, 10);
        assert_eq!((hits[0].start, hits[0].end, hits[0].forward), (500, 800, false));
        assert_eq!(hits.len(), 1);
    }

    #[test]
    fn mismatch_insertion_and_deletion() {
        let t = random_dna(3000, 2);
        let mut q = t[1000..1300].to_vec();
        // The template bases next to each gap differ from what is inserted/deleted, so no gap can slide:
        // the positions checked below are the only right answer.
        assert_eq!((t[1149], t[1150], t[1249], t[1254], t[1250], t[1255]), (b'G', b'A', b'G', b'A', b'T', b'G'));
        // Edit from the back so earlier positions stay valid: delete 5 bases, insert 3, mutate 1.
        q.drain(250..255);
        for (i, b) in b"TTT".iter().enumerate() { q.insert(150 + i, *b); }
        let wrong = if q[50] == b'A' { b'C' } else { b'A' };
        q[50] = wrong;
        let hits = align(&t, &q, false, 10);
        let h = &hits[0];
        assert!(h.forward && h.start == 1000 && h.end == 1300, "{h:?}");
        assert_eq!(h.mismatches, vec![(1050, wrong)]);
        assert_eq!(h.insertions.iter().map(|x| x.1).sum::<usize>(), 3);
        assert_eq!(h.deletions.iter().map(|x| x.1 - x.0).sum::<usize>(), 5);
        assert!(h.identity > 0.95);
        // --- exact places: 3 bases inserted before template 1150, template 1250..1255 deleted ---
        assert_eq!(hits.len(), 1);
        assert_eq!((h.qstart, h.qend), (0, 298)); // 300 - 5 + 3
        assert_eq!(h.insertions, vec![(1150, 3)]);
        assert_eq!(h.deletions, vec![(1250, 1255)]);
        // Same on the other strand, every difference reported in top-strand terms.
        let hits = align(&t, &revcomp(&q), false, 10);
        assert!(!hits[0].forward && hits[0].start == 1000 && hits[0].end == 1300);
        let r = &hits[0];
        assert_eq!((r.qstart, r.qend), (0, 298));
        assert_eq!((&r.mismatches, &r.insertions, &r.deletions), (&vec![(1050, wrong)], &vec![(1150, 3)], &vec![(1250, 1255)]));
    }

    #[test]
    fn a_mismatch_near_either_end_of_a_long_query_is_reported_not_trimmed_away() {
        let t = random_dna(3000, 12);
        for at in [0, 3, 6, 11, 288, 293, 296, 299] {
            let mut q = t[1000..1300].to_vec();
            let wrong = if q[at] == b'A' { b'C' } else { b'A' };
            q[at] = wrong;
            let h = &align(&t, &q, false, 10)[0];
            assert_eq!((h.start, h.end, h.qstart, h.qend), (1000, 1300, 0, 300), "mismatch at {at}");
            assert_eq!(h.mismatches, vec![(1000 + at, wrong)]);
        }
    }

    #[test]
    fn a_short_query_with_mismatches_near_both_ends_is_still_one_whole_hit() {
        let t = random_dna(3000, 13);
        let mut q = t[1000..1045].to_vec();
        for at in [9, 13, 35] {
            q[at] = if q[at] == b'A' { b'C' } else { b'A' };
        }
        let hits = align(&t, &q, false, 10);
        assert_eq!((hits[0].start, hits[0].end, hits[0].mismatches.len()), (1000, 1045, 3), "{hits:?}");
    }

    #[test]
    fn every_limit_reached_is_a_note() {
        let notes = |t: &[u8], q: &[u8], circular: bool| super::align(t, q, circular, 10, 0.0, 0).1.join(" ");
        let t = random_dna(1000, 30);
        assert_eq!(notes(&t, &t[200..500], true), "", "an ordinary query: no note");
        // A query longer than the circle.
        let mut q = t.clone();
        q.extend_from_slice(&t[..200]);
        assert!(notes(&t, &q, true).contains("longer than the circular template (1000 bp)"));
        // A query k-mer repeated more than 100 times in the query.
        let mut q = vec![b'A'; 300];
        q.extend_from_slice(&t[0..100]);
        assert!(notes(&t, &q, false).contains("more than 100 times in the query"));
        // 300 copies of a unit, each its own place: all aligned and found (a fixed cap of 200 places used to cut them), no note.
        let unit = random_dna(50, 31);
        let mut many = Vec::new();
        for i in 0..300 {
            many.extend_from_slice(&unit);
            many.extend(random_dna(100, 100 + i));
        }
        let (hits, n) = super::align(&many, &unit, false, 1000, 0.0, 0);
        assert_eq!((hits.len(), n.len()), (300, 0), "{n:?}");
        assert_eq!(super::align(&many, &unit, false, 1000, 0.0, 1).0, hits, "one thread or all: the same answer");
        // A time limit already over when the places are ready: none is aligned, and the note says how many were left out.
        let (hits, n) = super::align(&many, &unit, false, 1000, 1e-9, 0);
        assert!(hits.is_empty() && n[0].contains("were not aligned, so hits may be missing"), "{n:?}");
        // Seed limit: 96 copies of the 9-mer AAAAAAAAA in the query times 60,000 template places.
        let n = notes(&vec![b'A'; 60_000], &vec![b'A'; 104], false);
        assert!(n.contains("seed limit (5000000) was reached on the forward strand at template base"), "{n}");
    }

    #[test]
    fn no_hit_for_random_query() {
        let t = random_dna(5000, 4);
        let q = random_dna(500, 99);
        assert!(align(&t, &q, false, 10).is_empty());
    }

    // ---- a query that hangs off the end of the template -------------------------------------------------
    // (`random_dna(n, seed)` is a Vec<u8>; `&t[a..b]` is a slice, `.to_vec()` copies it, `extend_from_slice` appends.)

    #[test]
    fn query_hanging_off_the_template_end_aligns_only_the_overlap() {
        let t = random_dna(1000, 11);
        let overhang = random_dna(60, 12); // bases that are nowhere in the template
        // The query starts with the template's last 40 bases, then hangs off the end.
        let mut q = t[960..].to_vec();
        q.extend_from_slice(&overhang);
        let hits = align(&t, &q, false, 10);
        assert_eq!(hits.len(), 1);
        // The hit stops exactly at the template end (nothing invented past it) and covers query bases 0..40.
        assert_eq!((hits[0].start, hits[0].end, hits[0].forward), (960, 1000, true));
        assert_eq!((hits[0].qstart, hits[0].qend), (0, 40));
        // Reverse complement: same template place, strand false, the aligned query part now at the query's end.
        let hits = align(&t, &revcomp(&q), false, 10);
        assert_eq!((hits[0].start, hits[0].end, hits[0].forward), (960, 1000, false));
        assert_eq!((hits[0].qstart, hits[0].qend), (60, 100));
        // Mirror: the query ENDS with the template's first 40 bases and hangs off the start.
        let mut q = overhang.clone();
        q.extend_from_slice(&t[..40]);
        let hits = align(&t, &q, false, 10);
        assert_eq!((hits[0].start, hits[0].end, hits[0].qstart, hits[0].qend), (0, 40, 60, 100));
    }

    #[test]
    fn overhanging_query_wraps_on_a_circle_and_a_tiny_overlap_is_no_hit() {
        let t = random_dna(1000, 11);
        // End of the template followed by its start: one hit across the origin on a circle, two halves on a line.
        let mut wrap = t[960..].to_vec();
        wrap.extend_from_slice(&t[..60]);
        let hits = align(&t, &wrap, true, 10);
        assert_eq!(hits.len(), 1);
        assert_eq!((hits[0].start, hits[0].end, hits[0].qstart, hits[0].qend), (960, 1060, 0, 100));
        assert_eq!(align(&t, &wrap, false, 10).len(), 2);
        // A 60 base query sharing only 15 bases with the template end: below the score threshold
        // (min(1.2 * 75, 40) = 40 points = 20 matching bases), so nothing is reported.
        let mut q = t[985..].to_vec();
        q.extend_from_slice(&random_dna(60, 12));
        assert!(align(&t, &q, false, 10).is_empty());
        assert!(align(&t, &q, true, 10).is_empty());
    }

    // ---- a query made of several pieces (a structural variant) -------------------------------------------------
    // `(a, b)` below is a tuple (a small fixed group of values); `Vec<(usize, usize, bool, usize, usize)>` is a list of
    // 5-value tuples; `.iter().map(..).collect()` turns each hit into such a tuple, `sort()` orders the list.

    /// The hits as (query start, query end, forward?, template start, template end), in query order.
    fn pieces(hits: &[Hit]) -> Vec<(usize, usize, bool, usize, usize)> {
        let mut v: Vec<_> = hits.iter().map(|h| (h.qstart, h.qend, h.forward, h.start, h.end)).collect();
        v.sort();
        v
    }

    /// Same list, but every number may be off by 2 (a chance matching base next to a piece is fair).
    fn assert_pieces(hits: &[Hit], want: &[(usize, usize, bool, usize, usize)]) {
        let got = pieces(hits);
        assert_eq!(got.len(), want.len(), "number of pieces: {got:?}");
        let near = |x: usize, y: usize| x.abs_diff(y) <= 2; // absolute difference of two usize numbers
        for (g, w) in got.iter().zip(want) {
            assert!(g.2 == w.2 && near(g.0, w.0) && near(g.1, w.1) && near(g.3, w.3) && near(g.4, w.4), "{g:?} vs {w:?}");
        }
    }

    /// 100 bases of `t` from `at`, optionally reverse-complemented.
    fn piece(t: &[u8], at: usize, inverted: bool) -> Vec<u8> {
        let p = t[at..at + 100].to_vec();
        if inverted { revcomp(&p) } else { p }
    }

    /// A (forward) + 200 unmatched bases + B (forward) + C (INVERTED) + D (forward), as a genome inversion plus an
    /// insertion would give. Returns the query and the pieces we expect to be reported.
    fn sv_query(t: &[u8], a: usize, b: usize, c: usize, d: usize, seed: u64) -> (Vec<u8>, Vec<(usize, usize, bool, usize, usize)>) {
        let mut q = piece(t, a, false);
        q.extend(random_dna(200, seed));
        for (at, inverted) in [(b, false), (c, true), (d, false)] {
            q.extend(piece(t, at, inverted));
        }
        let want = vec![(0, 100, true, a, a + 100), (300, 400, true, b, b + 100), (400, 500, false, c, c + 100), (500, 600, true, d, d + 100)];
        (q, want)
    }

    #[test]
    fn every_piece_of_a_patchwork_query_is_reported_including_the_inverted_one() {
        let t = random_dna(2000, 21);
        // Pieces far apart, in the same order as in the query.
        let (q, want) = sv_query(&t, 100, 600, 1100, 1600, 5);
        assert_pieces(&align(&t, &q, false, 10), &want);
        assert_pieces(&align(&t, &q, true, 10), &want); // and on a circular template
        // Pieces far apart in another order (the query reorders the template).
        let (q, want) = sv_query(&t, 1400, 900, 500, 100, 6);
        assert_pieces(&align(&t, &q, false, 10), &want);
        // Pieces close together.
        let (q, want) = sv_query(&t, 500, 650, 800, 950, 7);
        assert_pieces(&align(&t, &q, false, 10), &want);
    }

    #[test]
    fn pieces_on_one_diagonal_with_a_long_unmatched_middle_are_not_merged() {
        // A, B and D sit on the SAME diagonal (the template has the same spacing as the query), with 200 and
        // 100 unmatched bases between them. This used to be returned as one long, low-identity hit (or only
        // its best piece); each piece must now be reported on its own and none may span a gap.
        let t = random_dna(2000, 22);
        let (q, want) = sv_query(&t, 500, 800, 1500, 1000, 8);
        let hits = align(&t, &q, false, 10);
        assert_pieces(&hits, &want);
        assert!(hits.iter().all(|h| h.end - h.start <= 102 && h.identity > 0.97));
    }

    #[test]
    fn max_hits_still_limits_the_number_of_pieces_and_keeps_the_best() {
        let t = random_dna(2000, 23);
        let (q, _) = sv_query(&t, 100, 600, 1100, 1600, 9);
        assert_eq!(align(&t, &q, false, 2).len(), 2);
        assert_eq!(align(&t, &q, false, 1).len(), 1);
        assert!(align(&t, &q, false, 0).is_empty());
        let all = align(&t, &q, false, 10);
        assert!(all.windows(2).all(|w| w[0].score >= w[1].score), "best first"); // windows(2): each neighbouring pair
        // --- the cut keeps the BEST pieces, not the first found (the pieces above all score the same): pieces
        // of 80, 100 and 150 bases, the longest one inverted so it is found last (second strand pass) ---
        let mut q = t[100..180].to_vec();
        q.extend(random_dna(60, 31));
        q.extend_from_slice(&t[600..700]);
        q.extend(random_dna(60, 32));
        q.extend(revcomp(&t[1100..1250]));
        let all = align(&t, &q, false, 10);
        let one = align(&t, &q, false, 1);
        assert_eq!((one[0].start, one[0].end, one[0].forward), (1100, 1250, false));
        assert_eq!(align(&t, &q, false, 2), all[..2].to_vec());
    }

    #[test]
    fn a_piece_repeated_in_the_query_is_found_twice_on_the_same_template_place() {
        // A tandem duplication: the same 100 bases twice in the query, once in the template.
        let t = random_dna(2000, 24);
        let mut q = random_dna(60, 10);
        q.extend_from_slice(&t[700..800]);
        q.extend(random_dna(60, 11));
        q.extend_from_slice(&t[700..800]);
        q.extend(random_dna(60, 12));
        assert_pieces(&align(&t, &q, false, 10), &[(60, 160, true, 700, 800), (220, 320, true, 700, 800)]);
    }

    #[test]
    fn an_unmatched_stretch_is_cut_out_but_a_short_one_is_still_bridged() {
        let t = random_dna(2000, 25);
        // 8 unmatched bases between two 60 base pieces on one diagonal (a small cluster of changes): cheap
        // enough, so ONE bridged hit over the whole query, with the changes shown as mismatches.
        let mut q = t[300..360].to_vec();
        q.extend(random_dna(8, 13));
        q.extend_from_slice(&t[368..428]);
        let hits = align(&t, &q, false, 10);
        assert_eq!(hits.len(), 1);
        assert_eq!((hits[0].start, hits[0].end), (300, 428));
        // --- the 8 changed bases are mismatches at their own template places, not gaps ---
        let model: Vec<(usize, u8)> = (0..8).filter(|&i| q[60 + i] != t[360 + i]).map(|i| (360 + i, q[60 + i])).collect();
        assert_eq!((hits[0].qstart, hits[0].qend), (0, 128));
        assert_eq!(hits[0].mismatches, model);
        assert!(hits[0].insertions.is_empty() && hits[0].deletions.is_empty());
        // 40 unmatched bases cost more than the rule allows (MAX_DRAWDOWN): two pieces.
        let mut q = t[300..360].to_vec();
        q.extend(random_dna(40, 15));
        q.extend_from_slice(&t[400..460]);
        assert_pieces(&align(&t, &q, false, 10), &[(0, 60, true, 300, 360), (100, 160, true, 400, 460)]);
        // 120 unmatched bases between them: two hits, the stretch is not bridged.
        let mut q = t[300..360].to_vec();
        q.extend(random_dna(120, 14));
        q.extend_from_slice(&t[480..540]);
        assert_pieces(&align(&t, &q, false, 10), &[(0, 60, true, 300, 360), (180, 240, true, 480, 540)]);
    }
}
