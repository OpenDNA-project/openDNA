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

use crate::primers::revcomp;

// Scores: +2 match, -3 mismatch, a gap of n bases costs 5 + 2 * (n - 1).
const MATCH: i32 = 2;
const MISMATCH: i32 = -3;
const GAP_OPEN: i32 = -5;
const GAP_EXTEND: i32 = -2;
/// Diagonals searched on each side of a seed cluster.
const BAND: i64 = 16;
/// Clusters aligned per strand (the ones with most seeds).
const MAX_CLUSTERS: usize = 200;
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

/// Uppercase, U -> T.
fn clean(seq: &[u8]) -> Vec<u8> {
    seq.iter().map(|b| match b.to_ascii_uppercase() { b'U' => b'T', c => c }).collect()
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

/// Seeds as (diagonal, template position) for one strand of the query.
fn seeds(text: &[u8], query: &[u8], k: usize) -> Vec<(i64, usize)> {
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
    // Very repeated k-mers would flood the seed list: ignore them.
    index.retain(|_, v| v.len() <= 100);
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
                break;
            }
        }
    }
    out
}

/// Groups seeds on nearby diagonals: (lowest diagonal, highest diagonal, number of seeds).
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
    out.truncate(MAX_CLUSTERS);
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
    // `extra` counts these added exact columns (they all match).
    let mut extra = 0usize;
    while qend < m && end < text.len() && query[qend] == text[end] && code(query[qend]).is_some() {
        end += 1; qend += 1; extra += 1;
    }
    while qstart > 0 && start > 0 && query[qstart - 1] == text[start - 1] && code(query[qstart - 1]).is_some() {
        start -= 1; qstart -= 1; extra += 1;
    }
    let (mut mismatches, mut insertions, mut deletions) = (Vec::new(), Vec::new(), Vec::new());
    let (mut matches, mut columns) = (0usize, 0usize);
    // Running score along the alignment, its best value so far, the query position reached, and
    // what the previous column was (0 = match/mismatch, 1 = deletion, 2 = insertion: a gap column
    // right after a gap column of the same kind is a cheap "extend", otherwise a dear "open").
    let (mut running, mut peak, mut qpos, mut prev) = (0i32, 0i32, i, 0u8);
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
        score: running + MATCH * extra as i32, // the score of the (possibly trimmed) alignment
        identity: (matches + extra) as f32 / (columns + extra).max(1) as f32,
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
fn pieces_of_cluster(text: &[u8], query: &[u8], lo: i64, hi: i64, min_score: i32) -> Vec<Hit> {
    let mut found: Vec<Hit> = Vec::new();
    // A to-do list of query ranges (from, to) still to be searched. `vec![(0, query.len())]`
    // builds a Vec holding that one tuple; `pop()` takes the last item (or None when empty).
    let mut todo = vec![(0usize, query.len())];
    while let Some((from, to)) = todo.pop() {
        if to - from < 4 || found.len() >= MAX_PIECES_PER_CLUSTER {
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
    found
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
pub fn align(template: &[u8], query: &[u8], circular: bool, max_hits: usize) -> Vec<Hit> {
    let (template, query) = (clean(template), clean(query));
    let (n, m) = (template.len(), query.len());
    if n == 0 || m < 4 || max_hits == 0 {
        return Vec::new();
    }
    // Circular: append the start of the template so hits can run through the origin.
    let mut text = template.clone();
    if circular {
        text.extend_from_slice(&template[..n.min(m + 2 * BAND as usize)]);
    }
    // Short queries get short k-mers: a 60-mer with 3 mismatches often has no two shared 11-mers.
    let k = if m >= 200 { 11 } else { (m / 4).clamp(4, 9) };
    let min_score = (m as i32 * 6 / 5).min(40); // random matches score ~20-30 on long queries
    let mut hits: Vec<Hit> = Vec::new();
    for forward in [true, false] {
        let q = if forward { query.clone() } else { revcomp(&query) };
        for (lo, hi, _) in clusters(seeds(&text, &q, k), 2) {
            for mut hit in pieces_of_cluster(&text, &q, lo - BAND, hi + BAND, min_score) {
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
    }
    hits.sort_by(|a, b| b.score.cmp(&a.score));
    let mut kept: Vec<Hit> = Vec::new();
    for h in hits {
        if kept.len() < max_hits && !kept.iter().any(|k| redundant(k, &h, if circular { n } else { 0 })) {
            kept.push(h);
        }
    }
    kept
}

#[cfg(test)]
mod tests {
    use super::*;

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
        // Edit from the back so earlier positions stay valid: delete 5 bases, insert 3, mutate 1.
        q.drain(250..255);
        for (i, b) in b"GGG".iter().enumerate() { q.insert(150 + i, *b); }
        let wrong = if q[50] == b'A' { b'C' } else { b'A' };
        q[50] = wrong;
        let hits = align(&t, &q, false, 10);
        let h = &hits[0];
        assert!(h.forward && h.start == 1000 && h.end == 1300, "{h:?}");
        assert_eq!(h.mismatches, vec![(1050, wrong)]);
        assert_eq!(h.insertions.iter().map(|x| x.1).sum::<usize>(), 3);
        assert_eq!(h.deletions.iter().map(|x| x.1 - x.0).sum::<usize>(), 5);
        assert!(h.identity > 0.95);
        // Same on the other strand.
        let hits = align(&t, &revcomp(&q), false, 10);
        assert!(!hits[0].forward && hits[0].start == 1000 && hits[0].end == 1300);
    }

    #[test]
    fn a_mismatch_near_either_end_of_a_long_query_is_reported_not_trimmed_away() {
        let t = random_dna(3000, 12);
        for at in [6, 11, 288, 293] {
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
