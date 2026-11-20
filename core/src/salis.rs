//! Two pieces of the Synthesis Success Calculator that are too slow in Python.
//!
//! Halper SM, Hossain A, Salis HM. "Synthesis success calculator: predicting the rapid synthesis of
//! DNA fragments with machine learning." ACS Synth Biol 2020, 9(7):1563-1571, doi:10.1021/acssynbio.9b00460.
//! Code: https://github.com/hsalis/SalisLabCode (Synthesis_Success_Calculator), GPLv3, (c) the authors.
//!
//! This is a port of `RepeatFinder.py` (Ayaan Hossain) and of the rule-based part of `seq_eval.py`
//! (`_is_hairpin_pass` / `fast_hairpins`), keeping their rules and quirks. The one change is speed:
//! the repeat finder looks k-mers up in a hash index instead of re-scanning the sequence with KMP,
//! which returns the same positions. Sequences must be A/C/G/T only (the caller checks).

use std::collections::{BTreeSet, HashMap, HashSet};

fn comp(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'G' => b'C',
        b'C' => b'G',
        o => o,
    }
}

/// One maximal repeat: a k-mer of `len` bases found at `starts` (0-based, sorted, no duplicates).
/// `invert` is the kind the original filed it under ("direct" or "invert").
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Repeat {
    pub len: usize,
    pub invert: bool,
    pub starts: Vec<usize>,
}

/// All maximal direct and inverted repeats of length `k_low` or more, as `RepeatFinder.get_repeat_dict`
/// returns them for one sequence. Quirk kept on purpose: every other occurrence of a k-mer is filed
/// as "direct", even when the k-mer was found as an inverted repeat.
pub fn maximal_repeats(seq: &[u8], k_low: usize) -> Vec<Repeat> {
    let n = seq.len();
    let mut found: HashMap<(Vec<u8>, bool), BTreeSet<usize>> = HashMap::new();
    // The original loops k over range(k_low, n - 1) and stops at the first k with no repeated k-mer.
    for k in k_low..n.saturating_sub(1) {
        let mut index: HashMap<&[u8], Vec<usize>> = HashMap::new();
        for s in 0..=n - k {
            index.entry(&seq[s..s + k]).or_default().push(s);
        }
        let mut covered: HashSet<Vec<u8>> = HashSet::new();
        for k_start in 0..=n - k {
            let k_end = k_start + k;
            let direct: Vec<u8> = seq[k_start..k_end].to_vec();
            let inverted: Vec<u8> = direct.iter().rev().map(|&b| comp(b)).collect();
            for (kmer, invert) in [(direct, false), (inverted, true)] {
                if covered.contains(&kmer) {
                    continue;
                }
                // Every other occurrence at or after k_start + 1, overlapping ones included.
                let entry: Vec<usize> = match index.get(kmer.as_slice()) {
                    Some(v) => v.iter().copied().filter(|&p| p > k_start).collect(),
                    None => Vec::new(),
                };
                if entry.is_empty() {
                    continue;
                }
                covered.insert(kmer.clone());
                if !left_maximal(seq, k_start, k, &entry, invert) || !right_maximal(seq, k_start, k, &entry, invert) {
                    continue;
                }
                let at = found.entry((kmer.clone(), invert)).or_default();
                at.insert(k_start);
                let direct_set = found.entry((kmer, false)).or_default();
                direct_set.extend(entry);
            }
        }
        if covered.is_empty() {
            break;
        }
    }
    let mut out: Vec<Repeat> = found
        .into_iter()
        .map(|((kmer, invert), set)| Repeat { len: kmer.len(), invert, starts: set.into_iter().collect() })
        .collect();
    out.sort_by(|a, b| (a.len, a.invert, &a.starts).cmp(&(b.len, b.invert, &b.starts)));
    out
}

fn left_maximal(seq: &[u8], k_start: usize, k: usize, entry: &[usize], invert: bool) -> bool {
    let k_end = k_start + k;
    if (!invert && k_start == 0) || (invert && k_end == seq.len()) {
        return true;
    }
    entry.iter().any(|&j| {
        if invert {
            comp(seq[k_end]) != seq[j - 1]
        } else {
            seq[k_start - 1] != seq[j - 1]
        }
    })
}

fn right_maximal(seq: &[u8], k_start: usize, k: usize, entry: &[usize], invert: bool) -> bool {
    let k_end = k_start + k;
    if (!invert && k_end == seq.len()) || (invert && k_start == 0) {
        return true;
    }
    if entry[entry.len() - 1] + k == seq.len() {
        return true;
    }
    entry.iter().any(|&j| {
        if invert {
            comp(seq[k_start - 1]) != seq[j + k]
        } else {
            seq[k_end] != seq[j + k]
        }
    })
}

fn is_gc(b: u8) -> bool {
    b == b'G' || b == b'C'
}

/// One pass of `_is_hairpin_pass`: every stem/loop/stem window whose stem has at most `max_mismatch`
/// mismatches; `best` rises to the highest GC fraction of a stem that reaches `gc_high`.
fn hairpin_pass(seq: &[u8], stem: usize, loop_: usize, max_mismatch: usize, gc_high: f64, best: &mut f64) {
    let total = stem + loop_ + stem;
    if seq.len() < total {
        return;
    }
    for i in 0..=seq.len() - total {
        let (mut j, mut mismatch, mut gc) = (0usize, 0usize, 0.0f64);
        while j < stem {
            if is_gc(seq[i + j]) {
                gc += 1.0;
            }
            if is_gc(seq[i + total - j - 1]) {
                gc += 1.0;
            }
            if seq[i + j] != comp(seq[i + total - j - 1]) {
                mismatch += 1;
            }
            if mismatch > max_mismatch {
                break;
            }
            j += 1;
        }
        if j == stem {
            let content = (gc / 2.0) / stem as f64;
            if content >= gc_high && content > *best {
                *best = content;
            }
        }
    }
}

/// The highest stem GC fraction among the original's rule-based hairpins (11 bp stems with up to 2
/// mismatches and GC >= 0.8, 17 bp with up to 3, 22 bp with up to 5 for sequences over 500 bp, and
/// 11 bp palindromes with up to 1), or 0.0 when there is none.
pub fn hairpin_rule_best_gc(seq: &[u8]) -> f64 {
    let mut best = 0.0;
    for loop_ in 3..=48 {
        hairpin_pass(seq, 11, loop_, 2, 0.80, &mut best);
    }
    for loop_ in 3..=100 {
        hairpin_pass(seq, 17, loop_, 3, 0.0, &mut best);
    }
    if seq.len() > 500 {
        for loop_ in 100..=500 {
            hairpin_pass(seq, 22, loop_, 5, 0.0, &mut best);
        }
    }
    hairpin_pass(seq, 11, 0, 1, 0.0, &mut best);
    best
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_direct_repeat_is_found_once_at_its_full_length_with_both_starts() {
        // Hand-built: the same 12-mer twice with a 20-base spacer; nothing else is 9 bases long.
        let x = "GATTACAGGCTA";
        let y = "CCTTGGAACCTTGGTTCCAA";
        let seq = format!("{x}{y}{x}");
        let r = maximal_repeats(seq.as_bytes(), 9);
        let long: Vec<&Repeat> = r.iter().filter(|r| r.len == 12 && !r.invert).collect();
        assert_eq!(long.len(), 1);
        assert_eq!(long[0].starts, vec![0, 32]);
        assert!(r.iter().all(|r| r.len >= 9));
    }

    #[test]
    fn a_sequence_without_repeats_has_none() {
        assert!(maximal_repeats(b"ACGTTGCAAGCTTAGGCCATAT", 9).is_empty());
        assert!(maximal_repeats(b"ACG", 9).is_empty());
    }

    #[test]
    fn a_gc_rich_stem_scores_one_and_a_poly_a_run_scores_zero() {
        let seq = "GCGCGCGCGCGAAACGCGCGCGCGC";
        assert_eq!(hairpin_rule_best_gc(seq.as_bytes()), 1.0);
        assert_eq!(hairpin_rule_best_gc("A".repeat(60).as_bytes()), 0.0);
    }
}
