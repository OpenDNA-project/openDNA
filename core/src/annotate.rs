//! Find known parts (promoters, origins, resistance genes, fluorescent proteins…) in a
//! sequence, for "Detect common features".
//!
//! DNA level: every 12-letter word of the template is indexed once. For each part (and
//! its reverse complement) a few of its 12-letter words are looked up; each hit proposes
//! an alignment ("diagonal": part placed at template position p - i), which is then scored
//! over the whole part. Kept if identity >= min_dna.
//!
//! Protein level (coding parts only): the same on translated sequences, in the 3 frames
//! of both strands, with 5-letter words. Finds codon-optimized versions of a protein,
//! whose DNA differs a lot.
//!
//! shortcut: alignments are ungapped (substitutions only), so a part with an insertion
//! or deletion inside is missed or split; add a banded gapped extension if that matters.

use std::collections::{HashMap, HashSet};

use crate::code::{self, GeneticCode};
use crate::primers::revcomp;

const K_DNA: usize = 12;
const K_AA: usize = 5;
/// Look up one word every STEP positions of a part (any exact 12-mer finds the alignment;
/// at 90 % identity there are many, so sampling loses nothing in practice and is ~4x faster).
const STEP: usize = 4;

#[derive(Clone, Debug, PartialEq)]
pub struct Hit {
    pub part: usize,     // index into the parts list
    pub start: usize,    // 0-based on the top strand; start + len may pass the origin
    pub len: usize,
    pub forward: bool,
    pub identity: f32,   // 0..1, of bases (DNA hit) or amino acids (protein hit)
    pub protein: bool,
}

/// Amino acids of `seq` read from its first base ('X' for codons with N or other letters), standard code.
/// The genetic code itself lives in code.rs (`GeneticCode`); `translate_with` takes any table.
pub fn translate(seq: &[u8]) -> Vec<u8> {
    translate_with(seq, code::standard())
}

pub fn translate_with(seq: &[u8], table: &GeneticCode) -> Vec<u8> {
    table.translate(seq)
}

/// word -> positions where it occurs. `&[u8]` slices of `seq` are the keys (no copies).
fn index(seq: &[u8], k: usize) -> HashMap<&[u8], Vec<usize>> {
    let mut map: HashMap<&[u8], Vec<usize>> = HashMap::new();
    if seq.len() >= k {
        for i in 0..=seq.len() - k {
            map.entry(&seq[i..i + k]).or_default().push(i);
        }
    }
    map
}

/// Alignment offsets of `query` in `target` proposed by shared words (see module doc).
fn diagonals(query: &[u8], index: &HashMap<&[u8], Vec<usize>>, k: usize) -> HashSet<isize> {
    let mut out = HashSet::new();
    if query.len() < k {
        return out;
    }
    let last = query.len() - k;
    // Every STEP-th word, plus the last one so the part's end is always sampled.
    for i in (0..=last).step_by(STEP).chain(std::iter::once(last)) {
        if let Some(positions) = index.get(&query[i..i + k]) {
            out.extend(positions.iter().map(|&p| p as isize - i as isize));
        }
    }
    out
}

/// Identity of `query` placed at `target[d..]`, or None if below `min` (stops early).
fn score(query: &[u8], target: &[u8], d: isize, min: f32) -> Option<f32> {
    if d < 0 || d as usize + query.len() > target.len() || query.is_empty() {
        return None;
    }
    let allowed = ((1.0 - min) * query.len() as f32).floor() as usize;
    let mut mismatches = 0;
    for (a, b) in query.iter().zip(&target[d as usize..]) {
        if a != b {
            mismatches += 1;
            if mismatches > allowed {
                return None;
            }
        }
    }
    Some(1.0 - mismatches as f32 / query.len() as f32)
}

/// Where each part occurs in `template` (both strands; across the origin if `circular`).
/// `coding[i]`: part i is a CDS, also searched as protein. Identities are 0..1.
pub fn find_parts(template: &[u8], circular: bool, parts: &[Vec<u8>], coding: &[bool], min_dna: f32, min_protein: f32) -> Vec<Hit> {
    find_parts_with(template, circular, parts, coding, min_dna, min_protein, code::standard())
}

/// `find_parts` with a chosen genetic code, used to translate both the template and the coding parts (the protein
/// level then finds a part whose protein is the same under that code).
pub fn find_parts_with(
    template: &[u8],
    circular: bool,
    parts: &[Vec<u8>],
    coding: &[bool],
    min_dna: f32,
    min_protein: f32,
    table: &GeneticCode,
) -> Vec<Hit> {
    let n = template.len();
    let longest = parts.iter().map(|p| p.len()).max().unwrap_or(0);
    // Circular: append the start again, so a part across the origin is one stretch.
    let mut top: Vec<u8> = template.iter().map(|b| b.to_ascii_uppercase()).collect();
    if circular && n > 0 {
        let extra = longest.saturating_sub(1).min(n);
        top.extend_from_slice(&template[..extra].to_ascii_uppercase());
    }
    let bottom = revcomp(&top);
    let len2 = top.len();
    let dna_index = index(&top, K_DNA);
    // Translated frames: (forward?, frame, amino acids, word index), built only if needed.
    let frames: Vec<(bool, usize, Vec<u8>)> = if coding.iter().any(|&c| c) {
        [(true, &top), (false, &bottom)]
            .iter()
            .flat_map(|(fwd, s)| (0..3).map(move |f| (*fwd, f, translate_with(&s[f.min(s.len())..], table))))
            .collect()
    } else {
        Vec::new()
    };
    let aa_indexes: Vec<_> = frames.iter().map(|(_, _, aa)| index(aa, K_AA)).collect();

    let mut hits = Vec::new();
    for (pi, part) in parts.iter().enumerate() {
        let part: Vec<u8> = part.iter().map(|b| b.to_ascii_uppercase()).collect();
        let mut found: Vec<Hit> = Vec::new();
        // DNA: the part as is (forward hit), and its reverse complement (reverse hit).
        for (forward, q) in [(true, part.clone()), (false, revcomp(&part))] {
            let candidates: Vec<isize> = if q.len() < K_DNA {
                // Too short for the index: plain exact search.
                (0..=len2.saturating_sub(q.len())).filter(|&p| top[p..].starts_with(&q)).map(|p| p as isize).collect()
            } else {
                diagonals(&q, &dna_index, K_DNA).into_iter().collect()
            };
            for d in candidates {
                let min = if q.len() < K_DNA { 1.0 } else { min_dna };
                if let Some(identity) = score(&q, &top, d, min) {
                    if (d as usize) < n || n == 0 {
                        found.push(Hit { part: pi, start: d as usize, len: q.len(), forward, identity, protein: false });
                    }
                }
            }
        }
        // Protein (coding parts not already found as DNA).
        if coding.get(pi).copied().unwrap_or(false) && found.is_empty() {
            let mut prot = translate_with(&part, table);
            while prot.last() == Some(&b'*') {
                prot.pop(); // the stop codon is not part of the protein
            }
            for ((forward, frame, aa), idx) in frames.iter().zip(&aa_indexes) {
                for d in diagonals(&prot, idx, K_AA) {
                    let Some(identity) = score(&prot, aa, d, min_protein) else { continue };
                    let len = part.len().min(len2); // keep the part's own length (stop included)
                    let s = frame + 3 * d as usize; // start on this strand
                    // Bottom-strand coordinates back to the top strand.
                    let start = if *forward { s } else { len2.saturating_sub(s + len) };
                    if start < n {
                        found.push(Hit { part: pi, start, len, forward: *forward, identity, protein: true });
                    }
                }
            }
        }
        // Overlapping hits of the same part and strand: keep the best one.
        found.sort_by(|a, b| b.identity.total_cmp(&a.identity));
        let mut kept: Vec<Hit> = Vec::new();
        for h in found {
            let overlaps = kept.iter().any(|k| k.forward == h.forward && h.start < k.start + k.len && k.start < h.start + h.len);
            if !overlaps {
                kept.push(h);
            }
        }
        hits.extend(kept);
    }
    hits.sort_by_key(|h| (h.start, h.part));
    hits
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Deterministic pseudo-random DNA (no external crate needed).
    fn random_dna(n: usize, mut seed: u64) -> Vec<u8> {
        (0..n)
            .map(|_| {
                seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                b"ACGT"[(seed >> 33) as usize % 4]
            })
            .collect()
    }

    #[test]
    fn finds_parts_on_both_strands_across_the_origin_and_as_protein() {
        let backbone = random_dna(3000, 1);
        let gene = random_dna(300, 2); // pretend CDS (100 codons)
        let promoter = random_dna(40, 3);
        let mut plasmid = backbone[..500].to_vec();
        plasmid.extend(&promoter); // 500..540, forward
        plasmid.extend(&backbone[500..1000]);
        let mut mutated = revcomp(&gene); // reverse strand at 1040..1340, 2 % mismatches
        for i in (0..300).step_by(50) {
            mutated[i] = if mutated[i] == b'A' { b'C' } else { b'A' };
        }
        plasmid.extend(&mutated);
        plasmid.extend(&backbone[1000..2000]);
        // The promoter again, across the origin: its first 15 bases at the end, the rest at the start.
        plasmid.extend(&promoter[..15]);
        let mut circ = promoter[15..].to_vec();
        circ.extend(&plasmid);
        let n = circ.len();

        let parts = vec![promoter.clone(), gene.clone()];
        let hits = find_parts(&circ, true, &parts, &[false, true], 0.95, 0.9);
        let simple: Vec<_> = hits.iter().map(|h| (h.part, h.start, h.len, h.forward, h.protein)).collect();
        assert!(simple.contains(&(0, 25 + 500, 40, true, false)), "{simple:?}"); // the promoter, shifted by 25
        assert!(simple.contains(&(0, n - 15, 40, true, false)), "{simple:?}"); // across the origin
        let g = hits.iter().find(|h| h.part == 1).unwrap();
        assert_eq!((g.start, g.forward, g.protein), (25 + 1040, false, false));
        assert!((g.identity - 294.0 / 300.0).abs() < 1e-6);
        // Linear: no hit across the origin.
        assert_eq!(find_parts(&circ, false, &parts[..1], &[false], 0.95, 0.9).len(), 1);

        // A recoded gene (synonymous codons): the DNA differs a lot, the protein not at all.
        let cds = b"ATGCTGCTGCTGAGCAGCAGCCGTCGTCGTGGTGGTGGTCTGCTGCTGTAA".to_vec(); // M LLL SSS RRR GGG LLL *
        let recoded = b"ATGTTATTGCTTTCTTCCTCACGCAGAAGGGGAGGCGGGTTGCTACTCTAG".to_vec(); // same protein
        let mut host = random_dna(200, 4);
        host.extend(&recoded);
        host.extend(random_dna(200, 5));
        let hits = find_parts(&host, false, &[cds.clone()], &[true], 0.95, 0.9);
        assert_eq!(hits.len(), 1, "{hits:?}");
        assert_eq!((hits[0].start, hits[0].len, hits[0].forward, hits[0].protein), (200, cds.len(), true, true));
        assert_eq!(hits[0].identity, 1.0);
    }
}
