//! The search bar's engine: a DNA motif (with IUPAC letters and a number of allowed mismatches, on either strand) or a
//! protein motif (in the three frames of both strands) in a sequence.
//!
//! Vocabulary, for readers who do not work with DNA:
//! - *IUPAC codes* are letters for "one of several bases": R is A or G, Y is C or T, N is any base, and so on. A
//!   query written with them matches every sequence they stand for.
//! - a *mismatch* is a position where the query and the sequence differ; "fuzzy" search allows a few.
//! - the *reverse strand*: DNA has two strands read in opposite directions, so a motif can also occur as its
//!   reverse complement (read backwards, with A<->T and C<->G swapped). Positions are always given on the top strand.
//! - a *frame*: a protein is read three bases at a time, and the start can be the 1st, 2nd or 3rd base: three
//!   frames on each strand, six in all.

use crate::code::{self, GeneticCode};
use crate::primers::revcomp;

/// One place a DNA motif occurs.
#[derive(Debug, Clone, PartialEq)]
pub struct MotifHit {
    /// First base on the top strand (0-based).
    pub start: usize,
    /// One past the last base; on a circular sequence this can pass the length (the hit wraps round the origin).
    pub end: usize,
    /// True if the motif reads along the top strand, false if its reverse complement does.
    pub forward: bool,
    /// How many positions differ.
    pub mismatches: usize,
}

/// Which bases a pattern position accepts, as bits: A = 1, C = 2, G = 4, T = 8.
/// `NOT_IUPAC_N` marks a literal N (used when IUPAC codes are off): it fits only an N of the sequence.
const NOT_IUPAC_N: u8 = 0x80;

/// The bits of the pattern letter `c`, or 0 when it is not a letter we can search for.
pub(crate) fn pattern_mask(c: u8, iupac: bool) -> u8 {
    match c.to_ascii_uppercase() {
        b'A' => 1,
        b'C' => 2,
        b'G' => 4,
        b'T' | b'U' => 8,
        b'N' if !iupac => NOT_IUPAC_N,
        b'R' if iupac => 1 | 4,
        b'Y' if iupac => 2 | 8,
        b'S' if iupac => 2 | 4,
        b'W' if iupac => 1 | 8,
        b'K' if iupac => 4 | 8,
        b'M' if iupac => 1 | 2,
        b'B' if iupac => 2 | 4 | 8,
        b'D' if iupac => 1 | 4 | 8,
        b'H' if iupac => 1 | 2 | 8,
        b'V' if iupac => 1 | 2 | 4,
        b'N' => 15,
        _ => 0,
    }
}

/// The bits of a base of the sequence: A, C, G, T one bit, N all four, anything else 0.
const fn base_table() -> [u8; 256] {
    let mut t = [0u8; 256];
    t[b'A' as usize] = 1;
    t[b'a' as usize] = 1;
    t[b'C' as usize] = 2;
    t[b'c' as usize] = 2;
    t[b'G' as usize] = 4;
    t[b'g' as usize] = 4;
    t[b'T' as usize] = 8;
    t[b't' as usize] = 8;
    t[b'N' as usize] = 15;
    t[b'n' as usize] = 15;
    t
}
static BASE: [u8; 256] = base_table();

/// Does the sequence base `b` fit the pattern position `mask`? A degenerate base of the sequence (N) fits only a
/// pattern position that accepts every base it could be.
#[inline]
fn fits(mask: u8, b: u8) -> bool {
    if mask == NOT_IUPAC_N {
        return b == b'N' || b == b'n';
    }
    let m = BASE[b as usize];
    m != 0 && m & !mask == 0
}

/// The pattern as it reads on the other strand: reversed, each position's bases complemented (A<->T, C<->G).
fn reverse_masks(masks: &[u8]) -> Vec<u8> {
    masks
        .iter()
        .rev()
        .map(|&m| {
            if m == NOT_IUPAC_N {
                m
            } else {
                // bit order A C G T = 1 2 4 8: complement swaps bit 0 with 3 and bit 1 with 2
                ((m & 1) << 3) | ((m & 2) << 1) | ((m & 4) >> 1) | ((m & 8) >> 3)
            }
        })
        .collect()
}

/// Every place `motif` (letters A C G T, and with `iupac` the ambiguity codes; U counts as T) occurs in `seq`, with at most
/// `max_mismatches` positions differing, on the top strand and, with `both`, as its reverse complement. A motif that is its
/// own reverse complement is reported once. On a `circular` sequence a hit can run across the origin. At most `limit` hits,
/// in order of position.
///
/// How it works: slide the motif along, comparing position by position and giving up on a place as soon as it has too
/// many mismatches (most places are given up on at the first or second base).
pub fn find_motif(
    seq: &[u8],
    motif: &[u8],
    max_mismatches: usize,
    iupac: bool,
    both: bool,
    circular: bool,
    limit: usize,
) -> Vec<MotifHit> {
    let m = motif.len();
    let n = seq.len();
    let masks: Vec<u8> = motif.iter().map(|&c| pattern_mask(c, iupac)).collect();
    if m == 0 || n == 0 || masks.contains(&0) || (!circular && m > n) {
        return Vec::new();
    }
    let reverse = reverse_masks(&masks);
    let mut patterns = vec![(masks.as_slice(), true)];
    if both && reverse != masks {
        patterns.push((reverse.as_slice(), false));
    }
    // A circular sequence is read with its start copied onto its end, so a hit can run across the origin.
    let wrapped: Vec<u8>;
    let text: &[u8] = if circular && m > 1 {
        wrapped = seq.iter().chain(seq[..(m - 1).min(n)].iter()).copied().collect();
        &wrapped
    } else {
        seq
    };
    let starts = if circular { n } else { n - m + 1 };
    let mut hits = Vec::new();
    'scan: for i in 0..starts {
        if i + m > text.len() {
            break;
        }
        let window = &text[i..i + m];
        for &(pattern, forward) in &patterns {
            let mut mismatches = 0;
            for (&mask, &b) in pattern.iter().zip(window) {
                if !fits(mask, b) {
                    mismatches += 1;
                    if mismatches > max_mismatches {
                        break;
                    }
                }
            }
            if mismatches <= max_mismatches {
                hits.push(MotifHit { start: i, end: i + m, forward, mismatches });
                if hits.len() >= limit {
                    break 'scan;
                }
            }
        }
    }
    hits
}

/// Every exact place of each of `motifs` (IUPAC letters) in `seq`, on both strands: for each motif, its hits in order of
/// position (`find_motif` with no mismatches and no limit, for many motifs in one pass: the enzymes' sites at every edit).
///
/// How it works: each pattern is filed under every 3-base word its first three positions accept, so at each place of the
/// sequence only the few patterns filed under the 3 bases found there are compared (all of them where those bases are not
/// plain A C G T).
pub fn find_motifs(seq: &[u8], motifs: &[&[u8]], circular: bool) -> Vec<Vec<MotifHit>> {
    const K: usize = 3;
    let n = seq.len();
    let mut out = vec![Vec::new(); motifs.len()];
    let mut patterns = Vec::new(); // (motif index, masks, forward)
    for (k, motif) in motifs.iter().enumerate() {
        let masks: Vec<u8> = motif.iter().map(|&c| pattern_mask(c, true)).collect();
        if masks.is_empty() || masks.contains(&0) || (!circular && masks.len() > n) {
            continue;
        }
        let reverse = reverse_masks(&masks);
        let palindrome = reverse == masks;
        patterns.push((k, masks, true));
        if !palindrome {
            patterns.push((k, reverse, false));
        }
    }
    let longest = patterns.iter().map(|p| p.1.len()).max().unwrap_or(0);
    if n == 0 || longest == 0 {
        return out;
    }
    let mut buckets: Vec<Vec<usize>> = vec![Vec::new(); 1 << (2 * K)];
    for (p, (_, masks, _)) in patterns.iter().enumerate() {
        for (word, bucket) in buckets.iter_mut().enumerate() {
            // position j of the word holds base (word >> 2j) & 3, as the bit 1 << that of BASE
            if (0..K.min(masks.len())).all(|j| masks[j] & (1 << ((word >> (2 * j)) & 3)) != 0) {
                bucket.push(p);
            }
        }
    }
    let all: Vec<usize> = (0..patterns.len()).collect();
    let wrapped: Vec<u8>;
    let text: &[u8] = if circular {
        wrapped = seq.iter().chain(&seq[..(longest - 1).min(n)]).copied().collect(); // (as find_motif reads a circle)
        &wrapped
    } else {
        seq
    };
    for i in 0..n {
        let mut word = 0;
        let mut plain = i + K <= text.len();
        for j in 0..K.min(text.len() - i) {
            match BASE[text[i + j] as usize] {
                1 => {}
                2 => word |= 1 << (2 * j),
                4 => word |= 2 << (2 * j),
                8 => word |= 3 << (2 * j),
                _ => plain = false,
            }
        }
        let checked = if plain { K } else { 0 }; // the bases the bucket already matched
        for &p in if plain { &buckets[word] } else { &all } {
            let (k, masks, forward) = &patterns[p];
            let m = masks.len();
            let from = checked.min(m);
            if i + m <= text.len() && masks[from..].iter().zip(&text[i + from..i + m]).all(|(&mask, &b)| fits(mask, b)) {
                out[*k].push(MotifHit { start: i, end: i + m, forward: *forward, mismatches: 0 });
            }
        }
    }
    out
}

/// Where each enzyme cuts `seq`: `enzymes` are (site in IUPAC letters, cuts as (top, bottom) offsets of the gap from the start of the
/// site, possibly outside it). For each, the sorted distinct (top, bottom) cut positions (the gap before that base on the top
/// strand's numbering): a site on the bottom strand cuts mirrored; on a `circular` sequence positions are taken round the
/// origin, on a linear one a cut on or past an end is no cut. (The app's `enzymes.find_cuts`.)
pub fn find_cuts(seq: &[u8], enzymes: &[(&[u8], Vec<(i64, i64)>)], circular: bool) -> Vec<Vec<(usize, usize)>> {
    let n = seq.len() as i64;
    let sites: Vec<&[u8]> = enzymes.iter().map(|e| e.0).collect();
    find_motifs(seq, &sites, circular)
        .into_iter()
        .zip(enzymes)
        .map(|(hits, (site, cuts))| {
            let l = site.len() as i64;
            let mut out = Vec::new();
            for h in hits {
                for &(top, bottom) in cuts {
                    // read the other way, the cut `t` from the start of the site is `l - t` from its end
                    let (t, b) = if h.forward { (top, bottom) } else { (l - bottom, l - top) };
                    let (t, b) = (h.start as i64 + t, h.start as i64 + b);
                    if circular {
                        out.push((t.rem_euclid(n) as usize, b.rem_euclid(n) as usize));
                    } else if 0 < t && t < n && 0 < b && b < n {
                        out.push((t as usize, b as usize));
                    }
                }
            }
            out.sort_unstable();
            out.dedup();
            out
        })
        .collect()
}

/// One place a protein motif occurs, read from DNA.
#[derive(Debug, Clone, PartialEq)]
pub struct ProteinHit {
    /// The DNA bases it is read from: top-strand coordinates, 0-based, end exclusive.
    pub start: usize,
    pub end: usize,
    /// True if read along the top strand.
    pub forward: bool,
    /// Which of the three frames of that strand (1, 2 or 3: the base the first codon starts at).
    pub frame: usize,
    pub mismatches: usize,
}

/// Every place the amino-acid `query` (one-letter codes, X = any, * = stop) occurs in the translations of `seq` in the
/// three frames of both strands, with at most `max_mismatches` amino acids differing. At most `limit` hits, in order of
/// position. Codons that contain an N or another letter translate to X, which only a query X matches.
pub fn find_protein(seq: &[u8], query: &[u8], max_mismatches: usize, limit: usize) -> Vec<ProteinHit> {
    find_protein_with(seq, query, max_mismatches, limit, code::standard())
}

/// `find_protein` reading the DNA with the genetic code `table` instead of the standard one.
pub fn find_protein_with(seq: &[u8], query: &[u8], max_mismatches: usize, limit: usize, table: &GeneticCode) -> Vec<ProteinHit> {
    let q: Vec<u8> = query.iter().map(|b| b.to_ascii_uppercase()).collect();
    let n = seq.len();
    if q.is_empty() || n < 3 * q.len() {
        return Vec::new();
    }
    let upper: Vec<u8> = seq.iter().map(|b| b.to_ascii_uppercase()).collect();
    let mut hits = Vec::new();
    let reverse = revcomp(&upper);
    for (forward, strand) in [(true, &upper), (false, &reverse)] {
        for frame in 0..3 {
            let protein = table.translate(&strand[frame..]);
            if protein.len() < q.len() {
                continue;
            }
            for i in 0..=protein.len() - q.len() {
                let mut mismatches = 0;
                for (&want, &got) in q.iter().zip(&protein[i..]) {
                    if want != b'X' && want != got {
                        mismatches += 1;
                        if mismatches > max_mismatches {
                            break;
                        }
                    }
                }
                if mismatches <= max_mismatches {
                    let (a, b) = (frame + 3 * i, frame + 3 * (i + q.len()));
                    // On the reverse strand the positions count from the other end.
                    let (start, end) = if forward { (a, b) } else { (n - b, n - a) };
                    hits.push(ProteinHit { start, end, forward, frame: frame + 1, mismatches });
                }
            }
        }
    }
    hits.sort_by_key(|h| (h.start, !h.forward));
    hits.truncate(limit);
    hits
}

#[cfg(test)]
mod tests {
    use super::*;

    fn starts(hits: &[MotifHit]) -> Vec<(usize, bool, usize)> {
        hits.iter().map(|h| (h.start, h.forward, h.mismatches)).collect()
    }

    #[test]
    fn finds_an_exact_motif_on_both_strands_in_position_order() {
        //           0         1         2
        //           012345678901234567890123
        let seq = b"TTGAATTCAAAAGGATCCAAGTTT";
        // GAATTC is its own reverse complement: once. GGATCC too.
        assert_eq!(starts(&find_motif(seq, b"GAATTC", 0, true, true, false, 100)), vec![(2, true, 0)]);
        // "AAGG" occurs at 10 on the top strand; CCTT, its reverse complement, is found at the same place on the bottom strand.
        assert_eq!(starts(&find_motif(seq, b"AAGG", 0, true, true, false, 100)), vec![(10, true, 0)]);
        assert_eq!(starts(&find_motif(seq, b"CCTT", 0, true, true, false, 100)), vec![(10, false, 0)], "found as the reverse complement");
        assert!(find_motif(seq, b"CCTT", 0, true, false, false, 100).is_empty(), "top strand only");
        // lower case in the sequence is fine
        assert_eq!(find_motif(b"ttgaattcaa", b"GAATTC", 0, true, true, false, 10).len(), 1);
    }

    #[test]
    fn iupac_letters_in_the_query_match_every_base_they_stand_for() {
        let seq = b"AAGAATTCAAGAGTTCAAGACTTC";
        // GARTTC: R = A or G. Hits: GAATTC at 2, GAGTTC at 10. And not GACTTC.
        let hits = find_motif(seq, b"GARTTC", 0, true, false, false, 100);
        assert_eq!(starts(&hits), vec![(2, true, 0), (10, true, 0)]);
        // With IUPAC off the R is not a base: no hit.
        assert!(find_motif(seq, b"GARTTC", 0, false, false, false, 100).is_empty());
        // N in the query matches anything; with IUPAC off it matches only an N of the sequence.
        assert_eq!(find_motif(b"ACGTACNT", b"GTAN", 0, true, false, false, 10).len(), 1);
        assert!(find_motif(b"ACGTACGT", b"GTAN", 0, false, false, false, 10).is_empty());
        assert_eq!(find_motif(b"ACGTACNT", b"CNT", 0, false, false, false, 10).len(), 1);
        // An N of the sequence matches only a query position that accepts every base.
        assert!(find_motif(b"ACGTNCGT", b"GTACGT", 0, true, false, false, 10).is_empty());
        assert_eq!(find_motif(b"ACGTNCGT", b"GTNCGT", 0, true, false, false, 10).len(), 1);
        // A letter that is not a base is no query at all.
        assert!(find_motif(seq, b"GAXTTC", 0, true, true, false, 100).is_empty());
    }

    #[test]
    fn allows_a_number_of_mismatches_and_reports_how_many() {
        let seq = b"TTTGAATTCTTTGAGTTCTTTGAGCTCTTT"; // GAATTC, then 1 off, then 2 off
        let hits = find_motif(seq, b"GAATTC", 1, true, true, false, 100);
        assert_eq!(starts(&hits), vec![(3, true, 0), (12, true, 1)]);
        let hits = find_motif(seq, b"GAATTC", 2, true, true, false, 100);
        assert_eq!(starts(&hits), vec![(3, true, 0), (12, true, 1), (21, true, 2)]);
        // and on the other strand: CCTA is the reverse complement of TAGG, one mismatch away from AAGG
        let hits = find_motif(b"GGCCTAGG", b"AAGG", 1, true, true, false, 100);
        assert_eq!(starts(&hits), vec![(2, false, 1), (4, true, 1)]);
    }

    #[test]
    fn a_hit_can_run_across_the_origin_of_a_circle() {
        // GAATTC split: GAA at the end, TTC at the start.
        let seq = b"TTCAAAAAAAAAAAAGAA";
        assert!(find_motif(seq, b"GAATTC", 0, true, true, false, 10).is_empty());
        let hits = find_motif(seq, b"GAATTC", 0, true, true, true, 10);
        assert_eq!(hits, vec![MotifHit { start: 15, end: 21, forward: true, mismatches: 0 }], "end passes the length: it wraps");
        // The reverse strand across the origin too: CCTAAG is the reverse complement of CTTAGG.
        let hits = find_motif(b"TAGGAAAAAAAACT", b"CCTAAG", 0, true, true, true, 10);
        assert_eq!(hits.len(), 1);
        assert!(!hits[0].forward && hits[0].start == 12);
    }

    #[test]
    fn many_motifs_in_one_pass_match_find_motif_one_by_one() {
        // A made-up sequence with lower case, N and a letter that is not a base; sites plain, palindromic, degenerate and short.
        let mut seq = Vec::new();
        let mut x: u32 = 7;
        for _ in 0..3000 {
            x = x.wrapping_mul(1_103_515_245).wrapping_add(12_345);
            seq.push(b"ACGTacgtNACGTACGTX"[(x >> 16) as usize % 18]);
        }
        let motifs: [&[u8]; 9] = [b"GAATTC", b"GGTCTC", b"GCNNNNNNNGC", b"RGATCY", b"AC", b"CCWGG", b"GAAGAC", b"NNN", b"GAXTTC"];
        for circular in [false, true] {
            let all = find_motifs(&seq, &motifs, circular);
            for (motif, hits) in motifs.iter().zip(&all) {
                assert_eq!(hits, &find_motif(&seq, motif, 0, true, true, circular, usize::MAX), "{}", String::from_utf8_lossy(motif));
            }
        }
        assert_eq!(find_motifs(b"TTCAAAGAA", &[b"GAATTC"], true)[0].len(), 1, "a site across the origin");
        // A site longer than a (tiny) circle: read once round, as find_motif does.
        let tiny = find_motifs(b"AATTCG", &[b"CGAATTC"], true);
        assert_eq!(tiny[0], find_motif(b"AATTCG", b"CGAATTC", 0, true, true, true, usize::MAX));
        assert!(!tiny[0].is_empty(), "a site longer than the circle");
        assert!(find_motifs(b"", &[b"GAATTC"], false)[0].is_empty());
    }

    #[test]
    fn cuts_mirror_on_the_bottom_strand_wrap_on_a_circle_and_stop_at_linear_ends() {
        let ecori: (&[u8], Vec<(i64, i64)>) = (b"GAATTC", vec![(1, 5)]);
        let bsai: (&[u8], Vec<(i64, i64)>) = (b"GGTCTC", vec![(7, 11)]); // cuts outside its site, not its own reverse complement
        // EcoRI at 2: cuts after G (3) and before the last C on the bottom strand (7).
        assert_eq!(find_cuts(b"AAGAATTCAA", &[ecori.clone()], false), vec![vec![(3, 7)]]);
        // BsaI forward at 0: (7, 11).
        assert_eq!(find_cuts(b"GGTCTCAAAAAAAAAA", &[bsai.clone()], false), vec![vec![(7, 11)]]);
        // BsaI on the bottom strand (GAGACC at 10): mirrored, (6 - 11, 6 - 7) from the site's start = (5, 9).
        assert_eq!(find_cuts(b"AAAAAAAAAAGAGACCAAAA", &[bsai.clone()], false), vec![vec![(5, 9)]]);
        // Linear: a cut past the end is no cut; circular: it is taken round the origin (11 % 8 = 3).
        assert_eq!(find_cuts(b"GGTCTCAA", &[bsai.clone()], false), vec![vec![]]);
        assert_eq!(find_cuts(b"GGTCTCAA", &[bsai, ecori], true), vec![vec![(7, 3)], vec![]]);
    }

    #[test]
    fn stops_at_the_limit_and_handles_tiny_and_empty_inputs() {
        let seq = b"ACACACACACACACAC";
        assert_eq!(find_motif(seq, b"AC", 0, true, false, false, 3).len(), 3);
        assert!(find_motif(b"", b"AC", 0, true, true, false, 10).is_empty());
        assert!(find_motif(b"AC", b"", 0, true, true, false, 10).is_empty());
        assert!(find_motif(b"AC", b"ACGT", 0, true, true, false, 10).is_empty());
        assert_eq!(find_motif(b"A", b"A", 0, true, true, true, 10).len(), 1);
    }

    #[test]
    fn finds_a_protein_motif_in_all_six_frames() {
        // ATG AAA GGG TTT = M K G F, starting at base 3 (0-based), so frame 1. On the reverse strand: the same four codons reverse-complemented.
        let fwd = b"CCCATGAAAGGGTTTCC";
        let hits = find_protein(fwd, b"MKGF", 0, 10);
        assert_eq!(hits, vec![ProteinHit { start: 3, end: 15, forward: true, frame: 1, mismatches: 0 }]);
        let rev = revcomp(fwd);
        let hits = find_protein(&rev, b"MKGF", 0, 10);
        assert_eq!(hits.len(), 1);
        assert!(!hits[0].forward);
        assert_eq!((hits[0].start, hits[0].end), (2, 14), "top-strand coordinates of the codons");
        // the other frames
        assert_eq!(find_protein(b"ATGAAAGGGTTT", b"KGF", 0, 10)[0].frame, 1);
        assert_eq!(find_protein(b"CATGAAAGGGTTT", b"KGF", 0, 10)[0].frame, 2);
        assert_eq!(find_protein(b"CCATGAAAGGGTTT", b"KGF", 0, 10)[0].frame, 3);
        // --- the three reverse frames: MKGF on the bottom strand at frame - 1, always top-strand 0..12 ---
        for (fwd, frame) in [(&b"ATGAAAGGGTTT"[..], 1usize), (b"CATGAAAGGGTTT", 2), (b"CCATGAAAGGGTTT", 3)] {
            let rc = revcomp(fwd);
            assert_eq!(find_protein(&rc, b"MKGF", 0, 10), vec![ProteinHit { start: 0, end: 12, forward: false, frame, mismatches: 0 }], "frame {frame}");
        }
        // X in the query is any amino acid; a mismatch is allowed when asked for
        assert_eq!(find_protein(fwd, b"MXGF", 0, 10).len(), 1);
        assert!(find_protein(fwd, b"MLGF", 0, 10).is_empty());
        assert_eq!(find_protein(fwd, b"MLGF", 1, 10)[0].mismatches, 1);
        // a stop codon is *
        assert_eq!(find_protein(b"ATGTAAGGG", b"M*G", 0, 10).len(), 1);
        assert!(find_protein(b"AT", b"M", 0, 10).is_empty());
    }
}
