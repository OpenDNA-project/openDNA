//! Small sequence utilities: open reading frames.

use crate::code::{self, GeneticCode};
use crate::primers::revcomp;

/// One open reading frame, top-strand coordinates `[start, end)` including the stop codon.
#[derive(Debug, PartialEq)]
pub struct Orf {
    pub start: usize,
    pub end: usize,
    pub forward: bool,
}

/// ORFs of at least `min_codons` codons (stop included), ATG to stop, in all 6 frames.
/// For each stop codon the longest ORF is reported (the first ATG after the previous stop).
/// `circular`: the sequence is a ring, so an ORF may start near the end and finish after the origin;
/// its `end` is then larger than the sequence length (the same convention as a feature that wraps).
pub fn find_orfs(seq: &[u8], min_codons: usize, circular: bool) -> Vec<Orf> {
    find_orfs_with(seq, min_codons, circular, code::standard(), false)
}

/// `find_orfs` for any genetic code: the stop codons are the table's, and the start codons are ATG only
/// (`alternative_starts` false, the default of the app) or every start codon the table lists (GTG, TTG, CTG, ATT...).
/// Reading from the first start after the previous stop means the longest ORF is reported for each stop.
pub fn find_orfs_with(seq: &[u8], min_codons: usize, circular: bool, table: &GeneticCode, alternative_starts: bool) -> Vec<Orf> {
    let top: Vec<u8> = seq.iter().map(|b| b.to_ascii_uppercase()).collect();
    let bottom = revcomp(&top);
    let n = top.len();
    let mut out = Vec::new();
    for (strand, s) in [(true, &top), (false, &bottom)] {
        // On a ring, read the sequence twice in a row: a start in the first copy can then run on past
        // the origin into the second. Only starts in the first copy count (the second copy would find
        // every ORF again), and an ORF is never longer than one lap.
        let text: Vec<u8> = if circular { [s.as_slice(), s.as_slice()].concat() } else { s.clone() };
        for frame in 0..3 {
            let mut start: Option<usize> = None;
            let mut i = frame;
            while i + 3 <= text.len() {
                if start.is_none() && i >= n {
                    break; // nothing more can start here: the rest is the second copy
                }
                if matches!(start, Some(a) if i + 3 > a + n) {
                    start = None; // a full lap without a stop: not an ORF
                }
                let codon = &text[i..i + 3];
                if start.is_none() && table.is_start(codon, alternative_starts) {
                    start = Some(i);
                } else if table.is_stop(codon) {
                    if let Some(a) = start.take() {
                        let end = i + 3;
                        if (end - a) / 3 >= min_codons {
                            // Bottom-strand coordinates map back to the top strand (modulo the length,
                            // because on a ring the ORF may start before the origin).
                            let (s0, e0) = if strand {
                                (a, end)
                            } else {
                                let s0 = (n as i64 - end as i64).rem_euclid(n as i64) as usize;
                                (s0, s0 + (end - a))
                            };
                            out.push(Orf { start: s0, end: e0, forward: strand });
                        }
                    }
                }
                i += 3;
            }
        }
    }
    if circular && n > 0 {
        // Each frame is read from base 0, unaware of an ORF that started before the origin and is still open there: an ATG
        // near base 0 then starts a shorter ORF ending at the same stop. On a ring one stop has one ORF, the longest, so keep
        // only that one per strand and stop (the stop's place taken modulo the length; in top-strand coordinates a bottom-strand
        // ORF's stop is at its start).
        let stop = |o: &Orf| if o.forward { o.end % n } else { o.start % n };
        out.sort_by_key(|o| (o.forward, stop(o), std::cmp::Reverse(o.end - o.start)));
        out.dedup_by_key(|o| (o.forward, stop(o)));
    }
    out.sort_by_key(|o| o.start);
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn finds_orfs_on_both_strands() {
        //                 ATG AAA AAA TAA  (forward, 4 codons at 2..14)
        let fwd = b"CCATGAAAAAATAACC";
        assert_eq!(find_orfs(fwd, 4, false), vec![Orf { start: 2, end: 14, forward: true }]);
        assert!(find_orfs(fwd, 5, false).is_empty());
        // The same ORF written on the bottom strand.
        let rev = revcomp(fwd);
        assert_eq!(find_orfs(&rev, 4, false), vec![Orf { start: 2, end: 14, forward: false }]);
        // --- the ORF above is centred, so its mirror lands on the same place; with unequal flanks (19 bp,
        // ORF 2..14) top position p is 18 - p on the reverse complement, so the ORF covers 5..17 there ---
        let fwd = b"CCATGAAAAAATAACCCCC";
        assert_eq!(find_orfs(&revcomp(fwd), 4, false), vec![Orf { start: 5, end: 17, forward: false }]);
    }

    #[test]
    fn alternative_starts_and_other_tables() {
        // GTG AAA AAA TAA: a GTG start (table 11) is only an ORF when alternative starts are allowed.
        let seq = b"CCGTGAAAAAATAACC";
        let t11 = code::by_id(11).unwrap();
        assert!(find_orfs_with(seq, 4, false, t11, false).is_empty(), "ATG only by default");
        assert_eq!(find_orfs_with(seq, 4, false, t11, true), vec![Orf { start: 2, end: 14, forward: true }]);
        // Table 1 has no GTG start: nothing even with alternatives allowed.
        assert!(find_orfs_with(seq, 4, false, code::standard(), true).is_empty());
        // TTG is a start in table 1; the longest ORF starts at the first start codon after the previous stop.
        let ttg = b"CCTTGAAAGTGAAATAACC";
        let found = find_orfs_with(ttg, 4, false, code::standard(), true);
        assert_eq!(found, vec![Orf { start: 2, end: 17, forward: true }]);
        // Table 4: TGA is tryptophan, so ATG TGA AAA TAA is one ORF of 4 codons (TGA does not stop it).
        let t4 = code::by_id(4).unwrap();
        assert_eq!(find_orfs_with(b"ATGTGAAAATAA", 4, false, t4, false).len(), 1);
        assert!(find_orfs_with(b"ATGTGAAAATAA", 4, false, code::standard(), false).is_empty());
    }

    #[test]
    fn finds_an_orf_that_crosses_the_origin_of_a_ring() {
        // ATG AAA | AAA TAA: the sequence starts in the middle of the ORF and the ATG is at the end.
        let ring = format!("AAATAA{}ATGAAA", "C".repeat(30));
        assert!(find_orfs(ring.as_bytes(), 4, false).is_empty(), "a line has no such ORF");
        let found = find_orfs(ring.as_bytes(), 4, true);
        assert_eq!(found, vec![Orf { start: 36, end: 36 + 12, forward: true }]);
        // The reverse complement of the ring has its ORF on the bottom strand. (This ORF has 6 bases on each
        // side of the origin, so its mirrored place is 36 again: the asymmetric ring below pins the mapping.)
        let rc = revcomp(ring.as_bytes());
        assert_eq!(find_orfs(&rc, 4, true), vec![Orf { start: 36, end: 48, forward: false }]);
        // --- asymmetric ring (45 bp): ATG AAA at 39..45, then AAA AAA TAA at 0..9, so the ORF is 39..54; on the
        // reverse complement p becomes 44 - p, which puts it at 36..51 on the bottom strand ---
        let ring = format!("AAAAAATAA{}ATGAAA", "C".repeat(30));
        assert_eq!(find_orfs(ring.as_bytes(), 4, true), vec![Orf { start: 39, end: 54, forward: true }]);
        assert_eq!(find_orfs(&revcomp(ring.as_bytes()), 4, true), vec![Orf { start: 36, end: 51, forward: false }]);
        // An ORF inside a ring is reported once, not twice.
        assert_eq!(find_orfs(b"CCATGAAAAAATAACC", 4, true).len(), 1);
        // --- an ATG just after the origin, inside an ORF that started before it, is not a second ORF: ATG AAA TAA at 0..9 and
        // ATG AAA at 39..45 read round the ring as one ORF 39..54 (its stop is the TAA at 6..9); the shorter 0..9 shares that stop
        // and used to be reported too. Same on the bottom strand (p becomes 44 - p: 36..51) ---
        let ring = format!("ATGAAATAA{}ATGAAA", "C".repeat(30));
        assert_eq!(find_orfs(ring.as_bytes(), 3, true), vec![Orf { start: 39, end: 54, forward: true }]);
        assert_eq!(find_orfs(&revcomp(ring.as_bytes()), 3, true), vec![Orf { start: 36, end: 51, forward: false }]);
        // as a line, both ATG..stop stretches that fit are separate ORFs: only 0..9 here (39..45 has no stop)
        assert_eq!(find_orfs(ring.as_bytes(), 3, false), vec![Orf { start: 0, end: 9, forward: true }]);
    }
}
