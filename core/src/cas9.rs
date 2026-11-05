//! Cas9 guide search: where a guide RNA's target (the 20 bases it pairs with, then the PAM) occurs in a
//! sequence, allowing mismatches. Used to count the other places a guide could cut, in the open sequence
//! and in a reference genome.
//!
//! Vocabulary, for readers who do not work with CRISPR:
//! - the *guide* (protospacer) is the ~20 bases the guide RNA copies; Cas9 cuts where the DNA matches it;
//! - the *PAM* is a short pattern (NGG for the usual SpCas9) that must sit right after the protospacer
//!   on the same strand, or Cas9 ignores the site;
//! - a *mismatch* is a base of the guide that differs from the DNA; Cas9 tolerates a few, most of all far
//!   from the PAM.
//!
//! A "site" is protospacer + PAM, on either strand of the sequence.

/// One place a guide matches.
#[derive(Debug, Clone, PartialEq)]
pub struct GuideSite {
    /// First base of the site on the top strand (0-based).
    pub start: usize,
    /// One past the last base. On a circular sequence this can pass the sequence length (it wraps).
    pub end: usize,
    /// True if the guide reads along the top strand (protospacer then PAM, left to right).
    /// False if it reads along the bottom strand: on the top strand the PAM is then on the left.
    pub forward: bool,
    /// How many of the guide's bases differ from the DNA.
    pub mismatches: usize,
    /// Which bases differ: bit `i` is set when guide base `i` (0 = the 5' end, far from the PAM) mismatches.
    pub mask: u32,
}

/// The complement of one base; anything that is not A, C, G or T is left alone.
fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        other => other,
    }
}

/// The other strand, read in its own 5'->3' direction (upper case).
fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|b| complement(b.to_ascii_uppercase())).collect()
}

/// Does the DNA base `base` (upper case) fit the IUPAC letter `code` of a PAM pattern?
/// "N" is any base, "R" is A or G, and so on. A base that is not A, C, G or T fits nothing.
fn iupac_fits(code: u8, base: u8) -> bool {
    // `matches!` is Rust's "does this value look like one of these patterns?" test.
    if !matches!(base, b'A' | b'C' | b'G' | b'T') {
        return false;
    }
    match code.to_ascii_uppercase() {
        b'A' | b'C' | b'G' | b'T' => code.to_ascii_uppercase() == base,
        b'N' => true,
        b'R' => matches!(base, b'A' | b'G'),
        b'Y' => matches!(base, b'C' | b'T'),
        b'S' => matches!(base, b'C' | b'G'),
        b'W' => matches!(base, b'A' | b'T'),
        b'K' => matches!(base, b'G' | b'T'),
        b'M' => matches!(base, b'A' | b'C'),
        b'B' => base != b'A',
        b'D' => base != b'C',
        b'H' => base != b'G',
        b'V' => base != b'T',
        _ => false,
    }
}

/// Every site of each guide in `template`, on both strands, with at most `max_mismatches` mismatches in the
/// protospacer and a PAM that fits `pam` exactly. All guides must have the same length (the first one's).
/// Returns one list per guide, in the order given.
///
/// How it works: per strand, first find every place where the PAM pattern fits (about one base in 16 for
/// NGG), then compare each guide with the bases in front of those places, giving up on a place as soon as it
/// has too many mismatches. That makes a whole bacterial genome take a few tens of milliseconds per guide.
pub fn find_guide_sites(
    template: &[u8],
    circular: bool,
    guides: &[Vec<u8>],
    pam: &[u8],
    max_mismatches: usize,
) -> Vec<Vec<GuideSite>> {
    find_sites(template, circular, guides, pam, max_mismatches, false)
}

/// The same search for either kind of nuclease. `pam_first` says where the PAM sits on the strand the guide reads
/// along: false (Cas9: protospacer, then PAM) or true (Cas12a: PAM, then protospacer, 5' of it). The site is the
/// same length either way and is reported in top-strand coordinates, so only the two offsets below differ.
/// The mismatch mask counts guide bases from the guide's own 5' end, which is next to the PAM for `pam_first`.
pub fn find_sites(
    template: &[u8],
    circular: bool,
    guides: &[Vec<u8>],
    pam: &[u8],
    max_mismatches: usize,
    pam_first: bool,
) -> Vec<Vec<GuideSite>> {
    let mut out: Vec<Vec<GuideSite>> = vec![Vec::new(); guides.len()];
    let n = template.len();
    if n == 0 || guides.is_empty() || pam.is_empty() {
        return out;
    }
    let len = guides[0].len();
    let span = len + pam.len(); // bases in a site
    // Where, counted from the start of a site, the PAM and the protospacer begin.
    let (pam_at, guide_at) = if pam_first { (0, pam.len()) } else { (len, 0) };
    if len == 0 || len > 32 || n + 1 < span {
        return out; // the mask holds 32 bases; a sequence shorter than a site has none
    }
    // The guides as upper-case bytes, so a lower-case guide still matches.
    let guides_up: Vec<Vec<u8>> = guides.iter().map(|g| g.iter().map(|b| b.to_ascii_uppercase()).collect()).collect();
    let top: Vec<u8> = template.iter().map(|b| b.to_ascii_uppercase()).collect();

    for forward in [true, false] {
        // Scan the strand the guide reads along: the top strand as it is, or its reverse complement.
        let mut strand = if forward { top.clone() } else { revcomp(&top) };
        if circular {
            // Copy the start onto the end so a site can run across the origin.
            let head: Vec<u8> = strand[..(span - 1).min(n)].to_vec();
            strand.extend(head);
        }
        // Where a site could start on this strand.
        let starts = if circular { n } else { n + 1 - span };
        // First pass: where does the PAM fit?
        let mut pam_here = Vec::new();
        for i in 0..starts {
            if i + span > strand.len() {
                break;
            }
            if pam.iter().enumerate().all(|(k, &code)| iupac_fits(code, strand[i + pam_at + k])) {
                pam_here.push(i);
            }
        }
        // Second pass: compare each guide with the bases before each PAM.
        for (g, guide) in guides_up.iter().enumerate() {
            if guide.len() != len {
                continue; // a guide of another length cannot be compared position by position
            }
            for &i in &pam_here {
                let (mut mismatches, mut mask) = (0usize, 0u32);
                for k in 0..len {
                    // A base that is not A/C/G/T never equals a guide base, so it counts as a mismatch.
                    if strand[i + guide_at + k] != guide[k] {
                        mismatches += 1;
                        mask |= 1 << k;
                        if mismatches > max_mismatches {
                            break;
                        }
                    }
                }
                if mismatches <= max_mismatches {
                    // Back to top-strand coordinates. A reverse-strand site that starts at `i` on the
                    // reverse complement occupies [n - i - span, n - i) on the top strand.
                    let start = if forward { i } else { (n as i64 - (i + span) as i64).rem_euclid(n as i64) as usize };
                    out[g].push(GuideSite { start, end: start + span, forward, mismatches, mask });
                }
            }
        }
    }
    // Same order whatever the strand: by position, then forward before reverse.
    for sites in &mut out {
        sites.sort_by_key(|s| (s.start, !s.forward));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn guide(s: &str) -> Vec<Vec<u8>> {
        vec![s.as_bytes().to_vec()]
    }

    #[test]
    fn finds_a_site_with_its_pam_on_the_top_strand() {
        // 5 filler bases, a 20-base protospacer, TGG (a NGG PAM), 5 filler bases.
        let g = "ACGTTGCAAGTCCATGCATG";
        let t = format!("TTTTT{g}TGGCCCCC");
        let hits = find_guide_sites(t.as_bytes(), false, &guide(g), b"NGG", 0);
        assert_eq!(hits[0].len(), 1);
        let h = &hits[0][0];
        assert_eq!((h.start, h.end, h.forward, h.mismatches, h.mask), (5, 28, true, 0, 0));
    }

    #[test]
    fn finds_a_site_on_the_bottom_strand() {
        let g = "ACGTTGCAAGTCCATGCATG";
        // The bottom strand reads protospacer + PAM, so the top strand holds the reverse complement.
        let site = String::from_utf8(revcomp(format!("{g}AGG").as_bytes())).unwrap();
        let t = format!("TTTTT{site}CCCCC");
        let hits = find_guide_sites(t.as_bytes(), false, &guide(g), b"NGG", 0);
        assert_eq!(hits[0].len(), 1);
        let h = &hits[0][0];
        assert_eq!((h.start, h.end, h.forward), (5, 28, false), "PAM is on the left of the top strand");
        assert_eq!(&t.as_bytes()[h.start..h.start + 3], b"CCT", "reverse complement of AGG");
    }

    #[test]
    fn needs_the_pam() {
        let g = "ACGTTGCAAGTCCATGCATG";
        let t = format!("TTTTT{g}TAACCCCC"); // TAA is no NGG
        assert!(find_guide_sites(t.as_bytes(), false, &guide(g), b"NGG", 3)[0].is_empty());
        // ...but it is a perfect NNGRRT site for another nuclease when the PAM is written that way.
        let t2 = format!("TTTTT{g}ATGAAT");
        assert_eq!(find_guide_sites(t2.as_bytes(), false, &guide(g), b"NNGRRT", 0)[0].len(), 1);
    }

    #[test]
    fn counts_mismatches_and_records_where() {
        let g = "ACGTTGCAAGTCCATGCATG";
        // Change guide base 0 (far from the PAM) and base 19 (next to it).
        let mut bent: Vec<u8> = g.bytes().collect();
        bent[0] = b'T';
        bent[19] = b'C';
        let t = format!("TTTTT{}CGGCCCCC", String::from_utf8(bent).unwrap());
        assert!(find_guide_sites(t.as_bytes(), false, &guide(g), b"NGG", 1)[0].is_empty(), "two mismatches, one allowed");
        let hits = find_guide_sites(t.as_bytes(), false, &guide(g), b"NGG", 2);
        assert_eq!(hits[0].len(), 1);
        assert_eq!(hits[0][0].mismatches, 2);
        assert_eq!(hits[0][0].mask, 1 | (1 << 19));
    }

    #[test]
    fn a_site_can_run_across_the_origin_of_a_circle() {
        let g = "ACGTTGCAAGTCCATGCATG";
        let site = format!("{g}TGG"); // 23 bases
        // Put the last 5 bases of the site at the start of the sequence and the rest at the end.
        let (tail, head) = site.split_at(18);
        let t = format!("{head}CCCCCCCCCCCCCCCCCCCC{tail}");
        let circle = find_guide_sites(t.as_bytes(), true, &guide(g), b"NGG", 0);
        assert_eq!(circle[0].len(), 1);
        assert_eq!(circle[0][0].start, t.len() - 18);
        assert!(circle[0][0].end > t.len(), "the site wraps past the end");
        assert!(find_guide_sites(t.as_bytes(), false, &guide(g), b"NGG", 0)[0].is_empty(), "a line has no such site");
    }

    #[test]
    fn several_guides_and_lower_case() {
        let a = "ACGTTGCAAGTCCATGCATG";
        let b = "GGATCCTTAAGGCCTTAAGC";
        let t = format!("tttttt{}TGGaaaaaa{}CGGtttttt", a.to_lowercase(), b);
        let hits = find_guide_sites(t.as_bytes(), false, &[a.as_bytes().to_vec(), b.as_bytes().to_vec()], b"NGG", 0);
        assert_eq!((hits[0].len(), hits[1].len()), (1, 1));
        assert!(hits[0][0].start < hits[1][0].start);
    }

    #[test]
    fn finds_a_site_whose_pam_comes_first() {
        // Cas12a: the PAM (TTTV) is on the 5' side. On the top strand: TTTA, then the 20-base protospacer.
        let g = "ACGTTGCAAGTCCATGCATG";
        let t = format!("CCCCCTTTA{g}CCCCC");
        let hits = find_sites(t.as_bytes(), false, &guide(g), b"TTTV", 0, true);
        assert_eq!(hits[0].len(), 1);
        let h = &hits[0][0];
        assert_eq!((h.start, h.end, h.forward, h.mismatches), (5, 29, true, 0));
        // The same text is no site when the PAM must come after the protospacer.
        assert!(find_sites(t.as_bytes(), false, &guide(g), b"TTTV", 0, false)[0].is_empty());
        // On the bottom strand the site reads reverse-complemented: protospacer first on the top strand, then the PAM's complement.
        let bottom = String::from_utf8(revcomp(format!("TTTC{g}").as_bytes())).unwrap();
        let t2 = format!("CCCCC{bottom}CCCCC");
        let hits = find_sites(t2.as_bytes(), false, &guide(g), b"TTTV", 0, true);
        assert_eq!(hits[0].len(), 1);
        assert_eq!((hits[0][0].start, hits[0][0].end, hits[0][0].forward), (5, 29, false));
    }

    #[test]
    fn mismatch_mask_counts_from_the_guides_5_prime_end_next_to_a_first_pam() {
        let g = "ACGTTGCAAGTCCATGCATG";
        let mut bent: Vec<u8> = g.bytes().collect();
        bent[0] = b'T'; // the base right after the PAM
        let t = format!("CCCCCTTTG{}CCCCC", String::from_utf8(bent).unwrap());
        let hits = find_sites(t.as_bytes(), false, &guide(g), b"TTTV", 1, true);
        assert_eq!((hits[0][0].mismatches, hits[0][0].mask), (1, 1));
    }
}
