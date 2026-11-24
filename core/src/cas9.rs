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
    pub mask: u64,
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

/// What a search found: the sites of each guide (in the order given), and for each guide how many more would have been
/// sites if an unknown base N of the DNA in the PAM were any base. Those are not reported (we cannot know), only counted.
pub struct Found {
    pub sites: Vec<Vec<GuideSite>>,
    pub skipped: Vec<usize>,
}

/// The search for either kind of nuclease, on both strands, with at most `max_mismatches` mismatches in the
/// protospacer and a PAM that fits `pam` exactly. `pam_first` says where the PAM sits on the strand the guide reads
/// along: false (Cas9: protospacer, then PAM) or true (Cas12a: PAM, then protospacer, 5' of it). The site is the
/// same length either way and is reported in top-strand coordinates, so only the two offsets below differ.
/// The mismatch mask counts guide bases from the guide's own 5' end, which is next to the PAM for `pam_first`.
/// All guides must have the same length, 1 to 64 (the mask is 64 bits), or the call fails: an empty answer
/// would read as "no site". A sequence shorter than one site has none, which is a true empty answer
/// (shortcut: a circular sequence shorter than a site could wrap onto itself and is not searched; plasmids never are).
///
/// How it works: per strand, first find every place where the PAM pattern fits (about one base in 16 for
/// NGG), then compare each guide with the bases in front of those places, giving up on a place as soon as it
/// has too many mismatches. That makes a whole bacterial genome take a few tens of milliseconds per guide.
pub fn find_sites(
    template: &[u8],
    circular: bool,
    guides: &[Vec<u8>],
    pam: &[u8],
    max_mismatches: usize,
    pam_first: bool,
) -> Result<Found, String> {
    let mut out = Found { sites: vec![Vec::new(); guides.len()], skipped: vec![0; guides.len()] };
    let n = template.len();
    if n == 0 || guides.is_empty() || pam.is_empty() {
        return Ok(out);
    }
    let len = guides[0].len();
    if len == 0 || len > 64 {
        return Err(format!("a guide must be 1 to 64 bases long, this one has {len}"));
    }
    if guides.iter().any(|g| g.len() != len) {
        return Err("all guides of one search must have the same length".to_string());
    }
    let span = len + pam.len(); // bases in a site
    // Where, counted from the start of a site, the PAM and the protospacer begin.
    let (pam_at, guide_at) = if pam_first { (0, pam.len()) } else { (len, 0) };
    if n + 1 < span {
        return Ok(out);
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
        // First pass: where does the PAM fit? And where would it fit if the N bases of the DNA were any base?
        let (mut pam_here, mut pam_maybe) = (Vec::new(), Vec::new());
        for i in 0..starts {
            if i + span > strand.len() {
                break;
            }
            let at = &strand[i + pam_at..i + pam_at + pam.len()];
            if pam.iter().zip(at).all(|(&code, &b)| iupac_fits(code, b)) {
                pam_here.push(i);
            } else if pam.iter().zip(at).all(|(&code, &b)| b == b'N' || iupac_fits(code, b)) {
                pam_maybe.push(i);
            }
        }
        // Second pass: compare each guide with the bases before each PAM. A base that is not A/C/G/T never
        // equals a guide base, so it counts as a mismatch. Some((mismatches, mask)) when within the limit.
        let compare = |i: usize, guide: &[u8]| {
            let (mut mismatches, mut mask) = (0usize, 0u64);
            for k in 0..len {
                if strand[i + guide_at + k] != guide[k] {
                    mismatches += 1;
                    mask |= 1 << k;
                    if mismatches > max_mismatches {
                        return None;
                    }
                }
            }
            Some((mismatches, mask))
        };
        for (g, guide) in guides_up.iter().enumerate() {
            for &i in &pam_here {
                if let Some((mismatches, mask)) = compare(i, guide) {
                    // Back to top-strand coordinates. A reverse-strand site that starts at `i` on the
                    // reverse complement occupies [n - i - span, n - i) on the top strand.
                    let start = if forward { i } else { (n as i64 - (i + span) as i64).rem_euclid(n as i64) as usize };
                    out.sites[g].push(GuideSite { start, end: start + span, forward, mismatches, mask });
                }
            }
            out.skipped[g] += pam_maybe.iter().filter(|&&i| compare(i, guide).is_some()).count();
        }
    }
    // Same order whatever the strand: by position, then forward before reverse.
    for sites in &mut out.sites {
        sites.sort_by_key(|s| (s.start, !s.forward));
    }
    Ok(out)
}

/// `find_sites` for a Cas9-type nuclease (protospacer, then PAM), sites only: for the tests.
#[cfg(test)]
fn find_guide_sites(template: &[u8], circular: bool, guides: &[Vec<u8>], pam: &[u8], max_mismatches: usize) -> Vec<Vec<GuideSite>> {
    find_sites(template, circular, guides, pam, max_mismatches, false).unwrap().sites
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
        // --- with 5 bases on both sides the reverse-strand index equals the top-strand start; unequal flanks
        // (5 / 10) tell them apart: top start 5, reverse-strand index 38 - 28 = 10 ---
        let t = format!("TTTTT{site}CCCCCCCCCC");
        let hits = find_guide_sites(t.as_bytes(), false, &guide(g), b"NGG", 0);
        assert_eq!(hits[0].iter().map(|s| (s.start, s.end, s.forward)).collect::<Vec<_>>(), vec![(5, 28, false)]);
        assert_eq!(&t.as_bytes()[5..8], b"CCT");
    }

    #[test]
    fn needs_the_pam() {
        let g = "ACGTTGCAAGTCCATGCATG";
        let t = format!("TTTTT{g}TAACCCCC"); // TAA is no NGG
        assert!(find_guide_sites(t.as_bytes(), false, &guide(g), b"NGG", 3)[0].is_empty());
        // ...but it is a perfect NNGRRT site for another nuclease when the PAM is written that way.
        let t2 = format!("TTTTT{g}ATGAAT");
        assert_eq!(find_guide_sites(t2.as_bytes(), false, &guide(g), b"NNGRRT", 0)[0].len(), 1);
        // --- R (A or G) refuses a C, and takes a G: ATGCAT is no NNGRRT, ATGGGT is one, at 5..31 ---
        let bad = format!("TTTTT{g}ATGCAT");
        assert!(find_guide_sites(bad.as_bytes(), false, &guide(g), b"NNGRRT", 0)[0].is_empty());
        let good = format!("TTTTT{g}ATGGGT");
        let hits = find_guide_sites(good.as_bytes(), false, &guide(g), b"NNGRRT", 0);
        assert_eq!(hits[0].iter().map(|s| (s.start, s.end, s.forward)).collect::<Vec<_>>(), vec![(5, 31, true)]);
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
        let hits = find_sites(t.as_bytes(), false, &guide(g), b"TTTV", 0, true).unwrap().sites;
        assert_eq!(hits[0].len(), 1);
        let h = &hits[0][0];
        assert_eq!((h.start, h.end, h.forward, h.mismatches), (5, 29, true, 0));
        // The same text is no site when the PAM must come after the protospacer.
        assert!(find_sites(t.as_bytes(), false, &guide(g), b"TTTV", 0, false).unwrap().sites[0].is_empty());
        // On the bottom strand the site reads reverse-complemented: protospacer first on the top strand, then the PAM's complement.
        let bottom = String::from_utf8(revcomp(format!("TTTC{g}").as_bytes())).unwrap();
        let t2 = format!("CCCCC{bottom}CCCCC");
        let hits = find_sites(t2.as_bytes(), false, &guide(g), b"TTTV", 0, true).unwrap().sites;
        assert_eq!(hits[0].len(), 1);
        assert_eq!((hits[0][0].start, hits[0][0].end, hits[0][0].forward), (5, 29, false));
        // --- unequal flanks (5 / 10): top start 5, reverse-strand index 39 - 29 = 10 ---
        let t2 = format!("CCCCC{bottom}CCCCCCCCCC");
        let hits = find_sites(t2.as_bytes(), false, &guide(g), b"TTTV", 0, true).unwrap().sites;
        assert_eq!(hits[0].iter().map(|s| (s.start, s.end, s.forward, s.mismatches)).collect::<Vec<_>>(), vec![(5, 29, false, 0)]);
    }

    #[test]
    fn mismatch_mask_counts_from_the_guides_5_prime_end_next_to_a_first_pam() {
        let g = "ACGTTGCAAGTCCATGCATG";
        let mut bent: Vec<u8> = g.bytes().collect();
        bent[0] = b'T'; // the base right after the PAM
        let t = format!("CCCCCTTTG{}CCCCC", String::from_utf8(bent).unwrap());
        let hits = find_sites(t.as_bytes(), false, &guide(g), b"TTTV", 1, true).unwrap().sites;
        assert_eq!((hits[0][0].mismatches, hits[0][0].mask), (1, 1));
    }

    #[test]
    fn a_guide_the_mask_cannot_hold_is_an_error_not_an_empty_answer() {
        let long = "A".repeat(40);
        // 40 bases now fit (a 64-bit mask): the site with its PAM is found and its mask sees base 39.
        let t = format!("CC{long}TGGCC");
        let found = find_sites(t.as_bytes(), false, &[long.as_bytes().to_vec()], b"NGG", 0, false).unwrap();
        assert_eq!(found.sites[0].len(), 1);
        let bent = format!("CC{}CTGGCC", &long[..39]);
        let found = find_sites(bent.as_bytes(), false, &[long.as_bytes().to_vec()], b"NGG", 1, false).unwrap();
        assert_eq!(found.sites[0][0].mask, 1 << 39);
        // 65 bases and guides of different lengths are refused, not answered with "no sites".
        assert!(find_sites(b"ACGT", false, &[vec![b'A'; 65]], b"NGG", 0, false).is_err());
        assert!(find_sites(t.as_bytes(), false, &[vec![b'A'; 20], vec![b'A'; 21]], b"NGG", 0, false).is_err());
        // A sequence shorter than one site is a true "none", not an error.
        assert!(find_sites(b"ACGT", false, &[vec![b'A'; 20]], b"NGG", 0, false).unwrap().sites[0].is_empty());
    }

    #[test]
    fn a_site_that_only_an_unknown_base_in_the_pam_hides_is_counted() {
        let g = "ACGTTGCAAGTCCATGCATG";
        // The PAM is A?G with ? unknown: NGG would fit if N were a G, so it is counted, not reported.
        let t = format!("TTTTT{g}ANGCCCCC");
        let found = find_sites(t.as_bytes(), false, &guide(g), b"NGG", 0, false).unwrap();
        assert_eq!((found.sites[0].len(), found.skipped[0]), (0, 1));
        // A known PAM is a site and nothing is skipped.
        let t = format!("TTTTT{g}AGGCCCCC");
        let found = find_sites(t.as_bytes(), false, &guide(g), b"NGG", 0, false).unwrap();
        assert_eq!((found.sites[0].len(), found.skipped[0]), (1, 0));
        // A PAM where a known base contradicts the pattern is no near-site.
        let t = format!("TTTTT{g}ANACCCCC");
        assert_eq!(find_sites(t.as_bytes(), false, &guide(g), b"NGG", 0, false).unwrap().skipped[0], 0);
    }
}
