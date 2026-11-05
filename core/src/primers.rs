//! Primer maths: melting temperature and where a primer binds a template.
//!
//! Tm uses the nearest-neighbour model with SantaLucia (1998) "unified" parameters, the
//! same model as Primer3 and Biopython's `Tm_NN` (DNA_NN3 table, salt correction 5):
//! - every adjacent base pair (e.g. "AG") adds a fixed ΔH and ΔS,
//! - each end adds an initiation term (different for A·T and G·C ends),
//! - salt shifts ΔS; Mg²⁺ is converted to an equivalent Na⁺ (von Ahsen 2001),
//! - Tm = ΔH / (ΔS + R·ln(C)) with C = primer concentration / 4 (or / 1 if the primer
//!   is its own reverse complement).

/// Reaction conditions for Tm. Defaults are Primer3's.
#[derive(Clone, Copy, Debug)]
pub struct Conditions {
    pub na_mm: f64,     // monovalent cations (Na⁺, K⁺), mM
    pub mg_mm: f64,     // Mg²⁺, mM
    pub dntp_mm: f64,   // dNTPs, mM (they bind Mg²⁺)
    pub primer_nm: f64, // primer concentration, nM
}

// `impl Default` = what `Conditions::default()` returns.
impl Default for Conditions {
    fn default() -> Self {
        Conditions { na_mm: 50.0, mg_mm: 1.5, dntp_mm: 0.6, primer_nm: 50.0 }
    }
}

const R: f64 = 1.987; // gas constant, cal/(K·mol)

/// ΔH (kcal/mol) and ΔS (cal/(K·mol)) of one nearest-neighbour pair, read 5'->3' on the
/// top strand. "CA" and "TG" are the same duplex step seen from either strand.
fn nn(pair: &[u8]) -> (f64, f64) {
    match pair {
        b"AA" | b"TT" => (-7.9, -22.2),
        b"AT" => (-7.2, -20.4),
        b"TA" => (-7.2, -21.3),
        b"CA" | b"TG" => (-8.5, -22.7),
        b"GT" | b"AC" => (-8.4, -22.4),
        b"CT" | b"AG" => (-7.8, -21.0),
        b"GA" | b"TC" => (-8.2, -22.2),
        b"CG" => (-10.6, -27.2),
        b"GC" => (-9.8, -24.4),
        b"GG" | b"CC" => (-8.0, -19.9),
        // `unreachable!` = crash with a message; tm() only passes A/C/G/T pairs.
        _ => unreachable!("not an ACGT pair"),
    }
}

/// Melting temperature in °C, or `None` if the sequence is shorter than 2 bases or
/// contains anything but A, C, G, T (U is read as T).
pub fn tm(seq: &[u8], c: &Conditions) -> Option<f64> {
    let s: Vec<u8> = seq.iter().map(|b| match b.to_ascii_uppercase() { b'U' => b'T', x => x }).collect();
    if s.len() < 2 || s.iter().any(|b| !b"ACGT".contains(b)) {
        return None;
    }
    let (mut dh, mut ds) = (0.0, 0.0);
    for end in [s[0], s[s.len() - 1]] {
        let (h, e) = if end == b'A' || end == b'T' { (2.3, 4.1) } else { (0.1, -2.8) };
        dh += h;
        ds += e;
    }
    // `windows(2)` yields every overlapping pair: "ACGT" -> "AC", "CG", "GT".
    for pair in s.windows(2) {
        let (h, e) = nn(pair);
        dh += h;
        ds += e;
    }
    let symmetric = s == revcomp(&s);
    if symmetric {
        ds += -1.4;
    }
    let na_eq = c.na_mm + 120.0 * (c.mg_mm - c.dntp_mm).max(0.0).sqrt();
    ds += 0.368 * (s.len() - 1) as f64 * (na_eq / 1000.0).ln();
    let ct = c.primer_nm * 1e-9;
    let k = if symmetric { ct } else { ct / 4.0 };
    Some(dh * 1000.0 / (ds + R * k.ln()) - 273.15)
}

/// Where a primer anneals on the template. Coordinates are on the top strand, 0-based:
/// bases `[start, start + len)`; on a circular template `start + len` may pass the end
/// (the site wraps through the origin).
#[derive(Debug, PartialEq)]
pub struct Site {
    pub start: usize,
    pub len: usize,
    /// true: the primer has the top-strand sequence (it extends to the right).
    pub forward: bool,
    pub mismatches: usize,
}

/// Default number of bases at the primer's 3' end that must match (with up to
/// `max_mismatches`) for a binding site. The annealed region then extends toward the 5'
/// end while bases keep matching; anything left over is a 5' tail (restriction site,
/// overlap, ...) that does not bind. A shorter seed (e.g. 12) finds weaker 3' homology:
/// places where the primer could mis-prime.
pub const SEED: usize = 15;
/// The last 3 bases at the 3' end must match exactly: the polymerase starts there.
const EXACT_3P: usize = 3;

/// All binding sites of `primer` on either strand.
/// shortcut: brute-force scan, ~0.2 s on a 5 Mb genome; use an index (FM/k-mer) if
/// whole-genome primer lists become slow.
pub fn binding_sites(template: &[u8], primer: &[u8], circular: bool, seed: usize, max_mismatches: usize) -> Vec<Site> {
    let p: Vec<u8> = primer.iter().map(|b| match b.to_ascii_uppercase() { b'U' => b'T', x => x }).collect();
    let n = template.len();
    if p.len() < 8 || n == 0 || seed < 4 {
        return Vec::new(); // too short to bind specifically
    }
    let seed = p.len().min(seed).min(n);
    // Template base at position i (wrapping on circular templates); None outside.
    let base = |i: i64| -> Option<u8> {
        if circular {
            Some(template[i.rem_euclid(n as i64) as usize].to_ascii_uppercase())
        } else if i >= 0 && (i as usize) < n {
            Some(template[i as usize].to_ascii_uppercase())
        } else {
            None
        }
    };
    let last_start = if circular { n } else { n + 1 - seed };
    let mut sites = Vec::new();
    // Case says nothing here: a .dna file writes annealing bases in lowercase too. A tail is
    // whatever the extension below does not reach.
    let reach = p.len() - seed; // bases the extension may cover

    // Forward primer: its 3' seed appears as-is on the top strand; extend leftwards.
    let fwd_seed = &p[p.len() - seed..];
    // Reverse primer: its reverse complement appears on the top strand, 3' end first;
    // extend rightwards.
    let rc = revcomp(&p);
    let rev_seed = &rc[..seed];
    for i in 0..last_start {
        // Forward: mismatches allowed except in the last EXACT_3P bases of the seed.
        if let Some(mm) = seed_mismatches(&base, i as i64, fwd_seed, max_mismatches, |k| k + EXACT_3P >= seed) {
            let (ext, ext_mm) = extend(reach, |j| Some((base(i as i64 - 1 - j as i64)?, p[p.len() - seed - 1 - j])));
            let len = (seed + ext).min(n);
            let start = (i as i64 - ext as i64).rem_euclid(n as i64) as usize;
            sites.push(Site { start, len, forward: true, mismatches: mm + ext_mm });
        }
        // Reverse: the primer's 3' end is the first EXACT_3P bases of `rev_seed`.
        if let Some(mm) = seed_mismatches(&base, i as i64, rev_seed, max_mismatches, |k| k < EXACT_3P) {
            let (ext, ext_mm) = extend(reach, |j| Some((base((i + seed + j) as i64)?, rc[seed + j])));
            sites.push(Site { start: i, len: (seed + ext).min(n), forward: false, mismatches: mm + ext_mm });
        }
    }
    sites
}

/// How far the annealed region goes past the seed toward the primer's 5' end, and the
/// mismatches in it. X-drop extension: +1 per matching base, -2 per mismatch, ending at
/// the best score and stopping once the score is 3 below it. So a lone mismatch followed
/// by matching bases is crossed (they still anneal), two in a row end the site (the rest
/// is a 5' tail). `pair(j)` = (template base, primer base) j steps past the seed, None
/// past the template's end.
fn extend(reach: usize, pair: impl Fn(usize) -> Option<(u8, u8)>) -> (usize, usize) {
    let (mut score, mut best, mut best_len, mut mm, mut best_mm) = (0i64, 0i64, 0, 0, 0);
    for j in 0..reach {
        let Some((t, q)) = pair(j) else { break };
        if t == q {
            score += 1;
        } else {
            score -= 2;
            mm += 1;
        }
        if score > best {
            (best, best_len, best_mm) = (score, j + 1, mm);
        }
        if best - score > 3 {
            break;
        }
    }
    (best_len, best_mm)
}

/// Mismatch count between `seed` and the template at `at`, or None if it is over the
/// limit or a mismatch falls where `exact(k)` requires a match.
/// `impl Fn(...)` = "any function or closure with this signature".
fn seed_mismatches(
    base: &impl Fn(i64) -> Option<u8>,
    at: i64,
    seed: &[u8],
    max: usize,
    exact: impl Fn(usize) -> bool,
) -> Option<usize> {
    let mut mm = 0;
    for (k, &b) in seed.iter().enumerate() {
        if base(at + k as i64)? != b {
            if exact(k) || mm == max {
                return None;
            }
            mm += 1;
        }
    }
    Some(mm)
}

/// ΔG at 37 °C (kcal/mol) of a perfectly paired duplex formed by `s` and its complement:
/// nearest-neighbour ΔH/ΔS with the same salt correction as `tm`. `init`: add the duplex
/// initiation terms (two separate strands); hairpin stems use a loop term instead.
fn stem_dg(s: &[u8], c: &Conditions, init: bool) -> f64 {
    let (mut dh, mut ds) = (0.0, 0.0);
    if init {
        for end in [s[0], s[s.len() - 1]] {
            let (h, e) = if end == b'A' || end == b'T' { (2.3, 4.1) } else { (0.1, -2.8) };
            dh += h;
            ds += e;
        }
    }
    for pair in s.windows(2) {
        let (h, e) = nn(pair);
        dh += h;
        ds += e;
    }
    let na_eq = c.na_mm + 120.0 * (c.mg_mm - c.dntp_mm).max(0.0).sqrt();
    ds += 0.368 * (s.len() - 1) as f64 * (na_eq / 1000.0).ln();
    dh - 310.15 * ds / 1000.0
}

fn pairs(a: u8, b: u8) -> bool {
    matches!((a, b), (b'A', b'T') | (b'T', b'A') | (b'C', b'G') | (b'G', b'C'))
}

fn clean(seq: &[u8]) -> Vec<u8> {
    seq.iter().map(|b| match b.to_ascii_uppercase() { b'U' => b'T', x => x }).collect()
}

/// Most stable primer-dimer between `a` and `b` (both 5'->3'; pass the same sequence twice
/// for a self-dimer). Returns (ΔG of the most stable dimer, ΔG of the most stable dimer
/// that pairs a 3' end), in kcal/mol; 0.0 when nothing pairs. 3'-end dimers matter most:
/// the polymerase can extend them.
///
/// shortcut: only perfectly paired stretches of ≥ 2 bp are scored (no internal loops,
/// mismatches or dangling ends), like quick primer checkers; replace with a full
/// thermodynamic alignment (Primer3's thal) if borderline cases matter.
pub fn dimer_dg(a: &[u8], b: &[u8], c: &Conditions) -> (f64, f64) {
    let (a, b) = (clean(a), clean(b));
    // Antiparallel: a (5'->3') against b read 3'->5'.
    let br: Vec<u8> = b.iter().rev().copied().collect();
    let (mut best, mut best_3p) = (0.0f64, 0.0f64);
    // a[i] faces br[i - shift]; try every overlap.
    for shift in -(br.len() as i64) + 1..a.len() as i64 {
        let mut run_start: Option<usize> = None;
        for i in 0..=a.len() {
            let j = i as i64 - shift;
            let ok = i < a.len() && j >= 0 && (j as usize) < br.len() && pairs(a[i], br[j as usize]);
            match (ok, run_start) {
                (true, None) => run_start = Some(i),
                (false, Some(start)) => {
                    run_start = None;
                    if i - start >= 2 {
                        let dg = stem_dg(&a[start..i], c, true);
                        best = best.min(dg);
                        // a's 3' end is a[len-1]; b's 3' end is br[0].
                        let touches_3p = i == a.len() || start as i64 - shift == 0;
                        if touches_3p {
                            best_3p = best_3p.min(dg);
                        }
                    }
                }
                _ => {}
            }
        }
    }
    (best, best_3p)
}

/// Hairpin loop initiation ΔG37 (kcal/mol), SantaLucia & Hicks 2004; longer loops grow
/// logarithmically.
fn loop_dg(n: usize) -> f64 {
    const TABLE: [f64; 8] = [3.5, 3.5, 3.3, 4.0, 4.2, 4.3, 4.5, 4.6]; // loops of 3..=10
    if n <= 10 { TABLE[n - 3] } else { 4.6 + 2.44 * R * 310.15 / 1000.0 * (n as f64 / 10.0).ln() }
}

/// Most stable hairpin: (ΔG, ΔG of the most stable hairpin whose stem includes the 3'-end
/// base), kcal/mol; 0.0 when none is stable. Loops of ≥ 3 bases, stems of ≥ 2 bp.
/// shortcut: perfectly paired stems only, as `dimer_dg`.
pub fn hairpin_dg(seq: &[u8], c: &Conditions) -> (f64, f64) {
    let s = clean(seq);
    let n = s.len();
    let (mut best, mut best_3p) = (0.0f64, 0.0f64);
    // Stem: s[i + t] pairs with s[j - t] for t = 0..k; the loop is s[i + k ..= j - k].
    for i in 0..n {
        for j in (i + 5)..n {
            // Only maximal stems: skip if the pair just outside also pairs.
            if !pairs(s[i], s[j]) || (i > 0 && j + 1 < n && pairs(s[i - 1], s[j + 1])) {
                continue;
            }
            let mut k = 0;
            while i + k < j - k && pairs(s[i + k], s[j - k]) && (j - k) - (i + k) - 1 >= 3 {
                k += 1;
            }
            if k < 2 {
                continue;
            }
            let dg = stem_dg(&s[i..i + k], c, false) + loop_dg((j - k + 1) - (i + k));
            best = best.min(dg);
            if j == n - 1 {
                best_3p = best_3p.min(dg);
            }
        }
    }
    (best, best_3p)
}

pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|b| match b.to_ascii_uppercase() {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            x => x,
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tm_is_sane() {
        // Exact values are checked against Biopython's Tm_NN (same model) in
        // tests/test_primers.py; here only plausibility and input checks.
        let c = Conditions::default();
        let t = tm(b"AGCGGATAACAATTTCACACAGGA", &c).unwrap();
        assert!(t > 55.0 && t < 65.0, "{t}");
        assert!(tm(b"ACGTN", &c).is_none());
        assert!(tm(b"A", &c).is_none());
        // More salt or more primer = higher Tm.
        let salty = Conditions { na_mm: 200.0, ..c };
        assert!(tm(b"AGCGGATAACAATTTCACACAGGA", &salty).unwrap() > t);
    }

    #[test]
    fn annealing_extends_across_a_lone_mismatch() {
        //                 0         1         2         3         4
        //                 0123456789012345678901234567890123456789012345
        let template = b"TTTTTTTTTTACGTACGGATCCAGTTGACCATGGTTTTTTTTTTTT";
        // Forward primer = top strand 10..34 but with base 14 wrong: the 15-nt seed is
        // 19..34; past it, the lone mismatch is crossed because the 4 bases after it
        // (10..14) still anneal (-2 for the mismatch, +4 after it).
        let mut fwd = template[10..34].to_vec();
        fwd[4] = b'C';
        assert_eq!(binding_sites(template, &fwd, false, SEED, 0),
                   vec![Site { start: 10, len: 24, forward: true, mismatches: 1 }]);
        // Only 2 good bases after a mismatch: not worth it, they stay a tail.
        let mut short = template[10..34].to_vec();
        short[2] = b'C';
        assert_eq!(binding_sites(template, &short, false, SEED, 0)[0].start, 13);
        // Reverse strand: its seed is on the left (10..25), the extension goes right.
        let mut top = template[10..34].to_vec();
        top[19] = b'A'; // base 29
        let rev = revcomp(&top);
        assert_eq!(binding_sites(template, &rev, false, SEED, 0),
                   vec![Site { start: 10, len: 24, forward: false, mismatches: 1 }]);
        // Two wrong bases in a row stop it: they and what lies beyond are a tail.
        let mut two = template[10..34].to_vec();
        two[5] = b'A';
        two[6] = b'A';
        let sites = binding_sites(template, &two, false, SEED, 0);
        assert_eq!((sites[0].start, sites[0].len, sites[0].mismatches), (17, 17, 0));
    }

    #[test]
    fn finds_sites_on_both_strands_with_tails_and_mismatches() {
        //                 0         1         2         3         4
        //                 0123456789012345678901234567890123456789012345
        let template = b"TTTTTTTTTTACGTACGGATCCAGTTGACCATGGTTTTTTTTTTTT";
        // Forward primer = top strand 10..34 with a 6-base 5' tail.
        let fwd = b"GAATTCACGTACGGATCCAGTTGACCATGG";
        let sites = binding_sites(template, fwd, false, SEED, 0);
        assert_eq!(sites, vec![Site { start: 10, len: 24, forward: true, mismatches: 0 }]);
        // Reverse primer: reverse complement of 14..34.
        let rev = revcomp(&template[14..34]);
        let sites = binding_sites(template, &rev, false, SEED, 0);
        assert_eq!(sites, vec![Site { start: 14, len: 20, forward: false, mismatches: 0 }]);
        // One mismatch in the seed (not at the 3' end): only found when allowed.
        let mut one_off = fwd.to_vec();
        one_off[20] = b'A'; // was T
        assert!(binding_sites(template, &one_off, false, SEED, 0).is_empty());
        assert_eq!(binding_sites(template, &one_off, false, SEED, 1)[0].mismatches, 1);
        // A mismatch at the 3'-terminal base never binds.
        let mut bad_end = fwd.to_vec();
        *bad_end.last_mut().unwrap() = b'A';
        assert!(binding_sites(template, &bad_end, false, SEED, 2).is_empty());
    }

    #[test]
    fn site_across_the_origin() {
        let template = b"GGATCCAGTTTTTTTTTTTTTTTTTTACGTACGGATCC";
        let primer = b"ACGTACGGATCCGGATCCAG"; // last 8 bases continue at the start
        let sites = binding_sites(template, primer, true, SEED, 0);
        assert_eq!(sites, vec![Site { start: 26, len: 20, forward: true, mismatches: 0 }]);
    }

    #[test]
    fn dimers_and_hairpins() {
        let c = Conditions::default();
        // A palindromic 3' end makes a strong 3' self-dimer.
        let (any, three) = dimer_dg(b"ACGTTTAGCGAATTCGC", b"ACGTTTAGCGAATTCGC", &c);
        assert!(three < -5.0 && any <= three, "{any} {three}");
        // Poly-A cannot pair with itself.
        assert_eq!(dimer_dg(b"AAAAAAAAAAAAAAAA", b"AAAAAAAAAAAAAAAA", &c), (0.0, 0.0));
        // A primer and its exact complement: the full duplex, strongest possible.
        let p = b"AGCGGATAACAATTTCACACAGGA";
        let full = dimer_dg(p, &revcomp(p), &c).0;
        assert!(full < -20.0, "{full}"); // ~ -1 kcal/mol per base pair
        // Hairpin with a GCGCGC / GCGCGC stem and a TTTT loop: the stem reaches the 3' end, so both values are the same.
        let (hp, hp3) = hairpin_dg(b"GCGCGCTTTTGCGCGC", &c);
        assert!(hp < -3.0 && hp3 == hp, "{hp} {hp3}");
        assert_eq!(hairpin_dg(b"AAAAAAAAAAAAAAA", &c), (0.0, 0.0));
    }
}
