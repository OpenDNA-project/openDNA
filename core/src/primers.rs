//! Primer maths: melting temperature and where a primer binds a template.
//!
//! Tm uses the nearest-neighbour model with SantaLucia (1998) "unified" parameters, the
//! same model as Primer3 and Biopython's `Tm_NN` (DNA_NN3 table, salt correction 5):
//! - every adjacent base pair (e.g. "AG") adds a fixed ΔH and ΔS,
//! - each end adds an initiation term (different for A·T and G·C ends),
//! - salt shifts ΔS; Mg²⁺ is converted to an equivalent Na⁺ (von Ahsen 2001),
//! - Tm = ΔH / (ΔS + R·ln(C)) with C = primer concentration / 4 (or / 1 if the primer
//!   is its own reverse complement).

use crate::search::pattern_mask;

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

/// Initiation (ΔH kcal/mol, ΔS cal/(K·mol)) of a duplex end: A·T ends differ from G·C ends.
fn end_term(b: u8) -> (f64, f64) {
    if b == b'A' || b == b'T' { (2.3, 4.1) } else { (0.1, -2.8) }
}

/// Salt correction of ΔS for a duplex of `len` bases (Mg²⁺ converted to equivalent Na⁺).
fn salt_ds(len: usize, c: &Conditions) -> f64 {
    let na_eq = c.na_mm + 120.0 * (c.mg_mm - c.dntp_mm).max(0.0).sqrt();
    0.368 * (len - 1) as f64 * (na_eq / 1000.0).ln()
}

/// Melting temperature in °C, or `None` if the sequence is shorter than 2 bases or
/// contains anything but A, C, G, T (U is read as T).
pub fn tm(seq: &[u8], c: &Conditions) -> Option<f64> {
    let s = clean(seq);
    if s.len() < 2 || s.iter().any(|b| !b"ACGT".contains(b)) {
        return None;
    }
    let (mut dh, mut ds) = (0.0, 0.0);
    for end in [s[0], s[s.len() - 1]] {
        let (h, e) = end_term(end);
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
    Some(melt(dh, ds, s.len(), symmetric, c))
}

/// Tm in °C of a duplex of `len` pairs with these nearest-neighbour sums (salt and strand concentration applied here).
fn melt(dh: f64, ds: f64, len: usize, symmetric: bool, c: &Conditions) -> f64 {
    let ct = c.primer_nm * 1e-9;
    let k = if symmetric { ct } else { ct / 4.0 };
    dh * 1000.0 / (ds + salt_ds(len, c) + R * k.ln()) - 273.15
}

/// Nearest-neighbour stacks with one mismatched pair, "XY/AB" = top 5'-XY-3' over bottom 3'-AB-5' (ΔH kcal/mol, ΔS cal/(K·mol)):
/// Allawi & SantaLucia 1997-98, Peyret et al. 1999 (Biopython's DNA_IMM1, without inosine). A stack is also found turned
/// around (the whole key reversed), as Biopython looks it up; together they cover all 96 stacks with one mismatch.
const MISMATCH_NN: [([u8; 5], f64, f64); 51] = [
    (*b"AG/TT", 1.0, 0.9),
    (*b"AT/TG", -2.5, -8.3),
    (*b"CG/GT", -4.1, -11.7),
    (*b"CT/GG", -2.8, -8.0),
    (*b"GG/CT", 3.3, 10.4),
    (*b"GG/TT", 5.8, 16.3),
    (*b"GT/CG", -4.4, -12.3),
    (*b"GT/TG", 4.1, 9.5),
    (*b"TG/AT", -0.1, -1.7),
    (*b"TG/GT", -1.4, -6.2),
    (*b"TT/AG", -1.3, -5.3),
    (*b"AA/TG", -0.6, -2.3),
    (*b"AG/TA", -0.7, -2.3),
    (*b"CA/GG", -0.7, -2.3),
    (*b"CG/GA", -4.0, -13.2),
    (*b"GA/CG", -0.6, -1.0),
    (*b"GG/CA", 0.5, 3.2),
    (*b"TA/AG", 0.7, 0.7),
    (*b"TG/AA", 3.0, 7.4),
    (*b"AC/TT", 0.7, 0.2),
    (*b"AT/TC", -1.2, -6.2),
    (*b"CC/GT", -0.8, -4.5),
    (*b"CT/GC", -1.5, -6.1),
    (*b"GC/CT", 2.3, 5.4),
    (*b"GT/CC", 5.2, 13.5),
    (*b"TC/AT", 1.2, 0.7),
    (*b"TT/AC", 1.0, 0.7),
    (*b"AA/TC", 2.3, 4.6),
    (*b"AC/TA", 5.3, 14.6),
    (*b"CA/GC", 1.9, 3.7),
    (*b"CC/GA", 0.6, -0.6),
    (*b"GA/CC", 5.2, 14.2),
    (*b"GC/CA", -0.7, -3.8),
    (*b"TA/AC", 3.4, 8.0),
    (*b"TC/AA", 7.6, 20.2),
    (*b"AA/TA", 1.2, 1.7),
    (*b"CA/GA", -0.9, -4.2),
    (*b"GA/CA", -2.9, -9.8),
    (*b"TA/AA", 4.7, 12.9),
    (*b"AC/TC", 0.0, -4.4),
    (*b"CC/GC", -1.5, -7.2),
    (*b"GC/CC", 3.6, 8.9),
    (*b"TC/AC", 6.1, 16.4),
    (*b"AG/TG", -3.1, -9.5),
    (*b"CG/GG", -4.9, -15.3),
    (*b"GG/CG", -6.0, -15.8),
    (*b"TG/AG", 1.6, 3.6),
    (*b"AT/TT", -2.7, -10.8),
    (*b"CT/GT", -5.0, -15.8),
    (*b"GT/CT", -2.2, -8.4),
    (*b"TT/AT", 0.2, -1.5),
];

fn mismatch_nn(top: &[u8], bottom: &[u8]) -> Option<(f64, f64)> {
    let key = [top[0], top[1], b'/', bottom[0], bottom[1]];
    let turned = [bottom[1], bottom[0], b'/', top[1], top[0]];
    MISMATCH_NN.iter().find(|(k, _, _)| *k == key || *k == turned).map(|&(_, h, s)| (h, s))
}

fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'G' => b'C',
        b'C' => b'G',
        x => x,
    }
}

/// Tm in °C of `primer` annealed at a site whose bases are `target`: the primer's sense, the same length, equal to the primer where
/// every base pairs (so a perfect site is `tm(primer)`). A mismatched pair inside the duplex takes its stacks from MISMATCH_NN.
/// `None`: lengths differ, a base is not A/C/G/T, or fewer than 2 pairs are left.
///
/// shortcut: mismatched pairs at the ends are left out of the duplex (no terminal-mismatch or dangling-end terms; a site found by
/// `binding_sites` starts and ends on a pair), and a stack of two mismatches in a row adds nothing (no data; Biopython skips it
/// too). Add those terms if sites with such ends or runs need a closer Tm.
pub fn site_tm(primer: &[u8], target: &[u8], c: &Conditions) -> Option<f64> {
    let (top, target) = (clean(primer), clean(target));
    if top.len() != target.len() || top.iter().chain(&target).any(|b| !b"ACGT".contains(b)) {
        return None;
    }
    if top == target {
        return tm(&top, c);
    }
    let pairs = |i: usize| top[i] == target[i];
    let (Some(first), Some(last)) = ((0..top.len()).find(|&i| pairs(i)), (0..top.len()).rposition(|i| pairs(i))) else {
        return None;
    };
    if last <= first {
        return None;
    }
    let (mut dh, mut ds) = (0.0, 0.0);
    for end in [top[first], top[last]] {
        let (h, e) = end_term(end);
        dh += h;
        ds += e;
    }
    for i in first..last {
        let (h, e) = match (pairs(i), pairs(i + 1)) {
            (true, true) => nn(&top[i..i + 2]),
            (false, false) => (0.0, 0.0),
            _ => mismatch_nn(&top[i..i + 2], &[complement(target[i]), complement(target[i + 1])]).unwrap_or((0.0, 0.0)),
        };
        dh += h;
        ds += e;
    }
    Some(melt(dh, ds, last + 1 - first, false, c))
}

/// The sites of `binding_sites(template, primer, circular, ...)` that melt at `min_tm` °C or above (`site_tm`: the primer's
/// annealed bases on the template's, mismatches counted). A Tm that cannot be computed (an ambiguity code) keeps the site;
/// `min_tm` <= 0 keeps every site. The one minimum-Tm rule of every primer site the app finds.
pub fn warm_sites(sites: Vec<Site>, template: &[u8], primer: &[u8], circular: bool, min_tm: f64, c: &Conditions) -> Vec<Site> {
    if min_tm <= 0.0 || template.is_empty() {
        return sites;
    }
    let n = template.len();
    sites
        .into_iter()
        .filter(|s| {
            let annealed = &primer[primer.len() - s.len.min(primer.len())..];
            let top: Vec<u8> = (s.start..s.start + s.len).map(|i| template[if circular { i % n } else { i.min(n - 1) }]).collect();
            let target = if s.forward { top } else { revcomp(&top) };
            site_tm(annealed, &target, c).is_none_or(|t| t >= min_tm)
        })
        .collect()
}

/// Lowest and highest Tm (°C) of a primer that may contain IUPAC codes: it is a mixture of every oligo the codes
/// stand for, and each one melts differently. `None` under the same conditions as `tm`, or for a letter that is
/// not a base or an IUPAC code. Not found by expanding the combinations (8 N's are 65,536 oligos): Tm = ΔH/(ΔS + c)
/// is optimised exactly over the choices with Dinkelbach's method, each round a small dynamic programme along the primer.
/// shortcut: the symmetric (self-complementary) correction is applied to the whole bracket when some expansion could be
/// self-complementary, so then the range may be a little wider than the true one; it always still contains every oligo.
pub fn tm_range(seq: &[u8], c: &Conditions) -> Option<(f64, f64)> {
    let s = clean(seq);
    let masks: Vec<u8> = s.iter().map(|&b| pattern_mask(b, true)).collect();
    if s.len() < 2 || masks.contains(&0) {
        return None;
    }
    let n = s.len();
    let salt = salt_ds(n, c);
    let ct = c.primer_nm * 1e-9;
    let options = |m: u8| -> Vec<u8> { b"ACGT".iter().enumerate().filter(|(i, _)| m >> i & 1 == 1).map(|(_, &b)| b).collect() };
    // Best (ΔH, ΔS) of any expansion for the objective ΔH·1000 − t·ΔS: smallest when `up`, else largest.
    let best = |t: f64, up: bool| -> (f64, f64) {
        let sign = if up { 1.0 } else { -1.0 };
        let cost = |(h, e): (f64, f64)| sign * (h * 1000.0 - t * e);
        let mut layer: Vec<(u8, f64, f64)> = options(masks[0]).into_iter().map(|b| (b, end_term(b).0, end_term(b).1)).collect();
        for &m in &masks[1..] {
            layer = options(m)
                .into_iter()
                .map(|b| {
                    let (_, h, e) = layer
                        .iter()
                        .map(|&(a, h, e)| { let (dh, de) = nn(&[a, b]); (a, h + dh, e + de) })
                        .min_by(|x, y| cost((x.1, x.2)).total_cmp(&cost((y.1, y.2))))
                        .unwrap();
                    (b, h, e)
                })
                .collect();
        }
        let (_, h, e) = layer
            .into_iter()
            .map(|(b, h, e)| (b, h + end_term(b).0, e + end_term(b).1))
            .min_by(|x, y| cost((x.1, x.2)).total_cmp(&cost((y.1, y.2))))
            .unwrap();
        (h, e)
    };
    let extreme = |up: bool, extra_ds: f64, k: f64| -> f64 {
        let tm_of = |(h, e): (f64, f64)| h * 1000.0 / (e + extra_ds + salt + R * k.ln());
        let mut t = 0.0;
        for _ in 0..100 {
            let next = tm_of(best(t, up));
            if (next - t).abs() < 1e-9 {
                break;
            }
            t = next;
        }
        t - 273.15
    };
    let (mut lo, mut hi) = (extreme(false, 0.0, ct / 4.0), extreme(true, 0.0, ct / 4.0));
    // Could some expansion be its own reverse complement? Then position i must allow the complement of position n-1-i.
    let comp = |m: u8| ((m & 1) << 3) | ((m & 2) << 1) | ((m & 4) >> 1) | ((m & 8) >> 3);
    if n % 2 == 0 && (0..n).all(|i| masks[i] & comp(masks[n - 1 - i]) != 0) {
        lo = lo.min(extreme(false, -1.4, ct));
        hi = hi.max(extreme(true, -1.4, ct));
    }
    Some((lo, hi))
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
/// Default number of bases at the 3' end that must match exactly: the polymerase starts there.
pub const EXACT_3P: usize = 3;

/// All binding sites of `primer` on either strand.
/// shortcut: brute-force scan, ~0.2 s on a 5 Mb genome; use an index (FM/k-mer) if
/// whole-genome primer lists become slow.
/// `exact_3p`: bases at the 3' end that must match exactly (at most the seed).
pub fn binding_sites(template: &[u8], primer: &[u8], circular: bool, seed: usize, max_mismatches: usize, exact_3p: usize) -> Vec<Site> {
    let p = clean(primer);
    // A degenerate base in the primer (R, N, ...) matches every base it stands for; the template is compared letter for letter.
    let same = |t: u8, q: u8| t == q || (b"ACGT".iter().position(|&x| x == t).map_or(0, |i| 1u8 << i) & pattern_mask(q, true)) != 0;
    let n = template.len();
    if p.len() < 8 || n == 0 || seed < 4 {
        return Vec::new(); // too short to bind specifically
    }
    let seed = p.len().min(seed).min(n);
    let exact_3p = exact_3p.min(seed);
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
        // Forward: mismatches allowed except in the last `exact_3p` bases of the seed.
        if let Some(mm) = seed_mismatches(&base, i as i64, fwd_seed, max_mismatches, &same, |k| k + exact_3p >= seed) {
            let (ext, ext_mm) = extend(reach, &same, |j| Some((base(i as i64 - 1 - j as i64)?, p[p.len() - seed - 1 - j])));
            let len = (seed + ext).min(n);
            let start = (i as i64 - ext as i64).rem_euclid(n as i64) as usize;
            sites.push(Site { start, len, forward: true, mismatches: mm + ext_mm });
        }
        // Reverse: the primer's 3' end is the first `exact_3p` bases of `rev_seed`.
        if let Some(mm) = seed_mismatches(&base, i as i64, rev_seed, max_mismatches, &same, |k| k < exact_3p) {
            let (ext, ext_mm) = extend(reach, &same, |j| Some((base((i + seed + j) as i64)?, rc[seed + j])));
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
fn extend(reach: usize, same: &impl Fn(u8, u8) -> bool, pair: impl Fn(usize) -> Option<(u8, u8)>) -> (usize, usize) {
    let (mut score, mut best, mut best_len, mut mm, mut best_mm) = (0i64, 0i64, 0, 0, 0);
    for j in 0..reach {
        let Some((t, q)) = pair(j) else { break };
        if same(t, q) {
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
    same: &impl Fn(u8, u8) -> bool,
    exact: impl Fn(usize) -> bool,
) -> Option<usize> {
    let mut mm = 0;
    for (k, &b) in seed.iter().enumerate() {
        if !same(base(at + k as i64)?, b) {
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
            let (h, e) = end_term(end);
            dh += h;
            ds += e;
        }
    }
    for pair in s.windows(2) {
        let (h, e) = nn(pair);
        dh += h;
        ds += e;
    }
    ds += salt_ds(s.len(), c);
    dh - 310.15 * ds / 1000.0
}

fn pairs(a: u8, b: u8) -> bool {
    matches!((a, b), (b'A', b'T') | (b'T', b'A') | (b'C', b'G') | (b'G', b'C'))
}

/// Uppercase, U -> T.
pub fn clean(seq: &[u8]) -> Vec<u8> {
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
            b'R' => b'Y',
            b'Y' => b'R',
            b'K' => b'M',
            b'M' => b'K',
            b'B' => b'V',
            b'V' => b'B',
            b'D' => b'H',
            b'H' => b'D',
            x => x, // S, W, N pair with themselves
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn warm_sites_keep_what_melts_high_enough() {
        // A primer and a copy with 2 mismatches, on either strand: the minimum Tm keeps or drops each by its own Tm.
        let c = Conditions::default();
        let p = b"AGCGGATAACAATTTCACACAGGA";
        let copy = b"AGCGGCTAACAATTTCACGCAGGA";
        let (perfect, lowered) = (tm(p, &c).unwrap(), site_tm(p, copy, &c).unwrap());
        let spacer = b"TTTTTTTTTT";
        let template = [&spacer[..], &p[..], &spacer[..], &revcomp(copy), &spacer[..]].concat();
        let found = |min: f64| -> Vec<(usize, bool)> {
            let sites = binding_sites(&template, p, false, SEED, 2, EXACT_3P);
            warm_sites(sites, &template, p, false, min, &c).iter().map(|s| (s.start, s.forward)).collect()
        };
        assert_eq!(found(0.0), vec![(10, true), (44, false)]);
        assert_eq!(found(lowered - 0.1), vec![(10, true), (44, false)]);
        assert_eq!(found(lowered + 0.1), vec![(10, true)]);
        assert!(found(perfect + 0.1).is_empty());
        // across the origin of a circle: the site's bases are read through it
        let circle = [&p[12..], &[b'T'; 20][..], &p[..12]].concat();
        let sites = binding_sites(&circle, p, true, SEED, 0, EXACT_3P);
        assert_eq!(sites.len(), 1);
        assert_eq!(warm_sites(sites, &circle, p, true, perfect - 0.1, &c).len(), 1);
        let sites = binding_sites(&circle, p, true, SEED, 0, EXACT_3P);
        assert!(warm_sites(sites, &circle, p, true, perfect + 0.1, &c).is_empty());
    }

    #[test]
    fn site_tm_with_mismatches() {
        // Against Biopython's Tm_NN (c_seq = the complement of the site, DNA_IMM1, saltcorr=5, dnac1 = dnac2 = primer/2).
        let c = Conditions::default();
        let p = b"AGCGGATAACAATTTCACACAGGA";
        assert_eq!(site_tm(p, p, &c), tm(p, &c));   // a perfect site is the primer's own Tm
        let one = site_tm(p, b"AGCGGATAACATTTTCACACAGGA", &c).unwrap();   // A·A at 12
        assert!((one - 56.58200646427804).abs() < 0.01, "{one}");
        let two = site_tm(p, b"AGCGGCTAACAATTTCACGCAGGA", &c).unwrap();   // A·G at 6, A·C at 19
        assert!((two - 52.154985731372165).abs() < 0.01, "{two}");
        assert!(two < one && one < tm(p, &c).unwrap());
        // ends that do not pair are left out: the same as the shorter duplex
        assert_eq!(site_tm(b"TAGCGGATAACAATTTCACACAGGA", b"GAGCGGATAACAATTTCACACAGGA", &c), tm(p, &c));
        assert!(site_tm(p, b"AGC", &c).is_none() && site_tm(b"ACGN", b"ACGT", &c).is_none());
    }

    #[test]
    fn tm_and_tm_range_of_plain_and_degenerate_primers() {
        // More exact values are checked against Biopython's Tm_NN (same model) in tests/test_core_sequence.py;
        // here a few exact values, input checks and tm_range.
        let c = Conditions::default();
        let t = tm(b"AGCGGATAACAATTTCACACAGGA", &c).unwrap();
        assert!(t > 55.0 && t < 65.0, "{t}");
        assert!(tm(b"ACGTN", &c).is_none());
        assert!(tm(b"A", &c).is_none());
        // More salt or more primer = higher Tm.
        let salty = Conditions { na_mm: 200.0, ..c };
        assert!(tm(b"AGCGGATAACAATTTCACACAGGA", &salty).unwrap() > t);

        // --- exact Tm (Biopython Tm_NN, DNA_NN3, saltcorr=5, Na 50, Mg 1.5, dNTPs 0.6) ---
        // plain primer: dnac1 = dnac2 = primer/2
        assert!((t - 60.80594739470723).abs() < 0.01, "{t}");
        let more = Conditions { primer_nm: 500.0, ..c };
        let t5 = tm(b"AGCGGATAACAATTTCACACAGGA", &more).unwrap();
        assert!((t5 - 63.596206495371575).abs() < 0.01, "{t5}");
        // self-complementary oligo: selfcomp=True, dnac1 = primer (C, not C/4, plus the symmetry term)
        let s = tm(b"GAATTCGAATTC", &c).unwrap();
        assert!((s - 30.94832964729642).abs() < 0.01, "{s}");

        // --- tm_range: a plain primer is its own Tm; bad letters give none (tm itself stays plain-ACGT only, above) ---
        let (lo, hi) = tm_range(b"AGCGGATAACAATTTCACACAGGA", &c).unwrap();
        assert!((lo - t).abs() < 1e-6 && (hi - t).abs() < 1e-6);
        assert!(tm_range(b"ACGTX", &c).is_none() && tm_range(b"A", &c).is_none());
        // --- tm_range of a degenerate primer brackets the Tm of every oligo it stands for ---
        for degenerate in ["AGCGGRTAACAATYTCACACNGGA", "ACGTSACGTWACGTKACGT", "GAATTCNNGAATTC", "RYGCGCRY"] {
            let tms: Vec<f64> = expand(degenerate).iter().map(|o| tm(o.as_bytes(), &c).unwrap()).collect();
            let (lo, hi) = tm_range(degenerate.as_bytes(), &c).unwrap();
            let (min, max) = (tms.iter().cloned().fold(f64::MAX, f64::min), tms.iter().cloned().fold(f64::MIN, f64::max));
            assert!(lo <= min + 1e-9 && hi >= max - 1e-9, "{degenerate}: [{lo}, {hi}] must bracket [{min}, {max}]");
            // tight unless the self-complementary widening applied
            if degenerate != "GAATTCNNGAATTC" && degenerate != "RYGCGCRY" {
                assert!((lo - min).abs() < 1e-6 && (hi - max).abs() < 1e-6, "{degenerate}: [{lo}, {hi}] vs [{min}, {max}]");
            }
        }
    }

    /// Every plain oligo a degenerate primer stands for (test helper; only for short ones).
    fn expand(seq: &str) -> Vec<String> {
        let mut out = vec![String::new()];
        for ch in seq.bytes() {
            let m = pattern_mask(ch, true);
            out = out
                .iter()
                .flat_map(|p| b"ACGT".iter().enumerate().filter(move |(i, _)| m >> i & 1 == 1).map(move |(_, &b)| format!("{p}{}", b as char)))
                .collect();
        }
        out
    }

    #[test]
    fn a_degenerate_primer_binds_every_site_it_matches_on_either_strand() {
        let flank = "T".repeat(20);
        let (a, g) = ("ACGTTGCAGGCTAGCTAGGA", "ACGTTGCAGGCTAGCTAGGA");
        // three degenerate positions: R (A/G), N, Y (C/T); the template carries three variants of the site
        let variant = |r: char, n: char, y: char| {
            let mut v: Vec<char> = a.chars().collect();
            (v[2], v[8], v[12]) = (r, n, y);
            v.into_iter().collect::<String>()
        };
        let template = format!("{flank}{}{flank}{}{flank}{}{flank}", variant('A', 'C', 'C'), variant('G', 'T', 'T'), variant('G', 'A', 'C'));
        let primer = {
            let mut v: Vec<char> = g.chars().collect();
            (v[2], v[8], v[12]) = ('R', 'N', 'Y');
            v.into_iter().collect::<String>()
        };
        let sites = binding_sites(template.as_bytes(), primer.as_bytes(), false, SEED, 0, EXACT_3P);
        let starts: Vec<usize> = sites.iter().filter(|s| s.forward).map(|s| s.start).collect();
        assert_eq!(starts, vec![20, 60, 100]);
        assert!(sites.iter().all(|s| s.len == 20 && s.mismatches == 0));
        // a base outside the code is still a mismatch (template has A where Y allows only C/T)
        let off = template.replacen(&variant('A', 'C', 'C'), &variant('A', 'C', 'A'), 1);
        let starts: Vec<usize> = binding_sites(off.as_bytes(), primer.as_bytes(), false, SEED, 0, EXACT_3P).iter().map(|s| s.start).collect();
        assert_eq!(starts, vec![60, 100]);
        // a degenerate base in the template still matches only the identical letter
        let n_template = template.replacen(&variant('A', 'C', 'C'), &variant('A', 'C', 'N'), 1);
        let starts: Vec<usize> = binding_sites(n_template.as_bytes(), primer.as_bytes(), false, SEED, 0, EXACT_3P).iter().map(|s| s.start).collect();
        assert_eq!(starts, vec![60, 100]);
        // --- a degenerate reverse primer binds the bottom strand through complemented codes ---
        let site2 = "ACGTTGCAGGCTAGCTAGGA";
        let rev_template = format!("{0}{site2}{0}", "T".repeat(20));
        // reverse primer = reverse complement of the site, with R in place of the base opposite the site's Y-compatible C/T
        let mut rc = revcomp(site2.as_bytes());
        rc[5] = b'R'; // faces a C on the top strand; R = A/G pairs with T/C
        let rev_sites = binding_sites(rev_template.as_bytes(), &rc, false, SEED, 0, EXACT_3P);
        assert_eq!(rev_sites.iter().filter(|s| !s.forward).map(|s| (s.start, s.len, s.mismatches)).collect::<Vec<_>>(), vec![(20, 20, 0)]);
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
        assert_eq!(binding_sites(template, &fwd, false, SEED, 0, EXACT_3P),
                   vec![Site { start: 10, len: 24, forward: true, mismatches: 1 }]);
        // Only 2 good bases after a mismatch: not worth it, they stay a tail.
        let mut short = template[10..34].to_vec();
        short[2] = b'C';
        assert_eq!(binding_sites(template, &short, false, SEED, 0, EXACT_3P)[0].start, 13);
        // Reverse strand: its seed is on the left (10..25), the extension goes right.
        let mut top = template[10..34].to_vec();
        top[19] = b'A'; // base 29
        let rev = revcomp(&top);
        assert_eq!(binding_sites(template, &rev, false, SEED, 0, EXACT_3P),
                   vec![Site { start: 10, len: 24, forward: false, mismatches: 1 }]);
        // Two wrong bases in a row stop it: they and what lies beyond are a tail.
        let mut two = template[10..34].to_vec();
        two[5] = b'A';
        two[6] = b'A';
        let sites = binding_sites(template, &two, false, SEED, 0, EXACT_3P);
        assert_eq!((sites[0].start, sites[0].len, sites[0].mismatches), (17, 17, 0));
    }

    #[test]
    fn finds_sites_on_both_strands_with_tails_and_mismatches() {
        //                 0         1         2         3         4
        //                 0123456789012345678901234567890123456789012345
        let template = b"TTTTTTTTTTACGTACGGATCCAGTTGACCATGGTTTTTTTTTTTT";
        // Forward primer = top strand 10..34 with a 6-base 5' tail.
        let fwd = b"GAATTCACGTACGGATCCAGTTGACCATGG";
        let sites = binding_sites(template, fwd, false, SEED, 0, EXACT_3P);
        assert_eq!(sites, vec![Site { start: 10, len: 24, forward: true, mismatches: 0 }]);
        // Reverse primer: reverse complement of 14..34.
        let rev = revcomp(&template[14..34]);
        let sites = binding_sites(template, &rev, false, SEED, 0, EXACT_3P);
        assert_eq!(sites, vec![Site { start: 14, len: 20, forward: false, mismatches: 0 }]);
        // One mismatch in the seed (not at the 3' end): only found when allowed.
        let mut one_off = fwd.to_vec();
        one_off[20] = b'A'; // was T
        assert!(binding_sites(template, &one_off, false, SEED, 0, EXACT_3P).is_empty());
        assert_eq!(binding_sites(template, &one_off, false, SEED, 1, EXACT_3P)[0].mismatches, 1);
        // A mismatch at the 3'-terminal base never binds.
        let mut bad_end = fwd.to_vec();
        *bad_end.last_mut().unwrap() = b'A';
        assert!(binding_sites(template, &bad_end, false, SEED, 2, EXACT_3P).is_empty());
        // A mismatch 5 bases from the 3' end: allowed with the default 3 exact bases, not when the last 5 must match exactly.
        let mut near_end = fwd.to_vec();
        near_end[25] = b'A'; // was C
        assert_eq!(binding_sites(template, &near_end, false, SEED, 1, EXACT_3P)[0].mismatches, 1);
        assert!(binding_sites(template, &near_end, false, SEED, 1, 5).is_empty());
    }

    #[test]
    fn site_across_the_origin() {
        let template = b"GGATCCAGTTTTTTTTTTTTTTTTTTACGTACGGATCC";
        let primer = b"ACGTACGGATCCGGATCCAG"; // last 8 bases continue at the start
        let sites = binding_sites(template, primer, true, SEED, 0, EXACT_3P);
        assert_eq!(sites, vec![Site { start: 26, len: 20, forward: true, mismatches: 0 }]);
    }

}
