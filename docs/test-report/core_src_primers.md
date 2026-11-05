# `core/src/primers.rs`

[← all test files](README.md) · source: [`core/src/primers.rs`](../../core/src/primers.rs)

## `tm_is_sane`

**Result:** ✅ passed · [line 322](../../core/src/primers.rs#L322) · **Area:** Primer Tm · **Checked against:** property/invariant · **Strength:** 🔴 weak

Nearest-neighbour Tm is in a plausible range, rejects non-ACGT and too-short input, and rises with salt.

| # | Input | Expected output |
|---|---|---|
| 1 | 'AGCGGATAACAATTTCACACAGGA', default conditions (50 mM Na, 1.5 mM Mg, 0.6 mM dNTP, 50 nM) | 55 < Tm < 65 °C |
| 2 | 'ACGTN'; 'A' | None; None |
| 3 | same primer, Na 200 mM | higher Tm |

**Review notes:**

- 10 °C window only; exact values are deferred to the Python test (tests/test_primers.py, per comment line 323).
- Comment line 330 says 'more salt or more primer = higher Tm' but primer concentration is never varied; the self-complementary (symmetric) branch is untested.

## `annealing_extends_across_a_lone_mismatch`

**Result:** ✅ passed · [line 336](../../core/src/primers.rs#L336) · **Area:** Primer binding · **Checked against:** hand-computed · **Strength:** 🟢 strong

Binding-site extension past the 15-nt 3' seed crosses a lone mismatch when enough matching bases follow, stops otherwise, and works on both strands.

| # | Input | Expected output |
|---|---|---|
| 1 | template 'TTTTTTTTTTACGTACGGATCCAGTTGACCATGGTTTTTTTTTTTT' (46 bp); fwd primer = t[10..34] with primer base 4 A->C; seed 15, 0 mm | [Site{start 10, len 24, fwd, mismatches 1}] |
| 2 | primer base 2 G->C | first site start 13 (2 good bases left as tail) |
| 3 | revcomp of t[10..34] with base 29 C->A | [Site{10, 24, reverse, 1 mismatch}] |
| 4 | primer bases 5 and 6 -> 'AA' | first site (start 17, len 17, 0 mismatches) |

## `finds_sites_on_both_strands_with_tails_and_mismatches`

**Result:** ✅ passed · [line 366](../../core/src/primers.rs#L366) · **Area:** Primer binding · **Checked against:** hand-computed · **Strength:** 🟢 strong

A primer with a 5' tail binds only over its annealing part; reverse primers are found; seed mismatches count only when allowed; a 3'-terminal mismatch never binds.

| # | Input | Expected output |
|---|---|---|
| 1 | same 46-bp template; 'GAATTC' + t[10..34] | [Site{10, 24, fwd, 0}] |
| 2 | revcomp(t[14..34]) | [Site{14, 20, reverse, 0}] |
| 3 | tailed primer with base 20 T->A; max 0 / max 1 | no site / first site mismatches 1 |
| 4 | last base -> 'A', max 2 | no site |

## `site_across_the_origin`

**Result:** ✅ passed · [line 390](../../core/src/primers.rs#L390) · **Area:** Primer binding · **Checked against:** hand-computed · **Strength:** 🟢 strong

On a circular template a primer whose 3' end continues past the origin is found as one site.

| # | Input | Expected output |
|---|---|---|
| 1 | template 'GGATCCAG' + 18 T + 'ACGTACGGATCC' (38 bp), circular; primer 'ACGTACGGATCCGGATCCAG' | [Site{start 26, len 20, fwd, 0}] |

**Review notes:**

- No linear counterpart (that the site is not found / is shorter on a linear template).

## `dimers_and_hairpins`

**Result:** ✅ passed · [line 398](../../core/src/primers.rs#L398) · **Area:** Primer dimers and hairpins · **Checked against:** property/invariant · **Strength:** 🔴 weak

Primer-dimer and hairpin ΔG are strongly negative for palindromic 3' ends, full complements and GC stems, and zero for poly-A.

| # | Input | Expected output |
|---|---|---|
| 1 | self-dimer 'ACGTTTAGCGAATTCGC' | 3' dimer ΔG < -5 and overall <= 3' value |
| 2 | poly-A 16 self-dimer | (0.0, 0.0) |
| 3 | 'AGCGGATAACAATTTCACACAGGA' vs its revcomp | ΔG < -20 kcal/mol |
| 4 | hairpin 'GCGCGCTTTTGCGCGC'; poly-A 15 | ΔG < -3 and 3' value == overall; (0.0, 0.0) |

**Review notes:**

- Only loose thresholds (< -5, < -20, < -3): wrong NN parameters, salt correction or loop penalties would still pass; no exact ΔG is checked against a reference (e.g. Primer3).
