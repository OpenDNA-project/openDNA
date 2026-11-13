# `core/src/primers.rs`

[← all test files](README.md) · source: [`core/src/primers.rs`](../../core/src/primers.rs)

## `tm_and_tm_range_of_plain_and_degenerate_primers`

**Result:** ✅ passed · [line 400](../../core/src/primers.rs#L400) · **Area:** Primer Tm · **Checked against:** property/invariant · **Strength:** 🔴 weak

**Why it matters:** Without it, the Rust Tm and the Tm range shown for degenerate primers could be implausible or not bracket the oligos they stand for.

Nearest-neighbour Tm is in a plausible range, rejects non-ACGT and too-short input, and rises with salt. For a primer with no degenerate letters tm_range returns (Tm, Tm); a non-IUPAC letter or a 1-base input gives None; tm() itself still refuses degenerate letters. tm_range of a degenerate primer brackets the Tm of every plain oligo it expands to, and equals the exact min and max when no expansion is self-complementary.

| # | Input | Expected output |
|---|---|---|
| 1 | 'AGCGGATAACAATTTCACACAGGA', default conditions (50 mM Na, 1.5 mM Mg, 0.6 mM dNTP, 50 nM) | 55 < Tm < 65 °C |
| 2 | 'ACGTN'; 'A' | None; None |
| 3 | same primer, Na 200 mM | higher Tm |
| 4 | tm_range('AGCGGATAACAATTTCACACAGGA'), default conditions | lo == hi == tm(same), within 1e-6 |
| 5 | tm_range('ACGTX'); tm_range('A') | None; None |
| 6 | tm('ACGTN') | None (tm stays plain-ACGT only) |
| 7 | 'AGCGGRTAACAATYTCACACNGGA' (16 oligos), default conditions | (lo, hi) == (min, max) of tm() over the 16 expansions, within 1e-6 |
| 8 | 'ACGTSACGTWACGTKACGT' (8 oligos) | (lo, hi) == exact (min, max), within 1e-6 |
| 9 | 'GAATTCNNGAATTC' (16 oligos, some self-complementary) | lo <= min and hi >= max (bracket only) |
| 10 | 'RYGCGCRY' (16 oligos, some self-complementary) | lo <= min and hi >= max (bracket only) |

**Review notes:**

- (was tm_is_sane) 10 °C window only; exact values are checked in Python by tests/test_core_sequence.py::test_melting_temperature_matches_santalucia_and_biopython. The code comment in the test still points to tests/test_primers.py, which no longer exists.
- (was tm_is_sane) Comment says 'more salt or more primer = higher Tm' but primer concentration is never varied; the self-complementary (symmetric) branch is untested.
- (was tm_is_sane) Superseded by tests/test_core_sequence.py::test_melting_temperature_matches_santalucia_and_biopython (exact Biopython values, the same input checks, salt/Mg/dNTP/primer direction); it only remains as a fast Rust smoke test.
- (was tm_range_of_a_plain_primer_is_its_tm_and_bad_letters_give_none) Duplicates the end of tests/test_core_sequence.py::test_melting_temperature_matches_santalucia_and_biopython (tm(degenerate) is None, tm_range(plain) == (tm, tm), tm_range('ACGTX') is None); the only extra is tm_range('A'), and tm_is_sane already checks tm('ACGTN') is None.
- (was tm_range_is_exactly_the_lowest_and_highest_tm_of_the_expanded_oligos) For the two self-complementary primers only bracketing is checked, so an arbitrarily wide range would pass there.
- (was tm_range_is_exactly_the_lowest_and_highest_tm_of_the_expanded_oligos) Overlaps tests/test_core_sequence.py::test_melting_temperature_matches_santalucia_and_biopython, which checks bracketing for 'AGCGGRTAACAATYTCACACNGGA' only; this Rust test is the stronger one (exact bounds, 4 primers).

**Fixed by the audit of November 2026:**

- 13 Nov 2026: gap closed. The proposed checks were added, pass on the current code and fail under the mutation that slipped through before. Was: tm() is only checked to lie in a 10 degree window (55..65) and to rise with salt; the comment's 'more primer = higher Tm' is never tested, and the self-complementary branch of tm() (C instead of C/4, -1.4 cal/K/mol) is never compared with anything: the bracket check compares tm_range with tm(), so a wrong tm() symmetric branch moves both sides together. \| The two self-complementary degenerate primers are checked only for bracketing (by design: the documented shortcut widens the range; on 'GAATTCNNGAATTC' lo = 32.17 vs true min 33.42).

## `a_degenerate_primer_binds_every_site_it_matches_on_either_strand`

**Result:** ✅ passed · [line 453](../../core/src/primers.rs#L453) · **Area:** Primer binding · **Checked against:** hand-computed · **Strength:** 🟡 adequate

**Why it matters:** Without it, a degenerate primer could miss sites it matches, or bind sites it does not, on either strand.

A 20-nt primer with R, N and Y binds three template variants each with zero mismatches; a template base outside the code (A against Y) is a mismatch; a degenerate N in the template matches only an identical primer letter. A reverse primer (reverse complement of a site) with one R in place of a G binds the opposite strand at the site with zero mismatches.

| # | Input | Expected output |
|---|---|---|
| 1 | template = T20 + site(A,C,C) + T20 + site(G,T,T) + T20 + site(G,A,C) + T20, where site = 'ACGTTGCAGGCTAGCTAGGA' with positions 2, 8, 12 set; primer = site with R, N, Y there; seed 15, 0 mismatches allowed | forward starts [20, 60, 100]; every site len 20, 0 mismatches |
| 2 | first site changed to (A,C,A): A at position 12 where the primer has Y | starts [60, 100] |
| 3 | first site changed to (A,C,N): template N against primer Y | starts [60, 100] |
| 4 | template = T20 + 'ACGTTGCAGGCTAGCTAGGA' + T20; primer = revcomp(site) = 'TCCTAGCTAGCCTGCAACGT' with base 5 (G) replaced by R; seed 15, 0 mismatches | reverse sites == [(start 20, len 20, 0 mismatches)] |

**Review notes:**

- (was a_degenerate_primer_binds_every_site_it_matches_with_no_mismatch) Near-duplicate of the degenerate-primer block of tests/test_primer_binding.py::test_binding_origin_ends_repeats_and_letters (three variants with R/Y/N, a base outside the code is a mismatch, N in the template matches only N), which also checks the reverse strand; one of the two is redundant.
- (was a_degenerate_primer_binds_every_site_it_matches_with_no_mismatch) The test binds the same string to two names (let (a, g) = (same, same)): harmless but confusing.
- (was a_degenerate_primer_binds_every_site_it_matches_with_no_mismatch) Degenerate position 2 lies outside the 15-nt 3' seed, so it only matters through the extension; the negative cases test position 12 (inside the seed) only.
- (was a_degenerate_reverse_primer_binds_the_opposite_strand_using_complemented_codes) One degenerate letter, one positive case and no negative case (an R facing a base it does not pair with), so code that let every IUPAC letter match anything would pass.
- (was a_degenerate_reverse_primer_binds_the_opposite_strand_using_complemented_codes) The code comment 'opposite a T on the top strand?' (with a question mark) is wrong: primer base 5 faces top-strand position 34, a C, which pairs with G (allowed by R).
- (was a_degenerate_reverse_primer_binds_the_opposite_strand_using_complemented_codes) Redundant with tests/test_primer_binding.py::test_binding_origin_ends_repeats_and_letters, which binds a three-letter degenerate primer to the reverse-complemented template (three sites, all reverse strand, 0 mismatches).

## `annealing_extends_across_a_lone_mismatch`

**Result:** ✅ passed · [line 491](../../core/src/primers.rs#L491) · **Area:** Primer binding · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** A primer's displayed binding site could stop at a lone mismatch (showing too much tail) or run on past two mismatches, on either strand.

Binding-site extension past the 15-nt 3' seed crosses a lone mismatch when enough matching bases follow, stops otherwise, and works on both strands.

| # | Input | Expected output |
|---|---|---|
| 1 | template 'TTTTTTTTTTACGTACGGATCCAGTTGACCATGGTTTTTTTTTTTT' (46 bp); fwd primer = t[10..34] with primer base 4 A->C; seed 15, 0 mm | [Site{start 10, len 24, fwd, mismatches 1}] |
| 2 | primer base 2 G->C | first site start 13 (2 good bases left as tail) |
| 3 | revcomp of t[10..34] with base 29 C->A | [Site{10, 24, reverse, 1 mismatch}] |
| 4 | primer bases 5 and 6 -> 'AA' | first site (start 17, len 17, 0 mismatches) |

## `finds_sites_on_both_strands_with_tails_and_mismatches`

**Result:** ✅ passed · [line 521](../../core/src/primers.rs#L521) · **Area:** Primer binding · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** Primers with 5' tails, reverse primers, or primers with mismatches could be shown binding wrongly, including a 3'-mismatched primer shown as binding.

A primer with a 5' tail binds only over its annealing part; reverse primers are found; seed mismatches count only when allowed; a 3'-terminal mismatch never binds.

| # | Input | Expected output |
|---|---|---|
| 1 | same 46-bp template; 'GAATTC' + t[10..34] | [Site{10, 24, fwd, 0}] |
| 2 | revcomp(t[14..34]) | [Site{14, 20, reverse, 0}] |
| 3 | tailed primer with base 20 T->A; max 0 / max 1 | no site / first site mismatches 1 |
| 4 | last base -> 'A', max 2 | no site |

## `site_across_the_origin`

**Result:** ✅ passed · [line 545](../../core/src/primers.rs#L545) · **Area:** Primer binding · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** A primer binding across the origin of a plasmid could be missed.

On a circular template a primer whose 3' end continues past the origin is found as one site.

| # | Input | Expected output |
|---|---|---|
| 1 | template 'GGATCCAG' + 18 T + 'ACGTACGGATCC' (38 bp), circular; primer 'ACGTACGGATCCGGATCCAG' | [Site{start 26, len 20, fwd, 0}] |

**Review notes:**

- No linear counterpart (that the site is not found / is shorter on a linear template).
