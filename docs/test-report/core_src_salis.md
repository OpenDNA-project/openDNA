# `core/src/salis.rs`

[← all test files](README.md) · source: [`core/src/salis.rs`](../../core/src/salis.rs)

## `a_direct_repeat_is_found_once_at_its_full_length_with_both_starts`

**Result:** ✅ passed · [line 177](../../core/src/salis.rs#L177) · **Area:** Synthesis Success Calculator (Salis) · **Checked against:** hand-computed · **Strength:** 🟡 adequate

**Why it matters:** The repeat finder could report a copy at the wrong place or twice.

A 12-mer placed twice with a 20-base spacer is one direct repeat of 12 with starts 0 and 32.

| # | Input | Expected output |
|---|---|---|
| 1 | GATTACAGGCTA + 20 spacer + GATTACAGGCTA | (12, direct) at [0, 32]; every repeat is 9 or longer |

**Review notes:**

- The spacer has repeats of its own (9 bases), so only the 12-mer is asserted exactly.

## `a_sequence_without_repeats_has_none`

**Result:** ✅ passed · [line 190](../../core/src/salis.rs#L190) · **Area:** Synthesis Success Calculator (Salis) · **Checked against:** hand-computed · **Strength:** 🟡 adequate

**Why it matters:** A repeat-free sequence must not be flagged.

Two sequences with no 9-base repeat give an empty list; one shorter than 9 as well.

| # | Input | Expected output |
|---|---|---|
| 1 | ACGTTGCAAGCTTAGGCCATAT | no repeat |
| 2 | ACG | no repeat |

## `a_gc_rich_stem_scores_one_and_a_poly_a_run_scores_zero`

**Result:** ✅ passed · [line 196](../../core/src/salis.rs#L196) · **Area:** Synthesis Success Calculator (Salis) · **Checked against:** hand-computed · **Strength:** 🟡 adequate

**Why it matters:** The rule-based hairpin pass feeds richest_hairpin.

A GC-only 11 bp stem with a 3-base loop scores 1.0; 60 A score 0.0.

| # | Input | Expected output |
|---|---|---|
| 1 | GCGCGCGCGCGAAACGCGCGCGCGC | 1.0 |
| 2 | 60 x A | 0.0 |
