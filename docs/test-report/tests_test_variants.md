# `tests/test_variants.py`

[← all test files](README.md) · source: [`tests/test_variants.py`](../../tests/test_variants.py)

## `test_change_text_and_which_variants_can_be_applied`

**Result:** ✅ passed · [line 43](../../tests/test_variants.py#L43) · **Area:** Variants from reads · **Checked against:** hand-computed · **Strength:** 🟢 strong

The label of each kind of variant, and which ones can be written into the sequence (an insertion whose bases the reads did not settle cannot).

| # | Input | Expected output |
|---|---|---|
| 1 | snp pos 4 A>G; del 'ACGTA'; del 11 x A; ins 'TT'; ins '' | 'A > G', 'deletion of 5 bp (ACGTA)', 'deletion of 11 bp', 'insertion of TT', 'insertion of ?' |
| 2 | can_apply on snp, del, ins 'TT', ins '', ins 'TNT' | True, True, True, False, False |

**Fixed by the audit of November 2026:**

- New test added by the audit: app/variants.py was 13 % covered and the GUI path to it is skipped without bwa/samtools.

## `test_each_kind_of_change_equals_the_string_edit`

**Result:** ✅ passed in 0.01 s · [line 55](../../tests/test_variants.py#L55) · **Area:** Variants from reads · **Checked against:** in-test reference model · **Strength:** 🟢 strong

Applying one variant changes the bases exactly like the same edit on a plain string, as one undo step; after an indel the user is told to map the reads again and the variant list is emptied, after a SNP it is rescanned.

| # | Input | Expected output |
|---|---|---|
| 1 | 32-base SEQ; snp at 5 C>T; at 0 A>G | SEQ with that base replaced; entry.variants None (rescan); no re-map message |
| 2 | del 'AAGG' at 8; del 'T' at 31 (last base) | SEQ[:8]+SEQ[12:]; SEQ[:31]; entry.variants []; 'map the reads again' message |
| 3 | ins 'CCC' at 8; ins 'GG' at 32 (after the end) | SEQ[:8]+'CCC'+SEQ[8:]; SEQ+'GG'; re-map message |

**Fixed by the audit of November 2026:**

- New test added by the audit: app/variants.py was 13 % covered and the GUI path to it is skipped without bwa/samtools.

## `test_several_changes_are_applied_last_position_first_and_each_undoes`

**Result:** ✅ passed in 0.00 s · [line 77](../../tests/test_variants.py#L77) · **Area:** Variants from reads · **Checked against:** in-test reference model · **Strength:** 🟢 strong

Several variants are written right to left so each position still refers to the original bases; the result does not depend on the order they were chosen in, and each change is its own undo step.

| # | Input | Expected output |
|---|---|---|
| 1 | ins TTT at 20, snp G>C at 2, del GGCC at 10, ins NN at 25; in this order and reversed | 3 applied (NN skipped); bases = SEQ[:2]+'C'+SEQ[3:10]+SEQ[14:20]+'TTT'+SEQ[20:] |
| 2 | undo until the original comes back | exactly 3 undo steps, then 3 redo steps available |

**Fixed by the audit of November 2026:**

- New test added by the audit: app/variants.py was 13 % covered and the GUI path to it is skipped without bwa/samtools.

## `test_nothing_is_applied_while_the_dna_is_locked_or_nothing_can_be`

**Result:** ✅ passed in 0.00 s · [line 95](../../tests/test_variants.py#L95) · **Area:** Variants from reads · **Checked against:** hand-computed · **Strength:** 🟢 strong

The DNA lock refuses every change and leaves the read entry untouched; an empty or unusable choice changes nothing.

| # | Input | Expected output |
|---|---|---|
| 1 | dna_locked = True; snp at 5 | 0 applied; bases and undo history unchanged; entry not reset; a 'locked' message |
| 2 | unlocked; [] and [ins ''] | 0 applied; nothing changes |

**Fixed by the audit of November 2026:**

- New test added by the audit: app/variants.py was 13 % covered and the GUI path to it is skipped without bwa/samtools.

## `test_variation_features_cover_the_right_bases`

**Result:** ✅ passed · [line 110](../../tests/test_variants.py#L110) · **Area:** Variants from reads · **Checked against:** hand-computed · **Strength:** 🟢 strong

The GenBank variation feature written for a variant covers the right bases (1-based), with the sample bases in /replace and the read counts in /note.

| # | Input | Expected output |
|---|---|---|
| 1 | snp at 0-based 5, C>T, 7 of 10 reads, sample.bam | ('variation', '6', label 'C > T', replace 'T', note '7 of 10 reads (70 %) in sample.bam') |
| 2 | del AAGG at 8; ins CCC at 8; ins at 0; ins at 32 (length 32) | '9..12' with replace ''; '8..9'; '1..1'; '32..32' |
| 3 | count 0, depth 0 | note '0 of 0 reads (0 %)' (no division by zero) |

**Fixed by the audit of November 2026:**

- New test added by the audit: app/variants.py was 13 % covered and the GUI path to it is skipped without bwa/samtools.

## `test_annotating_adds_one_feature_per_variant_in_one_undo_step_and_leaves_the_bases`

**Result:** ✅ passed in 0.00 s · [line 124](../../tests/test_variants.py#L124) · **Area:** Variants from reads · **Checked against:** hand-computed · **Strength:** 🟢 strong

Annotating adds one variation feature per variant on the reference bases of that variant, as one undo step, without changing the sequence.

| # | Input | Expected output |
|---|---|---|
| 1 | snp at 5 and del AAGG at 8 on SEQ | 2 features at '6' and '9..12' covering 'C' and 'AAGG'; bases unchanged |
| 2 | one undo; then annotate [] | both features gone; 0 added |

**Fixed by the audit of November 2026:**

- New test added by the audit: app/variants.py was 13 % covered and the GUI path to it is skipped without bwa/samtools.
