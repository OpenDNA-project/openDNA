# `core/src/cas9.rs`

[← all test files](README.md) · source: [`core/src/cas9.rs`](../../core/src/cas9.rs)

## `finds_a_site_with_its_pam_on_the_top_strand`

**Result:** ✅ passed · [line 177](../../core/src/cas9.rs#L177) · **Area:** CRISPR guide search · **Checked against:** hand-computed · **Strength:** 🟢 strong

A perfect 20-nt protospacer followed by an NGG PAM is found at the right coordinates on the top strand.

| # | Input | Expected output |
|---|---|---|
| 1 | guide 'ACGTTGCAAGTCCATGCATG'; template 'TTTTT' + guide + 'TGG' + 'CCCCC'; PAM NGG, 0 mismatches, linear | exactly 1 site: start 5, end 28, forward, 0 mismatches, mask 0 |

## `finds_a_site_on_the_bottom_strand`

**Result:** ✅ passed · [line 188](../../core/src/cas9.rs#L188) · **Area:** CRISPR guide search · **Checked against:** hand-computed · **Strength:** 🟢 strong

A site written on the bottom strand is found and reported in top-strand coordinates with the PAM complement on its left.

| # | Input | Expected output |
|---|---|---|
| 1 | template 'TTTTT' + revcomp(guide + 'AGG') + 'CCCCC'; NGG, 0 mm | exactly 1 site: (5, 28, forward=false); template[5..8] = 'CCT' |

## `needs_the_pam`

**Result:** ✅ passed · [line 201](../../core/src/cas9.rs#L201) · **Area:** CRISPR guide search · **Checked against:** hand-computed · **Strength:** 🟡 adequate

Without a fitting PAM there is no site even with mismatches allowed; other IUPAC PAMs (NNGRRT) are honoured.

| # | Input | Expected output |
|---|---|---|
| 1 | 'TTTTT' + guide + 'TAACCCCC', NGG, 3 mismatches allowed | no site |
| 2 | 'TTTTT' + guide + 'ATGAAT', PAM 'NNGRRT', 0 mm | 1 site |

**Review notes:**

- The NNGRRT case checks only the count, not position/strand.

## `counts_mismatches_and_records_where`

**Result:** ✅ passed · [line 211](../../core/src/cas9.rs#L211) · **Area:** CRISPR guide search · **Checked against:** hand-computed · **Strength:** 🟢 strong

Mismatches are counted against the limit and their guide positions are recorded in the bit mask.

| # | Input | Expected output |
|---|---|---|
| 1 | guide with base 0 A->T and base 19 G->C in template, PAM 'CGG'; max 1 mismatch | no site |
| 2 | same, max 2 | 1 site, mismatches 2, mask = 1 \| (1<<19) |

## `a_site_can_run_across_the_origin_of_a_circle`

**Result:** ✅ passed · [line 226](../../core/src/cas9.rs#L226) · **Area:** CRISPR guide search · **Checked against:** hand-computed · **Strength:** 🟡 adequate

On a circular template a site split across the origin is found once, with end past the length; not on a linear template.

| # | Input | Expected output |
|---|---|---|
| 1 | site = guide + 'TGG' (23 bp) split: last 5 bases at start, 20 C, first 18 bases at end (43 bp); circular | 1 site, start = 43-18 = 25, end > 43 |
| 2 | same, linear | no site |

**Review notes:**

- End is only checked as > length (line 235); exact end 48 and strand are not asserted.

## `several_guides_and_lower_case`

**Result:** ✅ passed · [line 240](../../core/src/cas9.rs#L240) · **Area:** CRISPR guide search · **Checked against:** hand-computed · **Strength:** 🟡 adequate

Several guides are searched in one call (one result list each) and a lower-case template matches.

| # | Input | Expected output |
|---|---|---|
| 1 | template 'tttttt' + lower(guide a) + 'TGGaaaaaa' + guide b + 'CGGtttttt'; guides [a, b], NGG, 0 mm | 1 site each; site of a starts before site of b |

**Review notes:**

- Only relative order is checked (line 246), not exact starts 6 and 35; lower-case guides (as opposed to template) are not tested.

## `finds_a_site_whose_pam_comes_first`

**Result:** ✅ passed · [line 250](../../core/src/cas9.rs#L250) · **Area:** CRISPR guide search (Cas12a) · **Checked against:** hand-computed · **Strength:** 🟢 strong

With pam_first (Cas12a, TTTV) the PAM is matched 5' of the protospacer, on both strands, and is not matched in Cas9 orientation.

| # | Input | Expected output |
|---|---|---|
| 1 | 'CCCCCTTTA' + guide + 'CCCCC', PAM TTTV, pam_first=true | 1 site (5, 29, fwd, 0 mm) |
| 2 | same text, pam_first=false | no site |
| 3 | 'CCCCC' + revcomp('TTTC' + guide) + 'CCCCC', pam_first=true | 1 site (5, 29, forward=false) |

## `mismatch_mask_counts_from_the_guides_5_prime_end_next_to_a_first_pam`

**Result:** ✅ passed · [line 269](../../core/src/cas9.rs#L269) · **Area:** CRISPR guide search (Cas12a) · **Checked against:** hand-computed · **Strength:** 🟡 adequate

For PAM-first nucleases bit 0 of the mismatch mask is the guide base next to the PAM.

| # | Input | Expected output |
|---|---|---|
| 1 | 'CCCCCTTTG' + guide with base 0 A->T + 'CCCCC', TTTV, 1 mm allowed, pam_first | first site: mismatches 1, mask 1 |

**Review notes:**

- Does not assert there is exactly one site or its position (line 275 indexes hits[0][0] directly).
