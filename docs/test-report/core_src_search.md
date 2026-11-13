# `core/src/search.rs`

[← all test files](README.md) · source: [`core/src/search.rs`](../../core/src/search.rs)

## `finds_an_exact_motif_on_both_strands_in_position_order`

**Result:** ✅ passed · [line 232](../../core/src/search.rs#L232) · **Area:** Sequence search (DNA motif) · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** Sequence search could miss bottom-strand hits, report palindromes twice, ignore the top-strand-only option, or be case-sensitive.

Exact motif search finds palindromes once, finds a motif via its reverse complement on the bottom strand, honours top-strand-only, and ignores case.

| # | Input | Expected output |
|---|---|---|
| 1 | 'TTGAATTCAAAAGGATCCAAGTTT'; 'GAATTC' both strands | [(2, fwd, 0)] (palindrome once) |
| 2 | 'AAGG'; 'CCTT' both; 'CCTT' top only | [(10,fwd,0)]; [(10,rev,0)]; empty |
| 3 | 'ttgaattcaa', 'GAATTC' | 1 hit |

**Review notes:**

- Comment line 236 says 'GGATCC too' but GGATCC is never searched.

## `iupac_letters_in_the_query_match_every_base_they_stand_for`

**Result:** ✅ passed · [line 247](../../core/src/search.rs#L247) · **Area:** Sequence search (DNA motif) · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** IUPAC search queries could miss bases they stand for, or an N in the sequence could match everything, giving false or missing hits.

IUPAC codes in the query match all bases they stand for; with IUPAC off, R is invalid and N matches only N; an N in the sequence matches only a query position accepting every base; non-base letters give no result.

| # | Input | Expected output |
|---|---|---|
| 1 | 'AAGAATTCAAGAGTTCAAGACTTC', 'GARTTC', IUPAC on / off | [(2),(10)] forward / empty |
| 2 | 'ACGTACNT' 'GTAN' IUPAC on; 'ACGTACGT' 'GTAN' off; 'ACGTACNT' 'CNT' off | 1 hit; empty; 1 hit |
| 3 | 'ACGTNCGT' 'GTACGT' / 'GTNCGT' | empty / 1 hit |
| 4 | 'GAXTTC' | empty |

## `allows_a_number_of_mismatches_and_reports_how_many`

**Result:** ✅ passed · [line 266](../../core/src/search.rs#L266) · **Area:** Sequence search (DNA motif) · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** Search with mismatches allowed could miss near-matches or report a wrong mismatch count.

Fuzzy search returns hits with up to N mismatches and reports the count, on both strands.

| # | Input | Expected output |
|---|---|---|
| 1 | 'TTTGAATTCTTTGAGTTCTTTGAGCTCTTT', 'GAATTC', max 1 / max 2 | [(3,f,0),(12,f,1)] / plus (21,f,2) |
| 2 | 'GGCCTAGG', 'AAGG', max 1, both | [(2,rev,1),(4,fwd,1)] |

## `a_hit_can_run_across_the_origin_of_a_circle`

**Result:** ✅ passed · [line 278](../../core/src/search.rs#L278) · **Area:** Sequence search (DNA motif) · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** A motif crossing the origin of a plasmid could be missed by the search.

On a circular sequence a motif split across the origin is found (end > length), on either strand; not on a linear one.

| # | Input | Expected output |
|---|---|---|
| 1 | 'TTCAAAAAAAAAAAAGAA' (18 bp), 'GAATTC', linear / circular | empty / [MotifHit{15, 21, fwd, 0}] |
| 2 | 'TAGGAAAAAAAACT', 'CCTAAG' circular | 1 hit, reverse, start 12 |

## `stops_at_the_limit_and_handles_tiny_and_empty_inputs`

**Result:** ✅ passed · [line 291](../../core/src/search.rs#L291) · **Area:** Sequence search (DNA motif) · **Checked against:** hand-computed · **Strength:** 🟡 adequate

**Why it matters:** A search could run past its hit limit or crash on an empty or very short sequence or query.

The hit limit is honoured and empty or too-short inputs return nothing without panicking; a 1-base circle works.

| # | Input | Expected output |
|---|---|---|
| 1 | 'ACACACACACACACAC', 'AC', limit 3 | 3 hits |
| 2 | empty seq; empty motif; motif longer than seq | empty each |
| 3 | 'A' circular, motif 'A' | 1 hit |

**Review notes:**

- Limit test does not check that the first 3 hits (0, 2, 4) are the ones kept.

## `finds_a_protein_motif_in_all_six_frames`

**Result:** ✅ passed · [line 301](../../core/src/search.rs#L301) · **Area:** Sequence search (protein motif) · **Checked against:** hand-computed · **Strength:** 🟡 adequate

**Why it matters:** Protein motif search could miss hits in some frames, report wrong DNA coordinates, or mishandle X wildcards and stops.

Protein motif search translates frames on both strands, reports top-strand coordinates and frame, supports X wildcards, mismatches and '*' stops.

| # | Input | Expected output |
|---|---|---|
| 1 | 'CCCATGAAAGGGTTTCC', 'MKGF' | [ProteinHit{3, 15, fwd, frame 1, 0}] |
| 2 | revcomp of it | 1 hit, reverse, (2, 14) |
| 3 | 'KGF' in 'ATGAAAGGGTTT', 'CATG...', 'CCATG...' | frames 1, 2, 3 |
| 4 | 'MXGF'; 'MLGF' 0 mm; 'MLGF' 1 mm | 1 hit; empty; mismatches 1 |
| 5 | 'ATGTAAGGG' 'M*G'; 'AT' 'M' | 1 hit; empty |

**Review notes:**

- Name says 'all six frames' but only the three forward frames are checked; the reverse-strand frame number is never asserted (line 309).

**Fixed by the audit of November 2026:**

- 13 Nov 2026: gap closed. The proposed checks were added, pass on the current code and fail under the mutation that slipped through before. Was: Name says 'all six frames' but only the three forward frames are checked for their frame number; the one reverse case asserts strand and coordinates but not .frame, and only reverse frame 1 is ever exercised, so reverse frames 2 and 3 (and their coordinate mapping) are untested.
