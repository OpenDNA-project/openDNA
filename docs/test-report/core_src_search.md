# `core/src/search.rs`

[← all test files](README.md) · source: [`core/src/search.rs`](../../core/src/search.rs)

## `finds_an_exact_motif_on_both_strands_in_position_order`

**Result:** ✅ passed · [line 347](../../core/src/search.rs#L347) · **Area:** Sequence search (DNA motif) · **Checked against:** hand-computed · **Strength:** 🟢 strong

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

**Result:** ✅ passed · [line 362](../../core/src/search.rs#L362) · **Area:** Sequence search (DNA motif) · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** IUPAC search queries could miss bases they stand for, or an N in the sequence could match everything, giving false or missing hits.

IUPAC codes in the query match all bases they stand for; with IUPAC off, R is invalid and N matches only N; an N in the sequence matches only a query position accepting every base; non-base letters give no result.

| # | Input | Expected output |
|---|---|---|
| 1 | 'AAGAATTCAAGAGTTCAAGACTTC', 'GARTTC', IUPAC on / off | [(2),(10)] forward / empty |
| 2 | 'ACGTACNT' 'GTAN' IUPAC on; 'ACGTACGT' 'GTAN' off; 'ACGTACNT' 'CNT' off | 1 hit; empty; 1 hit |
| 3 | 'ACGTNCGT' 'GTACGT' / 'GTNCGT' | empty / 1 hit |
| 4 | 'GAXTTC' | empty |

## `an_ambiguity_code_in_the_sequence_fits_a_query_that_accepts_all_it_stands_for`

**Result:** ✅ passed · [line 381](../../core/src/search.rs#L381) · **Area:** Motif search · **Checked against:** hand-computed · **Strength:** 🟡 adequate

**Why it matters:** A hit on a sequence holding R, Y or U could be lost without a word, so a search or a digest of a draft assembly reports fewer places than exist.

A sequence ambiguity code fits a query position that accepts every base it stands for (R fits R and N, not A, which costs one mismatch); U counts as T; lower case works; the reverse strand and IUPAC-off queries behave; check_motif names a letter that cannot be searched.

| # | Input | Expected output |
|---|---|---|
| 1 | sequence AAGARTTCAA, queries GARTTC / GANTTC / GAATTC (IUPAC on) | hit at 2; hit; none (1 mismatch allowed: hit with 1 mismatch) |
| 2 | 'ggaauucc', GAATTC | 1 hit (U is T, lower case) |
| 3 | AAGARTTCAA, GAGTTC with IUPAC off | none |
| 4 | 'TTCCYGGTT', CCRGG, both strands | one hit at 2 on the reverse strand |
| 5 | check_motif on GARTTC (on and off), GAXTTC, GANTTC (off), empty | ok; error naming 'R' when off; error naming 'X'; ok; ok |

**Review notes:**

- find_motifs (the enzyme pass) gets the same table; many_motifs_in_one_pass_match_find_motif_one_by_one compares it with find_motif over a sequence holding N, R, y and U.

## `allows_a_number_of_mismatches_and_reports_how_many`

**Result:** ✅ passed · [line 401](../../core/src/search.rs#L401) · **Area:** Sequence search (DNA motif) · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** Search with mismatches allowed could miss near-matches or report a wrong mismatch count.

Fuzzy search returns hits with up to N mismatches and reports the count, on both strands.

| # | Input | Expected output |
|---|---|---|
| 1 | 'TTTGAATTCTTTGAGTTCTTTGAGCTCTTT', 'GAATTC', max 1 / max 2 | [(3,f,0),(12,f,1)] / plus (21,f,2) |
| 2 | 'GGCCTAGG', 'AAGG', max 1, both | [(2,rev,1),(4,fwd,1)] |

## `a_hit_can_run_across_the_origin_of_a_circle`

**Result:** ✅ passed · [line 413](../../core/src/search.rs#L413) · **Area:** Sequence search (DNA motif) · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** A motif crossing the origin of a plasmid could be missed by the search.

On a circular sequence a motif split across the origin is found (end > length), on either strand; not on a linear one.

| # | Input | Expected output |
|---|---|---|
| 1 | 'TTCAAAAAAAAAAAAGAA' (18 bp), 'GAATTC', linear / circular | empty / [MotifHit{15, 21, fwd, 0}] |
| 2 | 'TAGGAAAAAAAACT', 'CCTAAG' circular | 1 hit, reverse, start 12 |

## `many_motifs_in_one_pass_match_find_motif_one_by_one`

**Result:** ✅ passed · [line 426](../../core/src/search.rs#L426) · **Area:** Sequence search (DNA motif) · **Checked against:** independent implementation · **Strength:** 🟢 strong

**Why it matters:** Restriction sites (map, sequence, Enzymes table, digests, cloning) could be missed or invented by the one-pass search that replaced the Python one.

find_motifs (all the enzymes' sites in one pass, used by every restriction-site display and digest) gives, for each motif, exactly what find_motif gives with no mismatches and no limit, on a linear and a circular sequence.

| # | Input | Expected output |
|---|---|---|
| 1 | 3000 made-up bases (A C G T, lower case, N, X); GAATTC, GGTCTC, GCNNNNNNNGC, RGATCY, AC, CCWGG, GAAGAC, NNN, GAXTTC; linear and circular | each motif's hits == find_motif(seq, motif, 0, iupac, both strands, circular, no limit) |
| 2 | 'TTCAAAGAA' circular, GAATTC | 1 hit (across the origin) |
| 3 | 'AATTCG' circular, CGAATTC (longer than the circle) | == find_motif's answer, not empty (read once round) |
| 4 | '' linear, GAATTC | none |

**Review notes:**

- find_motif is the reference, so a bug both share (the IUPAC masks of pattern_mask) would not be seen; those are checked by iupac_letters_in_the_query_match_every_base_they_stand_for.

## `cuts_mirror_on_the_bottom_strand_wrap_on_a_circle_and_stop_at_linear_ends`

**Result:** ✅ passed · [line 450](../../core/src/search.rs#L450) · **Area:** Restriction enzymes · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** Enzyme cut positions (map, sequence, digests, cloning) could be wrong for Type IIS enzymes on the bottom strand, through the origin or at linear ends.

find_cuts turns site hits into (top, bottom) cut positions: as given on the top strand, mirrored for a site on the bottom strand (an enzyme that cuts outside its site), dropped past a linear end, taken round the origin on a circle.

| # | Input | Expected output |
|---|---|---|
| 1 | 'AAGAATTCAA' linear, EcoRI GAATTC (1, 5) | [(3, 7)] |
| 2 | 'GGTCTCAAAAAAAAAA' linear, BsaI GGTCTC (7, 11) | [(7, 11)] |
| 3 | 'AAAAAAAAAAGAGACCAAAA' linear, BsaI (site on the bottom strand at 10) | [(5, 9)] |
| 4 | 'GGTCTCAA' linear / circular, BsaI and EcoRI | [] / [(7, 3)] and [] |

## `stops_at_the_limit_and_handles_tiny_and_empty_inputs`

**Result:** ✅ passed · [line 465](../../core/src/search.rs#L465) · **Area:** Sequence search (DNA motif) · **Checked against:** hand-computed · **Strength:** 🟡 adequate

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

**Result:** ✅ passed · [line 475](../../core/src/search.rs#L475) · **Area:** Sequence search (protein motif) · **Checked against:** hand-computed · **Strength:** 🟡 adequate

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
