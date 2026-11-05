# `tests/test_primer_binding.py`

[← all test files](README.md) · source: [`tests/test_primer_binding.py`](../../tests/test_primer_binding.py)

## `test_the_rust_search_agrees_with_the_rules_written_in_python`

**Result:** ✅ passed in 0.14 s · [line 153](../../tests/test_primer_binding.py#L153) · **Area:** Primer binding vs reference model · **Checked against:** independent implementation · **Strength:** 🟢 strong

The Rust primer search produces exactly the same site list (same order) as an independent Python reimplementation of the documented binding rules, over many randomized templates.

| # | Input | Expected output |
|---|---|---|
| 1 | 1500 random cases (seed 5): templates 12-90 bp with inserted 2-5 bp repeats, primers 6-30 bp cut out with ~6% noise, optional revcomp and random 5' tail, circular ~50%, limit 0-3, seed in {4,8,12,15,20} | find_primers output == reference_sites output for every case |
| 2 | total sites found | > 500 (cases really contain sites) |

## `test_document_binding_sites_and_find_primers_agree`

**Result:** ✅ passed in 0.00 s · [line 177](../../tests/test_primer_binding.py#L177) · **Area:** Primer binding API agreement · **Checked against:** round-trip/self-consistency · **Strength:** 🟡 adequate

Document.binding_sites and the module-level find_primers return identical sites for the same primer, so the two entry points cannot drift apart.

| # | Input | Expected output |
|---|---|---|
| 1 | 300 bp template, linear and circular; primers primer/revcomp/one-mismatch/tailed; limits 0,1,2; default and seed=12 | doc.binding_sites(p, limit[, seed]) == sites(template, p, ...) for all combinations |

**Review notes:**

- both sides ultimately call the same Rust search (sites() wraps find_primers), so this checks the two Rust entry points agree, not correctness against an independent oracle; no docstring

## `test_single_mismatch_positions_on_both_strands`

**Result:** ✅ passed in 0.00 s · [line 203](../../tests/test_primer_binding.py#L203) · **Area:** Primer binding: single mismatch · **Checked against:** hand-computed · **Strength:** 🟢 strong

A primer with a single changed base is reported at the right coordinates and count on both strands, with a mismatch becoming a 5' tail when too near the 5' end and the site vanishing when within the last 3 bases.

| # | Input | Expected output |
|---|---|---|
| 1 | lengths 16/20/24/30; mismatch at each position k from the 5' end; both strands | site == expected_single_mismatch(k, length, 100, 100+length, strand) (seed/extension/tail/no-bind rules) |
| 2 | 24 bp core, mismatch at k in {4,9,15,20}; forward primer and its reverse-complement partner | both give coords (100,124) with count 1; strands (1, -1) |

## `test_the_limit_counts_mismatches_and_decides_what_binds`

**Result:** ✅ passed in 0.00 s · [line 228](../../tests/test_primer_binding.py#L228) · **Area:** Primer binding: mismatch limit · **Checked against:** hand-computed · **Strength:** 🟢 strong

Mismatches are counted and the max_mismatches limit decides whether a site binds, including exactly-at-limit seed cases on both strands.

| # | Input | Expected output |
|---|---|---|
| 1 | 24 bp core, mismatch sets [12],[10,14],[9,15,19],[10,12,14,16]; limits 0-4; both strands | site with count len(positions) when len<=limit else no site |
| 2 | 30 bp core, seeds 8/12/15/20, limits 1/2/3; exactly `limit` mismatches in the seed vs one more | exactly-limit binds with count=limit; one more gives [] |

## `test_the_seed_and_the_3_prime_end_and_too_short_primers`

**Result:** ✅ passed in 0.01 s · [line 260](../../tests/test_primer_binding.py#L260) · **Area:** Primer binding: seed and 3' end · **Checked against:** hand-computed · **Strength:** 🟢 strong

A mismatch anywhere in the 3' seed is rejected (and the last 3 bases must always match), a mismatch just outside the seed is crossed by the extension, and too-short primers or seeds never bind.

| # | Input | Expected output |
|---|---|---|
| 1 | 30 bp core, seeds 8-26; limit 0, mismatch at each seed base up to last 3 | no site |
| 2 | mismatch in last 3 bases, limit 5 | no site |
| 3 | mismatch just outside the seed (k=length-seed-1), limit 0 | site (200,230,strand,1) |
| 4 | seed 3 or 0; primer lengths 0-7; empty template; empty primer | [] |

## `test_tails_are_not_counted_and_a_matching_tail_extends_the_site`

**Result:** ✅ passed in 0.00 s · [line 293](../../tests/test_primer_binding.py#L293) · **Area:** Primer binding: 5' tails · **Checked against:** hand-computed · **Strength:** 🟢 strong

A non-matching 5' tail of any length is ignored and does not move the site; a tail that happens to match the template extends the annealed site and keeps the count correct.

| # | Input | Expected output |
|---|---|---|
| 1 | non-matching tail length 0-30, both strands, linear and circular, limits 0 and 2 | always (40,64,strand,0) |
| 2 | matching upstream tail of length 0-30 (forward and reverse) | site grows to (100-extra,124,1,0) / (100,124+extra,-1,0) |
| 3 | one mismatch inside the matching tail crossed by following matches | (60,124,1,1) / (100,164,-1,1) |

## `test_repeats_tandem_repeats_and_palindromes_are_returned_once_in_order`

**Result:** ✅ passed in 0.00 s · [line 326](../../tests/test_primer_binding.py#L326) · **Area:** Primer binding: repeats and order · **Checked against:** hand-computed · **Strength:** 🟢 strong

Repeated and tandem-repeat binding sites are each returned once and in scan order, and a palindromic primer binds both strands at the same place once each.

| # | Input | Expected output |
|---|---|---|
| 1 | template with 3 forward copies + 1 reverse copy of a 24 bp core; primer = core, limit 0 | 4 sites at the known offsets with correct strands; no duplicates; sorted by start; revcomp primer gives same places with strands swapped |
| 2 | tandem 'GATTACA'*8, primer 'GATTACA'*3, limit 0 | 6 overlapping sites at 30+7i..51+7i; no duplicates; == reference_sites |
| 3 | 22 bp palindrome, limit 0; also random 2000 bp with seed 8 limit 3 | palindrome at (30,52) both strands; random case: no duplicates, anchors in scan order |

## `test_sites_across_the_origin_and_mismatches_on_each_side_of_it`

**Result:** ✅ passed in 0.00 s · [line 373](../../tests/test_primer_binding.py#L373) · **Area:** Primer binding: origin · **Checked against:** hand-computed · **Strength:** 🟢 strong

A primer spanning a circular template's origin is found (on both strands and after rotation) but not on a linear template, and mismatches on either side of the origin are counted correctly.

| # | Input | Expected output |
|---|---|---|
| 1 | 200 bp, 24 bp primer straddling origin, both strands, circular vs linear | circular (188,212,strand,0); linear []; rotated template shifts the coords |
| 2 | mismatches at (left,right) pairs both inside the 15 bp seed, limit 2 then 1 | limit 2: (188,212,strand,2); limit 1: []; single-side mismatch, limit 1: count 1 |

## `test_template_ends_linear_versus_circular_and_primers_longer_than_the_template`

**Result:** ✅ passed in 0.02 s · [line 404](../../tests/test_primer_binding.py#L404) · **Area:** Primer binding: template ends · **Checked against:** independent implementation · **Strength:** 🟢 strong

A site crossing a linear template's end is not found but is on a circle; matching vs non-matching run-off tails behave; and primers longer than (or far shorter than) the template never panic and match the reference model.

| # | Input | Expected output |
|---|---|---|
| 1 | 200 bp, primer across the join, both strands, linear vs circular | linear []; circular (190,214,1,0)/(190,214,-1,0) |
| 2 | non-matching run-off tail vs matching run-off tail, forward and reverse | non-matching: trimmed to in-bounds site; matching: extends through origin only on circle |
| 3 | templates n in {1,2,5,8,12,20}, primer lengths 0-40, circular/linear, seeds {4,8,15,30}, limit 1 | == reference_sites; every site in-bounds (0<=s<n, 0<e-s<=n); length<8 gives [] |
| 4 | primer longer than a circular template; primers length 8..SEED | site <= template length; whole-primer seed, one mismatch allowed except last 3 bases |

## `test_lower_case_u_and_other_letters_in_primer_or_template`

**Result:** ✅ passed in 0.00 s · [line 457](../../tests/test_primer_binding.py#L457) · **Area:** Primer binding: alphabet · **Checked against:** hand-computed · **Strength:** 🟢 strong

Case is ignored, U reads as T, and any non-ACGT letter (N or IUPAC) in primer or template is treated as a plain differing base, counted as a mismatch and never matching.

| # | Input | Expected output |
|---|---|---|
| 1 | 20 bp core flanked by T's; core, lowercase core, lowercase template, U-for-T core | all give (20,40,1,0) |
| 2 | N or R or Y at a mid position, limit 0 vs 1 | limit 0 []; limit 1 (20,40,1,1) |
| 3 | N or R at the 3' end, limit 3 | [] |
| 4 | N in the template, limit 0 vs 1 | limit 0 []; limit 1 (20,40,1,1) |

**Review notes:**

- no docstring (minor)

## `test_the_extension_stops_three_below_the_best_and_follows_the_ends_and_the_mirror`

**Result:** ✅ passed in 0.04 s · [line 485](../../tests/test_primer_binding.py#L485) · **Area:** Primer binding: extension · **Checked against:** hand-computed · **Strength:** 🟢 strong

The annealing extension scores +1/-2, keeps the best-scoring length and stops 3 below the best; it respects template ends and wraps on a circle; and reverse-complementing the template mirrors every site onto the other strand.

| # | Input | Expected output |
|---|---|---|
| 1 | 10 match/mismatch patterns (e.g. mmmxmmm->ext 7 mm 1; mmxmm->ext 2; xmm->ext 0; mxmxmmmmm->ext 9 mm 2), both strands | site == (computed start/end, strand, mismatches) for each |
| 2 | extension reaching before start / past end, linear vs circular; 200 random sub-primers | linear clipped to (0,20)/(n-20,n); circular wraps; all sites 0<=s<e<=n |
| 3 | 600 random templates vs their revcomp (seeds 23), circular and linear, limits 0-2 | mirrored coords/count/length match; >300 sites checked |
| 4 | lone mismatch at k=4 vs k=2 both strands | k=4 crossed count 1; k=2 becomes a tail (103,124,1,0)/(100,121,-1,0) |

## `test_mismatch_positions_are_template_coordinates`

**Result:** ✅ passed in 0.01 s · [line 582](../../tests/test_primer_binding.py#L582) · **Area:** Primer binding: mismatch display · **Checked against:** hand-computed · **Strength:** 🟢 strong

The mismatch positions the Editor exposes for drawing are template coordinates on both strands, with the tail length recorded and positions correct across the origin.

| # | Input | Expected output |
|---|---|---|
| 1 | 400 bp, 30 bp core, tails ''/'GGATCC'/12 bp, mismatch at template base 100+k for k in {3,5,12,20,25}, forward and reverse primers | f.mismatch_at==[100+k], mismatches 1, start 100, tail len; r same with end 130 |
| 2 | exact primers | mismatch_at [] |
| 3 | circular 300 bp, primer across origin with mismatches at across[4] and [20] (template n-11 and 5), both strands | spans [(285,300),(0,15)], mismatch_at [289,5], count 2; exact primer count 0 |

## `test_primer_pair_products_linear_circular_and_across_the_origin`

**Result:** ✅ passed in 0.01 s · [line 630](../../tests/test_primer_binding.py#L630) · **Area:** PCR product · **Checked against:** hand-computed · **Strength:** 🟢 strong

A primer pair's PCR product is computed with correct range, sequence (tails included), and circularity for linear, circular, upstream-reverse and across-origin cases, mismatched primers alter the product, and the shortest valid product is auto-picked.

| # | Input | Expected output |
|---|---|---|
| 1 | 500 bp linear, F at 100..124 tail 'GGAATTCC', R at 300..324 tail 'CCAAGCTTGG' | range (100,324); product = tail+template[100:324]+rc(tail), length 242, linear; dialog shows '242 bp' and same product |
| 2 | R upstream of F, linear vs circular | linear range None ('not downstream'); circular range (350,574), product length 228 |
| 3 | two forward primers only | dialog.pair() is None |
| 4 | F across origin, R downstream with non-matching tail | range (n-12,n+124); product template[n-12:]+template[:124]+rc(tail), length 142 |
| 5 | reverse primer with one mismatch at template base 108 | product carries the mutated base (template[50:100]+mutated) |
| 6 | four primers F1/F2/R1/R2 | pair picks F2,R1 (shortest, 124 bp) |

## `test_genome_search_in_a_fasta_file_finds_the_same_sites_as_the_search_in_memory`

**Result:** ✅ passed in 0.04 s · [line 708](../../tests/test_primer_binding.py#L708) · **Area:** Primer binding: FASTA genome · **Checked against:** independent implementation · **Strength:** 🟢 strong

Primer binding read from a multi-record FASTA file matches the in-memory search and the reference model, including a cross-strand copy, partial copies that are not sites, and regardless of line width or CRLF endings.

| # | Input | Expected output |
|---|---|---|
| 1 | genome.fa with chr1 (900 bp) and chr2 (1400 bp, a site at 300 plus a reverse copy at 700, a 12 bp partial copy in chr1); primers primer/revcomp/mismatch/untailed; limits 0-2 | doc.binding_sites == sites == reference_sites for every case |
| 2 | fa.search_motif for the annealing part | exact 0-mismatch sites match binding_sites of full-length hits |
| 3 | chr2 vs chr1 with the tailed primer, limit 0 | chr2 has 2 sites; chr1 (only a partial copy) has [] |
| 4 | same genome written at width 70 with CRLF | binding_sites identical to the LF/width-60 file |

## `test_non_specific_priming_in_the_sequence_and_in_the_genome`

**Result:** ✅ passed in 0.01 s · [line 742](../../tests/test_primer_binding.py#L742) · **Area:** Primers: non-specific priming check · **Checked against:** hand-computed · **Strength:** 🟢 strong

The mis-priming check warns (yellow) when a primer's 3' end matches elsewhere and fails (red) when the whole primer does, both in the open sequence and against a set genome, excluding the primer's own site and the sequence/genome being identical.

| # | Input | Expected output |
|---|---|---|
| 1 | homology_settings (12,0); 24 bp primer unique in 700 bp | 'Non-specific (sequence)' == 'ok' |
| 2 | last 12 bases copied elsewhere after 3 differing bases, forward and reverse | 'warn', text contains '404-415 ->'; its own 'Binds here' stays 'ok'; warn colour in html |
| 3 | whole primer copied a second time | 'bad', text '2 places'; bad colour '#c62828'; worst() 'bad' |
| 4 | lowercase 5' tail; reverse primer | tail excluded, unique stays 'ok'; reverse judged the same |
| 5 | genome unset; then genome ok/warn (3' end)/bad (whole revcomp); genome == the sequence itself | unset 'info'; ok/warn/bad respectively; self-genome 'info' |

## `test_our_sites_are_the_binding_sites_stored_in_each_real_dna_file`

**Result:** ✅ passed in 0.00 s · [line 800](../../tests/test_primer_binding.py#L800) · **Area:** Primer binding: real file regression · **Checked against:** regression value · **Strength:** 🟢 strong

For a real sample DNA file, every computed binding site (0 mismatches) corresponds to a primer_bind feature stored in the file, ending at the same 3' base on the same strand.

| # | Input | Expected output |
|---|---|---|
| 1 | sample DNA_FILE primers, binding_sites(sequence, 0) | each site matches a stored primer_bind by name, strand, covered bases (subset) and 3' base |

**Review notes:**

- two runtime pytest.skip at lines 804 and 807 (file missing, or no primers) silently skip in CI if the sample file is absent or primer-free, so the test may never actually run
- no docstring

**Fixed by the audit of November 2026:**

- The test now counts the sites it compares and asserts at least one per stored primer (15 on the demo .dna file), so it cannot pass by comparing nothing. The skips stay: they are reported as skips, and the demo file is in the repository.

## `test_primer_sites_written_on_save_do_not_change_the_document_or_its_undo_stack`

**Result:** ✅ passed in 1.61 s · [line 827](../../tests/test_primer_binding.py#L827) · **Area:** Primers: save without mutating doc · **Checked against:** round-trip/self-consistency · **Strength:** 🟢 strong

Saving writes primer_bind features to the GenBank file but never adds them to the live document, leaves sequence/features/history/modified untouched, writes identical bytes on re-save, and the written sites follow later edits.

| # | Input | Expected output |
|---|---|---|
| 1 | circular 600 bp, a feature, an insert (undo+redo), two added primers pFwd/pRev; snapshot then save | doc sequence/features/history/circular/length unchanged; modified False; no primer_bind in the live doc |
| 2 | reload saved file | primer_bind [('pFwd',((100,124),),1),('pRev',((400,424),),-1)]; sequence unchanged |
| 3 | save again | identical bytes; snapshot unchanged |
| 4 | undo stack | 4 steps (feature, insert, two primers); save is not a step; full undo restores edited sequence and empties features/db/sites |
| 5 | add primer, insert 10 bp at 0, save t2.gb | written site moved to (110,134) |
