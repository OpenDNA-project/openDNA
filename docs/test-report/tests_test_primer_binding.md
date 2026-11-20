# `tests/test_primer_binding.py`

[← all test files](README.md) · source: [`tests/test_primer_binding.py`](../../tests/test_primer_binding.py)

## `test_binding_agrees_with_the_documented_rules_written_in_python`

**Result:** ✅ passed in 0.12 s · [line 173](../../tests/test_primer_binding.py#L173) · **Area:** Primer binding · **Checked against:** in-test reference model · **Strength:** 🟢 strong

**Why it matters:** Without it, the primer binding sites shown on the map could deviate from the documented seed/extension rules on any of thousands of combinations of mismatches, tails, strands and topologies.

Compares find_primers with a Python re-implementation of the documented rules on 1500 random template/primer combinations; checks Document.binding_sites equals find_primers; checks binding on a genome FASTA file (any line width, CRLF) against the in-memory search and the file's own motif search.

| # | Input | Expected output |
|---|---|---|
| 1 | 1500 random templates 12-90 bp with tandem repeats, circular or not, primers 6-30 bp with 6% mutations, optional rc and random tail, limit 0-3, seed 4/8/12/15/20 | identical site list and order to reference_sites; > 500 sites in total |
| 2 | 300 bp template, primer 100..124 in 4 variants, limits 0-2, seed 15 and 12, linear and circular | Document.binding_sites == find_primers |
| 3 | FASTA genome chr1/chr2 with a tailed primer site at chr2:300, an rc copy at 700, a 12-base partial copy in chr1 | Document.from_fasta sequence equals input; binding_sites == find_primers == reference for 4 primer variants x 3 limits |
| 4 | exact sites of the 24-base annealing part | equal the file's search_motif hits |
| 5 | primer on chr2 / chr1, limit 0 | 2 sites / none |
| 6 | same genome written 70 columns, CRLF | same sites as the 60-column file |

## `test_binding_mismatches_seed_tails_and_extension`

**Result:** ✅ passed in 0.04 s · [line 234](../../tests/test_primer_binding.py#L234) · **Area:** Primer binding · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** Without it, specific binding rules (where a mismatch is allowed, how far the annealed part extends, how tails are handled) could change and primers would be shown binding where they do not, or not where they do.

Checks each documented binding rule with exact expected sites: every single-mismatch position, the mismatch limit in the seed, the exact 3' end, minimum seed/primer lengths, non-binding and matching tails, the extension scoring patterns, template ends, and mirror symmetry.

| # | Input | Expected output |
|---|---|---|
| 1 | primers 16/20/24/30 bp with one mismatch at each position k, both strands | site as in expected_single_mismatch: none in last 3 bases; counted in seed; crossed in extension; 5'-most 3 bases become a tail |
| 2 | mutations at 4, 9, 15, 20 forward vs reverse primer | same (100,124,1 mismatch), strands +1 / -1 |
| 3 | 1-4 mismatches at 12/10,14/... with limits 0-4 | found with count iff count <= limit |
| 4 | exactly `limit` seed mismatches vs limit+1, seeds 8/12/15/20 | found with count = limit; one more: not found |
| 5 | seed 8..26: mismatch inside seed at limit 0; in last 3 at limit 5; just outside seed at limit 0 | none; none; found with 1 mismatch |
| 6 | seed 3 or 0; primers shorter than 8; empty template/primer | no sites |
| 7 | non-matching 5' tail of 0-30 bases, both strands, linear/circular | site stays (40,64,strand,0) at limits 0 and 2 |
| 8 | matching upstream tail of 0-30 bases; one mismatch in it | site extended to include it; mismatch counted as 1 |
| 9 | extension patterns mmm, mmmxmmm, mmxmm, mmxxmmmm, xmmm, xmm, mxmxmmmmm, mxmxmxm, mmmmmxxmmmmmmmm, '' | extensions 3,7,2,2,4,0,9,1,5,0 with 0/1/0/0/1/0/2/0/0/0 mismatches |
| 10 | primer reaching 10 bases past template start/end | linear: clipped at 0 or n; circular: through the origin; 200 random primers never outside 0..n |
| 11 | 600 random cases per topology: template vs its reverse complement | mirrored site set; > 300 sites checked |
| 12 | mismatch 4 vs 2 bases from 5' end, both strands | (100,124,1 mm) vs tail: (103,124,0) / (100,121,0) |

**Review notes:**

- Partly redundant with the random comparison in test_binding_agrees_with_the_documented_rules_written_in_python, but it pins each rule by name, which makes failures easier to read.

## `test_binding_origin_ends_repeats_and_letters`

**Result:** ✅ passed in 0.02 s · [line 419](../../tests/test_primer_binding.py#L419) · **Area:** Primer binding · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** Without it, primers on circular plasmids across the origin, repeated or palindromic sites, lower-case/RNA/IUPAC primers or N in the template could produce missing, duplicated or wrong sites.

Checks site lists for repeated, tandem and palindromic sites, sites across the origin (with mismatches on each side), linear vs circular tails, primers longer than the template, and letter handling (case, U, X, N, IUPAC codes in primer and template).

| # | Input | Expected output |
|---|---|---|
| 1 | 24-mer placed at 4 places (one rc) in AT filler | 4 sites in scan order with right strands; rc primer gives the same with strands swapped |
| 2 | GATTACA x 8, primer = unit x 3 | 6 sites every 7 bases, equal to reference_sites |
| 3 | 22-base palindrome | two sites at (30,52), one per strand |
| 4 | random 2000 bp, seed 8, limit 3 | no duplicate sites, ordered by seed anchor |
| 5 | 24 bases across origin, both strands; rotated template | circular (n-12,n+12); linear none; rotation shifts by 40 |
| 6 | 2 mismatches either side of the origin | found with 2 at limit 2, not at limit 1; each single one found with 1 |
| 7 | non-matching / matching tail over a template end | linear always (0,24); circular extends through origin only when it matches |
| 8 | templates 1-20 bp, primers 0-40, seeds 4/8/15/30 | equal reference_sites; coordinates in range; < 8 bases none |
| 9 | primer longer than circular 12 bp template | site at most 12 bases |
| 10 | primers 8..15 bases (whole primer is seed) | exact found; 5' mismatch needs limit 1; last base never differs |
| 11 | lower case, U, X in primer | same site; X is a mismatch |
| 12 | N, Y, R in primer; N/Y at 3' end | N and Y match their bases with 0 mismatches; R vs T is a mismatch; Y vs A rejected at 3' end |
| 13 | N in template | a mismatch; template N matches only primer N |
| 14 | primer with R, Y, N at 3 positions vs 3 variants; rc template; disallowed base | 3 sites 0 mm; 3 reverse sites; the disallowed variant dropped |

## `test_binding_matches_the_primer_sites_stored_in_the_real_dna_file`

**Result:** ✅ passed in 0.00 s · [line 583](../../tests/test_primer_binding.py#L583) · **Area:** Primer binding, real file · **Checked against:** independent implementation · **Strength:** 🟢 strong

**Why it matters:** Without it, the primer sites the app computes could disagree with those SnapGene stored in a real lab file, which is the most direct check that users see the same sites as in their other tool.

Opens the sample .dna file with 15 primers and their stored primer_bind features, and checks in both directions that computed sites (limit 0) match the stored sites with the same 3' end, the only known differences being oAC812/oAC813 (SnapGene adds a chance-matching 5' tail).

| # | Input | Expected output |
|---|---|---|
| 1 | sample .dna file | 15 primers and 15 stored primer_bind features |
| 2 | each computed site | inside a stored site of the same primer and strand, sharing its 3' base; at least 15 compared |
| 3 | each stored site | found on its strand with the same 3' end |
| 4 | sites shorter than stored | exactly oAC812 and oAC813 |

**Review notes:**

- Skips if the sample file is missing; depends on one specific private sample file.

## `test_editor_primer_sites_pcr_product_specificity_and_save`

**Result:** ✅ passed in 0.03 s · [line 643](../../tests/test_primer_binding.py#L643) · **Area:** Primer sites in the editor, PCR product, specificity · **Checked against:** hand-computed · **Strength:** 🟡 adequate

**Why it matters:** Without it, the editor could colour the wrong mismatched bases, build a wrong PCR product (tails, origin, mismatches), give wrong non-specific priming warnings, or alter the document and undo history when saving primer sites.

Uses Editor objects with primers to check mismatch positions in template coordinates, PCR product ranges and sequences, the non-specific priming checks against the sequence and a genome, and that primer sites are written to the saved file without touching the document or undo stack.

| # | Input | Expected output |
|---|---|---|
| 1 | mismatch_indices('ACGTNR','ACCTGA') | [2] (N and R pair with what they stand for) |
| 2 | 30 bp core at 100 with mismatch at k in 3/5/12/20/25 and tails of 0/6/12 bases, fwd and rev | mismatch_at [100+k], 1 mismatch, start 100/end 130, tail length reported |
| 3 | 30 bases over origin with mismatches at 4 and 20 | spans [(n-15,n),(0,15)], mismatch_at [n-11, 5], 2 mismatches |
| 4 | F at 100..124 with 8-base tail, R at 300..324 with 10-base tail | product_range (100,324); product = tail + template[100:324] + rc(tail), length 242, linear |
| 5 | R upstream of F: linear; circular with GG/CC tails | None; (350,574) and product wrapping the origin, length 228 |
| 6 | F over origin, R at 100..124 with 6-base tail | product_range (n-12, n+124), length 142 |
| 7 | reverse primer with 1 mismatch | product carries the primer's base |
| 8 | non-specific check, sequence: unique; last 12 bases elsewhere; whole primer twice; tailed and rc | ok; warn with '404–415 →'; bad with '2 places'; tail ignored; rc same |
| 9 | non-specific check, genome: none; other genome; partial; rc copy; genome = the sequence | info; ok; warn; bad; info |
| 10 | editor with feature, edit, undo/redo, two primers, save twice | document, features, history unchanged; modified False; no primer_bind in doc; file holds pFwd (100,124,+1) and pRev (400,424,-1); identical bytes on second save |
| 11 | undo all steps | 4 steps; sequence/features/primers back to empty |
| 12 | edit after adding primer, then save | stored site moved to (110,134) |

**Review notes:**

- Bundles four unrelated features (mismatch display, PCR product, non-specific priming check, save) into one test; would read better split.
- Leftover comment '# the dialog (built, not run) says the same size...' with no code under it.
- Uses a module-level RNG shared across tests, so the random data depends on which tests ran before.

**Fixed by the audit of November 2026:**

- 13 Nov 2026: gap closed. The proposed checks were added, pass on the current code and fail under the mutation that slipped through before. Was: Order-dependent: the PCR product section's exact coordinates hold only for the rng state left by earlier tests; alone it fails (shown above), so it does not pin the rule it claims, just one lucky draw. \| 'A lower-case 5' tail is not part of what must match elsewhere': the only tailed check is on a unique template, where 'ok' comes out whether or not the tail is stripped (proven). \| Tails in the PCR/mismatch sections are random letters, not built to differ from the neighbouring template base, hence the chance extension.
