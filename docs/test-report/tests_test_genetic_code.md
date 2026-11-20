# `tests/test_genetic_code.py`

[← all test files](README.md) · source: [`tests/test_genetic_code.py`](../../tests/test_genetic_code.py)

## `test_translation_tables_match_biopython`

**Result:** ✅ passed in 0.00 s · [line 64](../../tests/test_genetic_code.py#L64) · **Area:** Genetic codes and translation · **Checked against:** independent implementation · **Strength:** 🟢 strong

**Why it matters:** Without it, translations, protein search and ORFs for non-standard organisms (mitochondria, mycoplasma, bacteria with alternative starts) could use wrong codon assignments.

Compares every codon of every embedded NCBI table (except 27, 28, 31) with Biopython on both strands, checks start/stop sets, Rust/Python table agreement, protein search over all 64 codons, incomplete codons and N, and table-dependent ORF starts and stops.

| # | Input | Expected output |
|---|---|---|
| 1 | embedded tables | same IDs as Biopython minus 27, 28, 31 |
| 2 | 64 codons x each table: translate, codon_map, reverse-strand codons(), whole string | equal Biopython; start and stop sets equal |
| 3 | TGA in 1/4/25; AGA, ATA in 1/2; CTG in 12/1 | * / W / G; R vs * ; I vs M; S vs L |
| 4 | table 11 vs 1; unknown table 99 | same codons, more starts; falls back to standard |
| 5 | genetic_codes() from Rust | same IDs, 64 entries per table |
| 6 | each codon and its rc with find_protein | (0,3,+1,0,1) and (0,3,-1,0,1); a different amino acid does not match |
| 7 | all 64 codons as one protein, upper and lower case | [(0,192,1,0,1)] |
| 8 | ATGGCCAAAGG; frame 2 WPK; ATGNNNGCC | MAK found, MAKG not; WPK at (1,10,frame 2); NNN is X: MXA found, MAA only with 1 mismatch |
| 9 | CC GTG AAA AAA TAA CC; TTG variant; rc | ORF only with table 11 and alt starts on; TTG with table 1 alt; reverse strand (2,14,-1) |
| 10 | ATG TGA AAA TAA; ATG AGA AAA TAA | table 4 ORF (0,12); table 2 stops at AGA (0,6) |
| 11 | find_protein 'WK' in AAATGAAAAAA | found with table 4, not table 1 |

**Review notes:**

- The line asserting str(Seq(codon).translate()) == aa only checks Biopython against itself.
- The ORF table/start section partly overlaps test_core_sequence.py::test_orfs_match_a_brute_force_scan (that one covers table 1 only).

## `test_cds_translation_and_translation_refresh`

**Result:** ✅ passed in 0.32 s · [line 178](../../tests/test_genetic_code.py#L178) · **Area:** CDS translation and /translation refresh · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** Without it, CDS amino acids shown in the views and exported, and the /translation qualifiers rewritten by the refresh command, could ignore /transl_table, /codon_start or /transl_except, or overwrite pseudo/partial genes.

Checks cds_codons, exports.protein_of/translate and Editor.refresh_translations on small hand-made CDSs under different table settings, including the undo step and agreement between the Rust refresh and the Python drawing/export.

| # | Input | Expected output |
|---|---|---|
| 1 | ATG TGA GTG TAA with /transl_table=4 under settings 1 and 11; plain CDS under 4 and 1 | MWV* regardless; plain follows setting: MWV* / M*V* |
| 2 | /codon_start=2; TTG start complete; partial <13..21 | CEC; MK*; LK* (no M for partial) |
| 3 | /transl_except Sec at 4..6; TERM at 4..5 | MU*; note says 'not reproducible', drawn M** |
| 4 | exports under table 1 / 4 | protein_of 'M*K' / 'MWK'; translate('ATGTGA') 'M*' / 'MW' |
| 5 | refresh on 7 CDSs (stale, codon_start 2, pseudo, transl_except, none, table 4, partial) | statuses updated/updated/pseudo/updated/updated/partial; values MK, MK, MU, MWK; pseudo and partial keep MQ; none stays without /translation |
| 6 | other qualifiers; second run; undo | unchanged; second run adds no step; one undo restores everything |
| 7 | setting 4 with CDS without own table; then setting 1 | 'updated' to MWK; then 'other' (internal stop) and left alone |
| 8 | table 11: reverse CDS with Sec, codon_start 2, TTG start | refresh equals exports.protein_of; rev 'MU', ttg 'MK' |
| 9 | summarize report | 'Updated /translation of 3 CDSs, skipped 1 pseudo'; '... of 1 CDS' |

**Review notes:**

- Changes shared tool settings during the test (restored by the fixture).
