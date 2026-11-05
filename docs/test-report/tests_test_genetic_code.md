# `tests/test_genetic_code.py`

[← all test files](README.md) · source: [`tests/test_genetic_code.py`](../../tests/test_genetic_code.py)

## `test_the_embedded_tables_match_biopython_and_differ_where_they_should`

**Result:** ✅ passed in 0.01 s · [line 64](../../tests/test_genetic_code.py#L64) · **Area:** Genetic codes · **Checked against:** independent implementation · **Strength:** 🟢 strong

Checks every embedded NCBI translation table (except 27, 28, 31) codon by codon, including start and stop sets, on both strands against Biopython. Rust and Python must hold the same tables.

| # | Input | Expected output |
|---|---|---|
| 1 | table ids | genetic_code.TABLES ids == Biopython ids minus {27,28,31} |
| 2 | every table x 64 codons via translate(), codon_map(), bottom-strand codons() and whole-string translate | equal Biopython amino acid or '*'; start/stop sets equal Biopython |
| 3 | TGA in 1/4/25, AGA in 1/2, ATA in 1/2, CTG in 12/1 | '*'/'W'/'G', 'R'/'*', 'I'/'M', 'S'/'L' |
| 4 | table 11 vs 1; table 99 | same codons, strictly more starts; TGA -> '*' fallback |
| 5 | opendna_core.genetic_codes(); citations | ids equal; 64-long strings; 'Elzanowski' and 'wprintgc' cited |

## `test_orfs_start_and_stop_follow_the_table`

**Result:** ✅ passed · [line 113](../../tests/test_genetic_code.py#L113) · **Area:** Genetic codes: ORFs / protein search · **Checked against:** hand-computed · **Strength:** 🟢 strong

Checks that the ORF finder uses alternative start codons only when allowed and stops from the chosen table, and that protein search translates with the table.

| # | Input | Expected output |
|---|---|---|
| 1 | 'CCGTGAAAAAATAACC' orfs(4), (4,11), (4,11,True), (4,1,True) | [], [], [(2,14,1)], [] |
| 2 | 'CCTTGAAAAAATAACC' table 1 alt on/off; rc of GTG seq table 11 alt | [(2,14,1)] / []; [(2,14,-1)] |
| 3 | 'ATGTGAAAATAA' table 1 / 4; 'ATGAGAAAATAA' table 2 min 1 | [] / [(0,12,1)]; [(0,6,1)] |
| 4 | find_protein('AAATGAAAAAA','WK', table 4 / 1) | truthy / empty |

**Review notes:**

- Lines 130-132: protein search result is only checked for truthiness, not position; the comment 'frame +3 starts at base 2' is wrong (the TGA AAA hit is in frame +1, 0-based index 3).

## `test_a_cds_is_read_with_its_own_table_frame_and_exceptions`

**Result:** ✅ passed in 0.03 s · [line 136](../../tests/test_genetic_code.py#L136) · **Area:** Genetic codes: CDS translation · **Checked against:** hand-computed · **Strength:** 🟢 strong

Checks that a CDS uses its own /transl_table, /codon_start and /transl_except, shows a start codon as M, and flags a /transl_except that cannot be reproduced.

| # | Input | Expected output |
|---|---|---|
| 1 | 'ATGTGAGTGTAA' CDS 1..12 with transl_table=4, under settings 1 and 11; plain CDS under 4 / 1 | 'MWV*' both; 'MWV*' / 'M*V*' |
| 2 | codon_start=2 on 1..12; 'TTGAAATAA' complete; '<13..21' partial | 'CEC'; 'MK*'; 'LK*' |
| 3 | 'ATGTGATAA' with (pos:4..6,aa:Sec); with (pos:4..5,aa:TERM) | 'MU*'; cds_note says 'not reproducible', drawn 'M**' |
| 4 | hover describe() of table-4 CDS and bad transl_except CDS | contains 'Table 4'; 'translation not reproducible' |

## `test_the_setting_persists_shows_in_the_dialog_and_drives_views_and_exports`

**Result:** ✅ passed in 0.28 s · [line 175](../../tests/test_genetic_code.py#L175) · **Area:** Genetic codes: settings / GUI / export · **Checked against:** GUI state · **Strength:** 🟢 strong

Checks that the chosen genetic code persists, is shown with its stops and starts in the settings dialog, and drives the protein view, exports and the sequence view's CDS lane.

| # | Input | Expected output |
|---|---|---|
| 1 | LABELS[0]; table set to 'x'; choose 11 alt, clear cache | '1 - Standard'; table_id 1; 11 and alt True after re-reading settings.ini |
| 2 | SettingsDialog with table 1, then table 2 | combo has len(IDS) items; label 'Stop codons: TAA, TAG, TGA . Start codons: TTG, CTG, ATG'; then 'TAA, TAG, AGA, AGG . ATT, ATC, ATA, ATG, GTG'; alt starts unchecked |
| 3 | CDS 3..14 of 'CCATGTGAAAATAACC', table 1 then 4 | protein_of 'M*K' then 'MWK'; translate('ATGTGA') 'M*' then 'MW'; ProteinView title gets 'Table 4' |
| 4 | non-CDS selection 2..14 under table 4 | src kind 'selection', codons 'MWK*' |
| 5 | MainWindow CDS lane, table 4 then 1 | 'MWK*' then 'M*K*' |

**Review notes:**

- Line 22 imports helpers (open_doc, close) from test_cas9.py: a test module depends on another test module. Settings are shared on disk across modules; the code_setting fixture restores them.

## `test_refresh_updates_stale_translations_and_leaves_the_rest_alone`

**Result:** ✅ passed in 0.02 s · [line 255](../../tests/test_genetic_code.py#L255) · **Area:** Genetic codes: /translation refresh · **Checked against:** independent implementation · **Strength:** 🟢 strong

Checks the /translation refresh: which CDSs it updates or skips, that only /translation changes, that it is one undo step, and that the Rust refresh matches what the Python views draw.

| # | Input | Expected output |
|---|---|---|
| 1 | 7 CDSs (stale 'MQ', codon_start 2, pseudo, Sec transl_except, no /translation, table 4, partial), table 1 | status g/cs2/exc/t4 'updated', pseudo 'pseudo', partial 'partial'; g 'MK', cs2 'MK', t4 'MWK', exc 'MU'; pseudo/partial stay 'MQ'; none stays None |
| 2 | other qualifiers; history depth; second run; undo | unchanged; +1 step; second run unchanged/no step; undo restores all qualifiers |
| 3 | 'ATGTGAAAATAA' translation 'M', table 4 then 1 | 'updated' -> 'MWK'; then 'other' (internal stop), still 'MWK' |
| 4 | table 11: complement CDS with Sec, codon_start 2, TTG start | each /translation equals exports.protein_of; rev 'MU', ttg 'MK' |
| 5 | summarize() of report lists | 'Updated /translation of 3 CDSs, skipped 1 pseudo'; '... of 1 CDS' |

## `test_features_menu_toggle_save_autosave_and_update_now`

**Result:** ✅ passed in 0.41 s · [line 312](../../tests/test_genetic_code.py#L312) · **Area:** Genetic codes: Feature menu / save · **Checked against:** GUI state · **Strength:** 🟢 strong

Checks the 'Keep /translation' menu toggle persists, that an explicit save refreshes translations but autosave never does, and that 'Update now' acts on the selected CDS or all as undo steps.

| # | Input | Expected output |
|---|---|---|
| 1 | new window; trigger the toggle; open another window | checked and checkable by default, text 'Keep /translation'; then setting False and new window unchecked |
| 2 | refresh_doc, modified, autosave() | version and history unchanged, g still 'MQ', '.autosave.gb' file contains 'MQ' |
| 3 | _write(path) twice; reread file; undo | g 'MK', not modified, file has 'MK' and pseudo 'MQ', status 'Updated /translation of 4 CDSs, skipped 1 partial, 1 pseudo'; no new step on 2nd save; undo gives 'MQ' and modified |
| 4 | keep off, _write | g stays 'MQ' |
| 5 | select feature 0, update_translations(); then none selected | only g 'MK'; then cs2 'MK', t4 'MWK'; +2 steps; status '... 3 CDSs ...'; undo reverts cs2 only |
