# `tests/test_app_logic.py`

[← all test files](README.md) · source: [`tests/test_app_logic.py`](../../tests/test_app_logic.py)

## `test_gene_groups`

**Result:** ✅ passed in 0.03 s · [line 111](../../tests/test_app_logic.py#L111) · **Area:** Gene groups · **Checked against:** hand-computed · **Strength:** 🟢 strong

Features sharing a locus_tag or /gene fold under their gene in table, map and display; lone features stay alone; and edits keep every feature on the same bases.

| # | Input | Expected output |
|---|---|---|
| 1 | 1500 bp (seed 11); gene gA 101..400, mRNA mA, CDS cA (all L1); promoter pA; misc solo; CDS orphan (L9, no gene); gene abc + CDS abcCDS (shared /gene); gene g2 (L2, no parts) | groups {L1, abc, L2}; L1 master gA, parts [mA, cA]; abc parts [abcCDS]; L2 parts [] |
| 2 | part_of and lone features | part_of keys = {mA, cA, abcCDS}; pA, solo, orphan neither part nor master |
| 3 | FeatureTable top rows | top = all non-part features; gA row carries 'L1' |
| 4 | default display (cds_only) | map {gA, pA, solo, orphan, abcCDS, g2}; abc not shown nor visible; abcCDS shown |
| 5 | cds_only=False; then collapse_genes=False | map {gA, pA, solo, orphan, abc, g2}; then all features |
| 6 | insert 40 A at 0, delete 300..310 | groups unchanged; cA bases = before[:110]+before[120:]; gA = before[:160]+before[170:]; mA 10 shorter; others unchanged |
| 7 | gene span vs its parts | each part inside the gene span |
| 8 | insert 'CCCCC' at 165 (inside gene, outside CDS) | cA bases unchanged; gA 5 longer |
| 9 | undo_steps(3) | returns 3; all bases and groups restored to start |

**Review notes:**

- comment at line 159 says the mRNA grows too, but mA length is not checked after the 'CCCCC' insert

## `test_display_filter_hides_exactly_the_types_switched_off`

**Result:** ✅ passed in 0.30 s · [line 173](../../tests/test_app_logic.py#L173) · **Area:** GUI: display filter · **Checked against:** GUI state · **Strength:** 🟢 strong

Unchecking a feature type in the Display panel hides exactly that type in all views but keeps it in the file and table; None/All, edits and per-document scope behave.

| # | Input | Expected output |
|---|---|---|
| 1 | MainWindow, gene_features(), cds_only=False; DisplayPanel | checkbox counts == Counter of kinds; all checked |
| 2 | uncheck promoter | hidden_kinds {promoter}; shown/visible/map exclude promoters; map_lanes == map features; all features still in doc and table |
| 3 | uncheck gene too | hidden {promoter, gene}; map has mA, cA, abcCDS and none of gA, abc, g2, pA; map == shown |
| 4 | recheck both | no hidden; map {gA, pA, solo, orphan, abc, g2} |
| 5 | None button, then All button | None: all kinds hidden, shown and map empty; All: none hidden |
| 6 | hide CDS, type 'A' at 5, undo; add terminator 1400..1420 | still {CDS} hidden; new terminator box checked; CDS box present but unchecked |
| 7 | second window, same features | its hidden_kinds empty and shown non-empty (per-document filter) |

**Review notes:**

- no docstring (minor)

## `test_history_texts`

**Result:** ✅ passed in 0.02 s · [line 230](../../tests/test_app_logic.py#L230)

*No written description yet: the docstring and the checks below are taken from the code as they are.*

History texts for each kind of edit (insert, delete, replace, paste, features, reverse complement, redone rows
after the marker) and for a circular sequence (rotation, linearising, edits across the origin).

Checks (input → expected output, as written in the test):

```python
assert texts[-2:] == [NOW, 'Reverse-complemented the whole sequence']
assert texts[0] == 'Moved base 1 (broke the circle elsewhere)'
assert texts[1] == 'Made the sequence linear'
assert texts[2] == 'Made the sequence circular'
assert texts[3] == 'Made the sequence linear (broke the circle elsewhere)'
assert texts[4] == 'Made the sequence circular'
assert texts[5] == 'Replaced 15 bp across the origin from 191 with GG'
assert last.startswith('Replaced 15 bp across the origin from ') and last.endswith(' (deleted)')
assert [t for _, _, t in group(ed.doc.history()[0])] == expected
assert [panel.list.item(i).text() for i in range(panel.list.count())] == ['Opened t.gb'] + expected + [NOW]
```

## `test_sidewinder`

**Result:** ✅ passed in 0.01 s · [line 329](../../tests/test_app_logic.py#L329) · **Area:** Assembly: Sidewinder · **Checked against:** hand-computed · **Strength:** 🟢 strong

Sidewinder design builds scarless 10 nt toeholds at each junction, coding/barcode oligos that tile the product, a simulated product carrying moved features, and stays linear even when circular is requested.

| # | Input | Expected output |
|---|---|---|
| 1 | 3 fragments (100/120/90 bp, seed loop 0..59 picks a clean one), analyse(parts, 2, False) | d.ok; product length 310; two toehold joints (a,b),(b,c) |
| 2 | joint seqs (TOEHOLD=10, half=5) | each 10 nt, == product[b-5:b+5] at bounds 100 and 220, unique, pos b-5 |
| 3 | coding rows | a/b/c-coding == product[:105], [105:225], [225:]; concatenate to product; a,b end with their toehold |
| 4 | barcode rows (reverse complement + 18 nt arm, none on last) | each starts with revcomp(coding); arm lengths [18,18,0]; first two arms differ; all oligos <= MAX_OLIGO (165) |
| 5 | simulate(parts, d) features | pro at (0,50); cds1 at (109,160); a/b/c at (0,100),(100,220),(220,310); sequence == product |
| 6 | analyse circular=True vs linear (seed 3) | both not circular, 2 joints, a 'Linear product' check; same product() |
| 7 | Gibson 2 fragments with auto_primers, circular | circular with 2 joints (loop closed) |
| 8 | analyse(one fragment, 2, False) | raises ValueError 'at least two' |

**Review notes:**

- seed-search loop (line 333) relies on one of 60 random draws giving d.ok; if assumptions changed it could fail non-deterministically, though the final assert catches it with a message

## `test_assembly_bad_input`

**Result:** ✅ passed in 0.08 s · [line 392](../../tests/test_app_logic.py#L392) · **Area:** Assembly: bad input · **Checked against:** hand-computed · **Strength:** 🟢 strong

Tiny or empty fragments produce a 'bad' check rather than an exception for every method and orientation, an empty source stretch yields a zero-base fragment, and a repeated Golden Gate overhang is flagged and blocks assembly.

| # | Input | Expected output |
|---|---|---|
| 1 | methods gibson/golden gate/sidewinder x circular/linear x middle fragment size 0/5/18 bp | analyse: not ok, has bad checks with non-empty text; auto_primers does not raise |
| 2 | after auto_primers, size 0 or sidewinder | not ok with bad checks (empty fragment can't join) |
| 3 | after auto_primers, non-sidewinder size 5/18 | 'Fragment length' == 'warn' naming fragment B |
| 4 | size 0 | fragment length 0 and seq '' |
| 5 | Frag on source with (5,5),(0,0),(50,20) ranges | seq '', core_seq '', tails (0,0) |
| 6 | AssemblyPanel with sources 300/0/250 bp | design not ok with bad checks; status starts '✗'; perform button disabled |
| 7 | Golden Gate ring with repeated overhang ['GGAG','GGAG','TCGC'], BsaI | a 'unique' bad check saying 'same overhang'; not ok |
| 8 | ring with three distinct overhangs ['GGAG','TCGC','CAGT'] | ok |

## `test_each_local_align_hit_is_its_own_row_in_the_alignment_view`

**Result:** ✅ passed in 0.19 s · [line 462](../../tests/test_app_logic.py#L462) · **Area:** GUI: local align view · **Checked against:** GUI state · **Strength:** 🟢 strong

Overlapping or repeated Local Align hits each get their own row in the genome view, labelled and numbered, and selecting a hit in the table or canvas selects it in the other.

| # | Input | Expected output |
|---|---|---|
| 1 | seq with a 60 bp unit repeated twice forward and once inverted; align_hits(unit) | 4 hits; ed.query_hits is dlg.hits |
| 2 | canvas 'Query hits' track layout | track height HITS_HEAD + 4*lane; 4 bars indexed 0..3, each on a distinct row at top+HITS_HEAD+i*lane |
| 3 | hit labels | label 1 starts '#1 '; one label carries the inversion mark; exactly one hit strand < 0 |
| 4 | click names column at row 3 y | hit_row_at(y)==2, hit_row_at(top+2)==-1; editor selection = hit 2 span; current_query_hit 2; table row 2 |
| 5 | select table row 0 | current_query_hit 0; editor selection = hit 0 span |
| 6 | condensed/expanded/short track styles; clear | no paint crash; cleared: no query_hits, current_query_hit -1 |

## `test_local_align_hits`

**Result:** ✅ passed in 0.37 s · [line 508](../../tests/test_app_logic.py#L508) · **Area:** Local align · **Checked against:** hand-computed · **Strength:** 🟢 strong

Local Align finds exact copies on either strand at the right coordinates, reports one mismatch for a single SNP, returns nothing for absent sequences, rejects bad queries, feeds a genome-view track, and handles a hit across a circle's origin.

| # | Input | Expected output |
|---|---|---|
| 1 | 3000 bp (seed 21); query ALIGN_SEQ[1000:1100] | one hit (1000, 1100, strand 1, identity 1.0); table '1,001-1,100', FORWARD_TEXT, '0 mismatches, 0 ins, 0 del', status '1 hit' |
| 2 | query revcomp of that piece | one hit (1000, 1100, -1); INVERTED_TEXT |
| 3 | sequence with one forward and one reverse copy | two hits sorted [(1000,1100,1),(1500,1600,-1)] |
| 4 | piece with base 40 mutated | one hit (1000,1100); mismatch at 1040 with the new base; no ins/del; identity ~0.99; cell '1 mismatch, 0 ins, 0 del' |
| 5 | genome view after a hit | query_hits length 1, 'Query hits' in canvas.kinds with a qhit at (1000,1100) |
| 6 | random queries absent from the sequence (seeds 31,32,33) | no hits, status 'No hit', table empty; no 'Query hits' track |
| 7 | query 'AC'; query 'ACGT'+'Z'*10 | 'Query too short' / a different non-'No hit' status; no hits |
| 8 | circular seq, query = last 30 + first 30 bp | one hit (2970, 3030, 1); table coord '2,971-30' (end past origin) |

## `test_primer_settings_and_check_defaults`

**Result:** ✅ passed in 0.01 s · [line 590](../../tests/test_app_logic.py#L590) · **Area:** Primers: settings and band · **Checked against:** hand-computed · **Strength:** 🟢 strong

Primer settings persist in the data folder (not the registry), the primer-check default limits equal the documented values, and the band function returns ok/warn/bad at its thresholds.

| # | Input | Expected output |
|---|---|---|
| 1 | primer_db.settings() | IniFormat at data_dir()/'primers.ini'; PP.settings is that function |
| 2 | set pcr/mg_mm = 3.0, read back via conditions() | conditions()['mg_mm'] == 3.0 |
| 3 | tool_settings.TOOLS[PP.CHECKS] defaults | equal DEFAULT_LIMITS |
| 4 | band(v,52,65,7) for v in 51.9,52,65,65.1,72,72.1,45,44.9 | ['warn','ok','ok','warn','warn','bad','warn','bad'] |

## `test_primer_check_levels`

**Result:** ✅ passed in 0.01 s · [line 618](../../tests/test_app_logic.py#L618) · **Area:** Primers: check thresholds · **Checked against:** hand-computed · **Strength:** 🟢 strong

The Tm, GC, 3' clamp, hairpin and self-dimer primer checks flip between ok/warn/bad at their configured limits, verdicts roll up, and too-short primers get no verdict.

| # | Input | Expected output |
|---|---|---|
| 1 | Tm values 51.9..20.0 with a fixed primer (Tm monkeypatched) | warn/ok/ok/ok/warn/warn/bad/warn/bad/bad; after moving tm_low=55,margin=2: 53.0 warn, 52.9 bad |
| 2 | GC content from gc_primer(k) for k in 4..20 of 20 bases | {4:bad,5:bad,6:warn,7:warn,8:ok,12:ok,13:warn,14:warn,15:bad,20:bad}; moved gc_low=45,margin=5: k 8/9/10 warn/ok/ok; mixed-case 50% ok |
| 3 | 3' clamp with 0..5 trailing G/C | {0:warn,1:ok,2:ok,3:ok,4:warn,5:warn}; with clamp_max=4 the 4-GC end becomes ok |
| 4 | hairpin/self-dimer on clean/loop/strong primers | clean ok/ok; loop Hairpin warn, strong Hairpin bad; two Self-dimer cases warn/bad |
| 5 | edge: set limit to structure(primer) +/- 0.01 | verdict flips warn<->bad exactly at the primer's own dG for hairpin and self-dimer (3' and anywhere limits) |
| 6 | primer_checks on 'ACGT' (4 bases); worst(...) rollups; checks_html | [('Primer','info','Type at least 8 bases.')]; worst gives warn/bad/info/info; html includes ok and bad marks and '60.0 °C' |

**Review notes:**

- Tm is monkeypatched to a fixed value, so the real melting_temp calculation is not exercised here (by design; covered elsewhere)

## `test_feature_table`

**Result:** ✅ passed in 0.01 s · [line 690](../../tests/test_app_logic.py#L690) · **Area:** GUI: feature table · **Checked against:** hand-computed · **Strength:** 🟢 strong

The feature table sorts by each column and filters by text, and editing a feature via the dialog updates document and table in one undo step while an impossible location changes nothing.

| # | Input | Expected output |
|---|---|---|
| 1 | 4 features; sort column 0 asc/desc | ['Alpha','beta','mid','zeta'] / reversed (case-insensitive) |
| 2 | sort columns 1 (type), 2 (position), 3 (length), 4 (strand) | type order; position order; lengths ['21','31','51','151']; strand: reverse-strand '←' first |
| 3 | filter text 'cds','ALPHA','complement','term','a','nothing like this','' | respective name subsets; features count stays 4 (filter deletes nothing) |
| 4 | double-click a row | edit_feature_requested emits that feature index |
| 5 | FeatureDialog set name Omega, kind gene, location 120..260, accept | history +1; feature is gene 120..260; no 'Alpha'; table shows Omega/gene/120..260; undo record (['Omega'],['Alpha']) |
| 6 | undo then redo | undo: Alpha 100..250, steps back; redo: Omega 120..260 |
| 7 | edit location to '5..999999' (beyond sequence) | Omega stays 120..260; no extra step |

## `test_settings`

**Result:** ✅ passed in 0.07 s · [line 789](../../tests/test_app_logic.py#L789) · **Area:** GUI: settings dialog · **Checked against:** hand-computed · **Strength:** 🟢 strong

Tool settings persist in settings.ini and fall back to declared defaults on bad values; the dialog saves only on change, emits one change notice per tool, and Restore Defaults repopulates the form without saving until OK.

| # | Input | Expected output |
|---|---|---|
| 1 | registered test tool defaults ratio/count/flag/mode | [1.5, 3, True, 'fast'] |
| 2 | set ratio 2.5, count 7, flag False, mode slow; clear cache | ini has 'ratio=2.5'; values read back as [2.5, 7, False, 'slow'] typed |
| 3 | store count='abc', mode='gone' | fall back to defaults 3 and 'fast' |
| 4 | dialog accept with no change | heard == [] (nothing saved) |
| 5 | dialog change all four, accept | heard == [tool]; values [4.5,9,False,'slow']; persisted in file |
| 6 | Restore Defaults then (a) before OK, (b) accept | form shows defaults but stored still 4.5 and heard []; after accept heard==[tool] and values back to defaults |
| 7 | change primer-checks tm_low to 40.0, accept | heard == [CHECKS]; PP.limit('tm_low') == 40.0 |

## `test_self_test_panel_counts_passes_and_failures_and_shows_the_message`

**Result:** ✅ passed in 0.01 s · [line 853](../../tests/test_app_logic.py#L853) · **Area:** GUI: self-test panel · **Checked against:** hand-computed · **Strength:** 🟢 strong

The self-test panel runs the registered checks, counts passes and failures, and shows each result's mark, group, name, timing and the failure or exception message.

| # | Input | Expected output |
|---|---|---|
| 1 | 4 checks: pass, fail with message, fail with no message, crash (ValueError) | status 'Not run' then '1 passed, 3 failed'; 4 rows |
| 2 | marks and messages | marks [ok, bad, bad, bad]; messages '', 'one plus one is not three', 'assertion failed', 'ValueError: boom' |
| 3 | group/name columns and timing | groups [Area,Area,Other,Other]; name row 1 'fails'; every row has a digit time |
| 4 | two passing checks | status '2 passed, 0 failed' |
