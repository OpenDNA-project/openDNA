# `tests/test_db_tabs.py`

[← all test files](README.md) · source: [`tests/test_db_tabs.py`](../../tests/test_db_tabs.py)

## `test_features_database_tab`

**Result:** ✅ passed in 0.21 s · [line 70](../../tests/test_db_tabs.py#L70) · **Area:** GUI: Features Database tab · **Checked against:** hand-computed · **Strength:** 🟢 strong

The Features Database opens as a tab on a document made of the database, shows one row per feature with source+file, displays a clicked feature alone, filters by Source, survives library changes, deletes only Local features, and follows the window across file opens/reopens.

| # | Input | Expected output |
|---|---|---|
| 1 | trigger 'Features Database…' | DatabaseTab tab titled 'Features Database'; its editor != main editor; hint shown, no sequence; columns ['Name','Type','Source','Length'] |
| 2 | rows | 'my promoter' source 'Local · pMine.gb'; 'AmpR' source 'Cloning vectors (NCBI records) · pBR322 (J01749)' |
| 3 | select 'my promoter' row | ReadOnlyView; sequence 'TTGACAATTAATCATCGGCTCG', dna_locked; feature (promoter,'my promoter',[(0,22)]); caption names it + file |
| 4 | filter Source to hide 'Cloning vectors' | value counts {Cloning vectors:2, LOCAL:1}; shown rows ['my promoter']; header has MARK |
| 5 | add 'rbs' part, reprocess | shown rows ['my promoter','rbs']; on-view feature still 'TTGACAA...' |
| 6 | delete 'my promoter' and 'AmpR' (AmpR not Local) | remaining parts {AmpR, lacZ alpha, rbs} (AmpR survives); sequence None |
| 7 | help dialog; open other file; close & reopen tab | help status/title 'Features Database'; tab follows window; reopened tab is the same object |

## `test_primer_database_tab`

**Result:** ✅ passed in 0.34 s · [line 122](../../tests/test_db_tabs.py#L122) · **Area:** GUI: Primer Database tab · **Checked against:** hand-computed · **Strength:** 🟢 strong

The Primer Database tab lists every primer (even ones binding nowhere) with note/added columns, shows a selected primer alone, sorts Length numerically with a stable sort through refresh, deletes a primer, and keeps the selection.

| # | Input | Expected output |
|---|---|---|
| 1 | add 3 primers, open tab | title 'Primer Database', PrimerPanel; tree items == all db names; tab_fwd note 'for the tab' and non-empty added; columns end ['Added','Note'] |
| 2 | select tab_fwd | entry sequence 'ACGTTGCAGGCTAGCATCGA' |
| 3 | click Length header once then again | lengths ascending then descending (numeric) |
| 4 | delete 'tab_amb' | not in names; db().find None; sort still descending after rebuild; current item still 'tab_fwd' |
| 5 | help dialog | status.title 'Primer Database' |

**Review notes:**

- Uses the shared global primer_db.db() with no private-db fixture (lines 123-126), so 'tab_fwd/tab_amb/tab_long' persist in the module's data dir; within-module ordering matters.

## `test_databases_menu`

**Result:** ✅ passed in 0.05 s · [line 158](../../tests/test_db_tabs.py#L158) · **Area:** GUI: Databases menu contents · **Checked against:** hand-computed · **Strength:** 🟡 adequate

The Databases menu lists exactly Features Database, Primer Database, Database Status (no enzyme or reference-genome entries).

| # | Input | Expected output |
|---|---|---|
| 1 | non-separator action texts of Databases menu | == ['&Features Database…','&Primer Database…','Database &Status'] |

## `test_value_filter_on_a_table_and_a_tree`

**Result:** ✅ passed in 0.01 s · [line 164](../../tests/test_db_tabs.py#L164) · **Area:** GUI: value filter widget · **Checked against:** hand-computed · **Strength:** 🟢 strong

The value-filter on a QTableWidget keeps its filter when the table is refilled and shows newly-seen values, the ValuePicker search+tick path hides a value, and the same works on a QTreeWidget's top rows.

| # | Input | Expected output |
|---|---|---|
| 1 | fill rows strands →,←,→,(empty); vf.values(1) | [('←',1),('→',2),(EMPTY,1)] |
| 2 | hide {←, EMPTY} | shown_rows == [0,2]; header 'Strand'+MARK |
| 3 | refill with ←,→,— ; process | shown_rows == [1,2] (still filtered, new '—' shown) |
| 4 | clear(); then ValuePicker search '←', _tick(False), _ok() | shown_rows == [0,1,2]; header back to 'Strand'; vf.hidden == {'Strand':{'←'}} |
| 5 | tree with p1/p2/p3 (match '1 site'/'no site'/'1 site'); hide {'no site'} | shown_rows == [0,2]; header(1) 'Match'+MARK |

## `test_reference_box_lists_the_last_three`

**Result:** ✅ passed in 0.46 s · [line 202](../../tests/test_db_tabs.py#L202) · **Area:** GUI: reference genome picker · **Checked against:** hand-computed · **Strength:** 🟢 strong

The toolbar reference-genome combobox remembers the last three browsed genomes (most recent first) plus Off/Add entries, and selecting an entry or Off updates the current reference.

| # | Input | Expected output |
|---|---|---|
| 1 | browse genome0..3 via stubbed getOpenFileName | items == ['genome3.fa','genome2.fa','genome1.fa', OFF, '', ADD]; index 0; current()==genome3 |
| 2 | activate index 2 | current()==genome1 file; text 'genome1.fa' |
| 3 | activate OFF | current() is None; text == OFF |

**Review notes:**

- Mutates the shared primer_db settings (primers/genome etc.) but restores them in a finally block (lines 223-227), so state is cleaned up.

## `test_primer_sites_grouped_with_mismatches_in_red`

**Result:** ✅ passed in 0.44 s · [line 230](../../tests/test_db_tabs.py#L230) · **Area:** GUI: primer panel site grouping/mismatches · **Checked against:** hand-computed · **Strength:** 🟢 strong

A primer's binding sites are grouped (this sequence, then reference genome once opened); mismatching primer bases are marked for red drawing; a genome hit on the - strand is reported with coordinates and 'whole primer, exact'.

| # | Input | Expected output |
|---|---|---|
| 1 | mismatch_indices('ACGTNR','ACCTGA') | == [2] (N,R pair ambiguously) |
| 2 | seeded 600bp seq; primer = seq[100:124] with 1 mismatch at index5; genome FASTA with rc(primer) on - strand; add 'red_test' | 'On this sequence (1)' group, 1 child, genome only once opened; site MISMATCHES==[5], text '1 mismatch', primer_item is top |
| 3 | expand top -> background genome searched on worker | genome group starts 'On the background genome, chrT'; hit text starts 'chrT:151–', ends '←', MISMATCHES==[], text 'whole primer, exact' |
| 4 | select the genome hit | no error (nothing to select on this sequence) |

**Review notes:**

- Opening the genome search uses a wall-clock busy loop of up to 200*(processEvents + 20ms sleep) (lines 262-268); on a slow machine the genome group may not appear in time and 'genome = top.child(1)' (line 269) would then IndexError rather than fail with a clear message.
