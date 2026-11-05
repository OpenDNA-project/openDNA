# `tests/test_primer_origin.py`

[← all test files](README.md) · source: [`tests/test_primer_origin.py`](../../tests/test_primer_origin.py)

## `test_a_primer_is_drawn_by_where_it_comes_from`

**Result:** ✅ passed in 6.26 s · [line 58](../../tests/test_primer_origin.py#L58) · **Area:** Primers: origin and colour · **Checked against:** GUI state · **Strength:** 🟢 strong

A primer's origin (database only, file, both, designed/temporary) sets its flags and colour, the views paint each origin in its colour, and the display filter hides file or database primers.

| # | Input | Expected output |
|---|---|---|
| 1 | circular 600 bp (seed 21); file has 'mine' 151..174, 'both' 251..274; DB has 'dbonly' SEQ[50:74], 'both'; temp 'temp' SEQ[350:374] | sites {dbonly, mine, both, temp} |
| 2 | flags (from_db, from_file, provisional) | dbonly (T,F,T); mine (F,T,F); both (T,T,F); temp temporary and provisional |
| 3 | primer_color of each | dbonly == temp == TEMP_COLOR; mine == both == PRIMER_COLOR |
| 4 | SequenceView 1000x700, count pixels within 12 RGB of each colour | both colours > 0; with file_primers off: PRIMER_COLOR 0, TEMP_COLOR > 0 |
| 5 | file_primers on, db_primers off | drawn_primers {mine, both, temp} |
| 6 | MapView grab | no crash |
| 7 | DisplayPanel | file/db swatches are PrimerSwatch with PRIMER_COLOR / TEMP_COLOR; tooltips non-empty, db swatch tooltip == checkbox tooltip; unchecking db_primers leaves {mine, both, temp} |

**Review notes:**

- map view only grabbed (line 93); comment claims same colours but nothing is checked

## `test_import_to_this_sequence_menu_entries_and_one_undo_step`

**Result:** ✅ passed in 0.74 s · [line 105](../../tests/test_primer_origin.py#L105) · **Area:** Primers: import to sequence · **Checked against:** GUI state · **Strength:** 🟢 strong

The primer context menu offers the right actions per origin, and Import to This Sequence turns a database primer into a file annotation in one undoable step.

| # | Input | Expected output |
|---|---|---|
| 1 | primer_actions menu for dbonly, temp, mine, both | Import enabled for dbonly and temp; Edit enabled for dbonly not temp; Import disabled for mine and both |
| 2 | primers_to_import(['dbonly','temp','both']) | [('dbonly', DBONLY), ('temp', TEMP)] |
| 3 | trigger Import on dbonly | history +1; last text 'Imported primer dbonly to this sequence'; site from_file, not provisional, PRIMER_COLOR; hover says 'in this file'; Import now disabled |
| 4 | undo, redo | undo: provisional again, hover 'not in this file'; redo: not provisional |

## `test_imported_and_made_primers_survive_save_and_reload_as_on_sequence_primers`

**Result:** ✅ passed in 1.79 s · [line 139](../../tests/test_primer_origin.py#L139) · **Area:** Primers: save and reload · **Checked against:** round-trip/self-consistency · **Strength:** 🟡 adequate

Imported and temporary primers are saved as primer_bind features and reload as regular file primers even if removed from the database; a primer made in the editor belongs to the file.

| # | Input | Expected output |
|---|---|---|
| 1 | import_primers [dbonly, temp]; save o.gb | temp not provisional; saved primer_bind labels {mine, both, dbonly, temp} |
| 2 | remove dbonly from DB; reopen o.gb | dbonly, temp, mine, both all from_file and not provisional |
| 3 | fresh editor on a new unsaved doc (p.gb) | no 'temp' site that is non-provisional |
| 4 | add_primer 'madehere' SEQ[450:474]; save m.gb | not provisional; 'madehere' in saved primer_bind labels |

**Review notes:**

- line 153: the check passes if 'temp' is missing altogether; it should assert temp exists and is provisional
