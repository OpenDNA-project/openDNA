# `tests/test_history.py`

[← all test files](README.md) · source: [`tests/test_history.py`](../../tests/test_history.py)

## `test_history_steps`

**Result:** ✅ passed in 8.20 s · [line 71](../../tests/test_history.py#L71) · **Area:** History: annotations and primers · **Checked against:** hand-computed · **Strength:** 🟢 strong

Annotation and primer changes are undoable History steps that interleave correctly with base edits, a primer-only change does not mark the document modified, importing records one step, a primer changed elsewhere is left alone on undo, and the History panel navigates the steps.

| # | Input | Expected output |
|---|---|---|
| 1 | add/edit(location)/edit(label)/delete feature | texts == ['Added feature alpha','Changed feature alpha','Added feature beta; Deleted feature alpha','Deleted feature beta'] |
| 2 | saved doc; add_primer pFwd, type 'A' at 11, add_primer pRev | primer add leaves doc unmodified; texts ['Added primer pFwd','Inserted 1 bp at 11: A','Added primer pRev']; undo removes pRev then base then pFwd in order; redo restores in order; a new edit after undo drops redo incl primer steps |
| 3 | edit/rename/delete primers p1->p1b etc. | 'Renamed primer p1 to p1b'; 'Edited primer p1b'; renaming to existing name refused (not a step, len==4); delete->'Deleted primers p1b, p2'; undo restores original row incl date; undo_steps(3)/(-3) work |
| 4 | import 3 primers a1,a2,a3; re-import a1; update a2 elsewhere then undo/redo | 'Imported 3 primers: a1, a2, a3'; re-import exists==1 nothing recorded; undo leaves a2 (changed since) with error 'changed since'; redo restores a1,a3 keeping a2; undoing a delete whose name was retaken does not overwrite ('used now' error) |
| 5 | HistoryPanel: add pA, type G, delete pA; setCurrentRow(1) go; last row go | items[1:4] == ['Added primer pA','Inserted 1 bp at 1: G','Deleted primer pA']; go to row1 undoes all after; go to last redoes all |
| 6 | PrimerDialog add 'viaDialog' then edit note | texts 'Added primer viaDialog' then 'Edited primer viaDialog'; note 'a note'; two undos clear names |

**Review notes:**

- Very long merged test (~100 lines, 6 sections); it does re-seed the primer db between sections via fresh_db and asserts exact strings, so a failure is diagnosable, but the single function makes the first failing section hide later ones.

## `test_primer_menu`

**Result:** ✅ passed in 2.57 s · [line 182](../../tests/test_history.py#L182) · **Area:** GUI: primer right-click menu · **Checked against:** hand-computed · **Strength:** 🟢 strong

The primer menu (Edit/Delete/Import/Hide) enables actions appropriately for database vs temporary primers, and right-clicking a primer in the sequence view, the linear/circular map, and the primer panel opens the correct menu and selects the primer.

| # | Input | Expected output |
|---|---|---|
| 1 | primer_actions for db primer 'menuP' | Import disabled (already in file); Edit/Delete/Hide enabled; Edit triggers edit request; Delete deletes with 'Undo brings it back' + step 'Deleted primer menuP'; undo restores |
| 2 | primer_actions for temporary 'temp1' | Edit/Delete disabled; Hide enabled |
| 3 | sequence view right-click on primer's 4th base | primer_menu called with site.key; real menu's first action starts 'Edit'; separator label == ['viewP']; current_primer set |
| 4 | map view linear & circular: hover/click primer mark | _primer_at returns key at the mark, None 40px away; cursor PointingHand on hover; click sets current_primer |
| 5 | primer panel right-click row | separator label ['rowP']; Delete action deletes (names empty) |
