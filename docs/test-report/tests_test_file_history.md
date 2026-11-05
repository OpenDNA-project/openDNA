# `tests/test_file_history.py`

[← all test files](README.md) · source: [`tests/test_file_history.py`](../../tests/test_file_history.py)

## `test_log`

**Result:** ✅ passed · [line 26](../../tests/test_file_history.py#L26) · **Area:** Multi-record undo log · **Checked against:** hand-computed · **Strength:** 🟢 strong

The file-wide undo log orders edits across records and follows undo/redo (a new edit drops the redo branch); steps() interleaves each record's own summaries in file order.

| # | Input | Expected output |
|---|---|---|
| 1 | observe A0,B0,A1,B1,A2 (undo counts) | undo_log ['A','B','A'], next_undo 'A' |
| 2 | observe A(1,1) undo in A; B(0,1) undo in B | undo ['A','B'], redo ['A'], next_redo 'A', next_undo 'B'; then undo ['A'], redo ['A','B'], next_redo 'B' |
| 3 | observe B(1,0) redo; A(2,0) new edit | undo ['A','B'], redo ['A']; then undo ['A','B','A'], redo [] |
| 4 | reset() | both logs empty, next_undo None |
| 5 | undo_log ['A','B','A'], redo_log ['B','A']; steps({A:(['a1','a2'],['a3']), B:(['b1'],['b2'])}) | undo [('A','a1'),('B','b1'),('A','a2')]; redo [('A','a3'),('B','b2')] |

## `test_file_as_one_document`

**Result:** ✅ passed in 2.19 s · [line 80](../../tests/test_file_history.py#L80) · **Area:** Multi-record undo / display / Save As · **Checked against:** hand-computed · **Strength:** 🟢 strong

A multi-record file behaves as one document: undo/redo walk back across records in edit order, the History panel lists the whole file and can jump back, display settings carry over between records, and Save As writes every edited record.

| # | Input | Expected output |
|---|---|---|
| 1 | two-record file: 'AAAA' at 10 in first, 'CCCC' at 20 in second, 'GGGG' at 30 in first | undo_log [first, second, first] |
| 2 | undo x3 | first without GGGG; second == original; first == original; undo_log [], 3 redo |
| 3 | redo x3 | first with AAAA at 10; second with CCCC at 20; first with GGGG at 30; logs restored |
| 4 | History panel after 2 edits; go to row 1 | first-record rows before second-record rows; first shown, undo_log [], second's history empty |
| 5 | set_display(enzymes, hidden {promoter}, ...), letter size 14, DNA locked; switch record | new editor carries all settings; size 14 |
| 6 | edits in both records, Save As copy.gb | 2 records with both edits at 10 and 20; original untouched; path copy.gb, nothing modified |

**Review notes:**

- Line 98: after undoing edit 3 only `seq[30:34] != 'GGGG'` is checked, not the exact expected sequence (base[:10]+'AAAA'+base[10:]).
- Lines 85, 112, 130, 141: four MainWindows with unsaved edits are opened and never closed.

## `test_theme_switch_repaints_a_multi_record_window_and_keeps_the_file_history`

**Result:** ✅ passed in 0.74 s · [line 161](../../tests/test_file_history.py#L161) · **Area:** GUI: theme switch with multi-record file · **Checked against:** GUI state · **Strength:** 🟡 adequate

Switching Light/Dark theme repaints a multi-record window without losing the History panel rows or the file-wide undo log.

| # | Input | Expected output |
|---|---|---|
| 1 | edits in first (AAAA at 10) and second (CCCC at 20); HistoryPanel; theme.refresh_views(); win.grab() | same panel row count; undo_log == [first, second] |

**Review notes:**

- No docstring; the repaint itself is smoke only (no colour or widget state checked after refresh_views).
