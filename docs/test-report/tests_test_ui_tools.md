# `tests/test_ui_tools.py`

[← all test files](README.md) · source: [`tests/test_ui_tools.py`](../../tests/test_ui_tools.py)

## `test_spacer_finder_closes_with_escape`

**Result:** ✅ passed in 0.74 s · [line 48](../../tests/test_ui_tools.py#L48) · **Area:** GUI: Cas9 dialog (regression) · **Checked against:** regression value · **Strength:** 🟢 strong

Regression: an instance attribute self.done once shadowed QDialog.done(), so Escape/close raised TypeError and the dialog never closed. Verifies done is still callable and Escape rejects the dialog.

| # | Input | Expected output |
|---|---|---|
| 1 | Cas9Dialog, assert dlg.done callable; QTimer press Escape at 300ms; exec() | seen['closed'] True (not visible after Escape), result == QDialog.Rejected |

## `test_spacer_finder_columns_fit_their_content`

**Result:** ✅ passed in 0.40 s · [line 65](../../tests/test_ui_tools.py#L65) · **Area:** GUI: Cas9 dialog table columns · **Checked against:** property/invariant · **Strength:** 🟡 adequate

The Cas9 guide table's columns are wide enough for header+cell text but not spread far wider, and it hides its vertical header.

| # | Input | Expected output |
|---|---|---|
| 1 | Cas9Dialog on demo at 1800x600, process events; rowCount>0; vertical header hidden | for each column, computed need (max header/cell advance + padding) fits via _fits: width <= max(2.2*need, need+120) |

**Review notes:**

- _fits (line 42-45) only bounds the width from above; a column too narrow to show its content would not be caught here.

## `test_enter_in_a_field_presses_no_stray_button`

**Result:** ✅ passed in 0.24 s · [line 83](../../tests/test_ui_tools.py#L83) · **Area:** GUI: ORF/Digest dialog default-button (regression) · **Checked against:** regression value · **Strength:** 🟢 strong

Regression: with no button box the first push button became the Enter default, so Enter in a text field added CDS features (ORF dialog) or a lane (Digest). Verifies Enter in a field presses nothing.

| # | Input | Expected output |
|---|---|---|
| 1 | OrfDialog min_aa=30, select a new ORF row (add enabled), focus min_aa, Enter | feature count unchanged; clearing selection disables add button |
| 2 | DigestDialog, focus filter box, Enter | len(lanes)==1 (no lane added) |

**Review notes:**

- Line 92 relies on the demo file having an unannotated ORF ('new' non-empty is asserted, so not vacuous), but the exact ORF content is not pinned; acceptable since the behaviour under test is Enter handling, not ORF detection.

## `test_one_enter_is_one_ncbi_search`

**Result:** ✅ passed in 0.02 s · [line 111](../../tests/test_ui_tools.py#L111) · **Area:** GUI: NCBI dialog (regression) · **Checked against:** regression value · **Strength:** 🟢 strong

Pressing Enter in the NCBI/Genome dialog query box triggers exactly one search request, not zero or two.

| # | Input | Expected output |
|---|---|---|
| 1 | for NcbiDialog and GenomeDialog: stub _run to append; set query 'pUC19', focus, Enter | len(sent)==1 for each class |

## `test_ncbi_error_text_is_not_read_as_markup`

**Result:** ✅ passed in 0.00 s · [line 125](../../tests/test_ui_tools.py#L125) · **Area:** GUI: NCBI error escaping (security/display) · **Checked against:** hand-computed · **Strength:** 🟢 strong

An NCBI error reply containing HTML-like text is shown escaped, not rendered as rich-text markup.

| # | Input | Expected output |
|---|---|---|
| 1 | NcbiDialog._done(token, None, 'bad <b>reply</b>') | '&lt;b&gt;' in status text |

## `test_the_digest_plus_sits_right_after_the_last_lane`

**Result:** ✅ passed in 0.29 s · [line 134](../../tests/test_ui_tools.py#L134) · **Area:** GUI: Digest dialog + button placement · **Checked against:** GUI state · **Strength:** 🟡 adequate

The Digest dialog's '+' add button sits immediately right of the last tab (lane) and keeps its tooltip as lanes are added.

| # | Input | Expected output |
|---|---|---|
| 1 | DigestDialog at 1200x600; twice: check + button then click it | add.parentWidget() is tab bar, visible; 0 <= add.left - lastTab.right <= 12; tooltip == 'Add another digest (lane)'; after two clicks len(lanes)==3 |

## `test_primer_database_columns_fit_and_rows_have_no_wasted_indent`

**Result:** ✅ passed in 2.16 s · [line 164](../../tests/test_ui_tools.py#L164) · **Area:** GUI: primer database tab columns · **Checked against:** property/invariant · **Strength:** 🟡 adequate

The Primer Database tree sizes columns to their content, wastes no indent column (rootIsDecorated off), and the sequence column is never narrower than its longest sequence.

| # | Input | Expected output |
|---|---|---|
| 1 | open primers tab with 3 imported primers at 1800x850; items non-empty; not rootIsDecorated | for each column _fits(need); sectionSize(1) >= max sequence advance (not cut off) |

**Review notes:**

- add_primers uses the shared global primer_db.db() (line 155) without a private-db fixture, so these test-only names ('tt_short' etc.) persist in the module's shared data dir and could interact with other tests in this file; _fits only bounds width above.

## `test_primer_database_delete_key_and_buttons`

**Result:** ✅ passed in 1.77 s · [line 179](../../tests/test_ui_tools.py#L179) · **Area:** GUI: primer database tab actions · **Checked against:** GUI state · **Strength:** 🟡 adequate

The Primer Database tab enables edit/delete only with a selection, Enter does not open the edit dialog, and the Delete key asks for confirmation then deletes the selected primer.

| # | Input | Expected output |
|---|---|---|
| 1 | clear selection | edit_button and delete_button disabled |
| 2 | select 'tt_x' | both buttons enabled |
| 3 | focus tree, press Enter (PrimerDialog.exec stubbed) | no dialog opened |
| 4 | press Delete (QMessageBox.question stubbed Yes) | question asked; db().find('tt_x') is None |

**Review notes:**

- Uses QTest.qWaitForWindowActive with a 2000ms timeout (line 198) and window activation that the comment notes is unreliable on Windows; the Delete-key path needs focus in the active window, so on a machine that won't activate the window the final assertion could behave unpredictably (best-effort by design).

## `test_delete_questions_default_to_no`

**Result:** ✅ passed in 0.00 s · [line 208](../../tests/test_ui_tools.py#L208) · **Area:** GUI: destructive-question default button · **Checked against:** hand-computed · **Strength:** 🟢 strong

A destructive delete confirmation defaults to No, so Enter or a stray click on the default keeps the primers.

| # | Input | Expected output |
|---|---|---|
| 1 | stub QMessageBox.question to return its last arg (default button); PrimersDatabase().delete(None, ['anything']) | question called; default button arg == QMessageBox.No |

## `test_a_running_mapping_is_not_hidden_by_escape`

**Result:** ✅ passed in 0.59 s · [line 233](../../tests/test_ui_tools.py#L233) · **Area:** GUI: run dialog Escape handling (regression) · **Checked against:** regression value · **Strength:** 🟡 adequate

Regression: Escape went to reject() (skipping closeEvent) so the window hid while bwa/minimap2 kept running. Verifies Escape with 'No' keeps running, with 'Yes' stops and closes.

| # | Input | Expected output |
|---|---|---|
| 1 | RunDialog running python sleep 30; qWait 300; Escape with question->No | still visible and running |
| 2 | Escape with question->Yes | not visible and not running |

**Review notes:**

- Spawns a real subprocess and uses QTest.qWait(300) to let it start (line 239); a slow start could leave running() False at the first assertion. Flaky risk under load.

## `test_a_failed_step_says_its_exit_code`

**Result:** ✅ passed in 0.15 s · [line 249](../../tests/test_ui_tools.py#L249) · **Area:** GUI: run dialog failure reporting · **Checked against:** hand-computed · **Strength:** 🟡 adequate

A run step that exits non-zero reports 'Failed (exit code N).' with the real exit code.

| # | Input | Expected output |
|---|---|---|
| 1 | RunDialog running python 'exit 3'; wait up to ~5s for not running | status text == 'Failed (exit code 3).' |

**Review notes:**

- Real subprocess + poll loop (up to 100*50ms); if the process hasn't finished, running() stays True and the final assertion may read a stale status. Timing-dependent.
