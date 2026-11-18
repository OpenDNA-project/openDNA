# `tests/test_file_tabs.py`

[← all test files](README.md) · source: [`tests/test_file_tabs.py`](../../tests/test_file_tabs.py)

## `test_file_set_keeps_order_dedupes_and_picks_the_next_page_on_close`

**Result:** ✅ passed in 0.01 s · [line 27](../../tests/test_file_tabs.py#L27) · **Area:** File tabs · **Checked against:** hand-computed · **Strength:** 🟡 adequate

**Why it matters:** Without it, opening a file twice could give two editors on one file, closing a file tab could leave a closed file shown or nothing shown, and a window could close without asking about a file with unsaved changes.

The bookkeeping of the files of one window (FileSet), with stand-in pages: tab order and dragging, stepping through the tabs, the shown page, finding a file by path whatever its spelling, which pages are unsaved, and which page is shown after a close.

| # | Input | Expected output |
|---|---|---|
| 1 | add a.gb, b.gb (unsaved), c.gb, one untitled unsaved page | the untitled one is shown; unsaved() = [b, untitled] |
| 2 | find b.gb; a missing file; sub/../a.gb; C.GB in capitals | b's page; None; a's page; c's page on Windows only (case-insensitive paths there) |
| 3 | drag a.gb to the third place; step(+1) and step(-1) from the shown last page | order b, c, a, untitled; next wraps round to b; previous is a |
| 4 | close the shown last page | the page to its left is shown |
| 5 | close the shown page b; then close a (not shown); then close c | c shown; c still shown; None and no pages |

**Review notes:**

- Uses stand-in pages, not real windows: that the host really shows the page FileSet says is shown is only checked by hand (no UI tests).

## `test_two_open_files_share_no_editor_state`

**Result:** ✅ passed in 0.00 s · [line 51](../../tests/test_file_tabs.py#L51) · **Area:** File tabs · **Checked against:** hand-computed · **Strength:** 🟡 adequate

**Why it matters:** Without it, an edit or cursor move in one open file could silently land in another open file.

Two Editors on two documents: typing into one changes only that one (its length, version and cursor), the other keeps its bases, version and cursor.

| # | Input | Expected output |
|---|---|---|
| 1 | Editor one (ACGT x50) and Editor two (TTTT x50); cursor at 0 in one, type_bases('GGGG') | one: 204 bp, version up, cursor 4; two: 200 bp, same version, cursor 0, still TTTT x50 |

**Review notes:**

- Covers the Editor objects only; the widgets of a page (search box, range box, docks) are separate objects per page by construction and are not exercised here (no UI tests).

## `test_a_second_launch_hands_its_files_to_the_running_program`

**Result:** ✅ passed in 0.27 s · [line 63](../../tests/test_file_tabs.py#L63) · **Area:** File tabs · **Checked against:** hand-computed · **Strength:** 🟡 adequate

**Why it matters:** Without it, double-clicking a file while openDNA is open would start a second copy of the program instead of opening the file in the one running, or the file names could be lost or garbled on the way.

The hand-over between launches: with nothing running the launch is told so (and goes on to start); with a listener running, a second process sends three paths (spaces, accents) and gets the answer, and the listener receives them whole and in order.

| # | Input | Expected output |
|---|---|---|
| 1 | send_to_running with no listener | False |
| 2 | a listener here; a second process sends ['<tmp>/a.gb', '<tmp>/b c.gb', 'ünï.gb'] | the sender exits 0 (it got the answer); the listener's files signal gave exactly that list |

**Review notes:**

- Does not open the files in a window (that is main.open_handed_over, checked by hand: no UI tests), and not a launch while the running program is busy or hung.
