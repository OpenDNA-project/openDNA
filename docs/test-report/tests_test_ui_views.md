# `tests/test_ui_views.py`

[← all test files](README.md) · source: [`tests/test_ui_views.py`](../../tests/test_ui_views.py)

## `test_ctrl_home_end_and_page_keys`

**Result:** ✅ passed in 0.02 s · [line 69](../../tests/test_ui_views.py#L69) · **Area:** GUI: sequence view keyboard navigation · **Checked against:** hand-computed · **Strength:** 🟢 strong

Ctrl+Home/End, End, PageUp/PageDown and their Shift variants move the cursor/selection by the right amount and scroll the view, so keyboard navigation of a 6,000 bp sequence works as expected.

| # | Input | Expected output |
|---|---|---|
| 1 | 6000 bp view, cursor at 1234, press Ctrl+End | start==end==length and vertical scrollbar at maximum |
| 2 | Ctrl+Shift+Home then Ctrl+Home | (start,end)==(0,length) then start==end==0 and scroll==0 |
| 3 | End alone | start == bases_per_line() (end of the row) |
| 4 | cursor 0, anchor 0, PageDown | moved > 2*bases_per_line, moved % bases_per_line == 0, scroll>0 |
| 5 | Shift+PageDown; then two PageUp | Shift extends (start==moved, end>moved); PageUp collapses to start<moved |

## `test_the_cursor_is_faint_in_a_view_without_the_keyboard`

**Result:** ✅ passed in 0.02 s · [line 94](../../tests/test_ui_views.py#L94) · **Area:** GUI: sequence view cursor rendering · **Checked against:** GUI state · **Strength:** 🟡 adequate

The text cursor is drawn solid in the focused view and fainter when focus leaves, so a user with several panes can tell where typing goes.

| # | Input | Expected output |
|---|---|---|
| 1 | 700x300 view, cursor at 5; grab a 5px-wide strip around the cursor x over one row; measure ink (sum \|pixel-background\| over RGB) with focus and without focus | hasFocus True then False; faint < solid |

**Review notes:**

- Relies on offscreen pixel rendering; faint<solid is a relative inequality (line 123) with no absolute magnitude, so a near-invisible 'solid' cursor could still pass as long as faint is even smaller.

## `test_wheel_zoom_steps_once_per_notch`

**Result:** ✅ passed in 0.03 s · [line 129](../../tests/test_ui_views.py#L129) · **Area:** GUI: sequence/map view wheel zoom · **Checked against:** hand-computed · **Strength:** 🟢 strong

Ctrl+wheel zooms one font-size step per full notch (two touchpad half-notches), ignores sideways scrolls, and the map view obeys the same rule.

| # | Input | Expected output |
|---|---|---|
| 1 | SequenceView, Ctrl+wheel dy=60 (half notch) twice | size unchanged after first, size+1 after second |
| 2 | Ctrl+wheel sideways dx=120, dy=0 | size unchanged (no zoom) |
| 3 | Ctrl+wheel dy=-240 (two notches) | size == size+1-2 == size-1 |
| 4 | MapView wheel dy=60 twice at centre, then sideways dx=0 | view None after first, not None after second; span unchanged by sideways |

## `test_a_map_click_that_wobbles_is_not_a_selection`

**Result:** ✅ passed in 0.03 s · [line 157](../../tests/test_ui_views.py#L157) · **Area:** GUI: map view click vs drag · **Checked against:** property/invariant · **Strength:** 🟡 adequate

On the map, a 1px wobble after pressing stays a click (cursor move only), but a real drag beyond startDragDistance makes a selection, so a plain click does not accidentally select a range.

| # | Input | Expected output |
|---|---|---|
| 1 | press at linear midpoint; move +1px | start==end==cursor (no selection) |
| 2 | move +(startDragDistance+30)px | end > start (selection made) |

## `test_dragging_past_the_bottom_edge_scrolls_and_extends_the_selection`

**Result:** ✅ passed in 0.02 s · [line 174](../../tests/test_ui_views.py#L174) · **Area:** GUI: sequence view autoscroll drag · **Checked against:** GUI state · **Strength:** 🟡 adequate

Dragging below the viewport with the button held starts an autoscroll timer that both scrolls down and extends the selection, so a selection can span past the visible rows.

| # | Input | Expected output |
|---|---|---|
| 1 | press inside row 0, move to y=viewport.height()+40 | scroll_timer.isActive() |
| 2 | call _autoscroll() 5 times | scrollbar value increased AND ed.end increased |
| 3 | release button | scroll_timer not active |

## `test_right_click_on_empty_space_offers_cut_copy_paste_select_all`

**Result:** ✅ passed in 0.03 s · [line 192](../../tests/test_ui_views.py#L192) · **Area:** GUI: context menu (sequence/map/read-only views) · **Checked against:** GUI state · **Strength:** 🟢 strong

Right-click / Menu key on empty space offers Cut/Copy/Paste/Select All with correct enabled states, greys Paste when DNA is locked, is refused by a read-only view, and the map shows the same menu.

| # | Input | Expected output |
|---|---|---|
| 1 | clipboard 'ACGT', cursor 10, right-click empty space (QMenu stubbed by Recorder) | items == {Cut:False, Copy:False, Paste:True, Select All:True} |
| 2 | trigger 'Select All' | (start,end)==(0,length) |
| 3 | dna_locked=True, Menu-key contextMenuEvent | 'Paste (DNA locked)' disabled, 'Select All' enabled |
| 4 | assembly_panel.ReadOnlyView Menu-key event | no new menu shown (count unchanged) |
| 5 | MapView right-click empty | one new menu shown containing 'Select All' |

## `test_the_alignment_view_scrolls_instead_of_growing_past_the_screen`

**Result:** ✅ passed in 0.42 s · [line 229](../../tests/test_ui_views.py#L229) · **Area:** GUI: genome/alignment view layout · **Checked against:** GUI state · **Strength:** 🟢 strong

Very tall tracks (many BAM files) scroll inside the genome view rather than making the view or window taller than the screen; a long info line does not widen the window; the canvas never steals the keyboard.

| # | Input | Expected output |
|---|---|---|
| 1 | load demo pFG051_Chlor.gb, show genome view, force GC%/Features track height=3000, refresh | canvas minimumHeight>6000; gv.minimumSizeHint and win.minimumSizeHint height < screen; tracks scrollbar maximum>0; gv.height() <= screen*0.6+80 |
| 2 | set info label text to 'read_1234 '*300 | window width unchanged; text starts 'read_1234'; tooltip == text |
| 3 | focus policies | canvas and tracks_area focusPolicy == NoFocus |

**Review notes:**

- Depends on the real screen's availableGeometry height; on a very small/large virtual screen the 0.6*screen+80 bound (line 250) could behave differently, though offscreen default is typically fine.

## `test_checkup_tables_keep_their_row_and_scroll_and_fit_their_columns`

**Result:** ✅ passed in 0.09 s · [line 265](../../tests/test_ui_views.py#L265) · **Area:** GUI: checkup dialog table refresh · **Checked against:** GUI state · **Strength:** 🟡 adequate

Refreshing the Checkup table (after an edit) keeps the selected row and the scroll position, and columns are sized to their content (details column not spread to 1000px).

| # | Input | Expected output |
|---|---|---|
| 1 | open CheckupDialog on demo file, select cell (5,0), scroll to bottom, record name+scroll, call refresh() | currentRow()==5, item(5,0) text unchanged, scroll unchanged |
| 2 | for each non-last column | sectionSize(c) >= min(need,260)-2 where need=max(sizeHintForColumn, header sectionSizeHint) |

**Review notes:**

- Column-fit check (lines 286-288) only asserts a lower bound (not too narrow); the docstring 'details column is not spread over 1000 px' (upper bound) is not asserted, so a column spread far too wide would pass.

## `test_protein_view_keeps_its_place_through_an_edit`

**Result:** ✅ passed in 0.04 s · [line 292](../../tests/test_ui_views.py#L292) · **Area:** GUI: protein view state preservation · **Checked against:** round-trip/self-consistency · **Strength:** 🟡 adequate

An edit elsewhere redraws the protein view showing the same protein with scroll position and text selection preserved; title is one line with full text in tooltip; empty feature list is hidden.

| # | Input | Expected output |
|---|---|---|
| 1 | seeded (Random(5)) 600-codon CDS at 11..length, select CDS, refresh | features list not visible; title tooltip set and 'a CDS with a name' in title text |
| 2 | scroll to max//2, select text positions 200-205, record scroll+selection; add a 'note' qualifier (an edit); refresh | scrollbar value and selected text unchanged |

## `test_export_is_light_on_a_dark_program_and_leaves_the_controls_out`

**Result:** ✅ passed in 6.93 s · [line 329](../../tests/test_ui_views.py#L329) · **Area:** GUI: image export · **Checked against:** hand-computed · **Strength:** 🟢 strong

Image export of a sequence view, while the app palette is dark, produces a light-background image when light=True (dark when False), excludes the zoom bar and scroll bar, and restores the live view afterwards.

| # | Input | Expected output |
|---|---|---|
| 1 | set dark palette; assert Base lightness<60; export light.png and dark.png | light.png pixel(2,2) lightness>200, dark.png pixel(2,2) lightness<60 |
| 2 | check the zoom-bar x and the right-edge scroll-bar x columns (y 2..130) of light.png | every pixel == '#ffffff' (no control strip in export) |
| 3 | after export | zoom_bar and scrollbar parent still visible; Base lightness still <60 (restored) |

## `test_the_export_dialog_does_not_grow_when_it_warns`

**Result:** ✅ passed in 1.02 s · [line 360](../../tests/test_ui_views.py#L360) · **Area:** GUI: image export dialog layout · **Checked against:** GUI state · **Strength:** 🟡 adequate

The export dialog reserves room for its 'Too large' warning so the dialog does not change height when the warning appears, and disables OK; light background defaults on.

| # | Input | Expected output |
|---|---|---|
| 1 | ImageExportDialog on 1400x800 view; record sizeHint height; set dpi=2400, scale=1000 | 'Too large' in warning text, OK disabled, sizeHint height unchanged, light checkbox checked |

## `test_a_long_local_structure_command_shows_progress_and_can_be_cancelled`

**Result:** ✅ passed in 0.54 s · [line 379](../../tests/test_ui_views.py#L379) · **Area:** GUI: structure dialog local command · **Checked against:** GUI state · **Strength:** 🟡 adequate

A long local structure-prediction command shows a busy/cancel UI, can be cancelled (killing the subprocess), reports 'Cancelled', and a late worker reply is not shown.

| # | Input | Expected output |
|---|---|---|
| 1 | StructureDialog, source 'Local command' running python sleeping 60s, start(); poll up to 10s for job.proc | proc running (poll() is None), busy and cancel button visible |
| 2 | cancel(); wait for proc | proc.poll() not None (stopped); get button enabled, busy hidden, 'Cancelled' in status |
| 3 | process events for 0.5s | status still 'Cancelled' (late reply ignored) |

**Review notes:**

- Spawns a real subprocess (sys.executable sleeping 60s) and uses wall-clock wait loops with 10s timeouts and time.sleep (lines 388-398); on a slow/loaded CI machine the proc may not appear within 10s, making it flaky and slow.
