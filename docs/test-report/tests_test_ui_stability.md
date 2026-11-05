# `tests/test_ui_stability.py`

[← all test files](README.md) · source: [`tests/test_ui_stability.py`](../../tests/test_ui_stability.py)

## `test_no_button_changes_size_when_pressed_toggled_or_hovered`

**Result:** ✅ passed in 0.41 s · [line 68](../../tests/test_ui_stability.py#L68) · **Area:** GUI: layout stability of buttons · **Checked against:** GUI state · **Strength:** 🟢 strong

No visible button may change size (and push neighbours) when hovered, pressed or toggled; protects against the old checked-state style sheet that made Lock/Filters/History/Circular jump by 5 px.

| # | Input | Expected output |
|---|---|---|
| 1 | Main window on demo/pFG051_Chlor.gb at 1400x850: every visible button hovered (WA_UnderMouse), setDown, toggled | _unstable_buttons == [] (sizeHint and sibling geometries unchanged) |
| 2 | Display dock shown | [] |
| 3 | History dock shown | [] |

## `test_feature_table_columns_fit_their_content_in_a_tab_of_its_own`

**Result:** ✅ passed in 0.28 s · [line 79](../../tests/test_ui_stability.py#L79) · **Area:** GUI: feature table column widths · **Checked against:** property/invariant · **Strength:** 🟡 adequate

A feature table alone in a wide tab must not stretch any column far beyond its content (Strand once took 800 px) and must not cut off names.

| # | Input | Expected output |
|---|---|---|
| 1 | New tab with only Features, pFG051_Chlor.gb | each column width <= max(2.2*need, need+120), need = widest header/cell text + padding |
| 2 | Name column | width >= widest name text |

**Review notes:**

- L95: tolerance max(2.2*need, need+120) is generous; a column 120 px wider than needed passes.

## `test_counted_words`

**Result:** ✅ passed · [line 100](../../tests/test_ui_stability.py#L100) · **Area:** UI text: plural() · **Checked against:** hand-computed · **Strength:** 🟢 strong

Counted words read correctly (singular/plural, -y to -ies, thousands separator, invariant 'bp').

| # | Input | Expected output |
|---|---|---|
| 1 | plural(1,'primer'), (3,'primer'), (1234,'feature'), (2,'copy'), (0,'mismatch'), (2,'entry'), (5,'bp') | '1 primer','3 primers','1,234 features','2 copies','0 mismatches','2 entries','5 bp' |

## `test_the_new_tab_button_stays_beside_the_tabs_and_in_sight`

**Result:** ✅ passed in 0.43 s · [line 106](../../tests/test_ui_stability.py#L106) · **Area:** GUI: tab bar '+' button · **Checked against:** GUI state · **Strength:** 🟡 adequate

The '+' new-tab button sits right after the last tab and stays visible inside the bar when tabs overflow; quick filters take only their own width.

| # | Input | Expected output |
|---|---|---|
| 1 | Default window | 0 <= add.left - lastTab.right <= 8 |
| 2 | Add one tab per workspace.VIEWS entry | bar.rect() contains '+' geometry; quick_filters width <= sizeHint+4 |

## `test_the_toolbar_fits_without_overflow_down_to_1024_px`

**Result:** ✅ passed in 0.17 s · [line 123](../../tests/test_ui_stability.py#L123) · **Area:** GUI: toolbar width · **Checked against:** GUI state · **Strength:** 🟡 adequate

At a 1024x768 window every toolbar action stays visible (nothing pushed to the overflow menu).

| # | Input | Expected output |
|---|---|---|
| 1 | Window 1024x768 | no visible non-separator action whose widget is hidden |

**Review notes:**

- No docstring. L126 checks only the first QToolBar (currently the only one, main.py:453).

## `test_menus_in_the_usual_order_each_item_its_own_letter`

**Result:** ✅ passed in 0.21 s · [line 133](../../tests/test_ui_stability.py#L133) · **Area:** GUI: menu bar order and mnemonics · **Checked against:** property/invariant · **Strength:** 🟢 strong

Menus follow the usual order (File, Edit, View, Import, Export … Help) and every Alt+letter mnemonic is unique within its menu, recursively.

| # | Input | Expected output |
|---|---|---|
| 1 | Menu bar titles | first 3 'File','Edit','View'; then 'Import','Export'; last 'Help' |
| 2 | Every menu and submenu | no duplicate '&' letters among items |
| 3 | Top-level titles | mnemonic letters all distinct; first action has a QMenu |

## `test_the_view_picker_keeps_its_width_whatever_view_it_shows`

**Result:** ✅ passed in 0.33 s · [line 154](../../tests/test_ui_stability.py#L154) · **Area:** GUI: view picker width · **Checked against:** GUI state · **Strength:** 🟡 adequate

The pane's view picker button keeps one width whatever view name it shows, so the pane header does not jump.

| # | Input | Expected output |
|---|---|---|
| 1 | picker.set_current(name) for each workspace.VIEWS entry | exactly one distinct picker width |

**Review notes:**

- No docstring. L160 uses picker.set_current, not pane.set_view, so a width change caused by actually switching the pane's view would not be seen.

## `test_the_window_remembers_its_size_view_switches_and_column_widths`

**Result:** ✅ passed in 0.34 s · [line 166](../../tests/test_ui_stability.py#L166) · **Area:** GUI: persisted window settings · **Checked against:** round-trip/self-consistency · **Strength:** 🟢 strong

Window size, View switches (translation) and a hand-dragged column width come back in the next window; which side panel was open does not.

| # | Input | Expected output |
|---|---|---|
| 1 | Window 790x560, translation toggled, display dock shown, feature column 1 resized to 151, close; new MainWindow | size within 10 px; translation state as set; display dock not visible after load; column 1 width == 151 |

## `test_the_range_box_reads_the_usual_ways_of_writing_a_stretch`

**Result:** ✅ passed in 0.00 s · [line 190](../../tests/test_ui_stability.py#L190) · **Area:** Range box parsing (parse_range) · **Checked against:** hand-computed · **Strength:** 🟢 strong

The 'select from-to' box accepts the usual notations (1-based inclusive) and rejects bad input with ValueError.

| # | Input | Expected output |
|---|---|---|
| 1 | '101..200', 3780 bp, circular | (100, 200) |
| 2 | ' 101-200 ', '101 to 200' (circular), '200..101' (linear) | all equal to each other |
| 3 | '1,000..2,000' | (999, 2000) |
| 4 | '3700..50' circular / linear | (3699, 3830) through the origin / (49, 3700) |
| 5 | '42' | (42, 42) cursor after base 42 |
| 6 | '', 'abc', '1..2..3', '0..10', '100..4000' | raises ValueError |

**Review notes:**

- L193: the chained equality does not pin a value, so ' 101-200 ', '101 to 200' and linear '200..101' could all return the same wrong result.
- The single '0' special case accepted at app/main.py:123 is not tested.
