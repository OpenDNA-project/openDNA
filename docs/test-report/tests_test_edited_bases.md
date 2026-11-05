# `tests/test_edited_bases.py`

[← all test files](README.md) · source: [`tests/test_edited_bases.py`](../../tests/test_edited_bases.py)

## `test_ranges_when_opened_typed_replaced_pasted_deleted_and_shifted_by_later_edits`

**Result:** ✅ passed in 0.01 s · [line 31](../../tests/test_edited_bases.py#L31) · **Area:** Edited-base marks · **Checked against:** hand-computed · **Strength:** 🟢 strong

Checks which base ranges the document marks as edited after typing, replacing, pasting and deleting, and how later edits shift them.

| # | Input | Expected output |
|---|---|---|
| 1 | fresh 200 bp doc (seed 5) | edited_ranges() == [] |
| 2 | type A, C, G at 10; replace 50..54 by 'TT'; paste 'GGGGG' at 100 | [(10,13)]; [(10,13),(50,52)]; last range (100,105) |
| 3 | type 'AAAA' at 100; delete 0..10; type 'CC' at 0; delete 93..95 | [(100,104)]; [(90,94)]; [(0,2),(92,96)]; [(0,2),(92,94)] |

**Review notes:**

- Line 49 checks only the last range after the paste, not that the two earlier ranges are unchanged.

## `test_undo_redo_save_and_the_origin_keep_the_ranges_consistent`

**Result:** ✅ passed in 0.01 s · [line 67](../../tests/test_edited_bases.py#L67) · **Area:** Edited-base marks · **Checked against:** hand-computed · **Strength:** 🟡 adequate

Checks that edited ranges follow undo/redo, are cleared by a real save but not by autosave, and survive edits across the origin and rotation.

| # | Input | Expected output |
|---|---|---|
| 1 | type 'GG' at 20, delete 5..8, undo, undo, redo, redo | [(17,19)], [(20,22)], [], [(20,22)], [(17,19)] |
| 2 | type 'TTT' at 30; write_genbank (autosave path); save | [(30,33)] kept; [] after save |
| 3 | type 'A' at 0; undo x2; redo; modified=False | [(0,1)]; []; []; [] |
| 4 | circular: replace 190..205 (over origin) by 'GGGGGG'; undo; redo; type 'A' at 0 and rotate(50); undo x3 | ranges total 6 bases; []; same ranges; total 7 bases; [] |

**Review notes:**

- Lines 110 and 118: across the origin and after rotate(50) only the total marked length (6, 7) is checked, not the coordinates, so ranges shifted to the wrong place would pass.

## `test_the_view_paints_edited_bases_blue_in_both_themes_and_a_search_hit_stays_a_hit`

**Result:** ✅ passed in 4.92 s · [line 142](../../tests/test_edited_bases.py#L142) · **Area:** GUI: sequence view · **Checked against:** GUI state · **Strength:** 🟡 adequate

Checks that edited bases are drawn in the blue 'edited' colour in light and dark themes, without a yellow background, and that a search hit over them is still visible.

| # | Input | Expected output |
|---|---|---|
| 1 | SequenceView 900x400 per theme; no edit | 0 pixels of theme edited_fg |
| 2 | type 'ACGTACGT' at 10, cursor to 0 | > 0 edited_fg pixels; 0 yellow pixels |
| 3 | set_search_hits([(10,18)]) | > 100 yellow pixels; edited_fg pixels still > 0 |
| 4 | clear hits, undo | 0 edited_fg pixels |
| 5 | theme.ROLES['edited_fg'] light/dark | light lightness < 150, dark > 128; blue > red+60 and > green+40 |

**Review notes:**

- Line 158 comment says 'no wide block of blue either', but nothing checks it; a blue background would pass since only yellow is counted.
