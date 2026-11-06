# `tests/test_map_view.py`

[← all test files](README.md) · source: [`tests/test_map_view.py`](../../tests/test_map_view.py)

## `test_map_layout_crowding_clicks_scale_and_geometry`

**Result:** ✅ passed in 0.28 s · [line 40](../../tests/test_map_view.py#L40) · **Area:** GUI: plasmid map layout · **Checked against:** property/invariant · **Strength:** 🟡 adequate

On a crowded map only the largest features that fit are drawn (a bigger map holds more), clicks hit only drawn features, the position scale picks readable 1/2/5 x 10^k ticks, and the circle/linear geometry is sized and centred.

| # | Input | Expected output |
|---|---|---|
| 1 | 6000 bp random seq (seed 1), circular and linear, CDS 'big' 101..1600 + 300 random misc_features (6-150 bp); MapView 500x500 | left_out > 0, drawn+left_out == all; 'big' drawn; smallest drawn >= largest left-out; band >= MIN_BAND; lanes keyed by drawn features |
| 2 | same views resized to 1100x1100 | fewer features left out than at 500x500 |
| 3 | circular, only 5 extra features, 500x500 | left_out == 0, drawn is ed.map_features, lanes is ed.map_lanes |
| 4 | crowded circular map: point at mid-angle of each of first 80 features, centre of lanes 0-3 | _hit returns None or a drawn feature index |
| 5 | ticks(length, px) for lengths 100, 999, 5936, 48502, 4641652, 10^7 at 300 and 1400 px | step in {1,2,5}x10^k; 2 <= length//step <= 40; spacing >= label width + 8 px; minor divides step and >= 4 px; 1400 px step <= 300 px step |
| 6 | tick_label (750,250),(2500,500),(15000,5000),(2500000,500000),(1,1000000) | ['750', '2,500', '15 kb', '2.5 Mb', '1'] |
| 7 | circular map resized 500x500, 900x600, 400x800; linear map 800x600 | 2*radius/min(w,h) == 0.8 +/- 0.01; linear backbone y in (150, 300) |
| 8 | both shapes, 10 features, 250 and 900 px, grab() | paints without error |

**Review notes:**

- Line 58: comment says it paints 'with the note' (how many are left out) but the note text is never checked.

**Fixed by the audit of November 2026:**

- The click check now also counts the clicks that hit a feature and asserts more than 40 of the 320 do, so a _hit() that finds nothing fails.

## `test_a_multi_fasta_lists_one_line_per_chromosome_and_a_click_opens_it`

**Result:** ✅ passed in 0.03 s · [line 134](../../tests/test_map_view.py#L134) · **Area:** GUI: map chromosome list (multi-FASTA) · **Checked against:** GUI state · **Strength:** 🟡 adequate

A FASTA with several sequences shows one clickable row per chromosome above the map; a click asks to open that record; a long list scrolls, opens on the current record and stops at the end.

| # | Input | Expected output |
|---|---|---|
| 1 | multi-FASTA chr1..chr4 of 300, 200, 100, 50 bp (seed 3), current = chr2, view 600x500 | _rows()==4, _top()>0, _room()==500-_top(); row i at y=2+17i+3 maps to chr(i+1); _first_row()==0 |
| 2 | left mouse press at (20, 2+2*17+3) | record_requested emits ['chr3'] |
| 3 | single-sequence FASTA (100 bp) | _top()==0 (no list) |
| 4 | 80 sequences (60..139 bp), current = index 70, view 600x400 | rows < 80, list height <= 0.45*400; first_row <= 70 < first_row+rows |
| 5 | scroll=0 then scroll=10000 | row at y=5 is chr1; _first_row() == 80 - rows (clamped) |

**Review notes:**

- Lines 140 and 144 hard-code a 17 px row height and 2 px offset, so a harmless row-height change breaks the test while the click logic itself is only checked against that constant.

## `test_human_genome_chromosomes_listed_drawn_and_zoomed`

**Result:** ⏭️ skipped in 0.00 s · [line 170](../../tests/test_map_view.py#L170) · **Area:** GUI: map on a real genome (hg38) · **Checked against:** regression value · **Strength:** 🔴 weak

> no human genome at V:/HM_Softwares/TestData/human_genome/hg38.fa

On the real human genome the chromosome list is read (chr1 length correct), chr21 opens as a linear sequence, its row is visible, and the wheel can zoom chr21 down to 50 bp.

| # | Input | Expected output |
|---|---|---|
| 1 | hg38.fa from OPENDNA_HUMAN_FASTA or V:/HM_Softwares/TestData/human_genome/hg38.fa | > 24 records, chr1 == 248,956,422 bp |
| 2 | open chr21, MapView 1000x700 | editor length == chr21 length, linear; chr21 row within visible rows; backbone y below the list |
| 3 | 80 wheel-up events at the centre of the chr21 backbone (1000x600) | window width rounds to 50 bp, 0 <= a < b <= length |

**Review notes:**

- Lines 170-171: runtime pytest.skip when the file is absent; the default path is a Windows drive (V:/...) outside git, so the test is skipped in CI and on Linux (it was skipped in the lead's run).

## `test_the_wheel_zooms_a_linear_map_at_the_pointer_and_shows_small_features`

**Result:** ✅ passed in 0.13 s · [line 209](../../tests/test_map_view.py#L209) · **Area:** GUI: linear map zoom and pan · **Checked against:** GUI state · **Strength:** 🟡 adequate

The mouse wheel zooms a linear map around the pointer down to 50 bp, showing all features in view; zooming out returns to the whole map; Shift+wheel and the slider pan/zoom; the map's zoom is independent of other views; circular maps do not zoom.

| # | Input | Expected output |
|---|---|---|
| 1 | linear crowded map (6000 bp, 300 features), 700x500, wheel up at 25% of backbone | view set; window < length; base under pointer moves < 1 bp |
| 2 | 30 more wheel-ups | window == 50 bp; no feature left out; every drawn feature overlaps [a,b) |
| 3 | 1 wheel-down, then 60 wheel-downs | _hit position == int(_pos_at(x)); finally view None and _win()==(0, length) |
| 4 | 10-feature map, 2 wheel-ups then Shift+wheel (mods 0x02000000) | window width unchanged (1e-6), start moved |
| 5 | zoom_slider 100 then 0 | window 50 bp; then view None |
| 6 | slider 80, ed.navigate.emit(5500 or 100 outside window) | window now contains that base |
| 7 | second MapView and a SequenceView on the same editor; slider -> 100 | other map view None; SequenceView.size unchanged |
| 8 | circular map 500x500, wheel at centre | view None and zoom bar hidden |

**Review notes:**

- Line 230: the hit position is compared with int(view._pos_at(x)), the same mapping _hit uses internally (map_view.py:808), so a wrong pixel-to-base mapping would go unnoticed here.
- Line 262: 'the sequence view is not moved' only checks SequenceView.size (letter size), not its scroll position or visible range.
