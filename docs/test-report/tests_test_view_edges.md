# `tests/test_view_edges.py`

[← all test files](README.md) · source: [`tests/test_view_edges.py`](../../tests/test_view_edges.py)

## `test_primer_tails_hanging_off_the_edges_are_drawn`

**Result:** ✅ passed in 0.21 s · [line 47](../../tests/test_view_edges.py#L47) · **Area:** GUI: sequence view, primer tails at sequence ends · **Checked against:** GUI state · **Strength:** 🟡 adequate

A primer's 5' tail that hangs off the start or the end of the sequence (or across a plasmid origin) is still drawn, and row width does not change with zoom.

| # | Input | Expected output |
|---|---|---|
| 1 | 120 random bases (seed 3); primer F = 'GGTCTCAAATGC' + SEQ[:22], linear | site start 0, strand 1, tail 12; left == base_left; orange (ALMOST_COLOR) pixels in the last 12 columns of the rows |
| 2 | Same primer, circular | site.extent[-1] == (108, 120); orange pixels in the 12 columns before the end of the last row |
| 3 | No primers; set_size 6, 9, 11, 16, 24 | one distinct row width (left + bpl*cw) |
| 4 | Primer R = 'CCGGTTAACC' + revcomp(SEQ[-22:]) | strand -1, tail 10, end 120; orange pixels right of the last base (or at row start if no room) |

**Review notes:**

- L38-44, L58, L67, L87: pixel checks only require > 0 orange pixels anywhere in the column range over the whole viewport height, so a single stray pixel (or another orange element, ALMOST_COLOR also marks mismatched sites) passes; the tail's length/position is not checked.

## `test_cuts_on_both_strands_are_known_drawn_and_shown_by_the_digest`

**Result:** ✅ passed in 0.04 s · [line 99](../../tests/test_view_edges.py#L99) · **Area:** Restriction cuts (both strands) and their drawing · **Checked against:** hand-computed · **Strength:** 🟢 strong

Enzymes report the top- and bottom-strand cut positions (sticky, blunt, Type IIS on either strand, across the origin), the sequence view draws each cut on its strand, and the Digest tool draws only ticked enzymes.

| # | Input | Expected output |
|---|---|---|
| 1 | find_cuts('TTGAATTCAA', linear) EcoRI; find_sites | [(3, 7)]; find_sites [3] |
| 2 | find_cuts('AAGATATCAA') EcoRV | [(5, 5)] blunt |
| 3 | BsaI on 'ATGGTCTCAAAAAAAAAAAAAA' and 'TTTTTTTTTTTTGAGACCAT' | (9, 13) and (7, 11); 4-nt overhang both |
| 4 | find_cuts('ATTCAAAAAAGA', circular) EcoRI | (11, 3) wrapping the origin |
| 5 | set_cut_marks [(11,15,'EcoRI'), (40,40,'blunt')] on 120-bp view | red at top half x=11 not bottom; red at bottom half x=15 not top; red dotted join between; blunt red on both halves at 40; cleared marks leave no red |
| 6 | DigestDialog on 'ACGT'*10+'GAATTC'+'ACGT'*10+'GGATCC'+'TTAA'*5; tick/untick EcoRI, BamHI | no marks before ticking; EcoRI mark (41, 45); labels ['BamHI','EcoRI']; ['BamHI'] after untick; no marks after both off |

**Review notes:**

- L125-126: dead branch, view.layout_seq_top does not exist anywhere in app/, so the hasattr test is always False.

## `test_clicking_the_alignment_view_and_hovering_features_in_the_sequence_view`

**Result:** ✅ passed in 0.28 s · [line 157](../../tests/test_view_edges.py#L157) · **Area:** GUI: alignment view click and sequence view cursor shape · **Checked against:** GUI state · **Strength:** 🟡 adequate

Clicking a base in the alignment view moves the sequence view to it (centred, no jump for nearby clicks), and the mouse pointer is a hand over a clickable feature and an I-beam elsewhere.

| # | Input | Expected output |
|---|---|---|
| 1 | 30,000 random bases; click alignment canvas at base 20000 (synthetic press/release) | editor.start within 2 of 20000; scrolled and row visible; row at 20-80% of viewport height |
| 2 | Click 30 bases further | scroll value unchanged |
| 3 | Click base 150 | row visible and scroll value smaller |
| 4 | CDS '10..40' on 120-bp view; hover at base 25 over the lane, over the bases, and at base 80 on the lane | PointingHandCursor; IBeamCursor; IBeamCursor |

**Review notes:**

- L179-181, L215: events are delivered by calling mousePressEvent/mouseMoveEvent directly, bypassing Qt's event dispatch (e.g. hit testing via the viewport).
