# `tests/test_cut_ui.py`

[← all test files](README.md) · source: [`tests/test_cut_ui.py`](../../tests/test_cut_ui.py)

## `test_pill_layouts`

**Result:** ✅ passed in 0.00 s · [line 46](../../tests/test_cut_ui.py#L46) · **Area:** Restriction enzymes: cut pill layout · **Checked against:** hand-computed · **Strength:** 🟢 strong

Checks the four pieces of the 'pill' diagram that shows how an enzyme cuts (top-left, bottom-left, top-right, bottom-right) for each kind of end, for reverse-strand and origin-spanning sites, and the ellipsis for a cut far from the site.

| # | Input | Expected output |
|---|---|---|
| 1 | EcoRI on 'TTGAATTCAA' | kind '5′', pieces ('G','CTTAA','AATTC','G') |
| 2 | KpnI on 'TTGGTACCAA' | kind '3′', pieces ('GGTAC','C','C','CATGG') |
| 3 | EcoRV on 'TTGATATCAA' | kind 'blunt', pieces ('GAT','CTA','ATC','TAG') |
| 4 | BsaI on 'TTGGTCTCAAGGAATTTT' | kind '5′', overhang 'AGGA', site strand +1 |
| 5 | BsaI on 'TTTTAAAAGGAAGAGACCTT' (site on bottom strand) | site strand -1, overhang 'AGGA' |
| 6 | EcoRI on circular 'ATTCAAAAAAGA' (site across origin) | pieces ('G','CTTAA','AATTC','G') |
| 7 | made-up Enzyme('FarI','GGTCTC',((20,24),)) on 'AA'+'GGTCTC'+'ACGT'*8, cuts at 22/26 | exactly one inner non-real column with '…'; first and last columns '…'; kind '5′', overhang length 4 |

## `test_pill_tooltips_and_pixmaps`

**Result:** ✅ passed in 1.24 s · [line 75](../../tests/test_cut_ui.py#L75) · **Area:** Restriction enzymes: cut tooltip · **Checked against:** GUI state · **Strength:** 🟡 adequate

Checks that hovering a cut shows a tooltip with the drawn pill image and a description, that it is reachable from unique cutters and Digest cut marks, and that the pill images are drawn in both themes and differ by kind of end.

| # | Input | Expected output |
|---|---|---|
| 1 | enzyme_mark_tip for EcoRI on 'TTGAATTCAA' | an <img src=data:image/png;base64,...> of 200+ chars; text has 'G^AATTC', '5′ overhang 4 nt (AATT)', 'after base 3' |
| 2 | QLabel rendering that HTML | more than 500 pixels in the theme's pill_bg colour |
| 3 | hover.describe for the EcoRI unique-cutter site, and for Digest cut mark 0 with EcoRI ticked | both contain 'data:image/png' |
| 4 | pill_pixmap for EcoRI, KpnI, EcoRV, BsaI, light and dark theme | non-null, >100x50 px, contains pill_bg and cut colours; all four images pairwise different |

## `test_guide_look_and_interface_highlight`

**Result:** ✅ passed in 1.11 s · [line 124](../../tests/test_cut_ui.py#L124) · **Area:** GUI: Cas9 guide drawing and theme · **Checked against:** GUI state · **Strength:** 🟡 adequate

Checks how a saved Cas9 guide is drawn (red upward cut arrow, black guide text that fits a one-letter lane) in both themes, and that the app palette's highlight colours come from the theme.

| # | Input | Expected output |
|---|---|---|
| 1 | guide_doc(1) in an 1800x300 sequence view, light and dark | cas9_arrow colour red>190, green<60, blue<80; >10 red pixels near the cut; first red row within top 2 px and widest row in top 70% of the lane |
| 2 | same view; pixels inside the guide boxes | theme.ROLES['cas9_text'] == ('#000000','#000000'); colour is black; zero pure-white pixels in the boxes |
| 3 | SequenceView.set_size(6, 8, 11, 16, 24, 32) | lane_h == lh for each |
| 4 | theme.apply('dark') and ('light') | palette Highlight, HighlightedText, Accent equal theme colours highlight, highlight_text, accent |

**Review notes:**

- Line 145 'widths[first] <= widths[widest]' and line 147 'widths[first] < widths[widest] or widest == first' are always true (widest is the argmax and rows before `first` are 0), so 'points up' rests only on 'first <= 2 and widest <= 0.7*h'.
- Docstring says the arrow is 'hollow'; nothing checks that.
- Lines 159-164: 'text is black' is tested as 'no white pixels'; if no glyphs are rendered (no fonts offscreen, as the comment admits) the check passes vacuously.

## `test_spacer_finder_automatic_check`

**Result:** ✅ passed in 0.64 s · [line 194](../../tests/test_cut_ui.py#L194)

*No written description yet: the docstring and the checks below are taken from the code as they are.*

The first 100 guides are checked automatically against the reference (and sorting keeps working); without a
reference nothing is checked and the header says so; the automatic check can be cancelled.

Checks (input → expected output, as written in the test):

```python
assert total > 100
assert all((s is not None for s in tool.specs[:100])) and all((s is None for s in tool.specs[100:]))
assert f'100 of {total} checked' in tool.status.text()
assert tool.table.item(0, 6).text() != '–' and tool.table.item(0, 7).text() != '–', 'here and genome columns are filled'
assert tool.table.item(100, 6).text() == '–'
assert scores == sorted(scores)
assert tool.table.item(r, 5).text() == str(tool.guides[i].score)
assert 'no reference' in tool.table.horizontalHeaderItem(7).text().lower()
assert 'No background genome' in tool.table.horizontalHeaderItem(7).toolTip()
assert not tool.busy.isVisible() and all((s is None for s in tool.specs))
assert 'not checked automatically' in tool.status.text()
assert tool.busy.isVisibleTo(tool) and (not tool.check.isEnabled())
assert done < 100 and '(cancelled)' in tool.status.text()
```

## `test_right_click_adds_the_same_feature_as_the_button`

**Result:** ✅ passed in 0.31 s · [line 252](../../tests/test_cut_ui.py#L252) · **Area:** GUI: Spacer Finder table menu · **Checked against:** round-trip/self-consistency · **Strength:** 🟡 adequate

The Spacer Finder's right-click 'Add as feature' saves exactly the same guide feature as the Add button, and adds several selected rows at once.

| # | Input | Expected output |
|---|---|---|
| 1 | linear seq dna(80)+GUIDE+'TGG'+dna(60)+GUIDE+'AGG'+dna(60) (A/T filler); genome check off; row 0 added with the button in window 1 | a feature is added (reference) |
| 2 | window 2, same seq: right-click on row 0, column 2 | one action 'Add as feature'; triggered, its (kind, location, qualifiers) equal the button's; has the CAS9_KEY qualifier |
| 3 | select all rows, right-click on row 1 | one action whose text contains '(2)'; triggered, feature count +2 |
| 4 | right-click below the last row (viewport bottom-right) | build_menu returns None |

**Review notes:**

- L285: adding both rows re-adds row 0's guide (already added at L274); only the count is checked, not that the two new features are the two guides.
- L257/L266: the two Cas9Dialog windows are never closed; L287-288 close only the main windows, and not at all if an assert fails.

## `test_clicking_a_digest_fragment`

**Result:** ✅ passed in 0.05 s · [line 303](../../tests/test_cut_ui.py#L303) · **Area:** GUI: Digest fragment click · **Checked against:** hand-computed · **Strength:** 🟢 strong

Clicking a fragment in the Digest table selects it in the sequence and pins its two end cuts (drawn and hoverable); the pins are cleared by another selection or a new listing; a fragment across a circle's origin is one wrapping selection.

| # | Input | Expected output |
|---|---|---|
| 1 | linear 'T'*20+'GAATTC'+'T'*30+'GGATCC'+'T'*20, EcoRI+BamHI; click row '22'-'57' | selection (21, 57); pinned marks EcoRI and BamHI, all pinned; hover text has 'data:image/png' |
| 2 | sequence view 1200x300 of that editor | row 0 in enzyme_rows; cut_at on the first pinned cut returns 'digest_ends' or 'digest'; some 'cut'-colour pixels in a 6x8 px box at the cut |
| 3 | ed.select(0, 3) | 'digest_ends' gone and no table row selected |
| 4 | click row '1'-'21' (left end of the sequence) | selection (0, 21); only 'EcoRI' pinned |
| 5 | tab 0 then reload(); dialog destroyed | 'digest_ends' gone; then 'digest' gone |
| 6 | circular 72 bp 'T'*10+'GAATTC'+'T'*40+'GGATCC'+'T'*10, EcoRI+BamHI; click row '58'-'11' | ed.wraps; (start, end) = (57, 83); selected_spans [(57,72),(0,11)]; pinned BamHI and EcoRI |
| 7 | same circle, EcoRI only; click row 0 | end-start = 72; pinned ['EcoRI','EcoRI'] |

**Review notes:**

- L321: cut_at accepts either 'digest_ends' or 'digest', so it does not prove the pinned mark is the one hit.
- L318: 'any(r in view.enzyme_rows for r in (0,))' is just '0 in enzyme_rows'; it does not show the names are those of the pinned cuts.

## `test_hovered_cut_site_enzyme_names_and_map_click`

**Result:** ✅ passed in 0.11 s · [line 356](../../tests/test_cut_ui.py#L356) · **Area:** GUI: cut hover, map click, enzyme names · **Checked against:** hand-computed · **Strength:** 🟡 adequate

Hovering a cut draws its recognition site; clicking an enzyme on the map selects its site (circle and line); enzyme names stay on the line above the bases even with a primer there.

| # | Input | Expected output |
|---|---|---|
| 1 | 'TTTTTTTTTTGAATTCAAAAAAAAAA', EcoRI digest mark; hover_site of mark 0 and of missing mark 5 | (10, 16, 11, 15); None |
| 2 | view._hover_changed(mark 0), then None | hovered_site (10,16,11,15) and the viewport image changes; then hovered_site None |
| 3 | mark (1,1,'Cas9') | hover_site None (no recognition site) |
| 4 | random 300 bp (seed 7) with GAATTC at 120; MapView 500x500, circular and linear; click the EcoRI label | EcoRI in site_hits; _site_at(centre) is its key; selection (120, 126) and reveal emitted [120] |
| 5 | 'TTGAATTCAA'+30 random+'GGATCC'+60 random (seed 5), enzymes shown, no primer | _name_y(0,3,5) == seq_top - lh == 0 |
| 6 | add forward primer pF = seq[0:20] (temp primer DB) | seq_top == primer_h; _name_y at pos 3 and 60 == max(0, seq_top-lh); tip anchor x at base 3 and below.y > above.y |

**Review notes:**

- L372: 'after != before' only shows that something was repainted, not that the site GAATTC is the part highlighted.
- L418: the tip-anchor check only compares x and the order of y; it does not check the tip is above the name line.

## `test_the_enzyme_filter`

**Result:** ✅ passed in 0.04 s · [line 433](../../tests/test_cut_ui.py#L433) · **Area:** Restriction enzyme filter · **Checked against:** hand-computed · **Strength:** 🟢 strong

The enzyme display filter combines independent criteria (site counts, site size, palindromic, Type IIS, end type) that must all hold; the Digest's own filter bar lists enzymes without changing the views' filter.

| # | Input | Expected output |
|---|---|---|
| 1 | filter_doc: C/T random body (seed 11) with EcoRI x1, BamHI x2, AluI x3; default EnzymeFilter() | EcoRI shown; BamHI and AluI not |
| 2 | EnzymeFilter(min_sites=1, max_sites=None, min_size=6) | EcoRI and BamHI shown, AluI not; BamHI shown 2 times |
| 3 | EnzymeFilter(1, None, 1, None) | EcoRI, BamHI, AluI shown; shown_sites sorted by position |
| 4 | EnzymeFilter(2, 3, 1, 4) | names == {AluI, Sau3AI} (Sau3AI GATC inside the 2 BamHI sites) |
| 5 | wanted(): palindromic, type_iis, ends, min_size=8 on EcoRI, BsaI, KpnI, EcoRV, NotI | EcoRI palindromic, BsaI not; BsaI Type IIS, EcoRI not; EcoRI 5', KpnI 3', EcoRV blunt; NotI >= 8 bp, EcoRI not |
| 6 | combined palindromic+ends '5'+min_size 6 | EcoRI yes; KpnI and EcoRV no |
| 7 | DigestDialog on filter_doc; max_sites=1 | list ['AluI','BamHI','EcoRI','Sau3AI']; then ['EcoRI'], editor filter still default |
| 8 | min_sites 0, max_sites 0; then max_sites -1 | PstI listed with count '0', EcoRI not; then every COMMON enzyme listed |
| 9 | dialog shown | 3 preset buttons, same y, 3 different x |

## `test_display_panel_and_hiding_features`

**Result:** ✅ passed in 0.07 s · [line 494](../../tests/test_cut_ui.py#L494) · **Area:** GUI: display panel, feature hiding · **Checked against:** GUI state · **Strength:** 🟡 adequate

The display panel edits and mirrors the enzyme filter, its presets and the inventory criterion work, its check boxes line up, and a feature can be hidden from the views and shown again.

| # | Input | Expected output |
|---|---|---|
| 1 | panel: max_sites 3, palindromic True, ends '5' | filter == EnzymeFilter(1, 3, 6, None, True, None, '5') |
| 2 | set filter EnzymeFilter(2, None, 6, 8, False, True, '3') | panel shows min 2, max -1 ('any'), max_size 8, False/True/'3' |
| 3 | default filter; AluI and EcoRI with count 1 | == EnzymeFilter(1,1,6); AluI not wanted, EcoRI wanted |
| 4 | click each PRESETS button (Unique 6+, Type IIS, All sites) | editor filter equals that preset |
| 5 | inventory-only filter, empty inventory; tick EcoRI's Owned box in the Digest list, then untick | no sites; then EcoRI owned and the only shown site, digest tick unchecked, bar 'owned' lists ['EcoRI']; then empty again |
| 6 | panel 320x900 shown, CDS and promoter features | criteria folded; all squares at one x; expander shows criteria; all_primers off turns file and DB primers off, on turns them on |
| 7 | 'ACGT'*100 with misc_feature 10..40; menu 'Hide', then 'Show on Sequence'; set_features_hidden on/off | hidden: not visible, ed.shown == [], still in features, greyed in table; shown again: visible, shown == [f] |

**Review notes:**

- L523-526: the first preset (Unique 6+) equals the filter already set, so that click cannot fail; only the 2 others test anything.
