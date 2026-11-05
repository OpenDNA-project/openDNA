# oDNA UI/UX pass, round 2 (on main at 3278bc0)

Proposal only; nothing changed in the repository.

## Method

New harness (`scratchpad/harness2.py`). It renders, at 1400x850:
- the empty window
- each database tab, with and without a file open
- every one of the 13 pane views alone in its own tab
- Assembly, the display panel, both Database Status dialogs and Settings

Renders are in `scratchpad/shots2/`.

Automated checks, run on every render:
1. Every button: does hover, press or toggle change its size or move its neighbours?
2. Clipped labels and button text.
3. **Tables:**
   - columns more than 2.2x wider than their widest content (and at least 120 px extra)
   - columns narrower than their content
   - columns overflowing into a horizontal scroll bar
4. **The "+" button:** its distance from the last tab, and whether it is still inside the tab bar.
5. Toolbar overflow at 1920, 1600, 1400, 1280 and 1024 px.

New criterion added this round: **ease of use**. Controls sit next to what they act on; no space is wasted while something else is cut off; things that do the same job look and sit the same everywhere.

---

## 1. Placement: controls far from what they act on

| # | Finding | Evidence | Proposal |
|---|---------|----------|----------|
| P1 | **"+" new-tab button disappears once tabs overflow.** With 14 tabs it sits 151 px past the end of the bar, behind the scroll arrows. It is placed by hand after the last tab and the tab bar doesn't know it exists. | harness: `inside bar: False` | Make "+" part of the bar: reserve its width in the tab bar's size hint so it always follows the last visible tab. When tabs scroll, pin it just left of the scroll arrows. |
| P2 | Tab bar capped at **3/4 of the window** (`QuickFilters.SHARE`): the five display icons always take 350 px at 1400 px, so tabs scroll early. The icons only matter for document tabs but also show on the database and Assembly tabs. | workspace.py:250; renders db_*.png | Give the quick filters their natural width (about 5 x 30 px) as a tab-bar corner widget; hide them on tabs where they do nothing. |
| P3 | Digest: **"+ lane" is at the far right** of the window, while the lane tabs ("Digest 1") start at the far left. That is about 1300 px from the thing it adds to. | v_Digest.png | Put "+ lane" (or a "+" like the main tabs) right after the last lane tab. |
| P4 | Toolbar: **"Reference [BW25113_2026.gb] +"**. The bare "+" doesn't say what it adds, and the reference picker is a document setting sitting among view toggles. | 01_plasmid.png | Make it a combo whose last item is "Add reference genome…", and drop the "+". |
| P5 | Primer database vs Features database tab headers: the primer tab has New / Edit / Delete / Import / Export / Options / Help; the features tab has only Sources / Help. The Primers *view* has yet another layout (a second row of left-aligned New / Edit / Database / Options). | db_primers.png, db_features.png, v_Primers.png | One header pattern for both database tabs and the Primers view: title and count on the left; the actions that apply (New, Edit, Delete, Import, Export) as one group on the right; Help last. |
| P6 | PCR in a pane shows **"Cancel"**, which does nothing useful there. The primer combo boxes stretch across all 1400 px. | v_PCR.png | Hide Cancel when embedded; cap the combos at about 400 px; left-align the form. |

## 2. Tables too wide for nothing, and tables cut off

The common cause: the **last column stretches** to fill the width (Qt's `stretchLastSection`). In a narrow pane that's fine; in a tab of its own the last column takes most of the screen, while the column that needs room is cut off.

| # | Table | Measured (1400 px) |
|---|-------|--------------------|
| T1 | Features (view) | **Strand 828 px for 61 px of content**, while Name is cut off ("distal 34-nt of of natural FR…") |
| T2 | Primers (view) | Structure 658 px for 115 px; Sequence cut off at 230 px; Name 210 px for 62 px |
| T3 | Primer Database tab | Sequence cut off; columns overflow into a horizontal scroll ("Added" cut to "Add") **while the right half (630 px) is empty** |
| T4 | Features Database tab | Source cut off at 280 px (needs 502) **while the right half is empty** |
| T5 | Spacer Finder | Verdict 747 px for 155 px |
| T6 | Detect ORFs | Annotation 1068 px for 299 px |
| T7 | Detect Features | Source 692 px for 166 px |
| T8 | Plasmid Checkup | Details 997 px for 278 px; Notes 1090 px for 408 px; "Replicates in" and "Unit / CDS" cut off |
| T9 | Digest | Fragments table: 4 columns use 200 px of a 450 px box. Enzyme table: an empty band after "Owned". |
| T10 | Local Align | "Before this piece" stretched; # column next to Qt's row numbers (two number columns) |

**Proposal T-all:** one shared helper, `fit_columns(view, stretch=<column>)`, used by every table.
- Columns size to their content (header included), with a cap: long text gets elided, and the full text is in a tooltip.
- **Only the one descriptive text column stretches**: Name, Annotation, Verdict text, Source, Notes. Short columns (Strand, Length, GC, Ends, Tm) never do.
- If even the descriptive column has nothing to fill, the space stays empty at the right instead of pushing columns apart.
- Numbers are right-aligned everywhere. Today Length is left-aligned in Features and the Features database but right-aligned in Primers.
- One row height for every table, from the font. Today: 30 px in Spacer Finder, 22 px in Features, 15 px in the primer tree.
- No Qt row-number column where the table already has "#".
- Database tabs: the list gets about 65% of the width by default. The right side is empty until something is clicked, so it starts narrow.
- Primer tree rows that have no children don't need the 20 px expand-arrow indent.

## 3. Layout stability: still open from round 1, verified again

- B1. **The toolbar toggles still resize when pressed**: Lock 58x22 to 53x18, Filters and History 32x31 to 27x27, Circular 54x18 to 59x22. The cause is still theme.py:39 `PRESSED_QSS` (checked state only). The new Theme button and the database tabs don't add new resizing buttons.
- New: **Digest filter spin boxes** are fixed so narrow that their up/down arrows are a sliver, and "to" touches the boxes. "Palindr." is abbreviated to fit. The "Unique 6+ / Type IIS / All sites" buttons use a 9 px font, and the grey "29 of the 53…" hint is low-contrast. Fix: let the filter column size from its content (about 230 px), use normal-size spin boxes, and spell out "Palindromic".
- New: **menu-button arrow drawn over the text** in "Sources▾" (database tab, Detect Features) and "Options▾" (toolbar). Fix: a QToolButton with `MenuButtonPopup` / `InstantPopup` and proper padding, or the same "▾"-free style used for the pane header picker.
- B2-B13 from round 1 still apply (search counter, Cas9 status labels, Assembly method columns, updater buttons, and so on).

## 4. Toolbar overcrowding (worse since the Reference picker and Theme button)

| Window width | Hidden in the » overflow |
|--------------|--------------------------|
| 1920 | nothing |
| 1600 | nothing |
| **1400** | Circular |
| 1280 | Theme, Circular |
| 1024 | Reference, Filters, History, Theme, Circular |

Proposal, in order of impact:
1. **Theme button:** move it to View > Theme. It is set once, not used daily. That saves about 50 px.
2. **Reference:** move it into the Primers / Spacer Finder options it serves, or show it as a compact button "Ref: BW25113 ▾". That saves about 200 px.
3. **"Select from [ ] to [ ] Select":** the two spin boxes plus a button take about 350 px. Edit > Go to Position / Select Range (Ctrl+G) already exists. Proposal: one compact range box ("1..3780", Enter selects), as in the Alignment view's "Go to" box.
4. Use real toolbar separators (`addSeparator`) instead of the grey "  |  " labels.
5. If something still has to overflow, it should be the least-used item, not Circular. Order the toolbar so the rarely used items come last.

With 1 to 3, everything fits at 1024 px.

## 5. Still open from round 1 (unchanged by the new commits)

| Item | Status |
|------|--------|
| A1 Spacer Finder cannot be closed (cas9.py:276 `self.done = _Done(self)`) | **open** |
| A2 Enter fires the first button (ORF, Detect, Digest add-lane, primer dialogs) | open; the new database tabs are QWidgets, so they don't have the problem |
| A3 NCBI searches twice per Enter | open (ncbi_dialog.py:46, 234) |
| A4 to A10 (origin-bar click, OSError on save, Yes-default on destructive prompts, shortcut clashes, run_busy) | open |
| C1 window size and position not remembered; opens at 1400x850 | open |
| D2 Open Recent | **done** in ebfb9fe |
| D7 menu tooltips | one menu has them now; the others still don't |
| D8 Digest and Spacer Finder still open as modal dialogs from the Tools menu, while also being pane views | open. Proposal: the Tools menu opens them as panes or tabs (the pattern the database tabs now use), so no modal tool windows remain. |
| D11 "(s)" plurals | 62 places (was 53) |
| E1 to E11, F1 to F7, G1 to G4 | open |

## 6. Consistency between the new tabs and the rest

- **Empty states:** the database tabs show a centred grey hint ("Click a feature to see its sequence."), which is good. Do the same in Local Align ("Paste a sequence and press Align"), Digest ("Tick an enzyme to see fragments"), and the empty Assembly lower area.
- **Tool tabs vs dialogs:** Assembly and the databases are tabs; Digest, Cas9, ORF, Local Align and Detect are modal dialogs from the menu but panes from the view picker. One rule: everything opens in the workspace.
- **Instruction prose inside the UI** ("Tick enzymes for the selected digest; click a header to sort, the ladder to change it:", "Click a gray box to choose its source · drag a fragment to reorder · …"): move it to tooltips, the "?" help, or empty-state hints. Professional tools rarely carry permanent instruction sentences.
- **Database Status dialogs:** key/value pairs in an editable-looking grid table, with the first cell selected on open. A plain form layout reads better; or at least no selection, no grid, and the column fitted to its content.
- **Settings:** still 9 tabs scrolling in 697 px ("Plasmid Checki…" cut off). Use a sidebar list.

## Suggested order (unchanged in spirit)

1. Bugs: A1 to A10.
2. Stability: B1 plus the Digest spin boxes and menu-button arrows, with the permanent "nothing resizes or moves" test.
3. **Placement and tables: P1 to P6, plus the shared `fit_columns` helper (T-all), with tests that check the "+" is inside the bar and that no column is over 2.2x its content.**
4. Toolbar declutter (section 4).
5. Window state and sizing (C).
6. Conventions, look, responsiveness, wording (D, E, F, G).
