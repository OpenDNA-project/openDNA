# oDNA / openDNA: UI/UX pass, findings and proposals

Status: proposal only. Nothing in the repository has been changed.

## How this was tested

- Built the Rust core and ran the app headless (Qt offscreen, Fusion style) on `demo/pFG051_Chlor.gb` and `demo/BW25113_2026.gb`.
- Took screenshots of the main window at 1400x850, 1024x640 and 800x600, at 100% and 150% scale, in a light palette and a forced dark palette.
  Also took screenshots of 15 dialogs, the Assembly tab and the Plasmid Checkup tab (the Cas9 dialog failed to open in the harness; see A1).
  The screenshots are in `scratchpad/shots/`.
- Ran an automated stability check on every visible button: hover, press and toggle. For each, did the size hint change, and did siblings move?
  It also checked every button and label for clipping at the dialog's minimum size.
- Three independent read-only code reviews covered the shell and panels, the painted views and the tool dialogs.
  I re-checked the highest-severity claims by running the code.
- Baseline test suite: 214 passed, 1 failed, 3 skipped.
  The failure is `test_cut_sites.py::test_a_saved_guide_site_is_drawn`, a pixel count that depends on Linux font metrics. It is not related to the UI pass.

---

## A. Real bugs found during the pass (recommended regardless)

| # | What | Where | Verified |
|---|------|-------|----------|
| A1 | **Spacer Finder (Cas9) cannot be closed.** `self.done = _Done(self)` hides `QDialog.done()`, so Escape and the window X raise TypeError and the modal stays open. It has no Close button either. | cas9.py:274 | ran it |
| A2 | **Enter fires a button nobody chose.** With no button box, the first QPushButton becomes the default.<br>- Primer DB filter: Enter opens New…<br>- ORF min length: Enter adds CDS features<br>- Detect table: Enter adds every ticked hit<br>- Digest filter: Enter adds a lane | primer_panel.py:626, tools.py:173, annotation.py:186, tools.py:470 | ran it |
| A3 | NCBI search sends **two requests per Enter** (`returnPressed` plus the default button). | ncbi_dialog.py:43, 216 | reviewer ran it |
| A4 | Clicking just after the last base of a circular sequence (where you click to append) **toggles endless scrolling** instead of placing the cursor. | sequence_view.py:938, 1177 | read code |
| A5 | Unhandled OSError on save. A read-only or full disk gives a traceback.<br>The assembly save also overwrites existing .gb files without asking. | ncbi_dialog.py:198, assembly_panel.py:1062, structure_dialog.py:105 | read code |
| A6 | `Yes` is the default on destructive questions: "Delete primers", and "written into the file, no undo". | primer_panel.py:682, annotation.py:255 | read code |
| A7 | Genome-check thread swallows errors and stays on "searching…" forever. | primer_panel.py:492 | read code |
| A8 | Map Reads: Escape hides the dialog while bwa keeps running, and the QProcess is destroyed while it is still running. | runner.py:354-445 | read code |
| A9 | Shortcut clashes on Linux and macOS. FindNext is Ctrl/Cmd+G, which clashes with Go to Position. FindPrevious clashes with Alignment View. | main.py:271-273, 326 | Qt docs |
| A10 | `run_busy`: clicks pass to the main window during the first 400 ms, and Esc hides the busy dialog while the work continues. Together these allow Save or Close in the middle of an operation. | main.py:119-142 | read code |

## B. Layout stability ("buttons resizing when clicked")

| # | What | Where | Verified |
|---|------|-------|----------|
| B1 | **Toolbar toggles change size when pressed and push their neighbours.** Lock 58x22 to 53x18, Circular 54x18 to 59x22, Filters and History 32x31 to 27x27.<br>Cause: the app-wide `QToolButton:checked, QPushButton:checked {…}` stylesheet styles only the checked state, so the box model changes.<br>Fix: one rule that gives every state the same border and padding, as `QuickFilters` already does. | theme.py:39, 65 | harness |
| B2 | Search hit counter has only a minimum width; its text ranges from "" to "12 hits · 3 / 24…". In a narrow window everything to its right shifts. Fix: fixed width. | search_bar.py:64 | harness (only narrow) |
| B3 | Pane header button text follows the view name ("Map ▾" vs "Plasmid Checkup ▾"), so ⚙ and ? jump. | view_picker.py:33, 41 | read code |
| B4 | Cas9 status and genome labels have no word-wrap and texts of about 230 characters, so the dialog or pane grows wider. Local Align and primer DB status labels likewise. | cas9.py:232, 264; local_align.py:154; primer_panel.py:641 | read code |
| B5 | Cas9 progress row shows and hides between the table and the buttons, so the table jumps. | cas9.py:263 | read code |
| B6 | Assembly: changing method hides grid columns, so "Annealing Tm" and "Circular" slide sideways. The map's `setFixedHeight` changes as cards appear, so the sequence view below jumps. | assembly_panel.py:617, 199, 982 | read code |
| B7 | Labels that rewrap on every progress tick or row click resize the tables above them: primer check report, design details. | primer_panel.py:433, 507; primer_design.py:150 | read code |
| B8 | Updater dialog: button text "Later" becomes "Cancel" and "Install" becomes "Try again", and the progress bar shows and hides, so the buttons move. | updater.py:584-633 | read code |
| B9 | Display panel: the "hidden features" row shows and hides and pushes the type list. | display_panel.py:199 | read code |
| B10 | Alignment view header info label (stretch 3, text changes on every pan) can force the pane wider. | genome_view.py:1014 | read code |
| B11 | Protein view resets scroll and selection on every edit. The title wraps to 1 to 3 lines, so the text jumps. | protein.py:121, 88 | read code |
| B12 | Plasmid Checkup tables lose the current row and scroll position after every edit. | checkup.py:951 | read code |
| B13 | Image export warning label appears and grows the dialog. | image_export.py:123 | read code |

## C. Window, sizing, DPI

| # | What | Where |
|---|------|-------|
| C1 | **Nothing is remembered between sessions.** The window always opens at 1400x850, taller than a 1366x768 laptop screen (or 1080p at 150%). Fix: save and restore geometry, dock state and table column widths; on first run use about 80% of the screen, centred. | main.py:176 |
| C2 | Digest dialog minimum width is **1132 px**, because the long hint label doesn't wrap. That is too wide for a 1366x768 laptop at 125%. | tools.py (hint label) |
| C3 | At 1024 px the toolbar overflows (») and hides **Circular**, while a stretch gap is still left after the search arrows. | main.py toolbar |
| C4 | Fixed pixel metrics that clip with larger fonts or DPI:<br>- feature table rows 22 px<br>- map chromosome list 17 px<br>- Alignment-view lanes and ruler<br>- 74 px track-name column<br>- assembly card boxes and footer<br>- gel header<br>- display panel `FIELD = 64` and `font-size: 9px` preset buttons<br>- 250 px fixed assembly side columns | feature_table.py:169, map_view.py:36, genome_view.py:39/52/340, assembly_panel.py:67/239/528, tools.py:368, display_panel.py:158/276 |
| C5 | Settings dialog: 9 tabs in 697 px, so the tabs scroll and "Plasmid Checkup" is cut off. Fix: a sidebar list with stacked pages (the usual Preferences layout). | tool_settings.py |
| C6 | The Alignment view has no vertical scroll. Its minimum height is the sum of its tracks, so a few BAMs push the window past the screen. | genome_view.py:1062 |
| C7 | Tables stretch the wrong column: NCBI Title fixed at 380 px while "Updated" stretches; in the primer DB "Added" stretches instead of "Note". | ncbi_dialog.py:50, 225; primer_panel.py:633 |
| C8 | Search placeholder is cut off ("…or enzyme (Ent…"). | search_bar.py:58 |
| C9 | Cut pill uses a fixed `DPR = 2`, so it is blurry at 250-300%. | cut_pill.py:34 |

## D. Professional conventions (menus, dialogs, keyboard)

| # | What | Where |
|---|------|-------|
| D1 | Edit menu has no Cut, Copy, Paste or Select All (they exist only as keys in the view). | main.py:261 |
| D2 | No **Open Recent** submenu (the list exists), no Close Window (Ctrl+W), and Quit closes only one window. | main.py:198-206 |
| D3 | Save As has no shortcut on Windows (`QKeySequence.SaveAs` is unbound there). Add Ctrl+Shift+S. | main.py:202 |
| D4 | Menu order File, Import, Export, Edit, Feature, Primers, Tools, View, Help. The convention is File, Edit, View, …, Tools, Help. Proposal: File, Edit, View, Feature, Primers, Tools, Help, with Import and Export as File submenus, or straight after File. **Needs your call.** | main.py:198-361 |
| D5 | Undo and Redo are always enabled and unnamed. Proposal: grey them when empty, and label them "Undo Insert 5 bp" from the history text. | main.py:262 |
| D6 | Duplicate and missing mnemonics (Edit R twice; View C twice and T three times; Tools D, C and S twice; Help A twice; five items with none). | main.py |
| D7 | Menu tooltips never show (no `setToolTipsVisible`), and there are no status tips anywhere. | main.py:289, 332 |
| D8 | Tool windows are inconsistent.<br>- Digest, ORF, Local Align, Detect and Cas9 open **modal**, with **no Close button**, so you cannot look at the sequence while they are open.<br>- Design primers is non-modal; PCR uses OK/Cancel; Report and Junction dialogs have no Close.<br>Proposal: every tool window non-modal, with a standard Close button. | main.py:303-310, tools.py, assembly_panel.py:423 |
| D9 | Window titles mix Title Case ("Local Align") and sentence case ("Restriction digest", "New primer", "Primer options"). One tool has three names (Detect Common Features / Detect common features / Detect Features). Proposal: Title Case everywhere, one name per tool. | many |
| D10 | Ellipsis rules: "..." in about 8 places vs "…"; missing on items that open dialogs (Keyboard Shortcuts, Welcome, Check Mapping Tools, Add Primer); present on "Assembly…", which opens a tab; double ellipsis on "Check thresholds (…)…". | many |
| D11 | "N primer(s)", "(1 parts)", "1 mismatches": 50+ places. Proposal: a `plural()` helper. | many |
| D12 | No copy or export from result tables (Cas9 guides, designed primer pairs, assembly oligos, digest fragments, ORFs, detected features, NCBI). Proposal: one shared helper with Ctrl+C as TSV, a context menu with "Copy" and "Export CSV…", and sortable columns. | tools/* |
| D13 | Buttons enabled when they can do nothing:<br>- Edit/Delete with no selection<br>- "Save pair" with no pairs<br>- "Remove Tracks" with no tracks (it also has no confirmation)<br>- the pane ✕ on the last pane and the tab ✕ on the last tab | primer_panel.py:637/730, primer_design.py:233, main.py:204, workspace.py:121/420 |
| D14 | Context menus: right-click on empty space shows nothing; no Menu-key or Shift+F10 support. Ctrl+Home and Ctrl+End act like Home and End. PageUp and PageDown don't move the cursor. | sequence_view.py:1155, 1246; map_view.py:826 |
| D15 | Message boxes with Yes/No/Cancel whose meaning is written in the text: "Delete gene" ("Yes: gene and parts. No: only the gene") and "Recover changes?", where No silently deletes the autosave. Proposal: verb buttons ("Delete Gene and Parts" / "Delete Gene Only"; "Recover" / "Discard Recovered Changes"). | feature_table.py:385, main.py:555 |
| D16 | Selection readout in the status bar disappears on every `showMessage()`. Fix: `addPermanentWidget`. | main.py:158 |
| D17 | No drag-and-drop to open files. | main.py |
| D18 | Raw error text ("Unexpected error: {e}", unescaped HTML) in the NCBI and Cas9 dialogs. | ncbi_dialog.py:102/118/268/288, cas9.py:434/497 |
| D19 | Welcome page opens window-modal every day, even when the app was started by double-clicking a file. Recent files open on a single click; Enter does nothing. | sponsors.py:136, 349 |

## E. Visual identity and theming

| # | What | Where |
|---|------|-------|
| E1 | Toolbar separators are `QLabel("  |  ", color: gray)`. Use `addSeparator()`. | main.py:372, 383 |
| E2 | Lock button is an emoji ("🔓 DNA"). It shows as empty boxes on systems without an emoji font, and the project's own CONTRIBUTING says no emojis. Pane headers use font glyphs (⚙ ◫ ⬓ ✕) whose size and baseline depend on font fallback. Fix: draw them in icons.py like the funnel and hourglass. | main.py:1141, workspace.py:111-121 |
| E3 | Tool status colours (ok/warn/bad) are fixed hex: `#f9a825` yellow text on white is about 1.9:1 contrast, and green and red on dark are about 2.5:1. The two places also disagree on the warning colour. Fix: theme roles with light and dark pairs. | primer_panel.py:189, assembly_panel.py:37, checkup.py:78 |
| E4 | Assembly fragment boxes are fixed light greys, which glare in dark mode. | assembly_panel.py:250-279 |
| E5 | Other fixed colours: search hit yellow, primer purple (about 2.3:1 on dark), `gray` and `#c62828` inline in about 10 places, the tooltip tail `#bf360c`, protein highlight text. | sequence_view.py:42-53, hover.py:39, protein.py:231, … |
| E6 | **Image export in dark mode gives dark figures** (wrong for papers and slides). The export also includes the scrollbar, zoom bar and Alignment header. Proposal: a "Light background" option, on by default, and hide the UI controls while exporting. | image_export.py:54-57 |
| E7 | Tab close buttons are bright red-orange squares on every tab. They draw attention and differ from the platform's look. Proposal: a plain grey ✕ drawn in icons.py, red only on hover. | workspace.py |
| E8 | Empty state is one centred line of text. Proposal: a start page in the same window with Open / New / Import from NCBI / Demo and the recent files, plus "drop a file here". | main.py:154 |
| E9 | Main sequence painter has no antialiasing, so primer arrows, cut crosses and arrowheads look jagged. Row numbers can run into the first base on very long sequences. | sequence_view.py:541, 228 |
| E10 | Sequence view gives no indication of which pane has keyboard focus: the selection and cursor look the same with or without focus. | sequence_view.py:546, 655 |
| E11 | Crowded labels overlap: the "N small features not shown" note sits over a map label; the primer name is drawn over the primer's own bases in the sequence view (oAC734, oAC640 in the demo). | map_view.py, sequence_view.py |

## F. Feedback and responsiveness

| # | What | Where |
|---|------|-------|
| F1 | Heavy work on the GUI thread with no busy cursor:<br>- primer design (window blank until done)<br>- Local Align<br>- Cas9 candidates (up to 2 Mb)<br>- Detect features<br>- ORFs (re-run on every spin-box step)<br>- self-test<br>Fix: route through the existing `run_busy`, and debounce the ORF spin box. | primer_design.py:182/267, local_align.py:210, cas9.py:326, annotation.py:212, tools.py:216, selftest.py:47 |
| F2 | Autosave writes GenBank on the GUI thread every 60 s, which freezes on chromosomes. | main.py:1309 |
| F3 | Genome download and structure lookup have no Cancel and no progress bar. | ncbi_dialog.py:325, structure_dialog.py:70 |
| F4 | Drag-select in the sequence view does not autoscroll at the edge. | sequence_view.py:1204 |
| F5 | Wheel: every event is one zoom step, so touchpads zoom far too fast; horizontal Ctrl+scroll zooms out. Map drag selects after 1 px movement (no drag threshold). | sequence_view.py:201, map_view.py:193/870 |
| F6 | Hover tips are recomputed and moved on every pause over the same item; the origin tooltip lingers for about 10 s. | hover.py:105-143, sequence_view.py:1212 |
| F7 | Clicking the Alignment view steals keyboard focus without using any keys, so typed bases go nowhere. | genome_view.py:108 |

## G. Wording and tone (need your call)

- G1. The Welcome/About page says, in bold, "If you paid for it, you were scammed." and "Early test version: keep copies of your files." For a 1.0 feel I propose "openDNA is free software. The official download is github.com/OpenDNA-project/openDNA." The test-version warning goes when you declare 1.0.
- G2. Developer tools shown to users: View > Performance Timings, the Settings "Self-test" tab, the primer DB window title showing the file path, citations mentioning `scripts/build_parts.py`, and live counts in the ORF window title. Proposal: move them to Help > Troubleshooting and clean up the titles.
- G3. "Cursor after base 0" in the status bar; proposal: "Cursor at start".
- G4. "Unique Cutters" (View menu) toggles every drawn enzyme site, while the panel calls it "Restriction enzymes". Proposal: one name, "Restriction Sites".

## What I recommend NOT changing (already good)

- QuickFilters' stylesheet (stable in every state). It is the model for B1.
- Theme roles in theme.py, icons drawn in palette colours, and redraw on palette change.
- Threaded load/save/export with a delayed busy dialog; stale-result tokens in search, Cas9 and NCBI.
- Autosave and recovery that never touch the user's file; the Save/Discard/Cancel prompt that lists other edited records.
- The F1 shortcut list built from the menus.
- View caches that allow for the device pixel ratio; the Alignment view's sliding pan cache.
- Validation that keeps dialogs open and focuses the bad field.
- Spin boxes mostly with units and ranges; thousands separators; en-dash ranges.

## Proposed order of work (each a separate commit, tests added and the suite run)

1. Bugs: A1-A10.
2. Layout stability: B1-B13, plus a **permanent regression test**. It is the harness check: every checkable or pressed button keeps its size hint, and no sibling moves.
3. Window state and sizing: C1-C9.
4. Conventions: D1-D19, after your answers on D4, D8 and D9.
5. Visual identity and theming: E1-E11.
6. Feedback and responsiveness: F1-F7.
7. Wording: G1-G4, as you decide.
