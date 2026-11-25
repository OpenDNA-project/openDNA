# Changelog

What changed in each version, newest first. The start page shows the section of the running version; the same text goes into the release description.
Every feature or fix adds its line under "Unreleased" in the change that makes it (see CLAUDE.md).

## Unreleased

## 0.9.0 (2026-11-25)

**Features**
- Added: the Alignment view's Tracks panel shows, hides, reorders, groups and restyles many tracks at once (filter, multi-select, drag and drop, colour ranges, Group by Name); a project keeps the tracks' order, groups and look, and is now saved as `.oDNA` (File > Save as Project; the old `.odp` files still open).
- Added: signal tracks (bigWig, bedGraph, wig) drawn as bars, a line or a heatmap, on their own range, their group's or a fixed one, several in one row, the value under the mouse; bigBed annotation tracks; Import > Track Hub reads a UCSC track hub and loads its tracks in the hub's groups, labels and colours.
- Added: highlighted regions across the Alignment view's tracks (Regions tab, BED import and export), a recap table of every track over a selection or the highlights (CSV), a region exported as an image, and track items copied into the sequence as features (Copy to Features, Regions Above Threshold of a signal, next / previous item); commands `recap`, `print`, `track-features` and `signal-regions`.

**Bugfixes**
- Two track files with the same name (from two folders) load as two tracks instead of sharing one display.
- Alignment view: a click in the Tracks panel outside a row no longer drags a selection over the tracks below it, and the names of pinned tracks no longer sit under the colour stripe.
- Alignment view: dragging the tracks sideways no longer drags the names column and the overview along with them.
- Loading many track files at once no longer redraws the view after each file: 200 files load in under a second instead of most of a minute, and their files are read side by side (120 BED files: 0.5 s instead of 2.6 s).

**UI/UX**
- Highlights: right-clicking one in the view offers the same Remove, Colour and Recap as the Regions tab, and several can be handled at once.
- The start page's opening message has GitHub Sponsors and Ko-fi links to support the program.
- Alignment view: hidden tracks take no room, a new one-line Dense display, collapsible group lines, and a Fit button that squeezes the tracks into the pane.
- Alignment view: only the tracks in sight are drawn, and dense rows sit closer together, so hundreds of tracks scroll smoothly.
- Loading many track files at once gives one message for all their warnings instead of one per file.
- Alignment view: the Tracks panel opens by default and the names of the selected tracks are highlighted; clicking a name no longer hides the track.
- Alignment view: the sequence and the features stay pinned in view when the tracks scroll, and any track can be pinned or unpinned (right-click); the bases come first and the GC % is hidden until shown.
- Windows file association also covers `.oDNA` projects; whoever accepted the other file types is asked once whether to add it.
- Alignment view: a narrower names column with wrapped names and the file type, a bar between the names and the tracks, and thin vertical guides (Guides button) at 100 bp to 100 kb or more.
- Alignment view: taller lanes for BED and other annotation tracks, rows as tall as the items in view need (Expanded shows every item at once), the bases drawn at every zoom (a line when far out), right-click a highlight to remove it, smaller panel buttons, a track menu regrouped (Display, Scale, Signal, Colour, Items, Group), both ends of a fixed range in one dialog, a click on a track in the panel scrolls to it, an icon toolbar and a status bar (range, pointer, selection), eye toggles and file-type tags in the panel, colour stripes on track names, a line under the pointer showing its position and signal value, a draggable overview, the wheel scrolls the tracks and Ctrl+wheel zooms.
- Alignment view: right-click the ruler or a selection for the region actions, an item of a track to add it as a feature or highlight it; the Tracks panel has a Regions tab.

## 0.8.8 (2026-11-24)

**Bugfixes**
- Renaming a primer of the file that is not in the Primer Database renames it instead of adding a copy beside it; an assembly's primers not saved yet keep their name until Save puts them in the database.
- Two primers never share a name: the database refuses a name a primer of the file has, and a file (or primer list) that brings a name twice opens with name_1, name_2.
- Features hidden by hand stay hidden after the file is saved and opened again (kept with the hidden primers).

**UI/UX**
- The features table has a Visible column: untick a feature to stop drawing it, tick it to bring it back.

## 0.8.7 (2026-11-24)

**Features**
- Added: the Primer Database's Match column says where each primer binds the file open in the window, those primers listed first; its In file box (or right-click > Import) adds a primer to that file or takes it out.
- Changed: a file opens with no background genome; pick one in the toolbar's genome box. Primer binding has a minimum Tm of a site (40 °C by default), on the sequence and in the background genome; a site's Tm now counts its mismatches.

**Bugfixes**
- Primers and binding sites hidden in a file stay hidden after it is saved and opened again (kept in the GenBank file as a comment other programs ignore).

**UI/UX**
- The background genome list opens at its full width and place instead of growing after its animation.
- Best site only is renamed *Show only the best binding site by default*: it sets which sites are ticked until you tick or untick one yourself.

## 0.8.6 (2026-11-24)

**Features**
- Changed: a FASTA file's topology is read from its header ([topology=circular] or [topology=linear]); with no tag, one sequence is circular and several are linear.
- Changed: importing a GFF3/GTF track asks which feature types to keep (with their number of lines) and offers each line's colour and note; projects keep the choice.
- Added: where a primer binds is set in Settings (3′ bases compared, how many must match exactly, mismatches); Best site only, and a single binding site can be hidden without hiding the primer. Hidden primers and sites now last until the file is closed (a project keeps them).

**Bugfixes**
- A GFF3 track made only of `region` lines is no longer read as empty: only a line covering the whole sequence is left unticked.
- A GFF3 feature across the origin (end before start) is one feature on a circular sequence instead of an unreadable line; on a linear sequence openDNA offers to make it circular.
- Restore Defaults on the Primers settings page now restores 2 mismatches, the number used before any change.
- The NCBI API key is no longer kept as plain text: it goes to the system's secure store, like the synthesis logins.
- The start page and the full changelog fit on the screen again: what is new and the sponsor column scroll instead of growing past it.

**UI/UX**
- Start page: larger, the logo top left with Full changelog beside the welcome lines, the sponsors over the whole height, the links in the bottom bar beside Continue; the Assembly button is gone.
- The start page opens at the first start of each day; *Show at every startup* (off by default) replaces *Don't show at startup*.
- The Display filter has a Best site only switch under the primer sources; right-click a primer site for Hide This Site, or untick a primer or one of its sites in the Primers panel's Visible column.
- Primers panel: a primer's sites are listed straight under it (a primer with one site is a single row), Match says whether each is on the sequence or the background genome; the panel is no wider than the Features view and opening a primer keeps the column widths.
- No background genome is set at first (it was the example E. coli genome): pick one in the toolbar's genome box or Settings → Primers. The box keeps the last 8 genomes used, however they were set, and its list opens wider to show their full names.
- Settings → Primers is built like the other pages: PCR conditions, Primer binding (with the background genome) and Primer checks.

## 0.8.5 (2026-11-20)

**Features**
- Added: Synthesis tool (synthesis success check with the problem stretches highlighted, IDT complexity screen) and pLannotate in Detect Common Features.
- Changed: Local Align aligns every candidate place on several CPU cores and finds hits it used to miss; an optional time limit warns when it cuts the search.
- Added: Settings > General > Performance (CPU cores openDNA may use) and Help > Run Setup Again; new `oDNA-cli` commands split, merge, synthesis and wsltools.

**Bugfixes**
- The window no longer freezes on analyses of long sequences: they run in the background and the status bar says which.
- Restriction sites are shown on whole genomes (the 250 kb limit is gone; a 4.6 Mb genome takes 1.2 s instead of 7.8 s).
- A digest with more than 1,000 cuts no longer freezes openDNA: it says so and does not draw them.
- Saving a large GenBank file, autosave and the search bar no longer block editing.
- A difference in the last bases of a Local Align query is listed as a mismatch instead of being trimmed off the hit.
- Variants, Golden Gate, recombineering, find and repeats no longer ignore what lies across the origin of a circular sequence.
- Every search limit that cuts an answer is now said in the status line instead of cutting silently.

**UI/UX**
- Detect Common Features shows the sequence beside the features it finds, with a "In the file" column; a feature annotated the other way round is shown and Add turns it round.
- Tables sort by any column, filter by value, and right-click a cell to show only or hide.
- Local Align has an "Aligned (bp)" column; Synthesis shows its nine measures in three columns with no scroll.
- The start page has up to 4 sponsor slots and a tip of the day.

## 0.8.4 (2026-11-19)

**Features**
- Changed: about 250 restriction enzymes instead of 47, with a filter for enzymes that cut only methylated DNA (DpnI).
- Changed: an Open DNA Project (.odp) also keeps each digest's ticked enzymes and filters, and the sequence view's letter size and spacing.

**Bugfixes**
- The circular map's title shrinks and wraps to fit inside the circle.

**UI/UX**
- File tabs keep at least 12 letters of a name and go on a second row when they do not fit.
- A slider under the zoom adds space between the letters of the sequence view; the search box sits against the range boxes.

## 0.8.3 (2026-11-18)

**Features**
- Added: opens Benchling GenBank files (its primers are read as primers), APE (.ape) and EMBL (.embl) files.
- Added: a command line for scripts and AI agents (oDNA-cli.exe): KLD, PCR primer and recombineering designs for a table of mutations or regions.

**Bugfixes**
- Changes through the origin of a plasmid are supported in the command line designs.

## 0.8.2 (2026-11-17)

**UI/UX**
- The README has a feature list and new pictures of the 0.8 interface.

## 0.8.1 (2026-11-17)

**Features**
- Changed: the PCR tab has three panes (template, primers, product); choose primers, click them on the template, or select a region and press Design. Designed pairs get a heterodimer check.
- Added: Edit > Set Origin to Selection Start.

**Bugfixes**
- Ctrl+E edits the clicked feature again.

**UI/UX**
- The menus are shorter and grouped: items that are on the toolbar, the right-click menus or the keys left them.
- Right-click on a feature: Delete feature, Hide feature, Edit feature, then its translation items.

## 0.8.0 (2026-11-17)

**Features**
- Added: several files in one window, each as a file tab with its own view tabs, side panels, search and undo history (View > Open Files in Tabs goes back to one window per file).
- Changed: double-clicking a file while openDNA is open sends it to that copy; autosave recovery copies are kept in the program's data folder.

**Bugfixes**
- The mouse wheel over a range box steps the right amount.

**UI/UX**
- The range box is two boxes, From and To, that follow the selection (wheel: Shift by 10, Ctrl by 100).
- The genetic code and the background genome are labelled in the second row.

## 0.7.5 (2026-11-17)

**Features**
- Changed: the multiple sequence alignment view is now the Multiple alignment tab of Local Align.
- Removed: the portable folder zip and tar.gz of releases (the installer and the AppImage remain).

**UI/UX**
- Hover Assembly in a pane's picker (or Tools) to choose the method and open the panel on it.

## 0.7.3 (2026-11-17)

**Features**
- Added: KLD method in the Assembly panel: edit a plasmid in place and get the primers for point mutations, insertions and deletions.

**UI/UX**
- Removed bases leave a small inverted ^ above the gap until the file is saved.
- Right-click on selected bases is grouped and has Translate; right-click on a CDS has Copy Translation; the search box is wider.

## 0.7.2 (2026-11-17)

**Features**
- Changed: releases carry one SHA256SUMS file instead of a .sha256 beside every download.

**Bugfixes**
- The selection, its brackets and the search marks are always blue instead of a system accent colour that can be hard to read.

## 0.7.1 (2026-11-17)

**Features**
- Added: Windows installs with a setup program (no administrator rights) and updates itself in place; macOS and Linux get a .dmg and an AppImage that update themselves.

## 0.7.0 (2026-11-14)

**Features**
- Added: Terminal tool: a bash session with pipelines to align reads, assemble with Unicycler and run MAFFT or Clustal Omega, with presets and a QC summary.

**UI/UX**
- The Alignment view labels every variant and marks where reads break and stretches with no reads.

## 0.6.4 (2026-11-13)

**Bugfixes**
- Circular sequences no longer list a second, shorter ORF inside one that crosses the origin.
- Error reports and the log hide file paths that contain spaces.
- Updates are only accepted with a checksum from openDNA's own releases.

## 0.6.3 (2026-11-13)

**Features**
- Added: forward / reverse primer from the selected bases with Ctrl+R / Ctrl+Shift+R.

**UI/UX**
- Local Align is laid out as Query, Search in and results; the buttons that run a tool share one blue.

## 0.6.2 (2026-11-12)

**Bugfixes**
- The updater finds releases whose files carry the version in their name.

## 0.6.1 (2026-11-12)

**Features**
- Changed: database primers are hidden by default when a sequence opens (Display > Primers > From the database).
- Changed: released files carry the version in their name.

**Bugfixes**
- Old openDNA file-association entries are cleaned up and follow the program when it moves.

## 0.6.0 (2026-11-11)

**Features**
- Added: restriction and ligation cloning, codon tools, degenerate primers and oligo order export.
- Added: multiple sequence alignment view, repeats dot plot, library and folder search.
- Added: a settings dialog; files are written more safely.

## 0.5.0 (2026-11-06)

**Features**
- Changed: Cut / Copy / Paste in Edit; Undo and Redo say what they undo; images export with a light background by default.

**UI/UX**
- Buttons and toolbars no longer move or change size; tables fit what they hold.
- Menus in the usual order (File, Edit, View...); one range box in the toolbar (101..200); View > Theme; Help > Troubleshooting.
- The window remembers its size, the View switches and the column widths you set; status colours are readable in light and dark.
