# Changelog

What changed in each version, newest first. The start page shows the section of the running version; the same text goes into the release description.
Every feature or fix adds its line under "Unreleased" in the change that makes it (see CLAUDE.md).

## Unreleased

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
