# Alignment view tracks: plan and hand-over

The Alignment view (`app/genome_view.py`) must stay easy to use with 100+ tracks: manage them in bulk, read signal
(bigWig) tracks, highlight and recap regions, turn track items into genome features, and edit track files. The owner
approved this design on 2026-11-24; every decision below is settled. Work it as five stacked PRs, in order, each branched
from the previous one (`tracks-1-panel`, `tracks-2-signal`, ...). Tick a PR's box when it is opened, and keep the
"Where it stands" notes current so another agent can take over mid-PR.

## Status

- [x] PR 1: Tracks panel (branch `tracks-1-panel`)
- [x] PR 2: Signal tracks and track hub import
- [ ] PR 3: Regions, recap and print
- [ ] PR 4: Tracks into the genome
- [ ] PR 5: Editing tracks and keeping them in sync

Where it stands: PR 1 is open for review (OpenDNA-project/openDNA#34). PR 2 (branch `tracks-2-signal`, made from
`tracks-1-panel`) is open for review on top of PR 1. PR 3
(regions, recap, print) is next, on branch `tracks-3-regions` made from `tracks-2-signal`.

What PR 1 built, for the next PRs to build on:

- `app/track_layout.py`: the layout model (no Qt). Every change goes through a `Layout` method, which calls
  `on_change(edited)` (`Editor._layout_changed`: `tracks_changed`, and a hand-made change marks an `.odp` modified); kinds of
  rows that take a scale, a colour or nothing (`SCALED`, `COLORED`, `FIXED`) live here.
- `Editor.row_table`: `{row key: (kind, loaded track or None)}`, the one place that knows a file's rows. A loaded track has
  `kind` ("track", "signal" or "bam": `editor.ENTRY_ROWS`) and a name unique against every row (`Editor.add_track`);
  `Editor.loaded()` lists them all. Projects save each track's name so its rows keep their place and look.
- `app/track_panel.py`: the panel and `fill_menu`, the one menu of actions on tracks (also the view's right-click).
- `main.load_many` loads many files with the Alignment views' drawing paused (one redraw at the end), one message for all.

What PR 2 built:

- `core/src/signals.rs` (`opendna_core.Signal`): bigWig through the `bigtools` crate (read only; its zoom levels when a bin
  spans many bases), bedGraph and wig parsed here with zoom levels built the same way. `bins(seq, start, end, n)` gives
  (mean, highest) per bin. bigBed files are read by `Track` (by their contents) as the BED lines their entries stand for.
- Signal rows: `Style.render` (bars, line, heatmap), `summary` (mean, max), `range` (auto, group, fixed with `low`/`high`),
  `overlay` (the row a signal is drawn in). The canvas reads each signal once per view (`entry.cache`), in a worker thread (`_read_signals`).
- The canvas draws only the rows in sight, a screen above and below (`GenomeCanvas.band`, `scrolled`); dense rows have a 1 px gap.
- `app/track_hub.py` (no Qt): hub.txt, genomes file, trackDb with includes, groups from containers (views skipped), hidden at
  first (visibility hide, `parent X off`, a super track not shown), colours inherited; `copy_tracks`. `app/hub_dialog.py` and
  `main.import_hub` in the window; the `hub` command and recipe for agents; `demo/hub` is a small hub on the demo plasmid.

Measured on the Meso_v7 data (offscreen): a bigWig bins query 0.1-3 ms; 13 signals and a gene track redrawn in 15 ms.

## Test data

The owner's real use: the Meso_v7 UCSC track hub (Mesoplasma florum L1, 793,224 bp, circular):
`https://g-3e7620.6d81c.5898.data.globus.org/Meso_v7/Meso_v7.hub.txt` (genomes file `Meso_v7.genomes.txt`, track list
`L1/trackDb.txt`, sequence `L1/L1.2bit`). It holds 201 data files: 15 bigBed (gene sets, TSS, transcription units,
terminators, riboswitches, essentiality) and 186 bigWig (insertion density, HDTM per sugar and passage, RNA-seq), most
of them in composite groups and 23 hidden by default.

On the owner's machine a copy is in `C:\Users\AC\Downloads\Meso_v7_L1`: every bigBed as `.bed`, every bigWig as `.bw`
(`bigwig/`) and `.bedGraph` (`bedgraph/`), the genome as `L1.fasta` and `L1.2bit`. The reference genome to test with is
`L1_RefSeq2022.gb` (686 RefSeq 2022 genes). `L1_all_tracks_as_features.gb` holds every BED track as features and is a
mistake of the download, not a normal annotation: do not use it as a reference. Elsewhere, the `hub` command of PR 2
rebuilds the same copy.

## Rules for all five PRs

- Each PR follows `CLAUDE.md`: changelog lines, `tests.json` and the rebuilt report for every test added or changed,
  `/simplify` then `/code-review` before pushing, a recipe and a `test_<command>_...` test for every new command.
- Per-file state (track order, groups, styles, regions, sync shifts) is an attribute of the `Editor`, saved in the `.odp`
  by `exports.project_meta` (see "Where a setting goes" in `CLAUDE.md`). Panel width and similar go in `window_state`.
- No widget tests exist. Put the logic (layout model, parsers, recap, conversions, shifting) in modules without widgets
  and test those; check the panel and the drawing by running the app on the test data, and say in each PR which parts
  were checked that way only.
- The view must stay smooth with 200 tracks: draw from what is visible, keep work off the GUI thread
  (see the module docstring of `genome_view.py`).

## Facts that shaped the plan

- A loaded track (`core/src/tracks.rs`, `Item`) keeps start, end, strand, name, type, colour and note only: not the line
  it came from, not BED columns past the colour. PR 4 (extra columns as notes) and PR 5 (writing back) need both kept.
- Track styles (`GenomeCanvas.styles`: mode, height, scale) live in the view, keyed by track name, and are lost when the
  file closes. PR 1 moves them to the `Editor` and the project.
- Tracks do not follow sequence edits: inserting bases leaves track items after the edit off by that many bp. Features
  do follow them, by a rule in the core. PR 5 reuses that rule.
- Export > Image of the current tab already writes the Alignment view as SVG (`app/image_export.py`, vector while
  `image_export.ACTIVE`). PR 3 reuses it for `print` and "Export this region as image".

## PR 1: Tracks panel

- A Tracks button in the view's header opens a panel on the left of that Alignment view (one per view; its width is
  remembered in `window_state`).
- The panel is a tree of groups and tracks: per track a checkbox (shown / hidden), a colour swatch and its display mode.
  The built-in tracks (GC %, Bases, Features, Primers, Query hits) sit in a "Sequence" group and can be hidden and
  reordered like the others.
- A filter box ("A6S" keeps the 21 matching tracks) with "select all that match". Ctrl- and Shift-click select several.
- A bulk bar above the list (Show, Hide, Heatmap once PR 2 exists, Colour, More...) and the right-click menu apply one
  change to every selected track: display mode, colour (one colour, or a light-to-dark gradient over the selection),
  height, scale, range, move to group, rename (the label only, never the file), remove.
- Display modes: expanded, condensed, hidden, plus "dense" (one line).
- Drag and drop: reorder tracks, and move tracks between groups (an easy, obvious drop target: the owner asked for it).
- Groups: several files loaded at once from one folder form a group named after the folder; a single file goes at the
  top level; "Group by name" (split on `_`) is a manual action in the panel; a hub import (PR 2) keeps the hub's groups.
- In the view, each group is a header line; a click collapses all its tracks into that one line.
- A Fit button squeezes the shown tracks into the pane's height.
- Order, groups, colours, labels and display modes are saved in the `.odp`.

## PR 2: Signal tracks and track hub import

- Add the `bigtools` Rust crate to the core: read bigWig and bigBed; bedGraph and wig are read as text.
- Draw a signal as bars, a line or a heatmap. Heatmap colour: the track's own colour from light (low) to dark (high);
  tracks sharing a range share the scale. Range: automatic, shared by the group, or fixed min / max; mean or max when
  zoomed out; several signals can overlay in one row; the value under the mouse shows.
- File > Import > Track Hub...: a hub URL or a local folder. Reads `hub.txt`, the genomes file and `trackDb.txt`, lists
  every track in the hub's own groups, all ticked (hidden and grouped ones too), and downloads them into a folder the user
  picks (default: a folder named after the hub beside the open file). The hub's hidden tracks load hidden; each track
  shows the hub's `shortLabel`; `color` and `itemRgb` are kept. Never read tracks remotely: edits must land on a copy.
- Command `hub` (list a hub's tracks, download them); in tests, canned replies, never the network.

## PR 3: Regions, recap and print

- Highlights: right-click a selection > Highlight region... (name, colour), drawn as a band across every track.
  A Regions tab in the Tracks panel lists them: jump to, rename, recolour, delete, import and export as BED. Saved in the
  `.odp`.
- Recap: right-click a selection > Recap tracks here, and Recap highlighted regions. A table, one row per track (and per
  region when several), plus a summary row per group (the mean of its tracks); sortable, copy and export as CSV.
  Content per track: annotation items overlapping (name, start, end, strand, share of the region covered, extra
  columns); signals: mean, min, max, sum, and with per-base on (regions up to 1,000 bp) the value at every position;
  BAM: mean depth and the variants in the region; also the genome's own features in the region and the GC %.
- Command `recap START..END` (one JSON answer, `--tsv` one row per track, `--regions places.bed`, `--per-base`). The
  window and the command use the same function, so the numbers agree.
- Print: right-click a selection > Export this region as image..., and command `print <project.odp> START..END out.svg`
  (`.png` by the extension, default width 1600 px). It draws the view as the project shows it (visible tracks, order,
  styles, highlights). Without a project: `print genome.gb --track a.bed --track b.bw ...` with default styles. Reuse
  the existing SVG export.

## PR 4: Tracks into the genome

- Right-click an item > Add as feature.
- Right-click one or more tracks > Copy to features...: scope (whole sequence, view, selection, highlighted regions),
  feature type (a track named `*_CDS` defaults to CDS), name field. Each new feature gets a note "from track <name>";
  BED extra columns become notes, one per column; an item matching an existing feature (same type, span and strand) is
  skipped and counted in the summary. The whole copy is one undo step.
- On a signal: Regions above threshold... (value, minimum length, gap to merge) to highlights or features.
- Next / previous item of a track; export a track's items (in view, or all) as BED.
- Commands `track-features` (a track's items into a GenBank file as features) and `signal-regions` (stretches above a
  threshold, as BED).

## PR 5: Editing tracks and keeping them in sync

The riskiest part; its own PR for its own review.

- Editing, for BED, GFF, GTF, bigBed, bedGraph, wig and bigWig (VCF stays read-only: it is evidence from reads):
  - right-click an item > Edit...: a form with every field of its line (name, start, end, strand, score, colour, each
    extra column under its own name);
  - right-click a signal > Set value...: one value over the selection, or over the stretch under the mouse;
  - right-click an item > Delete; right-click a selection > Add as item to track... (pick the track, then the form).
- Track edits share the undo history with sequence edits (Ctrl+Z undoes the last change, a base or a value).
- Writing: edits stay in memory; the track shows a "•" in the panel; files are written on save (Ctrl+S, or Save track)
  after one dialog listing the files and the edits. The first overwrite of a file keeps its original as `<name>.orig`.
  A binary file (bigWig, bigBed) is rewritten whole. No command-line editing.
- Sync: the first edit that changes the sequence length while tracks are loaded asks "Keep the N tracks in sync with
  your edits?". Kept in sync: items after the edit move with it, by the same rule as features (an insertion inside an
  item lengthens it, a deletion shortens it, an item deleted entirely is removed and named in the status bar; a signal
  gets no value over inserted bases and loses those of deleted ones); undoing the edit undoes the shift. Left as they
  are: a warning says the tracks are off by X bp after position Y. The choice is remembered for the file and saved in the
  `.odp`; the panel shows each track's state ("in sync", or "off by +100 bp after 12,345").
- Synced tracks are never written on their own. On save, files shifted by sync need a double confirmation: first a
  dialog listing each file and what changes in it ("RefSeq_2022_CDS.bed: 412 items moved, 1 removed") with Write shifted
  files / Keep the files as they are; then "This overwrites N files. Originals are kept as .orig. Continue?" with Cancel
  as the default. A file with both hand edits and sync shifts gets the double confirmation.
- Kept as they are on disk: the `.odp` records the edits that apply to each synced track (position, bp added or removed)
  and replays them when the project opens, so the tracks are in sync again without touching the files.
- Record what stays unsolved in `docs/KNOWN_ISSUES.md`.
