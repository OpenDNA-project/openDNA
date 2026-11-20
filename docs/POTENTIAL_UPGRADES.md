# Potential upgrades

Improvements that are wanted or worth considering and are not built. Not bugs: what is wrong is in [KNOWN_ISSUES.md](KNOWN_ISSUES.md).
An item leaves this file when it is built (or when it is decided against: say so in one line in the change).

Last full triage of the docs folder: 2026-11-20. The old plan files (TODO, file-tabs plan, UI passes v1 and v2, publishing audit, test plan and
catalogue, recombineering spec, synthesis research) were read point by point and deleted; what was still open is here.

## 1. Command line and agents

- `ncbi get --region`: fetch only a stretch of a record.
- `repeats.md`: recommend `--min-len` for large inputs.
- `codon`: wording of the GC warning.
- A single Undo step for an overwrite across the origin in the Rust core (today two).
- Origin-wrap helpers in one place: `Document.sequence(start, end)` and `Document.overwrite(pos, bases)` accepting `end > length` on a circular
  document, replacing `codon_tools.circular_bases` and the copies in `editor.py`, `cas9.py`, `cli_tools.py`.
- One `mutation.span(feature_or_range, doc)` returning `{start, end, crosses_origin}` for every CLI row and window location string.
- `variants`: a separate `wraps` flag instead of the text "across the origin" stored in `alt`.
- `annotation.note_for`: keep the Features Database's own source on the hit instead of undoing a string concatenation.
- One `uikit.Job` (token, cancel, progress) for the background
  jobs of `synthesis_panel`, `ncbi_dialog`, `cas9`, `folder_search`, `structure_dialog`; one `DataTable.add_filter_actions` for the context menus.
- `annotate --no-plannotate`: express as an inverted flag in `tool_settings_block` instead of `use_plannotate = None`.
- `synthesis`: ViennaRNA install route for the command line like `wsltools install`; move the region-to-sequence orchestration into `synthesis_salis`.

## 2. Speed (measured by reading, not timed)

- `core/src/salis.rs` `maximal_repeats`: O(n * L^2); visit only the starts whose (k-1)-mer repeated, and test `covered` before allocating.
- `core/src/alignments.rs`: a full clip vector per clipped read (only the 24-base tail is used, and not at all on linear references); the circle-junction loop re-scans every sample per junction.
- `app/synthesis_salis.py` `_vienna_richest`: run the windows' folds in a process pool (about 10 s per kb today).
- `app/cli_tools.py` `records_of` and `cmd_split`: the whole multi-record file is scanned again for each record (N+1 scans); let `C.load` take the names already read.
- `app/search.py` `protein_hits`: search only the junction window instead of a second full rotated copy.
- `app/editor.py` `Editor.find` on a circle: search only the join instead of copying the sequence twice per find.
- `app/annotation.py` `merge`: sort by start and bisect instead of a scan per pLannotate hit.
- `app/value_filter.py` `cell_actions`: read one cell instead of building the text of every row.
- `app/plannotate_tool.py`: skip the cold WSL start when neither WSL nor bash exists.
- `idt.post` could call `structure.fetch` (and `structure.fetch` could return the HTTP error body) instead of a second urllib copy.

## 2b. Multithreading and background work (review of 2026-11-20; sizes estimated, not timed)

Goal: the window never freezes. Every analysis runs on at least one background thread (with the busy dialog or a progress bar), and the
heavy ones use the CPU cores allowed in Settings > General > Performance (`tool_settings.cores()`, as Local Align's `align` does). Rust work
must release the GIL (`py.allow_threads`) for a Python thread to help; pure-Python loops need a process pool or a port to Rust.

- `app/synthesis_salis.py` `_vienna_richest`: about 6 folds per base (18,000 for a 3 kb fragment); fold the windows in parallel (process pool
  unless ViennaRNA releases the GIL), then the order-dependent dedup in one pass.
- `core/src/annotate.rs` `find_parts_with` (Detect Features, Plasmid Checkup): library parts are independent; split them over `std::thread::scope`, hits kept in part order.
- `app/enzymes.py` `find_cuts`: about 250 enzymes x 2 orientations in Python; move to Rust (`search::find_motif` reads IUPAC) with the enzymes over threads;
  lifts the 250 kb limit (see KNOWN_ISSUES section 4).
- `app/folder_search.py` `search_folder`: files one after another; a thread pool over files (the parsing and `align` release the GIL), shared progress and cancel.
- `core/src/primers.rs` `binding_sites` via `lib.rs` `find_primers`: every database primer scans the whole template in turn; split the primers over threads.
  `editor.py` runs it on the GUI thread below 300 kb.
- `app/reference.py` `Reference.search` (genome checks of primers): records and 1 Mb windows in turn; a pool over them.
- `core/src/cas9.rs` `find_sites` (Spacer Finder specificity): strands x guides x PAMs over every chromosome; threads per strand or guide range. The
  automatic check calls it per 10 guides, rescanning the genome's PAMs each time: one call for all. `cas9.py` `candidates` is a Python regex over up to 2 Mb.
- `app/search_bar.py`: the search bar runs `search.search_sequence` on the GUI thread (about 1 s per chromosome for a motif or protein search); use its
  file-search thread pattern, and chunk `find_motif` / `find_protein_with` in Rust.
- `app/primer_design.py` `design`: 50 separate genome scans for the candidates; one `find_primers` call for all.
- `app/cli.py` `run_rows` (kld, pcr, recombineering, codon batches): rows in turn, mostly pure Python; a process pool keeping row order for large batches.
- `core/src/alignments.rs` `Bam::variants`: 1 Mb windows in turn; pile them in parallel, then merge junctions across window edges in one pass.
- `app/repeats.py` (pair repeats, dot plot) and `app/msa.py`: pure Python, threads give nothing; port to Rust first.
- Small: `salis::hairpin_rule_best_gc` (about 550 independent passes); annotation of the other records of a multi-record file (one at a time).
- Rust calls that keep the GIL, so `run_busy` cannot keep the window alive during them: `write_genbank`, `save_genbank`, `Track.search`, `Track.query`,
  `Fasta.fetch`; add `py.allow_threads`.

## 3. Features not built

- Enzymes: written permission from NEB or REBASE for their catalogues; CNNR-type enzymes (MspJI, FspEI, LpnPI, AbaSI); list every isoschizomer
  (drop the `by_kind` de-duplication in `scripts/build_enzymes.py`); owner decision that the Biopython-derived data (itself from REBASE) is acceptable (see THIRD_PARTY_NOTICES.md).
- Structure: a ColabFold provider and a native 3D viewer. Documents: a protein document type; ENA, UniProt and Addgene import; NCBI Datasets for genomes over 100 Mb.
- Parts library (`docs/PARTS_LIBRARY_METHOD.md` worklist): add Flp, Bxb1 integrase, Tn5 and Tn7 transposases, T3 and SP6 RNA polymerases, StCas9, DsRed, Renilla luciferase, TEV protease;
  a short-DNA collection with `strength` filled in (promoters, terminators, loxP/FRT/attB, origins); restore the T-DNA right border, the araBAD operators O1/O2/I1/I2 and the
  pEX18Gm/pSEVA terminators once their boundaries are cited; a precedence key in the loader instead of "last file wins".
- Synthesis: SCP4ssd as an optional external process (MIT, Python, `auto-sklearn==0.14.6`, 31 features; not verified on Windows); clients for Twist, GenScript GenSmart and Thermo GeneArt
  (their APIs were never verified: the research tool refused the pages); S-index (no code found).
- Verification PCR and FLP/Cre excision for recombineering were declined on purpose: not planned.
- Window state: restore dock state and table column widths, not only the window geometry.

## 3b. Planned in the old roadmap (the old PLAN.md roadmap, phases 5 to 7, deleted), not built

- Mutagenesis primers: the overlapping (QuikChange-style) layout; only the back-to-back KLD layout exists.
- Primer design: a rule that the pairs offered differ (3' ends at least 10 bases from pairs already kept), and a minimum and maximum product size.
- A "what is this?" search of a sequence against NCBI BLAST (no BLAST code exists).
- Read tracks: CRAM, downsampling, sorting and grouping reads, tracks that follow edits of the sequence, coverage bars coloured by allele counts (not re-checked in the code).
- A plugin hook: Python scripts that receive the Document API (no plugin module found).
- Codon tools: silent editing of a site while keeping the protein (no separate function found); saving GFF/BED tracks with the document (not checked).
- `core/src/document.rs`: an interval tree instead of a linear scan for feature lookup.

## 4. Interface (from the UI passes)

- Menu order: View after Tools is the open question (File, Import, Export, Edit, Feature, Primers, Databases, Tools, View, Help); the owner decides.
- Tab bar: a "+" inside the bar, quick filters as a tab-bar corner widget; one header pattern for the database tabs and the Primers view; a Reference combo with "Add reference genome…";
  the PCR pane without Cancel; database Status dialogs as a plain form.
- Instruction prose inside tool windows moved to tooltips or a "?" help; status tips for all menus.
- Check by hand, once, that value filters work on every table: popup values and counts, filter survives a refill, two columns together, light and dark themes, the 4.6 Mb genome filtering by Type in under a second.

## 5. Release, security and project

- Code signing (SignPath Foundation or Azure Trusted Signing) and SmartScreen reputation; a signature over `SHA256SUMS` (minisign or Sigstore) checked by the updater.
- CI: `permissions: contents: read` and actions pinned by hash.
- A licence header (SPDX) in the source files; decide "AGPL-3.0 only" or "or later".
- A signed contributor agreement instead of a checkbox; a registered trademark; building the PyInstaller bootloader from source if antivirus false positives continue.
- After a PyInstaller or PySide6 upgrade, repeat the by-hand Qt-swap test (`docs/PUBLISHING.md`).
- GitHub Discussions (the commented link in `config.yml`); one label vocabulary.
- A 3Dmol licence comment in the generated 3D page; a one-line note on the NCBI email and key in the NCBI dialog.
- Owner actions that cannot be checked from the code: "Enable private vulnerability reporting" and two-factor authentication on the GitHub account.

## 6. Tests

- Test the Sidewinder "too few barcodes" failure; compare hairpin and dimer free energies with a published ntthal value; check digest fragment sizes against a reference, not only their sum.
