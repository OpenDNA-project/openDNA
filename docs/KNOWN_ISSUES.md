# Known issues

The central list of problems that are known and not fixed: things that are wrong, missing a warning, or behave against what the text says.
Improvements nobody has reported as wrong (ideas, features not built) are in [POTENTIAL_UPGRADES.md](POTENTIAL_UPGRADES.md).

Agents and people add to it when they find one (one line: file, what is wrong, who sees it, how it was found). A fixed item is removed in
the change that fixes it, with its test; the ones fixed for the origin and "never ignore a location" rule are in `audit-ignored-locations.md`.

Rule behind most of the first two sections: a tool never silently ignores a location (near an end, across the origin of a circular plasmid,
past a cap, in a repeated stretch). It handles it or says what was left out.

Last full triage of the docs folder: 2026-11-20 (each older plan and audit was read point by point against the code; what was still true is here).

## 1. Silent drops (the tool leaves something out and does not say)

Core (Rust):
- `core/src/cas9.rs`: a guide over 32 bases or a sequence shorter than a site returns an empty list ("no sites"); a PAM base that is N never fits.
- `core/src/search.rs` `find_motif`: sequence letters other than A, C, G, T, N never fit (U, R, Y...). `fasta.rs` `search_motif` is linear only.
- `core/src/annotate.rs`: the protein pass runs only when no DNA hit was found; overlapping copies of a part are merged (the test is linear);
  parts cut by a linear end are never reported.
- `core/src/bio.rs`: an ORF that runs off a linear end is dropped. `document.rs`: `features_in` and `gc_bins` do not wrap a window.
  `write_genbank` drops a stored primer_bind not matched by the supplied sites.
- `core/src/dna_file.rs`, `embl.rs`, `tracks.rs`: features, binding sites or a final record without `//` dropped silently; a BED/GFF line
  with end < start fails the whole track.

App (Python):
- `app/editor.py`: `MAX_SITES` (50) and the track-search cap (500, also `app/search.py`) cut silently.
- FASTA records are searched as linear (`app/search.py`, `app/reference.py`, `app/folder_search.py`), not said.
- `app/reference.py` `split_genbank`: a last record without `//` is never searched.
- `app/msa.py`: no warning when the alignment touches the band edge, none for circular inputs with different origins.
- `app/primer_design.py`: candidates are cut to 25 per side before the uniqueness check.
- `app/assembly.py`: Type IIS sites near a linear end are dropped as "no site".
- Locations like `(10.20)..30` and a lone remote location (`J00194.1:100..202`) are kept in the file but have no extent, so they are not drawn
  (the test only checks the text round trip, `tests/test_core_formats.py`).

## 2. Circular origin not handled (a feature through the origin is treated as linear)

- `app/annotation.py`, `app/genes.py`: gene groups, `already_there`, `span_of` and `import_files` treat a part across the origin as linear
  (a part is shown as not annotated and added twice; a gene master covers the whole plasmid).
- `app/tools.py`: an ORF through the origin is translated only up to the origin; `annotated()` never matches it.
- `app/protein.py`: extent, pin, highlight and residue selection of a CDS through the origin are linear (window only; found by reading).
- `app/checkup.py` `_show`: a unit through the origin selects only the part before it.
- `app/mutation.py`: `c.-N` does not wrap on a circle; only the first CDS a change falls in is evaluated; an insertion exactly at the origin
  of a CDS through it is "non-coding".
- Inserting bases exactly at the origin, strictly inside a feature that runs through it, does not grow the feature (not re-run since the audit).
- `app/local_align.py` / the paired Compare view: the topology is ignored (linear); the junction gap of a strand change is not wrapped.
- Recode window: the Undo of a CDS through the origin takes two steps (two overwrites), not one.

## 3. Wrong or inconsistent results

- The `source` feature does not grow when bases are inserted at either end (`core/src/document.rs`; `demo/pBAD30_eilR.gb`: `source` 1..3559 stays
  after an append and shifts to 4..3562 after an insert at base 1).
- GFF3 export writes "." as the CDS phase (`app/exports.py:68`); the GFF3 specification needs 0, 1 or 2.
- Sidewinder assembly: the strand layout is marked "Simplified. Check the strand layout against the paper" (`app/assembly.py`, docstring near line 600);
  the barcode arm and the toehold single-strandedness are not settled. A "too few barcodes" failure after the 5,000-draw search has no test.
- The high-fidelity overhang set cannot fail the `pair_risk` orthogonality check (`app/overhangs.py` exempts the published pairs), so its test
  proves less than its name says.
- Assembly "Save primers to database" is not an undo step.
- Primer dimer and hairpin energies score only perfectly paired stretches of 2 bp or more (`core/src/primers.rs`, a `shortcut:` comment): no internal loops, mismatches or
  dangling ends, so they are not a full nearest-neighbour (thal) model.
- `core/src/document.rs`: feature lookup is a linear scan (a `shortcut:` comment says to switch to an interval tree for very large files).
- `app/structure.py` `fetch` does not refuse a non-https URL in general (every caller passes a fixed https URL; AlphaFold's is checked separately).
- Right-click Paste in a sequence field only goes through the validator (the old input audit, not re-checked after `clean_pasted`).
- The NCBI API key is stored as plain text in `ncbi.ini` (`app/ncbi_dialog.py`, `app/ncbi.py`); `app/secret_store.py` (system keyring) is not used for it.
- `app/idt.py` is not verified against the live IDT service (its header says so): endpoints, scope, error bodies and terms of use are as published by third parties.

## 4. Interface (found by the 2026 UI passes, still present)

- A click within 4 px of the red origin bar after the last base toggles endless scrolling instead of placing the cursor (deliberate in the code comment, but it was the original complaint).
- Tools opened from the Tools menu are modal `.exec()` dialogs while the same tools also exist as pane views.
- Tool names differ between the setting and the window (`Detect Features` vs `Detect Common Features`); some titles are in sentence case.
- About 56 "(s)" plurals remain where `uikit.plural()` exists (`app/assembly.py`, `cas9.py`, `enzymes.py`, `pcr.py`...).
- Result tables (guides, primer pairs, oligos, digest fragments, ORFs) have no shared Copy or Export (`DataTable` has none; only `kld.py` and `folder_search.py` copy).
- Hard-coded colours bypass the theme (`kld.py`, `assembly.py`, `assembly_panel.py`, `hover.py`, `genome_view.py` BASE_COLORS).
- Menu tooltips are on only two menus; the other menus have none.

## 5. Project and repository

- `CONTACT_EMAIL` is empty (`app/project.py`), and `SECURITY.md` and `PRIVACY.md` still say `<project email>`.
- The GitHub labels do not match `.github/labels.yml`: the `type:` and `status:` labels do not exist on the repository, so the
  issue forms apply labels that GitHub drops; there are two vocabularies (`status: ...` in `.github/labels.yml`, `needs-triage`... in `docs/agents/triage-labels.md` and CLAUDE.md).

## 6. Test gaps

- Hairpin and dimer free energies are not compared with a published ntthal value; digest fragment sizes are checked only by their sum.

## Found by agents (add new ones here, newest first)

- 2026-11-20, `app/cli_tools.py` `cmd_wsltools`: `P.status()` and `runner.installed_tools()` each start WSL (1 to 2 s each), a third time after an install.
