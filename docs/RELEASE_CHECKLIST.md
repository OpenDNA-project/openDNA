# Release checklist: changes since 0.8.4 and what to test by hand

Automated state when this was written: 159 Python + 68 Rust tests pass, pyflakes clean, test report rebuilt (227 described).
Nothing is committed. No window was opened for any change below: **every GUI line is checked by reading only** (no widget tests exist).
Open problems are in `docs/KNOWN_ISSUES.md`. Rebuild the core first (`maturin develop`), then restart the app.

## What changed

### 1. Command line for agents (`oDNA-cli`)
- New commands: `split` and `merge` (multi-FASTA / multi-record GenBank; `FILE#NAME` picks one record; `--names`, a clash uses the file name),
  `synthesis`, `wsltools` (pLannotate and Terminal tools in WSL).
- `msa` takes a multi-record file. `pcr` stays at window 0 unless asked; bad primers give a warning that proposes a wider window and never applies it.
- Every limit reached is a warning in every command (`answer(caps=...)`): guides `--top`, repeats, align pieces, search hits, NCBI, PCR pairs.
- `checkup` origins are structured fields (no leak of internal text). `codon` works for a CDS across the origin. `annotate` has `--no-plannotate`.
- Agent docs: `app/agent/index.md` (new rules), recipes updated, new `split.md`, `synthesis.md`, `wsltools.md`; user guide `docs/AGENT_GUIDE.md`.

### 2. Circular origin and "never ignore a location"
- `variants` (Rust core + Python): nothing near an end is dropped; clips continuing the circle are reads over the origin; a deletion through the
  origin is one deletion; an uncovered start and end are one gap; hard-clipped reads near an end are reported as "could not be checked". `--apply` writes it.
- Golden Gate closing junction on a circle; recombineering left homology through the origin and exact amplicon site; PCR narrowed window said;
  DNA find and "selection" search scope through the origin; protein motif through the origin; repeats (2.5 copies, palindromes, cut repeats counted);
  guides over N counted; a residue range split by the origin is an error; protein motif list cap said.
- Per-tool detail: `docs/audit-ignored-locations.md`.

### 3. New features (from other work in the tree; read-only reviewed)
- **Synthesis tab** (Tools menu, right-click "Test Synthesis…"): Synthesis Success Calculator (Rust `salis.rs`, forest `app/salis_forest.json.gz`,
  ViennaRNA install helper) and an IDT complexity check with credentials in the system keyring (`secret_store.py`; IDT code is unverified against the live service).
- **pLannotate** in Detect Common Features: "Found by" column, install dialog, credit line, `Also search with pLannotate` setting.
- **Tables**: shared `DataTable` (sort, filter by value, right-click a cell: show only / hide) in Detect, Checkup, Local Align, Folder Search and more.
- Welcome page: up to 4 slots, slots that expire (`until`), a heading, tips in the status bar every 5 minutes.
- Other: NCBI "more results" note, structure dialog and citations edits, CI and packaging files.

### 4. Cleanup from this review
- One origin-wrap helper reused, shared complement functions, `secret_store.save` keeps a copy in memory only without a secure store, the Synthesis
  tab imports lazily, `salis_repeats` and `salis_hairpin_gc` release the GIL.
- Code review found 10 items, **not fixed**: checkup `cds` start/end across the origin (`cli_tools.py:562`), `synthesis` region end across the origin,
  `sponsors.load` with a non-string `until`, plus efficiency items (`maximal_repeats`, per-read clip vector, sequential ViennaRNA folds, rescans per record).

## Test by hand before release

Setup: `demo/BW25113_2026.gb` plus a plasmid turned round the origin (File: rotate, or any GenBank with a feature across base 1).

### Window (nothing here has been run in a window)
1. **Start-up**: app opens, no import error (Synthesis, plannotate and data_table now load lazily or new). Status bar shows a small grey tip; welcome page shows slots.
2. **Tables** (Detect Common Features, Plasmid Checkup, Local Align, Folder Search): click a header to sort, then click a row: the selection must land on the
   *right* feature (rows are mapped by data number). Tick boxes in Detect and Checkup after sorting: "Add ticked features" must add the ticked ones.
   Right-click a header and a cell: show only / hide / show all rows. Local Align: click a row, the Alignment view follows, and back.
3. **Detect Common Features**: without pLannotate, "Install pLannotate…" appears and the info line says it is not installed; with it, "Found by" shows
   openDNA / pLannotate / both, the credit line shows, added features carry the right /note. Close the install dialog with Escape mid-run: the install stops.
4. **Synthesis tab**: open from Tools and from the right-click menu with a selection (bases land in the box, nothing runs). Run the calculator on a 500 bp and a
   2 kb piece: the window stays responsive (GIL fix), verdict and feature table show. Without ViennaRNA: a clear message and the install route works.
   IDT section: untick "Remember", enter fake credentials, the call fails with a readable error (the real service is untested); with a secure store the
   "Remember" box restores them next time.
5. **Circular origin in the window**: Find (Ctrl+F) a DNA query that runs through base 1; search scope "selection" over a selection through the origin;
   protein motif through the origin. Recode CDS across the origin: it works, and Undo takes **two** steps (known).
6. **Variants pane** with a BAM across a circular plasmid (demo BAM): a deletion through the origin is one row; "could not be checked" clips at the ends are listed.
7. **PCR panel** on a circular plasmid with a large window: "Search window narrowed…" line; more than the shown pairs: "Only the best N pairs are listed".
8. **Welcome page**: edit `sponsors.json` with an `until` in the past: the slot disappears; with a number instead of a date: **currently crashes** (review finding).

### Command line
9. Run every recipe example once on a turned plasmid: `skills show index`, then `inspect`, `region 900..50`, `sites`, `guides --top 1`, `repeats`, `search`, `align`,
   `pcr check/run` (window 0, then a bad pair: warning only), `kld`, `recombineering`, `gibson`, `codon` (CDS across origin), `variants --apply`, `checkup`
   (check the `cds` rows across the origin: known wrong), `split` then `merge` (names, clash), `msa` on a multi-FASTA, `synthesis --region 900..50` (end wrong: known),
   `wsltools` (with and without WSL), `ncbi` (network).
10. A cap warning shows in each command that has one (`--top 1`, a tiny `--max-hits`, NCBI more results).

### Core and packaging
11. `maturin develop --release`, restart: `import opendna_core` and `opendna_core.salis_repeats` exist. Build the installer/spec: `salis_forest.json.gz` (1.2 MB) and
    `sponsors.json` are bundled; `oDNA-cli.exe` runs `synthesis` on a fresh machine (ViennaRNA message, not a traceback).
12. Windows installer's optional WSL task runs `oDNA-cli wsltools install` (not exercised here).

## Before tagging
- Decide on the 10 open review findings (at least the three correctness ones), run `/simplify` again if code changes, rebuild the report, commit with `tests.json` and generated files (not `junit.xml`, `cargo.txt`).
- Version bump and release notes (`whats_new`), `THIRD_PARTY_NOTICES.md` for ViennaRNA, pLannotate, the Salis model; `PRIVACY.md` already edited: reread for IDT (credentials, network).
