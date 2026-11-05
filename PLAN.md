# OpenGenEdit — Development Plan

A local, offline editor for plasmid / construct work and large-genome browsing (plasmid editors and genome browsers cover similar ground).
Opens and edits `.gb/.gbk` and `.fa/.fasta` (+ `.fai`), opens `.dna` files read-only (saved as GenBank), with first-class annotations and primers.

---

## 1. Architecture decision

**Recommended: Rust core + Python (PySide6/Qt) GUI.**

```
┌──────────────────────────────────────────────────────────┐
│  GUI  (Python 3.12 + PySide6)                            │
│  linear view · circular map · track browser · tables     │
│  dialogs · primer design · undo/redo UI · plugins        │
└───────────────▲──────────────────────────────────────────┘
                │  PyO3 bindings (maturin wheel: genomeedit_core)
┌───────────────┴──────────────────────────────────────────┐
│  CORE (Rust, UI-agnostic)                                │
│  I/O: GenBank · .dna · FASTA/.fai · bgzip · GFF │
│  Sequence store: mmap + piece table (cheap edits)        │
│  Annotation index: interval tree w/ edit-shift rules     │
│  Edit engine: command log → undo/redo                    │
│  Search: IUPAC motifs, enzymes, primers (both strands)   │
│  Summaries: LOD tiles (GC%, feature density) for zoom    │
│  Bio utils: translate, revcomp, Tm, ORFs, digest         │
└──────────────────────────────────────────────────────────┘
```

**Why this split**
- Everything that touches every base (parsing, searching, tiling, editing a 250 Mb chromosome) lives in Rust → no Python strings of genome size, no GIL-bound loops.
- The GUI only ever asks for *what is visible* (`get_bases(chrom, start, end)`, `features_in(start, end)`, `tile(chrom, zoom, idx)`), so rendering cost is bounded by the screen, not the genome.
- Qt is much stronger than Godot for this domain's UI: tables, tree views, rich text, dialogs, clipboard, drag & drop, printing/SVG export, native file dialogs.

**Godot alternative (kept open, not recommended for v1):** because the core is UI-agnostic, a Godot 4 frontend could later bind the same crate via `godot-rust/gdext`. Godot's GPU canvas is excellent for zooming/panning, but you'd rebuild tables, forms and text editing by hand. Decide after Phase 2 if Qt rendering is ever the bottleneck (it shouldn't be with LOD tiles).

### Key libraries
| Concern | Choice |
|---|---|
| Rust ↔ Python | `pyo3` + `maturin` |
| GenBank | `gb-io` crate (fast, round-trips well) — wrap/extend |
| FASTA/.fai, BAM, VCF, GFF, bgzf | `noodles` crates |
| mmap | `memmap2` |
| Interval index | own augmented interval tree / `rust-lapper` for read-only tracks |
| Multi-pattern search | `aho-corasick`, `memchr` |
| `.dna` | own Rust reader only (format is a simple TLV packet stream; see §4) |
| Primer design / Tm | Own Rust code: SantaLucia NN Tm, hairpin/dimer ΔG, primer picking. **Not `primer3-py`** (GPL v2, would force the whole app to be GPL); optionally call an external `primer3_core` binary like bwa (§6b) |
| GUI | `PySide6` (QGraphicsView for maps, custom `QWidget.paintEvent` for sequence view) |
| Packaging | `maturin` wheel + PyInstaller (or Nuitka) → Windows installer |
| Validation in tests | `Biopython` as a reference parser |

---

## 2. Data model (core)

```
Document
 ├─ sequences: [Sequence]          # one for a plasmid, many for a genome FASTA
 │    ├─ id, topology (linear/circular), alphabet (DNA/RNA/protein)
 │    ├─ store: PieceTable(original=mmap|bytes, add_buffer)
 │    └─ features: IntervalIndex<Feature>
 ├─ primers: [Primer]              # name, seq, binding sites (computed), Tm, notes
 ├─ metadata: header fields, references, comments, .dna notes
 └─ history: CommandLog            # undo/redo stack, dirty flag
Feature: id, type, name, location (CompoundLocation: join/order, strand, fuzzy ends),
         qualifiers (ordered multimap), display (color, arrow, visible)
```

**Edit rules ( must be unit-tested exhaustively):**
- Insert inside a feature → feature grows. Insert at a boundary → configurable (default: doesn't grow).
- Delete fully covering a feature → feature removed (recorded for undo). Partial → truncated.
- Replace = delete + insert in one command.
- Circular sequences: features may wrap the origin; "set origin" rotates sequence + all coordinates.
- Primers are *not* stored as coordinates: binding sites are recomputed (debounced) after edits.

Coordinates: 0-based half-open internally everywhere; convert only at I/O and display (1-based).

---

## 3. Performance strategy for large genomes

Targets (on hg38, chr1 ≈ 248 Mb, plus a 3 Gb whole-genome FASTA):
- Open indexed FASTA: < 0.5 s (mmap + `.fai`; build `.fai` if missing, in background with progress).
- Pan/zoom at 60 fps at any zoom level.
- Motif search over chr1: < 1 s, streamed results.
- Single-base edit: O(log n) — no copying the chromosome.

Techniques:
1. **Never load whole sequences into Python.** The GUI requests windows.
2. **mmap + piece table**: original file is read-only mmap; edits go to an append buffer; save streams pieces out.
3. **Level-of-detail tiles**: precompute (lazily, cached on disk next to file as `.geidx`) per-bin summaries at zoom levels ×4 apart: GC%, N-runs, feature counts. Zoomed out = draw tiles; zoomed in = draw bases.
4. **Feature rendering**: query interval index for visible range; at low zoom, collapse to density histogram; pack rows (IGV-style "expanded/squished/collapsed").
5. **Background work** (indexing, searching, digest) runs in Rust threads (`rayon`), releases the GIL, reports progress via callback/polling; GUI stays responsive.
6. Small plasmid files (< ~1 Mb) are simply loaded fully into memory — same API.

---

## 4. File formats

| Format | Read | Write | Notes |
|---|---|---|---|
| GenBank `.gb/.gbk/.gbff` | ✅ | ✅ | Primary native save format. Preserve unknown qualifiers & header order for lossless round-trip. Multi-record files supported. |
| FASTA `.fa/.fasta/.fna` (+`.gz`/bgzip) | ✅ | ✅ | Indexed random access via `.fai`/`.gzi`. Annotations stored in a sidecar (GFF3) since FASTA can't hold them. |
| `.dna` | ✅ | ❌ by design | Binary TLV stream: `[type:u8][len:u32 BE][payload]`. Used packets: 0x09 cookie, 0x00 DNA (flags byte + seq), 0x0A features (XML), 0x05 primers (XML → `primer_bind` features), 0x06 notes. History, caches and display packets are ignored. Edits are saved as GenBank. |
| GFF3 / BED | ✅ | ✅ | Annotation tracks for genomes; sidecar for FASTA. |
| BAM/CRAM, VCF | ✅ (Phase 5) | – | IGV-like read/variant tracks via `noodles`. |
| `.ab1` Sanger traces | Phase 6 | – | Chromatogram view + alignment to construct. |

Golden rule: **open → save without edits must be byte-identical or semantically identical** (checked by test corpus).

---

## 5. GUI design (PySide6)

Main window = tabbed documents, each with switchable/splittable views sharing one selection model:

1. **Sequence view** (standard layout): wrapped lines of bases, complementary strand, 3/6-frame translation, features as arrows above/below, primers drawn at binding sites, enzyme sites on top. Click/drag selection, keyboard typing to insert/overwrite, copy (seq / revcomp / FASTA / GenBank fragment), paste with features.
2. **Map view**: circular plasmid map (QGraphicsView) for circular sequences; linear map otherwise. Click feature → selects in all views.
3. **Genome browser view** (IGV-like): chromosome dropdown, locus box (`chr1:1,000,000-1,050,000` or gene name search), ruler, tracks (sequence, GC, features, GFF/BED, later BAM/VCF), drag-pan, wheel-zoom.
4. **Side panels**: Features table (sortable, editable, filter by type), Primers table, Enzymes panel (choose set: common/6-cutters/custom, unique cutters highlight), History panel, Properties/notes.
5. **Status bar**: selection start–end, length, GC%, Tm of selection, cursor position.

Rendering rule: each view asks the core for its visible window only; repaint on `document_changed(range)` signal.

---

## 6. Feature set by phase

### Phase 0 — Skeleton (1 week)
- Repo layout, maturin build, PyO3 hello-world, CI (GitHub Actions: Windows + Linux), pytest + `cargo test`, test-file corpus folder.

### Phase 1 — Core I/O & model (3–4 weeks)
- Rust: FASTA/.fai (mmap), GenBank read/write, `.dna` read, Document/Feature/Location model, piece table, interval index, edit commands with feature-shift rules, undo/redo.
- ✅ **Done** (Oct 2026): GenBank read/write (gb-io, vendored with a writer fix), `.dna` read, edits with feature-shift rules, undo/redo, Python `Document` API, tests on real files (E. coli genome opens in ~0.1 s, saves in ~0.1 s, byte-stable re-saves).
  Deferred on purpose: piece table (a `Vec<u8>` splice is fine up to chromosome size; needed only for editing multi-GB FASTA) and interval index (linear scan is µs for 10k features; needed for Phase 4 tracks).
- Python API: `open()`, `save()`, `get_bases()`, `features_in()`, `apply(edit)`, `undo()`.
- Tests: round-trip corpus, property-based tests (`proptest`) for edit/shift invariants, benchmarks.
- **Milestone:** CLI can open a `.dna`, insert bases, save `.gb`, and Biopython agrees on the result.

### Phase 2 — Plasmid editor MVP (4–5 weeks) → *usable for daily work*
- Sequence view with selection, typing, copy/paste, revcomp, translation.
- Circular/linear map.
- Feature table + add/edit/delete feature dialog (type, name, strand, color, qualifiers).
- Restriction enzymes: a list of common enzymes, site search, unique/2-cutter filter, show on map (a larger catalogue needs a source we may redistribute: see docs/TODO.md).
- Find (IUPAC motif, both strands, mismatches), go-to position.
- Undo/redo, dirty tracking, autosave/recovery file.
- ✅ **Done** (part 1): sequence view (wrapped rows, complement, feature lanes with arrows, click/drag/shift selection, typing IUPAC bases, backspace/delete, copy/paste), circular + linear map (cached rendering, labels without overlap, click to select), feature table (sort, filter, color, synced selection), open/save as GenBank, undo/redo, unsaved-changes prompt, status bar with GC%.
  Genome: open + first paint 0.4 s; typing one base ~0.3 s (feature map re-render).
- ✅ **Done** (part 2): feature add/edit/delete (dialog: name, type, GenBank location, strand, color, qualifiers; undoable in Rust), 47 built-in enzymes with unique cutters on map + sequence (both strands, across the origin), find (DNA both strands + feature names), go-to, 3-frame translation, copy reverse complement, whole-sequence reverse-complement, autosave every 60 s to `<file>.autosave.gb` with recovery prompt, per-row height (rows only reserve the lanes they use).
  Genome: open + first paint 0.6 s; typing one base ~0.4 s.
- **Left from the Phase 2 list:** IUPAC/mismatch search, 2-cutter filter and enzyme picker, genome cut sites (skipped above 1 Mb), selection-only reverse-complement, faster genome typing (debounce map redraw / incremental caches).
- **Milestone:** open/edit/save your real plasmid files without another editor.

### Phase 3 — Primers & bio tools (3 weeks)
- Primer object: add from selection, import list (CSV/FASTA), binding-site detection with mismatches + 5′ tails, Tm (SantaLucia NN, adjustable salt/primer conc.), hairpin/dimer ΔG (own Rust implementation, nearest-neighbor tables from published papers).
- Primer design: pick primers flanking selection (own scorer: Tm window, GC clamp, ΔG, uniqueness), sequencing primer walker. Optional "Use primer3" via external `primer3_core` process.
- PCR product simulation → new document with features carried over.
- ORF finder, translate selection, digest simulation (fragment list + virtual gel).
- ✅ **Done** (part 1): Tm in Rust (SantaLucia 1998 NN, salt + Mg²⁺ correction; matches Biopython `Tm_NN` to 0.01 °C), adjustable PCR conditions, Tm of the selection in the status bar, binding-site search (both strands, origin, 15-nt 3′ seed with up to 3 mismatches except the last 3 bases, 5′ tails), Primers tab (name, sequence, length, Tm, GC of the annealed part), new primer from selection or sequence, import FASTA/CSV/TSV list. Primers are stored as `primer_bind` features with `/note="sequence: …"` (the usual GenBank convention), so a primer that binds nowhere is not kept.
  Also: copy/cut/paste carries features (clipped), typing limited to ACGTN, paste strips FASTA/GenBank formatting and refuses non-DNA, readable fonts and black/white labels by background. Fixed: primer sites read from .dna files were 1 base off (binding sites are 0-based in .dna).
- ✅ **Done** (part 1b): primers are names + sequences in a primer database (`primers.tsv` in the app data folder; unique names; TSV/CSV/FASTA import/export; primers found in opened files are added). Their binding sites are **computed live** on the open sequence after every edit (all database primers; exact sites and, per Options, sites with up to N mismatches in the 3′ 15-nt seed; >300 kb matched in a background thread). Drawn as their own track: thin line on the annealed bases, 5′ tail raised, 3′ half arrow, dashed orange + red marks for mismatches; forward above the bases, reverse below. Primers tab = tree grouped by primer, collapsed by default. Primer dialog shows live where a sequence binds plus 3′ mis-priming hits here and in a genome (default BW25113 file next to the app). `primer_bind` features stored in files are hidden; on save (and autosave) they are replaced by the current exact sites of all database primers (file primers whose name means another sequence in the database are kept as they were).
- ✅ **Done** (part 2): hairpin / self- and cross-dimer ΔG in Rust (SantaLucia NN; perfectly paired stems only, a quick check, not a full fold), shown in the primer dialog and Primers tab with the usual warning limits; primer-pair design around a selection (own penalty: Tm, GC, 3′ clamp, runs, ΔG, 3′ uniqueness, pair Tm difference and cross-dimer; ~0.2 s on a plasmid), saved to the database as <name>_F/_R; PCR product → new untitled window (primer tails added, mismatched primer bases written into the product, features carried over); ORF finder (6 frames, add as CDS); restriction digest with fragment table and virtual gel; primer sites on the map.
- **Left:** bundled/downloadable BW25113 reference (with packaging, Phase 6); full thermodynamic folding (thal-like) if the quick ΔG is not enough; sequencing-primer walker; ORFs across the origin.

### Phase 4 — Large-genome browser (3–4 weeks)
- Genome view with LOD tiles and `.geidx` cache, chromosome navigation, GFF3/BED tracks, gene-name search index.
- Background indexing with progress; streamed genome-wide search.
- Edits on big genomes (piece table) + streaming save.
- **Milestone:** hg38 chr1 opens instantly, smooth zoom from whole chromosome to bases.
- ✅ **Done**: IGV-style genome view above the editor (overview, ruler, GC %, bases at deep zoom, features, primers, loaded tracks; wheel zoom, drag pan, Shift+drag select, click to select; locus box takes a range, a position or a name), synced with the sequence view and map. FASTA references open in the editor (sequence picker for multi-FASTA; save streams the file, replacing only the edited sequence). GFF3/GTF/BED tracks in Rust with an interval index (sorted starts + running max end) and density plots when crowded. Find/locus box search gene, locus_tag, product… and track names. Genome-wide motif search streamed per sequence. Big loads/saves run in a background thread behind a busy dialog. Sequence-view layout is sparse (only taller rows stored), and typing reuses feature lanes/colors.
  Measured on a synthetic 250 Mb chromosome: open (index + load + first paint) 0.9 s, whole-chromosome genome view 0.13 s, zoom/jump 0.01 s, typing one base 0.17 s, streaming save 0.7 s. E. coli genome: typing 0.15–0.23 s (was 0.45 s).
  Not built, by measurement: piece table (a `Vec<u8>` edit on 250 Mb takes ~10 ms) and an on-disk `.geidx` cache (GC per 1 kb block is computed in 50 ms and kept in memory; .fai covers random access).
  Left: multi-record GenBank files (still refused), IUPAC/mismatch genome search, saving tracks with the document.

### Phase 5 — Cloning & IGV extras (4+ weeks, pick by need)
- Cloning simulation: restriction/ligation, Gibson/HiFi overlap, Golden Gate (Type IIS), with history lineage ("made from").
- BAM/VCF tracks (coverage + reads at high zoom, variants).
- Pairwise alignment of a sequence/read to the construct (Rust `bio` crate / parasail), `.ab1` trace view.
- Export: SVG/PNG/PDF of maps, FASTA/GenBank of selection.
- ✅ **Done** (part 1, IGV extras): BAM alignments in the genome view (Rust, `noodles`): a coverage row (exact under 30 kb; wider windows from per-128-bp block totals computed once per sequence in a background thread) and, under 30 kb (IGV's visibility window), reads packed in 20 rows, coloured by strand, MAPQ 0 faint, mismatches as coloured bases, insertions as purple bars, deletions as a thin line; click a read for its name/position/MAPQ. Unmapped, duplicate and QC-failed reads are ignored. A BAM without `.bai` is indexed on open (if coordinate-sorted). VCF (plain or bgzipped) loads as a variant track (`ID REF>ALT`). External mappers (§6b): `app/tools.toml` recipes (bwa mem, minimap2 short/ONT/HiFi), run in WSL on Windows with a live log and Cancel (stops the WSL processes too), one run folder with reference, commands and log; the BAM opens as a track when done. Tools → Install Mapping Tools creates a `genomeedit-tools` conda env (bioconda) in WSL; Tools → Check Mapping Tools lists what is found.
- ✅ **Done** (part 2): copy-paste keeps annotations across windows/app instances (system clipboard); workspace tabs (any view or tool in any split pane, layout saved); gene groups (gene feature + parts sharing `/locus_tag`, collapsible in the feature table, one arrow on the map); Detect Common Features (Rust k-mer + protein matcher, NCBI-built `app/parts.tsv` + user library); Tools → Settings (per-tool parameters). Rust for IGV-style reads: pairing info per read, `Bam.pileup`, `Bam.variants` (SNP/indel scan).
- ✅ **Done** (part 3, IGV styling): reads coloured by pair orientation and insert size (IGV rules: grey normal, red larger / blue smaller than expected, teal LL, indigo RR, green RL, orange mate on another sequence, purple outline mate unmapped; limits = 0.5–99.5th percentile of proper pairs in view, or fixed), or by strand. Coverage split by base where reads disagree above the variant threshold; consensus changes (SNP, deletion, insertion) from a background scan flagged with red marks at any zoom. Variants pane: list, go to, apply ticked changes to the sequence (e.g. correct a plasmid map after whole-plasmid sequencing). Alignments settings: colour mode, insert limits, reads window, variant 20 % / consensus 80 %, min depth, base and mapping quality. Fixed: mate reads were never passed to paired-end mappers.
  Left for IGV extras: coloured coverage bars at variant positions (allele counts), downsampling, paired-read/insert-size colouring, CRAM, read sorting/grouping, tracks that follow edits.

### Phase 6 — Polish & distribution
- Preferences, themes (light/dark), keyboard shortcuts, recent files, file associations on Windows (`.gb`, `.dna`, `.fasta`).
- Plugin hook: Python scripts get the Document API (batch annotation, custom tools).
- Installer (PyInstaller + Inno Setup), crash log, user docs.

Rough total to a strong daily-driver (Phases 0–4): **~3–4 months** solo part-time-heavy; Phase 2 alone is the first big payoff.

---

## 6b. External command-line tools (bwa, samtools, minimap2, …), planned for Phase 5

**Principle:** don't re-implement aligners. OpenGenEdit *runs* them as separate processes and *reads their outputs* (BAM/VCF/BED) with the Rust core, exactly like a shell pipeline.

```
GUI "Align reads…" dialog
   │ fills a tool recipe (inputs, params)
   ▼
Python ToolRunner (QProcess: async, live log, cancel button)
   │ launches:  wsl -e bash -lc "bwa mem ref.fa r1.fq r2.fq | samtools sort -o out.bam && samtools index out.bam"
   ▼
output files ──► Rust core opens out.bam (noodles) ──► new track in the genome view
```

1. **Tool recipes in a plain file** (`tools/tools.toml`), one entry per tool: command template, inputs, outputs, and what to do with the result. Adding a tool = adding a recipe, no code:
   ```toml
   [bwa_mem]
   label   = "Map reads (bwa mem)"
   needs   = ["bwa", "samtools"]
   prepare = "bwa index {ref}"            # run only if {ref}.bwt is missing
   run     = "bwa mem -t {threads} {ref} {reads1} {reads2} | samtools sort -o {out}.bam && samtools index {out}.bam"
   result  = { track = "{out}.bam" }
   ```
2. **Runner backends**: bwa/samtools don't build natively on Windows, so:
   - **WSL** (default on Windows; Ubuntu is already installed here): prefix `wsl -e bash -lc`, and translate paths `V:\data\x.fa` → `/mnt/v/data/x.fa`.
   - **Native** on Linux/macOS, or Windows builds where they exist (minimap2).
   - **Docker** optional later (biocontainers images) for pinned versions.
3. **Install tools** inside WSL with conda/mamba from bioconda (`mamba install -c bioconda bwa samtools minimap2`); a *Tools → Check setup* menu runs `which` for every `needs` entry and shows what's missing.
4. **Safety**: arguments are quoted (`shlex.quote`), never built from raw user text; outputs go to a per-run folder with the full command + log saved next to them (reproducibility).
5. **Why Python, not Rust, runs them**: it's orchestration (processes, logs, dialogs), not per-base work, and Qt's `QProcess` already gives async I/O + cancel. Rust's job starts when the output file needs to be read fast.

---

## 7. Repository layout

```
OpenGenEdit/
├─ core/                    # Rust crate (cdylib + rlib)
│  ├─ src/io/{fasta,genbank,dna_file,gff}.rs
│  ├─ src/model/{document,feature,location}.rs
│  ├─ src/store/{piece_table,interval_index,tiles}.rs
│  ├─ src/edit/{commands,history}.rs
│  ├─ src/bio/{translate,tm,enzymes,search,orf}.rs
│  ├─ src/py.rs             # PyO3 bindings only
│  └─ benches/
├─ app/                     # Python package: genomeedit
│  ├─ main.py
│  ├─ views/{sequence_view,map_view,genome_view}.py
│  ├─ panels/{features,primers,enzymes,history}.py
│  ├─ dialogs/
│  └─ resources/{enzymes.json,icons,themes}
├─ tests/
│  ├─ corpus/               # real .gb/.dna/.fasta samples (+ expected outputs)
│  ├─ test_roundtrip.py
│  └─ test_edits.py
├─ pyproject.toml           # maturin build backend
└─ PLAN.md
```

---

## 8. Risks & mitigations

| Risk | Mitigation |
|---|---|
| `.dna` is a proprietary, undocumented format | Read-only; edits always saved as GenBank, so a `.dna` file is never modified. |
| GenBank edge cases (fuzzy locations, `join`/`order`, wrap-around, odd qualifiers) | Use `gb-io`, keep unknown data verbatim, corpus from NCBI/Addgene. |
| Feature-shift rules subtly wrong after edits | Property tests: random edit sequences + undo must restore exact state. |
| Qt rendering too slow at base-level for huge windows | Draw only visible rows; cache glyph pixmaps per base; LOD switch threshold; fallback to QOpenGLWidget. |
| Scope creep (it's two big products) | Phases are ordered so each ends in something you use daily; Phase 5 items are opt-in. |

---

## 9. Immediate next steps

1. Collect a test corpus: 20–50 of your own `.dna`/`.gb` files + a reference genome FASTA (e.g. hg38 chr1, E. coli, yeast).
2. Phase 0 scaffold: Rust crate + PyO3 + PySide6 window that opens a FASTA and prints length.
3. Implement GenBank + `.dna` readers and the round-trip test harness before any GUI work.
