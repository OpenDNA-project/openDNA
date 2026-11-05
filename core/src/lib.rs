//! Python bindings: the only file that knows about Python.
//! All real logic lives in plain Rust modules (fasta.rs, ...) so it stays testable
//! and could later be bound to another frontend.

// `mod fasta;` = "there is a module in fasta.rs" (Rust does not discover files on its own).
pub mod align;
pub mod alignments;
pub mod annotate;
pub mod bio;
pub mod code;
pub mod cas9;
pub mod document;
pub mod fasta;
pub mod locations;
pub mod primers;
pub mod search;
pub mod dna_file;
pub mod tracks;

use gb_io::seq::Topology;
use pyo3::prelude::*;
use std::path::PathBuf;

/// Python class `genomeedit_core.Document(path)`: an editable GenBank or .dna file.
/// All positions are 0-based, end-exclusive (Python slicing convention).
#[pyclass]
struct Document {
    inner: document::Document,
}

/// Read-only snapshot of one feature, handed to Python.
/// `get_all` makes every field readable as a Python attribute (`f.kind`, `f.spans`, ...).
#[pyclass(get_all, frozen)]
struct Feature {
    /// Position in the document's feature list (stable until the next edit).
    index: usize,
    kind: String,
    label: String,
    /// GenBank location text, 1-based (e.g. "complement(join(10..20,30..40))").
    location: String,
    /// The base ranges in order, 0-based end-exclusive.
    spans: Vec<(i64, i64)>,
    /// +1 forward, -1 reverse.
    strand: i8,
    qualifiers: Vec<(String, Option<String>)>,
}

impl Feature {
    fn from_gb(index: usize, f: &gb_io::seq::Feature) -> Self {
        let (spans, strand) = locations::spans(&f.location);
        // First available of these qualifiers names the feature; fall back to its type.
        let label = ["label", "gene", "locus_tag", "product", "note"]
            .iter()
            .find_map(|k| f.qualifier_values(k).next())
            .unwrap_or(&f.kind)
            .to_string();
        Feature {
            index,
            kind: f.kind.to_string(),
            label,
            location: f.location.to_gb_format(),
            spans,
            strand,
            qualifiers: f.qualifiers.iter().map(|(k, v)| (k.to_string(), v.clone())).collect(),
        }
    }
}

#[pymethods]
impl Feature {
    fn __repr__(&self) -> String {
        format!("<Feature {} '{}' {}>", self.kind, self.label, self.location)
    }
}

#[pymethods]
impl Document {
    #[new]
    fn new(py: Python<'_>, path: PathBuf) -> PyResult<Self> {
        let inner = py.allow_threads(|| document::Document::open(&path))?;
        Ok(Document { inner })
    }

    /// `Document.from_fasta(path, name)`: one sequence of a (multi-)FASTA file.
    /// File → New: `Document.new_sequence("pMine", "ACGT...", circular=True)`, unsaved.
    #[staticmethod]
    #[pyo3(signature = (name, bases="", circular=false))]
    fn new_sequence(name: &str, bases: &str, circular: bool) -> PyResult<Self> {
        Ok(Document { inner: document::Document::new_sequence(name, bases.as_bytes(), circular)? })
    }

    /// `#[staticmethod]` = called on the class, like Python's @staticmethod.
    #[staticmethod]
    fn from_fasta(py: Python<'_>, path: PathBuf, name: &str) -> PyResult<Self> {
        let inner = py.allow_threads(|| document::Document::from_fasta(&path, name))?;
        Ok(Document { inner })
    }

    /// GC fraction per bin over [start, end) (NaN = empty bin).
    fn gc_bins(&self, py: Python<'_>, start: usize, end: usize, bins: usize) -> Vec<f32> {
        let inner = &self.inner;
        py.allow_threads(|| inner.gc_bins(start, end, bins))
    }

    /// Write FASTA `dst`: a copy of `src` with this sequence replacing its record.
    /// Does not mark the document saved (set `doc.modified = False` after it succeeds):
    /// it only reads the document, so the window can repaint while a big file is written.
    fn write_fasta(&self, py: Python<'_>, src: PathBuf, dst: PathBuf) -> PyResult<()> {
        let inner = &self.inner;
        Ok(py.allow_threads(|| inner.write_fasta(&src, &dst))?)
    }

    // `#[getter]` = read-only Python property: `doc.name` (no parentheses).
    #[getter]
    fn name(&self) -> String {
        self.inner.record.name.clone().unwrap_or_default()
    }

    #[getter]
    fn definition(&self) -> Option<String> {
        self.inner.record.definition.clone()
    }

    #[getter]
    fn length(&self) -> usize {
        self.inner.len()
    }

    #[getter]
    fn circular(&self) -> bool {
        self.inner.record.topology == Topology::Circular
    }

    /// `doc.rotate(k)`: base k (0-based) becomes base 1, features along; one undo step.
    fn rotate(&mut self, k: usize) -> PyResult<()> {
        Ok(self.inner.rotate(k)?)
    }

    /// `doc.linearize(k)`: make a circular sequence linear, broken just before base k
    /// (0-based; 0 keeps the start). Rotation and topology change are ONE undo step.
    fn linearize(&mut self, k: usize) -> PyResult<()> {
        Ok(self.inner.linearize(k)?)
    }

    /// `doc.set_topology(circular)`: circular or linear as one undo step (unlike the
    /// `circular` property, which is not part of the history).
    fn set_topology(&mut self, circular: bool) {
        self.inner.set_topology(circular)
    }

    // `doc.circular = False` from Python: the setter pairs with the getter above.
    #[setter]
    fn set_circular(&mut self, circular: bool) {
        self.inner.set_circular(circular)
    }

    #[getter]
    fn modified(&self) -> bool {
        self.inner.is_modified()
    }

    /// `doc.edited_ranges()`: [(start, end)] (0-based, end exclusive) of the bases typed,
    /// pasted, replaced or inserted since the last save or load; the sequence view paints
    /// them yellow. Undo, redo and every later edit keep the list right; a save empties it.
    fn edited_ranges(&self) -> Vec<(usize, usize)> {
        self.inner.edited_ranges().to_vec() // `.to_vec()` copies the slice into a Vec for Python
    }

    // `#[setter]` = `doc.modified = True` from Python.
    #[setter]
    fn set_modified(&mut self, modified: bool) {
        self.inner.set_modified(modified)
    }

    /// `doc.sequence()` whole sequence, or `doc.sequence(start, end)`.
    /// `#[pyo3(signature = ...)]` declares Python default arguments.
    #[pyo3(signature = (start=0, end=None))]
    fn sequence(&self, start: usize, end: Option<usize>) -> String {
        let seq = &self.inner.record.seq;
        let end = end.unwrap_or(seq.len()).min(seq.len());
        let start = start.min(end);
        String::from_utf8_lossy(&seq[start..end]).into_owned()
    }

    fn features(&self) -> Vec<Feature> {
        self.inner.record.features.iter().enumerate().map(|(i, f)| Feature::from_gb(i, f)).collect()
    }

    /// Features overlapping [start, end).
    fn features_in(&self, start: i64, end: i64) -> Vec<Feature> {
        let all = &self.inner.record.features;
        self.inner.features_in(start, end).into_iter().map(|i| Feature::from_gb(i, &all[i])).collect()
    }

    /// Primers from a .dna file as (name, sequence). Their binding sites are also
    /// present as `primer_bind` features.
    fn primers(&self) -> Vec<(String, String)> {
        self.inner.primers.iter().map(|p| (p.name.clone(), p.sequence.clone())).collect()
    }

    fn insert(&mut self, pos: usize, bases: &str) -> PyResult<()> {
        Ok(self.inner.replace(pos, pos, bases.as_bytes())?)
    }

    fn delete(&mut self, start: usize, end: usize) -> PyResult<()> {
        Ok(self.inner.replace(start, end, b"")?)
    }

    fn replace(&mut self, start: usize, end: usize, bases: &str) -> PyResult<()> {
        Ok(self.inner.replace(start, end, bases.as_bytes())?)
    }

    /// Replace feature `index` (or add one when `index` is None) and return its index.
    /// `location` is GenBank text, 1-based: "10..200", "complement(join(1..5,9..20))".
    #[pyo3(signature = (index, kind, location, qualifiers))]
    fn set_feature(
        &mut self,
        index: Option<usize>,
        kind: String,
        location: &str,
        qualifiers: Vec<(String, Option<String>)>,
    ) -> PyResult<usize> {
        Ok(self.inner.set_feature(index, to_gb((kind, location.to_string(), qualifiers))?)?)
    }

    /// Add several features as one undo step: a list of (kind, location, qualifiers).
    /// Returns the index of the first one.
    fn add_features(&mut self, features: Vec<FeatureTuple>) -> PyResult<usize> {
        let features = features.into_iter().map(to_gb).collect::<PyResult<Vec<_>>>()?;
        Ok(self.inner.add_features(features)?)
    }

    /// Features overlapping [start, end), clipped to it, located relative to `start`
    /// (what Copy keeps so Paste can bring the annotations along).
    /// On circular records `end` may pass the origin (up to one full turn).
    fn copy_features(&self, start: usize, end: usize) -> Vec<Feature> {
        let n = self.inner.len();
        let found = if self.circular() && start < n && end > n {
            self.inner.features_across(start, end.min(start + n))
        } else {
            let end = end.min(n);
            self.inner.features_between(start.min(end), end)
        };
        found.iter().enumerate().map(|(i, f)| Feature::from_gb(i, f)).collect()
    }

    /// `copy_features` for "Copy Reverse Complement": the same clipped features, flipped onto
    /// the reverse-complemented copy (strand inverted, join parts reordered).
    fn copy_features_reversed(&self, start: usize, end: usize) -> Vec<Feature> {
        let n = self.inner.len();
        let end = if self.circular() && start < n { end.min(start + n) } else { end.min(n) };
        let start = start.min(end);
        self.inner
            .features_between_reversed(start, end)
            .iter()
            .enumerate()
            .map(|(i, f)| Feature::from_gb(i, f))
            .collect()
    }

    /// Replace [start, end) passing the origin of a circular record (length < end <=
    /// start + length) with `bases`, plus pasted `features` (as for `paste`) if given.
    /// One undo step. Returns where the new bases start.
    #[pyo3(signature = (start, end, bases, features=None))]
    fn edit_across_origin(&mut self, start: usize, end: usize, bases: &str, features: Option<Vec<FeatureTuple>>) -> PyResult<usize> {
        // Bases are checked in Rust by the edits themselves (`check_edit`), like `replace`.
        // `.map(...).transpose()?` turns Option<Result<..>> into Result<Option<..>>, then unpacks it.
        let features = features
            .map(|fs| fs.into_iter().map(to_gb).collect::<PyResult<Vec<_>>>())
            .transpose()?;
        Ok(self.inner.edit_across_origin(start, end, bases.as_bytes(), features)?)
    }

    /// Replace [start, end) with `bases` plus features (kind, location relative to the
    /// pasted bases, qualifiers), as one undo step.
    fn paste(&mut self, start: usize, end: usize, bases: &str, features: Vec<FeatureTuple>) -> PyResult<()> {
        let features = features.into_iter().map(to_gb).collect::<PyResult<Vec<_>>>()?;
        Ok(self.inner.paste(start, end, bases.as_bytes(), features)?)
    }

    /// Where `primer` anneals: list of (start, end, strand, mismatches), top-strand
    /// coordinates; on circular records `end` can be > length (wraps through the origin).
    /// `seed` = how many 3'-end bases must match; use ~12 to look for mis-priming.
    #[pyo3(signature = (primer, max_mismatches=1, seed=primers::SEED))]
    fn binding_sites(&self, py: Python<'_>, primer: &str, max_mismatches: usize, seed: usize) -> Vec<(usize, usize, i8, usize)> {
        let record = &self.inner.record;
        let circular = record.topology == Topology::Circular;
        let sites = py.allow_threads(|| primers::binding_sites(&record.seq, primer.as_bytes(), circular, seed, max_mismatches));
        sites
            .into_iter()
            .map(|s| (s.start, s.start + s.len, if s.forward { 1 } else { -1 }, s.mismatches))
            .collect()
    }

    fn delete_feature(&mut self, index: usize) -> PyResult<()> {
        Ok(self.inner.delete_feature(index)?)
    }

    /// Several feature changes as one undo step: a list of (index or None, feature tuple
    /// or None). (i, f) replaces #i, (None, f) adds f, (i, None) deletes #i. Indexes are
    /// those before the call.
    fn edit_features(&mut self, changes: Vec<(Option<usize>, Option<FeatureTuple>)>) -> PyResult<()> {
        let changes = changes
            .into_iter()
            .map(|(i, f)| Ok((i, f.map(to_gb).transpose()?))) // Option<PyResult> -> PyResult<Option>
            .collect::<PyResult<Vec<_>>>()?;
        Ok(self.inner.edit_features(changes)?)
    }

    /// Substitute bases at `pos` (same length; features do not move). Undoable.
    fn overwrite(&mut self, pos: usize, bases: &str) -> PyResult<()> {
        Ok(self.inner.overwrite(pos, bases.as_bytes())?)
    }

    /// PCR product as a new unsaved Document: prefix + bases [start, end) + suffix
    /// (end may pass the origin of a circular template). Features come along, clipped.
    fn pcr_product(&self, name: &str, start: usize, end: usize, prefix: &str, suffix: &str) -> PyResult<Document> {
        let inner = self.inner.pcr_product(name, start, end, prefix.as_bytes(), suffix.as_bytes())?;
        Ok(Document { inner })
    }

    /// Open reading frames of at least `min_codons` codons: list of (start, end, strand),
    /// end includes the stop codon; on a circular sequence an ORF across the origin has `end` past the length.
    /// `table`: NCBI genetic code (stops are its stops; unknown ids fall back to 1). `alternative_starts`: also start
    /// at the table's other start codons (GTG, TTG, ...), not only ATG.
    #[pyo3(signature = (min_codons=100, table=1, alternative_starts=false))]
    fn orfs(&self, py: Python<'_>, min_codons: usize, table: u8, alternative_starts: bool) -> Vec<(usize, usize, i8)> {
        let seq = &self.inner.record.seq;
        let circular = self.inner.record.topology == Topology::Circular;
        let code = code::by_id_or_standard(table);
        py.allow_threads(|| bio::find_orfs_with(seq, min_codons, circular, code, alternative_starts))
            .into_iter()
            .map(|o| (o.start, o.end, if o.forward { 1 } else { -1 }))
            .collect()
    }

    /// How CDS feature `index` is translated: (NCBI table id, frame 0-2 = /codon_start - 1, first codon may be read as
    /// Met when it is a start codon (the 5' end is complete and the frame is 0), [(codon number, amino acid letter or
    /// "*")] from /transl_except, why it cannot be (empty when fine)). `default_table` is used when the feature has
    /// no /transl_table. See code.rs `plan_cds`.
    #[pyo3(signature = (index, default_table=1))]
    fn cds_plan(&self, index: usize, default_table: u8) -> PyResult<(u8, usize, bool, Vec<(usize, String)>, String)> {
        let f = self
            .inner
            .record
            .features
            .get(index)
            .ok_or_else(|| pyo3::exceptions::PyIndexError::new_err(format!("no feature #{index}")))?;
        Ok(match code::plan_cds(f, self.inner.record.seq.len(), default_table) {
            Ok(p) => (
                p.table.id,
                p.frame,
                !p.partial5 && p.frame == 0,
                p.exceptions.iter().map(|&(c, aa)| (c, (aa as char).to_string())).collect(),
                String::new(),
            ),
            // A plan that cannot be made: draw with the default table so something sensible shows, and say why.
            Err(why) => (code::by_id_or_standard(default_table).id, 0, false, Vec::new(), why),
        })
    }

    /// Recompute the /translation of CDS features from their bases (one undo step; see
    /// document.rs `refresh_translations`). `indexes`: which features, or None for all. Returns one
    /// (feature index, "updated" | "unchanged" | "pseudo" | "partial" | "other", reason) per CDS that has a /translation.
    #[pyo3(signature = (indexes=None, default_table=1))]
    fn refresh_translations(&mut self, indexes: Option<Vec<usize>>, default_table: u8) -> Vec<(usize, String, String)> {
        use code::Skip;
        use document::Refresh;
        self.inner
            .refresh_translations(indexes.as_deref(), default_table)
            .into_iter()
            .map(|(i, r)| match r {
                Refresh::Updated => (i, "updated".to_string(), String::new()),
                Refresh::Unchanged => (i, "unchanged".to_string(), String::new()),
                Refresh::Skipped(kind, why) => {
                    let kind = match kind {
                        Skip::Pseudo => "pseudo",
                        Skip::Partial => "partial",
                        Skip::Other => "other",
                    };
                    (i, kind.to_string(), why)
                }
            })
            .collect()
    }

    /// Flip the whole record (sequence + features) to the other strand. Undoable.
    fn reverse_complement(&mut self) {
        self.inner.reverse_complement()
    }

    fn undo(&mut self) -> bool {
        self.inner.undo()
    }

    fn redo(&mut self) -> bool {
        self.inner.redo()
    }

    /// Add an undo step of the app's own that changes nothing in the document (a primer
    /// database change): in the history as `label` (its note is "external", its removed
    /// text is `token`), never makes the document modified. Undo / Redo pass over it like
    /// any step; the app asks `next_external` first and does the real work itself.
    fn external(&mut self, label: &str, token: &str) {
        self.inner.external(label, token)
    }

    /// The token of the next step to undo (`redo` false) or redo, if it is an external one.
    #[pyo3(signature = (redo=false))]
    fn next_external(&self, redo: bool) -> Option<String> {
        self.inner.next_external(redo)
    }

    /// The edit history for the History panel: (undo steps oldest first, redo steps next
    /// first), each (pos, removed bases (first 30), removed count, inserted bases (first 30),
    /// inserted count, added feature names, deleted feature names, note).
    fn history(&self) -> (Vec<HistoryStep>, Vec<HistoryStep>) {
        let (undo, redo) = self.inner.history();
        let tuple = |s: document::Summary| {
            (s.pos, s.removed, s.removed_len, s.inserted, s.inserted_len, s.added, s.deleted, s.note.to_string())
        };
        (undo.into_iter().map(tuple).collect(), redo.into_iter().map(tuple).collect())
    }

    /// `primers`: list of (kind, location, qualifiers) written instead of the document's
    /// own primer_bind features (the app passes the current primer sites).
    #[pyo3(signature = (path, primers=None))]
    fn save_genbank(&mut self, path: PathBuf, primers: Option<Vec<FeatureTuple>>) -> PyResult<()> {
        let primers = primers.map(|p| p.into_iter().map(to_gb).collect::<PyResult<Vec<_>>>()).transpose()?;
        Ok(self.inner.save_genbank(&path, primers)?)
    }

    /// Write a GenBank copy without marking the document saved (autosave).
    #[pyo3(signature = (path, primers=None))]
    fn write_genbank(&self, path: PathBuf, primers: Option<Vec<FeatureTuple>>) -> PyResult<()> {
        // `.transpose()` turns Option<Result<..>> into Result<Option<..>> so `?` works.
        let primers = primers.map(|p| p.into_iter().map(to_gb).collect::<PyResult<Vec<_>>>()).transpose()?;
        Ok(self.inner.write_genbank(&path, primers)?)
    }

    fn __repr__(&self) -> String {
        format!(
            "<Document '{}' {} bp {} {} features>",
            self.name(),
            self.length(),
            if self.circular() { "circular" } else { "linear" },
            self.inner.record.features.len()
        )
    }
}

/// Python class `genomeedit_core.Fasta(path)`.
///
/// `#[pyclass]` is a PyO3 macro: it generates the code that makes this struct a Python type.
#[pyclass]
struct Fasta {
    inner: fasta::IndexedFasta,
}

// `#[pymethods]` exposes the methods below to Python.
#[pymethods]
impl Fasta {
    /// Python: `Fasta("genome.fa")`. Builds `genome.fa.fai` on first open.
    #[new] // marks the constructor (Python's __init__)
    fn new(py: Python<'_>, path: PathBuf) -> PyResult<Self> {
        // allow_threads releases Python's GIL while we scan the file, so the GUI is not
        // blocked by Rust work (matters when indexing a multi-GB genome).
        // Rust I/O errors convert automatically into Python OSError via the `?`.
        let inner = py.allow_threads(|| fasta::IndexedFasta::open(&path))?;
        Ok(Fasta { inner })
    }

    /// Python: `fa.sequences()` -> list of (name, length) tuples.
    fn sequences(&self) -> Vec<(String, u64)> {
        // `.map(...).collect()` = Python's [(...) for e in entries]
        self.inner.entries.iter().map(|e| (e.name.clone(), e.length)).collect()
    }

    /// Python: `fa.search("chr1", "GAATTC")` -> [(start, strand)], both strands, exact.
    fn search(&mut self, py: Python<'_>, name: &str, motif: &str) -> PyResult<Vec<(u64, i8)>> {
        let inner = &mut self.inner;
        Ok(py.allow_threads(|| inner.search(name, motif.as_bytes()))?)
    }

    /// Python: `fa.search_motif("chr1", "GARTTC", max_mismatches=0, iupac=True, both=True, limit=100000)`
    /// -> [(start, end, strand, mismatches)], like `genomeedit_core.find_motif` on that sequence.
    #[pyo3(signature = (name, motif, max_mismatches=0, iupac=true, both=true, limit=100_000))]
    fn search_motif(
        &mut self,
        py: Python<'_>,
        name: &str,
        motif: &str,
        max_mismatches: usize,
        iupac: bool,
        both: bool,
        limit: usize,
    ) -> PyResult<Vec<(usize, usize, i8, usize)>> {
        let inner = &mut self.inner;
        let hits = py.allow_threads(|| inner.search_motif(name, motif.as_bytes(), max_mismatches, iupac, both, limit))?;
        Ok(hits.into_iter().map(|h| (h.start, h.end, if h.forward { 1 } else { -1 }, h.mismatches)).collect())
    }

    /// Python: `fa.fetch("chr1", start, end)` -> str of bases, 0-based half-open.
    fn fetch(&mut self, name: &str, start: u64, end: u64) -> PyResult<String> {
        let bytes = self.inner.fetch(name, start, end)?;
        Ok(String::from_utf8_lossy(&bytes).into_owned())
    }
}

/// Python class `genomeedit_core.Track(path)`: a GFF3/GTF/BED annotation file, indexed.
#[pyclass]
struct Track {
    inner: tracks::Track,
}

type TrackItem = (u64, u64, i8, String, String); // start, end, strand, name, kind

fn item_tuple(i: &tracks::Item) -> TrackItem {
    (i.start, i.end, i.strand, i.name.clone(), i.kind.clone())
}

#[pymethods]
impl Track {
    #[new]
    fn new(py: Python<'_>, path: PathBuf) -> PyResult<Self> {
        Ok(Track { inner: py.allow_threads(|| tracks::Track::open(&path))? })
    }

    /// [(sequence name, number of items)]
    fn sequences(&self) -> Vec<(String, usize)> {
        self.inner.sequences()
    }

    /// (items overlapping [start, end) as (start, end, strand, name, kind), more_than_limit)
    #[pyo3(signature = (seq, start, end, limit=5000))]
    fn query(&self, seq: &str, start: u64, end: u64, limit: usize) -> (Vec<TrackItem>, bool) {
        let (items, more) = self.inner.query(seq, start, end, limit);
        (items.into_iter().map(item_tuple).collect(), more)
    }

    /// Items per bin over [start, end) (for a density plot).
    fn density(&self, py: Python<'_>, seq: &str, start: u64, end: u64, bins: usize) -> Vec<u32> {
        let inner = &self.inner;
        py.allow_threads(|| inner.density(seq, start, end, bins))
    }

    /// [(sequence, start, end, strand, name, kind)] whose name contains `needle`.
    #[pyo3(signature = (needle, limit=200))]
    fn search(&self, needle: &str, limit: usize) -> Vec<(String, u64, u64, i8, String, String)> {
        self.inner
            .search(needle, limit)
            .into_iter()
            .map(|(s, i)| (s.to_string(), i.start, i.end, i.strand, i.name.clone(), i.kind.clone()))
            .collect()
    }
}

/// Python class `genomeedit_core.Bam(path)`: a sorted BAM file with its .bai/.csi index
/// (written next to it if missing).
#[pyclass]
struct Bam {
    inner: alignments::Bam,
}

/// (start, end, strand, mapq, name, mismatches [(pos, base)], insertions [(pos, len)],
/// deletions [(start, end)], flags, mate start or None, mate on another sequence, insert)
type ReadTuple = (u64, u64, i8, u8, String, Vec<(u64, char)>, Vec<(u64, u32)>, Vec<(u64, u64)>, u16, Option<u64>, bool, i64);

#[pymethods]
impl Bam {
    #[new]
    fn new(py: Python<'_>, path: PathBuf) -> PyResult<Self> {
        Ok(Bam { inner: py.allow_threads(|| alignments::Bam::open(&path))? })
    }

    /// [(sequence name, length)] from the BAM header.
    fn sequences(&self) -> Vec<(String, u64)> {
        self.inner.sequences()
    }

    /// (reads overlapping [start, end) as ReadTuple, more_than_limit). `reference` = the
    /// reference bases of [start, end), for mismatches.
    #[pyo3(signature = (seq, start, end, reference, limit=20000))]
    fn reads(&self, py: Python<'_>, seq: &str, start: u64, end: u64, reference: &str, limit: usize) -> PyResult<(Vec<ReadTuple>, bool)> {
        let inner = &self.inner;
        let (reads, more) = py.allow_threads(|| inner.reads(seq, start, end, reference.as_bytes(), limit))?;
        let tuples = reads
            .into_iter()
            .map(|r| {
                let mism = r.mismatches.into_iter().map(|(p, b)| (p, b as char)).collect();
                (r.start, r.end, if r.forward { 1 } else { -1 }, r.mapq, r.name, mism, r.insertions, r.deletions,
                 r.flags, r.mate_start, r.mate_elsewhere, r.insert)
            })
            .collect();
        Ok((tuples, more))
    }

    /// Base counts per position of [start, end): (A, C, G, T, N, deleted, insertion before).
    #[pyo3(signature = (seq, start, end, min_baseq=0, min_mapq=0))]
    fn pileup(&self, py: Python<'_>, seq: &str, start: u64, end: u64, min_baseq: u8, min_mapq: u8) -> PyResult<Vec<alignments::Counts>> {
        let inner = &self.inner;
        Ok(py.allow_threads(|| inner.pileup(seq, start, end, min_baseq, min_mapq))?)
    }

    /// [(pos, kind "snp"/"del"/"ins", reference bases, alt bases, reads, depth)] where at
    /// least `min_fraction` of `min_depth`+ reads differ from `reference` (the whole sequence).
    #[pyo3(signature = (seq, reference, min_depth=5, min_fraction=0.2, min_baseq=0, min_mapq=0))]
    fn variants(
        &self,
        py: Python<'_>,
        seq: &str,
        reference: &str,
        min_depth: u32,
        min_fraction: f32,
        min_baseq: u8,
        min_mapq: u8,
    ) -> PyResult<Vec<(u64, String, String, String, u32, u32)>> {
        let inner = &self.inner;
        let found = py.allow_threads(|| inner.variants(seq, reference.as_bytes(), min_depth, min_fraction, min_baseq, min_mapq))?;
        Ok(found.into_iter().map(|v| (v.pos, v.kind.to_string(), v.reference, v.alt, v.count, v.depth)).collect())
    }

    /// Mean read depth in each of `bins` bins of [start, end).
    fn coverage(&self, py: Python<'_>, seq: &str, start: u64, end: u64, bins: usize) -> PyResult<Vec<f32>> {
        let inner = &self.inner;
        Ok(py.allow_threads(|| inner.coverage(seq, start, end, bins))?)
    }

    /// Write the whole file as text SAM at `out` (releases Python's lock while it runs).
    fn write_sam(&self, py: Python<'_>, out: PathBuf) -> PyResult<()> {
        let inner = &self.inner;
        Ok(py.allow_threads(|| inner.write_sam(&out))?)
    }

    /// True when wide-window coverage of `seq` is ready (else `coverage` of a wide window
    /// first streams the whole sequence's reads: call it from a background thread).
    fn has_depth(&self, seq: &str) -> bool {
        self.inner.has_depth(seq)
    }
}

/// A feature as Python passes it in: (kind, GenBank location text, qualifiers).
type FeatureTuple = (String, String, Vec<(String, Option<String>)>);
// pos, removed, removed count, inserted, inserted count, added features, deleted features, note
type HistoryStep = (usize, String, usize, String, usize, Vec<String>, Vec<String>, String);

/// Python tuple -> gb-io feature; a location that does not parse becomes ValueError.
fn to_gb((kind, location, qualifiers): FeatureTuple) -> PyResult<gb_io::seq::Feature> {
    let parsed = gb_io::seq::Location::from_gb_format(location.trim()).map_err(|e| {
        pyo3::exceptions::PyValueError::new_err(format!("bad location '{location}': {e}"))
    })?;
    Ok(gb_io::seq::Feature {
        kind: kind.into(), // String -> Cow<'static, str> (gb-io's text type)
        location: parsed,
        qualifiers: qualifiers.into_iter().map(|(k, v)| (k.into(), v)).collect(),
    })
}

/// `genomeedit_core.tm("ACGT...")`: melting temperature in °C (None if not plain ACGT).
/// One quick-alignment hit: (start, end, strand, score, identity, query start, query end,
/// mismatches (pos, base), insertions (pos, length), deletions (start, end)). Positions are
/// on the template, 0-based; they may pass its length for a hit across the origin.
type AlignTuple = (usize, usize, i8, i32, f32, usize, usize, Vec<(usize, char)>, Vec<(usize, usize)>, Vec<(usize, usize)>);

/// `align(template, query, circular, max_hits=10)`: the best non-overlapping local
/// alignments of `query` on `template`, either strand (strand -1: the query's reverse
/// complement matches the top strand), best first. See align.rs.
#[pyfunction]
#[pyo3(name = "align")]
#[pyo3(signature = (template, query, circular=false, max_hits=10))]
fn quick_align(py: Python<'_>, template: &str, query: &str, circular: bool, max_hits: usize) -> Vec<AlignTuple> {
    // allow_threads releases Python's lock: scanning a genome must not freeze the GUI.
    let hits = py.allow_threads(|| align::align(template.as_bytes(), query.as_bytes(), circular, max_hits));
    hits.into_iter()
        .map(|h| {
            let mism = h.mismatches.into_iter().map(|(p, b)| (p, b as char)).collect();
            (h.start, h.end, if h.forward { 1 } else { -1 }, h.score, h.identity, h.qstart, h.qend, mism, h.insertions, h.deletions)
        })
        .collect()
}

/// `#[pyfunction]` = a module-level Python function rather than a method.
#[pyfunction]
#[pyo3(signature = (seq, na_mm=50.0, mg_mm=1.5, dntp_mm=0.6, primer_nm=50.0))]
fn tm(seq: &str, na_mm: f64, mg_mm: f64, dntp_mm: f64, primer_nm: f64) -> Option<f64> {
    primers::tm(seq.as_bytes(), &primers::Conditions { na_mm, mg_mm, dntp_mm, primer_nm })
}

/// `genomeedit_core.find_primers(template, circular, primers, max_mismatches, seed)`:
/// binding sites of every primer on a sequence string, as one list of
/// (start, end, strand, mismatches) per primer. Takes a copy of the sequence (not a
/// Document), so it can run in a background thread while the document is being edited.
#[pyfunction]
#[pyo3(signature = (template, circular, primers, max_mismatches=1, seed=primers::SEED))]
fn find_primers(
    py: Python<'_>,
    template: &str,
    circular: bool,
    primers: Vec<String>,
    max_mismatches: usize,
    seed: usize,
) -> Vec<Vec<(usize, usize, i8, usize)>> {
    // Uppercase once here instead of on every comparison.
    let upper = template.to_ascii_uppercase();
    py.allow_threads(|| {
        primers
            .iter()
            .map(|p| {
                primers::binding_sites(upper.as_bytes(), p.as_bytes(), circular, seed, max_mismatches)
                    .into_iter()
                    .map(|s| (s.start, s.start + s.len, if s.forward { 1 } else { -1 }, s.mismatches))
                    .collect()
            })
            .collect()
    })
}

/// `genomeedit_core.find_motif(seq, motif, max_mismatches=0, iupac=True, both=True, circular=False, limit=100000)`: where a DNA
/// motif occurs (IUPAC letters allowed, up to `max_mismatches` positions differing, on the top strand and with `both` as its
/// reverse complement), as (start, end, strand, mismatches) in order of position; see search.rs. Takes a copy of the sequence.
#[pyfunction]
#[pyo3(signature = (seq, motif, max_mismatches=0, iupac=true, both=true, circular=false, limit=100_000))]
fn find_motif(
    py: Python<'_>,
    seq: &str,
    motif: &str,
    max_mismatches: usize,
    iupac: bool,
    both: bool,
    circular: bool,
    limit: usize,
) -> Vec<(usize, usize, i8, usize)> {
    py.allow_threads(|| {
        search::find_motif(seq.as_bytes(), motif.as_bytes(), max_mismatches, iupac, both, circular, limit)
            .into_iter()
            .map(|h| (h.start, h.end, if h.forward { 1 } else { -1 }, h.mismatches))
            .collect()
    })
}

/// `genomeedit_core.find_protein(seq, query, max_mismatches=0, limit=100000)`: where an amino-acid motif occurs in the
/// translations of a DNA sequence in the three frames of both strands, as (start, end, strand, mismatches, frame) with the
/// DNA positions on the top strand, in order of position; see search.rs.
#[pyfunction]
#[pyo3(signature = (seq, query, max_mismatches=0, limit=100_000, table=1))]
fn find_protein(py: Python<'_>, seq: &str, query: &str, max_mismatches: usize, limit: usize, table: u8) -> Vec<(usize, usize, i8, usize, usize)> {
    let code = code::by_id_or_standard(table);
    py.allow_threads(|| {
        search::find_protein_with(seq.as_bytes(), query.as_bytes(), max_mismatches, limit, code)
            .into_iter()
            .map(|h| (h.start, h.end, if h.forward { 1 } else { -1 }, h.mismatches, h.frame))
            .collect()
    })
}

/// `genomeedit_core.genetic_codes()`: every embedded NCBI translation table as (id, name, amino acids, starts), the two
/// strings being NCBI's 64-letter ones in TCAG order (see code.rs).
#[pyfunction]
fn genetic_codes() -> Vec<(u8, String, String, String)> {
    code::TABLES.iter().map(|t| (t.id, t.name.to_string(), t.aas.to_string(), t.starts.to_string())).collect()
}

/// `genomeedit_core.translate(seq, table=1)`: amino acids of `seq` read from its first base ('X' for a codon with an N,
/// '*' for a stop). An unknown table id falls back to the standard code.
#[pyfunction]
#[pyo3(signature = (seq, table=1))]
fn translate(seq: &str, table: u8) -> String {
    String::from_utf8(code::by_id_or_standard(table).translate(seq.as_bytes())).unwrap_or_default()
}

/// `genomeedit_core.find_guide_sites(template, circular, guides, pam, max_mismatches, pam_first)`: every site of each
/// guide (protospacer + PAM, or PAM + protospacer for a Cas12a-type nuclease: `pam_first`) on either strand, with at most `max_mismatches` mismatches in the protospacer, as one
/// list per guide of (start, end, strand, mismatches, mask). `mask` bit i = guide base i (0 = 5' end) mismatches.
/// `end` may pass the length on a circular template. Takes a copy of the sequence, like `find_primers`, so it can
/// run in a background thread. All guides must have the same length (at most 32).
#[pyfunction]
#[pyo3(signature = (template, circular, guides, pam, max_mismatches=3, pam_first=false))]
fn find_guide_sites(
    py: Python<'_>,
    template: &str,
    circular: bool,
    guides: Vec<String>,
    pam: &str,
    max_mismatches: usize,
    pam_first: bool,
) -> Vec<Vec<(usize, usize, i8, usize, u32)>> {
    let guides: Vec<Vec<u8>> = guides.into_iter().map(String::into_bytes).collect();
    // allow_threads lets other Python threads run while Rust searches.
    let found = py.allow_threads(|| cas9::find_sites(template.as_bytes(), circular, &guides, pam.as_bytes(), max_mismatches, pam_first));
    found
        .into_iter()
        .map(|sites| {
            sites.into_iter().map(|s| (s.start, s.end, if s.forward { 1 } else { -1 }, s.mismatches, s.mask)).collect()
        })
        .collect()
}

/// `genomeedit_core.find_parts(template, circular, parts, coding, min_dna, min_protein)`:
/// where known parts occur (see annotate.rs), as (part index, start, end, strand,
/// identity 0..1, found as protein?). `end` may pass the origin of a circular template.
#[pyfunction]
#[pyo3(signature = (template, circular, parts, coding, min_dna=0.95, min_protein=0.9, table=1))]
fn find_parts(
    py: Python<'_>,
    template: &str,
    circular: bool,
    parts: Vec<String>,
    coding: Vec<bool>,
    min_dna: f32,
    min_protein: f32,
    table: u8,
) -> Vec<(usize, usize, usize, i8, f32, bool)> {
    let parts: Vec<Vec<u8>> = parts.into_iter().map(String::into_bytes).collect();
    let code = code::by_id_or_standard(table);
    py.allow_threads(|| {
        annotate::find_parts_with(template.as_bytes(), circular, &parts, &coding, min_dna, min_protein, code)
            .into_iter()
            .map(|h| (h.part, h.start, h.start + h.len, if h.forward { 1 } else { -1 }, h.identity, h.protein))
            .collect()
    })
}

/// `genomeedit_core.dimer(a, b)`: (ΔG of the most stable dimer, ΔG of the most stable
/// 3'-end dimer) in kcal/mol at 37 °C; same sequence twice = self-dimer.
#[pyfunction]
#[pyo3(signature = (a, b, na_mm=50.0, mg_mm=1.5, dntp_mm=0.6))]
fn dimer(a: &str, b: &str, na_mm: f64, mg_mm: f64, dntp_mm: f64) -> (f64, f64) {
    let c = primers::Conditions { na_mm, mg_mm, dntp_mm, ..Default::default() };
    primers::dimer_dg(a.as_bytes(), b.as_bytes(), &c)
}

/// `genomeedit_core.hairpin(seq)`: (ΔG of the most stable hairpin, ΔG of the most stable
/// hairpin involving the 3' end), kcal/mol at 37 °C.
#[pyfunction]
#[pyo3(signature = (seq, na_mm=50.0, mg_mm=1.5, dntp_mm=0.6))]
fn hairpin(seq: &str, na_mm: f64, mg_mm: f64, dntp_mm: f64) -> (f64, f64) {
    let c = primers::Conditions { na_mm, mg_mm, dntp_mm, ..Default::default() };
    primers::hairpin_dg(seq.as_bytes(), &c)
}

/// The module itself: `import genomeedit_core`. The function name must match the module name.
#[pymodule]
fn genomeedit_core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(tm, m)?)?;
    m.add_function(wrap_pyfunction!(quick_align, m)?)?;
    m.add_function(wrap_pyfunction!(find_primers, m)?)?;
    m.add_function(wrap_pyfunction!(find_parts, m)?)?;
    m.add_function(wrap_pyfunction!(find_guide_sites, m)?)?;
    m.add_function(wrap_pyfunction!(find_motif, m)?)?;
    m.add_function(wrap_pyfunction!(find_protein, m)?)?;
    m.add_function(wrap_pyfunction!(genetic_codes, m)?)?;
    m.add_function(wrap_pyfunction!(translate, m)?)?;
    m.add_function(wrap_pyfunction!(dimer, m)?)?;
    m.add_function(wrap_pyfunction!(hairpin, m)?)?;
    m.add_class::<Fasta>()?;
    m.add_class::<Document>()?;
    m.add_class::<Feature>()?;
    m.add_class::<Track>()?;
    m.add_class::<Bam>()?;
    Ok(())
}
