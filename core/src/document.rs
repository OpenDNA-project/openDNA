//! An open, editable sequence document: one record (sequence + features), its primers,
//! and the undo/redo history.
//!
//! Every edit goes through `replace(start, end, bases)`:
//!   insert = replace(pos, pos, "ACGT")      delete = replace(start, end, "")
//! so there is exactly one code path that changes the sequence and moves the features.

use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;
use std::sync::OnceLock;

use gb_io::seq::{Feature, Location, Seq};

use crate::locations;
use crate::dna_file::{self, Primer};

/// GenBank feature type for a primer binding site.
const PRIMER_KIND: &str = "primer_bind";

/// `note` of a step that is the app's own and changes nothing in the document (see `Document::external`).
const EXTERNAL: &str = "external";

/// One undoable edit. To undo, we put back `removed` and the feature list as it was.
///
/// shortcut: stores a full copy of the feature list per edit (~ms for 10k features);
/// switch to per-feature diffs if typing in very large annotated genomes feels slow.
struct Edit {
    pos: usize,
    removed: Vec<u8>,
    inserted: Vec<u8>,
    features_before: Vec<Feature>,
    /// `Some(list)`: the edit sets the feature list directly (feature edits, reverse
    /// complement). `None`: features just follow the bases (typing, deleting).
    features_after: Option<Vec<Feature>>,
    /// What kind of whole-sequence edit this is ("reverse complement", ...), "" otherwise.
    note: &'static str,
    /// For the history panel: filled when the edit is made (see `summarize`).
    summary: Option<Summary>,
    /// `Some((before, after))`: the edit also switches circular <-> linear (undo switches
    /// back). Making a plasmid linear is a rotation plus this, so it is one undo step.
    topology: Option<(gb_io::seq::Topology, gb_io::seq::Topology)>,
    /// How this edit moves the "edited bases" highlight (see `Remap`).
    remap: Remap,
    /// The edited stretches just before and after this edit (see `EditedMemo`).
    memo: EditedMemo,
}

/// A snapshot of the edited stretches around one edit, so that undo and redo restore them
/// EXACTLY (for instance an edited base that was deleted comes back still edited). It is only
/// trusted while no save happened since (`epoch` equal to the document's): a save empties
/// the stretches, and then undo / redo fall back to moving them by `remap_edited`.
#[derive(Default)]
struct EditedMemo {
    epoch: u64,
    before: Vec<(usize, usize)>,
    after: Vec<(usize, usize)>,
}

/// How an edit moves the list of edited stretches (`Document::edited`).
///
/// Most edits are "local": bases `[pos, pos + removed.len())` are swapped for `inserted`, and
/// everything after shifts. A few rewrite the whole sequence in a way that carries the bases
/// with it, so they need their own rule, or every base would look newly edited.
/// (`enum` = a type that is exactly one of several named variants, some carrying data;
/// `#[derive(Clone, Copy)]` lets it be copied around freely like a number.)
#[derive(Clone, Copy)]
enum Remap {
    /// Use the edit's own `pos` / `removed` / `inserted`.
    Local,
    /// Reverse complement: the stretches flip end to end.
    Mirror,
    /// Rotation so that old base `k` becomes base 0: the stretches slide round the circle.
    Rotate(usize),
    /// An edit across the origin of a circle, made of two local edits (see
    /// `edit_across_origin`): drop the first `head` bases, then replace `removed` bases at
    /// `pos - head` by `inserted` new ones.
    Across { head: usize, pos: usize, removed: usize, inserted: usize },
}

/// Sort stretches and join those that overlap or touch, dropping empty ones.
/// `&mut Vec<..>` = a borrowed list we may change in place (the caller keeps ownership).
fn normalize(ranges: &mut Vec<(usize, usize)>) {
    ranges.retain(|&(a, b)| a < b);
    ranges.sort_unstable();
    let mut out: Vec<(usize, usize)> = Vec::with_capacity(ranges.len());
    for &(a, b) in ranges.iter() {
        match out.last_mut() {
            Some(last) if a <= last.1 => last.1 = last.1.max(b), // touches the previous one: extend it
            _ => out.push((a, b)),
        }
    }
    *ranges = out;
}

/// Local edit going forward: `removed` bases at `pos` became `inserted` new ones. Stretches
/// lose the part inside the removed span, later ones shift, and the new bases are marked.
fn edited_forward(ranges: &mut Vec<(usize, usize)>, pos: usize, removed: usize, inserted: usize) {
    let cut_end = pos + removed;
    let mut out = Vec::with_capacity(ranges.len() + 1);
    for &(a, b) in ranges.iter() {
        if a < pos {
            out.push((a, b.min(pos))); // the part before the edit stays
        }
        if b > cut_end {
            // The part after the removed span moves by the change in length. Both ends are
            // at least `cut_end >= removed`, so subtracting `removed` last cannot underflow.
            out.push((a.max(cut_end) + inserted - removed, b + inserted - removed));
        }
    }
    out.push((pos, pos + inserted));
    normalize(&mut out);
    *ranges = out;
}

/// Local edit going backward (undo): the `inserted` bases at `pos` are gone again and
/// `removed` old ones are back. Stretches lose what was inside the inserted span and later
/// ones shift back. The bases that come back are not marked (they were not typed).
fn edited_backward(ranges: &mut Vec<(usize, usize)>, pos: usize, removed: usize, inserted: usize) {
    let ins_end = pos + inserted;
    let mut out = Vec::with_capacity(ranges.len());
    for &(a, b) in ranges.iter() {
        if a < pos {
            out.push((a, b.min(pos)));
        }
        if b > ins_end {
            out.push((a.max(ins_end) + removed - inserted, b + removed - inserted));
        }
    }
    normalize(&mut out);
    *ranges = out;
}

/// Slide stretches round a circle of length `n` so old base `k` becomes base 0.
fn edited_rotate(ranges: &mut Vec<(usize, usize)>, k: usize, n: usize) {
    let mut out = Vec::with_capacity(ranges.len() + 1);
    for &(a, b) in ranges.iter() {
        if a < k {
            out.push((a + n - k, b.min(k) + n - k)); // before k: they end up at the far end
        }
        if b > k {
            out.push((a.max(k) - k, b - k));
        }
    }
    normalize(&mut out);
    *ranges = out;
}

/// Move the stretches the way the edit's `remap` says (`forward`: the edit is being applied;
/// else it is being undone). `n` is the sequence length (the same before and after for
/// mirror and rotation, the only cases that use it).
fn remap_edited(ranges: &mut Vec<(usize, usize)>, edit: &Edit, forward: bool, n: usize) {
    match edit.remap {
        Remap::Local => {
            let (pos, rem, ins) = (edit.pos, edit.removed.len(), edit.inserted.len());
            if forward { edited_forward(ranges, pos, rem, ins) } else { edited_backward(ranges, pos, rem, ins) }
        }
        Remap::Mirror => {
            for r in ranges.iter_mut() {
                *r = (n - r.1, n - r.0);
            }
            normalize(ranges);
        }
        // Undoing a rotation by k is a rotation by n - k.
        Remap::Rotate(k) => edited_rotate(ranges, if forward { k } else { n - k }, n),
        Remap::Across { head, pos, removed, inserted } => {
            if forward {
                edited_forward(ranges, 0, head, 0);
                edited_forward(ranges, pos - head, removed, inserted);
            } else {
                edited_backward(ranges, pos - head, removed, inserted);
                edited_backward(ranges, 0, head, 0);
            }
        }
    }
}

/// One history step, short enough to list: where bases changed, the first 30 removed and
/// inserted bases, and the features the step added or deleted (a changed one is in both).
#[derive(Clone, Debug)]
pub struct Summary {
    pub pos: usize,
    pub removed: String,
    pub removed_len: usize,
    pub inserted: String,
    pub inserted_len: usize,
    pub added: Vec<String>,
    pub deleted: Vec<String>,
    pub note: &'static str,
}

/// History summary of an edit (see `Summary`).
fn summarize(edit: &Edit) -> Summary {
    // The whole replaced stretch (what was selected), not just the bases that differ:
    // "GAC > GAG" reads better than "C > G".
    let (r, i) = (&edit.removed, &edit.inserted);
    let short = |b: &[u8]| String::from_utf8_lossy(&b[..b.len().min(30)]).into_owned();
    let (added, deleted) = match &edit.features_after {
        Some(after) => feature_diff(&edit.features_before, after),
        None => (Vec::new(), Vec::new()),
    };
    Summary {
        pos: edit.pos,
        removed: short(r),
        removed_len: r.len(),
        inserted: short(i),
        inserted_len: i.len(),
        added,
        deleted,
        note: edit.note,
    }
}

/// Name shown for a feature: its label, gene, locus tag or product, else its type.
fn feature_name(f: &Feature) -> String {
    ["label", "gene", "locus_tag", "product"]
        .iter()
        .find_map(|k| f.qualifier_values(k).next())
        .unwrap_or(&f.kind)
        .to_string()
}

/// Features in `after` but not `before` (added), and in `before` but not `after` (deleted),
/// compared as whole features (type, location, qualifiers). O(n) with a hash map, so a
/// genome's 10k features cost a few ms per feature edit.
fn feature_diff(before: &[Feature], after: &[Feature]) -> (Vec<String>, Vec<String>) {
    let key = |f: &Feature| format!("{}\t{}\t{:?}", f.kind, f.location.to_gb_format(), f.qualifiers);
    let mut count: std::collections::HashMap<String, i64> = std::collections::HashMap::new();
    for f in before {
        *count.entry(key(f)).or_default() += 1;
    }
    let mut added = Vec::new();
    for f in after {
        let c = count.entry(key(f)).or_default();
        if *c > 0 {
            *c -= 1; // was there before: unchanged
        } else {
            added.push(feature_name(f));
        }
    }
    // What is left in `count` was deleted; walk `before` to keep its order.
    let deleted = before
        .iter()
        .filter(|f| count.get_mut(&key(f)).is_some_and(|c| if *c > 0 { *c -= 1; true } else { false }))
        .map(feature_name)
        .collect();
    (added, deleted)
}

/// What `Document::refresh_translations` did with one CDS.
#[derive(Debug, PartialEq)]
pub enum Refresh {
    Updated,
    Unchanged,
    Skipped(crate::code::Skip, String),
}

pub struct Document {
    pub record: Seq,
    pub primers: Vec<Primer>,
    undo_stack: Vec<Edit>,
    redo_stack: Vec<Edit>,
    /// Undo depth at the last save; `None` once that state can no longer be reached.
    saved_depth: Option<usize>,
    /// Stretches `[start, end)` of bases typed, pasted or replaced since the last save (or
    /// since load): the app paints them yellow. Sorted, no overlaps. Kept in step with every
    /// edit, undo and redo (`remap_edited`); a save empties it. Deletions leave nothing here.
    edited: Vec<(usize, usize)>,
    /// Counts the times `edited` was emptied (saves); see `EditedMemo`.
    edited_epoch: u64,
    /// G+C count per GC_BLOCK bases, computed on first use and dropped on every edit
    /// (zoomed-out GC plots of a 250 Mb chromosome read this, not the bases).
    /// `OnceLock` = a value filled at most once, safely, even through a shared `&self`.
    gc_blocks: OnceLock<Vec<u32>>,
}

const GC_BLOCK: usize = 1024;

/// 1 for G, C, S (strong) in either case, else 0: a table lookup per base is several
/// times faster than a `match` when counting a whole chromosome.
static GC_TABLE: [u32; 256] = {
    let mut t = [0u32; 256];
    t[b'G' as usize] = 1;
    t[b'C' as usize] = 1;
    t[b'S' as usize] = 1;
    t[b'g' as usize] = 1;
    t[b'c' as usize] = 1;
    t[b's' as usize] = 1;
    t
};

fn count_gc(bases: &[u8]) -> u32 {
    bases.iter().map(|&b| GC_TABLE[b as usize]).sum()
}

impl Document {
    pub fn new(mut record: Seq, primers: Vec<Primer>) -> Self {
        unwrap_text(&mut record);
        record.len = Some(record.seq.len()); // gb-io's revcomp() needs the length set
        Document {
            record,
            primers,
            undo_stack: Vec::new(),
            redo_stack: Vec::new(),
            saved_depth: Some(0),
            edited: Vec::new(),
            edited_epoch: 0,
            gc_blocks: OnceLock::new(),
        }
    }

    /// A new, unsaved document (File → New) holding `bases` (may be empty).
    pub fn new_sequence(name: &str, bases: &[u8], circular: bool) -> io::Result<Self> {
        let mut doc = Document::new(Seq { name: Some(name.to_string()), ..Seq::empty() }, Vec::new());
        doc.check_edit(0, 0, bases)?; // only sequence letters
        doc.record.seq = bases.to_vec();
        doc.record.len = Some(bases.len());
        doc.record.molecule_type = Some("DNA".into());
        doc.record.topology = if circular { gb_io::seq::Topology::Circular } else { gb_io::seq::Topology::Linear };
        doc.set_modified(true); // exists only in memory until saved
        Ok(doc)
    }

    /// Open a GenBank (.gb/.gbk/.gbff/.genbank) or .dna file.
    /// The format is detected from the first bytes, not the extension.
    pub fn open(path: &Path) -> io::Result<Self> {
        let head = {
            let mut buf = [0u8; 8];
            let n = io::Read::read(&mut File::open(path)?, &mut buf)?;
            buf[..n].to_vec()
        };
        if head.first() == Some(&0x09) {
            let (record, primers) = dna_file::read_file(path)?;
            return Ok(Document::new(record, primers));
        }
        // gb-io returns one `Seq` per record in the file.
        let mut records = gb_io::reader::parse_file(path)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, format!("GenBank: {e}")))?;
        match records.len() {
            1 => Ok(Document::new(records.pop().unwrap(), Vec::new())),
            0 => Err(io::Error::new(io::ErrorKind::InvalidData, "no GenBank record found")),
            // shortcut: refuse instead of silently dropping records on save; add multi-record
            // documents when a real file needs it.
            n => Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("file has {n} GenBank records; a document holds one record (split the file, or use it as the reference genome in Primers > Options, where every record is searched)"),
            )),
        }
    }

    /// Write the document as GenBank. Writes to a temporary file first and then renames
    /// it over the target, so a crash mid-save can never leave a half-written file.
    ///
    /// `primers`: when given, the file gets these `primer_bind` features instead of the
    /// ones in the document (primer sites are recomputed by the app, not stored). The
    /// document itself is not changed, so this adds no undo step.
    pub fn save_genbank(&mut self, path: &Path, primers: Option<Vec<Feature>>) -> io::Result<()> {
        self.write_genbank(path, primers)?;
        self.saved_depth = Some(self.real_depth());
        self.clear_edited(); // everything is on disk now
        Ok(())
    }

    /// Write into FASTA `dst` a copy of `src` where this document's sequence replaces the
    /// record of the same name (other sequences are streamed through untouched). FASTA
    /// has no features: the caller makes sure there are none to lose. Only reads the
    /// document (`&self`), so the GUI can keep drawing it during a long save; the caller
    /// marks it saved afterwards (`set_modified(false)`).
    pub fn write_fasta(&self, src: &Path, dst: &Path) -> io::Result<()> {
        let name = self.record.name.clone().unwrap_or_default();
        crate::fasta::write_replacing(src, dst, &name, &self.record.seq)
    }

    /// Write a GenBank copy without changing the "saved" state (used for autosave).
    pub fn write_genbank(&self, path: &Path, primers: Option<Vec<Feature>>) -> io::Result<()> {
        // `replaced` must live until the write is done, so it is declared out here.
        let replaced;
        let record = match primers {
            None => &self.record,
            Some(primers) => {
                let mut copy = self.record.clone();
                copy.features.retain(|f| f.kind != PRIMER_KIND);
                for mut f in primers {
                    check_feature(&mut f, self.len() as i64)?;
                    copy.features.push(f);
                }
                replaced = copy;
                &replaced
            }
        };
        let mut tmp = path.as_os_str().to_owned();
        tmp.push(".tmp");
        {
            // Inner block: the file is closed (flushed) when `out` goes out of scope.
            let mut out = BufWriter::new(File::create(&tmp)?);
            gb_io::writer::write(&mut out, record)?;
            out.flush()?;
        }
        let renamed = std::fs::rename(&tmp, path);
        if renamed.is_err() {
            let _ = std::fs::remove_file(&tmp); // e.g. the target is a folder: leave no stray .tmp file behind
        }
        renamed
    }

    /// Open one sequence of a FASTA file (through its .fai index) as a linear document.
    pub fn from_fasta(path: &Path, name: &str) -> io::Result<Self> {
        let mut fasta = crate::fasta::IndexedFasta::open(path)?;
        let length = fasta.entries.iter().find(|e| e.name == name).map(|e| e.length).unwrap_or(0);
        let seq = fasta.fetch(name, 0, length)?;
        let record = Seq {
            name: Some(name.to_string()),
            topology: gb_io::seq::Topology::Linear,
            molecule_type: Some("DNA".into()),
            len: Some(seq.len()),
            seq,
            ..Seq::empty()
        };
        Ok(Document::new(record, Vec::new()))
    }

    /// GC fraction in `bins` equal bins of [start, end) (NaN for an empty bin). Zoomed
    /// out (bins wider than 4 blocks) it sums cached per-block counts instead of reading
    /// every base, so a whole-chromosome plot stays instant.
    pub fn gc_bins(&self, start: usize, end: usize, bins: usize) -> Vec<f32> {
        let end = end.min(self.len());
        if bins == 0 || end <= start {
            return vec![f32::NAN; bins];
        }
        let seq = &self.record.seq;
        let width = (end - start) as f64 / bins as f64;
        let use_blocks = width >= 4.0 * GC_BLOCK as f64;
        let blocks = if use_blocks {
            Some(self.gc_blocks.get_or_init(|| seq.chunks(GC_BLOCK).map(count_gc).collect()))
        } else {
            None
        };
        (0..bins)
            .map(|i| {
                let a = start + (i as f64 * width) as usize;
                let b = (start + ((i + 1) as f64 * width) as usize).min(end);
                if b <= a {
                    return f32::NAN;
                }
                let gc = match blocks {
                    // Whole blocks only (bin edges rounded to blocks: invisible at this zoom).
                    Some(bl) => {
                        let (ba, bb) = (a / GC_BLOCK, (b / GC_BLOCK).max(a / GC_BLOCK + 1).min(bl.len()));
                        let total = ((bb - ba) * GC_BLOCK).min(seq.len() - ba * GC_BLOCK);
                        return bl[ba..bb].iter().sum::<u32>() as f32 / total as f32;
                    }
                    None => count_gc(&seq[a..b]) as usize,
                };
                gc as f32 / (b - a) as f32
            })
            .collect()
    }

    /// The edited stretches `[start, end)` (typed, pasted, replaced or inserted bases since the
    /// last save or load), for the yellow highlight. Sorted, never overlapping.
    pub fn edited_ranges(&self) -> &[(usize, usize)] {
        &self.edited
    }

    /// After an undo or redo that lands exactly on the saved state, nothing differs from the
    /// file, so nothing is highlighted (whatever the snapshots say about steps made before a save).
    fn edited_if_clean(&mut self) {
        if !self.is_modified() {
            self.edited.clear();
        }
    }

    /// Forget all highlights (the document was saved).
    fn clear_edited(&mut self) {
        self.edited.clear();
        self.edited_epoch += 1; // older snapshots no longer describe the document
    }

    pub fn len(&self) -> usize {
        self.record.seq.len()
    }

    pub fn is_modified(&self) -> bool {
        self.saved_depth != Some(self.real_depth())
    }

    /// Undo steps that changed the document. The app's own steps (`external`) are not
    /// counted: they never make the document modified, and undoing them never makes it clean.
    fn real_depth(&self) -> usize {
        self.undo_stack.iter().filter(|e| e.note != EXTERNAL).count()
    }

    /// Force the modified flag: `true` after recovering an autosave (nothing on disk
    /// matches it), `false` to treat the current state as saved.
    pub fn set_modified(&mut self, modified: bool) {
        self.saved_depth = if modified { None } else { Some(self.real_depth()) };
        if !modified {
            self.clear_edited(); // the app marked it saved (FASTA and project saves)
        }
    }

    /// An undo step of the app's own that changes nothing in the document (a change to the
    /// primer database, say): it sits in the history so that Undo / Redo reach it in order.
    /// `label` is the text for the History panel and `token` is what the app uses to find
    /// what to undo and redo; the app does that itself when `next_external` names the step.
    pub fn external(&mut self, label: &str, token: &str) {
        let summary = Summary {
            pos: 0,
            removed: token.to_string(),
            removed_len: 0,
            inserted: label.to_string(),
            inserted_len: 0,
            added: Vec::new(),
            deleted: Vec::new(),
            note: EXTERNAL,
        };
        self.push(Edit {
            pos: 0,
            removed: Vec::new(),
            inserted: Vec::new(),
            features_before: self.record.features.clone(),
            features_after: None,
            note: EXTERNAL,
            summary: Some(summary),
            topology: None,
            remap: Remap::Local,
            memo: EditedMemo::default(),
        });
    }

    /// The token of the next step to undo (`redo` false) or redo (`redo` true) when it is an
    /// `external` one, else `None`. Cheap: the app asks before every step.
    pub fn next_external(&self, redo: bool) -> Option<String> {
        let edit = if redo { self.redo_stack.last() } else { self.undo_stack.last() }?;
        (edit.note == EXTERNAL).then(|| edit.summary.as_ref().map(|s| s.removed.clone()).unwrap_or_default())
    }

    /// Circular (plasmid) or linear. Not an undo step (the history only records bases and
    /// features), so it marks the document modified directly.
    pub fn set_circular(&mut self, circular: bool) {
        let topology = if circular { gb_io::seq::Topology::Circular } else { gb_io::seq::Topology::Linear };
        if self.record.topology != topology {
            self.record.topology = topology;
            self.set_modified(true);
        }
    }

    /// Replace bases `[start, end)` with `bases`, moving features accordingly.
    /// This is the trust boundary for edits coming from the GUI, so it validates everything.
    pub fn replace(&mut self, start: usize, end: usize, bases: &[u8]) -> io::Result<()> {
        self.check_edit(start, end, bases)?;
        if start == end && bases.is_empty() {
            return Ok(()); // nothing to do; don't create an empty undo step
        }
        self.push(Edit {
            pos: start,
            removed: self.record.seq[start..end].to_vec(),
            inserted: bases.to_vec(),
            features_before: self.record.features.clone(),
            features_after: None,
            note: "",
            summary: None,
            topology: None,
            remap: Remap::Local,
            memo: EditedMemo::default(),
        });
        Ok(())
    }

    /// Paste: replace `[start, end)` with `bases` and add `features` (whose locations are
    /// relative to the pasted bases, as returned by `features_between`). One undo step.
    pub fn paste(&mut self, start: usize, end: usize, bases: &[u8], features: Vec<Feature>) -> io::Result<()> {
        self.check_edit(start, end, bases)?;
        let mut after = self.record.features.clone();
        move_features(&mut after, start as i64, end as i64, bases.len() as i64);
        for mut f in features {
            check_feature(&mut f, bases.len() as i64)?;
            // Inserting `start` bases at position 0 = shifting the location right by `start`.
            f.location = locations::after_insert(&f.location, 0, start as i64);
            after.push(f);
        }
        self.push(Edit {
            pos: start,
            removed: self.record.seq[start..end].to_vec(),
            inserted: bases.to_vec(),
            features_before: self.record.features.clone(),
            features_after: Some(after),
            note: "",
            summary: None,
            topology: None,
            remap: Remap::Local,
            memo: EditedMemo::default(),
        });
        Ok(())
    }

    /// The features overlapping `[start, end)`, clipped to it, with locations relative
    /// to `start` (what Copy puts on the clipboard). The whole-record `source` is skipped.
    pub fn features_between(&self, start: usize, end: usize) -> Vec<Feature> {
        let (s, e, len) = (start as i64, end as i64, self.len() as i64);
        self.record
            .features
            .iter()
            .filter(|f| f.kind != "source")
            .filter_map(|f| {
                // Clip = delete everything after the range, then everything before it.
                // A feature entirely outside vanishes (`?` returns None).
                let mut loc = f.location.clone();
                if e < len {
                    loc = locations::after_delete(&loc, e, len)?;
                }
                if s > 0 {
                    loc = locations::after_delete(&loc, 0, s)?;
                }
                Some(Feature { location: loc, ..f.clone() })
            })
            .collect()
    }

    /// `features_between` for a range that may pass the origin (`end` > length, circular
    /// records): the part before the origin, then the part after it. A feature crossing
    /// the origin stays ONE feature, its parts still in reading order (see
    /// `locations::clip_window`).
    pub fn features_across(&self, start: usize, end: usize) -> Vec<Feature> {
        let n = self.len();
        if end <= n {
            return self.features_between(start, end);
        }
        self.record
            .features
            .iter()
            .filter(|f| f.kind != "source")
            .filter_map(|f| {
                let loc = locations::clip_window(&f.location, start as i64, (end - start) as i64, n as i64)?;
                Some(Feature { location: loc, ..f.clone() })
            })
            .collect()
    }

    /// What "Copy Reverse Complement" keeps: the features of `[start, end)` (as `features_across`
    /// gives them, located relative to `start`) flipped onto the reverse-complemented copy.
    /// A feature on the top strand becomes one on the bottom strand, join parts come out in the
    /// new reading order, and every feature still covers exactly the same bases (now reversed).
    pub fn features_between_reversed(&self, start: usize, end: usize) -> Vec<Feature> {
        let width = end - start;
        // The flip is done by gb-io's own `Seq::revcomp`, which needs a record of the right
        // LENGTH to mirror positions against. So we hand it a small throw-away record holding
        // just the clipped features and `width` placeholder bases (the real bases are
        // reverse-complemented separately, by the caller). Linear: the copy is a plain stretch.
        let clipped = Seq {
            seq: vec![b'N'; width],
            features: self.features_across(start, end),
            ..Seq::empty()
        };
        clipped.revcomp().features
    }

    /// Replace `[start, end)` where `end` passes the origin of a circular record
    /// (start < length < end <= start + length) with `bases` (and `features` relative to
    /// them, like `paste`), as ONE undo step. Done as two edits (the part after the
    /// origin, then the rest) squashed into one. Returns where `bases` now start.
    pub fn edit_across_origin(&mut self, start: usize, end: usize, bases: &[u8], features: Option<Vec<Feature>>) -> io::Result<usize> {
        let n = self.len();
        if self.record.topology != gb_io::seq::Topology::Circular || start >= n || end <= n || end - start > n {
            return Err(io::Error::new(io::ErrorKind::InvalidInput, format!("bad range across the origin: {start}..{end}")));
        }
        let head = end - n; // bases after the origin
        let (depth, saved, real) = (self.undo_stack.len(), self.saved_depth, self.real_depth());
        let (seq_before, features_before) = (self.record.seq.clone(), self.record.features.clone());
        let edited_before = self.edited.clone();
        // `(|| { ... })()` = a block whose `?` errors land in `result` instead of returning.
        let result = (|| {
            self.replace(0, head, b"")?;
            let s = start - head;
            match features {
                Some(fs) => self.paste(s, n - head, bases, fs),
                None => self.replace(s, n - head, bases),
            }
        })();
        if let Err(e) = result {
            while self.undo_stack.len() > depth {
                self.undo(); // put back whatever part succeeded
            }
            self.redo_stack.clear();
            self.saved_depth = saved;
            return Err(e);
        }
        // Squash: the edits are applied; record them as one whole-sequence step.
        self.undo_stack.truncate(depth);
        self.saved_depth = saved.filter(|&d| d <= real);
        let mut edit = Edit {
            pos: 0,
            removed: seq_before,
            inserted: self.record.seq.clone(),
            features_before,
            features_after: Some(self.record.features.clone()),
            note: "across the origin",
            summary: None,
            topology: None,
            remap: Remap::Across { head, pos: start, removed: n - start, inserted: bases.len() },
            memo: EditedMemo { epoch: self.edited_epoch, before: edited_before, after: self.edited.clone() },
        };
        let mut summary = summarize(&edit);
        // The whole sequence was rewritten; say what was replaced, where (1-based in Python).
        summary.pos = start;
        summary.removed_len = end - start;
        summary.inserted_len = bases.len();
        summary.inserted = String::from_utf8_lossy(&bases[..bases.len().min(30)]).into_owned();
        summary.removed.clear();
        edit.summary = Some(summary);
        self.undo_stack.push(edit);
        Ok(start - head)
    }

    /// Overwrite bases at `pos` with `bases` (same length): substitutions only, so no
    /// feature moves or shrinks (a replace would trim features starting there).
    pub fn overwrite(&mut self, pos: usize, bases: &[u8]) -> io::Result<()> {
        self.check_edit(pos, pos + bases.len(), bases)?;
        if bases.is_empty() {
            return Ok(());
        }
        self.push(Edit {
            pos,
            removed: self.record.seq[pos..pos + bases.len()].to_vec(),
            inserted: bases.to_vec(),
            features_before: self.record.features.clone(),
            features_after: Some(self.record.features.clone()),
            note: "",
            summary: None,
            topology: None,
            remap: Remap::Local,
            memo: EditedMemo::default(),
        });
        Ok(())
    }

    /// A PCR product as a new, unsaved linear document: `prefix` + bases `[start, end)` +
    /// `suffix` (the primers' 5' tails). On a circular template `end` may pass the
    /// sequence end (the product spans the origin). Features in the range come along,
    /// clipped. A feature crossing the origin of the template becomes two pieces.
    pub fn pcr_product(&self, name: &str, start: usize, end: usize, prefix: &[u8], suffix: &[u8]) -> io::Result<Document> {
        let n = self.len();
        let circular = self.record.topology == gb_io::seq::Topology::Circular;
        if start >= n || end <= start || end > if circular { start + n } else { n } {
            return Err(io::Error::new(io::ErrorKind::InvalidInput, format!("bad product range {start}..{end}")));
        }
        self.check_edit(0, 0, prefix)?;
        self.check_edit(0, 0, suffix)?;
        // The template part, in one or two pieces.
        let pieces = if end <= n { vec![(start, end)] } else { vec![(start, n), (0, end - n)] };
        let mut seq = prefix.to_vec();
        let mut features = Vec::new();
        for (a, b) in pieces {
            let offset = seq.len() as i64; // where this piece starts in the product
            for mut f in self.features_between(a, b) {
                f.location = locations::after_insert(&f.location, 0, offset);
                features.push(f);
            }
            seq.extend_from_slice(&self.record.seq[a..b]);
        }
        seq.extend_from_slice(suffix);
        let record = Seq {
            name: Some(name.to_string()),
            topology: gb_io::seq::Topology::Linear,
            molecule_type: Some("DNA".into()),
            len: Some(seq.len()),
            seq,
            features,
            ..Seq::empty()
        };
        let mut doc = Document::new(record, Vec::new());
        doc.set_modified(true); // exists only in memory until saved
        Ok(doc)
    }

    /// Add several features as one undo step (e.g. a primer binding at several sites).
    /// Returns the index of the first one.
    pub fn add_features(&mut self, features: Vec<Feature>) -> io::Result<usize> {
        let first = self.record.features.len();
        let mut after = self.record.features.clone();
        for mut f in features {
            check_feature(&mut f, self.len() as i64)?;
            after.push(f);
        }
        self.push_features(after);
        Ok(first)
    }

    /// Range inside the sequence and only DNA letters.
    fn check_edit(&self, start: usize, end: usize, bases: &[u8]) -> io::Result<()> {
        if start > end || end > self.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("range {start}..{end} outside sequence of {} bp", self.len()),
            ));
        }
        if let Some(b) = bases.iter().find(|b| !is_iupac(**b)) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("'{}' is not a DNA letter (allowed: ACGTU + IUPAC codes, N, -)", *b as char),
            ));
        }
        Ok(())
    }

    /// Replace feature `index`, or append a new one when `index` is `None`.
    /// Returns the feature's index. Validated here: this is called with user input.
    pub fn set_feature(&mut self, index: Option<usize>, mut feature: Feature) -> io::Result<usize> {
        let bad = |msg: String| io::Error::new(io::ErrorKind::InvalidInput, msg);
        check_feature(&mut feature, self.len() as i64)?;
        let mut after = self.record.features.clone();
        let index = match index {
            Some(i) if i < after.len() => {
                after[i] = feature;
                i
            }
            Some(i) => return Err(bad(format!("no feature #{i}"))),
            None => {
                after.push(feature);
                after.len() - 1
            }
        };
        self.push_features(after);
        Ok(index)
    }

    pub fn delete_feature(&mut self, index: usize) -> io::Result<()> {
        if index >= self.record.features.len() {
            return Err(io::Error::new(io::ErrorKind::InvalidInput, format!("no feature #{index}")));
        }
        let mut after = self.record.features.clone();
        after.remove(index);
        self.push_features(after);
        Ok(())
    }

    /// Several feature changes as ONE undo step (e.g. grouping parts under a gene):
    /// `changes` = (Some(index), Some(f)) replace, (None, Some(f)) add at the end,
    /// (Some(index), None) delete. Indexes refer to the features before the call.
    /// Nothing changes if any of them is invalid.
    pub fn edit_features(&mut self, changes: Vec<(Option<usize>, Option<Feature>)>) -> io::Result<()> {
        let n = self.record.features.len();
        let bad = |msg: String| io::Error::new(io::ErrorKind::InvalidInput, msg);
        // `Option<Feature>` per old index: None = deleted. Built first, then compacted,
        // so deleting #2 does not shift what "#5" means for a later change.
        let mut slots: Vec<Option<Feature>> = self.record.features.iter().cloned().map(Some).collect();
        let mut added = Vec::new();
        for (index, feature) in changes {
            let feature = match feature {
                Some(mut f) => {
                    check_feature(&mut f, self.len() as i64)?;
                    Some(f)
                }
                None => None,
            };
            match (index, feature) {
                (Some(i), _) if i >= n => return Err(bad(format!("no feature #{i}"))),
                (Some(i), f) => slots[i] = f,
                (None, Some(f)) => added.push(f),
                (None, None) => {}
            }
        }
        // `flatten` drops the Nones (deleted features).
        let after: Vec<Feature> = slots.into_iter().flatten().chain(added).collect();
        self.push_features(after);
        Ok(())
    }

    /// Recompute the /translation qualifier of CDS features from their own bases (see `code::cds_protein` for the rules
    /// and for which features are left alone). `only`: the feature indexes to look at, or `None` for every CDS.
    /// A CDS that has no /translation never gets one. All the changes are ONE undo step; when nothing differs
    /// nothing is pushed (so running it again, or saving twice, never adds history or marks the file as edited).
    ///
    /// Returns one `(feature index, outcome)` per CDS that had a /translation: `Updated`, `Unchanged`, or
    /// `Skipped(kind, reason)`.
    pub fn refresh_translations(&mut self, only: Option<&[usize]>, default_table: u8) -> Vec<(usize, Refresh)> {
        let mut report = Vec::new();
        let mut features = self.record.features.clone();
        let mut changed = false;
        for (i, f) in features.iter_mut().enumerate() {
            if f.kind != "CDS" || only.is_some_and(|o| !o.contains(&i)) {
                continue;
            }
            // Position of the /translation qualifier (the first one, if somehow there are several).
            let Some(q) = f.qualifiers.iter().position(|(k, _)| k == "translation") else { continue };
            match crate::code::cds_protein(&self.record.seq, f, default_table) {
                Err((kind, why)) => report.push((i, Refresh::Skipped(kind, why))),
                Ok(protein) => {
                    // GenBank wraps long values over lines; compare without any white space.
                    let old: String = f.qualifiers[q].1.as_deref().unwrap_or("").split_whitespace().collect();
                    if old.eq_ignore_ascii_case(&protein) {
                        report.push((i, Refresh::Unchanged));
                    } else {
                        f.qualifiers[q].1 = Some(protein);
                        changed = true;
                        report.push((i, Refresh::Updated));
                    }
                }
            }
        }
        if changed {
            self.push_features(features);
        }
        report
    }

    /// Flip the whole record to the other strand: sequence and every feature location.
    pub fn reverse_complement(&mut self) {
        // gb-io does the work. It merges touching join parts (join(1..5,6) -> 1..6, same
        // bases) and would drop a feature whose location it cannot flip (none seen in
        // practice; undo brings it back).
        let flipped = self.record.revcomp();
        self.push(Edit {
            pos: 0,
            removed: self.record.seq.clone(),
            inserted: flipped.seq,
            features_before: self.record.features.clone(),
            features_after: Some(flipped.features),
            note: "reverse complement",
            summary: None,
            topology: None,
            remap: Remap::Mirror,
            memo: EditedMemo::default(),
        });
    }

    /// Rotate a circular sequence so base `k` (0-based) becomes base 1, features along
    /// (one crossing the new start becomes a join of its two pieces; `source` covers
    /// everything). One undo step.
    pub fn rotate(&mut self, k: usize) -> io::Result<()> {
        let edit = self.rotation(k)?;
        self.push(edit);
        Ok(())
    }

    /// Make a circular sequence linear, broken just before base `k` (0-based; 0 = keep the
    /// current start): the rotation of `rotate` plus the topology change, as ONE undo step
    /// (so undo gives back the circle with its old origin).
    pub fn linearize(&mut self, k: usize) -> io::Result<()> {
        if self.record.topology != gb_io::seq::Topology::Circular {
            return Err(io::Error::new(io::ErrorKind::InvalidInput, "the sequence is already linear"));
        }
        let mut edit = if k == 0 {
            self.topology_edit(false)
        } else {
            self.rotation(k)?
        };
        edit.topology = Some((gb_io::seq::Topology::Circular, gb_io::seq::Topology::Linear));
        edit.note = "made linear";
        self.push(edit);
        Ok(())
    }

    /// Circular <-> linear as one undo step (the Edit menu toggle). Nothing else changes.
    pub fn set_topology(&mut self, circular: bool) {
        if (self.record.topology == gb_io::seq::Topology::Circular) != circular {
            let edit = self.topology_edit(circular);
            self.push(edit);
        }
    }

    fn topology_edit(&self, circular: bool) -> Edit {
        use gb_io::seq::Topology::{Circular, Linear};
        let (before, after) = if circular { (Linear, Circular) } else { (Circular, Linear) };
        Edit {
            pos: 0,
            removed: Vec::new(),
            inserted: Vec::new(),
            features_before: self.record.features.clone(),
            features_after: None,
            note: if circular { "made circular" } else { "made linear" },
            summary: None,
            topology: Some((before, after)),
            remap: Remap::Local,
            memo: EditedMemo::default(),
        }
    }

    /// The edit for `rotate(k)`, not yet applied.
    fn rotation(&self, k: usize) -> io::Result<Edit> {
        let n = self.len();
        if k == 0 || k >= n {
            return Err(io::Error::new(io::ErrorKind::InvalidInput, format!("cannot start at base {}", k + 1)));
        }
        let mut features: Vec<Feature> = self
            .record
            .features
            .iter()
            .filter(|f| f.kind == "source")
            .map(|f| Feature { location: Location::simple_range(0, n as i64), ..f.clone() })
            .collect();
        features.extend(self.features_across(k, k + n));
        let mut seq = self.record.seq[k..].to_vec();
        seq.extend_from_slice(&self.record.seq[..k]);
        Ok(Edit {
            pos: 0,
            removed: self.record.seq.clone(),
            inserted: seq,
            features_before: self.record.features.clone(),
            features_after: Some(features),
            note: "rotated",
            summary: None,
            topology: None,
            remap: Remap::Rotate(k),
            memo: EditedMemo::default(),
        })
    }

    /// A feature-only edit: no bases change.
    fn push_features(&mut self, after: Vec<Feature>) {
        self.push(Edit {
            pos: 0,
            removed: Vec::new(),
            inserted: Vec::new(),
            features_before: self.record.features.clone(),
            features_after: Some(after),
            note: "",
            summary: None,
            topology: None,
            remap: Remap::Local,
            memo: EditedMemo::default(),
        });
    }

    /// The history: summaries of the undo steps (oldest first) and of the redo steps
    /// (the next one to redo first).
    pub fn history(&self) -> (Vec<Summary>, Vec<Summary>) {
        let summary = |e: &Edit| e.summary.clone().unwrap_or_else(|| summarize(e));
        (self.undo_stack.iter().map(summary).collect(), self.redo_stack.iter().rev().map(summary).collect())
    }

    /// Apply a new edit and record it for undo.
    fn push(&mut self, mut edit: Edit) {
        // A new edit below the saved point makes the saved state unreachable by undo/redo.
        if self.saved_depth.is_some_and(|d| d > self.real_depth()) {
            self.saved_depth = None;
        }
        if edit.summary.is_none() {
            edit.summary = Some(summarize(&edit)); // once, while it is cheap to say
        }
        let before = self.edited.clone();
        self.apply(&edit);
        edit.memo = EditedMemo { epoch: self.edited_epoch, before, after: self.edited.clone() };
        self.undo_stack.push(edit);
        self.redo_stack.clear();
    }

    /// Returns false when there is nothing to undo.
    pub fn undo(&mut self) -> bool {
        // `let Some(x) = ... else { return ... }` = unpack or leave early.
        let Some(mut edit) = self.undo_stack.pop() else { return false };
        let end = edit.pos + edit.inserted.len();
        // `splice` replaces a slice of the Vec in place (like Python's seq[a:b] = new).
        self.record.seq.splice(edit.pos..end, edit.removed.iter().copied());
        self.gc_blocks = OnceLock::new();
        if let Some((before, _)) = &edit.topology {
            self.record.topology = before.clone();
        }
        // Swap restores the old feature list and keeps the current one (unused) in `edit`.
        std::mem::swap(&mut self.record.features, &mut edit.features_before);
        self.record.len = Some(self.record.seq.len());
        let n = self.record.seq.len();
        if edit.memo.epoch == self.edited_epoch {
            self.edited = edit.memo.before.clone(); // no save since: restore exactly
        } else {
            remap_edited(&mut self.edited, &edit, false, n);
        }
        // Redo re-applies the edit to the restored features, so give it the restored list.
        edit.features_before = self.record.features.clone();
        self.redo_stack.push(edit);
        self.edited_if_clean();
        true
    }

    /// Returns false when there is nothing to redo.
    pub fn redo(&mut self) -> bool {
        let Some(edit) = self.redo_stack.pop() else { return false };
        self.apply(&edit);
        if edit.memo.epoch == self.edited_epoch {
            self.edited = edit.memo.after.clone(); // no save since: restore exactly
        }
        self.undo_stack.push(edit);
        self.edited_if_clean();
        true
    }

    /// Change sequence + features. `&Edit` = we only read the edit description.
    fn apply(&mut self, edit: &Edit) {
        let start = edit.pos;
        let end = start + edit.removed.len();
        self.record.seq.splice(start..end, edit.inserted.iter().copied());
        self.record.len = Some(self.record.seq.len());
        self.gc_blocks = OnceLock::new();
        remap_edited(&mut self.edited, edit, true, self.record.seq.len());
        if let Some((_, after)) = &edit.topology {
            self.record.topology = after.clone();
        }
        if let Some(features) = &edit.features_after {
            self.record.features = features.clone();
            return;
        }

        move_features(&mut self.record.features, start as i64, end as i64, edit.inserted.len() as i64);
    }

    /// Indices of features overlapping `[start, end)`.
    ///
    /// shortcut: linear scan, ~microseconds for 10k features; switch to an interval tree
    /// when loading tracks with millions of entries (GFF/BAM, Phase 4).
    pub fn features_in(&self, start: i64, end: i64) -> Vec<usize> {
        self.record
            .features
            .iter()
            .enumerate() // yields (index, feature) pairs, like Python's enumerate()
            .filter(|(_, f)| locations::spans(&f.location).0.iter().any(|&(a, b)| a < end && b > start))
            .map(|(i, _)| i)
            .collect()
    }
}

/// Move feature locations for "bases [s, e) replaced by n new bases".
fn move_features(features: &mut Vec<Feature>, s: i64, e: i64, n: i64) {
    // `retain_mut` walks the features, edits each in place, and drops the ones whose
    // closure returns false (fully deleted features).
    features.retain_mut(|f| {
        if e > s {
            match locations::after_delete(&f.location, s, e) {
                Some(l) => f.location = l,
                None => return false,
            }
        }
        if n > 0 {
            f.location = locations::after_insert(&f.location, s, n);
        }
        true
    });
}

/// Validate a feature coming from the user (type is one word, location inside a
/// sequence of `len` bases) and canonicalize its qualifier text.
fn check_feature(feature: &mut Feature, len: i64) -> io::Result<()> {
    let bad = |msg: String| io::Error::new(io::ErrorKind::InvalidInput, msg);
    if feature.kind.is_empty() || feature.kind.contains(char::is_whitespace) {
        return Err(bad(format!("feature type '{}' must be one word", feature.kind)));
    }
    let (spans, _) = locations::spans(&feature.location);
    if spans.is_empty() || spans.iter().any(|&(s, e)| s < 0 || e < s || e > len) {
        return Err(bad(format!("location {} is outside the {len} bp sequence", feature.location)));
    }
    for (key, value) in &mut feature.qualifiers {
        if let Some(v) = value {
            *v = canonical(key, v);
        }
    }
    Ok(())
}

/// GenBank wraps long values over several lines, and gb-io keeps those line breaks inside
/// the value; its writer then wraps again around them, so text drifted on every save.
/// We store values in canonical form: every run of whitespace (spaces, line breaks) becomes
/// one space, as Biopython does; `/translation` (a protein sequence) loses all whitespace.
/// The patched writer (vendor/gb-io) wraps only at single spaces and never cuts a word
/// (except in `/translation`), so save -> load gives back exactly the same text.
/// Cost: intentional double spaces in notes become single spaces.
fn unwrap_text(record: &mut Seq) {
    if let Some(d) = &mut record.definition {
        *d = canonical("", d);
    }
    for f in &mut record.features {
        // `&mut` in the pattern: we get mutable access to each (key, value) pair in place.
        for (key, value) in &mut f.qualifiers {
            if let Some(v) = value {
                *v = canonical(key, v);
            }
        }
    }
}

/// Canonical form of one qualifier value (see `unwrap_text`).
fn canonical(key: &str, value: &str) -> String {
    // `split_whitespace` splits on any run of spaces/tabs/newlines and drops empty pieces.
    let words = value.split_whitespace();
    if key == "translation" { words.collect() } else { words.collect::<Vec<_>>().join(" ") }
}

/// Letters accepted in a DNA sequence: ACGTU, IUPAC ambiguity codes, N and gap '-'.
fn is_iupac(b: u8) -> bool {
    matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T' | b'U' | b'R' | b'Y' | b'S' | b'W'
        | b'K' | b'M' | b'B' | b'D' | b'H' | b'V' | b'N' | b'-')
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dna_file::tests::sample_dna;
    use gb_io::seq::Location;
    use std::borrow::Cow;

    fn sample() -> Document {
        let (record, primers) = dna_file::parse(&sample_dna()).unwrap();
        Document::new(record, primers)
    }

    #[test]
    fn edit_moves_features_and_undo_restores_exactly() {
        let mut doc = sample();
        let original = doc.record.clone();
        doc.replace(4, 4, b"TTT").unwrap(); // insert inside lacZ (2..10)
        assert_eq!(doc.record.features[0].location.to_gb_format(), "2..13");
        assert_eq!(doc.len(), 23);
        assert!(doc.is_modified());
        assert!(doc.undo());
        assert_eq!(doc.record, original);
        assert!(!doc.is_modified());
        assert!(doc.redo());
        assert_eq!(doc.record.features[0].location.to_gb_format(), "2..13");
    }

    #[test]
    fn rejects_bad_edits() {
        let mut doc = sample();
        assert!(doc.replace(5, 4, b"").is_err());
        assert!(doc.replace(0, 999, b"").is_err());
        assert!(doc.replace(0, 0, b"AXG").is_err());
        assert!(!doc.undo()); // failed edits leave no history
    }

    #[test]
    fn wrapped_qualifiers_are_stable_across_saves() {
        // r#"..."# = raw string: backslashes and quotes are taken literally, spaces kept.
        let gb = r#"LOCUS       t                         70 bp    DNA     linear   UNK 01-JAN-1980
FEATURES             Location/Qualifiers
     CDS             1..9
                     /note="a long note that the original writer wrapped
                     over two lines.
                     "
                     /translation="MKVLAAGIVALLLAAGCSSSKEETPAPAAAPAPAAEAPAAPAAEAP
                     AAPAAE"
ORIGIN
        1 atgaaagtgc tggcagcagg tatagtagca ctgctactgg cagcaggttg cagcagcagc
       61 aaagaagaaa
//
"#;
        let path = std::env::temp_dir().join("ge_test_wrapped.gb");
        std::fs::write(&path, gb).unwrap();
        let mut doc = Document::open(&path).unwrap();
        let f = &doc.record.features[0];
        assert_eq!(f.qualifier_values("note").next(), Some("a long note that the original writer wrapped over two lines."));
        assert_eq!(
            f.qualifier_values("translation").next(),
            Some("MKVLAAGIVALLLAAGCSSSKEETPAPAAAPAPAAEAPAAPAAEAPAAPAAE")
        );
        // Save twice: the text must not drift.
        doc.save_genbank(&path, None).unwrap();
        let mut again = Document::open(&path).unwrap();
        assert_eq!(again.record.features, doc.record.features);
        again.save_genbank(&path, None).unwrap();
        assert_eq!(Document::open(&path).unwrap().record.features, doc.record.features);
    }

    #[test]
    fn quotes_survive_at_every_wrap_position() {
        // Regression test for the patched gb-io writer: slide a quoted word across the
        // wrap column; upstream dropped the quote or panicked at some offsets.
        let path = std::env::temp_dir().join("ge_test_quotes.gb");
        for k in 0..120 {
            let mut doc = sample();
            let note = format!("w{} \"quoted\" end", "x".repeat(k));
            doc.record.features[0].qualifiers.push((Cow::Borrowed("note"), Some(note.clone())));
            doc.save_genbank(&path, None).unwrap();
            let back = Document::open(&path).unwrap();
            let got: Vec<&str> = back.record.features[0].qualifier_values("note").collect();
            assert!(got.contains(&note.as_str()), "offset {k}: {got:?}");
        }
    }

    #[test]
    fn feature_edits_are_undoable() {
        let mut doc = sample();
        let original = doc.record.features.clone();
        let location = Location::from_gb_format("complement(3..8)").unwrap();
        let note = vec![("note".into(), Some("a\n  b".to_string()))];
        let new = Feature { kind: "misc_feature".into(), location, qualifiers: note };
        let i = doc.set_feature(None, new.clone()).unwrap();
        assert_eq!(i, original.len());
        assert_eq!(doc.record.features[i].qualifier_values("note").next(), Some("a b"));
        doc.delete_feature(0).unwrap();
        assert_eq!(doc.record.features.len(), original.len());
        // `..new.clone()` = "every other field copied from new".
        let outside = Feature { location: Location::simple_range(0, 999), ..new.clone() };
        assert!(doc.set_feature(None, outside).is_err());
        assert!(doc.set_feature(Some(99), new).is_err());
        assert!(doc.undo() && doc.undo() && !doc.undo());
        assert_eq!(doc.record.features, original);
        assert!(doc.redo());
        assert_eq!(doc.record.features.len(), original.len() + 1);
    }

    #[test]
    fn several_feature_changes_are_one_undo_step() {
        let mut doc = sample();
        let original = doc.record.features.clone();
        let gene = Feature { kind: "gene".into(), location: Location::simple_range(1, 10), qualifiers: vec![] };
        let renamed = Feature { kind: "CDS".into(), ..original[1].clone() };
        // Delete #0, replace #1 (still "#1" although #0 goes), add one.
        doc.edit_features(vec![(Some(0), None), (Some(1), Some(renamed.clone())), (None, Some(gene.clone()))]).unwrap();
        // Kept in their old order (#1 replaced in place), the added one last.
        let expected = [vec![renamed], original[2..].to_vec(), vec![gene]].concat();
        assert_eq!(doc.record.features, expected);
        assert!(doc.undo());
        assert_eq!(doc.record.features, original);
        // An invalid change applies nothing.
        assert!(doc.edit_features(vec![(Some(0), None), (Some(99), None)]).is_err());
        assert_eq!(doc.record.features, original);
    }

    #[test]
    fn refreshing_translations_is_one_undo_step_and_only_touches_stale_cds() {
        // ATG AAA TAA = M K stop, as a CDS with a wrong /translation; a second CDS without one; a pseudogene.
        let mut doc = Document::new_sequence("t", b"ATGAAATAAATGAAATAAATGAAATAA", false).unwrap();
        let cds = |loc: &str, quals: Vec<(&str, Option<&str>)>| Feature {
            kind: "CDS".into(),
            location: Location::from_gb_format(loc).unwrap(),
            qualifiers: quals.into_iter().map(|(k, v)| (k.to_string().into(), v.map(String::from))).collect(),
        };
        doc.add_features(vec![
            cds("1..9", vec![("translation", Some("MQ"))]),
            cds("10..18", vec![]),
            cds("19..27", vec![("pseudo", None), ("translation", Some("MQ"))]),
        ])
        .unwrap();
        let before = doc.record.features.clone();
        let depth = doc.history().0.len();
        let report = doc.refresh_translations(None, 1);
        // Only the CDS that has a /translation is looked at: #0 updated, #2 skipped; #1 never gets one.
        assert_eq!(report.len(), 2);
        assert_eq!(report[0], (0, Refresh::Updated));
        assert!(matches!(report[1], (2, Refresh::Skipped(crate::code::Skip::Pseudo, _))));
        assert_eq!(doc.record.features[0].qualifier_values("translation").next(), Some("MK"));
        assert_eq!(doc.record.features[1].qualifier_values("translation").next(), None);
        assert_eq!(doc.record.features[2].qualifier_values("translation").next(), Some("MQ"));
        // One undo step, and a second run finds nothing to do (so adds no step).
        assert_eq!(doc.history().0.len(), depth + 1);
        assert_eq!(doc.refresh_translations(None, 1)[0], (0, Refresh::Unchanged));
        assert_eq!(doc.history().0.len(), depth + 1);
        assert!(doc.undo());
        assert_eq!(doc.record.features, before);
    }

    #[test]
    fn copy_paste_carries_features() {
        let mut doc = sample(); // lacZ 2..10 forward, ori 18..3 across the origin
        let copied = doc.features_between(4, 12); // half of lacZ
        let lacz = copied.iter().find(|f| f.qualifier_values("label").next() == Some("lacZ")).unwrap();
        assert_eq!(lacz.location.to_gb_format(), "1..6"); // bases 5..10 -> relative 1..6
        let n = doc.record.features.len();
        let bases = doc.record.seq[4..12].to_vec();
        doc.paste(0, 0, &bases, copied.clone()).unwrap();
        assert_eq!(doc.len(), 28);
        assert_eq!(doc.record.features.len(), n + copied.len());
        assert_eq!(doc.record.features[n].location.to_gb_format(), "1..6");
        assert!(doc.undo());
        assert_eq!(doc.record.features.len(), n);
        // Features that do not fit the pasted bases are refused.
        assert!(doc.paste(0, 0, b"AC", copied).is_err());
    }

    #[test]
    fn reverse_complement_twice_is_identity() {
        let mut doc = sample();
        let original = doc.record.clone();
        doc.reverse_complement();
        // lacZ was 2..10 (forward) in 20 bp -> complement(11..19)
        assert_eq!(doc.record.features[0].location.to_gb_format(), "complement(11..19)");
        assert_eq!(doc.record.seq, b"ACGTAAAACCCCGGGGTTTT"); // sample is AAAACCCCGGGGTTTTACGT
        doc.reverse_complement();
        assert_eq!(doc.record, original);
        assert!(doc.undo() && doc.undo());
        assert_eq!(doc.record, original);
    }

    #[test]
    fn pcr_product_with_tails_and_across_the_origin() {
        let doc = sample(); // 20 bp circular: lacZ 2..10, ori complement(join(18..20,1..3))
        let p = doc.pcr_product("prod", 1, 10, b"GG", b"CC").unwrap();
        assert_eq!(p.record.seq, b"GGAAACCCCGGCC"); // GG + AAACCCCGG + CC
        let lacz = &p.record.features[0];
        assert_eq!(lacz.location.to_gb_format(), "3..11"); // shifted by the 2-base tail
        assert!(p.is_modified() && p.record.topology == gb_io::seq::Topology::Linear);
        // 16..24 on a 20 bp circle = bases 17-20 then 1-4.
        let wrap = doc.pcr_product("w", 16, 24, b"", b"").unwrap();
        assert_eq!(wrap.record.seq, b"ACGTAAAA");
        assert!(doc.pcr_product("bad", 5, 5, b"", b"").is_err());
        // Substitution keeps a feature that starts at the changed base intact.
        let mut q = doc.pcr_product("q", 1, 10, b"", b"").unwrap();
        let before = q.record.features[0].location.clone();
        q.overwrite(0, b"T").unwrap();
        assert_eq!(q.record.seq[0], b'T');
        assert_eq!(q.record.features[0].location, before);
        assert!(doc.pcr_product("bad", 0, 4, b"XX", b"").is_err());
    }

    #[test]
    fn saving_replaces_primer_sites_without_touching_the_document() {
        let mut doc = sample(); // has one primer_bind feature (p1)
        let before = doc.record.clone();
        let site = Feature {
            kind: "primer_bind".into(),
            location: Location::from_gb_format("2..12").unwrap(),
            qualifiers: vec![("label".into(), Some("new".to_string()))],
        };
        let path = std::env::temp_dir().join("ge_test_primers.gb");
        doc.save_genbank(&path, Some(vec![site])).unwrap();
        assert_eq!(doc.record, before); // in memory: unchanged
        assert!(!doc.is_modified());
        let back = Document::open(&path).unwrap();
        let primers: Vec<_> = back.record.features.iter().filter(|f| f.kind == "primer_bind").collect();
        assert_eq!(primers.len(), 1);
        assert_eq!(primers[0].qualifier_values("label").next(), Some("new"));
        assert_eq!(back.record.features.len(), before.features.len()); // others kept
        // A site outside the sequence is refused, nothing written.
        let bad = Feature { location: Location::simple_range(0, 999), ..primers[0].clone() };
        assert!(doc.save_genbank(&path, Some(vec![bad])).is_err());
    }

    #[test]
    fn gc_bins_exact_and_from_blocks() {
        let mut doc = sample(); // AAAACCCCGGGGTTTTACGT: GC 10/20
        assert_eq!(doc.gc_bins(0, 20, 2), vec![0.6, 0.4]); // AAAACCCCGG | GGTTTTACGT
        assert_eq!(doc.gc_bins(0, 8, 1), vec![0.5]);
        // A long sequence read through the block cache, then edited (cache dropped).
        doc.replace(0, 20, &b"GC".repeat(20_000)).unwrap();
        assert_eq!(doc.gc_bins(0, 40_000, 2), vec![1.0, 1.0]);
        doc.replace(0, 40_000, &b"AT".repeat(20_000)).unwrap();
        assert_eq!(doc.gc_bins(0, 40_000, 2), vec![0.0, 0.0]);
        assert!(doc.undo());
        assert_eq!(doc.gc_bins(0, 40_000, 1), vec![1.0]);
    }

    #[test]
    fn genbank_round_trip() {
        let mut doc = sample();
        doc.replace(0, 2, b"GG").unwrap();
        let path = std::env::temp_dir().join("ge_test_roundtrip.gb");
        doc.save_genbank(&path, None).unwrap();
        assert!(!doc.is_modified());
        let back = Document::open(&path).unwrap();
        assert_eq!(back.record.seq, doc.record.seq);
        assert_eq!(back.record.features, doc.record.features);
    }

    #[test]
    fn copy_and_edit_across_the_origin() {
        let mut doc = sample(); // 20 bp circular, ori = complement(join(18..20,1..3))
        let original = doc.record.clone();
        // Copy 16..25 (passes the origin): ori comes along as ONE feature, its 3 + 3 bases
        // now side by side, so one plain range.
        let copied = doc.features_across(15, 25);
        let ori: Vec<_> = copied.iter().filter(|f| locations::spans(&f.location).1 < 0).collect();
        assert_eq!(ori.len(), 1);
        let (spans, strand) = locations::spans(&ori[0].location);
        assert_eq!((spans, strand), (vec![(2, 8)], -1));
        // Replace 18..23 (bases 18-20 and 1-3) with GG: one undo step, undo is exact.
        let at = doc.edit_across_origin(17, 23, b"GG", None).unwrap();
        assert_eq!((at, doc.len()), (14, 16));
        assert_eq!(&doc.record.seq[14..16], b"GG");
        assert!(doc.undo());
        assert_eq!(doc.record, original);
        assert!(!doc.undo() && doc.redo() && doc.len() == 16);
        assert!(doc.edit_across_origin(5, 10, b"", None).is_err()); // does not pass the origin
    }

    #[test]
    fn rotate_moves_bases_and_features() {
        let mut doc = sample(); // 20 bp circular, lacZ 2..10, ori complement(join(18..20,1..3))
        let original = doc.record.clone();
        doc.rotate(5).unwrap(); // base 6 becomes base 1
        assert_eq!(&doc.record.seq[..15], &original.seq[5..]);
        assert_eq!(&doc.record.seq[15..], &original.seq[..5]);
        // ori (bases 18-20 and 1-3) no longer crosses the start: one piece, 13..18.
        let ori = doc.record.features.iter().find(|f| locations::spans(&f.location).1 < 0).unwrap();
        assert_eq!(locations::spans(&ori.location).0, vec![(12, 18)]);
        assert!(doc.undo());
        assert_eq!(doc.record, original);
        assert!(doc.rotate(0).is_err());
    }

    #[test]
    fn history_says_what_each_step_did() {
        let mut doc = sample(); // 20 bp, lacZ 2..10
        doc.replace(4, 7, b"GAG").unwrap(); // replace 3 bases
        doc.replace(0, 0, b"TT").unwrap(); // insert 2
        doc.delete_feature(0).unwrap();
        let (undo, redo) = doc.history();
        assert_eq!(undo.len(), 3);
        let s = &undo[0];
        assert_eq!((s.pos, s.removed_len, s.inserted.as_str()), (4, 3, "GAG"));
        assert_eq!((undo[1].pos, undo[1].inserted.as_str(), undo[1].removed_len), (0, "TT", 0));
        assert_eq!((undo[2].deleted.len(), undo[2].added.len()), (1, 0));
        assert!(redo.is_empty());
        assert!(doc.undo());
        let (undo, redo) = doc.history();
        assert_eq!((undo.len(), redo.len(), redo[0].deleted.len()), (2, 1, 1));
    }

    // ================================================================================
    // Sequence editing (docs/TESTING.md section 1) and annotation drift (section 2).
    // ================================================================================

    /// Tiny deterministic random numbers (no extra crate): the same seed always gives the
    /// same numbers, so a failing test can be replayed exactly.
    struct Lcg(u64);

    impl Lcg {
        /// A number in `0..m`.
        fn below(&mut self, m: usize) -> usize {
            self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            ((self.0 >> 33) as usize) % m.max(1)
        }


        fn bases(&mut self, n: usize) -> Vec<u8> {
            (0..n).map(|_| b"ACGT"[self.below(4)]).collect()
        }

        /// 0 to m - 1 random bases.
        fn some_bases(&mut self, m: usize) -> Vec<u8> {
            let n = self.below(m);
            self.bases(n)
        }
    }

    /// A feature of `kind` at GenBank (1-based) location text `loc`, with a `/label`.
    fn feat(kind: &str, loc: &str, label: &str) -> Feature {
        Feature {
            kind: kind.to_string().into(),
            location: Location::from_gb_format(loc).unwrap(),
            qualifiers: vec![("label".into(), Some(label.to_string()))],
        }
    }

    /// A document holding `seq` and exactly these features. They are set directly, so the
    /// document has no undo history and is not "modified".
    fn doc_with(seq: &[u8], circular: bool, features: Vec<Feature>) -> Document {
        let mut doc = Document::new_sequence("t", seq, circular).unwrap();
        doc.record.features = features;
        doc.set_modified(false);
        doc
    }

    fn rc(bases: &[u8]) -> Vec<u8> {
        bases
            .iter()
            .rev()
            .map(|&c| match c.to_ascii_uppercase() {
                b'A' => b'T',
                b'T' => b'A',
                b'C' => b'G',
                b'G' => b'C',
                other => other,
            })
            .collect()
    }

    fn index_of(doc: &Document, label: &str) -> usize {
        doc.record
            .features
            .iter()
            .position(|f| f.qualifier_values("label").next() == Some(label))
            .unwrap_or_else(|| panic!("feature {label} is gone"))
    }

    fn has(doc: &Document, label: &str) -> bool {
        doc.record.features.iter().any(|f| f.qualifier_values("label").next() == Some(label))
    }

    fn spans_of(doc: &Document, label: &str) -> Vec<(i64, i64)> {
        locations::spans(&doc.record.features[index_of(doc, label)].location).0
    }

    fn strand_of(doc: &Document, label: &str) -> i8 {
        locations::spans(&doc.record.features[index_of(doc, label)].location).1
    }

    /// The top-strand bases a feature covers, part after part in location order.
    fn covered(doc: &Document, label: &str) -> Vec<u8> {
        spans_of(doc, label)
            .iter()
            .flat_map(|&(a, b)| doc.record.seq[a as usize..b as usize].iter().copied())
            .collect()
    }

    /// The feature's own sequence (reverse-complemented when it is on the bottom strand).
    fn bio(doc: &Document, label: &str) -> Vec<u8> {
        let c = covered(doc, label);
        if strand_of(doc, label) < 0 { rc(&c) } else { c }
    }

    // ---------------------------------- section 1 -----------------------------------

    #[test]
    fn empty_edits_change_nothing() {
        let mut doc = doc_with(&Lcg(1).bases(50), true, vec![feat("gene", "11..20", "g")]);
        let before = doc.record.clone();
        for i in [0, 7, 50] {
            doc.replace(i, i, b"").unwrap(); // delete(i, i) and an empty insert
            doc.overwrite(i, b"").unwrap();
        }
        assert!(!doc.is_modified());
        assert!(!doc.undo(), "an empty edit left an undo step");
        assert_eq!(doc.record, before);
    }

    /// The big one: random edits of every kind, undo, redo and save, checked at every step
    /// against (a) a plain `Vec<u8>` model of the bases and (b) a list of whole-record
    /// snapshots (one per undo step) that undo / redo must walk exactly.
    #[test]
    fn random_edits_undo_redo_follow_a_snapshot_and_string_model() {
        for circular in [false, true] {
            let mut rng = Lcg(if circular { 11 } else { 12 });
            let seq = rng.bases(120);
            let mut features = vec![
                feat("gene", "10..40", "a"),
                feat("CDS", "complement(join(50..60,70..80))", "b"),
                feat("misc_feature", "<90..>110", "c"),
            ];
            if circular {
                features.push(feat("rep_origin", "join(115..120,1..5)", "o"));
            }
            let mut doc = doc_with(&seq, circular, features);
            let mut states = vec![doc.record.clone()]; // states[i] = the record after i steps
            let mut pos = 0usize; // which state we are in
            let mut saved = Some(0usize); // state number at the last save
            let mut applied = 0;
            for step in 0..600 {
                let cur = doc.record.seq.clone();
                let n = cur.len();
                let mut model: Option<Vec<u8>> = None; // what the bases must be afterwards
                let mut pushed = false; // did this make a new undo step?
                let op = rng.below(15);
                let what;
                match op {
                    0..=4 => {
                        let s = rng.below(n + 1);
                        let e = (s + rng.below(6)).min(n);
                        let ins = rng.some_bases(6);
                        what = format!("replace({s},{e},{})", String::from_utf8_lossy(&ins));
                        doc.replace(s, e, &ins).unwrap();
                        pushed = !(s == e && ins.is_empty());
                        model = Some([&cur[..s], &ins[..], &cur[e..]].concat());
                    }
                    5 | 6 => {
                        let s = rng.below(n + 1);
                        let k = rng.below(4).min(n - s);
                        let ins = rng.bases(k);
                        what = format!("overwrite({s},{})", String::from_utf8_lossy(&ins));
                        doc.overwrite(s, &ins).unwrap();
                        pushed = k > 0;
                        model = Some([&cur[..s], &ins[..], &cur[s + k..]].concat());
                    }
                    7 | 8 if n >= 4 => {
                        // Copy a stretch of this very document and paste it over a selection.
                        let a = rng.below(n - 1);
                        let b = a + 1 + rng.below((n - a - 1).min(30));
                        let copied = doc.features_between(a, b);
                        let s = rng.below(n + 1);
                        let e = (s + rng.below(8)).min(n);
                        what = format!("paste({s},{e}) of {a}..{b} with {} features", copied.len());
                        doc.paste(s, e, &cur[a..b], copied).unwrap();
                        pushed = true;
                        model = Some([&cur[..s], &cur[a..b], &cur[e..]].concat());
                    }
                    9 if circular && n >= 2 => {
                        let k = 1 + rng.below(n - 1);
                        what = format!("rotate({k})");
                        doc.rotate(k).unwrap();
                        pushed = true;
                        model = Some([&cur[k..], &cur[..k]].concat());
                    }
                    10 if n >= 1 => {
                        what = "reverse_complement".to_string();
                        doc.reverse_complement();
                        pushed = true;
                        model = Some(rc(&cur));
                    }
                    11 if n >= 2 => {
                        let s = rng.below(n);
                        let e = s + 1 + rng.below(n - s);
                        what = format!("add feature {}..{}", s + 1, e);
                        doc.set_feature(None, feat("misc_feature", &format!("{}..{}", s + 1, e), "x")).unwrap();
                        pushed = true;
                        model = Some(cur.clone());
                    }
                    12 if circular && n >= 12 => {
                        let start = n - 1 - rng.below(5);
                        let end = n + 1 + rng.below(5);
                        let bases = rng.some_bases(5);
                        what = format!("edit_across_origin({start},{end},{})", String::from_utf8_lossy(&bases));
                        let at = doc.edit_across_origin(start, end, &bases, None).unwrap();
                        assert_eq!(at, start - (end - n), "{what}: where the new bases start");
                        pushed = true;
                        model = Some([&cur[end - n..start], &bases[..]].concat());
                    }
                    13 => {
                        what = "undo".to_string();
                        assert_eq!(doc.undo(), pos > 0, "step {step} {what}: undo availability");
                        pos = pos.saturating_sub(1);
                    }
                    14 => {
                        if rng.below(3) == 0 {
                            what = "save".to_string();
                            doc.set_modified(false);
                            saved = Some(pos);
                        } else {
                            what = "redo".to_string();
                            let can = pos + 1 < states.len();
                            assert_eq!(doc.redo(), can, "step {step} {what}: redo availability");
                            if can {
                                pos += 1;
                            }
                        }
                    }
                    _ => {
                        what = "(skipped)".to_string();
                    }
                }
                let ctx = format!("circular={circular} step {step} after {what}");
                if pushed {
                    if saved.is_some_and(|d| d > pos) {
                        saved = None; // the saved state can no longer be reached
                    }
                    states.truncate(pos + 1); // a new edit drops the redo stack
                    states.push(doc.record.clone());
                    pos += 1;
                    applied += 1;
                    assert!(!doc.redo(), "{ctx}: redo available after a new edit");
                }
                if let Some(m) = model {
                    assert_eq!(doc.record.seq, m, "{ctx}: bases differ from the string model");
                }
                assert_eq!(doc.record, states[pos], "{ctx}: record differs from snapshot {pos}");
                assert_eq!(doc.record.len, Some(doc.record.seq.len()), "{ctx}: stale length");
                assert_eq!(doc.is_modified(), saved != Some(pos), "{ctx}: modified flag");
            }
            assert!(applied > 100, "the random run made too few edits ({applied})");
            // Redo everything, then undo everything back to the original, then redo again.
            while doc.redo() {
                pos += 1;
            }
            assert_eq!(doc.record, states[pos]);
            while pos > 0 {
                assert!(doc.undo());
                pos -= 1;
                assert_eq!(doc.record, states[pos], "circular={circular}: undoing back to state {pos}");
            }
            assert!(!doc.undo());
            assert_eq!(doc.record.seq, seq);
            while doc.redo() {
                pos += 1;
            }
            assert_eq!(doc.record, states[pos], "redo of everything");
        }
    }

    #[test]
    fn circular_edit_across_the_origin_equals_the_rotated_string_edit() {
        let n = 60;
        let seq = Lcg(5).bases(n);
        let cases: [(usize, usize, &[u8]); 7] =
            [(55, 63, b"GGG"), (50, 61, b""), (59, 61, b"TTTTT"), (30, 90, b""), (30, 90, b"A"), (59, 61, b""), (10, 70, b"CC")];
        for (start, end, new) in cases {
            let mut doc = doc_with(&seq, true, vec![feat("gene", "11..30", "x"), feat("gene", "36..50", "y")]);
            let before = doc.record.clone();
            let at = doc.edit_across_origin(start, end, new, None).unwrap();
            // The string model: rotate so `start` is first, replace the first end-start
            // bases, which leaves new + the untouched rest (a rotation of what Rust holds).
            let rot = [&seq[start..], &seq[..start]].concat();
            let model = [new, &rot[end - start..]].concat();
            let got = doc.record.seq.clone();
            assert_eq!(got.len(), model.len(), "{start}..{end}: length");
            assert!(
                model.is_empty() || [&got[..], &got[..]].concat().windows(model.len()).any(|w| w == &model[..]),
                "{start}..{end}: {got:?} is not a rotation of {model:?}"
            );
            assert_eq!(&got[at..at + new.len()], new, "{start}..{end}: returned position");
            assert_eq!(doc.record.topology, gb_io::seq::Topology::Circular);
            assert_eq!(doc.history().0.len(), 1, "{start}..{end}: one undo step");
            assert!(doc.undo() && !doc.undo());
            assert_eq!(doc.record, before, "{start}..{end}: undo is exact");
        }
    }

    #[test]
    fn circular_inserts_at_zero_and_at_the_end_and_wrapping_deletes() {
        let n = 60;
        let seq = Lcg(6).bases(n);
        let mut doc = doc_with(&seq, true, vec![feat("gene", "11..20", "g"), feat("rep_origin", "join(56..60,1..5)", "o")]);
        let before = doc.record.clone();
        let o = [&seq[55..], &seq[..5]].concat();
        let g = seq[10..20].to_vec();
        for at in [0usize, n] {
            doc.replace(at, at, b"TTT").unwrap();
            assert_eq!(doc.len(), n + 3);
            assert_eq!(doc.record.topology, gb_io::seq::Topology::Circular);
            assert_eq!(covered(&doc, "o"), o, "insert at {at} moved the origin-crossing feature");
            assert_eq!(covered(&doc, "g"), g, "insert at {at} moved gene g");
            assert_eq!(doc.record.seq.len(), doc.record.len.unwrap());
            assert!(doc.undo());
            assert_eq!(doc.record, before);
        }
        // A delete wrapping the origin: 0-based 57..60 and 0..3.
        doc.edit_across_origin(57, 63, b"", None).unwrap();
        assert_eq!(doc.len(), n - 6);
        assert_eq!(doc.record.seq, seq[3..57].to_vec());
        assert_eq!(doc.record.topology, gb_io::seq::Topology::Circular);
        assert_eq!(covered(&doc, "g"), g);
        // o lost 3 bases at each side of the origin, keeps the other 2 + 2.
        assert_eq!(covered(&doc, "o"), [&seq[55..57], &seq[3..5]].concat());
        assert!(doc.undo());
        assert_eq!(doc.record, before);
    }

    #[test]
    fn paste_at_every_boundary_matches_the_string_model_and_is_one_step() {
        let seq = Lcg(7).bases(100);
        let base = doc_with(&seq, false, vec![feat("gene", "11..40", "g")]);
        let copied = base.features_between(20, 35); // part of g, relative: 1..15
        let bases = seq[20..35].to_vec();
        for (s, e) in [(0usize, 0usize), (0, 10), (100, 100), (90, 100), (0, 100), (40, 50), (40, 40), (45, 45), (12, 14)] {
            let mut doc = doc_with(&seq, false, vec![feat("gene", "11..40", "g")]);
            let before = doc.record.clone();
            doc.paste(s, e, &bases, copied.clone()).unwrap();
            assert_eq!(doc.record.seq, [&seq[..s], &bases[..], &seq[e..]].concat(), "paste({s},{e}): bases");
            // The pasted feature (last) covers exactly the pasted bases.
            let last = doc.record.features.last().unwrap();
            let (sp, _) = locations::spans(&last.location);
            let got: Vec<u8> = sp.iter().flat_map(|&(a, b)| doc.record.seq[a as usize..b as usize].to_vec()).collect();
            assert_eq!(got, bases, "paste({s},{e}): pasted feature covers the wrong bases");
            assert_eq!(sp, vec![(s as i64, (s + bases.len()) as i64)]);
            // The old feature, if it survives, still covers its own bases.
            if (s, e) == (0, 100) {
                assert_eq!(doc.record.features.len(), 1, "g should be gone");
            } else if (s, e) == (12, 14) {
                // Deleting 2 bases inside g and pasting 15 there: g loses 2 and gains 15.
                assert_eq!(covered(&doc, "g").len(), 30 - 2 + bases.len());
            } else {
                assert_eq!(covered(&doc, "g"), seq[10..40].to_vec(), "paste({s},{e}): g drifted");
            }
            assert_eq!(doc.history().0.len(), 1);
            assert!(doc.undo() && !doc.undo());
            assert_eq!(doc.record, before, "paste({s},{e}): undo");
        }
    }

    #[test]
    fn overwrite_keeps_length_moves_nothing_and_is_one_step_per_call() {
        let seq = Lcg(8).bases(80);
        let mut doc = doc_with(&seq, true, vec![feat("gene", "11..40", "g"), feat("x", "complement(join(50..60,70..80))", "h")]);
        let before = doc.record.clone();
        doc.overwrite(10, b"TTTT").unwrap(); // starts exactly at g's first base
        doc.overwrite(35, b"CCCCCC").unwrap(); // runs over g's end
        doc.overwrite(76, b"GGGG").unwrap(); // to the very end
        assert_eq!(doc.len(), 80);
        assert_eq!(doc.record.features, before.features, "an overwrite moved or trimmed a feature");
        assert_eq!(doc.history().0.len(), 3);
        assert!(doc.overwrite(78, b"ACG").is_err(), "overwrite past the end");
        assert_eq!(doc.history().0.len(), 3);
        assert!(doc.undo() && doc.undo() && doc.undo() && !doc.undo());
        assert_eq!(doc.record, before);
    }

    #[test]
    fn letters_are_kept_as_typed_and_invalid_ones_are_refused_not_dropped() {
        let mut doc = doc_with(b"ACGTACGT", false, vec![feat("gene", "2..5", "g")]);
        let before = doc.record.clone();
        let bad: [&[u8]; 6] = [b"A C", b"A1", b"A*", b"\xc3\xa9", b"A\n", b"A.G"];
        for b in bad {
            assert!(doc.replace(2, 2, b).is_err(), "{b:?} accepted by replace");
            assert!(doc.overwrite(0, b).is_err(), "{b:?} accepted by overwrite");
            assert!(doc.paste(0, 0, b, vec![]).is_err(), "{b:?} accepted by paste");
            assert!(Document::new_sequence("n", b, false).is_err(), "{b:?} accepted by new_sequence");
        }
        assert_eq!(doc.record, before, "a refused edit changed the document");
        assert!(!doc.undo(), "a refused edit left an undo step");
        // Accepted exactly as typed: lower case, U, every IUPAC code, N and the gap.
        let ok = b"acguRYSWKMBDHVNn-";
        doc.replace(4, 4, ok).unwrap();
        assert_eq!(doc.len(), 8 + ok.len(), "bases were dropped");
        assert_eq!(&doc.record.seq[4..4 + ok.len()], &ok[..], "bases were changed");
        assert!(doc.undo());
        assert_eq!(doc.record, before);
    }

    #[test]
    fn editing_an_empty_sequence() {
        for circular in [false, true] {
            let mut doc = Document::new_sequence("e", b"", circular).unwrap();
            assert_eq!(doc.len(), 0);
            doc.replace(0, 0, b"").unwrap(); // nothing
            assert!(!doc.undo());
            assert!(doc.replace(0, 1, b"").is_err() && doc.replace(1, 1, b"A").is_err());
            doc.replace(0, 0, b"ACG").unwrap();
            assert_eq!(doc.record.seq, b"ACG");
            assert!(doc.undo());
            assert_eq!(doc.len(), 0);
            assert_eq!(doc.record.len, Some(0));
            assert!(doc.redo() && doc.record.seq == b"ACG");
            doc.replace(0, 3, b"").unwrap(); // delete everything again
            assert_eq!(doc.len(), 0);
            assert!(doc.rotate(1).is_err() && doc.edit_across_origin(0, 1, b"", None).is_err());
            doc.reverse_complement(); // must not panic
            assert_eq!(doc.gc_bins(0, 0, 3).len(), 3);
            assert!(doc.features_between(0, 0).is_empty());
        }
    }

    #[test]
    fn modified_flag_follows_edits_undo_and_save() {
        let mut doc = sample();
        assert!(!doc.is_modified());
        doc.replace(0, 1, b"G").unwrap();
        assert!(doc.is_modified());
        assert!(doc.undo());
        assert!(!doc.is_modified(), "undoing to the opened state must clear modified");
        assert!(doc.redo() && doc.is_modified());
        let path = std::env::temp_dir().join("ge_test_modified.gb");
        doc.save_genbank(&path, None).unwrap();
        assert!(!doc.is_modified(), "saving clears modified");
        doc.replace(0, 0, b"T").unwrap();
        assert!(doc.is_modified());
        assert!(doc.undo() && !doc.is_modified(), "undo back to the saved state");
        assert!(doc.undo() && doc.is_modified(), "undo past the saved state");
        assert!(doc.redo() && !doc.is_modified(), "redo to the saved state");
        // A new edit after undoing below the saved state makes that state unreachable.
        assert!(doc.undo());
        doc.replace(3, 3, b"AA").unwrap();
        assert!(doc.is_modified());
        assert!(doc.undo() && doc.is_modified(), "the saved state is gone from this history");
        let fresh = Document::new_sequence("n", b"ACGT", false).unwrap();
        assert!(fresh.is_modified(), "a new, unsaved sequence is modified");
    }

    #[test]
    fn external_steps_sit_in_the_history_in_order_and_never_change_the_modified_flag() {
        let mut doc = sample();
        let before = doc.record.seq.clone();
        doc.external("Added primer p1", "t1");
        assert!(!doc.is_modified(), "a primer database change does not make the document modified");
        doc.replace(0, 1, b"G").unwrap();
        doc.external("Deleted primer p2", "t2");
        assert!(doc.is_modified());
        let (undo, redo) = doc.history();
        assert_eq!(undo.len(), 3);
        assert_eq!((undo[0].note, undo[0].inserted.as_str(), undo[0].removed.as_str()), ("external", "Added primer p1", "t1"));
        assert_eq!(undo[1].note, "");
        assert!(redo.is_empty());
        // Undo walks back through them in order; the app is told about each external one first.
        assert_eq!(doc.next_external(false).as_deref(), Some("t2"));
        assert!(doc.undo() && doc.is_modified());
        assert_eq!(doc.next_external(false), None, "the next step is the base edit");
        assert_eq!(doc.next_external(true).as_deref(), Some("t2"));
        assert!(doc.undo() && !doc.is_modified(), "undoing the base edit gives back the saved state, with an external step still in the history");
        assert_eq!(doc.next_external(false).as_deref(), Some("t1"));
        assert!(doc.undo() && !doc.is_modified() && doc.record.seq == before);
        assert!(!doc.undo(), "nothing left");
        assert!(doc.redo() && doc.redo() && doc.redo());
        assert!(doc.is_modified() && doc.redo() == false);
        // A new edit drops the redo steps, external ones included.
        assert!(doc.undo() && doc.undo());
        doc.replace(2, 2, b"T").unwrap();
        assert_eq!(doc.next_external(true), None);
        // Saving counts the base edits only: an external step before or after it changes nothing.
        let path = std::env::temp_dir().join("ge_test_external.gb");
        doc.save_genbank(&path, None).unwrap();
        doc.external("Edited primer p3", "t3");
        assert!(!doc.is_modified());
        assert!(doc.undo() && !doc.is_modified(), "undoing only the external step");
    }

    // ---------------------------------- section 2 -----------------------------------
    // A feature must keep covering the same bases. Each test builds a 200 bp random
    // sequence `s` and states which bases a feature must cover afterwards.

    fn drift_doc(features: Vec<Feature>) -> (Document, Vec<u8>) {
        let s = Lcg(99).bases(200);
        (doc_with(&s, false, features), s)
    }

    #[test]
    fn insert_at_the_exact_start_or_end_of_a_feature_does_not_grow_it() {
        let (mut doc, s) = drift_doc(vec![feat("gene", "51..100", "f")]);
        let before = doc.record.clone();
        let own = s[50..100].to_vec();
        doc.replace(50, 50, b"TTT").unwrap(); // exactly at the first base
        assert_eq!(spans_of(&doc, "f"), vec![(53, 103)], "insert at the start must shift the feature right");
        assert_eq!(covered(&doc, "f"), own);
        assert!(doc.undo());
        doc.replace(100, 100, b"TTT").unwrap(); // exactly after the last base
        assert_eq!(spans_of(&doc, "f"), vec![(50, 100)], "insert at the end must not grow the feature");
        assert_eq!(covered(&doc, "f"), own);
        assert!(doc.undo());
        assert_eq!(doc.record, before);
    }

    #[test]
    fn delete_overlapping_the_start_the_end_or_the_middle_keeps_the_remaining_bases() {
        let cases: [(usize, usize); 6] = [(40, 60), (90, 110), (70, 80), (50, 60), (90, 100), (20, 70)];
        for (a, b) in cases {
            let (mut doc, s) = drift_doc(vec![feat("gene", "51..100", "f")]);
            let before = doc.record.clone();
            doc.replace(a, b, b"").unwrap();
            let want: Vec<u8> = (50..100).filter(|i| *i < a || *i >= b).map(|i| s[i]).collect();
            assert_eq!(covered(&doc, "f"), want, "delete({a},{b}): wrong bases left in the feature");
            assert_eq!(spans_of(&doc, "f").len(), 1);
            assert!(doc.undo());
            assert_eq!(doc.record, before, "delete({a},{b}): undo");
        }
    }

    #[test]
    fn deleting_exactly_the_feature_or_more_removes_it_and_undo_restores_it_whole() {
        for (a, b) in [(50usize, 100usize), (49, 101), (0, 200), (30, 100), (50, 130)] {
            let mut f = feat("CDS", "complement(join(51..70,81..100))", "f");
            f.qualifiers.push(("note".into(), Some("keep me".to_string())));
            f.qualifiers.push(("codon_start".into(), Some("1".to_string())));
            let (mut doc, _) = drift_doc(vec![f]);
            let before = doc.record.clone();
            doc.replace(a, b, b"").unwrap();
            if a <= 50 && b >= 100 {
                assert!(!has(&doc, "f") && doc.record.features.is_empty(), "delete({a},{b}) must remove the feature");
            } else {
                assert_eq!(doc.record.features.len(), 1, "delete({a},{b}) must keep what is left of it");
            }
            assert!(doc.undo());
            assert_eq!(doc.record, before, "delete({a},{b}): undo must restore location and qualifiers");
            assert!(doc.redo() && doc.undo());
            assert_eq!(doc.record.features, before.features);
        }
    }

    #[test]
    fn replace_over_part_of_a_feature_follows_the_documented_rules() {
        // Documented rule (locations.rs): a replace is "delete the range, then insert at
        // its start", so the new bases join the feature only when the replaced stretch is
        // strictly inside it; bases replaced at its very start or end leave the feature.
        let x: &[u8] = b"TTTTTTTTTTTTT"; // 13 bases
        let (mut doc, s) = drift_doc(vec![feat("gene", "51..100", "f")]);
        // Strictly inside, longer replacement: feature grows by 3.
        doc.replace(60, 70, x).unwrap();
        assert_eq!(covered(&doc, "f"), [&s[50..60], x, &s[70..100]].concat());
        assert!(doc.undo());
        // Strictly inside, same length: still the same length, replaced bases covered.
        doc.replace(60, 70, &x[..10]).unwrap();
        assert_eq!(spans_of(&doc, "f"), vec![(50, 100)]);
        assert_eq!(covered(&doc, "f"), [&s[50..60], &x[..10], &s[70..100]].concat());
        assert!(doc.undo());
        // Shorter replacement inside: feature shrinks.
        doc.replace(60, 70, b"GG").unwrap();
        assert_eq!(covered(&doc, "f"), [&s[50..60], b"GG".as_slice(), &s[70..100]].concat());
        assert!(doc.undo());
        // Covering the start: the overlapped bases leave the feature, the new ones sit before it.
        doc.replace(40, 60, x).unwrap();
        assert_eq!(covered(&doc, "f"), s[60..100].to_vec());
        assert!(doc.undo());
        // Covering the end: likewise.
        doc.replace(90, 110, x).unwrap();
        assert_eq!(covered(&doc, "f"), s[50..90].to_vec());
        assert!(doc.undo());
        // Covering all of it: the feature is gone.
        doc.replace(50, 100, x).unwrap();
        assert!(doc.record.features.is_empty());
    }

    #[test]
    fn complement_feature_keeps_its_bases_and_strand_through_edits() {
        let (mut doc, s) = drift_doc(vec![feat("gene", "complement(51..100)", "f")]);
        let own = s[50..100].to_vec();
        doc.replace(10, 10, b"ACGTA").unwrap(); // before: the feature moves to 55..105
        doc.replace(190, 190, b"CC").unwrap(); // after
        assert_eq!(strand_of(&doc, "f"), -1);
        assert_eq!(covered(&doc, "f"), own, "edits outside moved the feature off its bases");
        doc.replace(70, 75, b"").unwrap(); // inside: original bases 65..70 of the feature
        assert_eq!(covered(&doc, "f"), [&own[..15], &own[20..]].concat());
        assert_eq!(bio(&doc, "f"), rc(&covered(&doc, "f")));
        assert_eq!(doc.record.features[0].location.to_gb_format(), "complement(56..100)");
    }

    #[test]
    fn join_edit_between_parts_grows_the_gap_inside_a_part_grows_that_part() {
        let (mut doc, s) = drift_doc(vec![feat("gene", "join(21..40,81..100)", "j")]);
        let own = [&s[20..40], &s[80..100]].concat();
        let x: &[u8] = b"ACGT";
        let cases: [(usize, Vec<(i64, i64)>); 8] = [
            (60, vec![(20, 40), (84, 104)]),  // in the gap: the gap grows
            (30, vec![(20, 44), (84, 104)]),  // inside part 1
            (90, vec![(20, 40), (80, 104)]),  // inside part 2
            (40, vec![(20, 40), (84, 104)]),  // at the end of part 1: gap side
            (80, vec![(20, 40), (84, 104)]),  // at the start of part 2: gap side
            (20, vec![(24, 44), (84, 104)]),  // at the start of part 1: shifts
            (100, vec![(20, 40), (80, 100)]), // at the end of part 2: nothing
            (0, vec![(24, 44), (84, 104)]),   // before all of it
        ];
        for (at, want) in cases {
            doc.replace(at, at, x).unwrap();
            assert_eq!(spans_of(&doc, "j"), want, "insert at {at}");
            let got = covered(&doc, "j");
            let size: i64 = want.iter().map(|(a, b)| b - a).sum();
            assert_eq!(got.len() as i64, size);
            if size as usize == own.len() {
                assert_eq!(got, own, "insert at {at}: bases changed");
            }
            assert!(doc.undo());
        }
        // A delete spanning the gap trims both parts.
        doc.replace(30, 90, b"").unwrap();
        assert_eq!(covered(&doc, "j"), [&s[20..30], &s[90..100]].concat());
    }

    #[test]
    fn complement_join_order_survives_edits_and_reverse_complement() {
        for text in ["complement(join(21..40,81..100))", "join(21..40,81..100)"] {
            let (mut doc, s) = drift_doc(vec![feat("CDS", text, "cj")]);
            let reverse = text.starts_with("complement");
            let before_bio = bio(&doc, "cj");
            assert_eq!(before_bio, if reverse { rc(&[&s[20..40], &s[80..100]].concat()) } else { [&s[20..40], &s[80..100]].concat() });
            doc.replace(60, 60, b"GGGG").unwrap(); // in the gap
            doc.replace(0, 5, b"").unwrap(); // before everything
            assert_eq!(bio(&doc, "cj"), before_bio, "{text}: edits outside the parts changed the feature");
            let parts = spans_of(&doc, "cj");
            assert_eq!(parts.len(), 2);
            assert!(parts[0].1 <= parts[1].0, "{text}: part order lost: {parts:?}");
            // Reverse-complement the whole sequence: same molecule, other strand.
            let old_strand = strand_of(&doc, "cj");
            let old_parts = parts.clone();
            let n = doc.len() as i64;
            doc.reverse_complement();
            assert_eq!(strand_of(&doc, "cj"), -old_strand, "{text}: strand not flipped");
            assert_eq!(bio(&doc, "cj"), before_bio, "{text}: the feature's own sequence changed");
            let mut mirrored: Vec<(i64, i64)> = old_parts.iter().map(|&(a, b)| (n - b, n - a)).collect();
            mirrored.sort();
            let mut now = spans_of(&doc, "cj");
            now.sort();
            assert_eq!(now, mirrored, "{text}: mirrored parts");
            doc.reverse_complement(); // twice = identity
            assert_eq!(spans_of(&doc, "cj"), old_parts);
        }
    }

    #[test]
    fn feature_crossing_the_origin_keeps_its_bases() {
        let n = 200;
        let s = Lcg(98).bases(n);
        let mk = || doc_with(&s, true, vec![feat("rep_origin", "join(191..200,1..10)", "o"), feat("gene", "101..150", "g")]);
        let own = [&s[190..], &s[..10]].concat();
        let x: &[u8] = b"ACGT";
        // Outside the feature, on either side: unchanged bases.
        for (a, b) in [(0usize, 0usize), (200, 200), (100, 110), (160, 170), (10, 10), (190, 190)] {
            let mut doc = mk();
            doc.replace(a, b, x).unwrap();
            assert_eq!(covered(&doc, "o"), own, "insert at {a}: origin feature drifted");
            assert!(doc.undo());
        }
        // Inside either part: that part grows.
        let mut doc = mk();
        doc.replace(195, 195, x).unwrap();
        assert_eq!(covered(&doc, "o"), [&s[190..195], x, &s[195..], &s[..10]].concat());
        assert!(doc.undo());
        doc.replace(5, 5, x).unwrap();
        assert_eq!(covered(&doc, "o"), [&s[190..], &s[..5], x, &s[5..10]].concat());
        assert!(doc.undo());
        // Deleting on both sides of the origin.
        doc.replace(0, 5, b"").unwrap();
        assert_eq!(covered(&doc, "o"), [&s[190..], &s[5..10]].concat());
        assert!(doc.undo());
        doc.replace(195, 200, b"").unwrap();
        assert_eq!(covered(&doc, "o"), [&s[190..195], &s[..10]].concat());
        assert!(doc.undo());
        // Rotating: any new start, including ones that cut through either part of the
        // feature, leaves the feature covering the same bases, in the same order.
        for k in [1usize, 5, 9, 10, 11, 60, 100, 150, 190, 191, 195, 199] {
            let mut doc = mk();
            let before = doc.record.clone();
            doc.rotate(k).unwrap();
            assert_eq!(bio(&doc, "o"), own, "rotate({k}): the origin feature changed its bases or their order");
            assert_eq!(bio(&doc, "g"), s[100..150].to_vec(), "rotate({k}): gene g drifted");
            assert!(doc.undo());
            assert_eq!(doc.record, before);
        }
    }

    #[test]
    fn edit_across_the_origin_keeps_features_on_both_sides() {
        let n = 200;
        let s = Lcg(97).bases(n);
        let mut doc = doc_with(
            &s,
            true,
            vec![feat("gene", "5..30", "a"), feat("gene", "100..150", "b"), feat("gene", "170..195", "c")],
        );
        let before = doc.record.clone();
        // Replace 190..205 (bases 191-200 and 1-5) with 3 new bases.
        doc.edit_across_origin(190, 205, b"TTT", None).unwrap();
        // Of a (0-based 4..30) the edit deletes only base 4; c (169..195) loses its last 5 bases.
        assert_eq!(covered(&doc, "a"), s[5..30].to_vec(), "a lost its first base and nothing else");
        assert_eq!(covered(&doc, "b"), s[99..150].to_vec());
        assert_eq!(covered(&doc, "c"), s[169..190].to_vec());
        assert!(doc.undo());
        assert_eq!(doc.record, before);
        // Replacing with pasted features: the pasted one covers the pasted bases.
        let pasted = vec![feat("misc_feature", "1..3", "p")];
        let at = doc.edit_across_origin(190, 205, b"TTT", Some(pasted)).unwrap();
        assert_eq!(covered(&doc, "p"), b"TTT".to_vec());
        assert_eq!(spans_of(&doc, "p"), vec![(at as i64, at as i64 + 3)]);
        assert_eq!(covered(&doc, "b"), s[99..150].to_vec());
    }

    #[test]
    fn fuzzy_locations_survive_edits_and_saves_with_their_markers() {
        let (mut doc, _) = drift_doc(vec![feat("misc_feature", "<21..>60", "z"), feat("CDS", "complement(<101..>150)", "y")]);
        assert_eq!(doc.record.features[0].location.to_gb_format(), "<21..>60");
        doc.replace(0, 0, b"TTT").unwrap();
        assert_eq!(doc.record.features[0].location.to_gb_format(), "<24..>63", "marker lost or shifted wrongly");
        doc.replace(40, 40, b"AA").unwrap(); // inside z
        assert_eq!(doc.record.features[0].location.to_gb_format(), "<24..>65");
        doc.replace(150, 160, b"").unwrap();
        assert_eq!(doc.record.features[1].location.to_gb_format(), "complement(<106..>150)");
        let path = std::env::temp_dir().join("ge_test_fuzzy.gb");
        doc.save_genbank(&path, None).unwrap();
        let back = Document::open(&path).unwrap();
        assert_eq!(back.record.features, doc.record.features, "fuzzy markers changed on save");
        assert!(doc.undo() && doc.undo() && doc.undo());
        assert_eq!(doc.record.features[0].location.to_gb_format(), "<21..>60");
    }

    #[test]
    fn reverse_complement_flips_every_feature_onto_the_same_bases() {
        let feats = vec![
            feat("gene", "21..60", "f1"),
            feat("gene", "complement(81..120)", "f2"),
            feat("CDS", "join(131..140,151..160)", "f3"),
            feat("CDS", "complement(join(171..180,191..200))", "f4"),
        ];
        for circular in [false, true] {
            let s = Lcg(96).bases(200);
            let mut doc = doc_with(&s, circular, feats.clone());
            let before = doc.record.clone();
            let old: Vec<(Vec<u8>, i8, Vec<(i64, i64)>)> =
                ["f1", "f2", "f3", "f4"].iter().map(|l| (bio(&doc, l), strand_of(&doc, l), spans_of(&doc, l))).collect();
            doc.reverse_complement();
            assert_eq!(doc.record.seq, rc(&s));
            for (l, (b, strand, spans)) in ["f1", "f2", "f3", "f4"].iter().zip(old) {
                assert_eq!(bio(&doc, l), b, "{l}: the feature's own sequence changed");
                assert_eq!(strand_of(&doc, l), -strand, "{l}: strand not flipped");
                let mut want: Vec<(i64, i64)> = spans.iter().map(|&(a, e)| (200 - e, 200 - a)).collect();
                want.sort();
                let mut got = spans_of(&doc, l);
                got.sort();
                assert_eq!(got, want, "{l}: wrong bases covered");
                // Qualifiers untouched.
                assert_eq!(doc.record.features[index_of(&doc, l)].qualifiers.len(), 1);
            }
            assert!(doc.undo());
            assert_eq!(doc.record, before);
        }
    }

    #[test]
    fn copying_the_reverse_complement_keeps_every_feature_on_the_same_bases() {
        // 200 random bases; a copy of 41..160 (0-based 40..160) is 120 bases long.
        let feats = vec![
            feat("gene", "61..100", "top"),
            feat("gene", "complement(101..140)", "bottom"),
            feat("CDS", "join(51..60,111..120)", "split"),
        ];
        let s = Lcg(11).bases(200);
        let doc = doc_with(&s, false, feats);
        let copy = doc.features_between(40, 160);
        let flipped = doc.features_between_reversed(40, 160);
        assert_eq!(copy.len(), flipped.len());
        // Build the pasted-in-reverse molecule: reverse-complemented bases + flipped features.
        let rc_bases = rc(&s[40..160]);
        let target = doc_with(&rc_bases, false, flipped);
        for l in ["top", "bottom", "split"] {
            // The feature's own sequence, read along its strand, is unchanged by the flip...
            let before = doc_with(&s[40..160], false, copy.clone());
            assert_eq!(bio(&target, l), bio(&before, l), "{l}: its sequence changed");
            // ...and it is on the opposite strand.
            assert_eq!(strand_of(&target, l), -strand_of(&before, l), "{l}: strand not flipped");
        }
    }

    #[test]
    fn pasted_features_cover_the_pasted_bases_and_existing_ones_shift_correctly() {
        let (mut doc, s) = drift_doc(vec![feat("gene", "51..100", "g"), feat("gene", "121..150", "h"), feat("gene", "complement(join(11..20,31..40))", "k")]);
        let copied = doc.features_between(60, 130); // g (60..100) and h (120..130), relative
        assert_eq!(copied.len(), 2);
        let bases = s[60..130].to_vec();
        doc.paste(0, 0, &bases, copied).unwrap(); // paste at the very start
        // Old features moved right by 70, keeping their bases.
        let g_old = doc.record.features.iter().position(|f| f.qualifier_values("label").next() == Some("g")).unwrap();
        assert_eq!(locations::spans(&doc.record.features[g_old].location).0, vec![(120, 170)]);
        assert_eq!(doc.record.seq[120..170].to_vec(), s[50..100].to_vec());
        // The two pasted copies (appended at the end) sit on the pasted bases at 0..70.
        let n = doc.record.features.len();
        let (g2, h2) = (&doc.record.features[n - 2], &doc.record.features[n - 1]);
        assert_eq!(locations::spans(&g2.location).0, vec![(0, 40)]);
        assert_eq!(doc.record.seq[0..40].to_vec(), s[60..100].to_vec());
        assert_eq!(locations::spans(&h2.location).0, vec![(60, 70)]);
        assert_eq!(doc.record.seq[60..70].to_vec(), s[120..130].to_vec());
    }

    #[test]
    fn qualifiers_are_untouched_by_edits() {
        let mut f = feat("CDS", "join(21..40,81..100)", "f");
        f.qualifiers.push(("note".into(), Some("a long note with words".to_string())));
        f.qualifiers.push(("translation".into(), Some("MKVLAAGIVALLL".to_string())));
        f.qualifiers.push(("pseudo".into(), None));
        let (mut doc, _) = drift_doc(vec![f]);
        let q = doc.record.features[0].qualifiers.clone();
        let mut rng = Lcg(5);
        for _ in 0..60 {
            let n = doc.len();
            let a = rng.below(n + 1);
            let b = (a + rng.below(6)).min(n);
            let ins = rng.some_bases(5);
            doc.replace(a, b, &ins).unwrap();
            if doc.record.features.is_empty() {
                break;
            }
            assert_eq!(doc.record.features[0].qualifiers, q, "an edit changed the qualifiers");
        }
    }

    #[test]
    fn making_linear_is_one_undo_step_that_rotates_and_undo_gives_the_circle_back() {
        let seq = Lcg(3).bases(100);
        let mut doc = doc_with(
            &seq,
            true,
            vec![feat("gene", "21..40", "g"), feat("rep_origin", "join(91..100,1..10)", "o"), feat("misc_feature", "31..70", "cut")],
        );
        let before = doc.record.clone();
        let o = [&seq[90..], &seq[..10]].concat();
        doc.linearize(50).unwrap(); // break before base 51: cuts the feature "cut" (bases 31..70)
        assert_eq!(doc.record.topology, gb_io::seq::Topology::Linear);
        assert_eq!(doc.record.seq, [&seq[50..], &seq[..50]].concat());
        assert_eq!(bio(&doc, "g"), seq[20..40].to_vec());
        assert_eq!(bio(&doc, "o"), o, "the origin feature changed its bases");
        assert_eq!(doc.history().0.len(), 1, "make linear must be ONE undo step");
        assert!(doc.undo() && !doc.undo());
        assert_eq!(doc.record, before, "undo must give back the circle with its old origin");
        assert!(doc.redo());
        assert_eq!(doc.record.topology, gb_io::seq::Topology::Linear);
        assert!(doc.undo());
        // Keep the current start; already linear; and the toggle is a step too.
        doc.linearize(0).unwrap();
        assert_eq!(doc.record.seq, seq);
        assert!(doc.linearize(5).is_err(), "already linear");
        assert!(doc.undo());
        assert_eq!(doc.record.topology, gb_io::seq::Topology::Circular);
        doc.set_topology(false);
        doc.set_topology(false); // no change, no step
        assert_eq!(doc.history().0.len(), 1);
        assert!(doc.undo());
        assert_eq!(doc.record.topology, gb_io::seq::Topology::Circular);
        assert!(!doc.is_modified());
    }

    /// The yellow "edited bases" list, checked against a simple model: a Vec of one flag per
    /// base that is edited with the very same operations (a plain, obviously-right version).
    /// Random edits of every kind, with undo / redo, must keep the two in agreement.
    #[test]
    fn edited_ranges_follow_every_edit_undo_and_redo() {
        // Turn the flag-per-base model into the [start, end) list the document reports.
        fn runs(flags: &[bool]) -> Vec<(usize, usize)> {
            let mut out: Vec<(usize, usize)> = Vec::new();
            for (i, &f) in flags.iter().enumerate() {
                if !f {
                    continue;
                }
                match out.last_mut() {
                    Some(l) if l.1 == i => l.1 = i + 1,
                    _ => out.push((i, i + 1)),
                }
            }
            out
        }
        for circular in [false, true] {
            let mut rng = Lcg(77 + circular as u64);
            let mut doc = doc_with(&rng.bases(40), circular, vec![]);
            // One flag list per undo depth, so undo / redo can be compared with the model.
            let mut states: Vec<Vec<bool>> = vec![vec![false; 40]];
            let mut pos = 0usize;
            for step in 0..300 {
                let n = doc.len();
                let mut flags = states[pos].clone();
                match rng.below(6) {
                    0 | 1 | 2 => {
                        // replace / insert / delete in the middle
                        let a = rng.below(n + 1);
                        let b = (a + rng.below(6)).min(n);
                        let ins = rng.below(5);
                        if a == b && ins == 0 {
                            continue;
                        }
                        doc.replace(a, b, &rng.bases(ins)).unwrap();
                        flags.splice(a..b, vec![true; ins]);
                    }
                    3 => {
                        let a = rng.below(n);
                        doc.overwrite(a, b"A").unwrap();
                        flags[a] = true;
                    }
                    4 if circular => {
                        let k = 1 + rng.below(n - 1);
                        doc.rotate(k).unwrap();
                        flags.rotate_left(k);
                    }
                    4 => {
                        doc.reverse_complement();
                        flags.reverse();
                    }
                    _ => {
                        // go back one step (and sometimes forward again)
                        if pos > 0 && doc.undo() {
                            pos -= 1;
                            assert_eq!(doc.edited_ranges(), runs(&states[pos]).as_slice(), "step {step}: after undo");
                            if rng.below(2) == 0 {
                                assert!(doc.redo());
                                pos += 1;
                                assert_eq!(doc.edited_ranges(), runs(&states[pos]).as_slice(), "step {step}: after redo");
                            }
                        }
                        continue;
                    }
                }
                states.truncate(pos + 1);
                states.push(flags);
                pos += 1;
                assert_eq!(doc.edited_ranges(), runs(&states[pos]).as_slice(), "circular={circular} step {step}");
            }
        }
    }

    /// Typing a run joins into one stretch, a deletion leaves nothing, an edit across the
    /// origin marks the new bases on both sides of it, and a save empties the list.
    #[test]
    fn edited_ranges_basics_origin_and_save() {
        let mut doc = doc_with(b"ACGTACGTACGTACGTACGT", false, vec![]);
        assert!(doc.edited_ranges().is_empty(), "nothing is edited when opened");
        doc.replace(5, 5, b"G").unwrap();
        doc.replace(6, 6, b"G").unwrap();
        doc.replace(7, 7, b"G").unwrap();
        assert_eq!(doc.edited_ranges(), &[(5, 8)], "typed run");
        doc.replace(0, 2, b"").unwrap(); // delete before: the run moves left
        assert_eq!(doc.edited_ranges(), &[(3, 6)]);
        doc.replace(4, 5, b"").unwrap(); // delete inside: the run shrinks
        assert_eq!(doc.edited_ranges(), &[(3, 5)]);
        assert!(doc.undo() && doc.undo());
        assert_eq!(doc.edited_ranges(), &[(5, 8)]);
        assert!(doc.undo() && doc.undo() && doc.undo());
        assert!(doc.edited_ranges().is_empty(), "undoing an insertion removes its highlight");
        doc.replace(3, 4, b"TT").unwrap();
        let path = std::env::temp_dir().join("opendna_edited_ranges_test.gb");
        doc.save_genbank(&path, None).unwrap();
        let _ = std::fs::remove_file(&path);
        assert!(doc.edited_ranges().is_empty(), "save clears");
        assert!(doc.undo() && doc.edited_ranges().is_empty(), "undo after save marks nothing");

        // Circular, across the origin: replace bases 18..22 (18, 19, 0, 1) with "GGG".
        let mut c = doc_with(b"ACGTACGTACGTACGTACGT", true, vec![]);
        c.edit_across_origin(18, 22, b"GGG", None).unwrap();
        assert_eq!(c.edited_ranges().len(), 1);
        let n = c.len();
        let (a, b) = c.edited_ranges()[0];
        assert_eq!(b - a, 3, "the three new bases are marked");
        assert_eq!(n, 19);
        assert!(c.undo() && c.edited_ranges().is_empty());
        assert!(c.redo() && c.edited_ranges().len() == 1);
    }
}
