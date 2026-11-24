"""Shared state for one open document: the Rust Document, the selection, and cached
feature layout. Views never talk to each other: they call Editor methods and repaint
on its signals.

Positions are 0-based, end-exclusive. The selection is [start, end); when it is empty
it is just a cursor sitting before base `start`.
"""
from collections import Counter
import weakref
import json
from bisect import bisect_left, insort
from pathlib import Path
from types import SimpleNamespace

from PySide6.QtCore import QEvent, QMimeData, QObject, QRegularExpression, QThreadPool, Signal
from PySide6.QtGui import QColor, QFont, QFontDatabase, QKeySequence, QRegularExpressionValidator
from PySide6.QtWidgets import QApplication, QToolTip

from opendna_core import find_primers

import multigb

import genetic_code
import inventory
from enzymes import COMMON, MAX_DRAWN_CUTS, EnzymeFilter, find_sites, is_owned, wanted
from genes import GENE, redundant_genes
from perf import timed
from uikit import SYNC_LIMIT, Latest
from units import is_tu
import tool_settings
from primer_db import BINDING, PCR, binding, db, notifier, tm_rule

PRIMER_KIND = "primer_bind"
# `source` spans the whole record (.dna files hide it too). `primer_bind` features stored in
# files are not drawn: primers come from the primer database and their sites are
# recomputed on the current sequence (see Editor.match_primers).
HIDDEN_KINDS = {"source", PRIMER_KIND}
# A GenBank COMMENT holding the primers and sites hidden in the file (other programs show it as a line of text, nothing more)
DISPLAY_NOTE = "openDNA-display:"
CLIP_MIME = "application/x-opendna-clip"  # copied bases + features, JSON
# Qualifiers searched by Find / the genome view's locus box.
NAME_KEYS = {"label", "gene", "locus_tag", "old_locus_tag", "product", "gene_synonym", "protein_id", "note"}
MAX_SITES = 50               # per primer and sequence DRAWN and listed; more means it is not specific anyway (the primers panel and the status bar say so)
TRACK_HITS = 20_000          # track items a search stops at (the search bar's counter says so)
# Standard colour of each annotation type: drawn when a feature has no colour of its own,
# given to new features, and applied by Feature → Recolor Annotations by Type.
KIND_COLORS = {
    # coding
    "CDS": "#e8a33d", "gene": "#9ccc65", "mRNA": "#c5e1a5", "exon": "#aed581", "intron": "#e0e0e0",
    "sig_peptide": "#ffcc80", "mat_peptide": "#ffb74d", "transit_peptide": "#ffe0b2", "5'UTR": "#dcedc8",
    "3'UTR": "#dcedc8",
    # expression control
    "promoter": "#4fc3f7", "-35_signal": "#81d4fa", "-10_signal": "#29b6f6", "RBS": "#4db6ac",
    "regulatory": "#4dd0e1", "enhancer": "#80deea", "protein_bind": "#26a69a", "terminator": "#e57373",
    "polyA_signal": "#ff8a65",
    # replication and structure
    "rep_origin": "#ffd54f", "oriT": "#fff176", "LTR": "#a1887f", "repeat_region": "#bcaaa4",
    "stem_loop": "#d7ccc8", "mobile_element": "#ba68c8", "misc_recomb": "#ce93d8", "misc_binding": "#80cbc4",
    "centromere": "#b39ddb", "telomere": "#b39ddb",
    # RNA
    "tRNA": "#f06292", "rRNA": "#f06292", "ncRNA": "#f48fb1", "misc_RNA": "#f48fb1",
    # others
    "primer_bind": "#9575cd", "misc_feature": "#90a4ae",
}
DEFAULT_COLOR = "#90a4ae"
COLOR_KEYS = ("ApEinfo_fwdcolor", "ApEinfo_revcolor")  # what .dna files and ApE read and write


# Direction: a feature has an arrow (its strand) unless its type has none by default or it says otherwise. GenBank cannot
# say "no direction" in a location, so a /directional="no" (or "yes") qualifier overrides the type's default; the .dna file's
# own flag is read into it (core/src/dna_file.rs).
DIRECTIONLESS_KINDS = {"misc_feature"}
DIRECTION_KEY = "directional"


def is_directional(kind, qualifiers):
    value = next((v for k, v in qualifiers if k == DIRECTION_KEY), None)
    if value is not None and value.strip().lower() in ("yes", "no"):
        return value.strip().lower() == "yes"
    return kind not in DIRECTIONLESS_KINDS


def arrow(f):
    """+1 / -1: the arrow of feature `f` (its strand); 0: drawn without one."""
    return f.strand if is_directional(f.kind, f.qualifiers) else 0


def type_color(kind):
    """Standard colour of an annotation type (grey for unknown types)."""
    return QColor(KIND_COLORS.get(kind, DEFAULT_COLOR))
IUPAC = set("ACGTURYSWKMBDHVN-acgturyswkmbdhvn")
TYPABLE = set("ACGTNacgtn")  # keyboard input; ambiguity codes only arrive by paste
_COMPLEMENT = str.maketrans("ACGTURYSWKMBDHVNacgturyswkmbdhvn", "TGCAAYRSWMKVHDBNtgcaayrswmkvhdbn")


def complement(seq):
    return seq.translate(_COMPLEMENT)


def reverse_complement(seq):
    return seq.translate(_COMPLEMENT)[::-1]


def gc_percent(seq):
    return 100 * sum(c in "GC" for c in seq.upper()) / max(len(seq), 1)


def mono_font(size):
    """A clear fixed-width font for bases (Windows' default Courier New is thin)."""
    families = QFontDatabase.families()
    for family in ("Cascadia Mono", "Consolas", "DejaVu Sans Mono", "Menlo"):
        if family in families:
            return QFont(family, size)
    font = QFontDatabase.systemFont(QFontDatabase.FixedFont)
    font.setPointSize(size)
    return font


def label_font(base):
    """Feature labels drawn on colored boxes: the UI font, a bit bolder."""
    font = QFont(base)
    font.setPointSizeF(9)
    font.setWeight(QFont.DemiBold)
    return font


def luminance(color):
    """WCAG relative luminance of a colour: 0 (black) to 1 (white), on linear light, so that green weighs most and blue least."""
    def linear(v):
        v /= 255
        return v / 12.92 if v <= 0.04045 else ((v + 0.055) / 1.055) ** 2.4
    return 0.2126 * linear(color.red()) + 0.7152 * linear(color.green()) + 0.0722 * linear(color.blue())


def contrast(a, b):
    """WCAG contrast ratio of two colours, 1 (none) to 21 (black on white); 4.5 is the minimum for readable text."""
    hi, lo = sorted((luminance(a), luminance(b)), reverse=True)
    return (hi + 0.05) / (lo + 0.05)


def text_color_on(background):
    """Black or white, whichever has the higher contrast on `background`. Whatever the colour, that is at least 4.58 to 1 (the worst
    colour, a grey of luminance 0.18, gives it with either): a mid-tone that a brightness threshold sent the wrong way is read right."""
    black, white = QColor("black"), QColor("white")
    return black if contrast(black, background) >= contrast(white, background) else white


def feature_color(f):
    for key, value in f.qualifiers:
        if key == "ApEinfo_fwdcolor" and value:
            c = QColor(value)
            if c.isValid() and c.lightness() < 235:  # the .dna default white is invisible
                return c
    return type_color(f.kind)


def feature_bounds(f):
    """Selectable range of a feature. Origin-spanning features select their first part."""
    starts = [s for s, _ in f.spans]
    if starts != sorted(starts):  # wraps around the origin
        return f.spans[0]
    return f.spans[0][0], max(e for _, e in f.spans)


def primer_sequence(f):
    """Full 5'->3' sequence of a primer_bind feature from its /note="sequence: ...", or
    None (then only the annealed bases are known)."""
    for key, value in f.qualifiers:
        if key == "note" and value and value.startswith("sequence: "):
            return value[len("sequence: "):]
    return None


def site_location(start, end, strand, length):
    """GenBank location text for a binding site (end may pass the origin)."""
    loc = f"{start + 1}..{end}" if end <= length else f"join({start + 1}..{length},1..{end - length})"
    return f"complement({loc})" if strand < 0 else loc


def primer_extent(f, tail, length=None):
    """Spans of a primer including its 5' tail (which hangs off the template). `length`: of a circular sequence,
    whose tail goes on past the origin to the other end (a linear one is cut at its start)."""
    spans = list(f.spans)
    if f.strand > 0:  # 5' end on the left of the first part
        s, e = spans[0]
        spans[0] = (max(0, s - tail), e)
        if length and s - tail < 0:
            spans.append((max(0, length + s - tail), length))
    else:             # 5' end on the right of the last part
        s, e = spans[-1]
        spans[-1] = (s, e + tail)
        if length and e + tail > length:
            spans[-1] = (s, length)
            spans.append((0, min(e + tail - length, length)))
    return spans


def pack_lanes(features):
    """Greedy packing: each feature gets the first lane where none of its spans overlap
    another feature. Returns ({feature index: lane}, number of lanes)."""
    lanes = []  # per lane: sorted, non-overlapping (start, end) intervals
    result = {}

    def free(occupied, s, e):
        j = bisect_left(occupied, (s, s))
        return (j == 0 or occupied[j - 1][1] <= s) and (j == len(occupied) or occupied[j][0] >= e)

    for f in sorted(features, key=lambda f: min(s for s, _ in f.spans)):
        spans = [(s, max(e, s + 1)) for s, e in f.spans]
        lane = next((i for i, occ in enumerate(lanes) if all(free(occ, s, e) for s, e in spans)), None)
        if lane is None:
            lanes.append([])
            lane = len(lanes) - 1
        for span in spans:
            insort(lanes[lane], span)
        result[f.index] = lane
    return result, len(lanes)


class _PasteCleaner(QObject):
    """Ctrl+V (Shift+Insert) in a sequence field: the clipboard goes through `clean`, so a pasted FASTA record or a GenBank block
    (header line, spaces, line breaks, numbers) lands as bare letters; text that is not a sequence is refused with a tooltip
    that says why. (The validator alone would have refused the whole paste.) The right-click menu's Paste still goes
    through the validator only."""

    def __init__(self, field, clean):
        super().__init__(field)
        self.field, self.clean = field, clean
        field.installEventFilter(self)

    def eventFilter(self, obj, event):
        if event.type() == QEvent.KeyPress and event.matches(QKeySequence.Paste):
            try:
                self.field.insert(self.clean(QApplication.clipboard().text()))
            except ValueError as err:
                QToolTip.showText(self.field.mapToGlobal(self.field.rect().bottomLeft()), str(err), self.field)
            return True
        return False


def dna_field(field):
    """Let a text field take only DNA letters (IUPAC codes, either case, spaces): other
    keys do nothing, and a paste with anything else is refused (a pasted FASTA record is cleaned first)."""
    field.setValidator(QRegularExpressionValidator(QRegularExpression(r"[ACGTURYSWKMBDHVNacgturyswkmbdhvn\s]*"), field))
    _PasteCleaner(field, clean_pasted)
    return field


def clean_pasted_protein(text):
    """Amino acids from pasted text: drops a FASTA header line, spaces, line breaks and digits. Raises ValueError for other characters."""
    if text.startswith(">"):
        text = text.partition("\n")[2]
    letters = "".join(c for c in text if not c.isspace() and not c.isdigit())
    bad = sorted({c for c in letters if c.upper() not in "ACDEFGHIKLMNPQRSTVWYX*"})
    if bad:
        raise ValueError(f"Not a protein, nothing pasted: {' '.join(bad[:10])}")
    return letters


def protein_field(field):
    """Same for amino acids (the 20, X for any, * for stop)."""
    field.setValidator(QRegularExpressionValidator(QRegularExpression(r"[ACDEFGHIKLMNPQRSTVWYXacdefghiklmnpqrstvwyx*\s]*"),
                                                   field))
    _PasteCleaner(field, clean_pasted_protein)
    return field


def clean_pasted(text):
    """Bases from pasted text: drops a FASTA header line, spaces, line breaks and digits
    (GenBank ORIGIN blocks). Raises ValueError if anything else is left."""
    if text.startswith(">"):
        text = text.partition("\n")[2]
    bases = "".join(c for c in text if not c.isspace() and not c.isdigit())
    bad = sorted({c for c in bases if c not in IUPAC})
    if bad:
        raise ValueError(f"Not DNA, nothing pasted: {' '.join(bad[:10])}")
    return bases


class _Relay(QObject):
    """Carries the primer sites found on the other records of a file to the GUI thread. Lives forever, so a
    late result never hits a deleted Editor."""
    elsewhere_ready = Signal(object, object)      # token, binding sites on the other records of a multi-record file


_relay = _Relay()


def feature_key(f):
    """What names a feature in `Editor.hidden_features` (its index changes when features are added or removed)."""
    return (f.kind, f.label, f.location)


def site_key(site):
    """What names a primer site among the sites of its primer, whatever else binds: its strand and its 3′ end."""
    return site.strand, site.end if site.strand > 0 else site.start


def site_rank(site, n):
    """Sort key of a site, the best first: fewest mismatches in the last 5 bases at the 3′ end (they decide whether the polymerase
    starts), then fewest mismatches, then the longest annealed part. `n`: the length of the sequence (a site may cross the origin)."""
    end3 = (site.end - 1) % n if site.strand > 0 else site.start
    near = sum(((end3 - p) if site.strand > 0 else (p - end3)) % n < 5 for p in site.mismatch_at)
    return near, site.mismatches, site.start - site.end


def _connect_weakly(signal, method):
    """Connect a bound method to a signal without keeping its object alive: once the object is gone the link does nothing (it
    is a closure of a few bytes, left connected: disconnecting from inside the signal only makes PySide warn)."""
    ref = weakref.WeakMethod(method)

    def call(*args):
        target = ref()
        if target is not None:
            target(*args)

    signal.connect(call)


class Editor(QObject):
    changed = Signal()            # sequence or features changed (also after save: title/modified)
    selection_changed = Signal()
    display_changed = Signal()    # a View option was toggled (enzymes, translation)
    search_changed = Signal()     # the places the search bar found (highlighted in the views) changed
    reveal = Signal(int)          # ask views to scroll this position into view
    navigate = Signal(int)        # a view was used to go to this position: the others bring it to their middle
    error = Signal(str)
    edit_feature_requested = Signal(int)  # a view was double-clicked on this feature
    edit_primer_requested = Signal(str)   # ... or on a primer site (primer name)
    primers_changed = Signal()    # primer binding sites were recomputed
    modified_changed = Signal()   # the file has something new to save that is not an edit (a primer or site hidden): the title
    primer_clicked = Signal(str)  # a primer site was clicked (or shift-clicked) in a view (its site key)
    elsewhere_changed = Signal()      # the primer sites on the other records of the file were (re)computed
    record_hit_requested = Signal(str, int, int)   # open this record of the file and select these bases
    record_requested = Signal(str)   # the map's chromosome list was clicked: open this sequence of the FASTA file
    records_order_requested = Signal(object)   # the map's chromosome list was reordered (drag and drop): [names] in the new order
    record_remove_requested = Signal(str)      # the map's chromosome list: remove this record from the file
    record_export_requested = Signal(str)      # the map's chromosome list: export this record as a GenBank file
    query_hit_changed = Signal(int)   # the Local Align hit (row) picked in the Alignment view or the dialog's table: its index, -1 none

    # Last Copy in this app, shared by all open documents: (bases, features relative to
    # the copied bases). Paste uses it when the clipboard still holds those bases.
    clip = None

    def __init__(self, doc, path):
        super().__init__()
        self.doc = doc
        self.path = Path(path)
        self.untitled = False        # made in memory: `path` is only a suggested name
        self.fasta_record = None     # sequence name when `path` is a FASTA file
        self.elsewhere, self._elsewhere_key, self._elsewhere_token = [], None, 0   # primer sites on the other records: [SimpleNamespace(record, label, sequence, start, end, strand, mismatches)]
        self.multi_genbank = False   # `fasta_record` names a record of a multi-record GenBank file (Save puts it back in the file)
        self.fasta_records = []      # [(name, length)] of every sequence of that file (the map lists them, one line each)
        self.start = self.end = 0
        self.current_feature = None  # index of the feature last selected as a whole
        self.version = 0             # bumped on every change (autosave compares it)
        self.externals, self._external_count = {}, 0   # token -> (undo, redo) of the history's external steps (record_external)
        self.show_enzymes = False    # restriction sites on map + sequence (View → Restriction Enzymes)
        self.enzyme_filter = EnzymeFilter()   # which of them (Display filter): the unique cutters until the criteria are changed
        self.site_cuts, self._shown_sites = {}, {}   # enzyme name -> its cut positions; filter -> [(position, name)] (cache)
        self.show_translation = False  # 3-frame translation under the bases
        self.show_complement = True    # the second strand under the top one (off: the top strand alone, one line of bases)
        self.cds_only = True         # a gene at exactly the bases of a CDS is not drawn: its CDS says it all (every view)
        self.redundant_genes = set() # indices of those gene features
        # Display filter (side panel): feature types not drawn, primers on/off, and the
        # amino acids drawn inside CDS arrows (all CDSs, minus those switched off one by one).
        self.hidden_kinds = set()
        self.hidden_features = set()    # feature_key of the features hidden one by one (right-click > Hide)
        self.show_file_primers = True   # Display filter > Primers: the file's own (and designed, unsaved) primers
        self.show_db_primers = False    # ... and the primers of the database that bind here (hidden until switched on)
        self.show_units = True          # ... transcription units (units.py), whatever their feature type
        self.temporary_primers = set()   # {(name, upper-case sequence)} of primers matched here that are not in the database yet
        # Primers that belong to THIS file: {(name, upper-case sequence): (name, sequence)}. They are the ones Save writes as
        # primer_bind features; a primer that is only in the database is drawn as a temporary one until it is imported.
        self.file_primer_set = self._primers_stored_in_file()
        self.search_strand = {}   # (start, end) -> 1 / -1 / 0: the strand a search hit was found on (0: both, or not a strand)
        self.search_spans, self._search_reach = [], 0   # [(start, end)] the search bar found, sorted (see set_search_hits)
        self.cut_marks = {}          # owner ->[SimpleNamespace(top, bottom, label, tip)]: cuts the views draw (digest, Cas9...)
        self.show_cds_translation = True
        self.untranslated = set()    # (location, label) of CDSs whose translation is off
        self.primer_sites, self.site_by_key, self.primer_lanes = [], {}, {}
        # Primers of the database not drawn on this file (still in the database, still listed greyed in the Primers view).
        self.hidden_primers = set()   # (for this session only; an .odp project keeps them)
        # Sites shown or hidden by hand: {primer name: {"sequence": upper-case, "sites": {site_key: shown}}}, under the binding
        # rule (seed, exact 3′, mismatches) `site_rule` they were made with; another rule finds other sites, so it starts afresh.
        self.site_rule, self.site_choices = None, {}
        self._read_display()
        self.current_primer = None   # key of the selected primer site
        self.tracks = []             # loaded annotation tracks (SimpleNamespace name, track, seqid)
        self.alignments = []         # loaded BAM files (SimpleNamespace name, bam, seqid, length, ...)
        self.query_hits = []         # Local Align hits (SimpleNamespace), drawn by the genome view, one row each
        self.current_query_hit = -1  # index in query_hits of the picked one
        self.dna_locked = False      # lock button of the toolbar: bases cannot change (features etc. still can)
        self._primer_job = Latest(self, "Primer sites")
        self._sites_job = Latest(self, "Restriction sites")
        # Weak links: a global signal must not keep a closed window's editor alive (each one would go on re-matching the primers).
        _connect_weakly(_relay.elsewhere_ready, self._take_elsewhere)
        _connect_weakly(notifier.changed, self.match_primers)
        _connect_weakly(tool_settings.notifier.changed, self._settings_changed)
        _connect_weakly(inventory.notifier.changed, self._inventory_changed)
        self._refresh()

    # ---- cached layout -------------------------------------------------------------
    def _refresh(self, keep_lanes=False):
        # shortcut: rebuilds the whole feature cache on every edit (~0.1 s on a 10k-feature
        # genome, instant on plasmids); make incremental if typing in genomes feels slow.
        self.length = self.doc.length
        self.query_hits, self.current_query_hit = [], -1  # any edit makes Local Align positions stale
        everything = self.doc.features()
        self.features = [f for f in everything if f.kind not in HIDDEN_KINDS]
        self.file_primers = [f for f in everything if f.kind == PRIMER_KIND]  # as stored in the file
        self.by_index = {f.index: f for f in self.features}
        # Transcription units and the count of every other type, once per edit (the Display panel and the filter asked each
        # feature again on every edit and display switch: tens of thousands of qualifier lists read from the core on a genome)
        self.unit_indexes = {f.index for f in self.features if is_tu(f)}
        self.kind_counts = Counter(f.kind for f in self.features if f.index not in self.unit_indexes)
        self.redundant_genes = redundant_genes(self.features)
        self._lay_out(keep_lanes)
        self._search_text = None
        self.match_primers()
        self.find_cut_sites()

    def find_cut_sites(self):
        """Where the enzymes cut (`site_cuts`), while they are shown. Small sequences: now. Large ones: in a background thread, and
        none drawn until they are found (never the sites of the sequence as it was)."""
        self.site_cuts, self._shown_sites = {}, {}
        if not self.show_enzymes:
            self._sites_job.cancel()
            return
        seq, circular, now = self.doc.sequence(), self.doc.circular, self.length <= SYNC_LIMIT
        self._sites_job.run(lambda: find_sites(seq, circular), lambda sites: self._take_cut_sites(sites, now), now)

    def _take_cut_sites(self, sites, now):
        self.site_cuts, self._shown_sites = sites, {}
        if not now:   # (drawn after the edit that asked for them: draw again)
            self.display_changed.emit()

    def _lay_out(self, keep_lanes=False):
        """What the views draw and in which lanes (no primer matching: set_overlay calls this alone)."""
        # What the views draw: all features but the types filtered out in the side panel, then the tools' temporary
        # ones (set_overlay) under the same filter.
        self.shown = [f for f in self.features if self.visible(f)] + [o for o in self.overlay_features if self.visible(o)]
        # Typing/deleting bases never makes two features overlap that did not before, so
        # the packing stays valid while the features are the same ones (same count).
        if not (keep_lanes and len(self.lanes) == len(self.shown)):
            self.lanes, self.lane_count = pack_lanes(self.shown)
            self.colors = {f.index: feature_color(f) for f in self.features + self.overlay_features}

    # ---- primers (from the database, matched on the current sequence) -----------------
    extra = ()   # [(name, sequence)] matched in addition to the database's (see set_extra)

    def _primers_stored_in_file(self):
        """{(name, SEQUENCE): (name, sequence)} of the primers the opened file stores: its primer_bind features (the full
        sequence from the note, else the annealed bases) and, for a .dna file, its primer list."""
        out = {}
        for name, seq in self.doc.primers():
            out[(name, seq.upper())] = (name, seq)
        for f in self.doc.features():
            if f.kind == PRIMER_KIND:
                out.setdefault(self._stored_key(f), (f.label, self._stored_sequence(f)))
        return out

    def _stored_sequence(self, f):
        """A stored primer_bind's primer: the full sequence from its note, else the bases it covers (its strand's)."""
        seq = primer_sequence(f)
        if seq is None:
            bases = "".join(self.doc.sequence(a, b) for a, b in f.spans)
            seq = reverse_complement(bases) if f.strand < 0 else bases
        return seq

    def _stored_key(self, f):
        return f.label, self._stored_sequence(f).upper()

    def in_file(self, name, sequence):
        """Is the primer (name, sequence) one of this file's own?"""
        return (name, sequence.upper()) in self.file_primer_set

    def set_in_file(self, pairs, present):
        """Add / remove primers [(name, sequence)] among the file's own, without a history step (the caller records one)."""
        for n, q in pairs:
            if present:
                self.file_primer_set[(n, q.upper())] = (n, q)
            else:
                self.file_primer_set.pop((n, q.upper()), None)
        self.match_primers()

    def rename_in_file(self, before, after):
        """A database primer was edited (row dicts): the file's own copy follows it."""
        if self.in_file(before["name"], before["sequence"]):
            self.set_in_file([(before["name"], before["sequence"])], False)
            self.set_in_file([(after["name"], after["sequence"])], True)

    def import_primers(self, pairs, present=True):
        """Import to This Sequence: make the primers `pairs` [(name, sequence)] part of this file (the regular primer
        annotation Save writes), as one step in the history; `present` False takes them out of it again (they stay in the
        database). Returns how many changed."""
        new = {(n, s.upper()): (n, s) for n, s in pairs if ((n, s.upper()) in self.file_primer_set) != present}
        if not new:
            return 0

        def put(add):
            for key, pair in new.items():
                if add:
                    self.file_primer_set[key] = pair
                else:
                    self.file_primer_set.pop(key, None)
            self.match_primers()
        put(present)
        names = [n for n, _ in new.values()]
        what = f"primer {names[0]}" if len(names) == 1 else f"{len(names)} primers"
        label = f"Imported {what} to this sequence" if present else f"Removed {what} from this sequence"
        self.record_external(label, lambda: put(not present), lambda: put(present))
        self.doc.modified = True   # the file stores them now: Save has something to write
        self.changed.emit()
        return len(new)

    def primer_list(self):
        """[(name, sequence)] of the primers whose binding sites are looked for: the database's, the file's own, then `extra`."""
        out, seen = [], set()
        for name, seq in [(r["name"], r["sequence"]) for r in db().rows] + list(self.file_primer_set.values()) + list(self.extra):
            if (name, seq.upper()) not in seen:
                seen.add((name, seq.upper()))
                out.append((name, seq))
        return out

    def set_extra(self, primers, temporary=()):
        """The primers to look for besides the database's: `primers` (the file's own) and `temporary` ones (designed, not saved:
        they are drawn in their own colour and have no menu until they are in the database)."""
        primers = list(primers) + [p for p in temporary if p not in primers]
        keys = {(n, s.upper()) for n, s in temporary}
        if primers != list(self.extra) or keys != self.temporary_primers:
            self.extra, self.temporary_primers = primers, keys
            self.match_primers()

    def match_primers(self):
        """Recompute where the database primers bind. Small sequences: now. Large ones:
        in a background thread (the old sites stay drawn meanwhile)."""
        primers = self.primer_list()
        seq, circular, rule = self.doc.sequence(), self.doc.circular, binding()   # copies: the document may be edited meanwhile
        seed, exact, mm = rule
        warm = tm_rule()
        # Only the minimum Tm (or the PCR conditions) changed: the sites it hides keep their choices, for when they come back.
        prune = self._matched_with is None or self._matched_with[1] == warm
        self._matched_with = (rule, warm)
        seqs = [p for _, p in primers]
        self._primer_job.run(lambda: find_primers(seq, circular, seqs, mm, seed, exact, **warm),
                             lambda hits: self._set_sites(primers, hits, rule, prune), self.length <= SYNC_LIMIT)

    def _set_sites(self, primers, hits, rule, prune=True):
        n = self.length
        sites = []
        in_db = self._db_primer_keys()
        self.site_totals = {name: len(found) for (name, _), found in zip(primers, hits)}   # (the sites drawn are the first MAX_SITES of each)
        for (name, seq), found in zip(primers, hits):
            for i, (s, e, strand, mm) in enumerate(found[:MAX_SITES]):
                site = SimpleNamespace(key=f"{name}#{i}", label=name, sequence=seq, strand=strand, start=s, end=e,
                                       spans=[(s, e)] if e <= n else [(s, n), (0, e - n)],
                                       mismatches=mm, tail=max(0, len(seq) - (e - s)))
                key = (name, seq.upper())
                site.from_db = key in in_db      # in the primer database
                site.from_file = self.primer_from_file(name, seq, in_db)
                site.temporary = key in self.temporary_primers and key not in self.file_primer_set   # designed, not imported
                # Drawn in the temporary style: designed here, or only in the database (not in this file).
                site.provisional = site.temporary or not site.from_file
                site.extent = primer_extent(site, site.tail, n if self.circular else None)
                site.mismatch_at = self._mismatch_positions(site) if mm else []
                sites.append(site)
        best = {}
        for site in sites:
            site.rank = site_rank(site, n)
            best[site.label] = min(best.get(site.label, site.rank), site.rank)
        for site in sites:
            site.best = site.rank == best[site.label]
        self.primer_sites = sites
        self._keep_site_choices(list(rule), primers, prune)
        self._say_cut_primers()
        self.match_elsewhere()   # (no-op unless the file has several records, and the primers are the same as last time)
        self.site_by_key = {s.key: s for s in sites}
        self._pack_primers()
        if self.current_primer not in self.site_by_key:
            self.current_primer = None
        self.primers_changed.emit()

    _matched_with = None  # (binding(), tm_rule()) of the last match_primers
    site_totals = {}      # primer name -> every site found on this sequence (more than are drawn when above MAX_SITES)
    away_totals = {}      # primer name -> every site found on the other records of the file
    away_capped = set()   # primer names with more than MAX_SITES on one of those records
    _cut_said = None

    def cut_primers(self):
        """{primer name: sites on all the sequences} for the primers that bind more places than are drawn (MAX_SITES a sequence)."""
        cut = {n for n, c in self.site_totals.items() if c > MAX_SITES} | self.away_capped
        return {n: self.site_totals.get(n, 0) + self.away_totals.get(n, 0) for n in cut}

    def _say_cut_primers(self):
        """The status line says which primers have sites that are not drawn (once, until the answer changes)."""
        cut = self.cut_primers()
        if cut == self._cut_said:
            return
        self._cut_said = cut
        if len(cut) == 1:
            (name, count), = cut.items()
            self.error.emit(f"Primer {name} binds {count:,} places: the first {MAX_SITES} on each sequence are shown.")
        elif cut:
            self.error.emit(f"{len(cut)} primers bind more than {MAX_SITES} places ({', '.join(sorted(cut)[:4])}{'…' if len(cut) > 4 else ''}): "
                            f"the first {MAX_SITES} of each on each sequence are shown.")

    @staticmethod
    def _db_primer_keys():
        return {(r["name"], r["sequence"].upper()) for r in db().rows}

    def primer_from_file(self, name, seq, in_db=None):
        """Is this primer "From this sequence" (Display filter) rather than "From the database"? In this file: stored in it /
        imported to it, designed here, or not from the database at all. A primer that is in the file and in the database counts
        as the file's (and is drawn as one). `in_db`: the database's (name, SEQUENCE) pairs, when already at hand."""
        key = (name, seq.upper())
        if key in self.file_primer_set or any(key == (n, q.upper()) for n, q in self.extra):
            return True
        return key not in (self._db_primer_keys() if in_db is None else in_db)

    def primer_source_shown(self, name, seq, in_db=None):
        """Is this primer's source switched on in the Display filter (From this sequence / From the database)?"""
        return self.show_file_primers if self.primer_from_file(name, seq, in_db) else self.show_db_primers

    def match_elsewhere(self):
        """Multi-record file: where the primers bind on the other records, found on a worker thread (each record is read from the
        file in turn). The result is `self.elsewhere`; the Primers panel lists it beside the sites on the open record."""
        if not self.multi_genbank:
            self.elsewhere, self.away_totals, self.away_capped = [], {}, set()
            return
        primers, (seed, exact, mm), warm = self.primer_list(), binding(), tm_rule()
        key = (str(self.path), tuple(primers), seed, exact, mm, tuple(sorted(warm.items())), tuple(self.fasta_records))
        if key == self._elsewhere_key:
            return
        self._elsewhere_key = key
        self._elsewhere_token += 1
        token, path, current, records = (id(self), self._elsewhere_token), self.path, self.fasta_record, list(self.fasta_records)
        seqs = [p for _, p in primers]

        def work():
            found, totals, capped = [], {}, set()
            try:
                for name, _ in records:
                    if name == current:
                        continue
                    doc = multigb.read_record(path, name)
                    hits = find_primers(doc.sequence(), doc.circular, seqs, mm, seed, exact, **warm)
                    for (label, seq), per_primer in zip(primers, hits):
                        totals[label] = totals.get(label, 0) + len(per_primer)
                        if len(per_primer) > MAX_SITES:
                            capped.add(label)
                        for s, e, strand, m in per_primer[:MAX_SITES]:
                            n = doc.length   # the template's bases (a site may run through the origin of a circle)
                            bases = doc.sequence(s, min(e, n)) + (doc.sequence(0, e - n) if e > n else "")
                            found.append(SimpleNamespace(record=name, label=label, sequence=seq, start=s, end=e, strand=strand, mismatches=m,
                                                         bases=bases))
                _relay.elsewhere_ready.emit(token, (found, totals, capped))
            except (OSError, RuntimeError):   # an unreadable file, or the program is closing
                pass
        QThreadPool.globalInstance().start(work)

    def _take_elsewhere(self, token, result):
        if token == (id(self), self._elsewhere_token):   # ignore the result for an older primer list
            self.elsewhere, self.away_totals, self.away_capped = result
            self._say_cut_primers()
            self.elsewhere_changed.emit()

    def _pack_primers(self):
        """Own lanes per strand for the drawn sites: forward primers above the bases, reverse below."""
        self.primer_lanes = {}
        best_only = tool_settings.value(BINDING, "best_only")
        for strand in (1, -1):
            group = [SimpleNamespace(index=s.key, spans=s.extent) for s in self.primer_sites
                     if s.strand == strand and self.site_drawn(s, best_only)]
            self.primer_lanes.update(pack_lanes(group)[0])

    def set_primer_hidden(self, name, hidden):
        """Hide in This File / Show in This File: the primer stays in the database; only this file stops drawing it."""
        self.set_primers_hidden([name], hidden)

    def _settings_changed(self, tool):
        """Settings > Primers: another binding rule finds other sites (matched again); the rest (Show only the best binding site by default, the PCR
        conditions, the 3′ check, the genome) only changes what is drawn or checked."""
        if tool in (BINDING, PCR) and (binding(), tm_rule()) != self._matched_with:   # (the minimum Tm: PCR too)
            self.match_primers()
        elif tool in (BINDING, PCR):
            self._pack_primers()
            self.primers_changed.emit()
            self.display_changed.emit()

    @staticmethod
    def _site_default(site, best_only):
        """Shown unless chosen otherwise: every site, or with Show only the best binding site by default the best ones."""
        return site.best or not best_only

    def site_shown(self, site, best_only=None):
        """Is this site ticked? As chosen by hand, else its default. `best_only`: the setting, when the caller already read it."""
        choice = self.site_choices.get(site.label, {}).get("sites", {}).get(site_key(site))
        if choice is not None:
            return choice
        return self._site_default(site, tool_settings.value(BINDING, "best_only") if best_only is None else best_only)

    def site_drawn(self, site, best_only=None):
        return site.label not in self.hidden_primers and self.site_shown(site, best_only)

    def set_sites_shown(self, sites, shown):
        """Show or hide these sites (right-click > Hide This Site, the Primers view's checkboxes); the primer stays. A choice that
        is what the site would be anyway is not kept."""
        best_only = tool_settings.value(BINDING, "best_only")
        for site in sites:
            entry = self.site_choices.setdefault(site.label, {"sequence": site.sequence.upper(), "sites": {}})
            if shown == self._site_default(site, best_only):
                entry["sites"].pop(site_key(site), None)
            else:
                entry["sites"][site_key(site)] = shown
            if not entry["sites"]:
                del self.site_choices[site.label]
        if not shown and self.current_primer in self.site_by_key and not self.site_drawn(self.site_by_key[self.current_primer]):
            self.current_primer = None
        self._pack_primers()
        self.primers_changed.emit()
        self.display_changed.emit()
        self._choice_changed()   # (after the repack: the views redraw on it)

    def saved_site_choices(self):
        """The site choices as a project keeps them (JSON)."""
        return {"rule": self.site_rule, "primers": {name: {"sequence": e["sequence"], "sites": [[*k, v] for k, v in e["sites"].items()]}
                                                    for name, e in self.site_choices.items()}}

    def _choice_changed(self):
        """A primer or site was shown or hidden: the file saves it (_write_display), so it has something to save. Not a FASTA file
        (no place for it)."""
        if not self.fasta_record or self.multi_genbank:
            self.doc.modified = True
            self.modified_changed.emit()   # (not `changed`: nothing was edited, every view would redo its layout)

    def display_choices(self):
        """The primers and sites hidden in this file, as a project (exports.project_meta) and the file's DISPLAY_NOTE keep them."""
        return {"hidden_primers": sorted(self.hidden_primers), "primer_sites": self.saved_site_choices()}

    def apply_display_choices(self, saved):
        """Take the primers and sites hidden that a project or the file kept (display_choices); not a change to save."""
        self.hidden_primers = {str(n) for n in saved.get("hidden_primers") or []}
        if isinstance(saved.get("primer_sites"), dict):
            self.load_site_choices(saved["primer_sites"])
        self._pack_primers()
        self.primers_changed.emit()
        self.display_changed.emit()

    def _read_display(self):
        """The primers and sites hidden when the file was saved (its DISPLAY_NOTE comment)."""
        for comment in self.doc.comments:
            if comment.startswith(DISPLAY_NOTE):
                try:   # (the writer wraps long lines: the note has no space of its own, see _write_display)
                    self.apply_display_choices(json.loads("".join(comment[len(DISPLAY_NOTE):].split())))
                except (ValueError, TypeError, AttributeError):
                    pass   # a hand-edited note: nothing hidden

    def _write_display(self):
        """The primers and sites hidden, as the document's DISPLAY_NOTE comment (written with the file), or none when nothing is.
        Spaces in names are escaped, so the only spaces are the ones added after commas for the writer to wrap at."""
        kept = [c for c in self.doc.comments if not c.startswith(DISPLAY_NOTE)]
        if self.hidden_primers or self.site_choices:
            saved = json.dumps(self.display_choices(), separators=(",", ":"))
            kept.append(DISPLAY_NOTE + " " + saved.replace(" ", r"\u0020").replace(",", ", "))
        if kept != self.doc.comments:   # (autosave calls this on its thread, after the GUI thread did: then nothing to write)
            self.doc.comments = kept

    def load_site_choices(self, saved):
        """Take the choices of a project; dropped when they were made under another binding rule."""
        try:
            choices = {name: {"sequence": str(e["sequence"]).upper(), "sites": {(int(a), int(b)): bool(v) for a, b, v in e["sites"]}}
                       for name, e in (saved.get("primers") or {}).items()}
        except (TypeError, ValueError, KeyError, AttributeError):   # a hand-edited file
            choices = {}
        rule = list(binding())
        if saved.get("rule") == rule and choices:
            self.site_rule, self.site_choices = rule, choices
            if self.primer_sites:
                self._pack_primers()
                self.primers_changed.emit()

    def _keep_site_choices(self, rule, primers, prune=True):
        """After a new match: under another rule nothing is kept; a primer whose sequence changed loses its choices; a site that no
        longer binds (the sequence was edited) loses its own, unless not `prune` (only the minimum Tm changed: it may come back)."""
        if rule != self.site_rule:
            self.site_rule, self.site_choices = rule, {}
        sequences = {name: seq.upper() for name, seq in primers}
        keys = {}
        for site in self.primer_sites:
            keys.setdefault(site.label, set()).add(site_key(site))
        for name in list(self.site_choices):
            entry = self.site_choices[name]
            kept = {k: v for k, v in entry["sites"].items() if not prune or k in keys.get(name, ())}
            entry["sites"] = kept
            if sequences.get(name) != entry["sequence"] or not kept:
                del self.site_choices[name]

    def set_primers_hidden(self, names, hidden):
        """set_primer_hidden for several primers at once (one redraw)."""
        names = set(names)
        before = set(self.hidden_primers)
        self.hidden_primers = before | names if hidden else before - names
        if self.hidden_primers == before:
            return
        if hidden and self.current_primer and self.site_by_key.get(self.current_primer) is not None \
                and self.site_by_key[self.current_primer].label in names:
            self.current_primer = None
        self._pack_primers()
        self.primers_changed.emit()
        self.display_changed.emit()
        self._choice_changed()   # (after the repack: the views redraw on it)

    def _mismatch_positions(self, site):
        """Template positions where the annealed part of the primer does not match."""
        positions = [p for a, b in site.spans for p in range(a, b)]
        template = "".join(self.doc.sequence(a, b) for a, b in site.spans).upper()
        expected = site.sequence[-len(template):].upper()
        if site.strand < 0:
            expected = reverse_complement(expected)  # compare on the top strand
        return [positions[k] for k in range(len(template)) if template[k] != expected[k]]

    def select_primer(self, key, extend=False):
        """Select a primer site; with `extend` (shift-click) and a primer on the other strand
        already selected, select the stretch they amplify instead."""
        site = self.site_by_key.get(key)
        if site is None:
            return
        self.primer_clicked.emit(key)
        other = self.site_by_key.get(self.current_primer) if extend else None
        if other is not None and other.strand != site.strand:
            fwd, rev = (other, site) if other.strand > 0 else (site, other)
            s, e = fwd.spans[0][0], rev.spans[-1][1]
            if e <= s and self.circular:
                e += self.length  # the product runs through the origin
            if s < e:
                self.select(s, e)
                self.reveal.emit(self.start)
                return
        s, e = site.spans[0][0], site.spans[-1][1]
        if e < s:  # wraps through the origin: all of it on a circular sequence, else the first part
            e = e + self.length if self.circular else site.spans[0][1]
        self.select(s, e, primer=key)
        self.reveal.emit(self.start)

    def _changed(self, keep_lanes=False):
        with timed("editor: refresh after an edit"):
            self._refresh(keep_lanes)
        self.version += 1
        self.changed.emit()

    def set_display(self, enzymes=None, translation=None, hidden_kinds=None, primers=None,
                    cds_translation=None, file_primers=None, db_primers=None, units=None, enzyme_filter=None, cds_only=None,
                    complement=None):
        if enzymes is not None and enzymes != self.show_enzymes:
            self.show_enzymes = enzymes
            self.find_cut_sites()
        if enzyme_filter is not None:
            self.enzyme_filter = enzyme_filter
        if translation is not None:
            self.show_translation = translation
        if complement is not None:
            self.show_complement = complement
        if primers is not None:
            self.show_primers = primers
        if file_primers is not None:
            self.show_file_primers = file_primers
        if db_primers is not None:
            self.show_db_primers = db_primers
        if cds_translation is not None:
            self.show_cds_translation = cds_translation
        if hidden_kinds is not None or units is not None or cds_only is not None:
            if cds_only is not None:
                self.cds_only = cds_only
            if units is not None:
                self.show_units = units
            if hidden_kinds is not None:
                self.hidden_kinds = set(hidden_kinds)
            # only what is drawn changes: the features, their genes, primer sites and cut sites stay (_refresh read them all
            # again from the core, half a second on a genome, and dropped the Local Align hits)
            self._lay_out()
        self.display_changed.emit()

    def set_cut_marks(self, owner, marks):
        """Show cuts on the sequence: `marks` = [(top, bottom, label[, tip[, pinned]])], the gap between two bases where each
        strand is cut (equal for a blunt cut); `tip`: the hover text (rich text), or a function that makes it; `pinned`: drawn
        emphasised with its label (the name) written above, as if it were being hovered (the Digest's selected fragment).
        Each tool sets its own (`owner`); an empty list removes them."""
        n = self.length
        fixed = [SimpleNamespace(top=m[0] % n if self.circular and n else m[0], bottom=m[1] % n if self.circular and n else m[1],
                                 label=m[2], tip=m[3] if len(m) > 3 else None, pinned=bool(m[4]) if len(m) > 4 else False)
                 for m in marks]
        if fixed:
            self.cut_marks[owner] = fixed
        elif owner not in self.cut_marks:
            return
        else:
            del self.cut_marks[owner]
        try:
            self.display_changed.emit()
        except RuntimeError:
            pass   # the editor is already being deleted with its window (a tool clearing its cuts as it closes)

    def set_search_hits(self, spans):
        """The places the search bar found, which the views highlight: [(start, end)] or [(start, end, strand)] (end may pass
        the length on a circular sequence: it is split in two here). The strand (1 top, -1 bottom, 0 both) tells the sequence
        view which of the two strands to colour. An empty list clears them."""
        n = self.length
        out, strands = [], {}
        for hit in spans:
            s, e = hit[0], hit[1]
            strand = hit[2] if len(hit) > 2 else 0
            for piece in ([(s, n), (0, e - n)] if e > n and n else [(s, e)]):
                out.append(piece)
                # the same place found on both strands (a palindrome) is both
                strands[piece] = strand if strands.get(piece, strand) == strand else 0
        out.sort()
        changed_strands = strands != self.search_strand
        self.search_strand = strands
        if out != self.search_spans or changed_strands:
            self.search_spans = out
            self._search_reach = max((e - s for s, e in out), default=0)
            try:
                self.search_changed.emit()
            except RuntimeError:
                pass   # the editor is being deleted with its window

    def hits_between(self, a, b):
        """The highlighted search spans that overlap [a, b), in order (the views ask for each row they draw)."""
        spans = self.search_spans
        i = bisect_left(spans, (a - self._search_reach, -1))
        out = []
        while i < len(spans) and spans[i][0] < b:
            if spans[i][1] > a:
                out.append(spans[i])
            i += 1
        return out

    def shown_sites(self):
        """[(cut position, enzyme name)] to draw (the display filter's choice of enzymes), sorted by position. None when the filter
        lets through more than MAX_DRAWN_CUTS (drawing them froze the window): `withheld_sites` says how many."""
        return self._sites()[0]

    def withheld_sites(self):
        """How many sites the display filter lets through but are not drawn (more than MAX_DRAWN_CUTS), 0 when they all are."""
        return self._sites()[1]

    def _sites(self):
        if not self.show_enzymes:
            return [], 0
        if self.enzyme_filter not in self._shown_sites:
            enzymes, mine = {e.name: e for e in COMMON}, inventory.owned()
            sites = sorted((pos, name) for name, cuts in self.site_cuts.items()
                           if wanted(self.enzyme_filter, enzymes[name], len(cuts), is_owned(enzymes[name], mine)) for pos in cuts)
            self._shown_sites[self.enzyme_filter] = ([], len(sites)) if len(sites) > MAX_DRAWN_CUTS else (sites, 0)
        return self._shown_sites[self.enzyme_filter]

    def _inventory_changed(self):
        """An enzyme was tagged / untagged as in the user's inventory: the filter may draw other sites."""
        self._shown_sites = {}
        if self.enzyme_filter.inventory is not None:
            self.display_changed.emit()

    @property
    def show_primers(self):
        """Primers are drawn at all (either source switched on)."""
        return self.show_file_primers or self.show_db_primers

    @show_primers.setter
    def show_primers(self, on):
        self.show_file_primers = self.show_db_primers = bool(on)

    @property
    def drawn_primers(self):
        """Primer sites to draw: never those hidden in this file, and only from the sources switched on in the display
        filter: "From this sequence" for the primers in this file (also when the database has them too), "From the database"
        for the ones that are only there."""
        best_only = tool_settings.value(BINDING, "best_only")   # (read once: the views ask for this list on every paint)
        return [s for s in self.primer_sites if self.site_drawn(s, best_only)
                and (self.show_file_primers if s.from_file else self.show_db_primers)]

    def visible(self, f):
        """Drawn, as far as the display filter goes? Transcription units have a switch of their own (their type, mRNA,
        is not what they are)."""
        if feature_key(f) in self.hidden_features:
            return False
        if self.cds_only and f.kind == GENE and f.index in self.redundant_genes:
            return False
        return self.show_units if self.is_unit(f) else f.kind not in self.hidden_kinds

    def is_unit(self, f):
        """A transcription unit? Known for the document's own features (see _refresh); a tool's temporary one is asked."""
        if self.by_index.get(f.index) is f:
            return f.index in self.unit_indexes
        return is_tu(f)

    def set_features_hidden(self, keys, hidden):
        """Hide / show (right-click > Hide, Show) the features with these `feature_key`s: only the drawing changes."""
        keys = set(keys)
        before = set(self.hidden_features)
        self.hidden_features = before | keys if hidden else before - keys
        if self.hidden_features != before:
            if self.current_feature is not None and hidden and feature_key(self.by_index[self.current_feature]) in keys:
                self.current_feature = None
            self._lay_out()
            self.display_changed.emit()

    def visible_features(self, start, end, overlays=True):
        """Features drawn in [start, end) (not those of hidden types); with `overlays`, the tools' temporary ones too."""
        out = [f for f in self.doc.features_in(start, end) if f.kind not in HIDDEN_KINDS and self.visible(f)]
        if overlays:
            out += [o for o in self.overlay_features
                    if self.visible(o) and any(s < end and e > start for s, e in o.spans)]
        return out

    # ---- temporary features (a tool's results, shown without being added) ---------------------------------------
    overlays = {}           # owner -> [feature-like SimpleNamespace]; set per instance by set_overlay
    overlay_features = []   # all of them
    overlay_index = {}      # negative index -> feature

    def set_overlay(self, owner, features):
        """Show temporary features in every view of this document without adding them to it (the Plasmid Checkup's
        transcription units and found elements): `features` = [(kind, spans, strand, qualifiers)], spans as on the
        sequence (two parts across the origin). They are drawn and filtered like the others but are not in `features`,
        the feature table, the file or the history, and cannot be edited (their index is negative). [] removes them."""
        self.overlays = {k: v for k, v in self.overlays.items() if k != owner}
        made = []
        for kind, spans, strand, qualifiers in features:
            label = next((v for k, v in qualifiers if k == "label" and v), kind)
            parts = ",".join(f"{s + 1}..{e}" for s, e in spans)
            loc = parts if len(spans) == 1 else f"join({parts})"
            made.append(SimpleNamespace(index=0, kind=kind, label=label, spans=list(spans), strand=strand,
                                        qualifiers=list(qualifiers), location=f"complement({loc})" if strand < 0 else loc,
                                        temporary=True))
        if made:
            self.overlays[owner] = made
        self.overlay_features = [f for fs in self.overlays.values() for f in fs]
        for i, f in enumerate(self.overlay_features):
            f.index = -1 - i
        self.overlay_index = {f.index: f for f in self.overlay_features}
        self._lay_out()
        self.display_changed.emit()

    def feature_by_index(self, index):
        """A feature of the document, or a temporary one (negative index)."""
        return self.by_index.get(index) or self.overlay_index.get(index)

    def translated(self, f):
        """Draw this CDS's amino acids inside its arrow?"""
        return f.kind == "CDS" and self.show_cds_translation and (f.location, f.label) not in self.untranslated

    def toggle_translation(self, index):
        f = self.by_index.get(index)
        if f is not None:
            self.untranslated ^= {(f.location, f.label)}  # ^= : add if absent, remove if present
            self.display_changed.emit()

    # ---- selection -----------------------------------------------------------------
    def set_query_hits(self, hits):
        """Show Local Align hits (empty list = clear them) in the genome view."""
        self.query_hits, self.current_query_hit = hits, -1
        self.display_changed.emit()

    def select_site(self, name, pos):
        """Select the recognition site of enzyme `name` that cuts the top strand at gap `pos` (a unique cutter on the map), and
        bring the views to it; a cursor at the cut when the site cannot be found."""
        from cut_sites import hover_site
        site = hover_site(self, ("site", name, pos))
        if site is None:
            self.set_cursor(pos)
        else:
            self.select(site[0], site[1])
        self.reveal.emit(self.start)

    def select_query_hit(self, index):
        """Pick Local Align hit `index` (its row in the Alignment view and in the dialog's table): select its bases."""
        if not 0 <= index < len(self.query_hits):
            return
        h = self.query_hits[index]
        self.current_query_hit = index
        self.select(h.start, h.end)
        self.reveal.emit(h.start)
        self.display_changed.emit()
        self.query_hit_changed.emit(index)

    def select(self, start, end, feature=None, primer=None):
        """Select [start, end). On a circular sequence the selection may pass the origin:
        then end > length (start < length < end), like a GenBank join(s..n,1..e). Positions
        of any turn (e.g. the sequence view's second copy, or -1) are brought back to the
        first one, and a selection never covers more than one full copy."""
        n = self.length
        if self.circular and n:
            start, end = sorted((start, end))
            end = min(end, start + n)
            # Into 0 <= start <= n; start == n only for a cursor at the very end (so typing
            # there appends, not inserts at base 1).
            turn = (start - 1) // n * n if start > 0 else start // n * n
            if start - turn == n and end > start:
                turn += n
            start, end = start - turn, end - turn
        else:
            start, end = sorted((max(0, min(start, n)), max(0, min(end, n))))
        self.start, self.end, self.current_feature, self.current_primer = start, end, feature, primer
        self.selection_changed.emit()

    @property
    def circular(self):
        return self.doc.circular

    @property
    def wraps(self):
        """The selection passes the origin."""
        return self.end > self.length

    def selected_spans(self):
        """The selection as [(start, end)] on the sequence: two parts if it passes the origin."""
        n = self.length
        return [(self.start, n), (0, self.end - n)] if self.end > n else [(self.start, self.end)]

    def set_cursor(self, pos):
        self.select(pos, pos)

    def select_feature(self, index):
        f = self.feature_by_index(index)
        if f is None:
            return
        start, end = feature_bounds(f)
        starts = [s for s, _ in f.spans]
        if self.circular and starts != sorted(starts):  # crosses the origin: all of it
            start, end = f.spans[0][0], self.length + f.spans[-1][1]
        # A temporary feature selects its bases only: there is no feature in the document to edit or delete.
        self.select(start, end, feature=index if index in self.by_index else None)
        self.reveal.emit(self.start)

    def selected_sequence(self):
        return "".join(self.doc.sequence(a, b) for a, b in self.selected_spans())

    # ---- editing -------------------------------------------------------------------
    def dna_editable(self):
        """False (with a status message) while the DNA is locked: every change to the bases checks this first."""
        if self.dna_locked:
            self.error.emit("The DNA is locked: unlock it (lock button of the toolbar) to change bases.")
        return not self.dna_locked

    def _edit(self, action, cursor):
        try:
            action()
        except (OSError, ValueError) as e:
            self.error.emit(str(e))
            return False
        self._changed(keep_lanes=True)
        self.set_cursor(cursor)
        self.reveal.emit(cursor)
        return True

    def type_bases(self, text):
        """Keyboard input: insert at the cursor, replacing the selection. Only A, C, G, T
        and N are accepted; anything else is refused with a status message."""
        if not self.dna_editable():
            return False
        if not text or any(c not in TYPABLE for c in text):
            self.error.emit("Only A, C, G, T and N can be typed (paste for ambiguity codes)")
            return False
        return self._replace_selection(text)

    def _replace_selection(self, bases, features=None):
        """Replace the selection with `bases` (+ pasted `features`), cursor after them.
        A selection across the origin is one edit too (Rust squashes its two parts)."""
        if not self.dna_editable():
            return False
        s, e = self.start, self.end
        if e > self.length:
            at = []
            done = self._edit(lambda: at.append(self.doc.edit_across_origin(s, e, bases, features)), 0)
            if done:
                self.set_cursor(at[0] + len(bases))
                self.reveal.emit(self.start)
            return done
        if features is not None:
            return self._edit(lambda: self.doc.paste(s, e, bases, features), s + len(bases))
        return self._edit(lambda: self.doc.replace(s, e, bases), s + len(bases))

    def overwrite_bases(self, pos, bases):
        """Replace the bases from `pos` by `bases` of the same length (a recoded CDS): nothing moves, so no feature changes. One
        undo step. The new bases stay selected."""
        if not self.dna_editable():
            return False
        n = self.doc.length

        def write():
            if pos + len(bases) <= n:
                self.doc.overwrite(pos, bases)
            else:   # through the origin of a circular sequence: in two pieces (Undo takes two steps)
                self.doc.overwrite(pos, bases[:n - pos])
                self.doc.overwrite(0, bases[n - pos:])
        done = self._edit(write, pos)
        if done:
            self.select(pos, pos + len(bases))
        return done

    def copy(self):
        """Selected bases as text (for the clipboard), remembering their features."""
        bases = self.selected_sequence()
        features = [(f.kind, f.location, f.qualifiers) for f in self.doc.copy_features(self.start, self.end)]
        Editor.clip = (bases, features)
        return bases

    def copy_mime(self):
        """Clipboard content for a copy: the bases as plain text (for any other program),
        plus the bases and features as JSON under CLIP_MIME, so pasting in another
        openDNA window, even one in another running copy of the app, keeps them."""
        self.copy()
        return self._clip_mime()

    @staticmethod
    def _clip_mime():
        mime = QMimeData()
        mime.setText(Editor.clip[0])
        mime.setData(CLIP_MIME, json.dumps(Editor.clip).encode())
        return mime

    def copy_reverse_complement_mime(self):
        """Like `copy_mime`, but the clipboard holds the REVERSE COMPLEMENT of the selection, and its features come
        along flipped onto it (strand inverted, join parts reordered), so pasting it gives the same annotations on the
        other strand."""
        bases = reverse_complement(self.selected_sequence())
        features = [(f.kind, f.location, f.qualifiers) for f in self.doc.copy_features_reversed(self.start, self.end)]
        Editor.clip = (bases, features)
        return self._clip_mime()

    def paste_mime(self, mime):
        clip = None
        if mime.hasFormat(CLIP_MIME):
            try:
                bases, features = json.loads(bytes(mime.data(CLIP_MIME)).decode())
                clip = (bases, [(k, loc, [tuple(q) for q in quals]) for k, loc, quals in features])
            except (ValueError, TypeError):
                clip = None  # malformed: paste the text alone
        return self.paste(mime.text(), clip)

    def paste(self, text, clip=None):
        """Paste clipboard text over the selection. If it is what was copied (`clip`, or
        the last copy in this app), the copied features come along (one undo step)."""
        clip = clip or Editor.clip
        if clip and text == clip[0] and text:
            return self._replace_selection(text, clip[1])
        try:
            bases = clean_pasted(text)
        except ValueError as err:
            self.error.emit(str(err))
            return False
        if not bases:
            return False
        return self._replace_selection(bases)

    def backspace(self):
        if not self.dna_editable():
            return False
        s, e = self.start, self.end
        if e > self.length:
            return self._replace_selection("")
        if s == e:
            if s == 0:
                return False
            s -= 1
        return self._edit(lambda: self.doc.delete(s, e), s)

    def delete_forward(self):
        if not self.dna_editable():
            return False
        s, e = self.start, self.end
        if e > self.length:
            return self._replace_selection("")
        if s == e:
            if e >= self.length:
                return False
            e += 1
        return self._edit(lambda: self.doc.delete(s, e), s)

    def undo(self):
        return self._history(self._undo_one)

    def redo(self):
        return self._history(self._redo_one)

    def record_external(self, label, undo, redo):
        """A step in the history for something that is not in the document (a change to the primer database): `label` is its
        History text; Undo / Redo reach it in order with the edits and call `undo()` / `redo()`, which return a message
        (shown in the status bar) when they could not do it, else None. It never makes the document modified."""
        self._external_count += 1
        token = f"x{self._external_count}"
        self.externals[token] = (undo, redo)
        self.doc.external(label, token)
        self.version += 1
        self.changed.emit()   # the History panel; nothing in the sequence or features moved, so no full refresh

    def _step(self, undo):
        """One undo or redo step of the document; when it is an external one, also what it stands for."""
        token = self.doc.next_external(not undo)
        if not (self.doc.undo() if undo else self.doc.redo()):
            return False
        if token is not None:
            run = self.externals.get(token)
            message = run[0 if undo else 1]() if run else None
            if message:
                self.error.emit(message)
        return True

    def _undo_one(self):
        return self._step(True)

    def _redo_one(self):
        return self._step(False)

    def undo_steps(self, count):
        """Several undo (or redo, count < 0) steps at once, with one refresh (History panel)."""
        step = self._undo_one if count > 0 else self._redo_one
        before = self.doc.sequence() if self.dna_locked else None
        done = 0
        while done < abs(count) and step():
            done += 1
        if done and before is not None and self.doc.sequence() != before:
            back = self._redo_one if count > 0 else self._undo_one   # it would change the bases: step back
            for _ in range(done):
                back()
            self.dna_editable()
            return 0
        if done:
            self._changed()
            self.set_cursor(min(self.start, self.length))
        return done

    def _history(self, step):
        before = self.doc.sequence() if self.dna_locked else None
        if not step():
            return False
        if before is not None and self.doc.sequence() != before:   # locked: a step that changes bases is refused
            (self._redo_one if step == self._undo_one else self._undo_one)()
            self.dna_editable()
            return False
        self._changed()
        self.set_cursor(min(self.start, self.length))
        return True

    # ---- features and whole-sequence operations -------------------------------------
    def set_feature(self, index, kind, location, qualifiers):
        """Add (index None) or replace a feature, then select it.
        Raises OSError/ValueError on bad input, so a dialog can show it and stay open."""
        index = self.doc.set_feature(index, kind, location, qualifiers)
        self._changed()
        self.select_feature(index)
        return index

    def edit_features(self, changes):
        """Several feature changes as one undo step (see Document.edit_features)."""
        self.doc.edit_features(changes)
        self._changed()
        self.set_cursor(self.start)

    def add_features(self, features):
        """Add (kind, location, qualifiers) tuples as one undo step; selects the first."""
        first = self.doc.add_features(features)
        self._changed()
        self.select_feature(first)
        return first

    def delete_feature(self, index):
        self.doc.delete_feature(index)
        self._changed()
        self.set_cursor(self.start)

    def recolor_by_type(self):
        """Every annotation takes its type's standard colour (written as .dna/ApE colour
        qualifiers, so other programs show it too). One undo step; returns how many changed."""
        changes = []
        for f in self.features:
            color = type_color(f.kind).name()
            qualifiers = [(k, v) for k, v in f.qualifiers if k not in COLOR_KEYS] + [(k, color) for k in COLOR_KEYS]
            if sorted(qualifiers, key=str) != sorted(f.qualifiers, key=str):
                changes.append((f.index, (f.kind, f.location, qualifiers)))
        if changes:
            self.edit_features(changes)
        return len(changes)

    def refresh_translations(self, indexes=None):
        """Recompute the /translation of CDSs that have one, from their own bases, /codon_start and /transl_table (else the
        table chosen in Tools > Settings): Rust does the work (Document.refresh_translations), as ONE undo step, and only
        when something differs. `indexes`: feature indexes, or None for every CDS. Returns [(index, status, reason)]:
        status "updated", "unchanged", or why it was left alone ("pseudo", "partial", "other")."""
        report = self.doc.refresh_translations(indexes, genetic_code.table_id())
        if any(status == "updated" for _, status, _ in report):
            self._changed()
        return report

    def features_cut_at(self, pos):
        """Features that breaking a circular sequence just before base `pos` (0-based) would
        cut in two."""
        out = []
        for f in self.features:
            starts = [s for s, _ in f.spans]
            wraps = starts != sorted(starts)  # crosses the current origin
            if any(s < pos < e for s, e in f.spans) or (pos == 0 and wraps):
                out.append(f)
        return out

    def linearize(self, pos):
        """Make a circular sequence linear, broken just before base `pos` (0-based): that base
        becomes base 1 (a rotation) and the topology is set to linear: ONE undo step."""
        if not self.doc.circular or not self.dna_editable():
            return
        self.doc.linearize(pos % max(1, self.length))
        self._changed()
        self.set_cursor(0)

    def set_origin(self, pos):
        """Rotate a circular sequence so base `pos` (0-based) becomes base 1, features along: ONE undo step. False when nothing
        changed (linear, locked, or `pos` is already the origin)."""
        pos %= max(1, self.length)
        if not self.doc.circular or not self.dna_editable() or pos == 0:
            return False
        self.doc.rotate(pos)
        self._changed()
        self.set_cursor(0)
        return True

    def set_circular(self, circular):
        """Plasmid or linear: changes the map, cut sites and primer sites across the origin."""
        if circular != self.doc.circular:
            self.doc.set_topology(circular)  # an undo step, like any edit
            self._changed()

    def reverse_complement(self):
        if not self.dna_editable():
            return
        start, end = self.start, self.end
        self.doc.reverse_complement()
        self._changed()
        self.select(self.length - end, self.length - start)  # same bases, now flipped

    def find(self, query, after):
        """Next match after position `after` (wrapping to the top): (start, end, feature
        index or None), or None. Searches DNA on both strands when the query is only A, C,
        G, T, N (so a name like "SsrA" is looked up as a name only), and always names:
        feature label/gene/locus_tag/product/..., primers, and items of loaded annotation
        tracks. Protein motifs: the ORFs tool searches translations.
        shortcut: exact letters only (N in the query matches only N); expand IUPAC codes
        if primer-style searches are needed."""
        q = query.strip()
        if not q:
            return None
        hits = []
        if set(q.upper()) <= set("ACGTN"):
            seq, pat = self.doc.sequence().upper(), q.upper()
            text = seq + seq[:len(pat) - 1] if self.doc.circular else seq   # (a circle read on past its origin: a query may run through it)
            for p in {pat, reverse_complement(pat)}:
                for start in (after + 1, 0):  # next one after the cursor, else the first
                    i = text.find(p, start)
                    if i != -1 and i < len(seq):   # (a match wholly in the copied tail is the one at the front)
                        hits.append((i, i + len(p), None))
        needle = q.lower()
        if self._search_text is None:  # built on the first search after a change
            self._search_text = [(" ".join(v for k, v in f.qualifiers if k in NAME_KEYS and v) + " " + f.label).lower()
                                 for f in self.features]
        for f, text in zip(self.features, self._search_text):
            if needle in text:
                hits.append((*feature_bounds(f), f.index))
        for t in self.tracks:
            found = t.track.search(q, TRACK_HITS, t.seqid)
            if len(found) >= TRACK_HITS:
                self.error.emit(f"{t.name}: more than {TRACK_HITS:,} items match: Next and Previous see the first {TRACK_HITS:,} of them.")
            hits += [(s, e, None) for _, s, e, *_ in found]
        for site in self.primer_sites:
            if needle in site.label.lower():
                hits.append((site.spans[0][0], site.spans[0][1], None))
        if not hits:
            return None
        later = [h for h in hits if h[0] > after]
        return min(later or hits)

    def file_primer_features(self):
        """primer_bind features to write in the GenBank file: every exact site of every primer that belongs to this file
        (stored in it or imported to it; matched now, a background result could be stale). A primer only in the database is
        not written: Import to This Sequence makes it part of the file. A stored primer with no exact site (a few
        mismatches) is kept as it was. (Every save writes these: the primers and sites hidden go with them, _write_display.)"""
        self._write_display()
        primers = list(self.file_primer_set.values())
        hits = find_primers(self.doc.sequence(), self.doc.circular, [p for _, p in primers], 0)   # (every one: not only those drawn)
        out, written = [], set()
        for (name, seq), found in zip(primers, hits):
            for s, e, strand, _ in found:   # (every exact site: a file keeps what other programs read from it)
                written.add((name, seq.upper()))
                out.append((PRIMER_KIND, site_location(s, e, strand, self.length),
                            [("label", name), ("note", f"sequence: {seq}")]))
        for f in self.file_primers:
            key = self._stored_key(f)
            if key in self.file_primer_set and key not in written:   # (not one taken out of the file: In file unticked)
                out.append((f.kind, f.location, f.qualifiers))
        return out

    def lookup(self, name):
        """Best match for a name (the genome view's "go to gene"): a feature whose label,
        gene, locus_tag... is exactly `name` (any case), else a track item with that exact
        name, else the first Find hit. (start, end, feature index or None) or None."""
        needle = name.strip().lower()
        for f in self.features:
            if needle == f.label.lower() or any(k in NAME_KEYS and v and v.lower() == needle for k, v in f.qualifiers):
                return (*feature_bounds(f), f.index)
        for t in self.tracks:
            for _, s, e, _, item_name, _ in t.track.search(name, TRACK_HITS, t.seqid):
                if item_name.lower() == needle:
                    return (s, e, None)
        return self.find(name, -1)

    def save(self, path):
        self.doc.save_genbank(path, self.file_primer_features())
        self.untitled = False
        self.path = Path(path)
        self.changed.emit()

    def gc_percent(self):
        if self.end <= self.start:
            return None
        # Rust: fast on any size; a selection across the origin = its two parts, weighted.
        spans = self.selected_spans()
        return 100.0 * sum(self.doc.gc_bins(a, b, 1)[0] * (b - a) for a, b in spans) / (self.end - self.start)