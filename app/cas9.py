"""Spacer Finder (Cas9 / Cas12a guides): good places to cut, on the open sequence and checked against a background genome.

Engine (no GUI): `candidates` lists every guide (protospacer + PAM) of a nuclease in a stretch of sequence, `score_guide`
rates it with a few transparent rules, `scan` counts the other places a guide matches (in the open sequence and in the
genome) and `verdict` turns the counts into one mark. `Cas9Dialog` is the tool: a table of guides, a button that checks
their specificity, and a button that adds one as a feature (the first AUTO_CHECK are checked by themselves, in the background). The background genome is the one of the Primers tool
(Primers > Options): both tools read the same setting, so they always agree on it.

Not a trained model: the score is a plain checklist (GC content, no TTTT, no long runs), not Doench's Rule Set 2.
"""
import re
import threading
import value_filter
from collections import namedtuple

from PySide6.QtCore import QObject, Qt, Signal
from PySide6.QtGui import QColor
from PySide6.QtWidgets import (
    QAbstractItemView, QComboBox, QHBoxLayout, QLabel, QMenu, QProgressBar, QPushButton, QSizePolicy, QSpinBox, QTableWidget,
    QTableWidgetItem, QVBoxLayout, QWidget,
)

import opendna_core as core

from cut_sites import CAS9_KEY, NUCLEASES, cas9_mark_tip, cas9_qualifier, cut_places
from editor import reverse_complement, site_location
from primer_panel import MARKS, OptionsDialog, genome, genome_label_text
from reference import as_reference
from tools import ToolDialog
from uikit import fit_columns, no_auto_default, plural, right_aligned, run_busy

# Nuclease and NUCLEASES (and cut_places) live in cut_sites.py (the view and the tooltips need them without this tool)
IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "N": "[ACGT]", "R": "[AG]", "Y": "[CT]", "S": "[CG]", "W": "[AT]",
         "K": "[GT]", "M": "[AC]", "B": "[CGT]", "D": "[AGT]", "H": "[ACT]", "V": "[ACG]"}
MAX_CANDIDATES = 300           # rows shown: the best by score
MAX_REGION = 2_000_000         # bp: longer stretches must be narrowed with a selection first
SEED = 12                      # bases next to the PAM where a mismatch matters most
AUTO_CHECK = 100               # the first guides of the list are checked against the reference without being asked
CHUNK = 10                     # guides scanned per step of that check: a result lands in the table (and a cancel is seen) at each
INDEX = Qt.UserRole + 1        # item data: which guide (index in Cas9Dialog.guides) a table row shows (rows may be sorted)

# One guide: where it sits, how it reads, where the nuclease cuts it (the top-strand base the cut falls before; `bottom`: the same
# for the bottom strand, which differs from `cut` for Cas12a) and its score.
Guide = namedtuple("Guide", "start end strand guide pam cut gc score notes bottom", defaults=(None,))


def gc_percent(seq):
    return 100 * sum(c in "GC" for c in seq.upper()) / max(len(seq), 1)


def score_guide(guide):
    """(score 0-100, notes) from a checklist of common design advice: GC content 40-70 %, no TTTT (it stops transcription
    from the U6 promoter that usually drives a guide) and no run of five identical bases. 100 = nothing to say against it.
    This is a rule of thumb, not a measured activity."""
    score, notes = 100, []
    gc = gc_percent(guide)
    if gc < 25 or gc > 80:
        score -= 40
        notes.append(f"GC {gc:.0f} % is extreme")
    elif gc < 40 or gc > 70:
        score -= 15
        notes.append(f"GC {gc:.0f} % is outside 40-70 %")
    if "TTTT" in guide.upper():
        score -= 40
        notes.append("TTTT ends transcription from a U6 promoter")
    if re.search(r"(.)\1{4,}", guide.upper()):
        score -= 15
        notes.append("run of 5 or more identical bases")
    return max(score, 0), notes


def pam_regex(pam):
    return "".join(IUPAC[c] for c in pam)


def candidates(seq, circular, nuclease, offset=0, wrap=False):
    """Every guide of `nuclease` in `seq` (a stretch starting at `offset` on the full sequence), on both strands, best
    first (score, then position). `wrap`: `seq` is a whole circular sequence, so sites across its origin count too."""
    seq = seq.upper()
    n, L, P = len(seq), nuclease.length, len(nuclease.pam)
    span = L + P
    text = seq + seq[:span - 1] if wrap and n >= span else seq
    out = []
    pam, guide_re = f"({pam_regex(nuclease.pam)})", f"([ACGT]{{{L}}})"
    pattern = f"(?={pam + guide_re if nuclease.pam_first else guide_re + pam})"
    for strand in (1, -1):
        s = text if strand > 0 else reverse_complement(text)
        for m in re.finditer(pattern, s):
            i = m.start()
            if wrap and i >= n:
                continue   # the same site starts again at the other end of the copy
            guide, pam_found = (m.group(2), m.group(1)) if nuclease.pam_first else (m.group(1), m.group(2))
            # A site that starts at `i` on the reverse complement starts at len(s) - i - span on the top strand.
            start = i if strand > 0 else len(s) - i - span
            if not wrap and start < 0:
                continue
            top, bottom = cut_places(nuclease, start, strand)
            score, notes = score_guide(guide)
            out.append(Guide(start + offset, start + span + offset, strand, guide, pam_found, top + offset,
                             round(gc_percent(guide)), score, notes, bottom + offset))
    out.sort(key=lambda g: (-g.score, g.start))
    return out


def seed_mismatches(mask, length):
    """How many of the mismatches in `mask` (bit i = guide base i, 0 = the 5' end) are in the SEED bases next to the PAM."""
    return bin(mask >> max(length - SEED, 0)).count("1")


Specificity = namedtuple("Specificity", "here genome genome_name")   # each: counts of other sites by mismatch count, or None


def _counts(hits, length, max_mm, own=None):
    """[sites with 0 mismatches, 1, ..., max_mm] from find_guide_sites hits, leaving out the guide's own site `own`."""
    counts = [0] * (max_mm + 1)
    for start, end, strand, mm, mask in hits:
        if own is not None and (start, strand) == own and mm == 0:
            continue
        counts[mm] += 1
    return counts


class Cancelled(Exception):
    """A scan stopped because `cancel()` said so."""


def scan(text, circular, name, guides, nuclease, max_mm=3, ref=None, cancel=None):
    """Specificity of each of `guides` (Guide tuples found in the sequence `text`, circular or not, called `name`): the other
    sites it matches there and in the background genome `ref` (a Document, or None), as counts by number of mismatches.
    A site that is the guide itself is left out. If `ref` is that same sequence it is not scanned twice. `cancel()` -> True
    stops the scan (between two sequences of the reference) by raising Cancelled."""
    n, seqs = len(text), [g.guide for g in guides]
    here = core.find_guide_sites(text, circular, seqs, nuclease.pam, max_mm, nuclease.pam_first)
    ref = as_reference(ref)
    same = ref is not None and ref.single and ref.name == name and ref.length == n
    # Every sequence of the reference (all chromosomes of a multi-record FASTA), one at a time, counts add up.
    far = None
    if ref is not None and not same:
        far = [[0] * (max_mm + 1) for _ in guides]
        for rname, rtext, rcircular in ref.sequences():
            if cancel is not None and cancel():
                raise Cancelled()
            if not ref.single and (rname, len(rtext)) == (name, n):
                continue   # the open sequence is one of the genome's: not compared with itself
            for total, hits in zip(far, core.find_guide_sites(rtext, rcircular, seqs, nuclease.pam, max_mm, nuclease.pam_first)):
                total[:] = [a + b for a, b in zip(total, _counts(hits, nuclease.length, max_mm))]
    out = []
    for g, h, f in zip(guides, here, far or [None] * len(guides)):
        out.append(Specificity(_counts(h, nuclease.length, max_mm, (g.start % n, g.strand) if n else None), f,
                               ref.name if far is not None else None))
    return out


def verdict(spec, g):
    """(level, text): ok / warn / bad from the counts. Another perfect site is bad, a near one (up to 2 mismatches) a
    warning, and a poor checklist score a warning."""
    reasons, level = [], "ok"

    def worse(new):
        nonlocal level
        level = new if ["ok", "warn", "bad"].index(new) > ["ok", "warn", "bad"].index(level) else level

    if spec is None:
        return "info", "specificity not checked"
    if spec.here[0]:
        worse("bad")
        reasons.append(f"{plural(spec.here[0], 'other perfect site')} in this sequence")
    elif sum(spec.here[1:3]):
        worse("warn")
        reasons.append(f"{plural(sum(spec.here[1:3]), 'similar site')} here (up to 2 mismatches)")
    if spec.genome is not None:
        if spec.genome[0] > 1:
            worse("bad")
            reasons.append(f"{spec.genome[0]} perfect sites in {spec.genome_name}")
        elif sum(spec.genome[1:3]):
            worse("warn")
            reasons.append(f"{plural(sum(spec.genome[1:3]), 'similar site')} in {spec.genome_name} (up to 2 mismatches)")
    if g.score < 70:
        worse("warn")
        reasons.extend(g.notes)
    return level, "; ".join(reasons) or "no other site up to 2 mismatches"


def cut_text(nuclease, g, n):
    """Where the nuclease cuts a guide, in words (bases numbered from 1; a cut is before the base named)."""
    if nuclease.pam_first:
        return (f"cuts the top strand before base {g.cut % n + 1:,} and the bottom strand before base {g.bottom % n + 1:,} "
                f"(a {abs(g.bottom - g.cut)} nt 5′ overhang, away from the PAM)")
    return f"cuts {nuclease.cut_from_pam} bp from the PAM, before base {g.cut % n + 1:,}"


# ---- the tool -------------------------------------------------------------------------------------------------------


class _Done(QObject):
    """Brings a background scan's results back to the GUI thread (the scan's version number comes first: a result of an older
    list is ignored)."""
    finished = Signal(int, object)        # the manual check: all the results, or the exception
    chunk = Signal(int, object, object)   # the automatic check: guide indices and their Specificity, as they come
    ended = Signal(int, object)           # the automatic check is over: "" (complete), "cancelled" or an error text


class _Item(QTableWidgetItem):
    """A cell that sorts by a key (a number, a tuple) kept in its UserRole data instead of by its text."""

    def __lt__(self, other):
        a, b = self.data(Qt.UserRole), other.data(Qt.UserRole)
        return a < b if a is not None and b is not None else super().__lt__(other)


def counts_key(counts):
    """Sort key of a column of site counts: unchecked first, then by perfect sites, then by the near ones."""
    return (-1, 0) if counts is None else (counts[0], sum(counts[1:]))


class Cas9Dialog(ToolDialog):
    def __init__(self, parent, editor):
        super().__init__(parent)
        self.win, self.editor, self.guides, self.specs, self.version = parent, editor, [], [], 0
        self.stop = {"event": threading.Event()}   # the automatic check's cancel flag (a dict: the closing hook below sees the current one)
        self.setWindowTitle("Spacer Finder (Cas9 / Cas12a)")
        self.resize(900, 520)
        self.nuclease = QComboBox()
        self.nuclease.addItems(list(NUCLEASES))
        self.region = QComboBox()
        self.region.addItems(["Selection (else the whole sequence)", "Whole sequence"])
        self.max_mm = QSpinBox(minimum=1, maximum=4, value=3, suffix=" mismatches")
        self.max_mm.setToolTip("Other sites with up to this many mismatches are counted")
        self.max_mm.valueChanged.connect(lambda v: self.max_mm.setSuffix(" mismatch" if v == 1 else " mismatches"))   # not "1 mismatches"
        find = QPushButton("Find guides", clicked=self.reload)
        top = QHBoxLayout()
        for w in (QLabel("Nuclease"), self.nuclease, QLabel("Look in"), self.region, QLabel("Count sites up to"), self.max_mm, find):
            top.addWidget(w)
        top.addStretch()
        self.genome_label = QLabel(wordWrap=True)   # (a long genome name wraps instead of widening the window)
        self.genome_label.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Preferred)
        options = QPushButton("Genome options…", clicked=self.genome_options,
                              toolTip="The background genome is shared with the Primers tool (Primers > Options)")
        genome_row = QHBoxLayout()
        genome_row.addWidget(self.genome_label, 1)
        genome_row.addWidget(options)
        self.table = QTableWidget(0, 9)
        self.table.setHorizontalHeaderLabels(["Position", "Strand", "PAM and guide 5′→3′", "Cuts after base", "GC", "Score", "Here", "Genome", "Verdict"])
        value_filter.install(self.table, ["Strand", "Here", "Genome", "Verdict"])   # right-click a header: filter by values
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        fit_columns(self.table, stretch=8)
        self.table.horizontalHeader().setSortIndicator(-1, Qt.AscendingOrder)   # unsorted until a header is clicked
        self.table.currentCellChanged.connect(lambda row, *_: self._select(self._index(row)))
        self.table.setContextMenuPolicy(Qt.CustomContextMenu)
        self.table.customContextMenuRequested.connect(self._menu)
        self.check = QPushButton("Check specificity", clicked=self.start_check,
                                 toolTip="Count the other places each guide matches, here and in the background genome")
        add = QPushButton("Add selected as feature", clicked=lambda: self.add())
        buttons = QHBoxLayout()
        buttons.addWidget(self.check)
        buttons.addWidget(add)
        buttons.addStretch()
        # progress and cancel of the automatic check of the first AUTO_CHECK guides
        self.progress = QProgressBar(textVisible=True, format="%v of %m guides checked")
        self.cancel = QPushButton("Cancel check", clicked=self.stop_auto, toolTip="Stop the automatic specificity check")
        busy = QHBoxLayout()
        busy.addWidget(self.progress, 1)
        busy.addWidget(self.cancel)
        busy.setContentsMargins(0, 0, 0, 0)
        self.busy = QWidget()
        self.busy.setLayout(busy)
        keep = self.busy.sizePolicy()
        keep.setRetainSizeWhenHidden(True)   # its room stays when it hides: the table above does not jump at each check
        self.busy.setSizePolicy(keep)
        self.busy.hide()
        # Plain text (an error's text is never read as HTML), wrapped, and two lines tall from the start: a long message neither
        # widens the window nor pushes the table up when it appears.
        self.status = QLabel(wordWrap=True, textFormat=Qt.PlainText, alignment=Qt.AlignLeft | Qt.AlignTop)
        self.status.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Preferred)
        self.status.setMinimumHeight(2 * self.status.fontMetrics().lineSpacing() + 2)
        layout = QVBoxLayout(self)
        layout.addLayout(top)
        layout.addLayout(genome_row)
        layout.addWidget(self.table, 1)
        layout.addWidget(self.busy)
        layout.addLayout(buttons)
        layout.addWidget(self.status)
        # the cut goes with the tool, and so does a check still running
        self.destroyed.connect(lambda _=None, ed=editor, stop=self.stop: (ed.set_cut_marks("cas9", []), stop["event"].set()))
        # (not `self.done`: that would hide QDialog.done(), which Escape and the window's close button go through)
        self.relay = _Done(self)
        self.relay.finished.connect(self._checked)
        self.relay.chunk.connect(self._auto_chunk)
        self.relay.ended.connect(self._auto_ended)
        no_auto_default(self)
        find.setDefault(True)   # Enter in a field finds the guides, never adds one as a feature
        self.nuclease.currentIndexChanged.connect(self.reload)
        self.region.currentIndexChanged.connect(self.reload)
        self.update_genome_label()
        self.reload()

    # -- the reference genome: the Primers tool's
    def update_genome_label(self):
        g = genome()
        self.genome_label.setText(f"Background genome: {genome_label_text(g)}, shared with the Primers tool" if g is not None
                                  else "Background genome: none set (Primers > Options), so only this sequence is checked")
        header = self.table.horizontalHeaderItem(7)
        if g is None:
            header.setText("Genome (no reference)")
            header.setToolTip("No background genome is set (Primers > Options): the first guides are not checked automatically, "
                              "and this column stays empty. “Check specificity” still counts the other sites in this sequence.")
        else:
            header.setText("Genome")
            header.setToolTip(f"Other sites in {genome_label_text(g)} by number of mismatches (0/1/2/…). The first {AUTO_CHECK} "
                              "guides of the list are counted automatically; “Check specificity” does the rest.")

    def genome_options(self):
        if OptionsDialog(self).exec():
            self.update_genome_label()
            if self.win is not None and hasattr(self.win, "primer_options_changed"):
                self.win.primer_options_changed()

    def refresh(self):
        self.reload()   # after edits, once typing pauses

    # -- finding guides
    def _stretch(self):
        """(offset, sequence, whole circular sequence?) the guides are looked for in, or an error text."""
        ed = self.editor
        if self.region.currentIndex() == 0 and ed.end > ed.start:
            a, b = ed.start, min(ed.end, ed.length)
            return a, ed.doc.sequence(a, b), False, ""
        if ed.length > MAX_REGION:
            return 0, "", False, f"This sequence is {ed.length:,} bp: select a region first (up to {MAX_REGION:,} bp)."
        return 0, ed.doc.sequence(0, None), ed.circular, ""

    def reload(self):
        self.version += 1
        self.stop_auto(quiet=True)
        self.busy.hide()
        self.check.setEnabled(True)
        self.editor.set_cut_marks("cas9", [])
        nuc = NUCLEASES[self.nuclease.currentText()]
        offset, seq, wrap, error = self._stretch()
        circular = self.editor.circular   # (up to 2 Mb: in a thread, behind a busy dialog if it takes a while)
        found = run_busy(self, "Finding guides…", lambda: candidates(seq, circular, nuc, offset, wrap)) if seq else []
        self.guides, self.specs = found[:MAX_CANDIDATES], [None] * min(len(found), MAX_CANDIDATES)
        self._fill()
        base = error or f"{plural(len(found), 'guide')} in {len(seq):,} bp" + (f"; the best {MAX_CANDIDATES} are listed" if len(found) > MAX_CANDIDATES else "")
        self.status.setText(base + ".")
        if self.guides and not error:
            self._auto_check(base)

    def _cells(self, g, spec, nuc, n):
        """[(text, sort key)] of the nine cells of a guide's row."""
        level, text = verdict(spec, g)
        here = "–" if spec is None else "/".join(map(str, spec.here))
        far = "–" if spec is None or spec.genome is None else "/".join(map(str, spec.genome))
        pam_guide = f"{g.pam} {g.guide}" if nuc.pam_first else f"{g.guide} {g.pam}"
        return [(f"{g.start % n + 1:,}–{(g.end - 1) % n + 1:,}", g.start % n), ("→" if g.strand > 0 else "←", g.strand), (pam_guide, pam_guide),
                (f"{g.cut % n:,}", g.cut % n), (f"{g.gc} %", g.gc), (str(g.score), g.score),
                (here, counts_key(None if spec is None else spec.here)),
                (far, counts_key(None if spec is None else spec.genome)),
                (f"{MARKS[level][0]} {text}", (["info", "ok", "warn", "bad"].index(level), text))]

    def _set_row(self, r, i):
        """Write guide `i` into table row `r` (sorting must be off while rows are written)."""
        nuc, n, g, spec = NUCLEASES[self.nuclease.currentText()], self.editor.length, self.guides[i], self.specs[i]
        level = verdict(spec, g)[0]
        for c, (v, key) in enumerate(self._cells(g, spec, nuc, n)):
            item = self.table.item(r, c)
            if item is None:
                item = _Item()
                self.table.setItem(r, c, item)
            item.setText(v)
            if c in (0, 3, 4, 5):
                right_aligned(item)
            item.setData(Qt.UserRole, key)
            item.setToolTip(f"Sites by number of mismatches: 0/1/2/…; {nuc.name} {cut_text(nuc, g, n)}" if c in (6, 7) else "; ".join(g.notes))
            if c == 8:
                item.setForeground(QColor(MARKS[level][1]))
        self.table.item(r, 0).setData(INDEX, i)

    def _fill(self):
        header = self.table.horizontalHeader()
        section, order = header.sortIndicatorSection(), header.sortIndicatorOrder()   # a clicked header's order is kept
        self.table.setSortingEnabled(False)
        self.table.setRowCount(0)
        self.table.setRowCount(len(self.guides))
        for r in range(len(self.guides)):
            self._set_row(r, r)
        header.setSortIndicator(section, order)
        self.table.setSortingEnabled(section >= 0)

    def _index(self, row):
        """The guide (index into self.guides) shown in table row `row`, or None."""
        item = self.table.item(row, 0) if row >= 0 else None
        i = item.data(INDEX) if item is not None else None
        return i if i is not None and 0 <= i < len(self.guides) else None

    def _row_of(self, i):
        return next((r for r in range(self.table.rowCount()) if self._index(r) == i), None)

    def _select(self, i):
        if i is not None:
            g = self.guides[i]
            self.editor.select(g.start, g.end)
            self.editor.reveal.emit(g.start)
            nuc = NUCLEASES[self.nuclease.currentText()]
            n = self.editor.length
            ed = self.editor
            tip = lambda g=g, nuc=nuc: cas9_mark_tip(nuc, g.start % n, g.strand, n, ed.doc.sequence, f"gRNA {g.guide[:8]}")
            self.editor.set_cut_marks("cas9", [(g.cut, g.bottom, nuc.name, tip)])   # Cas9: the same place on both strands

    # -- specificity of the first guides, by itself, in a background thread (a genome takes a moment)
    def _auto_check(self, base):
        """Count the other sites of the first AUTO_CHECK guides, CHUNK at a time, so the table fills as results come in and a
        cancel is seen within a chunk. Needs a background genome (Primers > Options); without one, nothing is done (the Genome
        column's header says so)."""
        ref = genome()
        if ref is None:
            self.status.setText(base + ". No background genome is set (Primers > Options), so specificity is not checked "
                                "automatically: press “Check specificity” to count the other sites in this sequence.")
            return
        nuc, mm, version = NUCLEASES[self.nuclease.currentText()], self.max_mm.value(), self.version
        indices = list(range(min(AUTO_CHECK, len(self.guides))))
        guides = [self.guides[i] for i in indices]
        doc = self.editor.doc
        text, circular, name = doc.sequence(0, None), doc.circular, doc.name   # copied here, on the GUI thread
        stop = self.stop["event"] = threading.Event()
        self.progress.setRange(0, len(indices))
        self.progress.setValue(0)
        self.busy.show()
        self.check.setEnabled(False)
        self.auto_base = base
        self.status.setText(f"{base}. Checking the first {plural(len(indices), 'guide')} against this sequence and {ref.name}…")

        def emit(signal, *args):
            try:
                signal.emit(*args)
            except RuntimeError:
                pass   # the dialog is gone

        def work():
            try:
                for k in range(0, len(guides), CHUNK):
                    if stop.is_set():
                        raise Cancelled()
                    emit(self.relay.chunk, version, indices[k:k + CHUNK], scan(text, circular, name, guides[k:k + CHUNK], nuc, mm, ref, stop.is_set))
                emit(self.relay.ended, version, "")
            except Cancelled:
                emit(self.relay.ended, version, "cancelled")
            except Exception as e:   # shown in the status line, never lost in a thread
                emit(self.relay.ended, version, str(e) or type(e).__name__)

        threading.Thread(target=work, daemon=True).start()

    def stop_auto(self, quiet=False):
        """Cancel the automatic check (the button, or a new list of guides)."""
        self.stop["event"].set()
        if not quiet:
            self.status.setText("Cancelling the specificity check…")

    def _auto_chunk(self, version, indices, results):
        if version != self.version:
            return
        sorting = self.table.isSortingEnabled()
        self.table.setSortingEnabled(False)   # (writing a cell of the sorted column would move the row under us)
        for i, spec in zip(indices, results):
            self.specs[i] = spec
            row = self._row_of(i)
            if row is not None:
                self._set_row(row, i)
        self.table.setSortingEnabled(sorting)
        self.progress.setValue(self.progress.value() + len(indices))

    def _auto_ended(self, version, how):
        if version != self.version:
            return
        self.busy.hide()
        self.check.setEnabled(True)
        done = sum(s is not None for s in self.specs)
        total = len(self.guides)
        if how and how != "cancelled":
            self.status.setText(f"{self.auto_base}. The specificity check failed: {how}")
        else:
            self.status.setText(f"{self.auto_base}. {done} of {total} checked" + (" (cancelled)" if how else "")
                                + (f"; press “Check specificity” to count the other {total - done}." if done < total else "."))

    # -- specificity of all of them, on request
    def start_check(self):
        if not self.guides:
            return
        self.stop_auto(quiet=True)
        nuc, mm, version = NUCLEASES[self.nuclease.currentText()], self.max_mm.value(), self.version
        guides, doc, ref = list(self.guides), self.editor.doc, genome()
        self.busy.hide()
        self.check.setEnabled(False)
        self.status.setText(f"Checking {plural(len(guides), 'guide')} against this sequence" + (f" and {ref.name}…" if ref is not None else "…"))
        # shortcut: copies the sequence on the GUI thread, scans on another; an edit meanwhile only makes the result stale.
        text, circular, name = doc.sequence(0, None), doc.circular, doc.name

        def work():
            try:
                result = scan(text, circular, name, guides, nuc, mm, ref)
            except Exception as e:   # shown in the status line, never lost in a thread
                result = e
            try:
                self.relay.finished.emit(version, result)
            except RuntimeError:
                pass   # the dialog is gone

        threading.Thread(target=work, daemon=True).start()

    def _checked(self, version, result):
        self.check.setEnabled(True)
        if version != self.version:
            return   # the sequence or the options changed meanwhile
        if isinstance(result, Exception):
            self.status.setText(f"The specificity check failed: {str(result) or type(result).__name__}")
            return
        order = sorted(range(len(self.guides)), key=lambda i: (["ok", "warn", "bad"].index(verdict(result[i], self.guides[i])[0]),
                                                                -self.guides[i].score, self.guides[i].start))
        self.guides, self.specs = [self.guides[i] for i in order], [result[i] for i in order]
        self._fill()
        good = sum(verdict(s, g)[0] == "ok" for s, g in zip(self.specs, self.guides))
        self.status.setText(f"{good} of {plural(len(self.guides), 'guide')} {'has' if good == 1 else 'have'} no other site up to 2 mismatches. "
                            "Sites are counted as 0/1/2/… mismatches; best guides first.")

    # -- using a guide
    def selected_guides(self):
        """Indices (into self.guides) of the selected rows, in table order; the current row if nothing is selected."""
        rows = sorted({ix.row() for ix in self.table.selectedIndexes()}) or [self.table.currentRow()]
        return [i for i in (self._index(r) for r in rows) if i is not None]

    def feature_for(self, g):
        """The misc_feature (kind, location, qualifiers) that saves guide `g`: the button and the right-click menu both use this."""
        nuc = NUCLEASES[self.nuclease.currentText()]
        n = self.editor.length
        location = site_location(g.start % n, g.start % n + (g.end - g.start), g.strand, n)
        return ("misc_feature", location, [
            ("label", f"gRNA {g.guide[:8]}"),
            ("note", f"{nuc.name} guide {g.guide}, PAM {g.pam}; {cut_text(nuc, g, n)}; checklist score {g.score}"),
            (CAS9_KEY, cas9_qualifier(nuc))])   # the last one lets this program redraw protospacer, PAM and cut (see cut_sites.py)

    def add(self, indices=None):
        """Save guides as features: `indices` (into self.guides), else the selected rows."""
        indices = self.selected_guides() if indices is None else indices
        if not indices:
            self.status.setText("Select a guide in the table first.")
            return
        self.editor.add_features([self.feature_for(self.guides[i]) for i in indices])
        self.status.setText(f"Added guide {self.guides[indices[0]].guide} as a feature." if len(indices) == 1
                            else f"Added {len(indices)} guides as features.")

    def build_menu(self, pos):
        """The right-click menu of the table at viewport position `pos` (None over no row); the row under the pointer joins the
        selection if it was not in it, as in any list."""
        item = self.table.itemAt(pos)
        if item is None:
            return None
        if not item.isSelected():
            self.table.selectRow(item.row())
        chosen = self.selected_guides()
        menu = QMenu(self.table)
        menu.addAction("Add as feature" if len(chosen) == 1 else f"Add selected as features ({len(chosen)})",
                       lambda chosen=chosen: self.add(chosen))
        return menu

    def _menu(self, pos):
        menu = self.build_menu(pos)
        if menu is not None:
            menu.exec(self.table.viewport().mapToGlobal(pos))
