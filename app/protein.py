"""The Protein view: the amino acids of the selected CDS, or of the translation of the selected DNA, numbered in blocks of ten, with
the protein's length, mass, isoelectric point and charge, and the features that lie inside the CDS.

Select a CDS (click it) or a stretch of DNA (a frame is chosen beside it) and the view follows. Select letters of the protein to select
their codons in the sequence. Under each line of letters, tracks (protein_tracks.py) show the residues of the selection in the sequence,
the secondary structure and confidence of the 3D model, hydrophobicity, size, charge and property; the buttons beside the letters
switch them on and off. A structure got with "3D structure…" shows beside the letters (structure_view.py) for as long as that protein
is shown: the selection is highlighted on it too, and clicking a residue on it selects its codons. While a CDS is shown, a selection
inside it (a feature, a stretch) keeps it shown.

Mass: sum of the average residue masses + water. Isoelectric point and charge: Henderson-Hasselbalch with the pKa values of EMBOSS's
iep (Rice et al. 2000, see Help > Citations). Everything is computed up to the first stop codon; the numbers are estimates of the
free chain (no signal peptide cut, no modification).
"""
from html import escape

from PySide6.QtCore import QTimer, Qt
from PySide6.QtGui import QPainter, QPalette, QTextCursor, QTextDocument
from PySide6.QtWidgets import (
    QComboBox, QHBoxLayout, QLabel, QListWidget, QListWidgetItem, QPushButton, QSizePolicy, QSplitter, QVBoxLayout, QWidget,
)

from editor import HIDDEN_KINDS, mono_font
import genetic_code
import structure
import tool_settings
from sequence_view import cds_codons, codons
from structure_dialog import StructureDialog
from protein_tracks import TrackButtons, TrackText
from structure_view import StructureView

PER_LINE, BLOCK, MARGIN = 50, 10, 8          # amino acids per line, per block, and the width of the number column
MAX_BASES = 300_000                          # a longer selection is not translated here
FRAMES = ["+1", "+2", "+3", "−1", "−2", "−3"]
RESIDUE_MASS = {   # average mass of an amino acid in a chain (Da), as on ExPASy
    "G": 57.0519, "A": 71.0788, "S": 87.0782, "P": 97.1167, "V": 99.1326, "T": 101.1051, "C": 103.1388, "L": 113.1594, "I": 113.1594,
    "N": 114.1038, "D": 115.0886, "Q": 128.1307, "K": 128.1741, "E": 129.1155, "M": 131.1926, "H": 137.1411, "F": 147.1766,
    "R": 156.1875, "Y": 163.1760, "W": 186.2132,
}
WATER = 18.0153
PKA_POSITIVE = {"K": 10.8, "R": 12.5, "H": 6.5}      # side chains
PKA_NEGATIVE = {"D": 3.9, "E": 4.1, "C": 8.5, "Y": 10.1}
PKA_N_TERM, PKA_C_TERM = 8.6, 3.6


def until_stop(protein):
    return protein.split("*")[0]


def mass(protein):
    """Average mass in Da of the chain `protein` (letters it does not know, such as X, are left out)."""
    return WATER + sum(RESIDUE_MASS.get(a, 0.0) for a in protein)


def charge(protein, ph):
    """Net charge of the free chain at `ph`."""
    pos = 1 / (1 + 10 ** (ph - PKA_N_TERM)) + sum(1 / (1 + 10 ** (ph - PKA_POSITIVE[a])) for a in protein if a in PKA_POSITIVE)
    neg = 1 / (1 + 10 ** (PKA_C_TERM - ph)) + sum(1 / (1 + 10 ** (PKA_NEGATIVE[a] - ph)) for a in protein if a in PKA_NEGATIVE)
    return pos - neg


def isoelectric_point(protein):
    """The pH where the net charge is zero (bisection: the charge falls as the pH rises)."""
    low, high = 0.0, 14.0
    for _ in range(40):
        mid = (low + high) / 2
        low, high = (mid, high) if charge(protein, mid) > 0 else (low, mid)
    return (low + high) / 2


def index_at(offset, block_number):
    """The amino acid (0-based) at `offset` characters into the line `block_number` of the text."""
    col = max(0, offset - MARGIN)
    return block_number * PER_LINE + min(col - col // (BLOCK + 1), PER_LINE - 1)


def offset_of(i):
    """(line, characters into the line) of amino acid `i` in the text."""
    line, col = divmod(i, PER_LINE)
    return line, MARGIN + col + col // BLOCK


def layout(protein):
    """The text: numbered lines of PER_LINE amino acids in blocks of BLOCK."""
    lines = []
    for start in range(0, len(protein), PER_LINE):
        chunk = protein[start:start + PER_LINE]
        lines.append(f"{start + 1:>{MARGIN - 2}}  " + " ".join(chunk[i:i + BLOCK] for i in range(0, len(chunk), BLOCK)))
    return "\n".join(lines)


class TitleLabel(QLabel):
    """The view's title on one line beside the buttons: a title that wrapped onto two or three lines made the text under it jump
    as the selection changed. Too long for its room, it is written plain and cut short with "…", the whole title in its tooltip.
    text() is the whole title."""

    def __init__(self):
        super().__init__()
        self.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Preferred)
        self.setMinimumWidth(0)

    def _document(self):
        doc = QTextDocument()
        doc.setDocumentMargin(0)
        doc.setDefaultFont(self.font())
        doc.setHtml(self.text())
        return doc

    def _cut(self):
        """The plain title when it does not fit, else None."""
        doc = self._document()
        return doc.toPlainText() if doc.idealWidth() > self.contentsRect().width() else None

    def setText(self, text):
        super().setText(text)
        self.setToolTip(self._cut() or "")

    def resizeEvent(self, event):
        super().resizeEvent(event)
        self.setToolTip(self._cut() or "")

    def paintEvent(self, event):
        plain = self._cut()
        if plain is None:
            return super().paintEvent(event)
        p = QPainter(self)
        r = self.contentsRect()
        p.setPen(self.palette().color(QPalette.WindowText))
        p.drawText(r, Qt.AlignLeft | Qt.AlignVCenter, self.fontMetrics().elidedText(plain, Qt.ElideRight, r.width()))
        p.end()


class ProteinView(QWidget):
    def __init__(self, editor):
        super().__init__()
        self.editor, self.src, self.own, self.dialogs = editor, None, False, []
        self.structures = {}  # amino acids -> the Structure got for them this session (the few last ones)
        self.anchor = None    # the residue last clicked on the 3D view (Shift+click selects from it)
        self.protein = ""     # the amino acids in the text
        self.shown = None     # what the text shows (see _show)
        self.places = {}      # what was shown before -> the place the reader was in it (see _set_text)
        self.title = TitleLabel()
        self.frame = QComboBox()
        self.frame.addItems(FRAMES)
        self.frame.setToolTip("Reading frame of the selected DNA (a selected CDS has its own)")
        self.frame.currentIndexChanged.connect(lambda _: self.refresh(True))
        self.text = TrackText(PER_LINE, offset_of, index_at)
        self.text.setFont(mono_font(max(9, round(self.font().pointSizeF()))))   # (the system's fixed font is a thin Courier on Windows)
        self.text.selectionChanged.connect(self._text_selected)
        self.props = QLabel(wordWrap=True, textInteractionFlags=Qt.TextSelectableByMouse)
        self.features = QListWidget()
        self.features.setMaximumHeight(130)
        self.features.setToolTip("Features inside this CDS: click one to select it in the sequence")
        self.features.itemClicked.connect(lambda item: self.editor.select_feature(item.data(Qt.UserRole)))
        self.features.hide()   # shown when the CDS has features inside (see _features)
        self.fold = QPushButton("3D structure…", clicked=self.show_structure, enabled=False,
                                toolTip="Look up or predict the 3D structure of this protein")
        self.tracks = TrackButtons(self._update_tracks)
        letters = QWidget()   # the letters and, on their right, the buttons of the tracks
        row = QHBoxLayout(letters)
        row.setContentsMargins(0, 0, 0, 0)
        row.addWidget(self.text, 1)
        row.addWidget(self.tracks)
        top = QHBoxLayout()
        top.addWidget(self.title, 1)
        top.addWidget(self.fold)
        top.addWidget(QLabel("Frame"))
        top.addWidget(self.frame)
        self.view3d = StructureView()
        self.view3d.picked.connect(self._picked)
        self.view3d.hide()
        self.split = QSplitter(Qt.Vertical, childrenCollapsible=False)   # side by side in a wide pane: see resizeEvent
        self.split.addWidget(letters)
        self.split.addWidget(self.view3d)
        self.split.setStretchFactor(1, 1)
        layout_ = QVBoxLayout(self)
        layout_.addLayout(top)
        layout_.addWidget(self.split, 1)
        layout_.addWidget(self.props)
        layout_.addWidget(self.features)
        # Dragging a selection fires often: look again once it pauses.
        self.timer = QTimer(self, singleShot=True, interval=120, timeout=self.refresh)
        editor.selection_changed.connect(self._selection_changed)
        tool_settings.notifier.changed.connect(self._setting_changed)   # a bound method, like _edited
        editor.changed.connect(self._edited)   # a bound method, not a lambda: Qt drops it when this pane is closed
        self.refresh(True)

    def resizeEvent(self, event):
        super().resizeEvent(event)
        wide = Qt.Horizontal if self.width() > 1.5 * self.height() else Qt.Vertical
        if wide != self.split.orientation():
            self.split.setOrientation(wide)
            if wide == Qt.Horizontal:   # the letters get a whole line (number, 50 amino acids, spaces), the 3D view the rest
                line = (self.text.fontMetrics().horizontalAdvance("0" * (MARGIN + PER_LINE + PER_LINE // BLOCK)) + 30
                        + self.tracks.sizeHint().width())
                self.split.setSizes([line, max(self.width() - line, 1)])

    # ---- what to show ---------------------------------------------------------------------------------------------------------------
    def _edited(self):
        self.refresh(True)

    def _setting_changed(self, tool):
        if tool == genetic_code.TOOL:   # another genetic code: translate again
            self.refresh(True)

    def _selection_changed(self):
        if not self.own:
            self.timer.start()

    def _cds_source(self, f):
        ed = self.editor
        plan = genetic_code.cds_plan(ed.doc, f)     # /transl_table, /codon_start and /transl_except of this CDS
        cods = cds_codons(ed.doc, f)
        extent = (min(s for s, _ in f.spans), max(e for _, e in f.spans))
        return dict(kind="CDS", feature=f.index, name=f.label or "CDS", accession=structure.accession_of(f.qualifiers),
                    links=structure.pdb_links(f.qualifiers), cods=cods, extent=extent, version=ed.version,
                    title=f"<b>{escape(f.label or 'CDS')}</b> · CDS {extent[0] + 1:,}–{extent[1]:,}, {'reverse' if f.strand < 0 else 'forward'} strand"
                          + genetic_code.table_note(plan.table) + (" · translation not reproducible" if plan.why else ""))

    def _selection_source(self):
        ed = self.editor
        if not 0 < ed.end - ed.start <= MAX_BASES:
            return None
        k = self.frame.currentIndex()
        strand = 1 if k < 3 else -1
        table = genetic_code.table_id()
        cods = codons(ed.doc, ed.selected_spans(), strand, k % 3, table)
        return dict(kind="selection", feature=None, name=f"translation {ed.start + 1}", accession="", links=[], cods=cods, extent=(ed.start, ed.end), version=ed.version,
                    title=f"Translation of {ed.start + 1:,}–{ed.end - (ed.length if ed.wraps else 0):,} in frame {FRAMES[k]}" + genetic_code.table_note(table))

    def refresh(self, force=False):
        """Choose what to show for the current selection (keeping a CDS shown while the selection is inside it) and show it."""
        ed, src = self.editor, self.src
        f = ed.by_index.get(ed.current_feature) if ed.current_feature is not None else None
        pinned = (not force and src and src["kind"] == "CDS" and src["version"] == ed.version and
                  src["extent"][0] <= ed.start and ed.end <= src["extent"][1] and not (f and f.kind == "CDS" and f.index != src["feature"]))
        if pinned:
            self._highlight()
            return
        if f is not None and f.kind == "CDS":
            src = self._cds_source(f)
        elif ed.end - ed.start > MAX_BASES:
            src = None
            self._empty(f"The selection is longer than {MAX_BASES:,} bases: select a CDS or a shorter stretch.")
        else:
            src = self._selection_source()
            if src is None:
                self._empty("Select a CDS (click it) or a stretch of DNA to see its protein.")
        self.src = src
        self.frame.setEnabled(src is None or src["kind"] == "selection")
        if src is not None:
            self._show(src)

    def _empty(self, message):
        self._leave()
        self.src = self.shown = None
        self.protein = ""
        self.title.setText(message)
        self.fold.setEnabled(False)
        self.view3d.hide()
        self.text.protein = ""
        self.text.setPlainText("")
        self._update_tracks()
        self.props.setText("")
        self.features.clear()
        self.features.hide()

    def _show(self, src):
        protein = "".join(aa for _, aa, _, _ in src["cods"])
        # what is shown: a CDS, or a stretch read in a frame (another stretch is another text, read from its top)
        key = ("CDS", src["feature"]) if src["kind"] == "CDS" else ("selection", src["extent"], self.frame.currentIndex())
        self.protein = self.text.protein = protein   # (before the text: its tracks are read from it)
        self._set_text(key, layout(protein))
        chain = until_stop(protein)
        # A stop codon is not an amino acid: count the residues before it, like the properties line below.
        self.title.setText(f"{src['title']} · {len(chain):,} aa" + (" + stop" if len(chain) < len(protein) else ""))
        self.fold.setEnabled(bool(chain))
        self._show_3d(chain)
        if not chain:
            self.props.setText("No amino acids before the first stop codon.")
        else:
            stop = f" · stop codon at aa {len(chain) + 1:,}" if "*" in protein[:-1] else ""
            unknown = f" · {protein.count('X'):,} unknown (X) left out of the mass" if "X" in chain else ""
            self.props.setText(f"<b>{len(chain):,} aa</b> · {mass(chain) / 1000:,.2f} kDa · pI {isoelectric_point(chain):.2f} · "
                               f"charge {charge(chain, 7.0):+.1f} at pH 7{stop}{unknown}")
        self._features(src)
        self._highlight()

    def _place(self):
        """Where the reader is in the text: the scroll bars and the letters selected."""
        cursor = self.text.textCursor()
        return [self.text.verticalScrollBar().value(), self.text.horizontalScrollBar().value()], cursor.anchor(), cursor.position()

    def _leave(self):
        """Note the place in what is shown before something else replaces it (see _set_text)."""
        if self.shown is not None:
            self.places.pop(self.shown, None)
            self.places[self.shown] = self._place()
            while len(self.places) > 8:   # the few last ones are enough
                self.places.pop(next(iter(self.places)))

    def _set_text(self, key, text):
        """Show `text` for `key` (see _show). The same protein again (an edit elsewhere in the sequence shows it anew, often
        after a moment of showing something else) comes back at the place the reader was: same scroll position, same letters
        selected. It used to jump back to the top and drop the selection. The same text is not set again."""
        if key == self.shown and text == self.text.toPlainText():
            return
        if key != self.shown:
            self._leave()
            place = self.places.get(key)
        else:
            place = self._place()
        self.shown = key
        self.text.blockSignals(True)   # (putting the selection back is not the user selecting letters: see _text_selected)
        try:
            self.text.setPlainText(text)
            if place is not None:
                scroll, anchor, position = place
                cursor = QTextCursor(self.text.document())
                cursor.setPosition(min(anchor, len(text)))
                cursor.setPosition(min(position, len(text)), QTextCursor.KeepAnchor)
                self.text.setTextCursor(cursor)
        finally:
            self.text.blockSignals(False)
        if place is not None:
            for bar, value in zip((self.text.verticalScrollBar(), self.text.horizontalScrollBar()), place[0]):
                bar.setValue(value)

    def show_structure(self):
        """The 3D structure window for the protein shown (it stays open while you go on selecting)."""
        chain = until_stop("".join(aa for _, aa, _, _ in self.src["cods"]))
        dialog = StructureDialog(self.window(), self.src["name"], chain, self.src["accession"], self.src["links"])
        dialog.found.connect(self._got_structure)
        self.dialogs = [d for d in self.dialogs if d.isVisible()] + [dialog]
        dialog.show()

    def _got_structure(self, chain, found):
        self.structures.pop(chain, None)
        self.structures[chain] = found
        while len(self.structures) > 5:   # a model is a few hundred kB of text
            self.structures.pop(next(iter(self.structures)))
        self.view3d.protein = None        # (a new structure of the same protein replaces the one shown)
        if self.src is not None:
            self._show_3d(until_stop("".join(aa for _, aa, _, _ in self.src["cods"])))
            self._highlight()

    def _show_3d(self, chain):
        """The 3D view of the protein shown, when a structure was got for it; hidden otherwise."""
        found = self.structures.get(chain)
        if found is None:
            self.view3d.hide()
        else:
            if self.view3d.protein != chain:
                self.view3d.set_structure(found, chain)
                self.anchor = None
            self.view3d.show()
        self._update_tracks()

    def _update_tracks(self):
        """The tracks under the letters: the model's (structure, confidence) when this protein has one, and those switched on."""
        canvas, text = self.view3d.canvas, self.text
        chain = until_stop(self.protein)
        has_model = bool(chain) and self.view3d.protein == chain and chain in self.structures
        res = canvas.residues
        text.model = {i: (res[ks[0]].ss, res[ks[0]].b) for i, ks in canvas.to_model.items()} if has_model else {}
        text.plddt = has_model and canvas.plddt
        self.tracks.set_model(has_model)
        on = [n for n in self.tracks.on() if has_model or n not in ("Structure", "Confidence")]
        if not text.plddt and "Confidence" in on:   # an experimental structure has no pLDDT
            on.remove("Confidence")
        text.set_tracks(["Selection", *on])

    def _picked(self, i, extend):
        """A residue clicked on the 3D view: select its codons (with Shift, from the residue clicked before)."""
        lo, hi = (min(i, self.anchor), max(i, self.anchor)) if extend and self.anchor is not None else (i, i)
        if not extend:
            self.anchor = i
        self._select_residues(lo, hi)

    def _features(self, src):
        """The features inside the CDS, each with the amino acids it covers."""
        self.features.clear()
        if src["kind"] != "CDS":
            self.features.hide()
            return
        cods, ed = src["cods"], self.editor
        for f in ed.features:
            if f.index == src["feature"] or f.kind in HIDDEN_KINDS or f.kind in ("CDS", "gene"):
                continue
            lo, hi = min(s for s, _ in f.spans), max(e for _, e in f.spans)
            covered = [i for i, c in enumerate(cods) if c[2] < hi and c[3] > lo]
            if covered:
                item = QListWidgetItem(f"{f.label or f.kind} ({f.kind}) · aa {covered[0] + 1:,}–{covered[-1] + 1:,}")
                item.setData(Qt.UserRole, f.index)
                self.features.addItem(item)
        self.features.setVisible(self.features.count() > 0)   # an empty list is not left taking room

    def _highlight(self):
        """Mark the residues of the selection in the sequence: on the Selection track (the letters stay plain) and on the 3D view."""
        ed, src = self.editor, self.src
        hit = [i for i, c in enumerate(src["cods"]) if c[2] < ed.end and c[3] > ed.start] if src and ed.end > ed.start else []
        if set(hit) != self.text.marked:
            self.text.marked = set(hit)
            self.text.viewport().update()
        self.view3d.set_highlight(hit)

    # ---- selecting residues selects their codons ------------------------------------------------------------------------------
    def _text_selected(self):
        cur = self.text.textCursor()
        if not cur.hasSelection() or self.src is None:
            return
        doc = self.text.document()
        a, b = cur.selectionStart(), cur.selectionEnd() - 1
        ia = index_at(a - doc.findBlock(a).position(), doc.findBlock(a).blockNumber())
        ib = index_at(b - doc.findBlock(b).position(), doc.findBlock(b).blockNumber())
        cods = self.src["cods"]
        self._select_residues(min(ia, len(cods) - 1), min(max(ib, ia), len(cods) - 1))

    def _select_residues(self, ia, ib):
        """Select the codons of the residues `ia` to `ib` in the sequence (the CDS stays shown)."""
        part = self.src["cods"][ia:ib + 1]
        if part:
            self.own = True
            try:
                self.editor.select(min(c[2] for c in part), max(c[3] for c in part))
                self.editor.reveal.emit(part[0][2])
            finally:
                self.own = False
            self._highlight()

