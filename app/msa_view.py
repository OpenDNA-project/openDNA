"""The Multiple Alignment view: several sequences compared at once (clones against their reference, homologues of a gene, a part across
constructs). One row per sequence under its name, the columns coloured by how conserved they are, a consensus row, and beside each
name the percentage identity of that sequence to the reference (or to the consensus).

Which sequences: the records of the open file (a multi-record GenBank / FASTA file, read as in add_record.py), files you pick, or
the open sequence plus the bases selected in it. How they are aligned: MAFFT or Clustal Omega when installed (the recipes of
tools.toml, run by runner.py with progress and Cancel), else the built-in progressive aligner of msa.py.

A click on a column puts the open sequence's cursor at that column's position and shows, for every sequence, the position it has
there; dragging selects the columns' bases in the open sequence. Selecting in the other views marks the columns in this one.
"""
import tempfile
from pathlib import Path

from PySide6.QtCore import QRect, Qt, Signal
from PySide6.QtGui import QColor, QPainter
from PySide6.QtWidgets import (
    QAbstractScrollArea, QComboBox, QDialog, QFileDialog, QHBoxLayout, QLabel, QMessageBox, QPushButton, QVBoxLayout, QWidget,
)

import msa
import runner
from add_record import source_records
from editor import mono_font
from uikit import run_busy

SOURCES = ["Records of this file", "Files…", "Open sequence + selection"]
BUILT_IN = "Built-in (progressive)"
AUTO_LIMIT = 30_000   # bp in all: a file of records this small is aligned when the view opens, a bigger one waits for Align
NAME_MAX = 220        # px: widest name column
IDENTITY_W = 64       # px: the identity column
CONSERVED = ((1.0, QColor(46, 160, 67, 190)), (0.75, QColor(150, 190, 60, 170)), (0.5, QColor(235, 170, 60, 150)))


def shade(conservation):
    """The background of a residue in a column of that conservation (None: not conserved enough to colour)."""
    return next((c for limit, c in CONSERVED if conservation >= limit), None)


class MsaCanvas(QAbstractScrollArea):
    """Draws the alignment: a frozen left part (names, identities), a ruler on top, the rows, then the consensus row."""
    clicked = Signal(int)          # a column was clicked
    dragged = Signal(int, int)     # columns first..last were selected by dragging

    def __init__(self):
        super().__init__()
        self.alignment, self.marked, self.anchor = None, None, None   # marked: (first, last) columns highlighted
        self.consensus, self.cons = "", []
        self.setFont(mono_font(max(9, round(self.font().pointSizeF()))))
        self.viewport().setMouseTracking(False)

    def set_alignment(self, alignment):
        self.alignment, self.marked = alignment, None
        if alignment is not None:
            self.consensus, self.cons = alignment.consensus(), alignment.conservation()
        self._layout()
        self.viewport().update()

    # geometry: one extra row on top (the ruler) and one at the bottom (the consensus)
    def _metrics(self):
        fm = self.fontMetrics()
        return fm.horizontalAdvance("A"), fm.height() + 4

    def _left(self):
        fm = self.fontMetrics()
        names = max((fm.horizontalAdvance(n) for n in self.alignment.names), default=0) if self.alignment else 0
        return min(NAME_MAX, names + 12) + IDENTITY_W

    def _layout(self):
        cw, rh = self._metrics()
        a = self.alignment
        if a is None:
            return
        view = self.viewport()
        self.horizontalScrollBar().setRange(0, max(0, a.width * cw - (view.width() - self._left())))
        self.horizontalScrollBar().setPageStep(max(1, view.width() - self._left()))
        self.horizontalScrollBar().setSingleStep(cw * 5)
        self.verticalScrollBar().setRange(0, max(0, (len(a.rows) + 2) * rh - view.height()))
        self.verticalScrollBar().setPageStep(view.height())
        self.verticalScrollBar().setSingleStep(rh)

    def resizeEvent(self, event):
        super().resizeEvent(event)
        self._layout()

    def column_at(self, x):
        cw, _ = self._metrics()
        return max(0, min(self.alignment.width - 1, (x - self._left() + self.horizontalScrollBar().value()) // cw))

    def ensure_visible(self, col):
        cw, _ = self._metrics()
        bar, room = self.horizontalScrollBar(), self.viewport().width() - self._left()
        if not bar.value() <= col * cw <= bar.value() + room - cw:
            bar.setValue(col * cw - room // 2)

    def mark(self, first, last):
        self.marked = None if first is None else (min(first, last), max(first, last))
        self.viewport().update()

    def mousePressEvent(self, event):
        if self.alignment is None or event.button() != Qt.LeftButton or event.position().x() < self._left():
            return
        self.anchor = self.column_at(int(event.position().x()))
        self.mark(self.anchor, self.anchor)
        self.clicked.emit(self.anchor)

    def mouseMoveEvent(self, event):
        if self.anchor is not None and event.buttons() & Qt.LeftButton:
            col = self.column_at(int(event.position().x()))
            self.mark(self.anchor, col)
            self.dragged.emit(min(self.anchor, col), max(self.anchor, col))

    def mouseReleaseEvent(self, event):
        self.anchor = None

    def paintEvent(self, event):
        a = self.alignment
        p = QPainter(self.viewport())
        p.setFont(self.font())
        pal = self.palette()
        p.fillRect(self.viewport().rect(), pal.base())
        if a is None:
            return
        cw, rh = self._metrics()
        left, h, y0 = self._left(), self.horizontalScrollBar().value(), -self.verticalScrollBar().value()
        first, last = h // cw, min(a.width, (h + self.viewport().width() - left) // cw + 1)
        text, ascent = pal.text().color(), self.fontMetrics().ascent() + 2
        identities = a.identities()
        for r in range(len(a.rows) + 2):   # row 0 is the ruler, the last one the consensus
            y = y0 + r * rh
            if y + rh < 0 or y > self.viewport().height():
                continue
            seq = self.consensus if r == len(a.rows) + 1 else (a.rows[r - 1] if 1 <= r <= len(a.rows) else None)
            for c in range(first, last):
                x = left + c * cw - h
                if r == 0:
                    if (c + 1) % 10 == 0:
                        p.setPen(text)
                        p.drawText(x - 3 * cw, y, 4 * cw, rh, Qt.AlignRight | Qt.AlignVCenter, str(c + 1))
                    continue
                color = shade(self.cons[c]) if seq[c] != msa.GAP else None
                if color is not None:
                    p.fillRect(x, y, cw, rh, color)
                p.setPen(text if seq[c] != msa.GAP else QColor(150, 150, 150))
                p.drawText(x, y + ascent, seq[c])
            if self.marked and r:   # the clicked / selected columns: a frame down every row
                x0, x1 = left + self.marked[0] * cw - h, left + (self.marked[1] + 1) * cw - h
                p.setPen(pal.highlight().color())
                p.drawRect(QRect(x0, y, x1 - x0, rh))
        # the frozen part is drawn last, over the letters that scroll under it
        p.fillRect(0, 0, left, self.viewport().height(), pal.base())
        for r in range(1, len(a.rows) + 2):
            y = y0 + r * rh
            if r == len(a.rows) + 1:
                label, pct = "Consensus", ""
            else:
                label, pct = a.names[r - 1], f"{identities[r - 1]:.1f}%" + ("" if a.reference != r - 1 else " ref")
            p.setPen(text)
            p.drawText(QRect(4, y, left - IDENTITY_W - 8, rh), Qt.AlignLeft | Qt.AlignVCenter,
                       self.fontMetrics().elidedText(label, Qt.ElideRight, left - IDENTITY_W - 8))
            p.drawText(QRect(left - IDENTITY_W, y, IDENTITY_W - 4, rh), Qt.AlignRight | Qt.AlignVCenter, pct)
        p.drawText(QRect(4, 0, left - 8, rh), Qt.AlignLeft | Qt.AlignVCenter,
                   "identity to " + ("consensus" if a.reference is None else "reference"))


class MsaView(QWidget):
    def __init__(self, editor):
        super().__init__()
        self.editor, self.alignment, self.sources, self.own = editor, None, [], False
        self.source = QComboBox()
        self.source.addItems(SOURCES)
        self.source.setToolTip("Which sequences to align")
        self.method = QComboBox()
        self.method.addItem(BUILT_IN, None)
        for key, recipe in runner.recipes().items():
            if recipe.get("kind") == "msa":
                self.method.addItem(recipe.get("label", key), key)
        self.method.setToolTip("MAFFT and Clustal Omega must be installed (see Map Reads, the same tools environment)")
        self.reference = QComboBox()
        self.reference.setToolTip("Identity is measured against this sequence, or against the consensus")
        self.reference.activated.connect(self._reference_changed)
        self.status = QLabel()
        self.canvas = MsaCanvas()
        self.canvas.clicked.connect(self._column_clicked)
        self.canvas.dragged.connect(self._columns_dragged)
        top = QHBoxLayout()
        for label, widget in (("Sequences", self.source), ("Method", self.method)):
            top.addWidget(QLabel(label))
            top.addWidget(widget)
        top.addWidget(QPushButton("Align", clicked=self.align_now))
        top.addWidget(QLabel("Identity to"))
        top.addWidget(self.reference)
        top.addStretch()
        box = QVBoxLayout(self)
        box.addLayout(top)
        box.addWidget(self.canvas, 1)
        box.addWidget(self.status)
        editor.selection_changed.connect(self._selection_changed)
        self._say("Choose the sequences and press Align.")
        self._auto()

    # ---- gathering the sequences ---------------------------------------------------------------------------------------------------
    def _open_name(self):
        return self.editor.fasta_record or self.editor.doc.name or "open sequence"

    def gather(self):
        """[(name, bases)] of the sequences chosen in "Sequences", or None when the choice was cancelled."""
        ed, kind = self.editor, self.source.currentText()
        if kind == SOURCES[2]:
            out = [(self._open_name(), ed.doc.sequence())]
            if ed.end > ed.start:
                out.append((f"selection {ed.start + 1}-{ed.end}", ed.selected_sequence()))
            return out
        if kind == SOURCES[1]:
            paths, _ = QFileDialog.getOpenFileNames(self, "Sequence files to align", "", "Sequences (*.gb *.gbk *.gbff *.fa *.fasta *.fna *.dna);;All files (*)")
            return None if not paths else self._from_files(paths)
        if not ed.path or not Path(ed.path).exists():
            return [(self._open_name(), ed.doc.sequence())]
        out = self._from_files([ed.path])
        # the open record may hold edits not saved yet: align what is on screen
        return [(n, ed.doc.sequence() if n == self._open_name() else s) for n, s in out]

    def _from_files(self, paths):
        out, taken = [], set()
        for path in paths:
            for rec in source_records(path):
                name = rec.name
                while name in taken:   # two files with a record of the same name: told apart by the file
                    name = f"{rec.name} ({Path(path).stem})"
                    rec = type(rec)(name=name, load=rec.load)
                taken.add(name)
                out.append((name, rec.load().sequence()))
        return out

    def _auto(self):
        """Align the records of a small multi-record file at once, so the view opens showing something."""
        try:
            if len(self.editor.fasta_records) < 2:
                return
            seqs = self.gather()
        except (OSError, ValueError):
            return
        if seqs and len(seqs) > 1 and sum(len(s) for _, s in seqs) <= AUTO_LIMIT:
            self._align(seqs, None)

    # ---- aligning ------------------------------------------------------------------------------------------------------------------
    def align_now(self):
        try:
            seqs = self.gather()
        except (OSError, ValueError) as e:
            QMessageBox.warning(self, "Multiple alignment", str(e))
            return
        if seqs is None:
            return
        self._align(seqs, self.method.currentData())

    def _align(self, seqs, tool):
        names, bases = [n for n, _ in seqs], [s for _, s in seqs]
        if len(seqs) < 2:
            self._say("Pick at least two sequences to align.")
            return
        try:
            if tool is None:
                alignment = run_busy(self, "Aligning…", lambda: msa.align(names, bases))
            else:
                alignment = self._run_tool(tool, names, bases)
                if alignment is None:
                    return
        except ValueError as e:
            QMessageBox.warning(self, "Multiple alignment", str(e))
            return
        self.sources = seqs
        self.set_alignment(alignment)

    def _run_tool(self, key, names, bases):
        """The Alignment MAFFT / Clustal Omega gives for the sequences (None if the run failed or was cancelled; its log says why).
        Sequences go in under n0, n1… so that no name can confuse the tool; the names are put back from the order."""
        recipe = runner.recipes()[key]
        folder = Path(tempfile.mkdtemp(prefix="opendna_msa_"))
        files = {"ref": str(folder / "sequences.fa"), "out": str(folder / "aligned")}
        (folder / "sequences.fa").write_text(msa.to_fasta([f"n{i}" for i in range(len(names))], bases), encoding="utf-8")
        run = runner.RunDialog(self, f"Aligning with {recipe.get('label', key)}", runner.steps_for(recipe, files, 4), folder)
        if run.exec() != QDialog.Accepted or not run.ok:
            self._say(f"{recipe.get('label', key)} did not finish: see its log (Cancel, install it, or pick the built-in aligner).")
            return None
        got, rows = msa.parse_aligned_fasta(Path(recipe["result"].format(**files)).read_text(encoding="utf-8"))
        by_id = dict(zip(got, rows))   # Clustal Omega may reorder: back to the order given
        if sorted(by_id) != sorted(f"n{i}" for i in range(len(names))):
            raise ValueError("The aligner's output does not hold the sequences that were sent.")
        return msa.Alignment(names, [by_id[f"n{i}"] for i in range(len(names))])

    def set_alignment(self, alignment):
        self.alignment = alignment
        self.reference.clear()
        self.reference.addItems(["Consensus"] + alignment.names)
        self.reference.setCurrentIndex(1 if alignment.reference is not None else 0)
        self.canvas.set_alignment(alignment)
        self._summary()

    def _reference_changed(self, index):
        if self.alignment is not None:
            self.alignment.reference = None if index == 0 else index - 1
            self.canvas.viewport().update()

    def _summary(self):
        a = self.alignment
        identical = sum(c == 1.0 for c in self.canvas.cons)
        self._say(f"{len(a.rows)} sequences, {a.width:,} columns, {identical:,} identical ({100 * identical / a.width:.1f} %). "
                  "Click a column to go to it in every sequence.")

    def _say(self, text):
        self.status.setText(text)

    # ---- column <-> position ---------------------------------------------------------------------------------------------------------
    def open_row(self):
        """The row that is the open sequence (its record, or the whole sequence in "Open sequence + selection"), or None."""
        if self.alignment is None:
            return None
        name = self._open_name()
        return self.alignment.names.index(name) if name in self.alignment.names else None

    def _column_clicked(self, col):
        a, row = self.alignment, self.open_row()
        places = "  ".join(f"{n}: {a.position(i, col) + 1:,}" + ("" if a.rows[i][col] != msa.GAP else "(gap)") for i, n in enumerate(a.names))
        self._say(f"Column {col + 1:,} ({a.consensus()[col]}, {100 * self.canvas.cons[col]:.0f} % conserved) · position in  {places}")
        if row is not None:
            self.own = True
            try:
                self.editor.set_cursor(a.position(row, col))
            finally:
                self.own = False

    def _columns_dragged(self, first, last):
        a, row = self.alignment, self.open_row()
        if row is not None:
            self.own = True
            try:
                self.editor.select(a.position(row, first), a.position(row, last + 1))
            finally:
                self.own = False

    def _selection_changed(self):
        """Mark the columns of the bases selected in the open sequence (not when this view made the selection)."""
        row, ed = self.open_row(), self.editor
        if self.own or row is None:
            return
        a = self.alignment
        if ed.end > ed.start and ed.end <= len(a.rows[row].replace(msa.GAP, "")):
            first, last = a.column_of(row, ed.start), a.column_of(row, ed.end - 1)
        else:
            first = last = a.column_of(row, ed.start) if ed.start < len(a.rows[row].replace(msa.GAP, "")) else None
        self.canvas.mark(first, last)
        if first is not None:
            self.canvas.ensure_visible(first)
