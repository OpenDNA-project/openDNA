"""Local Align: paste a DNA sequence, see where each piece of it matches the open sequence (both
strands, mismatches and gaps allowed; Rust `align`), and show the hits in the genome view.

A query can be a patchwork of pieces (a structural variant: a deletion, an insertion, a piece in
another order, an INVERTED piece). Each matching piece is its own row, ordered by its place in the
query, and a summary line plus a "before this piece" column say how consecutive pieces are joined
(gap in the query, gap in the template, strand change).

The dialog reads top to bottom: Query (the box, Paste from clipboard, Clear, and the "Compare with a second sequence" check box),
Search in (This sequence or Folder, and the one Align button), the result tabs (each with the buttons that act on its own list) and
the status line. With the second box (the check box) Align aligns the two with each other (the Alignment tab) and compares them by
a dot plot and a list of the stretches they share (the Dot plot tab, repeats.py), the same sequence in both boxes finding its own
repeats. The Folder choice aligns the query against every file of a folder (folder_search.py)."""
from types import SimpleNamespace

from PySide6.QtCore import Qt
from PySide6.QtGui import QBrush, QColor, QFont
from PySide6.QtWidgets import (
    QAbstractItemView, QApplication, QCheckBox, QGroupBox, QHBoxLayout, QLabel, QPlainTextEdit, QPushButton, QRadioButton,
    QSizePolicy, QTabWidget, QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget,
)

from opendna_core import align

import msa
from alignments import INVERTED_COLOR
from editor import clean_pasted
from msa_view import MsaCanvas
from repeats import PairRepeats
from tools import ToolDialog, empty_hint
from uikit import fit_columns, no_auto_default, plural, right_aligned, run_busy
import theme
import value_filter

MIN_QUERY = 12        # bp: shorter text from the clipboard is not offered as a query
MAX_QUERY = 200_000   # bp: longer queries are refused (the alignment takes too long and too much memory)
MAX_HITS = 30         # pieces / hits listed
FORWARD_TEXT = "forward"
INVERTED_TEXT = "↺ reverse complement (inverted)"


class QueryError(ValueError):
    """The text is not usable as a query; str(e) says why, in words for the user."""


def clean_query(text):
    """The bases of a query typed or pasted as text. Takes plain bases or FASTA (a header line
    `>name` is dropped; of several records only the first is used), ignores spaces, line breaks
    and digits (GenBank ORIGIN blocks) and alignment gap dashes, turns lowercase into uppercase
    and RNA's U into T. Raises QueryError when nothing usable is left."""
    text = text.strip()
    if not text:
        raise QueryError("Nothing to align: the query is empty")
    # FASTA: keep the lines of the first record only (a record starts at a ">" line).
    lines, records = [], 0
    for line in text.splitlines():
        if line.lstrip().startswith(">"):
            records += 1
            if records > 1:
                break
            continue
        lines.append(line)
    body = "\n".join(lines)
    try:
        bases = clean_pasted(body).replace("-", "").upper().replace("U", "T")
    except ValueError as e:
        bad = str(e).split(":", 1)[-1].split()
        if set(bad) & set("EFILPQZJOX*"):  # letters that exist only in protein sequences
            raise QueryError("This looks like a protein sequence; Local Align needs DNA "
                             f"(unexpected letters: {' '.join(bad[:10])})") from None
        raise QueryError(f"Not DNA: unexpected characters {' '.join(bad[:10])}") from None
    if not bases:
        raise QueryError("No DNA found" + (" (only a FASTA header)" if records else ""))
    if len(bases) > MAX_QUERY:
        raise QueryError(f"Query too long: {len(bases):,} bp (the limit is {MAX_QUERY:,} bp)")
    return bases


def clipboard_query():
    """(bases, message) from the system clipboard: bases = "" and a message saying why when it
    holds nothing usable."""
    text = QApplication.clipboard().text()
    if not text.strip():
        return "", "The clipboard is empty"
    try:
        return clean_query(text), ""
    except QueryError as e:
        return "", f"The clipboard has no DNA: {e}"


def strand_text(strand):
    return FORWARD_TEXT if strand > 0 else INVERTED_TEXT


def _junction(prev, cur, n, circular):
    """(query gap, template gap, label) between two pieces that follow each other in the query.
    The query gap is the unmatched query between them (negative: they overlap); the template gap is
    the template skipped going the way the query goes (negative: the template jumps back)."""
    qgap = cur.qstart - prev.qend
    if prev.strand != cur.strand:
        label = "inversion breakpoint (strand changes)"
        # The template pieces meet at their near ends: whichever pair is closest.
        tgap = min(abs(cur.start - prev.end), abs(prev.start - cur.end), abs(cur.start - prev.start), abs(cur.end - prev.end))
        return qgap, tgap, label
    tgap = cur.start - prev.end if prev.strand > 0 else prev.start - cur.end
    if circular and n:
        tgap = (tgap + n // 2) % n - n // 2  # the short way round the circle
    if tgap < 0:
        label = "template jumps back (duplication or pieces out of order)"
    elif qgap > 0 and tgap == 0:
        label = "insertion in the query"
    elif tgap > 0 and qgap <= 0:
        label = "deletion from the query"
    elif qgap > 0:
        label = "unmatched query in place of template"
    else:
        label = "contiguous"
    return qgap, tgap, label


def junction_text(prev, cur, n, circular):
    """The "before this piece" cell: how `cur` follows `prev` in the query."""
    qgap, tgap, label = _junction(prev, cur, n, circular)
    q = f"{qgap:,} bp unmatched in the query" if qgap > 0 else (f"queries overlap by {-qgap:,} bp" if qgap < 0 else "")
    t = f"{tgap:,} bp of the template skipped" if tgap > 0 else (f"template back {-tgap:,} bp" if tgap < 0 else "")
    detail = "; ".join(x for x in (q, t) if x)
    return f"{label}" + (f": {detail}" if detail else "")


def summarize(hits, n, circular):
    """One line about all the pieces, in query order: how many, which strands, the unmatched query
    between them and the structural events. `hits` must be in query order."""
    if not hits:
        return ""
    forward = sum(h.strand > 0 for h in hits)
    inverted = len(hits) - forward
    parts = [f"{forward} forward"] if forward else []
    if inverted:
        parts.append(f"{inverted} inverted (reverse complement)")
    noun = "piece" if len(hits) == 1 else "pieces"
    text = f"{len(hits)} {noun}: {', '.join(parts)}"
    gaps, events = [], []
    for i in range(1, len(hits)):
        qgap, tgap, label = _junction(hits[i - 1], hits[i], n, circular)
        if qgap > 0:
            gaps.append(f"{qgap:,} bp unmatched (#{i}→#{i + 1})")
        if label != "contiguous":
            events.append(f"{label.split(' (')[0]} (#{i}→#{i + 1})")
    if len(hits) > 1:
        text += "; gaps in the query between pieces: " + (", ".join(gaps) if gaps else "none")
        if events:
            text += "; " + ", ".join(events)
    return text


class LocalAlignDialog(ToolDialog):
    COLUMNS = ["#", "Query", "Template", "Strand", "Identity", "Score", "Differences", "Before this piece"]

    def __init__(self, parent, editor):
        super().__init__(parent)
        from folder_search import FolderSearchPanel   # (it imports this module's strand_text)
        self.win, self.editor, self.hits = parent, editor, []
        self.setWindowTitle("Local Align")
        self.resize(900, 620)
        self.query = QPlainTextEdit(placeholderText="Paste the query here (FASTA or plain bases), or use the clipboard button")
        self.query.setMaximumHeight(110)
        self.query.setPlainText(self._clipboard_dna())
        self.query_b = QPlainTextEdit(placeholderText="Second sequence (B), aligned with the first (A)", visible=False)
        self.query_b.setMaximumHeight(110)
        self.paired = QCheckBox("Compare with a second sequence", toolTip="Add a second box and align the two with each other, with a "
                                "dot plot and the stretches they share (the same sequence in both finds its repeats)")
        self.paired.toggled.connect(self._pair_toggled)
        # (wrapped, and unable to widen the window: a long message beside the buttons used to stretch it)
        self.status = QLabel(wordWrap=True, textFormat=Qt.PlainText)
        self.status.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Preferred)
        self.summary = QLabel(wordWrap=True)
        self.table = QTableWidget(0, len(self.COLUMNS))
        self.table.setHorizontalHeaderLabels(self.COLUMNS)
        value_filter.install(self.table, ["Strand"])   # right-click a header: filter by values
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        fit_columns(self.table, stretch=7)   # one number column: "#" (Qt's own row numbers are hidden)
        self.table.setWordWrap(False)
        empty_hint(self.table, "Paste a query and press Align")
        self.table.currentCellChanged.connect(lambda row, *_: self._select(row))
        editor.query_hit_changed.connect(self._picked)   # a row clicked in the Alignment view picks the same row here
        self.length = QLabel()   # the query's size, as the box is typed in
        self.query.textChanged.connect(self._query_length)
        self._query_length()
        # The three sections read top to bottom: what to align (Query), against what and go (Search in), the answers (tabs, each with
        # the buttons that act on its own list), and the status line.
        self.align_button = theme.mark_primary(QPushButton("Align", clicked=self.run, toolTip="Align the query (Enter)"))
        self.scope_this = QRadioButton("This sequence", checked=True)
        self.scope_folder = QRadioButton("Folder", toolTip="Align the query against every sequence file in a folder")
        self.scope_folder.toggled.connect(self._scope_changed)
        self.folder = FolderSearchPanel(parent, lambda: self._query(self.query), self.status.setText,
                                        lambda: self.tabs.setCurrentWidget(self.folder_tab), self._busy)
        self.canvas, self.pair = MsaCanvas(), PairRepeats()

        query_group = QGroupBox("Query")
        query_box = QVBoxLayout(query_group)
        boxes = QHBoxLayout()
        boxes.addWidget(self.query)
        boxes.addWidget(self.query_b)
        query_box.addLayout(boxes)
        helpers = QHBoxLayout()
        helpers.addWidget(QPushButton("Paste from clipboard", clicked=self.paste_clipboard,
                                      toolTip="Replace the query with the DNA on the clipboard (FASTA header, spaces and digits are dropped)"))
        helpers.addWidget(QPushButton("Clear", clicked=self.clear_query, toolTip="Empty the query box"))
        helpers.addWidget(self.length)
        helpers.addStretch(1)
        helpers.addWidget(self.paired)
        query_box.addLayout(helpers)

        self.search_group = QGroupBox()
        search_box = QVBoxLayout(self.search_group)
        go = QHBoxLayout()
        for w in (self.scope_this, self.scope_folder):
            go.addWidget(w)
        go.addStretch(1)
        go.addWidget(self.folder.cancel_button)
        go.addWidget(self.align_button)
        search_box.addLayout(go)
        search_box.addWidget(self.folder)

        hits = self._tab([self.summary, "Copy results", self.copy_hits, "Clear hits", lambda: self.clear()], self.table)
        self.folder_tab = self._tab([1, "Copy results", self.folder.copy_results], self.folder.table)
        self.tabs = QTabWidget()
        for widget, title in ((hits, "Hits in this sequence"), (self.folder_tab, "Folder results"), (self.canvas, "Alignment"),
                              (self.pair, "Dot plot")):
            self.tabs.addTab(widget, title)
        self._pair_toggled(False)

        status = QHBoxLayout()
        status.addWidget(self.status, 1)
        status.addWidget(self.folder.bar)
        layout = QVBoxLayout(self)
        layout.addWidget(query_group)
        layout.addWidget(self.search_group)
        layout.addWidget(self.tabs, 1)
        layout.addLayout(status)
        no_auto_default(self)
        self.align_button.setDefault(True)   # Enter aligns (the query box itself takes Enter as a new line)

    @staticmethod
    def _tab(bar, content):
        """A results tab: a row of buttons acting on this list (`bar`: a label to stretch, or 1 for a stretch, then pairs of a button's
        text and its handler) above the list."""
        tab = QWidget()
        box = QVBoxLayout(tab)
        box.setContentsMargins(0, 4, 0, 0)
        row = QHBoxLayout()
        row.setContentsMargins(0, 0, 0, 0)
        items = iter(bar)
        first = next(items)
        row.addWidget(first, 1) if isinstance(first, QWidget) else row.addStretch(first)
        for text in items:
            row.addWidget(QPushButton(text, clicked=next(items)))
        box.addLayout(row)
        box.addWidget(content)
        return tab

    def _busy(self, on):
        """A folder search starts or ends: Align waits for it (Cancel and the progress bar are the panel's)."""
        self.align_button.setEnabled(not on)

    def _query_length(self):
        try:
            self.length.setText(f"{len(clean_query(self.query.toPlainText())):,} bp")
        except QueryError:
            self.length.setText("")   # (said when Align is pressed)

    def clear_query(self):
        self.query.clear()
        self.query_b.clear()

    def _scope_changed(self, folder):
        """The Folder choice shows the folder row and its results; This sequence shows the hits."""
        self.folder.setVisible(folder and not self.paired.isChecked())
        if not self.paired.isChecked():
            self.tabs.setCurrentWidget(self.folder_tab if folder else self.tabs.widget(0))

    def _pair_toggled(self, on):
        """Two sequences: the second box and the Alignment and Dot plot tabs take the place of the search in the open sequence or a
        folder (Align then compares A with B), whose choice and folder row are hidden."""
        self.query_b.setVisible(on)
        self.search_group.setTitle("Compare A with B" if on else "Search in")
        for w in (self.scope_this, self.scope_folder):
            w.setVisible(not on)
        for i, shown in enumerate((not on, not on, on, on)):
            self.tabs.setTabVisible(i, shown)
        folder = self.scope_folder.isChecked()
        self.folder.setVisible(folder and not on)
        self.tabs.setCurrentIndex(2 if on else 1 if folder else 0)

    def copy_hits(self):
        """The hits table as tab-separated text, to paste in a lab notebook or a spreadsheet."""
        table = self.table
        rows = [[table.horizontalHeaderItem(c).text() for c in range(table.columnCount())]]
        rows += [[table.item(r, c).text() for c in range(table.columnCount())] for r in range(table.rowCount())]
        QApplication.clipboard().setText("\n".join("\t".join(row) for row in rows))
        self.status.setText(f"Copied {plural(table.rowCount(), 'hit')}")

    @staticmethod
    def _query(box):
        """The bases of a query box; raises QueryError (a ValueError) saying why they are not usable."""
        query = clean_query(box.toPlainText())
        if len(query) < 4:
            raise QueryError("Query too short")
        return query

    def done(self, result):
        self.folder.stop()   # closing the dialog stops a folder search too
        super().done(result)

    @staticmethod
    def _clipboard_dna():
        """The clipboard text if it is DNA (or FASTA) of a useful length, else nothing."""
        bases, _ = clipboard_query()
        return bases if len(bases) >= MIN_QUERY else ""

    def paste_clipboard(self):
        """Put the DNA on the clipboard in the query box; say why when there is none. True if done."""
        bases, message = clipboard_query()
        if not bases:
            self.status.setText(message)
            return False
        self.query.setPlainText(bases)
        self.status.setText(f"Query: {len(bases):,} bp from the clipboard")
        return True

    @classmethod
    def align_clipboard(cls, parent, editor):
        """Tools > Align Clipboard: open the dialog with the clipboard as the query and align it."""
        dlg = cls(parent, editor)
        if dlg.paste_clipboard():
            dlg.run()
        dlg.exec()

    def run_pair(self):
        """Paired mode: A and B aligned with each other, then the dot plot and shared stretches."""
        try:
            a, b = self._query(self.query), self._query(self.query_b)
            alignment = run_busy(self, "Aligning…", lambda: msa.align(["A", "B"], [a, b]))
        except ValueError as e:
            self.status.setText(str(e))
            return
        self.canvas.set_alignment(alignment)
        self.pair.show_pair(a, b)
        self.tabs.setCurrentWidget(self.canvas)
        self.status.setText(f"A {len(a):,} bp, B {len(b):,} bp: {alignment.identity(1):.1f} % identity over {alignment.width:,} columns")

    def run(self):
        if self.paired.isChecked():
            return self.run_pair()
        if self.scope_folder.isChecked():
            return self.folder.start()
        ed = self.editor
        try:
            query = self._query(self.query)
        except QueryError as e:
            self.clear(str(e))
            return
        text, circular = ed.doc.sequence(), ed.doc.circular   # (up to 200 kb of query: in a thread, behind a busy dialog)
        found = run_busy(self, "Aligning…", lambda: align(text, query, circular, MAX_HITS))
        hits = [SimpleNamespace(start=h[0], end=h[1], strand=h[2], score=h[3], identity=h[4], qstart=h[5],
                                qend=h[6], mismatches=h[7], insertions=h[8], deletions=h[9]) for h in found]
        # One row per piece, in the order the pieces have in the query.
        self.hits = sorted(hits, key=lambda h: (h.qstart, h.qend, h.start))
        self.fill()
        ed.set_query_hits(self.hits)
        limit = " (the limit; weaker ones are left out)" if len(self.hits) >= MAX_HITS else ""
        self.status.setText(f"{plural(len(self.hits), 'hit')}{limit}" if self.hits else "No hit")
        if self.hits:
            self.win.set_genome_visible(True)  # show the Alignment view and frame the pieces
            n = ed.length
            best = max(self.hits, key=lambda h: h.score)
            lo, hi = min(h.start for h in self.hits), max(h.end for h in self.hits)
            if hi > n or hi - lo > n * 0.9:  # pieces all round a circle: frame the best one
                lo, hi = best.start, best.end
            margin = max(20, (hi - lo) // 10)
            self.win.genome_view.set_view(lo - margin, min(hi, n) + margin)
            self.table.selectRow(self.hits.index(best) if len(self.hits) == 1 else 0)

    def fill(self):
        n, circular = self.editor.length, self.editor.doc.circular
        self.summary.setText(summarize(self.hits, n, circular))
        self.table.setRowCount(len(self.hits))
        inverted_brush = QBrush(QColor(INVERTED_COLOR))
        bold = QFont(self.table.font())
        bold.setBold(True)
        for i, h in enumerate(self.hits):
            before = junction_text(self.hits[i - 1], h, n, circular) if i else "start of the query"
            values = [str(i + 1), f"{h.qstart + 1:,}–{h.qend:,}", f"{h.start + 1:,}–{(h.end - 1) % n + 1:,}",
                      strand_text(h.strand), f"{h.identity:.1%}", str(h.score),
                      f"{plural(len(h.mismatches), 'mismatch')}, {len(h.insertions)} ins, {len(h.deletions)} del", before]
            for c, v in enumerate(values):
                item = QTableWidgetItem(v)
                if c in (0, 4, 5):
                    right_aligned(item)
                if h.strand < 0 and c in (0, 3):  # the inverted pieces stand out
                    item.setBackground(inverted_brush)
                    item.setForeground(QBrush(QColor("white")))
                    item.setFont(bold)
                self.table.setItem(i, c, item)

    def clear(self, message=""):
        """Remove the hits (and their track); `message` goes in the status line."""
        self.hits = []
        self.fill()
        self.editor.set_query_hits([])
        self.status.setText(message)

    def refresh(self):
        if self.hits and not self.editor.query_hits:  # an edit made the hits stale
            self.clear()

    def _select(self, row):
        if 0 <= row < len(self.hits):
            self.editor.select_query_hit(row)   # the Alignment view marks the same row

    def _picked(self, index):
        """The Alignment view picked hit `index`: select its row in the table (same order, same numbers)."""
        if 0 <= index < self.table.rowCount() and self.table.currentRow() != index:
            self.table.selectRow(index)
