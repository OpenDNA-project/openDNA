"""Local Align: paste a DNA sequence, see where each piece of it matches the open sequence (both
strands, mismatches and gaps allowed; Rust `align`), and show the hits in the genome view.

A query can be a patchwork of pieces (a structural variant: a deletion, an insertion, a piece in
another order, an INVERTED piece). Each matching piece is its own row, ordered by its place in the
query, and a summary line plus a "before this piece" column say how consecutive pieces are joined
(gap in the query, gap in the template, strand change)."""
from types import SimpleNamespace

from PySide6.QtCore import Qt
from PySide6.QtGui import QBrush, QColor, QFont
from PySide6.QtWidgets import (
    QAbstractItemView, QApplication, QHBoxLayout, QLabel, QPlainTextEdit, QPushButton, QSizePolicy, QTableWidget,
    QTableWidgetItem, QVBoxLayout,
)

from opendna_core import align

from alignments import INVERTED_COLOR
from editor import clean_pasted
from tools import ToolDialog, empty_hint
from uikit import fit_columns, no_auto_default, plural, right_aligned, run_busy
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
        self.win, self.editor, self.hits = parent, editor, []
        self.setWindowTitle("Local Align")
        self.resize(760, 520)
        self.query = QPlainTextEdit(placeholderText="Paste the query here (FASTA or plain bases), or use the clipboard button")
        self.query.setMaximumHeight(110)
        self.query.setPlainText(self._clipboard_dna())
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
        buttons = QHBoxLayout()
        align_button = QPushButton("Align", clicked=self.run)
        buttons.addWidget(align_button)
        buttons.addWidget(QPushButton("Paste from clipboard", clicked=self.paste_clipboard,
                                      toolTip="Replace the query with the DNA on the clipboard (FASTA header, spaces and digits are dropped)"))
        buttons.addWidget(QPushButton("Clear hits", clicked=lambda: self.clear()))
        buttons.addWidget(self.status, 1)
        layout = QVBoxLayout(self)
        layout.addWidget(self.query)
        layout.addLayout(buttons)
        layout.addWidget(self.summary)
        layout.addWidget(self.table)
        no_auto_default(self)
        align_button.setDefault(True)   # Enter aligns (the query box itself takes Enter as a new line)

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

    def run(self):
        ed = self.editor
        try:
            query = clean_query(self.query.toPlainText())
        except QueryError as e:
            self.clear(str(e))
            return
        if len(query) < 4:
            self.clear("Query too short")
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
