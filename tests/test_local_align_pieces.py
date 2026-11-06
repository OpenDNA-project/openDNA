"""Local Align as a structural-variant reader (docs/TEST_CATALOGUE.md section 8): a query made of [A] + 200 unmatched bases +
[B] + [C, inverted] + [D] lists every piece as its own row, in query order, with the strand in words, the gaps between pieces
and a summary line; and the query can come from the clipboard (FASTA header, spaces, digits dropped; clear messages when there
is no DNA). The aligner itself is tested in test_align_structural.py."""
import os
import random
import sys
import tempfile
from pathlib import Path
from types import SimpleNamespace

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ["OPENDNA_DATA_DIR"] = tempfile.mkdtemp(prefix="opendna_test_")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "app"))

from PySide6.QtCore import Qt
from PySide6.QtGui import QAction
from PySide6.QtTest import QTest
from PySide6.QtWidgets import QApplication

from opendna_core import Document

import local_align
from editor import reverse_complement as rc
from local_align import (FORWARD_TEXT, INVERTED_TEXT, MAX_QUERY, LocalAlignDialog, QueryError, clean_query, junction_text,
                         summarize)
from main import MainWindow

app = QApplication.instance() or QApplication([])


@pytest.fixture
def errors(monkeypatch):
    caught = []
    monkeypatch.setattr(sys, "excepthook", lambda *exc: caught.append(exc))
    yield caught
    assert not caught, caught


def dna(n, seed):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def open_window(seq, circular=False):
    win = MainWindow()
    win.resize(1200, 800)
    win.show()
    win.show_document(Document.new_sequence("t", seq, circular), Path("t.gb"))
    app.processEvents()
    return win


def close(win):
    win.editor.doc.modified = False
    win.close()


T = dna(2000, 73)   # (seed chosen so no chance base extends a piece: the coordinates below are exact)
# A (forward) + 200 unmatched + B (forward) + C (INVERTED) + D (forward); the junction bases are fixed so no chance base extends a piece
A_AT, B_AT, C_AT, D_AT = 100, 600, 1100, 1600
QUERY = T[A_AT:A_AT + 100] + dna(200, 5) + T[B_AT:B_AT + 100] + rc(T[C_AT:C_AT + 100]) + T[D_AT:D_AT + 100]


def hit(qstart, qend, start, end, strand=1):
    return SimpleNamespace(qstart=qstart, qend=qend, start=start, end=end, strand=strand)


# ---- the query text ----------------------------------------------------------------------------------------------------

def test_query_text_and_junction_texts():
    """clean_query accepts FASTA, lowercase, RNA, digits and gaps and says what is wrong otherwise; the junction and summary texts
    name insertions, deletions, inversions and duplications."""
    # --- clean_query accepts ---
    assert clean_query(">seq1 description\nacgt acgt\n  ggcc\n") == "ACGTACGTGGCC"
    assert clean_query("1 acgtacgt 9\n11 ggcctt") == "ACGTACGTGGCCTT"            # a GenBank ORIGIN block
    assert clean_query("ACGU UGCA") == "ACGTTGCA"                                  # RNA: U becomes T
    assert clean_query("AC-GT--AC") == "ACGTAC"                                    # alignment gap dashes
    assert clean_query("ACGTNNRYACGT") == "ACGTNNRYACGT"                           # IUPAC codes are DNA
    # several FASTA records: the first one only
    assert clean_query(">a\nAAAACCCC\n>b\nGGGGTTTT\n") == "AAAACCCC"
    assert clean_query("\n\n  >a\r\nACGT\r\nACGT\r\n") == "ACGTACGT"
    # --- clean_query errors ---
    for text, words in (("", "empty"), ("  \n ", "empty"), (">only a header", "only a FASTA header"), ("123 456", "No DNA")):
        with pytest.raises(QueryError, match=words):
            clean_query(text)
    with pytest.raises(QueryError, match="protein"):
        clean_query("MKVLAAGIVGLLLAQPSEEF")                                          # amino acids
    with pytest.raises(QueryError, match="Not DNA"):
        clean_query("ACGT hello 12")                                                 # H, O: not protein-only letters
    with pytest.raises(QueryError, match=r"too long.*200,000"):
        clean_query("A" * (MAX_QUERY + 1))
    assert len(clean_query("ACGT" * (MAX_QUERY // 4))) == MAX_QUERY                 # exactly the limit is fine
    # --- junction and summary texts ---
    n = 5000
    a = hit(0, 100, 1000, 1100)
    assert junction_text(a, hit(100, 200, 1100, 1200), n, False) == "contiguous"
    assert junction_text(a, hit(300, 400, 1100, 1200), n, False).startswith("insertion in the query: 200 bp unmatched")
    assert junction_text(a, hit(100, 200, 1150, 1250), n, False).startswith("deletion from the query: 50 bp of the template skipped")
    assert junction_text(a, hit(100, 200, 900, 1000), n, False).startswith("template jumps back")
    assert junction_text(a, hit(100, 200, 1500, 1600, -1), n, False).startswith("inversion breakpoint")
    # two reverse pieces: the template runs backwards along the query
    r1, r2 = hit(0, 100, 1000, 1100, -1), hit(100, 200, 900, 1000, -1)
    assert junction_text(r1, r2, n, False) == "contiguous"
    # a circle: the short way round the origin is a small gap, not a jump back of nearly the whole template
    assert junction_text(hit(0, 100, 4900, 5000), hit(100, 200, 50, 150), n, True).startswith("deletion from the query: 50 bp")
    assert summarize([a], n, False) == "1 piece: 1 forward"
    s = summarize([a, hit(300, 400, 1100, 1200), hit(400, 500, 3000, 3100, -1)], n, False)
    assert s.startswith("3 pieces: 2 forward, 1 inverted (reverse complement); gaps in the query between pieces: 200 bp unmatched (#1→#2)")
    assert "insertion in the query (#1→#2)" in s and "inversion breakpoint (#2→#3)" in s
    assert summarize([], n, False) == ""


# ---- the table, the summary and the track -----------------------------------------------------------------------------

def test_structural_variant_query(errors):
    """A structural-variant query gives one row per piece in query order (table, summary, colours, row click, alignment bars),
    and on a circular template typed as its reverse complement the same four places come back."""
    # --- one row per piece ---
    win = open_window(T)
    ed = win.editor
    dlg = LocalAlignDialog(win, ed)
    dlg.query.setPlainText(QUERY)
    dlg.run()
    app.processEvents()
    assert [(h.qstart, h.qend, h.strand) for h in dlg.hits] == [(0, 100, 1), (300, 400, 1), (400, 500, -1), (500, 600, 1)]
    assert [(h.start, h.end) for h in dlg.hits] == [(100, 200), (600, 700), (1100, 1200), (1600, 1700)]
    t = dlg.table
    assert t.rowCount() == 4
    assert [t.item(r, 1).text() for r in range(4)] == ["1–100", "301–400", "401–500", "501–600"]
    assert [t.item(r, 2).text() for r in range(4)] == ["101–200", "601–700", "1,101–1,200", "1,601–1,700"]
    assert [t.item(r, 3).text() for r in range(4)] == [FORWARD_TEXT, FORWARD_TEXT, INVERTED_TEXT, FORWARD_TEXT]
    assert all(t.item(r, 4).text() == "100.0%" for r in range(4)) and dlg.status.text() == "4 hits"
    # how each piece follows the previous one
    assert t.item(0, 7).text() == "start of the query"
    assert t.item(1, 7).text() == ("unmatched query in place of template: 200 bp unmatched in the query; "
                                   "400 bp of the template skipped")
    assert "inversion breakpoint" in t.item(2, 7).text() and "inversion breakpoint" in t.item(3, 7).text()
    # the summary line above the table
    s = dlg.summary.text()
    assert s.startswith("4 pieces: 3 forward, 1 inverted (reverse complement)")
    assert "gaps in the query between pieces: 200 bp unmatched (#1→#2)" in s and "inversion breakpoint (#2→#3)" in s
    # the inverted row is marked with a colour of its own, the forward ones are not
    assert t.item(2, 3).background().color().name() == local_align.INVERTED_COLOR and t.item(1, 3).background().color().name() != local_align.INVERTED_COLOR
    # clicking a row selects its template span
    ed.set_cursor(0)
    QTest.mouseClick(t.viewport(), Qt.LeftButton, pos=t.visualItemRect(t.item(2, 2)).center())
    assert (ed.start, ed.end) == (1100, 1200)
    QTest.mouseClick(t.viewport(), Qt.LeftButton, pos=t.visualItemRect(t.item(3, 2)).center())
    assert (ed.start, ed.end) == (1600, 1700)
    # the alignment view: every piece is a bar, the inverted one says so
    canvas = win.genome_view.canvas
    win.grab()
    bars = {(ed.query_hits[a[1]].start, ed.query_hits[a[1]].end): a[2] for _, a in canvas.hits if a[0] == "qhit"}
    assert set(bars) == {(100, 200), (600, 700), (1100, 1200), (1600, 1700)}
    assert "inverted" in bars[(1100, 1200)] and "forward" in bars[(600, 700)] and "inverted" not in bars[(600, 700)]
    dlg.clear()
    assert dlg.summary.text() == "" and dlg.table.rowCount() == 0
    close(win)
    # --- circular template, reverse complement ---
    win = open_window(T, circular=True)
    dlg = LocalAlignDialog(win, win.editor)
    dlg.query.setPlainText(rc(QUERY))
    dlg.run()
    # the same four places; typed backwards the order of the pieces reverses and every strand flips
    assert [(h.start, h.end, h.strand) for h in dlg.hits] == [(1600, 1700, -1), (1100, 1200, 1), (600, 700, -1), (100, 200, -1)]
    assert [(h.qstart, h.qend) for h in dlg.hits] == [(0, 100), (100, 200), (200, 300), (500, 600)]
    assert dlg.summary.text().startswith("4 pieces: 1 forward, 3 inverted (reverse complement)")
    close(win)


# ---- query from the clipboard -----------------------------------------------------------------------------------------

def test_query_box_and_clipboard(errors, monkeypatch):
    """Paste from the clipboard cleans the text and explains when there is no DNA; the query box takes lowercase, RNA, protein
    and huge text; the Align Clipboard menu entry opens the dialog with the clipboard aligned."""
    # --- paste from the clipboard ---
    win = open_window(T)
    dlg = LocalAlignDialog(win, win.editor)
    clip = QApplication.clipboard()
    dlg.query.setPlainText("keep me")
    for text, message in (("", "The clipboard is empty"), ("  \n", "The clipboard is empty"),
                          ("Meeting at 3pm, bring coffee", "no DNA"), (">header only", "only a FASTA header"),
                          ("MKVLAAGIVGLLLAQPSEEF", "protein")):
        clip.setText(text)
        assert dlg.paste_clipboard() is False
        assert message in dlg.status.text() and dlg.query.toPlainText() == "keep me", text   # the query is left alone
    clip.setText(">my query\n" + "\n".join(f"{i * 10 + 1:>5} {T[i * 10:i * 10 + 10].lower()}" for i in range(6)))
    assert dlg.paste_clipboard() is True
    assert dlg.query.toPlainText() == T[:60] and dlg.status.text() == "Query: 60 bp from the clipboard"
    close(win)
    # --- query box input ---
    win = open_window(T)
    dlg = LocalAlignDialog(win, win.editor)
    # typed lowercase, with RNA's U in place of T: found like the DNA
    dlg.query.setPlainText(T[300:400].lower().replace("t", "u"))
    dlg.run()
    assert [(h.start, h.end, h.strand) for h in dlg.hits] == [(300, 400, 1)]
    # a protein, an empty box, a huge query: a message each, no hits, no exception
    for text, message in (("MKVLAAGIVGLLLAQPSEEF", "protein"), ("", "empty"), ("ACGT" * (MAX_QUERY // 4 + 1), "too long")):
        dlg.query.setPlainText(text)
        dlg.run()
        assert message in dlg.status.text() and dlg.hits == [] and dlg.table.rowCount() == 0, message
    close(win)
    # --- Align Clipboard menu entry ---
    win = open_window(T)
    entry = next(a for a in win.findChildren(QAction) if a.text() == "Align &Clipboard…")
    shown = []
    monkeypatch.setattr(LocalAlignDialog, "exec", lambda self: shown.append(self) or 0)   # do not block on a modal dialog
    QApplication.clipboard().setText(">q\n" + QUERY)
    entry.trigger()
    dlg = shown[0]
    assert dlg.query.toPlainText() == QUERY and len(dlg.hits) == 4 and dlg.table.rowCount() == 4
    assert [h.strand for h in dlg.hits] == [1, 1, -1, 1]
    # nothing usable on the clipboard: the dialog still opens and says why, without hits
    QApplication.clipboard().setText("")
    entry.trigger()
    assert "clipboard is empty" in shown[1].status.text() and shown[1].hits == []
    # the Tools menu offers both entries
    texts = [a.text() for a in win.findChildren(QAction)]
    assert "&Local Align…" in texts and "Align &Clipboard…" in texts
    close(win)


