"""Regression tests for bugs the TESTING.md tests found (docs/TODO.md, "Known problems"): each fails without its fix."""
import os
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ["GENOMEEDIT_DATA_DIR"] = tempfile.mkdtemp(prefix="genomeedit_test_")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "app"))

import pytest

import enzymes
import genomeedit_core as core
from genomeedit_core import Document, Fasta


def test_tracks_enzymes_and_files(tmp_path):
    """Regressions in data handling: a track over the limit keeps the first items of the window; a GFF line with an empty name
    falls back to the next attribute; a cut off the end of a linear sequence is no cut; a FASTA with a byte order mark is indexed
    and edited; replacing a record of a CRLF FASTA keeps CRLF; a failed save leaves no temporary file."""
    # --- track over the limit ---
    p = tmp_path / "t.bed"
    p.write_text("".join(f"c\t{s}\t{s + 5}\tf{s:04d}\n" for s in range(0, 1000, 10)))   # 100 items
    t = core.Track(str(p))
    items, more = t.query("c", 0, 1000, 10)
    assert more is True and [i[3] for i in items] == [f"f{s:04d}" for s in range(0, 100, 10)], "the first ten, in start order"
    items, more = t.query("c", 500, 1000, 3)
    assert [i[3] for i in items] == ["f0500", "f0510", "f0520"] and more is True
    assert t.query("c", 0, 1000, 100)[1] is False, "exactly the limit: nothing is left out"
    # an item that starts before the window but overlaps it counts, and comes first
    assert t.query("c", 503, 520, 5)[0][0][3] == "f0500"
    # --- GFF empty name ---
    p = tmp_path / "t.gff3"
    p.write_text("##gff-version 3\nc\ts\tgene\t10\t50\t.\t+\t.\tID=g1;Name=;gene=lacZ\n"
                 "c\ts\tgene\t60\t90\t.\t+\t.\tID=g2;Name=\n")
    t = core.Track(str(p))
    names = [i[3] for i in t.query("c", 0, 100, 10)[0]]
    assert names == ["lacZ", "g2"], "an empty Name= is no name: gene, then ID"
    # --- cut off the end of a linear sequence ---
    bsai = next(e for e in enzymes.COMMON if e.name == "BsaI")       # GGTCTC(1/5): the bottom strand is cut 5 bases on
    assert enzymes.find_cuts("GGTCTCA", False, [bsai]) == {}, "the bottom-strand cut would be past the end"
    assert enzymes.find_cuts("GGTCTCAAAAAAAAAA", False, [bsai])["BsaI"] == [(7, 11)]
    assert enzymes.find_cuts("GGTCTCA", True, [bsai])["BsaI"], "on a circle it wraps round and is a cut"
    # --- FASTA with a byte order mark ---
    path = tmp_path / "bom.fa"
    path.write_bytes(b"\xef\xbb\xbf" + FA.encode())
    fasta = Fasta(str(path))
    assert [(n, length) for n, length in fasta.sequences()] == [("chr1", 14), ("chr2", 12)]
    assert fasta.fetch("chr1", 0, 14) == "ACGTACGTACGTAC" and fasta.fetch("chr2", 0, 4) == "TTTT"
    doc = Document.from_fasta(str(path), "chr1")
    assert doc.sequence() == "ACGTACGTACGTAC"
    out = tmp_path / "out.fa"
    Document.new_sequence("chr1", "GGGG", False).write_fasta(str(path), str(out))      # the first record, whose header has the mark
    assert out.read_bytes() == b"\xef\xbb\xbf>chr1 first\nGGGG\n>chr2\nTTTTTTTTTT\nTT\n"
    # --- CRLF FASTA ---
    path = tmp_path / "crlf.fa"
    path.write_bytes(FA.replace("\n", "\r\n").encode())
    out = tmp_path / "out.fa"
    Document.new_sequence("chr2", "G" * 70, False).write_fasta(str(path), str(out))
    data = out.read_bytes()
    assert b"\n" not in data.replace(b"\r\n", b""), "every line ends with CRLF"
    assert data.startswith(b">chr1 first\r\nACGTACGTAC\r\nGTAC\r\n>chr2\r\n" + b"G" * 60 + b"\r\n" + b"G" * 10 + b"\r\n")
    # --- failed save ---
    doc = Document.new_sequence("t", "ACGT" * 5, False)
    target = tmp_path / "taken.gb"
    target.mkdir()                                                   # a folder in the way: the rename fails
    with pytest.raises(OSError):
        doc.save_genbank(str(target))
    assert not list(tmp_path.glob("*.tmp")), "no half-finished file is left beside it"
    fa = tmp_path / "g.fa"
    fa.write_text(FA)
    blocked = tmp_path / "blocked.fa"
    blocked.mkdir()
    with pytest.raises(OSError):
        Document.new_sequence("chr1", "ACGT", False).write_fasta(str(fa), str(blocked))
    assert not list(tmp_path.glob("*.tmp"))
    # (the fuzzy-location regression is checked, with a save, in test_core_functions.py::test_locations_in_genbank_files)


FA = ">chr1 first\nACGTACGTAC\nGTAC\n>chr2\nTTTTTTTTTT\nTT\n"


def test_closing_a_pane_does_not_leave_the_editor_calling_a_deleted_view(monkeypatch):
    """Editing after the Protein or Map pane was closed raised "Internal C++ object already deleted" (a lambda connected to the editor)."""
    from PySide6.QtCore import QEvent
    from PySide6.QtWidgets import QApplication
    from main import MainWindow
    from map_view import MapView
    from protein import ProteinView
    app = QApplication.instance() or QApplication([])
    errors = []
    monkeypatch.setattr(sys, "excepthook", lambda *exc: errors.append(exc))
    win = MainWindow()
    win.show()
    win.show_document(Document.new_sequence("t", "ACGT" * 50, True), Path("t.gb"))
    ed, ws = win.editor, win.workspace
    protein, map_view = ws.ensure(ProteinView), ws.ensure(MapView)
    for view in (protein, map_view):
        pane = next(p for p in ws.panes() if p.view is view)
        ws.close_pane(pane)
    QApplication.sendPostedEvents(None, QEvent.DeferredDelete)
    app.processEvents()
    ed.set_cursor(3)
    ed.type_bases("A")
    app.processEvents()
    win.editor.doc.modified = False
    win.close()
    assert not errors, errors
