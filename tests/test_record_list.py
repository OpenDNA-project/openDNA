"""The map's chromosome list of a file of several records (GenBank or FASTA): a line dragged to another place reorders the
records in the file (each copied byte for byte); right-click > Move Up / Down, Export as GenBank (one record as a file of its
own, with its unsaved edits) and Remove from File."""
import os
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ["OPENDNA_DATA_DIR"] = tempfile.mkdtemp(prefix="opendna_test_")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "app"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import pytest
from PySide6.QtCore import QEvent, QPointF, Qt
from PySide6.QtGui import QMouseEvent
from PySide6.QtWidgets import QApplication

from opendna_core import Document, Fasta

import multigb
from main import MainWindow
from sample_files import PLASMID, PLASMID2

app = QApplication.instance() or QApplication([])


def three_records(folder):
    """A GenBank file of three records: the two demo plasmids and a copy of the first named "copy"."""
    if not (PLASMID.exists() and PLASMID2.exists()):
        pytest.skip("demo plasmids missing")
    first = PLASMID.read_text(encoding="utf-8").rstrip("\n") + "\n"
    copy = first.replace(Document(str(PLASMID)).name, "copy", 1)
    path = folder / "three.gb"
    path.write_text(first + PLASMID2.read_text(encoding="utf-8").rstrip("\n") + "\n" + copy, encoding="utf-8")
    return path


def chunks(path):
    """{name: the bytes of the record} of a GenBank file."""
    data = path.read_bytes()
    _, spans = multigb.record_spans(path)
    return {n: data[s:e] for n, s, e in spans}


def names(path):
    return [n for n, *_ in multigb.records(path)]


def opened(path):
    win = MainWindow()
    win.resize(1200, 800)
    win.show()
    assert win.load(str(path))
    app.processEvents()
    return win


def test_rewrite_records_reorders_and_removes_byte_for_byte(tmp_path):
    path = three_records(tmp_path)
    before = chunks(path)
    a, b, c = names(path)
    multigb.rewrite_records(path, [c, a, b])
    assert names(path) == [c, a, b] and chunks(path) == before
    multigb.rewrite_records(path, [c, b])
    assert names(path) == [c, b] and chunks(path) == {c: before[c], b: before[b]}
    for bad in ([], [c, c], [c, "nope"]):
        with pytest.raises(OSError):
            multigb.rewrite_records(path, bad)
    assert names(path) == [c, b] and not list(tmp_path.glob("*.saving"))
    # FASTA: a byte order mark, Windows line endings, no line ending at the end, a stale index
    fa = tmp_path / "g.fa"
    fa.write_bytes(b"\xef\xbb\xbf>chr1 one\r\nACGT\r\nAC\r\n>chr2\r\nGGGG\r\n>chr3\r\nTT")
    Fasta(str(fa))
    multigb.rewrite_records(fa, ["chr3", "chr1", "chr2"])
    assert fa.read_bytes() == b"\xef\xbb\xbf>chr3\r\nTT\r\n>chr1 one\r\nACGT\r\nAC\r\n>chr2\r\nGGGG\r\n"
    assert [n for n, _ in Fasta(str(fa)).sequences()] == ["chr3", "chr1", "chr2"]
    assert Document.from_fasta(str(fa), "chr1").sequence() == "ACGTAC"
    multigb.rewrite_records(fa, ["chr2"])
    assert dict(Fasta(str(fa)).sequences()) == {"chr2": 4}


def test_a_line_dragged_in_the_map_reorders_the_file(tmp_path):
    """Drag the first line below the last: the file's records follow, the list and the picker too, the open record (with its
    unsaved edit) stays open, and a plain click still opens a record."""
    path = three_records(tmp_path)
    before = chunks(path)
    win = opened(path)
    a, b, c = names(path)
    assert [n for n, _ in win.editor.fasta_records] == [a, b, c] and win.editor.fasta_record == a, "listed in file order"
    win.editor.set_cursor(3)
    win.editor.type_bases("GG")
    view = win.map_view
    view.resize(600, 500)
    row = lambda i: QPointF(30, 2 + i * view.row_h() + view.row_h() / 2)
    press = lambda pt: QMouseEvent(QEvent.MouseButtonPress, pt, Qt.LeftButton, Qt.LeftButton, Qt.NoModifier)
    move = lambda pt: QMouseEvent(QEvent.MouseMove, pt, Qt.NoButton, Qt.LeftButton, Qt.NoModifier)
    release = lambda pt: QMouseEvent(QEvent.MouseButtonRelease, pt, Qt.LeftButton, Qt.NoButton, Qt.NoModifier)
    view.mousePressEvent(press(row(0)))
    view.mouseMoveEvent(move(row(1)))
    end = QPointF(30, 2 + 3 * view.row_h())   # the gap after the last line
    view.mouseMoveEvent(move(end))
    assert view.list_drop == 3
    view.grab()   # draws the drop line
    view.mouseReleaseEvent(release(end))
    app.processEvents()
    assert names(path) == [b, c, a] and chunks(path) == before, "reordered, every record unchanged"
    ed = win.editor
    assert [n for n, _ in ed.fasta_records] == [b, c, a] and [win.picker.itemText(i) for i in range(3)] == [b, c, a]
    assert ed.fasta_record == a and ed.doc.modified, "the open record and its edit are kept"
    assert win.save() and multigb.read_record(path, a).sequence()[3:5] == "GG" and names(path) == [b, c, a]
    # a click (no drag) opens the line
    view = win.map_view
    view.resize(600, 500)
    view.mousePressEvent(press(row(1)))
    view.mouseReleaseEvent(release(row(1)))
    app.processEvents()
    assert win.editor.fasta_record == c and names(path) == [b, c, a]


def test_right_click_menu_moves_exports_and_removes(tmp_path, monkeypatch):
    path = three_records(tmp_path)
    win = opened(path)
    a, b, c = names(path)
    menu = win.map_view.record_menu(a)
    actions = {x.text(): x for x in menu.actions() if x.text()}
    assert not actions[f"&Open {a}"].isEnabled() and not actions["Move &Up"].isEnabled() and actions["Move &Down"].isEnabled()
    actions["Move &Down"].trigger()
    app.processEvents()
    assert names(path) == [b, a, c]
    actions = {x.text(): x for x in win.map_view.record_menu(c).actions() if x.text()}
    assert not actions["Move &Down"].isEnabled()
    actions["Move &Up"].trigger()
    app.processEvents()
    assert names(path) == [b, c, a]
    # --- export: the open record with its unsaved edit, and one read from the file ---
    win.editor.set_cursor(0)
    win.editor.type_bases("TTTT")
    edited = win.editor.doc.sequence()
    out = tmp_path / "out"
    out.mkdir()
    assert win.export_record(a, out / "a.gb") and win.export_record(b, out / "b.gb")
    assert Document(str(out / "a.gb")).sequence() == edited and Document(str(out / "a.gb")).name == a
    assert Document(str(out / "b.gb")).sequence() == multigb.read_record(path, b).sequence()
    assert not multigb.is_multi_genbank(out / "b.gb") and win.editor.doc.modified, "the file is not saved by an export"
    warned = []
    monkeypatch.setattr("PySide6.QtWidgets.QMessageBox.warning", lambda *a, **k: warned.append(a))
    assert not win.export_record(b, path) and warned, "never over the file itself"
    monkeypatch.setattr("PySide6.QtWidgets.QFileDialog.getSaveFileName", lambda *a, **k: (str(out / "c.gb"), ""))
    actions = {x.text(): x for x in win.map_view.record_menu(c).actions() if x.text()}
    actions["&Export as GenBank (.gb)…"].trigger()
    app.processEvents()
    assert Document(str(out / "c.gb")).name == "copy"
    # --- remove: declined, then another record, then the open one with its edit (the next one opens) ---
    asked = []
    monkeypatch.setattr(MainWindow, "_confirm_remove", lambda self, name, edited: asked.append((name, edited)) or False)
    assert not win.remove_record(b) and names(path) == [b, c, a]
    monkeypatch.setattr(MainWindow, "_confirm_remove", lambda self, name, edited: asked.append((name, edited)) or True)
    actions = {x.text(): x for x in win.map_view.record_menu(b).actions() if x.text()}
    actions["&Remove from File…"].trigger()
    app.processEvents()
    assert names(path) == [c, a] and [n for n, _ in win.editor.fasta_records] == [c, a] and win.picker.count() == 2
    assert win.remove_record(a)
    assert asked[-1] == (a, True), "the question says the unsaved edits are lost"
    assert names(path) == [c] and win.editor.fasta_record == c and a not in win.record_docs
    assert not win.map_view._records() and not win.picker_action.isVisible(), "one record left: no list"
    assert not win.remove_record(c), "the last record stays"
    win.undo()   # the removed record's edit is not in the file's history any more
    assert win.editor.fasta_record == c


def test_fasta_records_are_reordered_exported_and_removed(tmp_path, monkeypatch):
    fa = tmp_path / "genome.fa"
    fa.write_text(">chrA\n" + "ACGT" * 30 + "\n>chrB\n" + "GGCC" * 10 + "\n>chrC\n" + "TTAA" * 5 + "\n")
    win = opened(fa)
    assert win.editor.fasta_record == "chrA" and not win.editor.multi_genbank
    assert win.reorder_records(["chrC", "chrA", "chrB"])
    assert [n for n, _ in Fasta(str(fa)).sequences()] == ["chrC", "chrA", "chrB"]
    assert win.export_record("chrB", tmp_path / "b.gb") and Document(str(tmp_path / "b.gb")).sequence() == "GGCC" * 10
    monkeypatch.setattr(MainWindow, "_confirm_remove", lambda self, name, edited: True)
    assert win.remove_record("chrA") and win.editor.fasta_record == "chrB"
    assert dict(Fasta(str(fa)).sequences()) == {"chrC": 20, "chrB": 40}
    win.editor.record_requested.emit("chrC")
    app.processEvents()
    assert win.editor.fasta_record == "chrC" and win.editor.doc.sequence() == "TTAA" * 5


def test_switching_record_does_not_read_the_file_through_again(tmp_path, monkeypatch):
    """Opening another record of the open file keeps the window's list of records (the whole file was read again to list them:
    seconds on a genome of gigabytes), and nothing is painted half-built meanwhile (the window waits, then shows the new record)."""
    path = three_records(tmp_path)
    win = opened(path)
    a, b, c = names(path)
    monkeypatch.setattr(multigb, "records", lambda *a: pytest.fail("the file was listed again"))
    monkeypatch.setattr(multigb, "is_multi_genbank", lambda *a: pytest.fail("the file was scanned again"))
    seen = []
    real = MainWindow.show_document
    monkeypatch.setattr(MainWindow, "show_document", lambda self, *a, **k: seen.append(self.updatesEnabled()) or real(self, *a, **k))
    win.editor.record_requested.emit(b)
    app.processEvents()
    assert win.editor.fasta_record == b and [n for n, _ in win.editor.fasta_records] == [a, b, c]
    assert seen == [False] and win.updatesEnabled(), "painting is held during the switch, then back on"
