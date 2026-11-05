"""A file of several records is one document to the user: one undo history across its records, the display settings kept when
another record is shown, and Save / Save As writing the whole file."""
import os
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ["OPENDNA_DATA_DIR"] = tempfile.mkdtemp(prefix="opendna_test_")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "app"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import pytest
from PySide6.QtWidgets import QApplication, QFileDialog

import multigb
from file_history import FileHistory
from history_panel import HistoryPanel
from main import MainWindow
from sample_files import PLASMID, PLASMID2

app = QApplication.instance() or QApplication([])


# ---- the log itself (no window) --------------------------------------------------------------------------------------------
def test_log():
    """The log orders edits across records and follows undo and redo; steps take each record's own summaries in file order."""
    # --- the log across records ---
    h = FileHistory()
    h.observe("A", 0, 0)
    h.observe("B", 0, 0)                 # starting points
    h.observe("A", 1, 0)
    h.observe("B", 1, 0)
    h.observe("A", 2, 0)
    assert h.undo_log == ["A", "B", "A"] and h.next_undo() == "A"
    h.observe("A", 1, 1)                 # undo in A
    assert h.undo_log == ["A", "B"] and h.redo_log == ["A"] and h.next_redo() == "A" and h.next_undo() == "B"
    h.observe("B", 0, 1)                 # undo in B
    assert h.undo_log == ["A"] and h.redo_log == ["A", "B"] and h.next_redo() == "B"
    h.observe("B", 1, 0)                 # redo in B
    assert h.undo_log == ["A", "B"] and h.redo_log == ["A"]
    h.observe("A", 2, 0)                 # a new edit in A: A's undone step cannot come back
    assert h.undo_log == ["A", "B", "A"] and h.redo_log == []
    h.reset()
    assert h.undo_log == h.redo_log == [] and h.next_undo() is None
    # --- steps in file order ---
    h = FileHistory()
    h.undo_log, h.redo_log = ["A", "B", "A"], ["B", "A"]
    undo, redo = h.steps({"A": (["a1", "a2"], ["a3"]), "B": (["b1"], ["b2"])})
    assert undo == [("A", "a1"), ("B", "b1"), ("A", "a2")]
    assert redo == [("A", "a3"), ("B", "b2")], "redo order: the next one first"


# ---- in the window ---------------------------------------------------------------------------------------------------------
def two_records(tmp_path):
    if not (PLASMID.exists() and PLASMID2.exists()):
        pytest.skip("demo plasmids missing")
    path = tmp_path / "two.gb"
    path.write_text(PLASMID.read_text(encoding="utf-8").rstrip("\n") + "\n" + PLASMID2.read_text(encoding="utf-8"), encoding="utf-8")
    return path


def opened(tmp_path):
    path = two_records(tmp_path)
    win = MainWindow()
    win.resize(1300, 800)
    win.show()
    assert win.load(str(path))
    first = win.editor.fasta_record
    second = [n for n, _ in win.editor.fasta_records if n != first][0]
    return win, path, first, second


def switch(win, name):
    win.editor.record_requested.emit(name)
    app.processEvents()
    assert win.editor.fasta_record == name


def test_file_as_one_document(tmp_path, monkeypatch):
    """One undo history across the records of a file (undo / redo walk back through the records in the order of the edits, the
    History panel lists the whole file), display settings and letter size kept when another record is shown, and Save As
    writing the whole file with every edit."""
    # --- undo and redo across records ---
    win, path, first, second = opened(tmp_path)
    base = {n: multigb.read_record(path, n).sequence() for n in (first, second)}
    win.editor.set_cursor(10)
    win.editor.type_bases("AAAA")          # edit 1, in the first record
    switch(win, second)
    win.editor.set_cursor(20)
    win.editor.type_bases("CCCC")          # edit 2, in the second
    switch(win, first)
    win.editor.set_cursor(30)
    win.editor.type_bases("GGGG")          # edit 3, in the first again
    fh = win.file_history
    assert fh.undo_log == [first, second, first]
    win.undo()                              # edit 3
    assert win.editor.fasta_record == first and win.editor.doc.sequence()[30:34] != "GGGG"
    win.undo()                              # edit 2: shows the second record, undone there
    assert win.editor.fasta_record == second and win.editor.doc.sequence() == base[second]
    win.undo()                              # edit 1
    assert win.editor.fasta_record == first and win.editor.doc.sequence() == base[first]
    assert fh.undo_log == [] and len(fh.redo_log) == 3
    win.redo()
    assert win.editor.fasta_record == first and win.editor.doc.sequence() == base[first][:10] + "AAAA" + base[first][10:]
    win.redo()
    assert win.editor.fasta_record == second and win.editor.doc.sequence()[20:24] == "CCCC"
    win.redo()
    assert win.editor.fasta_record == first and win.editor.doc.sequence()[30:34] == "GGGG"
    assert fh.undo_log == [first, second, first] and fh.redo_log == []
    # --- History panel ---
    win, path, first, second = opened(tmp_path)
    win.editor.set_cursor(10)
    win.editor.type_bases("AAAA")
    switch(win, second)
    win.editor.set_cursor(20)
    win.editor.type_bases("CCCC")
    panel = HistoryPanel(win.editor)
    texts = [panel.list.item(i).text() for i in range(panel.list.count())]
    mine = [i for i, t in enumerate(texts) if t.startswith(f"{first}:")]
    other = [i for i, t in enumerate(texts) if t.startswith(f"{second}:")]
    assert mine and other and mine[0] < other[0], texts
    panel.list.setCurrentRow(1)             # the first row after "Opened": back to before the first edit
    panel.go()
    app.processEvents()
    app.processEvents()
    assert win.editor.fasta_record == first and win.file_history.undo_log == []
    assert win.record_docs[second].history()[0] == []
    # --- display settings and letter size ---
    win, path, first, second = opened(tmp_path)
    ed = win.editor
    ed.set_display(enzymes=True, hidden_kinds={"promoter"}, collapse_genes=False, cds_only=False, translation=False)
    win.sequence_view.set_size(14)
    win.set_dna_locked(True)
    switch(win, second)
    new = win.editor
    assert new is not ed
    assert new.show_enzymes and new.hidden_kinds == {"promoter"} and not new.collapse_genes and not new.cds_only
    assert not new.show_translation and new.dna_locked and win.sequence_view.size == 14
    # --- Save As ---
    win, path, first, second = opened(tmp_path)
    before = {n: multigb.read_record(path, n).sequence() for n in (first, second)}
    win.editor.set_cursor(10)
    win.editor.type_bases("AAAA")
    switch(win, second)                      # the first stays in memory, edited
    win.editor.set_cursor(20)
    win.editor.type_bases("CCCC")
    target = tmp_path / "copy.gb"
    monkeypatch.setattr(QFileDialog, "getSaveFileName", lambda *a, **k: (str(target), ""))
    assert win.save_as()
    assert len(multigb.records(target)) == 2
    assert multigb.read_record(target, first).sequence() == before[first][:10] + "AAAA" + before[first][10:]
    assert multigb.read_record(target, second).sequence() == before[second][:20] + "CCCC" + before[second][20:]
    assert multigb.read_record(path, first).sequence() == before[first], "the original file is untouched"
    assert win.editor.path == target and win.editor.multi_genbank and not win.editor.doc.modified
    assert not any(d.modified for d in win.record_docs.values())




def test_theme_switch_repaints_a_multi_record_window_and_keeps_the_file_history(tmp_path):
    import theme
    win, path, first, second = opened(tmp_path)
    win.editor.set_cursor(10)
    win.editor.type_bases("AAAA")
    switch(win, second)
    win.editor.set_cursor(20)
    win.editor.type_bases("CCCC")
    panel = HistoryPanel(win.editor)
    rows = panel.list.count()
    theme.refresh_views()                    # what a Light / Dark switch does to every widget
    win.grab()
    assert panel.list.count() == rows and win.file_history.undo_log == [first, second]
