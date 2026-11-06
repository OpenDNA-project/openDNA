"""A GenBank file with several records opens like a multi-FASTA: the map lists the records, one is open at a time, Save puts it back
among the others unchanged, and the search covers every record (DNA and feature names)."""
import os
import sys
import tempfile
import time
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ["OPENDNA_DATA_DIR"] = tempfile.mkdtemp(prefix="opendna_test_")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "app"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import pytest
from PySide6.QtWidgets import QApplication

import multigb
from main import MainWindow
from sample_files import PLASMID, PLASMID2

app = QApplication.instance() or QApplication([])
VIBRIO = Path(os.environ.get("OPENDNA_VIBRIO_GB", "V:/HM_Softwares/TestData/vibrio_cholerae/vibrio_cholerae_N16961.gb"))   # outside git


def two_records(tmp_path):
    if not (PLASMID.exists() and PLASMID2.exists()):
        pytest.skip("demo plasmids missing")
    path = tmp_path / "two.gb"
    path.write_text(PLASMID.read_text(encoding="utf-8").rstrip("\n") + "\n" + PLASMID2.read_text(encoding="utf-8"), encoding="utf-8")
    return path


def wait_search(win):
    sb = win.search_bar
    for _ in range(400):
        app.processEvents()
        time.sleep(0.02)
        if sb.token is not None and sb.seen >= sb.total:
            break
    app.processEvents()
    return sb


def test_records_and_saving(tmp_path, monkeypatch):
    """The records are listed in file order and switched from the map; Save puts the record back among the others; Save As keeps
    the whole file; switching records keeps the tabs and the edits in memory (one Save writes both)."""
    # --- records listed and switched ---
    path = two_records(tmp_path)
    assert multigb.is_multi_genbank(path) and not multigb.is_multi_genbank(PLASMID)
    win = MainWindow()
    win.show()
    assert win.load(str(path))
    ed = win.editor
    sizes = [l for _, l in ed.fasta_records]
    assert len(sizes) == 2 and [n for n, _ in ed.fasta_records] == [n for n, *_ in multigb.records(path)] and ed.multi_genbank
    other = ed.fasta_records[1][0]
    ed.record_requested.emit(other)
    app.processEvents()
    assert win.editor.fasta_record == other and win.editor.length == sizes[1] and win.picker.currentText() == other
    win.grab()
    # --- save puts the record back ---
    path = two_records(tmp_path)
    win = MainWindow()
    win.show()
    assert win.load(str(path))
    ed = win.editor
    first, second = ed.fasta_record, ed.fasta_records[1][0]
    before = {n: multigb.read_record(path, n).sequence() for n, _ in ed.fasta_records}
    ed.set_cursor(10)
    ed.type_bases("ACGTACGT")
    assert win.save() and not ed.doc.modified
    assert path.read_text(encoding="utf-8").count("\nLOCUS") + path.read_text(encoding="utf-8").startswith("LOCUS") == 2
    assert multigb.read_record(path, first).sequence() == before[first][:10] + "ACGTACGT" + before[first][10:]
    assert multigb.read_record(path, second).sequence() == before[second], "the other record is untouched"
    assert [n for n, *_ in multigb.records(path)].count(first) == 1 and not list(tmp_path.glob("*.saving"))
    # --- save as keeps the whole file ---
    from PySide6.QtWidgets import QFileDialog
    path = two_records(tmp_path)
    win = MainWindow()
    win.show()
    assert win.load(str(path))
    target = tmp_path / "copy.gb"
    monkeypatch.setattr(QFileDialog, "getSaveFileName", lambda *a, **k: (str(target), ""))
    assert win.save_as()
    assert multigb.is_multi_genbank(target) and multigb.is_multi_genbank(path), "both files hold both records"
    assert win.editor.multi_genbank and win.editor.path == target
    # --- switching records keeps tabs and edits ---
    path = two_records(tmp_path)
    win = MainWindow()
    win.resize(1300, 800)
    win.show()
    assert win.load(str(path))
    first = win.editor.fasta_record
    second = [n for n, _ in win.editor.fasta_records if n != first][0]
    # a layout of the user's own: one more tab than the standard ones, and the second tab current
    win.workspace.setCurrentIndex(min(1, win.workspace.count() - 1))
    layout, current = shape(win.workspace.dump()), win.workspace.currentIndex()
    before = multigb.read_record(path, first).sequence()
    win.editor.set_cursor(5)
    win.editor.type_bases("GGGG")
    assert win.editor.doc.modified
    win.editor.record_requested.emit(second)
    app.processEvents()
    assert win.editor.fasta_record == second and shape(win.workspace.dump()) == layout and win.workspace.currentIndex() == current
    assert win.record_docs[first].modified, "the first record stays in memory with its unsaved edit (no 'discard changes?' question)"
    win.editor.record_requested.emit(first)
    app.processEvents()
    assert win.editor.fasta_record == first and win.editor.doc.modified and win.editor.doc.sequence() == before[:5] + "GGGG" + before[5:]
    assert second in win.record_docs and shape(win.workspace.dump()) == layout
    # edit the second too, then one Save writes both
    win.editor.record_requested.emit(second)
    app.processEvents()
    seq2 = win.editor.doc.sequence()
    win.editor.set_cursor(7)
    win.editor.type_bases("TT")
    assert win.save()
    assert multigb.read_record(path, second).sequence() == seq2[:7] + "TT" + seq2[7:]
    assert multigb.read_record(path, first).sequence() == before[:5] + "GGGG" + before[5:], "the record switched away from was saved too"
    assert not win.record_docs[first].modified


def test_search_primers_and_detect_cover_every_record(tmp_path):
    """The search covers every record (DNA and feature names); primers, the primer search and Detect Features cover the other
    records too."""
    # --- search covers every record ---
    path = two_records(tmp_path)
    win = MainWindow()
    win.show()
    assert win.load(str(path))
    ed = win.editor
    names = [n for n, _ in ed.fasta_records]
    seq = multigb.read_record(path, names[1]).sequence()
    probe = seq[1000:1024]   # a DNA stretch that only the second (open: no) record is sure to have
    sb = win.search_bar
    sb.box.setText(probe)
    sb.run(1)
    sb = wait_search(win)
    assert {h.record for h in sb.hits} >= {names[1]}
    assert win.editor.fasta_record == names[1], "the hit is in another record: it was opened"
    # a feature name of the record that is not open
    other = names[0]
    feature = next(f for f in multigb.read_record(path, other).features() if f.kind == "CDS" and f.label)
    win.search_bar.box.setText(feature.label)
    win.search_bar.run(1)
    sb = wait_search(win)
    assert any(h.record == other and h.kind == "feature" for h in sb.hits)
    # --- primers and Detect Features ---
    from annotation import DetectDialog
    path = two_records(tmp_path)
    win = MainWindow()
    win.show()
    assert win.load(str(path))
    ed = win.editor
    other = [n for n, _ in ed.fasta_records if n != ed.fasta_record][0]
    primer = multigb.read_record(path, other).sequence()[2000:2022]
    ed.set_extra([("probe", primer)])
    for _ in range(300):
        app.processEvents()
        time.sleep(0.02)
        if ed.elsewhere:
            break
    assert any(s.record == other and s.label == "probe" and s.start == 2000 for s in ed.elsewhere)
    tree = win.primer_panel.tree
    texts = [tree.topLevelItem(i).text(4) for i in range(tree.topLevelItemCount())]
    assert any("other sequences" in t for t in texts), texts
    # the primer name is found by the search on the record that is not open
    win.search_bar.box.setText("probe")
    win.search_bar.run(1)
    sb = wait_search(win)
    assert any(h.kind == "primer" and h.record == other for h in sb.hits)
    # Detect Features runs on each record
    dlg = DetectDialog(win, win.editor)
    now = win.editor.fasta_record   # (the search opened the other record)
    assert dlg.columns[-1] == "Record" and set(dlg.others) == {n for n, _ in win.editor.fasta_records if n != now}


def test_vibrio_cholerae_two_chromosomes():
    if not VIBRIO.exists():
        pytest.skip(f"no V. cholerae file at {VIBRIO}")
    win = MainWindow()
    win.resize(1400, 850)
    win.show()
    assert win.load(str(VIBRIO))
    ed = win.editor
    assert ed.fasta_records == [("NC_002505", 2961149), ("NC_002506", 1072315)] and ed.doc.circular
    assert len(ed.features) > 1000
    win.grab()
    ed.record_requested.emit("NC_002506")
    app.processEvents()
    assert win.editor.fasta_record == "NC_002506" and win.editor.length == 1072315
    win.search_bar.box.setText("toxin")
    win.search_bar.run(1)
    sb = wait_search(win)
    assert {h.record for h in sb.hits} == {"NC_002505", "NC_002506"}


def shape(layout):
    """A workspace layout without the splitter sizes (pixels): which tabs, which panes in them."""
    if isinstance(layout, dict):
        return {k: shape(v) for k, v in layout.items() if k != "sizes"}
    if isinstance(layout, list):
        return [shape(x) for x in layout]
    return layout


