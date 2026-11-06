"""Databases menu: the Features Database and the Primer Database as tabs of the window (the Features / Primers view beside the
sequence, on a document made of the database); filter by values on the tables' headers; the reference genome box of the toolbar."""
import os
import random
import sys
import tempfile
import time
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ["OPENDNA_DATA_DIR"] = tempfile.mkdtemp(prefix="opendna_test_")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "app"))

import pytest
from PySide6.QtCore import Qt
from PySide6.QtWidgets import QApplication, QMessageBox, QTableWidget, QTableWidgetItem, QTreeWidget, QTreeWidgetItem

from opendna_core import Document

import db_tabs
import parts_library as pl
import primer_db
import reference_picker
import value_filter
from assembly_panel import ReadOnlyView
from feature_table import FeatureTable
from primer_panel import PrimerPanel
from test_cas9 import close, dna, open_doc

app = QApplication.instance() or QApplication([])

PLAIN = """#title: Cloning vectors (NCBI records)
name\ttype\tsequence\torganism\tstrength\tunits\tsource
AmpR\tCDS\tATGAGTATTCAACATTTCCGTGTCG\t\t\t\tpBR322 (J01749)
lacZ alpha\tCDS\tATGACCATGATTACGGATTCACTGG\t\t\t\tpUC19 (L09137)
"""


@pytest.fixture
def library(tmp_path, monkeypatch):
    """A library of this test's own: one published source and your own (Local) features."""
    (tmp_path / "parts").mkdir()
    (tmp_path / "parts" / "vectors.tsv").write_text(PLAIN, encoding="utf-8")
    monkeypatch.setattr(pl, "BUNDLED_DIR", tmp_path / "parts")
    monkeypatch.setattr(pl, "user_library", lambda: tmp_path / "mine.tsv")
    pl.set_on("vectors", True)
    pl.add_parts([("my promoter", "promoter", "TTGACAATTAATCATCGGCTCG", "pMine.gb")])
    yield
    pl.set_on("vectors", True)


@pytest.fixture
def win(monkeypatch):
    w = open_doc(Document.new_sequence("t", dna(500), False))
    monkeypatch.setattr(QMessageBox, "question", staticmethod(lambda *a, **k: QMessageBox.Yes))
    yield w
    close(w)


def menu_action(w, menu, text):
    m = next(a.menu() for a in w.menuBar().actions() if a.text() == menu)
    return next(a for a in m.actions() if a.text() == text)


def app_process():
    for _ in range(3):
        app.processEvents()


def test_features_database_tab(win, library):
    """A tab of this window, the Features view beside the sequence of the feature clicked (alone); one row per feature with its
    source and the file it came from; filter by Source; deleting a Local feature; the filters stay when the library changes;
    the tab moves with the window when another file is opened; reopening shows the same tab."""
    menu_action(win, "Data&bases", "&Features Database…").trigger()
    tab = win.workspace.currentWidget()
    assert isinstance(tab, db_tabs.DatabaseTab) and win.workspace.tabText(win.workspace.currentIndex()) == "Features Database"
    assert isinstance(tab.table, FeatureTable) and tab.editor is not win.editor
    assert tab.sequence is None and tab.hint.isVisibleTo(tab), "nothing clicked yet: a hint, no sequence"
    model = tab.table.model
    assert model.columns == ["Name", "Type", "Source", "Length"], "no Location; no Measured when nothing is measured"
    rows = {r[0][0]: r[0] for r in model.shown}
    assert rows["my promoter"][2] == "Local · pMine.gb", "the source and the file it came from, in one field"
    assert rows["AmpR"][2] == "Cloning vectors (NCBI records) · pBR322 (J01749)"
    # Clicking a row shows that feature alone, its bases with the feature over them.
    row = next(i for i, r in enumerate(model.shown) if r[0][0] == "my promoter")
    tab.table.view.selectRow(row)
    shown = tab.entry_editor
    assert isinstance(tab.sequence, ReadOnlyView) and shown.doc.sequence(0, None) == "TTGACAATTAATCATCGGCTCG" and shown.dna_locked
    (feature,) = shown.doc.features()
    assert (feature.kind, feature.label, feature.spans) == ("promoter", "my promoter", [(0, 22)])
    assert "my promoter" in tab.caption.text() and "pMine.gb" in tab.caption.text()
    # Filter by values: only Local (the Source column is filtered by the source's name, not by each record).
    vf = tab.view().value_filter
    source = model.columns.index("Source")
    assert dict(vf.values(source)) == {"Cloning vectors (NCBI records)": 2, pl.LOCAL: 1}
    vf.set_hidden("Source", {"Cloning vectors (NCBI records)"})
    assert [model.shown[r][0][0] for r in vf.shown_rows()] == ["my promoter"]
    assert model.headerData(source, Qt.Horizontal).endswith(value_filter.MARK)
    # The library changes (a feature added): the list is made again, the filter and the feature on view stay.
    pl.add_parts([("rbs", "RBS", "AAGGAGGTAAAAAATGAAACC", "b.gb")])
    app_process()
    model = tab.table.model
    assert sorted(model.shown[r][0][0] for r in tab.view().value_filter.shown_rows()) == ["my promoter", "rbs"]
    assert tab.entry_editor.doc.sequence(0, None) == "TTGACAATTAATCATCGGCTCG"
    # The right-click menu is the database's: only Local features can be deleted; the one on view goes, the hint comes back.
    tab.view().value_filter.clear()
    by_name = {p.name: i for i, p in tab.editor.database.by_index.items()}
    tab.editor.database.delete(tab, [by_name["my promoter"], by_name["AmpR"]])
    app_process()
    assert {p.name for p in pl.all_parts()} == {"AmpR", "lacZ alpha", "rbs"} and tab.sequence is None
    # Help: how to use the tab, and the database's status.
    help = db_tabs.help_dialog(tab, db_tabs.FEATURES)
    assert help.status.title == "Features Database" and help.windowTitle() == "Features Database Help"
    # Another file opened in the window: the tab follows, as it was.
    win.show_document(Document.new_sequence("u", dna(300), False), Path(tempfile.mkdtemp()) / "u.gb", untitled=True)
    titles = [win.workspace.tabText(i) for i in range(win.workspace.count())]
    assert "Features Database" in titles and win.workspace.currentWidget() is tab
    win.workspace.close_tab(win.workspace.indexOf(tab))
    assert db_tabs.open_tab(win, db_tabs.FEATURES) is tab, "closed and opened again: the same tab"


def test_primer_database_tab(win):
    db = primer_db.db()
    db.add("tab_fwd", "ACGTTGCAGGCTAGCATCGA", "for the tab")
    db.add("tab_amb", "ACGTNNRYACGTACGTAC", "")
    db.add("tab_long", "GGCCTTAAGGCCTTAAGGCCTTAAGG", "")
    menu_action(win, "&Primers", "Primer &Database…").trigger()
    tab = win.workspace.currentWidget()
    assert win.workspace.tabText(win.workspace.currentIndex()) == "Primer Database" and isinstance(tab.table, PrimerPanel)
    tree = tab.table.tree
    root = tree.invisibleRootItem()
    items = {root.child(i).text(0): root.child(i) for i in range(root.childCount())}
    assert set(items) == {r["name"] for r in db.rows}, "every primer, also one that binds nowhere"
    assert items["tab_fwd"].text(7) == "for the tab" and items["tab_fwd"].text(6)
    assert tab.table.columns[-2:] == ["Added", "Note"]
    tree.setCurrentItem(items["tab_fwd"])
    assert tab.entry_editor.doc.sequence(0, None) == "ACGTTGCAGGCTAGCATCGA", "the primer alone on the right"
    # Sorting: a header click sorts (Length as a number), a second click the other way; it stays after a refresh.
    tree.header().sectionClicked.emit(2)
    lengths = [int(root.child(i).text(2)) for i in range(root.childCount())]
    assert lengths == sorted(lengths)
    tree.header().sectionClicked.emit(2)
    lengths = [int(root.child(i).text(2)) for i in range(root.childCount())]
    assert lengths == sorted(lengths, reverse=True)
    tab.editor.database.delete(tab, ["tab_amb"])
    app_process()
    tree = tab.table.tree
    root = tree.invisibleRootItem()
    names = [root.child(i).text(0) for i in range(root.childCount())]
    assert "tab_amb" not in names and db.find("tab_amb") is None
    lengths = [int(root.child(i).text(2)) for i in range(root.childCount())]
    assert lengths == sorted(lengths, reverse=True), "the sort survives the list being made again"
    current = tree.currentItem()
    assert current is not None and current.text(0) == "tab_fwd", "the primer picked stays picked"
    assert db_tabs.help_dialog(tab, db_tabs.PRIMERS).status.title == "Primer Database"


def test_databases_menu(win):
    databases = next(a.menu() for a in win.menuBar().actions() if a.text() == "Data&bases")
    texts = [a.text() for a in databases.actions() if not a.isSeparator()]
    assert texts == ["&Features Database…", "&Primer Database…", "Database &Status"], "no enzyme, no reference genome"


def test_value_filter_on_a_table_and_a_tree():
    """A table refilled keeps its filter; a value never seen before is shown; the tree filters its top rows."""
    table = QTableWidget(0, 2)
    table.setHorizontalHeaderLabels(["Name", "Strand"])
    vf = value_filter.install(table, ["Strand"])

    def fill(rows):
        table.setRowCount(0)
        table.setRowCount(len(rows))
        for i, row in enumerate(rows):
            for c, text in enumerate(row):
                table.setItem(i, c, QTableWidgetItem(text))

    fill([("a", "→"), ("b", "←"), ("c", "→"), ("d", "")])
    assert vf.values(1) == [("←", 1), ("→", 2), (value_filter.EMPTY, 1)]
    vf.set_hidden("Strand", {"←", value_filter.EMPTY})
    assert vf.shown_rows() == [0, 2] and table.horizontalHeaderItem(1).text() == "Strand" + value_filter.MARK
    fill([("e", "←"), ("f", "→"), ("g", "—")])
    app_process()
    assert vf.shown_rows() == [1, 2], "refilled: still filtered, and a new value shows"
    vf.clear()
    assert vf.shown_rows() == [0, 1, 2] and table.horizontalHeaderItem(1).text() == "Strand"
    picker = value_filter.ValuePicker(vf, "Strand", 1)
    picker.search.setText("←")
    picker._tick(False)   # None acts on the values listed
    picker._ok()
    assert vf.hidden == {"Strand": {"←"}}
    tree = QTreeWidget()
    tree.setHeaderLabels(["Primer", "Match"])
    for name, match in (("p1", "1 site"), ("p2", "no site"), ("p3", "1 site")):
        item = QTreeWidgetItem([name, match])
        item.addChild(QTreeWidgetItem(["site", "exact"]))
        tree.addTopLevelItem(item)
    tvf = value_filter.install(tree, ["Match"])
    tvf.set_hidden("Match", {"no site"})
    assert tvf.shown_rows() == [0, 2] and tree.headerItem().text(1) == "Match" + value_filter.MARK


def test_reference_box_lists_the_last_three(win, tmp_path, monkeypatch):
    files = []
    for i in range(4):
        f = tmp_path / f"genome{i}.fa"
        f.write_text(">c\nACGTACGTAC\n", encoding="utf-8")
        files.append(f)
    s = primer_db.settings()
    before = (s.value("primers/genome", ""), s.value("primers/genome_check", "true"), s.value(reference_picker.KEY, []))
    try:
        picker = win.reference_picker
        for f in files:
            monkeypatch.setattr(reference_picker.QFileDialog, "getOpenFileName", staticmethod(lambda *a, f=f, **k: (str(f), "")))
            picker.browse()
        box = picker.box
        assert [box.itemText(i) for i in range(box.count())] == ["genome3.fa", "genome2.fa", "genome1.fa", reference_picker.OFF, "",
                                                                 reference_picker.ADD]
        assert box.currentIndex() == 0 and reference_picker.current() == str(files[3])
        box.activated.emit(2)
        assert reference_picker.current() == str(files[1]) and box.currentText() == "genome1.fa"
        box.activated.emit(box.findText(reference_picker.OFF))   # no background genome: the check is off
        assert reference_picker.current() is None and box.currentText() == reference_picker.OFF
    finally:
        s.setValue("primers/genome", before[0])
        s.setValue("primers/genome_check", before[1])
        s.setValue(reference_picker.KEY, before[2])
        primer_db.notifier.changed.emit()


def test_primer_sites_grouped_with_mismatches_in_red(win, tmp_path):
    """A primer's sites are grouped: on this sequence, then (once the primer is opened) in the reference genome; in each site the
    primer's bases that do not pair with the template are marked (drawn red), on this sequence and in the genome alike."""
    from primer_panel import GROUP_ROLE, MISMATCHES, mismatch_indices, primer_item
    from editor import reverse_complement
    assert mismatch_indices("ACGTNR", "ACCTGA") == [2], "N and R pair with what they stand for"
    rnd = random.Random(5)
    win.show_document(Document.new_sequence("t2", "".join(rnd.choice("ACGT") for _ in range(600)), False),
                      Path(tempfile.mkdtemp()) / "t2.gb", untitled=True)
    ed = win.editor
    target = ed.doc.sequence(100, 124)
    primer = target[:5] + ("A" if target[5] != "A" else "C") + target[6:]   # one mismatch, far from the 3' end
    genome_file = tmp_path / "genome.fa"
    bases = "".join(random.Random(11).choice("ACGT") for _ in range(400))
    rc = reverse_complement(primer)
    genome_file.write_text(f">chrT\n{bases[:150]}{rc}{bases[150:]}\n", encoding="utf-8")   # the whole primer, on the - strand
    s = primer_db.settings()
    before = (s.value("primers/genome", ""), s.value("primers/genome_check", "true"))
    try:
        reference_picker.use(genome_file)
        primer_db.db().add("red_test", primer, "")
        primer_db.notifier.changed.emit()   # (what the primer dialogs do after a change)
        app_process()
        panel = win.primer_panel
        root = panel.tree.invisibleRootItem()
        top = next(root.child(i) for i in range(root.childCount()) if root.child(i).text(0) == "red_test")
        assert panel.tree.headerItem().text(0) == "Name"
        here = top.child(0)
        assert here.data(0, GROUP_ROLE) and here.text(0) == "On this sequence (1)" and top.childCount() == 1, "genome: only once opened"
        site = here.child(0)
        assert site.data(1, MISMATCHES) == [5] and site.text(4) == "1 mismatch" and primer_item(site) is top
        top.setExpanded(True)   # opened: the genome is searched (on a worker thread)
        for _ in range(200):
            app_process()
            root = panel.tree.invisibleRootItem()
            top = next(root.child(i) for i in range(root.childCount()) if root.child(i).text(0) == "red_test")
            if top.childCount() > 1 and "searching" not in top.child(1).text(0):
                break
            time.sleep(0.02)
        genome = top.child(1)
        assert genome.text(0).startswith("On the background genome, chrT"), genome.text(0)
        hit = genome.child(0)
        assert hit.text(0).startswith("chrT:151–") and hit.text(0).endswith("←") and hit.data(1, MISMATCHES) == []
        assert hit.text(4) == "whole primer, exact"
        panel.tree.setCurrentItem(hit)   # a genome site is not on this sequence: nothing to select, no error
    finally:
        primer_db.db().remove(["red_test"])
        s.setValue("primers/genome", before[0])
        s.setValue("primers/genome_check", before[1])
        primer_db.notifier.changed.emit()
