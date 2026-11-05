"""Headless GUI tests (Qt offscreen) on the real files in the repo root; skipped if absent."""
import os
import sys
import tempfile
import time
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
# Keep the tests' primer database away from the real one.
os.environ["GENOMEEDIT_DATA_DIR"] = tempfile.mkdtemp(prefix="genomeedit_test_")
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "app"))

import pytest
from PySide6.QtCore import QEvent, QPoint, Qt
from PySide6.QtTest import QTest
from PySide6.QtWidgets import QAbstractSlider, QApplication, QDialog, QFileDialog, QLabel, QMessageBox

from PySide6.QtGui import QColor

from editor import feature_bounds, reverse_complement, text_color_on
from feature_dialog import FeatureDialog
from genomeedit_core import Document, Fasta
from main import MainWindow, autosave_path
from primer_db import db
from primer_design import DesignDialog
from primer_panel import structure
from local_align import INVERTED_TEXT, LocalAlignDialog
from workspace import Pane
from tools import DigestDialog, OrfDialog, PcrDialog, fragments, pcr_product
from primer_panel import (OptionsDialog, PrimerDialog, binding_report, import_primers,
                          melting_temp, parse_primer_list)

from sample_files import GENOME as GENOME_DNA, PLASMID as PLASMID_DNA, DNA_FILE   # copies (see sample_files.py)
SHOTS = os.environ.get("GE_SHOTS")  # folder to save screenshots into, for eyeballing


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


@pytest.fixture
def errors(monkeypatch):
    """Qt swallows exceptions raised in paint/event handlers; collect them instead."""
    caught = []
    monkeypatch.setattr(sys, "excepthook", lambda *exc: caught.append(exc))
    yield caught
    assert not caught, caught


def open_window(app, path):
    if not path.exists():
        pytest.skip(f"missing {path.name}")
    win = MainWindow()
    win.resize(1400, 850)
    win.show()
    t = time.perf_counter()
    win.load(str(path))
    app.processEvents()
    win.grab()  # forces a full paint of every view
    print(f"\n  {path.name}: loaded + first paint in {time.perf_counter() - t:.2f} s")
    return win


def shot(win, name):
    if SHOTS:
        win.grab().save(str(Path(SHOTS) / f"{name}.png"))


def test_plasmid_edit_features_find_and_autosave(app, errors, tmp_path, monkeypatch):
    """Plasmid: select a table row, click / type / backspace / undo in the sequence view, Shift+Right and the map click, saving a .dna goes to GenBank; then a new feature through the dialog (a bad location changes nothing), Delete in the table, reverse-complement twice, Find, translation rows, autosave recovery."""
    # --- edit flow, selection, typing and saving a .dna ---
    win = open_window(app, PLASMID_DNA)
    ed, seq_view = win.editor, win.sequence_view
    shot(win, "plasmid")
    n = ed.length

    # Selecting a row in the table selects that feature everywhere.
    table = win.feature_table.view
    table.selectRow(3)
    app.processEvents()
    assert ed.current_feature is not None
    assert (ed.start, ed.end) == feature_bounds(ed.by_index[ed.current_feature])
    shot(win, "plasmid_feature_selected")

    # Click in the sequence (top row), type bases, backspace, undo.
    seq_view.verticalScrollBar().setValue(0)
    ed.set_display(primers=False)   # a primer site crosses the origin, drawn over the top row
    app.processEvents()
    QTest.mouseClick(seq_view.viewport(), Qt.LeftButton, pos=QPoint(int(seq_view.left + 5 * seq_view.cw), 5))
    assert ed.start == ed.end == 5
    ed.set_display(primers=True)
    QTest.keyClicks(seq_view, "gattaca")
    assert ed.length == n + 7 and ed.doc.sequence(5, 12) == "gattaca"
    assert ed.doc.modified and win.windowTitle().startswith("*")
    QTest.keyClick(seq_view, Qt.Key_Backspace)
    assert ed.length == n + 6
    QTest.keyClicks(seq_view, "qzx!1")  # not DNA letters (IUPAC): ignored
    assert ed.length == n + 6
    shot(win, "plasmid_edited")
    while ed.undo():
        pass
    assert ed.length == n and not ed.doc.modified

    # Shift+Right extends the selection from the cursor.
    ed.set_cursor(10)
    seq_view.anchor = 10
    for _ in range(4):
        QTest.keyClick(seq_view, Qt.Key_Right, Qt.ShiftModifier)
    assert (ed.start, ed.end) == (10, 14)

    # Clicking the map selects something; feature lanes or the cursor.
    QTest.mouseClick(win.map_view, Qt.LeftButton, pos=win.map_view.rect().center() + QPoint(0, -60))
    win.map_view.grab()

    # Saving a .dna never overwrites it: goes to a GenBank file instead.
    dna = open_window(app, DNA_FILE)
    original = DNA_FILE.read_bytes()
    dna.editor.type_bases("A")
    out = tmp_path / "saved.gb"
    monkeypatch.setattr(QFileDialog, "getSaveFileName", lambda *a, **k: (str(out), ""))
    assert dna.save()
    assert out.exists() and DNA_FILE.read_bytes() == original
    assert not dna.editor.doc.modified and dna.editor.path == out

    # --- features, find, translation and autosave ---
    if not PLASMID_DNA.exists():
        pytest.skip(f"missing {PLASMID_DNA.name}")
    work = tmp_path / "plasmid.dna"  # a copy: autosave files land next to it
    work.write_bytes(PLASMID_DNA.read_bytes())
    win = open_window(app, work)
    ed, seq_view = win.editor, win.sequence_view
    n_features = len(ed.features)

    # Rows are only as tall as their content.
    heights = {seq_view.offset(r + 1) - seq_view.offset(r) for r in range(seq_view.row_count())}
    assert all(seq_view.row_at(seq_view.offset(r) + 1) == r for r in range(seq_view.row_count()))
    assert len(heights) > 1, heights
    assert ed.unique_sites, "expected some unique cutters on a plasmid"

    # New feature through the dialog (filled in programmatically, then OK).
    ed.select(9, 30)
    dialog = FeatureDialog(win, ed)
    dialog.name.setText("my site")
    dialog.direction.setCurrentIndex(1)   # Reverse
    dialog._direction_chosen(1)
    dialog.qualifiers.setPlainText('/note="added in test"')
    dialog.accept()
    f = ed.by_index[ed.current_feature]
    assert (f.label, f.location, f.strand) == ("my site", "complement(10..30)", -1)
    assert ("note", "added in test") in f.qualifiers
    shot(win, "plasmid_new_feature")

    # A bad location keeps the dialog open and changes nothing.
    monkeypatch.setattr(QMessageBox, "warning", lambda *a, **k: None)
    bad = FeatureDialog(win, ed, f.index)
    bad.location.setText("5..999999")
    bad.accept()
    assert len(ed.features) == n_features + 1 and ed.by_index[f.index].location == "complement(10..30)"

    # Delete from the table with the Delete key; undo brings it back.
    table = win.feature_table
    table.view.selectRow(table.model.row_of(f.index))
    # the Delete shortcut needs the table focused in the active window; other tests leave windows open, and activateWindow()
    # is asynchronous: wait until this window is the active one
    win.raise_()
    win.activateWindow()
    QTest.qWaitForWindowActive(win, 2000)   # (best effort: the Windows test machine does not let a window become active)
    table.view.setFocus()
    QTest.keyClick(table.view, Qt.Key_Delete)
    assert len(ed.features) == n_features
    ed.undo()
    assert len(ed.features) == n_features + 1

    # Reverse-complement twice = same sequence and features. Compared by bases covered:
    # gb-io merges touching join parts, e.g. join(1804..1882,1883) -> 1804..1883.
    def snapshot():
        bases = lambda x: tuple(sorted({b for s, e in x.spans for b in range(s, e)}))
        return ed.doc.sequence(), sorted((x.kind, bases(x), x.strand) for x in ed.features)
    before = snapshot()
    ed.reverse_complement()
    assert ed.doc.sequence()[:20] == reverse_complement(before[0][-20:])
    ed.reverse_complement()
    assert snapshot() == before

    # Find: DNA on either strand, then a feature name.
    probe = ed.doc.sequence(500, 515)
    win.find_box.setText(reverse_complement(probe))
    ed.set_cursor(0)
    win.find_next()
    assert (ed.start, ed.end) == (500, 515)
    win.find_box.setText("my sit")
    win.find_next()
    assert ed.current_feature == f.index

    # Translation rows paint without errors and make rows taller.
    total = seq_view.offset(seq_view.row_count())
    win.translation_action.trigger()
    assert seq_view.offset(seq_view.row_count()) > total
    seq_view.viewport().grab()
    shot(win, "plasmid_translation")

    # Autosave writes a recovery copy; reopening offers to recover it.
    win.autosave()
    recovery = autosave_path(work)
    assert recovery.exists() and work.read_bytes() == PLASMID_DNA.read_bytes()
    ed.doc.modified = False  # pretend the user chose "Discard" ...
    win.autosave_file = None  # ... after a crash: the recovery file stays
    monkeypatch.setattr(MainWindow, "_ask_recover", lambda *a: True)   # "Recover Changes"
    win.load(str(work))
    assert win.editor.doc.modified and len(win.editor.features) == n_features + 1
    win.editor.doc.modified = False
    win.close()
    assert not recovery.exists()


def test_primers_clipboard_design_pcr_orfs_and_digest(app, errors, monkeypatch, tmp_path):
    """Clipboard copy / paste between windows and from text, typing, the primer database (dialog, sites, panel tab, import / export, binding report, checks, saving as primer_bind), then primer design, PCR product, ORFs, digest and secondary structure."""
    # --- clipboard, typing and primers ---
    win = open_window(app, PLASMID_DNA)
    ed, seq_view = win.editor, win.sequence_view
    n, n_features = ed.length, len(ed.features)
    clipboard = QApplication.clipboard()

    # Copy a feature, paste it at the start: the bases and the annotation come along.
    lacz = next(f for f in ed.features if f.kind == "CDS")
    ed.select_feature(lacz.index)
    size = ed.end - ed.start
    seq_view.setFocus()
    QTest.keyClick(seq_view, Qt.Key_C, Qt.ControlModifier)
    ed.set_cursor(0)
    QTest.keyClick(seq_view, Qt.Key_V, Qt.ControlModifier)
    assert ed.length == n + size
    copies = [f for f in ed.features if f.label == lacz.label]
    assert len(copies) == 2 and any(f.spans == [(0, size)] for f in copies)
    ed.undo()
    assert ed.length == n and len(ed.features) == n_features

    # Into another window, even of another running app (the in-memory copy is gone): the
    # features travel on the system clipboard.
    from editor import Editor
    other = MainWindow()
    assert other.load(str(PLASMID_DNA))
    Editor.clip = None
    other.editor.set_cursor(0)
    other.sequence_view.setFocus()
    QTest.keyClick(other.sequence_view, Qt.Key_V, Qt.ControlModifier)
    assert sum(f.label == lacz.label for f in other.editor.features) == 2
    other.editor.doc.modified = False
    other.close()

    # Text from elsewhere: FASTA/number formatting stripped, non-DNA refused.
    clipboard.setText(">x\n1 acgt acgt\n")
    ed.set_cursor(0)
    ed.paste(clipboard.text())
    assert ed.doc.sequence(0, 8) == "acgtacgt"
    ed.undo()
    assert not ed.paste("hello world") and ed.length == n

    # Typing: only ACGTN; ambiguity codes like R are refused.
    QTest.keyClicks(seq_view, "acgtnr")
    assert ed.length == n + 5
    while ed.undo():
        pass

    # Label text color adapts to the feature color.
    assert text_color_on(QColor("#1a237e")) == QColor("white")
    assert text_color_on(QColor("#ffd54f")) == QColor("black")

    # A primer is a name + sequence in the database. The dialog shows live where it binds.
    ed.select(2000, 2022)
    dialog = PrimerDialog(win, ed)
    dialog.name.setText("testF")
    dialog.sequence.setText("gaattc" + ed.doc.sequence(2000, 2022))  # with a 5' tail
    dialog._update_report()
    assert "2,001–2,022 → exact" in dialog.report.text() and "Non-specific (sequence)" in dialog.report.text()
    dialog.accept()
    assert db().find("testF")["sequence"].startswith("gaattc")

    def site_of(name):
        return next(s for s in ed.primer_sites if s.label == name)
    f = site_of("testF")
    assert (f.spans, f.strand, f.mismatches, f.tail) == ([(2000, 2022)], 1, 0, 6)
    assert all(x.kind != "primer_bind" for x in ed.features)
    assert any(s.label == "oAC1345" for s in ed.primer_sites)  # file primers went to the DB

    # Sites follow the sequence: 5 bases typed upstream move the site.
    ed.set_cursor(100)
    ed.type_bases("AAAAA")
    assert site_of("testF").spans == [(2005, 2027)]
    ed.undo()

    # "Almost binds": one mismatch in the 3' seed (not in the last 3 bases) is shown,
    # with the mismatched base marked.
    template = list(ed.doc.sequence(2300, 2325).upper())
    template[8] = "A" if template[8] != "A" else "C"
    assert db().add("testR", reverse_complement("".join(template))) == "added"
    r = site_of("testR")
    assert (r.spans, r.strand, r.mismatches, r.mismatch_at) == ([(2300, 2325)], -1, 1, [2308])

    # Primers tab: one group per primer, collapsed by default; a site row selects it.
    panel = win.primer_panel
    groups = {panel.tree.topLevelItem(i).text(0): panel.tree.topLevelItem(i)
              for i in range(panel.tree.topLevelItemCount())}
    g = groups["testF"]
    # Its sites are grouped by where they are: here first. (At the default of up to 2 mismatches a chance near-copy elsewhere in
    # the plasmid may add a row: only the exact one is pinned.)
    here = g.child(0)
    assert here.text(0).startswith("On this sequence")
    exact_rows = [here.child(i) for i in range(here.childCount()) if here.child(i).text(4) == "exact"]
    assert not g.isExpanded() and len(exact_rows) == 1 and here.child(0).text(4) == "exact"
    assert float(here.child(0).text(3)) == pytest.approx(melting_temp(ed.doc.sequence(2000, 2022)), abs=0.05)
    panel.tree.setCurrentItem(here.child(0))
    assert ed.current_primer == f.key and (ed.start, ed.end) == (2000, 2022)
    assert "Tm" in win.status.text()
    g.setExpanded(True)
    panel.reload()  # refreshed after an edit: stays open
    assert {panel.tree.topLevelItem(i).text(0): panel.tree.topLevelItem(i)
            for i in range(panel.tree.topLevelItemCount())}["testF"].isExpanded()

    # Clicking on the primer line (above the bases) selects that site.
    ed.set_cursor(0)
    row = 2000 // seq_view.bases_per_line()
    seq_view.verticalScrollBar().setValue(seq_view.offset(row))
    lane_top = seq_view._primer_lane_top(row, 1, ed.primer_lanes[f.key])
    x = seq_view.left + (2010 - row * seq_view.bases_per_line() + 0.5) * seq_view.cw
    QTest.mouseClick(seq_view.viewport(), Qt.LeftButton, pos=QPoint(int(x), int(lane_top + seq_view.primer_h - 4)))
    assert ed.current_primer == f.key
    ed.select(2290, 2330)
    shot(win, "plasmid_primer_lines")

    # Import / export the database.
    lst = tmp_path / "primers.csv"
    lst.write_text("name,sequence\ntestR2," + reverse_complement(ed.doc.sequence(2600, 2625)) +
                   "\njunk,ACGTACGTACGTACGTACGT\n")
    assert "2 new" in import_primers(str(lst))
    assert site_of("testR2").spans == [(2600, 2625)] and db().find("junk")
    assert not any(s.label == "junk" for s in ed.primer_sites)
    out = tmp_path / "export.fasta"
    db().export(out)
    assert ">testF\ngaattc" in out.read_text()
    assert parse_primer_list(out.read_text())[-1][0] == "junk"

    # The Primer Database tab: every primer of the database, in the Primers view on a document of its own.
    import db_tabs
    tab = db_tabs.open_tab(win, db_tabs.PRIMERS)
    root = tab.table.tree.invisibleRootItem()
    names = [root.child(i).text(0) for i in range(root.childCount())]
    assert {"testF", "testR", "junk", "oAC1345"} <= set(names) and len(names) == len(db().rows)
    assert win.workspace.currentWidget() is tab
    OptionsDialog(win)

    # Binding report: one exact site; the genome line when the BW25113 file is present.
    report = binding_report(ed, ed.doc.sequence(3000, 3020))
    assert "Binds here: 3,001–3,020 → exact" in report
    if GENOME_DNA.exists():
        assert "Non-specific (genome): BW25113_2026" in report

    # Checks with marks: a primer whose whole sequence occurs twice here is red (✖) for the
    # sequence; a unique one is green (✔); its 3′ end alone elsewhere would be yellow (?).
    from primer_panel import primer_checks
    unique = dict((c, lvl) for c, lvl, _ in primer_checks(ed, ed.doc.sequence(3000, 3020)))
    assert unique["Non-specific (sequence)"] == "ok" and unique["Binds here"] == "ok"
    ed.set_cursor(10)
    ed.type_bases(ed.doc.sequence(3000, 3020))   # a second copy of the same 20 bases
    twice = dict((c, lvl) for c, lvl, _ in primer_checks(ed, ed.doc.sequence(3020, 3040)))
    assert twice["Non-specific (sequence)"] == "bad"
    ed.undo()
    # The limits are settings (Tools → Settings → Primer checks).
    import tool_settings
    from primer_panel import CHECKS
    probe = ed.doc.sequence(3000, 3020)
    t = melting_temp(probe)
    tool_settings.set_value(CHECKS, "tm_low", t + 1)
    tool_settings.set_value(CHECKS, "tm_margin", 0.5)
    assert dict((c, lvl) for c, lvl, _ in primer_checks(ed, probe))["Tm"] == "bad"
    tool_settings.set_value(CHECKS, "tm_low", 52.0)
    tool_settings.set_value(CHECKS, "tm_margin", 7.0)

    # Saving writes the current exact sites as primer_bind features; a primer stored in
    # the file whose name now means another sequence in the database is kept as it was.
    db().remove(["oAC1346"])
    db().add("oAC1346", "ACGTTTGACCAGTTTACGGA")  # same name, different primer
    out = tmp_path / "with_primers.gb"
    ed.save(out)
    saved = [f for f in Document(str(out)).features() if f.kind == "primer_bind"]
    by_label = {}
    for x in saved:
        by_label.setdefault(x.label, []).append(x)
    assert [x.location for x in by_label["testF"]] == ["2001..2022"]
    assert ("note", "sequence: gaattc" + ed.doc.sequence(2000, 2022)) in by_label["testF"][0].qualifiers
    assert "testR" not in by_label  # 1 mismatch: shown in the app, not written
    assert len(by_label["oAC1346"]) == 1
    assert any(v and "GTAGAGGATCTGCTCATG" in v for _, v in by_label["oAC1346"][0].qualifiers)
    original = next(v for k, v in by_label["oAC1346"][0].qualifiers if k == "note").split("sequence: ")[1]
    db().remove(["oAC1346"])
    db().add("oAC1346", original)   # the database is shared by the tests of this file: put the real primer back
    assert not ed.doc.modified
    while ed.undo():
        pass

    # --- primer design, PCR product, ORFs, digest and structure ---
    win = open_window(app, PLASMID_DNA)
    ed = win.editor
    seq = ed.doc.sequence().upper()

    # Design a pair around 1,001-1,200: primers outside the target, Tm near 60 °C.
    t = time.perf_counter()
    ed.select(1000, 1200)
    dialog = DesignDialog(win, ed)
    dialog.run()
    print(f"\n  primer design: {time.perf_counter() - t:.2f} s, best penalty {dialog.pairs[0]['penalty']:.1f}")
    best = dialog.pairs[0]
    f, r = best["fwd"], best["rev"]
    assert f["end"] <= 1000 and r["start"] >= 1200 and abs(f["tm"] - 60) < 5 and abs(r["tm"] - 60) < 5
    assert f["seq"] == seq[f["start"]:f["end"]] and r["seq"] == reverse_complement(seq[r["start"]:r["end"]])
    # The selected pair is drawn as temporary primers (not in the database); saving makes them regular under their new names.
    temp = {s.label: s for s in ed.primer_sites if s.temporary}
    assert set(temp) == {"pair1_F", "pair1_R"} and (temp["pair1_F"].start, temp["pair1_R"].end) == (f["start"], r["end"])
    assert "pair1_F" not in [row["name"] for row in db().rows]
    dialog.table.selectRow(1)
    assert {s.label for s in ed.primer_sites if s.temporary} == {"pair2_F", "pair2_R"}
    dialog.table.selectRow(0)
    assert dialog.save("des") == ["des_F", "des_R"]
    assert not [s for s in ed.primer_sites if s.temporary] and not ed.extra
    sites = {s.label: s for s in ed.primer_sites}
    assert (sites["des_F"].start, sites["des_R"].end) == (f["start"], r["end"])

    # Right-click a selection > Design Primers for This Sequence: its own window, closed once a pair is saved.
    from primer_design import design_for_selection
    ed.select(1000, 1200)
    window = design_for_selection(win, ed)
    assert window.isVisible() and window.pairs and ed.extra and any(s.temporary for s in ed.primer_sites)
    window.table.selectRow(0)
    assert window.save("des2") == ["des2_F", "des2_R"] and not window.isVisible() and not ed.extra
    ed.select(1000, 1200)
    window = design_for_selection(win, ed)
    window.close()   # closed without saving: the pairs it drew go with it
    assert not ed.extra

    # PCR product from the designed pair opens in a new, untitled window.
    pcr = PcrDialog(win, ed)
    pcr.fwd.setCurrentIndex(pcr.fwd_sites.index(sites["des_F"]))
    pcr.rev.setCurrentIndex(pcr.rev_sites.index(sites["des_R"]))
    monkeypatch.setattr(PcrDialog, "exec", lambda self: QDialog.Accepted)
    monkeypatch.setattr(PcrDialog, "_pick_default", lambda self: None)
    monkeypatch.setattr(PcrDialog, "pair", lambda self: (sites["des_F"], sites["des_R"]))
    product_win = win.pcr_product()
    prod = product_win.editor
    assert prod.doc.sequence().upper() == seq[f["start"]:r["end"]] and prod.untitled and prod.doc.modified
    assert prod.features, "features in the amplified region come along"
    out = tmp_path / "product.gb"
    monkeypatch.setattr(QFileDialog, "getSaveFileName", lambda *a, **k: (str(out), ""))
    assert product_win.save() and out.exists() and not prod.untitled
    product_win.editor.doc.modified = False
    product_win.close()

    # A primer with a mismatch mutates the product (the primer's base wins).
    mutant = list(seq[1500:1525])
    mutant[12] = "A" if mutant[12] != "A" else "C"
    db().add("mutF", "".join(mutant))
    db().add("mutR", reverse_complement(seq[1700:1725]))
    sites = {s.label: s for s in ed.primer_sites}
    assert sites["mutF"].mismatches == 1
    product = pcr_product(ed, sites["mutF"], sites["mutR"])
    assert product.sequence(0, 25).upper() == "".join(mutant) and product.length == 225

    # ORFs: the plasmid's CDS features that start with ATG are found as ORFs.
    orfs = OrfDialog(win, ed)
    found = set(orfs.orfs)
    cds = [(x.spans[0][0], x.spans[0][1], x.strand) for x in ed.features if x.kind == "CDS" and len(x.spans) == 1
           and (x.spans[0][1] - x.spans[0][0]) // 3 > 100]
    assert any(c in found for c in cds), (cds, sorted(found)[:10])
    n_features = len(ed.features)
    orfs.table.selectRow(orfs.links.index(None))  # an ORF not annotated yet (annotated ones are skipped)
    orfs.add()
    assert len(ed.features) == n_features + 1 and ed.features[-1].label.startswith("ORF")
    ed.undo()

    # Digest: fragment sizes add up to the plasmid; the gel has ladder + digest lanes.
    digest = DigestDialog(win, ed)
    for name in ("BamHI", "PvuI"):   # one site each
        digest.enzymes.item_of(name).setCheckState(0, Qt.Checked)
    assert len(digest.frags) == 2 and sum(fr[0] for fr in digest.frags) == ed.length
    assert len(digest.gel.lanes) == 2
    digest.gel.resize(260, 400)
    digest.gel.grab()
    assert fragments(100, False, {30: ["X"]}) == [(30, 0, 30, [], ["X"]), (70, 30, 100, ["X"], [])]

    # Secondary structure warnings; the map draws primer sites.
    assert any("3′ self-dimer" in w for w in structure("ACGTTTAGCGAATTCGC")["warnings"])
    assert structure("AGCGGATAACAATTTCACACAGGA")["warnings"] == []
    win.map_view.grab()
    shot(win, "plasmid_design")


def test_fasta_file_tracks_genome_view_search_and_chromosome_list(app, errors, tmp_path, monkeypatch):
    """A multi-sequence FASTA: edit and save, features offer GenBank, picker, GFF / BED tracks in the genome view, track display, legend, search in the whole file; the map's chromosome list switches sequence; chromosomes are listed largest first."""
    # --- FASTA, tracks, genome view and search ---
    import random
    rng = random.Random(7)
    chr_a = "".join(rng.choice("ACGT") for _ in range(150_000))
    chr_b = "".join(rng.choice("ACGT") for _ in range(20_000)) + "GGTACCTTGGTACC"
    fa = tmp_path / "genome.fa"
    fa.write_text("".join(f">{n} test\n" + "\n".join(s[i:i + 60] for i in range(0, len(s), 60)) + "\n"
                          for n, s in (("chrA", chr_a), ("chrB", chr_b))))
    win = MainWindow()
    win.resize(1400, 850)
    win.show()
    assert win.load(str(fa))
    ed, gv = win.editor, win.genome_view
    assert ed.fasta_record == "chrA" and ed.length == 150_000
    assert win.picker_action.isVisible() and win.picker.count() == 2
    assert not gv.isHidden()  # > 100 kb: genome view opens by itself
    gv.canvas.grab()

    # Edit and save: the FASTA is rewritten with chrB untouched.
    ed.set_cursor(10)
    ed.type_bases("NNNN")
    assert win.save() and not ed.doc.modified
    back = Fasta(str(fa))
    assert back.fetch("chrA", 0, 20) == chr_a[:10] + "NNNN" + chr_a[10:16]
    assert back.fetch("chrB", 0, 30_000) == chr_b

    # A feature cannot live in a FASTA: saving offers GenBank instead.
    ed.select(100, 200)
    ed.add_features([("gene", "101..200", [("label", "x")])])
    shown = []
    monkeypatch.setattr(QMessageBox, "information", lambda *a, **k: shown.append(a))
    monkeypatch.setattr(QFileDialog, "getSaveFileName", lambda *a, **k: (str(tmp_path / "chrA.gb"), ""))
    assert win.save() and shown and (tmp_path / "chrA.gb").exists() and ed.fasta_record is None

    # Switch sequence with the picker.
    assert win.load(str(fa), "chrB") and win.editor.fasta_record == "chrB"
    assert win.load(str(fa))
    ed, gv = win.editor, win.genome_view

    # GFF track: drawn in the genome view, its names found by the locus box.
    gff = tmp_path / "ann.gff3"
    gff.write_text("##gff-version 3\n" + "".join(
        f"chrA\tt\tgene\t{s + 1}\t{s + 800}\t.\t{'+' if i % 2 else '-'}\t.\tID=g{i};Name=gene{i}\n"
        for i, s in enumerate(range(0, 140_000, 1000))) + "chrB\tt\tgene\t1\t50\t.\t+\t.\tName=other\n")
    track = win.load_track(str(gff))
    assert track.seqid == "chrA" and ed.tracks == [track]
    gv.show_all()
    gv.canvas.grab()  # 140 items: boxes; zoomed out on more would switch to density
    gv.locus.setText("gene42")
    gv.go()
    assert (ed.start, ed.end) == (42_000, 42_800) and gv.vstart <= 42_000 and gv.vend >= 42_800
    gv.locus.setText("1,000-2,000")
    gv.go()
    assert (ed.start, ed.end) == (999, 2000)
    gv.zoom(0.001)  # deepest zoom: bases drawn
    assert any(name == "Bases" for name, _, _ in gv.canvas.layout_rows()[0])
    gv.canvas.grab()
    shot(win, "genome_view_fasta")

    # Clicking a track item selects its range.
    gv.locus.setText("gene10")
    gv.go()
    gv.canvas.grab()
    rect, action = next((r, a) for r, a in gv.canvas.hits if a == ("range", 10_000, 10_800))
    ed.set_cursor(0)
    QTest.mouseClick(gv.canvas, Qt.LeftButton, pos=rect.center().toPoint())
    assert (ed.start, ed.end) == (10_000, 10_800)

    # Track display: click its name to hide / show, condensed is thinner, drag the line under
    # it to resize, log scale for bars. Several files at once from File → Load Tracks.
    from genome_view import GAP, HIDDEN_H
    c = gv.canvas

    def row(name):
        return next((top, h) for n, top, h in c.layout_rows()[0] if n == name)

    top, h = row(track.name)
    QTest.mouseClick(c, Qt.LeftButton, pos=QPoint(10, int(top) + 5))
    assert c.style(track.name).mode == "hidden" and row(track.name)[1] == HIDDEN_H
    c.grab()
    QTest.mouseClick(c, Qt.LeftButton, pos=QPoint(10, int(top) + 5))
    assert row(track.name)[1] == h
    c.set_style(track.name, mode="condensed")
    top, small = row(track.name)
    assert small < h
    y = int(top + small + GAP / 2)
    assert c.edge_at(y)[0] == track.name
    QTest.mousePress(c, Qt.LeftButton, pos=QPoint(400, y))
    QTest.mouseMove(c, QPoint(400, y + 40))
    QTest.mouseRelease(c, Qt.LeftButton, pos=QPoint(400, y + 40))
    assert row(track.name)[1] == small + 40 and c.minimumHeight() >= c.layout_rows()[1]
    c.set_style(track.name, scale="log10")
    gv.show_all()
    c.grab()
    bed = tmp_path / "more.bed"
    bed.write_text("chrA\t100\t200\tb1\n")
    monkeypatch.setattr(QFileDialog, "getOpenFileNames", lambda *a, **k: ([str(gff), str(bed)], ""))
    assert len(win.load_tracks()) == 2 and len(ed.tracks) == 3
    # Moving the cursor / selection redraws only the overlay: the tracks image is reused.
    import perf
    c.grab()
    perf.STATS.clear()
    ed.select(5, 50)
    c.grab()
    assert "genome view: tracks (total)" not in perf.STATS and perf.STATS["genome view: cursor + selection (on the cached tracks)"]
    win.show_timings()
    assert win.timings.isVisible() and win.timings.table.rowCount() >= 1
    win.timings.close()

    # Legend: its own window, toggled by the button or closed by its close box.
    gv.legend.click()
    window = gv.legend_window
    assert window.isVisible() and "mate unmapped" in window.label.text() and "80 %" in window.label.text()
    window.close()
    assert not gv.legend.isChecked() and not window.isVisible()
    gv.legend.click()
    assert window.isVisible()
    gv.legend.click()
    assert not window.isVisible()

    # Search in the whole file: hits stream in per sequence; a hit in chrB opens chrB.
    bar = win.search_bar
    bar.scopes["file"].setChecked(True)
    bar.box.setText("GGTACC")
    bar.step(1)
    t = time.perf_counter()
    while bar.seen < bar.total and time.perf_counter() - t < 10:
        QTest.qWait(20)
    expected = sum(s.count("GGTACC") for s in (chr_a, chr_b))  # palindrome: one strand
    assert len(bar.hits) == expected and bar.seen == bar.total == 2
    bar._go(next(i for i, h in enumerate(bar.hits) if (h.record, h.start) == ("chrB", 20_000)))
    assert win.editor.fasta_record == "chrB" and (win.editor.start, win.editor.end) == (20_000, 20_006)
    assert bar.counter.text().endswith(f"of {expected}") and win.editor.search_spans   # the hits of chrB are highlighted
    bar.scopes["sequence"].setChecked(True)
    win.editor.doc.modified = False
    win.close()

    # --- the map's list of chromosomes switches sequence ---
    fa = tmp_path / "two.fa"
    fa.write_text(">chrA\n" + "ACGT" * 500 + "\n>chrB\n" + "TTGCA" * 300 + "\n")
    win = MainWindow()
    win.resize(1200, 800)
    win.show()
    assert win.load(str(fa))
    assert win.editor.fasta_records == [("chrA", 2000), ("chrB", 1500)]
    win.editor.record_requested.emit("chrB")
    app.processEvents()
    assert win.editor.fasta_record == "chrB" and win.editor.length == 1500
    assert win.editor.fasta_records == [("chrA", 2000), ("chrB", 1500)] and win.picker.currentText() == "chrB"
    win.grab()

    # --- the chromosome list and picker are largest first ---
    fa = tmp_path / "unsorted.fa"
    fa.write_text(">small\n" + "ACGT" * 50 + "\n>big\n" + "ACGT" * 500 + "\n>mid\n" + "ACGT" * 200 + "\n")
    win = MainWindow()
    win.show()
    assert win.load(str(fa))
    assert [n for n, _ in win.editor.fasta_records] == ["big", "mid", "small"]
    assert [win.picker.itemText(i) for i in range(3)] == ["big", "mid", "small"] and win.editor.fasta_record == "big"


def test_workspace_tabs_and_view_picker(app, errors):
    """Workspace tabs and panes (ORF pane follows edits, split / close, "+", middle-click copy, view switching, standard tabs for a new window), and the view picker groups with automatic tab names."""
    # --- workspace tabs and panes ---
    from map_view import MapView
    from sequence_view import SequenceView
    from tools import OrfDialog, PcrDialog
    win = open_window(app, PLASMID_DNA)
    win.reset_layout()   # earlier tests may have saved other tabs
    ws, ed = win.workspace, win.editor
    from workspace import pane, split
    # One standard tab, named after the views in it (the first three, then a count).
    assert [ws.tabText(i) for i in range(ws.count())] == ["Sequence · Map · Features"]
    ws.add_tab(None, split("v", [500, 350], split("h", [850, 550], pane("Sequence"), pane("Detect ORFs")), pane("Digest")))
    assert ws.tabText(1) == "Sequence · Detect ORFs · Digest"
    # Panes of tabs never shown have no view yet.
    assert all(p.view is None for p in ws.panes(ws.widget(1)))

    # Tools tab: the ORF pane follows edits (after the typing pause) without reopening.
    ws.setCurrentIndex(1)
    app.processEvents()
    orfs = ws.views(OrfDialog)[0]
    before = len(orfs.orfs)
    ed.set_cursor(0)
    ed.type_bases("ATG" + "GCT" * 150 + "TAA")
    QTest.qWait(600)
    # the new ORF is listed (a plasmid is a ring: an insert at the origin can also join or split ORFs that cross it)
    assert (0, 456, 1) in orfs.orfs and len(orfs.orfs) > before
    ed.undo()

    # Split a pane, turn the new one into a map, close it again.
    seq_pane = next(p for p in ws.panes(ws.currentWidget()) if p.view_name == "Sequence")
    new = ws.split(seq_pane, Qt.Vertical)
    app.processEvents()
    new.set_view("Map")
    assert isinstance(new.view, MapView) and new.view.isVisible()
    ed.select(100, 200)  # every view of the window follows the same selection
    assert new.view.editor is ed and seq_pane.view.editor is ed
    shot(win, "workspace_split")
    ws.close_pane(new)
    app.processEvents()
    assert all(not isinstance(p.view, MapView) for p in ws.panes(ws.currentWidget()))

    # The "+" right after the last tab asks which view or tool to open (the picker's popup) and adds no tab until one is
    # chosen: closing the popup opens nothing; choosing a view opens a tab with that view only.
    from PySide6.QtWidgets import QPushButton
    plus = ws.tabBar().add
    assert plus.isVisible() and plus.x() >= ws.tabBar().tabRect(ws.count() - 1).right()
    count = ws.count()
    plus.click()
    popup = ws.tabBar().picker.popup
    assert popup.isVisible() and ws.count() == count, "the menu first, no tab yet"
    popup.close()
    assert ws.count() == count, "dismissed: nothing opened"
    plus.click()
    next(b for b in ws.tabBar().picker.popup.findChildren(QPushButton, "item") if b.text() == "Map").click()
    app.processEvents()
    assert ws.count() == count + 1 and ws.currentIndex() == count
    assert [p.view_name for p in ws.panes(ws.currentWidget())] == ["Map"], "only the chosen view"
    ws.close_tab(count)

    # Middle-click a tab: a copy of it right after (same panes).
    QTest.mouseClick(ws.tabBar(), Qt.MiddleButton, pos=ws.tabBar().tabRect(0).center())
    assert ws.tabText(1) == ws.tabText(0) and ws.currentIndex() == 1
    assert [p.view_name for p in ws.panes(ws.widget(1))] == [p.view_name for p in ws.panes(ws.widget(0))]
    ws.close_tab(1)

    # Switching a pane's view back and forth: the old views are deleted in between (crashed
    # when the genome view's link to the sequence view was undone).
    p = ws.panes(ws.currentWidget())[0]
    for name in ("Alignment", "Sequence", "Alignment", "Map", "Alignment"):
        p.set_view(name)
        QApplication.sendPostedEvents(None, QEvent.DeferredDelete)
    p.set_view("Sequence")

    # A tab added by hand lives as long as the window; the next window starts from the standard tabs again, so what
    # you get on opening a file never depends on what you left behind (only a .odp project carries its own layout).
    ws.add_tab("Mine", {"split": "h", "sizes": [500, 500],
                        "children": [{"view": "Sequence"}, {"view": "PCR"}]}, select=True)
    app.processEvents()
    assert isinstance(ws.views(PcrDialog)[0], PcrDialog)
    before = [ws.tabText(i) for i in range(ws.count())]
    assert before[-1] == "Mine" and "Sequence · Detect ORFs · Digest" in before
    win.close()
    other = open_window(app, PLASMID_DNA)
    names = [other.workspace.tabText(i) for i in range(other.workspace.count())]
    assert names != before and "Mine" not in names, "a new window starts from the standard tabs"
    assert names == ["Sequence · Map · Features"], "one standard tab, nothing more"
    other.reset_layout()
    assert other.workspace.count() == 1 and isinstance(other.sequence_view, SequenceView)
    other.close()

    # --- view picker groups and automatic tab names ---
    import theme
    import view_picker
    import workspace
    win = open_window(app, PLASMID_DNA)
    ws = win.workspace
    groups = workspace.grouped_views()
    assert [g for g, _ in groups] == ["Views", "Annotations", "Tools"]
    assert groups[0][1] == ["Alignment", "Map", "Protein", "Sequence"]
    assert groups[1][1] == ["Features", "Primers", workspace.SEPARATOR, "Detect Features", "Detect ORFs"]
    tools = groups[2][1]
    assert tools == sorted(tools) and {"Assembly", "Digest", "Local Align", "PCR", "Spacer Finder"} <= set(tools)
    assert "Display Filter" not in sum((n for _, n in groups), []) and "Quick Align" not in tools and "Primer Design" not in tools
    assert sum(len([n for n in names if n]) for _, names in groups) == len(workspace.VIEWS) + len(workspace.TAB_TOOLS), \
        "every view is in exactly one group"

    win.reset_layout()   # earlier tests may have saved other tabs
    ws = win.workspace
    pane = ws.ordered_panes(ws.widget(0))[0]   # the first pane in the layout, Sequence
    assert pane.picker.current() == "Sequence" and pane.picker.text().startswith("Sequence")
    # The picker opens a popup with a column per category: a heading, then its views, indented.
    from PySide6.QtWidgets import QLabel, QPushButton, QWidget
    popup = pane.picker.open_popup()
    columns = popup.findChildren(QWidget, "column")
    assert [c.findChild(QLabel, "heading").text() for c in columns] == ["Views", "Annotations", "Tools"], "one column per category"
    assert [[b.text() for b in c.findChildren(QPushButton, "item")] for c in columns] == [[x for x in n if x] for _, n in groups]
    assert len(columns[1].findChildren(QWidget, "divider")) == 1, "a thin line between Features / Primers and the Detect tools"
    style = popup.styleSheet()
    assert f"padding: 5px 16px 5px {view_picker.INDENT}px" in style and view_picker.INDENT == 16, \
        "the columns already group the views: no indentation, items line up with their heading"
    for role in ("group_header_bg", "group_item_bg", "group_item_hover_bg", "group_border"):
        assert theme.hex(role) in style, "the popup takes its colours from the theme"
    assert any(b.property("current") for b in popup.findChildren(QPushButton, "item")), "the current view is marked"
    # Choosing an item (the Map one) switches the pane and closes the popup.
    next(b for b in popup.findChildren(QPushButton, "item") if b.text() == "Map").click()
    assert pane.view_name == "Map" and pane.picker.current() == "Map" and not popup.isVisible()
    # Assembly is in the list too (under Tools); it is a tab, so choosing it opens that tab and the pane keeps its view.
    assert win.assembly is None
    popup = pane.picker.open_popup()
    next(b for b in popup.findChildren(QPushButton, "item") if b.text() == "Assembly").click()
    assert win.assembly is not None and ws.currentWidget() is win.assembly and ws.indexOf(win.assembly) >= 0
    assert pane.view_name == "Map" and pane.picker.current() == "Map", "the pane still shows what it showed"
    ws.setCurrentIndex(0)

    # The theme holds the colours: a darker band for a heading than for its items, in both themes.
    for dark in (False, True):
        light_hex, dark_hex = zip(*theme.ROLES.values())
        roles = dict(zip(theme.ROLES, dark_hex if dark else light_hex))
        assert QColor(roles["group_header_bg"]).lightness() < QColor(roles["group_item_bg"]).lightness()
    stub = theme.is_dark
    theme.is_dark = lambda: True
    assert theme.hex("group_item_bg") == theme.ROLES["group_item_bg"][1]
    theme.is_dark = lambda: False
    assert theme.hex("group_item_bg") == theme.ROLES["group_item_bg"][0]
    theme.is_dark = stub

    # A tab is named after its views, and follows when they change.
    assert "Map" in ws.tabText(0) and ws.tabToolTip(0).startswith(ws.tabText(0).split(" +")[0])
    pane.set_view("Protein")
    assert ws.tabText(0).startswith("Protein")
    # A name given by hand sticks and is saved; an empty name goes back to automatic.
    ws.set_tab_name(0, "My edit")
    pane.set_view("Map")
    assert ws.tabText(0) == "My edit" and ws.dump()[0]["name"] == "My edit"
    ws.set_tab_name(0, "")
    assert ws.tabText(0) != "My edit" and "name" not in ws.dump()[0]
    # Layouts saved by an older version with the old fixed names are named after their views now.
    old = workspace.Workspace(win, [{"name": "Edit", "layout": workspace.pane("Sequence")},
                                    {"name": "Tab 3", "layout": workspace.pane("Features")},
                                    {"name": "Mine", "layout": workspace.pane("Map")}])
    assert [old.tabText(i) for i in range(3)] == ["Sequence", "Features", "Mine"]
    old.deleteLater()
    win.reset_layout()   # closing saves the tabs: leave the default ones for the tests that follow
    win.close()


def test_primer_pair_selection_and_circular_origin(app, errors):
    """Click a primer, shift-click one on the other strand: the stretch they amplify is selected, and the circular button asks where to break; circular sequence view shows the sequence twice, selections pass the origin, copy / paste / delete treat them as one stretch."""
    # --- shift-click primer pair and the circular button ---
    win = open_window(app, PLASMID_DNA)
    ed = win.editor
    key = {s.label: k for k, s in ed.site_by_key.items()}
    fwd, rev = ed.site_by_key[key["Verif pBAD30_F"]], ed.site_by_key[key["oAC1346"]]
    for first, second in ((key["oAC1346"], key["Verif pBAD30_F"]), (key["Verif pBAD30_F"], key["oAC1346"])):
        ed.select_primer(first)
        ed.select_primer(second, extend=True)
        assert (ed.start, ed.end) == (fwd.spans[0][0], rev.spans[-1][1])
    # Same strand, or no shift: just the clicked primer.
    other_fwd = next(k for k, s in ed.site_by_key.items() if s.strand > 0 and k != key["Verif pBAD30_F"])
    ed.select_primer(key["Verif pBAD30_F"])
    ed.select_primer(other_fwd, extend=True)
    assert ed.current_primer == other_fwd

    # A pair whose product runs through the origin selects across it.
    ed.select_primer(key["oAC736"])                      # 5931..5936 + 1..23
    ed.select_primer(key["Verif pBAD30_R"], extend=True)  # complement(75..94)
    assert (ed.start, ed.end) == (ed.site_by_key[key["oAC736"]].spans[0][0], ed.length + 94) and ed.wraps

    # Circular button: unpressing asks where to break; that base becomes base 1. Cancel
    # keeps it circular; a break inside a feature asks first.
    from PySide6.QtWidgets import QInputDialog
    seq = ed.doc.sequence()
    answers = {"base": (1001, True), "confirm": QMessageBox.Yes}
    orig_get, orig_q = QInputDialog.getInt, QMessageBox.question
    QInputDialog.getInt = staticmethod(lambda *a, **k: answers["base"])
    QMessageBox.question = staticmethod(lambda *a, **k: answers["confirm"])
    try:
        assert ed.doc.circular and win.circular_action.isChecked()
        answers["base"] = (1001, False)
        win.circular_action.trigger()                      # cancelled
        assert ed.doc.circular and win.circular_action.isChecked()
        assert [f.label for f in ed.features_cut_at(1000)] == ["Ampicilin"]  # 521..1381
        answers.update(base=(1001, True), confirm=QMessageBox.No)
        win.circular_action.trigger()                      # cuts the bla CDS: declined
        assert ed.doc.circular and ed.doc.sequence() == seq
        answers["confirm"] = QMessageBox.Yes
        win.circular_action.trigger()
        assert not ed.doc.circular and not win.circular_action.isChecked()
        assert ed.doc.sequence() == seq[1000:] + seq[:1000] and ed.doc.modified
        amp = next(f for f in ed.features if f.label == "Ampicilin")
        assert len(amp.spans) == 2                         # cut in two by the break
        assert win.history_dock.widget() is not None
        ed.undo()                                          # one step: the circle and its old origin come back
        assert ed.doc.sequence() == seq and ed.doc.circular and win.circular_action.isChecked()
        ed.redo()
        assert not ed.doc.circular and not win.circular_action.isChecked()
        win.circular_action.trigger()                      # pressed again: circular, an undo step too
        assert ed.doc.circular
        ed.undo()
        assert not ed.doc.circular and not win.circular_action.isChecked()
    finally:
        QInputDialog.getInt, QMessageBox.question = orig_get, orig_q
    ed.doc.modified = False
    win.close()

    # --- circular selection through the origin ---
    win = open_window(app, PLASMID_DNA)
    ed, sv = win.editor, win.sequence_view
    n, seq = ed.length, ed.doc.sequence()
    # Endless scrolling is off until the red origin bar (end of the last row) is clicked.
    assert sv.copies == 1 and not sv.loop
    bpl, last = sv.bases_per_line(), sv.rpc - 1
    sv.verticalScrollBar().setValue(sv.offset(last))
    y = sv.offset(last) - sv.verticalScrollBar().value() + sv._parts(last)[1] + 5
    x = sv.left + (n - last * bpl) * sv.cw + 1
    assert sv.origin_at(x, y) and not sv.origin_at(x - 3 * sv.cw, y)
    QTest.mouseClick(sv.viewport(), Qt.LeftButton, pos=QPoint(int(x), int(y)))
    assert sv.loop and sv.copies >= 3 and sv.rows == sv.copies * sv.rpc

    # Endless scrolling: past the end comes base 1 again, above base 1 the end; the scroll
    # position stays in the middle copy (shifting by a whole copy changes nothing on screen).
    sb, one = sv.verticalScrollBar(), sv.offset(sv.rpc)
    assert one <= sb.value() < 2 * one and sb.maximum() >= 2 * one
    sb.setValue(2 * one + 30)
    assert sb.value() == one + 30
    sb.setValue(one - 30)
    assert sb.value() == 2 * one - 30
    for _ in range(3 * sv.rpc):  # scroll down three whole turns, one line at a time
        sb.triggerAction(QAbstractSlider.SliderSingleStepAdd)
    assert one <= sb.value() < 2 * one

    # Keyboard: from 2 bases before the end, Shift+Right x5 runs through the origin.
    ed.set_cursor(n - 2)
    sv.anchor = n - 2
    sv.setFocus()
    for _ in range(5):
        QTest.keyClick(sv, Qt.Key_Right, Qt.ShiftModifier)
    assert (ed.start, ed.end) == (n - 2, n + 3) and ed.selected_sequence() == seq[-2:] + seq[:3]
    assert "across the origin" in win.status.text()

    # A screen position in the second copy is that base again; at most one full copy.
    pos, row, _ = sv.pos_at(sv.left + 3 * sv.cw, sv.offset(sv.rpc) - sv.verticalScrollBar().value() + 5)
    assert row == sv.rpc and pos == n + 3
    ed.select(10, 10 + 3 * n)
    assert (ed.start, ed.end) == (10, 10 + n)
    ed.set_cursor(n)
    assert ed.start == n  # a cursor at the very end stays there (typing appends)

    # A feature across the origin: selected whole; copied as ONE feature; pasted elsewhere.
    ed.add_features([("misc_feature", f"join({n - 9}..{n},1..10)", [("label", "across")])])
    assert (ed.start, ed.end) == (n - 10, n + 10) and ed.wraps
    sv.grab()
    QTest.keyClick(sv, Qt.Key_C, Qt.ControlModifier)
    ed.set_cursor(1000)
    QTest.keyClick(sv, Qt.Key_V, Qt.ControlModifier)
    pasted = [f for f in ed.features if f.label == "across" and f.spans[0][0] == 1000]
    assert len(pasted) == 1 and sum(e - s for s, e in pasted[0].spans) == 20
    assert ed.doc.sequence(1000, 1020) == seq[-10:] + seq[:10]
    ed.undo()

    # Delete across the origin: one undo step brings it all back.
    ed.select_feature(next(f.index for f in ed.features if f.label == "across"))
    ed.backspace()
    assert ed.length == n - 20 and ed.doc.sequence() == seq[10:-10]
    ed.undo()
    assert ed.doc.sequence() == seq
    ed.undo()  # the added feature
    assert not ed.doc.modified

    # Linear: one copy, selections clamp at the ends.
    ed.set_circular(False)
    assert sv.copies == 1 and sv.rows == sv.rpc
    ed.select(n - 5, n + 5)
    assert (ed.start, ed.end) == (n - 5, n)
    ed.set_circular(True)
    ed.doc.modified = False
    win.close()


def test_display_filters_translation_menus_history_and_new_file(app, errors, monkeypatch):
    """Display filter, quick filters, standard colours, select-from-to, CDS translation, letter size, selection brackets, Find by name, ORFs and the shortcuts dialog; the history panel (changes in words, go back / forward, drag the marker); File > New and the selection menu."""
    # --- display filter, translation menus, zoom and shortcuts ---
    from PySide6.QtWidgets import QMenu
    from feature_table import feature_actions
    from shortcuts import ShortcutsDialog
    win = open_window(app, PLASMID_DNA)
    ed, sv, panel = win.editor, win.sequence_view, win.display_dock.widget()
    n_cds = sum(f.kind == "CDS" for f in ed.features)
    # The side panel starts folded away; the Filters button opens and closes it.
    assert not win.display_dock.isVisible()
    win.display_dock.toggleViewAction().trigger()
    assert win.display_dock.isVisible()

    # Display filter: a type switched off is not drawn (but stays in the table), primers off.
    panel.kinds["CDS"].setChecked(False)
    assert "CDS" in ed.hidden_kinds and not any(f.kind == "CDS" for f in ed.shown)
    assert sum(f.kind == "CDS" for f in ed.features) == n_cds and all(f.kind != "CDS" for f in ed.visible_features(0, ed.length))
    win.map_view.grab()
    panel.kinds["CDS"].setChecked(True)
    assert not ed.hidden_kinds
    panel.file_primers.setChecked(False)
    panel.db_primers.setChecked(False)
    assert ed.drawn_primers == [] and ed.primer_sites and not sv.row_fwd and not sv.row_rev
    panel.file_primers.setChecked(True)
    panel.db_primers.setChecked(True)
    assert not ed.show_enzymes and not win.enzymes_action.isChecked()  # off by default
    win.enzymes_action.trigger()  # the menu and the panel stay in step
    assert ed.show_enzymes and panel.enzymes.isChecked()
    win.enzymes_action.trigger()

    # Quick filters at the right of the tab row: the same switches, each way; as wide as their buttons (the tabs keep the rest).
    quick = win.workspace.quick_filters
    assert quick.width() <= quick.sizeHint().width() + 4 < win.workspace.width() // 4
    quick.buttons["enzymes"].click()
    assert ed.show_enzymes and panel.enzymes.isChecked() and win.enzymes_action.isChecked()
    win.enzymes_action.trigger()
    assert not quick.buttons["enzymes"].isChecked()
    quick.buttons["primers"].click()
    assert not ed.show_primers and not panel.file_primers.isChecked()
    quick.buttons["primers"].click()
    assert ed.show_file_primers and ed.show_db_primers
    panel.kinds["CDS"].setChecked(False)        # Annotations off then on gives back the types hidden one by one
    quick.buttons["annotations"].click()
    assert all(not b.isChecked() for b in panel.kinds.values())
    quick.buttons["annotations"].click()
    assert ed.hidden_kinds == {"CDS"} and quick.buttons["annotations"].isChecked()
    panel.kinds["CDS"].setChecked(True)
    if SHOTS:
        win.grab().save(os.path.join(SHOTS, "quick_filters.png"))

    # Standard colours: a new feature's colour follows its type; recolour all in one step.
    from editor import type_color
    from feature_dialog import FeatureDialog
    ed.select(100, 200)
    dialog = FeatureDialog(win, ed)
    dialog.kind.setCurrentText("terminator")
    assert dialog.color == type_color("terminator")
    n_changed = ed.recolor_by_type()
    assert n_changed and all(ed.colors[f.index] == type_color(f.kind) for f in ed.features)
    ed.undo()

    # The range box: 101..200 (1-based, inclusive); an end before the start goes through the origin.
    win.range_box.setText("101..200")
    win.range_box.returnPressed.emit()
    assert (ed.start, ed.end) == (100, 200)
    win.select_range(f"{ed.length - 9:,} - 10")
    assert (ed.start, ed.end) == (ed.length - 10, ed.length + 10)

    # CDS translation: eilR's codons, from its start codon; switched off per CDS.
    fos = next(f for f in ed.features if f.label == "eilR")
    codons = sv._codons(fos)
    assert codons[0][:2] == (fos.spans[0][0] + 1, "M") and len(codons) == (fos.spans[0][1] - fos.spans[0][0]) // 3
    assert codons[0][2:] == (fos.spans[0][0], fos.spans[0][0] + 3)  # its bases (for the ticks)
    assert ed.translated(fos)
    from hover import describe
    assert "eilR" in describe(ed, fos.index) and "aa" in describe(ed, fos.index)
    menu = QMenu()
    feature_actions(menu, win, ed, fos.index)
    texts = [a.text() for a in menu.actions()]
    assert {"Edit…", "Delete", "Show Translation"} <= set(texts)
    next(a for a in menu.actions() if a.text() == "Show Translation").trigger()
    assert not ed.translated(fos)
    ed.select_feature(fos.index)
    sv.grab()

    # Letter size: Ctrl+wheel / Ctrl+= / Ctrl+0; the same bases stay at the top.
    ed.set_cursor(2000)
    top = sv.real_row(sv.row_at(sv.verticalScrollBar().value())) * sv.bases_per_line()
    cw = sv.cw0
    sv.setFocus()
    QTest.keyClick(sv, Qt.Key_Equal, Qt.ControlModifier)
    QTest.keyClick(sv, Qt.Key_Equal, Qt.ControlModifier)
    assert sv.size == 13 and sv.cw0 > cw
    top2 = sv.real_row(sv.row_at(sv.verticalScrollBar().value())) * sv.bases_per_line()
    assert abs(top2 - top) < 2 * sv.bases_per_line()
    QTest.keyClick(sv, Qt.Key_0, Qt.ControlModifier)
    assert sv.size == 11

    # Selection brackets: drag the "]" two bases right to extend the selection.
    ed.select(100, 110)
    sv.ensure_visible(100)
    app.processEvents()
    bpl, sb = sv.bases_per_line(), sv.verticalScrollBar()
    row = next(r for r in range(sv.rows) if sv.real_row(r) == 110 // bpl and sv.offset(r) >= sb.value())
    y = sv.offset(row) - sb.value() + sv._parts(110 // bpl)[1] + 5
    x = sv.left + (110 % bpl) * sv.cw
    assert sv.bracket_at(x, y) is not None and sv.bracket_at(x + sv.cw / 2, y) is None
    QTest.mousePress(sv.viewport(), Qt.LeftButton, pos=QPoint(int(x), int(y)))
    QTest.mouseMove(sv.viewport(), QPoint(int(x + 2 * sv.cw), int(y)))
    QTest.mouseRelease(sv.viewport(), Qt.LeftButton, pos=QPoint(int(x + 2 * sv.cw), int(y)))
    assert (ed.start, ed.end) == (100, 112)

    # Find: a gene whose name is made of DNA letters (cat, the chloramphenicol gene) is looked up as a name, not only as
    # bases. "Go to" prefers the feature; plain Find still reports the places where those bases occur.
    ed.add_features([("CDS", "200..260", [("label", "cat")])])
    named = next(f for f in ed.features if f.label == "cat")
    assert ed.lookup("cat") == (199, 260, named.index), "an exact name wins over the DNA that spells it"
    assert ed.lookup("CAT")[2] == named.index, "any case"
    seq = ed.doc.sequence().upper()
    hit = ed.find("cat", 0)
    assert hit is not None and (hit[2] == named.index or seq[hit[0]:hit[1]] in ("CAT", "ATG")), "Find also looks for the bases"
    ed.undo()

    # ORFs: an ORF already annotated as a CDS links to it; protein motif search.
    from tools import OrfDialog
    orf = OrfDialog(win, ed)
    row = next(i for i, f in enumerate(orf.links) if f is not None and f.label == "eilR")
    assert "eilR (annotated)" in orf.table.item(row, 3).text()
    orf._select(row)
    assert ed.current_feature == fos.index
    orf.motif.setText("".join(aa for _, aa, _, _ in sv._codons(fos)[1:5]))  # 4 aa of eilR
    assert orf.found and orf.hits.item(0).text().startswith("eilR:")
    orf.hits.setCurrentRow(0)
    assert (ed.start, ed.end) == (fos.spans[0][0] + 3, fos.spans[0][0] + 15)

    # Help → Keyboard Shortcuts: menu shortcuts read from the menus, plus a tab per panel.
    dialog = ShortcutsDialog(win)
    menus = dialog.tabs.widget(0)
    keys = {menus.item(i, 0).text() for i in range(menus.rowCount())}
    assert "Ctrl+Shift+F" in keys and dialog.tabs.count() >= 5
    ed.doc.modified = False
    win.close()

    # --- history panel ---
    win = open_window(app, PLASMID_DNA)
    ed = win.editor
    seq = ed.doc.sequence()
    win.history_dock.toggleViewAction().trigger()
    assert win.history_dock.isVisible() and not win.display_dock.isVisible()  # one or the other
    panel = win.history_dock.widget()

    ed.set_cursor(100)
    for c in "GAG":
        ed.type_bases(c)                       # 3 undo steps, one row
    ed.select(200, 203)
    ed.type_bases("NNN")                       # replace
    ed.set_cursor(300)
    ed.backspace()
    ed.backspace()                             # 2 steps, one row
    ed.add_features([("misc_feature", "10..20", [("label", "mine")])])
    from history_panel import NOW
    texts = [panel.list.item(i).text() for i in range(panel.list.count())]
    assert texts[0] == f"Opened {PLASMID_DNA.name}"
    assert texts[1] == "Typed 3 bp at 101: GAG"
    assert texts[2] == f"{seq[197:200]} > NNN at 201"  # (GAG was inserted before it)
    assert texts[3].startswith("Deleted 2 bp at 299")
    assert texts[4] == "Added feature mine" and texts[5] == NOW

    panel.list.setCurrentRow(2)                # back to before the replacement
    panel.go()
    assert ed.doc.sequence() == seq[:100] + "GAG" + seq[100:]
    assert panel.list.item(2).text() == NOW and panel.button.text() == "Go back to before this change"
    panel.list.setCurrentRow(panel.list.count() - 1)  # the last undone change: redo it all
    assert panel.button.text() == "Redo up to this change"
    panel.go()
    assert len(ed.features) and any(f.label == "mine" for f in ed.features)

    # Dragging the current-sequence marker: up above "Typed…" = back to before it; then
    # down to the end = everything redone. (Qt's drop reports the old row it lands before.)
    panel._moved(None, panel.now_row, panel.now_row, None, 1)
    QTest.qWait(10)
    assert ed.doc.sequence() == seq and panel.list.item(1).text() == NOW
    panel._moved(None, panel.now_row, panel.now_row, None, panel.list.count())
    QTest.qWait(10)
    assert any(f.label == "mine" for f in ed.features) and panel.list.item(panel.list.count() - 1).text() == NOW

    panel.list.setCurrentRow(0)                # "Opened …": back to the file as opened
    assert panel.button.text() == "Back to the file as opened"
    panel.go()
    assert ed.doc.sequence() == seq
    win.display_dock.toggleViewAction().trigger()
    assert win.display_dock.isVisible() and not win.history_dock.isVisible()
    ed.doc.modified = False
    win.close()

    # --- new file and selection menu ---
    from PySide6.QtWidgets import QInputDialog, QPlainTextEdit
    from feature_table import features_inside, group_as_gene, reverse_primer_dialog

    # File → New: into this window when it is empty; pasted text is cleaned up.
    win = MainWindow()
    win.show()
    filled = {"text": "1 acgtac gtNN\n7 ACGT"}

    def fill(dialog):
        dialog.findChild(QPlainTextEdit).setPlainText(filled["text"])
        return QDialog.Accepted

    with monkeypatch.context() as m:
        m.setattr(QDialog, "exec", fill)
        assert win.new_file() is win
        ed = win.editor
        assert ed.doc.sequence() == "acgtacgtNNACGT" and ed.circular and ed.untitled and ed.doc.modified
        filled["text"] = ""
        other = win.new_file()                     # this window has a document: a new window
        assert other is not win and other.editor.length == 0
        other.editor.type_bases("ACGT")            # an empty sequence can be typed into
        assert other.editor.doc.sequence() == "ACGT"
        for w in (win, other):
            w.editor.doc.modified = False
            w.close()

    # Right-click menu on a selection: reverse primer prefilled; gene named after the largest CDS.
    win = open_window(app, PLASMID_DNA)
    ed = win.editor
    ed.select(944, 1004)
    assert reverse_primer_dialog(win, ed).sequence.text() == reverse_complement(ed.selected_sequence())
    eil = next(f for f in ed.features if f.label == "eilR")
    ed.select(eil.spans[0][0] - 20, eil.spans[-1][1] + 20)
    inside = features_inside(ed)
    assert eil.index in inside
    asked = {}

    def answer(parent, title, prompt, text=""):
        asked.update(prompt=prompt, text=text)
        return "", False                       # cancel: only the question is checked

    monkeypatch.setattr(QInputDialog, "getText", answer)
    group_as_gene(win, ed, inside)
    assert asked["text"] == "eilR" and "no CDS" not in asked["prompt"]
    group_as_gene(win, ed, [f.index for f in ed.features if f.kind == "RBS"])
    assert "no CDS" in asked["prompt"]

    # DNA-only and protein-only fields ignore other keys.
    from primer_panel import PrimerDialog
    from tools import OrfDialog
    dialog = PrimerDialog(win, ed)
    dialog.sequence.clear()
    QTest.keyClicks(dialog.sequence, "gaATCG1x-Z?T")
    assert dialog.sequence.text() == "gaATCGT"
    orf = OrfDialog(win, ed)
    QTest.keyClicks(orf.motif, "MK1L*j")
    assert orf.motif.text() == "MKL*"
    win.close()


def test_start_page_sponsors_and_updates(app, errors, monkeypatch, tmp_path):
    """Startup page: welcome message, shown once a day and never again once turned off; a click on a recent file opens it; the update check answers in place; pinned / removed files, news after an update, tips, quick actions; sponsor cards and contact address only when set."""
    # --- sponsors page ---
    import sponsors
    dialog = sponsors.SponsorsDialog()
    assert dialog.cards == [], "no sponsors: no empty 'Your ad here' boxes"
    welcome = " ".join(w.text() for w in dialog.findChildren(QLabel))
    assert "Antoine Champie" in welcome and "official download" in welcome and "github.com/Ant-Champ/GenomeEdit" in welcome
    assert sponsors.due()
    page = sponsors.show_at_startup(None)
    page.accept()
    assert not sponsors.due(), "already shown today"
    assert sponsors.show_at_startup(None) is None
    sponsors.settings().setValue("sponsors/last_shown", "2000-01-01")
    assert sponsors.due(), "a new day: shown again"
    dialog = sponsors.SponsorsDialog()
    dialog.snooze.setChecked(True)
    dialog.accept()
    assert not sponsors.due(), "turned off"
    dialog = sponsors.SponsorsDialog()
    assert dialog.snooze.isChecked()
    dialog.snooze.setChecked(False)
    dialog.accept()
    assert sponsors.due(), "turned back on from Help > Welcome"
    sponsors.settings().setValue("sponsors/never", True)
    assert sponsors.show_at_startup(None) is None, "turned off"
    forced = sponsors.show_at_startup(None, force=True)   # the --welcome option
    assert forced is not None and forced.snooze.isChecked()
    forced.accept()
    assert sponsors.hidden(), "opening it by force does not turn it back on"

    # --- recent files and the update check ---
    with monkeypatch.context() as m:
        from types import SimpleNamespace
        from PySide6.QtWidgets import QDialog, QWidget
        import project
        import sponsors
        import updater
        gb = tmp_path / "pX.gb"
        gb.write_text("x")
        m.setattr(sponsors, "recent_files", lambda: [str(gb)])
        window, opened = QWidget(), []
        window.open_path = opened.append
        dialog = sponsors.SponsorsDialog(window)
        assert dialog.recent.count() == 1 and dialog.recent.item(0).text().startswith("pX.gb")
        dialog.recent.itemClicked.emit(dialog.recent.item(0))
        assert opened == [str(gb)] and dialog.result() == QDialog.Accepted

        m.setattr(sponsors, "recent_files", lambda: [])
        dialog = sponsors.SponsorsDialog()
        assert not dialog.recent.isEnabled(), "nothing to open: the list says so and takes no click"
        m.setattr(project, "UPDATES_ENABLED", False)
        dialog.check_updates()
        assert "Not available yet" in dialog.update_status.text()
        m.setattr(project, "UPDATES_ENABLED", True)
        m.setattr(updater, "mark_checked", lambda *a: None)
        answers = iter([SimpleNamespace(update=None, error="", page_url=""),
                        SimpleNamespace(update=SimpleNamespace(version="9.9.9"), error="", page_url="")])
        m.setattr(updater, "check_in_background", lambda callback, **_: callback(next(answers)))
        dialog.check_updates()
        assert "Up to date" in dialog.update_status.text()
        dialog.check_updates()
        assert "9.9.9" in dialog.update_status.text() and "install:" in dialog.update_status.text()

    # --- pins, removal, news, tips and quick actions ---
    with monkeypatch.context() as m:
        from PySide6.QtWidgets import QWidget
        import project
        import sourcefiles
        import sponsors
        a, b = tmp_path / "a.gb", tmp_path / "b.gb"
        for f in (a, b):
            f.write_text("x")
        sourcefiles.add_recent(a)
        sourcefiles.add_recent(b)    # newest first: b, a
        dialog = sponsors.SponsorsDialog()
        names = lambda: [dialog.recent.item(i).data(Qt.UserRole) for i in range(dialog.recent.count())]
        assert names()[:2] == [str(b.resolve()), str(a.resolve())]
        dialog.pin(str(a.resolve()), True)
        assert names()[0] == str(a.resolve()) and dialog.recent.item(0).font().bold()
        dialog.remove(str(a.resolve()))
        assert str(a.resolve()) not in names() and a.exists() and str(a.resolve()) not in sourcefiles.pinned_files()
        sourcefiles.remove_recent(str(b.resolve()))

        sponsors.settings().remove("sponsors/seen_version")
        assert sponsors.whats_new() is None, "a first install has nothing 'new'"
        sponsors.settings().setValue("sponsors/seen_version", "0.0.1")
        m.setitem(project.WHATS_NEW, project.VERSION, ["Pinned files"])
        dialog = sponsors.SponsorsDialog()
        assert not dialog.news.isHidden() and "Pinned files" in dialog.news.text()
        assert sponsors.SponsorsDialog().news.isHidden(), "shown once"

        first = dialog.tip.text()
        dialog._link("app:next-tip")
        assert dialog.tip.text() != first and len(set(sponsors.TIPS)) == len(sponsors.TIPS) > 1

        assert not any(b.isEnabled() for b in dialog.quick.values()), "no main window behind: nothing to act on"
        assert sponsors.DEMO is not None and sponsors.DEMO.name == sponsors.DEMO_NAME
        window, opened = QWidget(), []
        window.open_untitled = lambda doc, path: opened.append((doc.name, path))
        dialog = sponsors.SponsorsDialog(window)
        assert dialog.quick["demo"].isEnabled() and not dialog.quick["new_file"].isEnabled()
        dialog.quick["demo"].click()
        assert opened and opened[0][1].name == sponsors.DEMO_NAME and dialog.result() == QDialog.Accepted
        sponsors.settings().remove("sponsors")

    # --- sponsor cards and the contact address ---
    with monkeypatch.context() as m:
        import project
        import sponsors
        m.setattr(sponsors, "load", lambda: ([{"name": "Acme", "text": "Reagents", "url": "https://example.org"}, {}, {}], "ads@example.org"))
        dialog = sponsors.SponsorsDialog()
        assert len(dialog.cards) == 3 and "Acme" in dialog.cards[0].findChild(QLabel).text()
        assert any(isinstance(w, QLabel) and w.text() == "SPONSORED" and not w.isHidden() for w in dialog.findChildren(QLabel))
        m.setattr(project, "CONTACT_EMAIL", "")
        assert "mailto" not in sponsors.welcome_html(), "no address set: none is shown"
        m.setattr(project, "CONTACT_EMAIL", "genomeedit@example.org")
        assert "mailto:genomeedit@example.org" in sponsors.welcome_html()


def test_gene_groups_and_detect_common_features(app, errors, tmp_path):
    """Gene groups (make, table, map, move a part, copy-paste, save / reopen, ungroup / delete) and Detect Common Features on a synthetic plasmid of library parts (one recoded), grouped adding, settings, your own parts."""
    # --- gene groups ---
    import genes
    win = open_window(app, PLASMID_DNA)
    ed, table = win.editor, win.feature_table
    n = len(ed.features)
    cds = max((f for f in ed.features if f.kind == "CDS"), key=lambda f: f.spans[0][1] - f.spans[0][0])
    s, e = cds.spans[0]
    # A promoter and an RBS just upstream (on the CDS strand), as typed by a user.
    up = (lambda a, b: f"{a + 1}..{b}") if cds.strand > 0 else (lambda a, b: f"complement({a + 1}..{b})")
    pro_at, rbs_at = ((s - 60, s - 20), (s - 15, s - 8)) if cds.strand > 0 else ((e + 20, e + 60), (e + 8, e + 15))
    ed.add_features([("promoter", up(*pro_at), [("label", "pTest")]), ("RBS", up(*rbs_at), [("label", "rbsTest")])])
    pro, rbs = (f for f in ed.features if f.label in ("pTest", "rbsTest"))
    key = genes.make_gene(ed, [pro.index, rbs.index, cds.index], "myGene")
    g = ed.groups[key]
    assert [p.kind for p in g.parts] == ["promoter", "RBS", "CDS"]
    assert genes.span_of([g.master]) == (min(pro_at[0], s), max(pro_at[1], e)) and g.master.strand == cds.strand
    assert len(ed.features) == n + 3
    ed.undo()  # grouping is one step
    assert not ed.groups and len(ed.features) == n + 2
    ed.redo()
    g = ed.groups[key]

    # Table: one bold gene row, parts hidden until opened; map: the gene, not its parts.
    rows = [r[3] for r in table.model.shown]
    assert g.master.index in rows and all(p.index not in rows for p in g.parts)
    table.model.toggle(table.model.row_of(g.master.index))
    rows = [r[3] for r in table.model.shown]
    i = rows.index(g.master.index)
    assert rows[i + 1:i + 4] == [p.index for p in g.parts]
    assert g.master.index in ed.map_lanes and all(p.index not in ed.map_lanes for p in g.parts)
    win.collapse_action.trigger()  # off: parts on the map too
    assert all(p.index in ed.map_lanes for p in g.parts)
    win.collapse_action.trigger()
    shot(win, "gene_group")

    # Moving a part outside the gene grows the gene, in the same undo step.
    pro = next(p for p in g.parts if p.kind == "promoter")
    far = (pro.spans[0][0] - 100, pro.spans[0][1]) if cds.strand > 0 else (pro.spans[0][0], pro.spans[0][1] + 100)
    ed.set_feature(pro.index, "promoter", up(*far), pro.qualifiers)
    assert genes.span_of([ed.groups[key].master])[0 if cds.strand > 0 else 1] == far[0 if cds.strand > 0 else 1]
    ed.undo()
    assert genes.span_of([ed.groups[key].master]) == genes.span_of([g.master])

    # Copy-paste the whole gene: the copy is a second, separate gene.
    ed.select_feature(ed.groups[key].master.index)
    win.sequence_view.setFocus()
    QTest.keyClick(win.sequence_view, Qt.Key_C, Qt.ControlModifier)
    ed.set_cursor(0)
    QTest.keyClick(win.sequence_view, Qt.Key_V, Qt.ControlModifier)
    assert len(ed.groups) == 2 and all(len(x.parts) == 3 for x in ed.groups.values())
    ed.undo()

    # Saved as plain GenBank (shared /locus_tag); groups come back on reopening.
    out = tmp_path / "grouped.gb"
    ed.save(out)
    assert out.read_text().count(f'/locus_tag="{key}"') == 4
    assert win.load(str(out))
    ed = win.editor
    assert [p.kind for p in ed.groups[key].parts] == ["promoter", "RBS", "CDS"]

    # Ungroup keeps the parts; delete takes them.
    genes.ungroup(ed, key)
    assert not ed.groups and len(ed.features) == n + 2
    ed.undo()
    genes.delete_gene(ed, key)
    assert not ed.groups and len(ed.features) == n - 1
    ed.doc.modified = False
    win.close()

    # --- detect common features ---
    import random
    import annotation
    import tool_settings
    from annotation import DetectDialog, load_parts
    parts = {p.name: p for p in load_parts()}
    # Canonical names: the library keeps one row per part (scripts/build_parts.py, docs/PARTS_LIBRARY_METHOD.md).
    promoter, gfp, term = parts["lac promoter"], parts["sfGFP"], parts["T7 terminator"]
    amp = parts["AmpR (bla, TEM-1)"]
    rng = random.Random(11)
    rand = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
    # sfGFP recoded: the third base of every Gly/Ala/Val/Leu codon (GGx/GCx/GTx/CTx) changed,
    # which keeps the protein but changes ~10 % of the bases (below the 95 % DNA threshold).
    codons = [gfp.sequence[i:i + 3] for i in range(0, len(gfp.sequence), 3)]
    recoded = "".join(c[:2] + ("C" if c[2] != "C" else "T") if c[:2] in ("GG", "GC", "GT", "CT") else c
                      for c in codons)
    assert sum(a != b for a, b in zip(recoded, gfp.sequence)) > 0.05 * len(recoded)
    gene_part = promoter.sequence + rand(20) + recoded + rand(30) + term.sequence
    seq = rand(400) + gene_part + rand(300) + reverse_complement(amp.sequence) + rand(400)
    gb = tmp_path / "synthetic.gb"
    fa = tmp_path / "synthetic.fa"
    fa.write_text(">synthetic\n" + seq + "\n")
    win = MainWindow()
    win.resize(1400, 850)
    win.show()
    assert win.load(str(fa))
    ed = win.editor
    dialog = DetectDialog(win, ed)
    found = {(h.name, h.strand, h.protein): h for h in dialog.hits}
    for want in (("AmpR (bla, TEM-1)", -1, False), ("lac promoter", 1, False), ("T7 terminator", 1, False)):
        assert want in found, want
    gfp = found[("sfGFP", 1, True)]  # found through its protein: the DNA differs too much
    assert gfp.start == 400 + len(promoter.sequence) + 20 and gfp.identity == 1.0
    assert all(not h.known for h in dialog.hits)

    # Add the ticked ones, grouped: Lac promoter + sfGFP + T7 terminator become one gene.
    n = len(ed.features)
    dialog.add()
    assert len(ed.features) > n
    g = next(g for g in ed.groups.values() if g.master.label == "sfGFP")
    assert sorted(p.kind for p in g.parts) == ["CDS", "promoter", "terminator"]
    assert all(h.known for h in dialog.hits)  # refreshed: everything is annotated now
    ed.undo()
    assert len(ed.features) == n  # one undo step

    # Settings: with grouping off, no gene is made.
    tool_settings.set_value(annotation.TOOL, "group_genes", False)
    tool_settings.notifier.changed.emit(annotation.TOOL)
    dialog.add()
    assert not ed.groups
    tool_settings.set_value(annotation.TOOL, "group_genes", True)
    ed.undo()
    settings = tool_settings.SettingsDialog(win, annotation.TOOL)
    assert settings.tabs.tabText(settings.tabs.currentIndex()) == tool_settings.LABELS.get(annotation.TOOL, annotation.TOOL)
    p, w = settings.widgets[(annotation.TOOL, "min_length")]
    w.setValue(31)  # the 30-bp Lac promoter is now too short to look for
    settings.accept()  # saves and tells open tools: the dialog re-runs
    assert tool_settings.value(annotation.TOOL, "min_length") == 31
    assert not any(h.name == "Lac promoter" for h in dialog.hits)
    tool_settings.set_value(annotation.TOOL, "min_length", 15)

    # Your own parts: import a file's features, and they are found next time.
    ed.select(100, 160)
    ed.add_features([("misc_feature", "101..160", [("label", "myTag")])])
    ed.save(gb)
    assert annotation.import_files([gb]) >= 1
    ed.delete_feature(next(f.index for f in ed.features if f.label == "myTag"))
    assert any(h.name == "myTag" and not h.known for h in DetectDialog(win, ed).hits)
    ed.doc.modified = False
    win.close()


def mapping_tools_present():
    import subprocess
    import runner
    try:
        out = subprocess.run(runner.command_line(runner.script(runner.CHECK)), capture_output=True, timeout=60)
    except (OSError, subprocess.TimeoutExpired):
        return False
    return out.returncode == 0 and b"MISSING" not in out.stdout


def test_read_colours_local_align_and_map_reads(app, errors, tmp_path, monkeypatch):
    """IGV's pair rules on hand-made reads; Local Align (hits, the alignment view track, old layout names); then bwa mem on simulated reads (skipped without the mapping tools): reads, variants, the variants track, coverage, cancel."""
    # --- read colour categories ---
    from alignments import classify
    pair = lambda start, strand, mate, mate_rev, insert, flags=0x1 | 0x2, elsewhere=False: (
        start, start + 100, strand, 60, "r", [], [], [], flags | (0x20 if mate_rev else 0), mate, elsewhere, insert)
    mode, limits = "Pair orientation & insert size", (50, 1000)
    assert classify(pair(100, 1, 300, True, 300), limits, mode) == "normal"      # → ←
    assert classify(pair(300, -1, 100, False, -300), limits, mode) == "normal"   # its mate's view
    assert classify(pair(100, 1, 4000, True, 4000), limits, mode) == "larger"
    assert classify(pair(100, 1, 110, True, 40), limits, mode) == "smaller"
    assert classify(pair(100, 1, 300, False, 300), limits, mode) == "LL"         # → →
    assert classify(pair(100, -1, 300, True, 300), limits, mode) == "RR"         # ← ←
    assert classify(pair(100, -1, 300, False, 300), limits, mode) == "RL"        # ← →
    assert classify(pair(100, 1, 300, True, 0, elsewhere=True), limits, mode) == "elsewhere"
    assert classify(pair(100, 1, None, True, 0, flags=0x1 | 0x8), limits, mode) == "normal"  # mate unmapped
    assert classify(pair(100, -1, 300, True, 300), limits, "Strand") == "reverse"

    # --- local align ---
    win = open_window(app, PLASMID_DNA)
    ed = win.editor
    seq = ed.doc.sequence()
    # 300 bp of the plasmid with a mismatch, an insertion and a deletion, reverse-complemented.
    piece = list(seq[1000:1300].upper())
    del piece[250:255]
    piece[150:150] = "GGG"
    piece[50] = "A" if piece[50] != "A" else "C"
    query = reverse_complement("".join(piece))
    QApplication.clipboard().setText(">query\n" + query)
    dlg = LocalAlignDialog(win, ed)
    assert dlg.query.toPlainText().replace("\n", "").endswith(query)  # prefilled from the clipboard
    dlg.run()
    app.processEvents()
    top = dlg.hits[0]
    assert (top.start, top.end, top.strand) == (1000, 1300, -1)
    assert len(top.mismatches) == 1 and sum(n for _, n in top.insertions) == 3
    assert sum(e - s for s, e in top.deletions) == 5
    assert dlg.table.item(0, 2).text() == "1,001–1,300" and dlg.table.item(0, 3).text() == INVERTED_TEXT
    # Clicking the row selects the hit.
    ed.set_cursor(0)
    QTest.mouseClick(dlg.table.viewport(), Qt.LeftButton, pos=dlg.table.visualItemRect(dlg.table.item(0, 2)).center())
    assert (ed.start, ed.end) == (1000, 1300)
    # The alignment view draws it as a track.
    canvas = win.genome_view.canvas
    win.grab()
    assert len(ed.query_hits) == len(dlg.hits) and "Query hits" in canvas.kinds
    assert any(a[0] == "qhit" and (ed.query_hits[a[1]].start, ed.query_hits[a[1]].end) == (1000, 1300) for _, a in canvas.hits)
    shot(win, "local_align")
    dlg.clear()
    win.grab()
    assert not ed.query_hits and "Query hits" not in canvas.kinds
    # Saved layouts from before the rename still open.
    assert Pane(win.workspace, "Genome").view_name == "Alignment"
    win.close()

    # --- map reads and the BAM track ---
    import random
    import alignments
    import runner
    import tool_settings
    from variants import annotate_variants, apply_variants
    if not mapping_tools_present():
        pytest.skip("bwa/samtools/minimap2 not installed (Tools → Install Mapping Tools)")
    rng = random.Random(3)
    ref = "".join(rng.choice("ACGT") for _ in range(50_000))
    snp, dele = 20_000, 30_000
    alt = "A" if ref[snp] != "A" else "C"
    sample = ref[:snp] + alt + ref[snp + 1:dele] + ref[dele + 5:]
    fa = tmp_path / "ref one.fa"  # a space: paths must be quoted for the shell
    fa.write_text(">refA\n" + "\n".join(ref[i:i + 60] for i in range(0, len(ref), 60)) + "\n")
    r1, r2 = tmp_path / "reads_1.fastq", tmp_path / "reads_2.fastq"
    with open(r1, "w") as f1, open(r2, "w") as f2:
        pairs = [(sample[i:i + 100], sample[i + 200:i + 300]) for i in range(0, len(sample) - 300, 15)]
        pairs.append((sample[5000:5100], sample[9000:9100]))  # a 4.1-kb insert: "larger than expected"
        for n, (a, b) in enumerate(pairs):
            f1.write(f"@p{n}/1\n{a}\n+\n{'I' * 100}\n")
            f2.write(f"@p{n}/2\n{reverse_complement(b)}\n+\n{'I' * 100}\n")

    win = MainWindow()
    win.resize(1400, 900)
    win.show()
    assert win.load(str(fa))
    ed, gv = win.editor, win.genome_view
    out = tmp_path / "mapped.bam"

    def fill(dialog):
        dialog.reads1.setText(str(r1))
        dialog.reads2.setText(str(r2))
        dialog.out.setText(str(out))
        return QDialog.Accepted

    def wait(dialog):
        t = time.perf_counter()
        while dialog.button.text() != "Close" and time.perf_counter() - t < 120:
            QTest.qWait(50)
        print("\n" + dialog.log.toPlainText()[-600:])
        return QDialog.Accepted if dialog.ok else QDialog.Rejected

    monkeypatch.setattr(runner.MapReadsDialog, "exec", fill)
    monkeypatch.setattr(runner.RunDialog, "exec", wait)
    t = time.perf_counter()
    entry = win.map_reads()
    print(f"  bwa index + mem + sort + index: {time.perf_counter() - t:.1f} s")
    assert entry is not None and out.exists() and Path(str(out) + ".bai").exists()
    assert ed.alignments == [entry] and entry.seqid == "refA"
    folder = tmp_path / "mapped_run"
    assert (folder / "ref.fa").exists() and "bwa mem" in (folder / "command.txt").read_text()

    # Zoomed on the SNP: every read over it shows it; normal pairs are grey.
    gv.set_view(snp - 150, snp + 150)
    gv.canvas.grab()
    laid, more = entry.cache[3], entry.cache[4]
    covering = [r for _, r in laid if r[0] <= snp < r[1]]
    assert len(covering) >= 10 and not more
    assert all((snp, alt) in r[5] for r in covering)
    assert set(gv.canvas.read_kinds) == {"normal"}
    pile = entry.bam.pileup("refA", snp, snp + 1)[0]
    assert pile["ACGT".index(alt)] == len(covering) and alignments.flag_level(len(covering), sum(pile[:5])) == "consensus"
    shot(win, "bam_reads")
    # Condensed reads: more rows fit in the same height; log scale for the coverage.
    gv.canvas.set_style(entry.name + " reads", mode="condensed")
    gv.canvas.set_style(entry.name, scale="log2")
    gv.canvas.grab()
    assert entry.cache[6] > 20
    gv.canvas.set_style(entry.name + " reads", mode="expanded")
    gv.canvas.set_style(entry.name, scale="linear")

    # The variant scan flags the SNP and the deletion as consensus changes (red marks).
    t = time.perf_counter()
    while entry.variants is None and time.perf_counter() - t < 20:
        QTest.qWait(20)
    gv.canvas.grab()
    found = {(v.kind, v.pos, v.level) for v in entry.variants}
    assert ("snp", snp, "consensus") in found
    del_at = next(v for v in entry.variants if v.kind == "del")
    assert abs(del_at.pos - dele) <= 5 and len(del_at.ref) == 5 and del_at.level == "consensus"
    assert len(entry.variants) == 2
    assert any(a[0] == "read" and a[3].startswith(f"{snp + 1:,} ") for _, a in gv.canvas.hits)  # its red mark

    # Clicking a read selects it and shows its name and insert size.
    rect, action = next((r, a) for r, a in gv.canvas.hits if a[0] == "read" and "MAPQ" in a[3])
    QTest.mouseClick(gv.canvas, Qt.LeftButton, pos=rect.center().toPoint())
    assert (ed.start, ed.end) == action[1:3] and "MAPQ" in gv.info.text() and "insert" in gv.info.text()

    # The odd pair: red (larger than expected); by strand once the setting says so.
    gv.set_view(4900, 5300)
    gv.canvas.grab()
    assert gv.canvas.read_kinds.get("larger") == 1
    tool_settings.set_value(alignments.TOOL, "color_by", "Strand")
    tool_settings.notifier.changed.emit(alignments.TOOL)
    gv.canvas.grab()
    assert set(gv.canvas.read_kinds) == {"forward", "reverse"}
    tool_settings.set_value(alignments.TOOL, "color_by", alignments.COLOR_MODES[0])
    tool_settings.notifier.changed.emit(alignments.TOOL)

    # Variants track: each change is a box under the coverage; annotating, then applying both consensus changes.
    gv.set_view(snp - 200, snp + 200)
    t = time.perf_counter()   # the colour setting above started a new background scan: wait until its boxes are drawn
    while not any(a[0] == "variant" for _, a in gv.canvas.hits) and time.perf_counter() - t < 20:
        gv.canvas.invalidate()
        gv.canvas.grab()
        QTest.qWait(20)
    assert any(a[0] == "variant" for _, a in gv.canvas.hits)
    consensus = [v for v in entry.variants if v.level == "consensus"]
    assert annotate_variants(ed, entry, consensus) == 2
    assert sum(f.kind == "variation" for f in ed.features) == 2
    ed.undo()
    ed.dna_locked = True
    assert apply_variants(ed, entry, consensus) == 0 and ed.doc.sequence() == ref   # locked: refused
    ed.dna_locked = False
    assert apply_variants(ed, entry, consensus) == 2
    assert ed.length == 50_000 - 5 and ed.doc.sequence(snp, snp + 1) == alt
    assert ed.doc.sequence() == sample  # the plasmid map now matches the sample
    ed.undo()
    ed.undo()
    assert ed.doc.sequence() == ref

    # Whole sequence (wider than the reads window): coverage from one background pass.
    gv.show_all()
    gv.canvas.grab()
    t = time.perf_counter()
    while not entry.bam.has_depth("refA") and time.perf_counter() - t < 20:
        QTest.qWait(20)
    gv.canvas.grab()
    shot(win, "bam_coverage")
    mean = sum(entry.bam.coverage("refA", 0, 50_000, 10)) / 10
    assert 12 < mean < 14, mean  # 2 x 100 bases every 15 bp

    # Cancel stops the command inside WSL too.
    run = runner.RunDialog(win, "sleep", ["sleep 30"])
    t = time.perf_counter()
    while run.pid is None and time.perf_counter() - t < 20:
        QTest.qWait(20)
    run.cancel()
    while run.button.text() != "Close" and time.perf_counter() - t < 20:
        QTest.qWait(20)
    assert not run.ok and time.perf_counter() - t < 20
    run.close()
    win.close()


def test_genome_is_responsive(app, errors):
    win = open_window(app, GENOME_DNA)
    ed, seq_view = win.editor, win.sequence_view
    shot(win, "genome")
    print(f"  lanes: {ed.lane_count}, rows: {seq_view.row_count():,}")

    sb = seq_view.verticalScrollBar()
    t = time.perf_counter()
    for value in (sb.maximum() // 3, sb.maximum() // 2, sb.maximum()):
        sb.setValue(value)
        seq_view.viewport().grab()
    print(f"  3 jumps + repaints: {(time.perf_counter() - t) * 1000:.0f} ms")
    shot(win, "genome_end")

    ed.set_cursor(1000)
    t = time.perf_counter()
    QTest.keyClicks(seq_view, "A")
    app.processEvents()
    print(f"  one typed base (edit + relayout): {(time.perf_counter() - t) * 1000:.0f} ms")
    assert ed.length == 4_631_470
    ed.undo()

    # Genome-sized: primers are matched in a background thread after a short pause.
    t = time.perf_counter()
    db().add("genomeP", ed.doc.sequence(500_000, 500_022))
    while not any(s.label == "genomeP" for s in ed.primer_sites) and time.perf_counter() - t < 60:   # (about 5 s alone; slower when other tests ran first: it only has to arrive, not be quick)
        QTest.qWait(50)
    print(f"  {len(db().rows)} primers matched on the genome in {time.perf_counter() - t:.2f} s (background)")
    # (near matches of its 3′ end elsewhere may be listed too, in any order)
    assert any(s.spans == [(500_000, 500_022)] for s in ed.primer_sites if s.label == "genomeP")
