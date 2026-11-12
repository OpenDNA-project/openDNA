"""Import / Export menus, the Open DNA Project file (.odp), the DNA lock and reference genomes from NCBI Datasets (network faked)."""
import io
import json
import os
import sys
import tempfile
import time
import zipfile
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ["OPENDNA_DATA_DIR"] = tempfile.mkdtemp(prefix="opendna_test_")
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "app"))

import pytest
from PySide6.QtCore import Qt
from PySide6.QtTest import QTest
from PySide6.QtWidgets import QApplication, QFileDialog, QMessageBox

import exports
import ncbi
from main import MainWindow

app = QApplication.instance() or QApplication([])
DEMO = ROOT / "demo" / "demo_vector.gb"


@pytest.fixture
def win():
    w = MainWindow()
    assert w.load(DEMO)
    yield w
    w.editor.doc.modified = False
    w.close()


def test_exports(win, tmp_path, monkeypatch):
    """Text exports (GFF3, BED, feature table, proteins, FASTA), the submission pair, the FASTA and GenBank copies, and the
    Export menu greying what is missing."""
    # --- text exports ---
    ed = win.editor
    gff = exports.gff3_text(ed)
    assert gff.startswith("##gff-version 3") and gff.count("\n") >= len(ed.features) + 2
    cds = exports.cds_features(ed)
    assert cds, "the demo vector has a CDS"
    f = cds[0]
    line = next(l for l in gff.splitlines() if l.split("\t")[2:3] == ["CDS"])
    assert line.split("\t")[6] == ("+" if f.strand > 0 else "-")
    assert len(exports.bed_text(ed).splitlines()) == sum(len(f.spans) for f in ed.features)
    tbl = exports.tbl_text(ed).splitlines()
    assert tbl[0] == f">Feature {exports.seq_id(ed)}"
    first = next(l for l in tbl if l.endswith("\tCDS"))
    a, b = (int(x.strip("<>")) for x in first.split("\t")[:2])
    assert (a < b) == (f.strand > 0), "the table lists minus-strand parts as start > stop"
    assert not any("ApEinfo" in l or "\tlabel\t" in l for l in tbl)
    assert exports.proteins_fasta(ed).startswith(">")
    # A real CDS on each strand (the demo files are random sequences).
    from editor import Editor, reverse_complement
    from opendna_core import Document
    gene = "ATGAAACCCGGGTAA"
    built = Editor(Document.new_sequence("t", "TT" + gene + "TT" + reverse_complement(gene) + "TT", False), "t.gb")
    built.add_features([("CDS", "3..17", [("label", "fwd")]), ("CDS", "complement(20..34)", [("label", "rev")])])
    assert [exports.protein_of(built, f) for f in exports.cds_features(built)] == ["MKPG", "MKPG"]
    lines = exports.tbl_text(built).splitlines()
    assert "3\t17\tCDS" in lines and "34\t20\tCDS" in lines and "\t\t\tnote\trev" in lines
    assert exports.translate("ATGAAATAA") == "MK*"
    assert exports.fasta_text("x", "A" * 75).splitlines()[1:] == ["A" * 70, "A" * 5]
    # --- submission pair ---
    fsa, tbl = exports.write_submission(win.editor, tmp_path / "sub.tbl")
    assert fsa.read_text().startswith(f">{exports.seq_id(win.editor)}") and tbl.exists()
    # --- FASTA and GenBank copy ---
    target = tmp_path / "out.fasta"
    monkeypatch.setattr(QFileDialog, "getSaveFileName", lambda *a, **k: (str(target), ""))
    assert win.export_fasta(False) == target
    assert target.read_text().replace("\n", "").endswith(win.editor.doc.sequence()[-20:])
    target = tmp_path / "copy.gb"
    assert win.export_genbank() == target and target.read_text().startswith("LOCUS")
    assert win.editor.path == DEMO, "an export is a copy: the open file keeps its name"
    # --- the Export menu greys what is missing ---
    win._update_export_menu()
    assert win.export_whole.menuAction().isEnabled() and win.export_features.menuAction().isEnabled()
    assert not win.export_alignments.menuAction().isEnabled(), "no BAM loaded"
    assert not win.export_selection.isEnabled()
    win.editor.select(0, 30)
    win._update_export_menu()
    assert win.export_selection.isEnabled() and win.export_translation.isEnabled()
    empty = MainWindow()
    empty._update_export_menu()
    assert not any(m.menuAction().isEnabled() for m in (empty.export_whole, empty.export_dna, empty.export_features,
                                                         empty.export_protein, empty.export_alignments))
    empty.close()


def test_dna_lock_and_project_file(win, tmp_path):
    """The DNA lock refuses base edits (not feature edits, no linearising) and unlocks again; a .odp project round-trips
    the sequence, tracks, display settings and the lock."""
    # --- DNA lock ---
    ed = win.editor
    seq, n = ed.doc.sequence(), len(ed.features)
    ed.set_cursor(5)
    assert ed.type_bases("ACGT")       # unlocked: a base edit
    win.set_dna_locked(True)
    messages = []
    ed.error.connect(messages.append)
    assert not ed.type_bases("A") and not ed.backspace() and not ed.paste("ACGT")
    assert messages and "locked" in messages[0]
    ed.select(0, 20)
    ed.add_features([("misc_feature", "1..20", [("label", "ok")])])   # features still change
    assert len(ed.features) == n + 1
    assert ed.undo()                    # the feature: allowed
    assert not ed.undo()                # the typed bases: refused while locked
    assert ed.doc.sequence() != seq
    assert not win.set_circular(False) and ed.circular == ed.doc.circular   # linearising moves base 1
    win.set_dna_locked(False)
    assert ed.undo() and ed.doc.sequence() == seq
    # --- project round trip ---
    ed = win.editor
    gff = tmp_path / "extra.gff3"
    gff.write_text(f"{ed.doc.name}\tx\tgene\t10\t50\t.\t+\t.\tName=g1\n")
    assert win.load_track(str(gff)) is not None
    ed.set_display(enzymes=True, hidden_kinds={"misc_feature"})
    win.set_dna_locked(True)
    project = tmp_path / "work.odp"
    assert win._write(project) and ed.path == project and not ed.doc.modified
    with zipfile.ZipFile(project) as z:
        meta = json.loads(z.read("project.json"))
        assert set(z.namelist()) == {"sequence.gb", "project.json"}
    assert meta["dna_locked"] and meta["tracks"][0]["relative"] == "extra.gff3"

    other = MainWindow()
    assert other.load(project)
    o = other.editor
    assert o.doc.sequence() == ed.doc.sequence() and len(o.features) == len(ed.features)
    assert o.dna_locked and other.lock_action.isChecked() and other.lock_action.iconText() == "DNA"
    assert [t.name for t in o.tracks] == ["extra.gff3"] and o.show_enzymes and o.hidden_kinds == {"misc_feature"}
    other.close()


# ---- reference genomes ---------------------------------------------------------------------------------------------------
REPORT = {"reports": [{"accession": "GCF_000000001.1", "organism": {"organism_name": "Madeup coli"},
                       "assembly_info": {"assembly_name": "ASM1v1", "assembly_level": "Complete Genome", "release_date": "2026-01-01",
                                         "refseq_category": "reference genome"},
                       "assembly_stats": {"total_sequence_length": "60", "number_of_component_sequences": 1}}], "total_count": 1}


def fake_zip():
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as z:
        z.writestr("README.md", "x")
        z.writestr("ncbi_dataset/data/GCF_000000001.1/GCF_000000001.1_ASM1v1_genomic.fna", ">chr1 Madeup\n" + "ACGT" * 15 + "\n")
        z.writestr("ncbi_dataset/data/GCF_000000001.1/genomic.gff", "chr1\tRefSeq\tgene\t5\t20\t.\t+\t.\tName=abcA\n")
    return buf.getvalue()


def test_reference_genomes(monkeypatch):
    """Searching NCBI Datasets, downloading (and keeping) a genome, and the dialog opening it with its GFF track and using it
    as the off-target reference (network faked)."""
    # --- search and download ---
    urls = []
    monkeypatch.setattr(ncbi, "fetch", lambda url, timeout=60: urls.append(url) or json.dumps(REPORT))
    found = ncbi.search_genomes("Madeup coli")
    assert "taxon/Madeup%20coli/dataset_report" in urls[0] and "reference_only=true" in urls[0]
    a = found[0]
    assert (a.accession, a.length, a.sequences, a.level) == ("GCF_000000001.1", 60, 1, "Complete Genome")
    ncbi.search_genomes("gcf_000000001.1")
    assert "accession/GCF_000000001.1/dataset_report" in urls[1]
    with pytest.raises(ncbi.NcbiError):
        ncbi.search_genomes("  ")

    calls, sizes = [], []
    monkeypatch.setattr(ncbi, "open_url", lambda url, headers, timeout=120: calls.append(url) or io.BytesIO(fake_zip()))
    fasta, gff, cached = ncbi.download_genome(a, progress=sizes.append)
    assert not cached and sizes and fasta.read_text().startswith(">chr1") and "abcA" in gff.read_text()
    assert "GENOME_FASTA" in calls[0] and "GENOME_GFF" in calls[0]
    assert not list(fasta.parent.glob("*.part"))
    assert ncbi.download_genome(a)[2] and len(calls) == 1, "kept: not downloaded twice"
    # --- dialog opens the genome with its track ---
    from ncbi_dialog import GenomeDialog
    from primer_db import settings
    # Setting a reference genome writes a shared setting that the primer tests of this run also read: put it back.
    before = settings().value("primers/genome", "")
    monkeypatch.setattr(ncbi, "fetch", lambda url, timeout=60: json.dumps(REPORT))
    monkeypatch.setattr(ncbi, "open_url", lambda url, headers, timeout=120: io.BytesIO(fake_zip()))
    monkeypatch.setattr(QMessageBox, "question", lambda *a, **k: QMessageBox.Yes)
    w = MainWindow()
    dialog = GenomeDialog(w)
    dialog.query.setText("Madeup coli")
    dialog.search()
    t = time.perf_counter()
    while not dialog.assemblies and time.perf_counter() - t < 10:
        QTest.qWait(20)
    assert dialog.reference_button.isEnabled(), "one sequence: usable for the off-target checks"
    dialog.open()
    t = time.perf_counter()
    while dialog.opened is None and time.perf_counter() - t < 10:
        QTest.qWait(20)
    assert dialog.opened is w and w.editor.length == 60 and [t.name for t in w.editor.tracks] == ["GCF_000000001.1.gff"]
    dialog.use_as_reference()
    t = time.perf_counter()
    while "reference genome" not in dialog.status.text() and time.perf_counter() - t < 10:
        QTest.qWait(20)
    from primer_panel import genome_path
    try:
        assert genome_path().endswith("GCF_000000001.1.fna") and Path(genome_path()).exists()
    finally:
        settings().setValue("primers/genome", before)
        settings().sync()
    w.close()


def test_example_files_and_feature_direction(win):
    """The example plasmids and genome bundled in demo/ open (read only); a .dna file marks features without direction and
    they come in with no arrow; a new misc_feature has no arrow by default and the direction can be chosen."""
    # --- example files ---
    from opendna_core import Document
    for name, length in [("pBAD30_eilR.gb", 5936), ("pFG051_Chlor.gb", 3780), ("BW25113_2026.gb", None)]:
        doc = Document(str(ROOT / "demo" / name))
        assert doc.features() and (doc.length == length if length else doc.length > 4_000_000), name
    # --- direction read from a .dna file ---
    from sample_files import DNA_FILE
    from editor import arrow
    from opendna_core import Document
    if DNA_FILE.exists():   # (the file is optional)
        feats = Document(str(DNA_FILE)).features()
        cds = [f for f in feats if f.kind == "CDS"]
        assert cds and all(arrow(f) != 0 for f in cds), "CDSs keep their arrow"
        assert any(arrow(f) == 0 for f in feats if f.kind == "misc_feature")
    # --- features without direction ---
    from editor import arrow, is_directional
    from feature_dialog import FeatureDialog
    assert not is_directional("misc_feature", []) and is_directional("CDS", [])
    assert is_directional("misc_feature", [("directional", "yes")]) and not is_directional("promoter", [("directional", "no")])
    ed = win.editor
    ed.select(10, 40)
    new = FeatureDialog(win, ed)
    assert new.direction.currentText().startswith("None"), "a new misc_feature has no arrow by default"
    new.name.setText("plain")
    new.accept()
    plain = next(f for f in ed.features if f.label == "plain")
    assert arrow(plain) == 0 and not any(k == "directional" for k, _ in plain.qualifiers), "the default is not written"
    edit = FeatureDialog(win, ed, plain.index)
    edit.direction.setCurrentIndex(1)
    edit._direction_chosen(1)
    edit.accept()
    plain = next(f for f in ed.features if f.label == "plain")
    assert arrow(plain) == -1 and ("directional", "yes") in plain.qualifiers and plain.location.startswith("complement(")
    ed.select(50, 80)
    cds = FeatureDialog(win, ed)
    cds.kind.setCurrentText("CDS")
    assert cds.direction.currentText().startswith("Forward")
    cds.direction.setCurrentIndex(2)
    cds._direction_chosen(2)
    cds.name.setText("arrowless")
    cds.accept()
    assert arrow(next(f for f in ed.features if f.label == "arrowless")) == 0
    win.feature_table.view.grab()   # drawn with "—"
    win.sequence_view.grab()
    win.map_view.grab()


# ---- primers hidden in one file ----------------------------------------------------------------------------------------
def test_primer_display_and_hiding(win):
    """'Hide in This File' hides a primer from the drawing but keeps it found, listed greyed and remembered per file; the
    Display filter > Primers 'From this sequence' and 'From the database' switch independently (a primer in both shows while
    either is on) and a file-only primer can be hidden too."""
    # --- hide a primer in this file ---
    from editor import Editor
    from primer_db import db
    from primer_panel import PrimerPanel
    ed = win.editor
    ed.set_display(db_primers=True)   # hidden by default
    db().add("hideme", ed.doc.sequence(100, 122))
    ed.match_primers()
    ed._match_in_background()
    t = time.perf_counter()
    while not any(s.label == "hideme" for s in ed.primer_sites) and time.perf_counter() - t < 10:
        QTest.qWait(20)
    assert any(s.label == "hideme" for s in ed.drawn_primers)
    ed.set_primer_hidden("hideme", True)
    assert not any(s.label == "hideme" for s in ed.drawn_primers), "not drawn"
    assert any(s.label == "hideme" for s in ed.primer_sites) and db().find("hideme"), "still found and still in the database"
    panel = PrimerPanel(ed)
    item = next(panel.tree.topLevelItem(i) for i in range(panel.tree.topLevelItemCount())
                if panel.tree.topLevelItem(i).data(0, Qt.UserRole) == "hideme")
    assert item.text(0) == "hideme (hidden)", "listed greyed in the Primers view"
    assert "hideme" in Editor(ed.doc, ed.path).hidden_primers, "remembered for this file"
    ed.set_primer_hidden("hideme", False)
    assert any(s.label == "hideme" for s in ed.drawn_primers)
    db().remove(["hideme"])
    # --- filter by source ---
    from PySide6.QtWidgets import QMenu
    from display_panel import DisplayPanel
    from feature_table import primer_actions
    from primer_db import db
    ed = win.editor
    own, shared = ed.doc.sequence(300, 322), ed.doc.sequence(500, 522)
    db().add("fromdb", ed.doc.sequence(100, 122))
    db().add("both", shared)
    ed.set_extra([("fromfile", own), ("both", shared)])
    ed.match_primers()
    ed._match_in_background()
    t = time.perf_counter()
    while {"fromdb", "fromfile", "both"} - {s.label for s in ed.primer_sites} and time.perf_counter() - t < 10:
        QTest.qWait(20)

    def drawn():
        return {s.label for s in ed.drawn_primers}
    try:
        assert {"fromdb", "fromfile", "both"} <= drawn()
        ed.set_display(db_primers=False)
        assert {"fromfile", "both"} <= drawn() and "fromdb" not in drawn()
        ed.set_display(db_primers=True, file_primers=False)
        assert drawn() == {"fromdb"}, "a primer in the file and in the database counts as the file's: only 'this sequence' shows it"
        ed.set_display(db_primers=False)
        assert not drawn() and not ed.show_primers
        panel = DisplayPanel(ed)
        panel.file_primers.setChecked(True)
        assert ed.show_file_primers and "fromfile" in drawn()
        ed.set_display(primers=True)
        assert ed.show_file_primers and ed.show_db_primers
        menu = QMenu()
        primer_actions(menu, ed, "fromfile")
        hide = next(a for a in menu.actions() if a.text() == "Hide in This File")
        assert hide.isEnabled(), "a primer that is only in this file can be hidden"
    finally:
        ed.set_display(primers=True)
        ed.set_extra([])
        db().remove(["fromdb", "both"])


# ---- tabs and layout ------------------------------------------------------------------------------------------------------


def test_tab_names_and_standard_layout(win, tmp_path):
    """A tab is named after the panes you can see (hiding the Alignment view renames it, showing it restores the name); every
    file opens with the standard 70 / 30 layout and one tab, only an Open DNA Project keeps a layout of its own."""
    # --- tab named after the visible panes ---
    ws = win.workspace
    i = ws.currentIndex()
    assert ws.tabText(i) == "Sequence · Map · Features"
    win.set_genome_visible(True)
    assert ws.tabText(i).startswith("Alignment") and "Alignment" in ws.tabToolTip(i)
    pane = next(p for p in ws.ordered_panes(ws.widget(i)) if p.view_name == "Alignment")
    win.set_genome_visible(False)
    assert pane.view is not None and not pane.isVisibleTo(ws.widget(i)), "hidden, but still there"
    assert "Alignment" not in ws.tabText(i) and "Alignment" not in ws.tabToolTip(i)
    # --- standard tabs ---
    w = MainWindow()
    w.resize(1400, 850)
    w.show()
    assert w.load(DEMO)
    app.processEvents()
    root = w.workspace.widget(0)
    sizes = root.sizes()
    assert [round(100 * x / sum(sizes)) for x in sizes] == [70, 30]
    inner = root.widget(1)
    assert [round(100 * x / sum(inner.sizes())) for x in inner.sizes()] == [50, 50]
    assert [p.view_name for p in w.workspace.ordered_panes(root)] == ["Sequence", "Map", "Features"]
    # One tab and nothing more: no extra tabs open with a file.
    assert w.workspace.count() == 1

    # Rearranged, saved as a project, reopened: the project keeps its tabs, a plain file does not.
    w.set_genome_visible(True)
    project = tmp_path / "work.odp"
    assert w._write(project)
    again = MainWindow()
    again.resize(1400, 850)
    assert again.load(project)
    assert again.workspace.tabText(0).startswith("Alignment"), "the project opens as it was saved"
    plain = MainWindow()
    plain.resize(1400, 850)
    assert plain.load(DEMO)
    assert plain.workspace.tabText(0) == "Sequence · Map · Features", "a plain file is unaffected"
    for x in (w, again, plain):
        x.editor.doc.modified = False
        x.close()


# ---- Export > Image of the Current Tab ---------------------------------------------------------------------------------
def test_tab_image_export(win, tmp_path, monkeypatch):
    """The image dialog keeps scale and DPI in step. Every pane of the current tab is drawn into the file, at the size and
    DPI chosen, as png, tiff and svg; in an SVG the views are vector paths (the map and alignment views skip their cached
    screen pixmaps while exporting)."""
    import re
    import image_export
    from PySide6.QtGui import QImage
    # --- the dialog follows the DPI ---
    d = image_export.ImageExportDialog(win, win.workspace.currentWidget())
    d.dpi.setValue(600)
    assert d.scale.value() == round(600 / 96 * 100)
    d.scale.setValue(100)
    assert d.options()[1:3] == (1.0, 600)
    assert "cm" in d.info.text()

    # --- exports ---
    win.resize(1000, 700)
    win.show()
    QTest.qWait(50)
    root = win.workspace.currentWidget()
    monkeypatch.setattr(image_export.ImageExportDialog, "exec", lambda self: True)
    monkeypatch.setattr(QFileDialog, "getSaveFileName", lambda *a, **k: (str(path), ""))
    for fmt in ("png", "tiff", "svg"):
        path = tmp_path / f"tab.{fmt}"
        monkeypatch.setattr(image_export.ImageExportDialog, "options", lambda self, fmt=fmt: (fmt, 2.0, 300, False))
        assert win.export_tab_image() == path and path.stat().st_size > 1000, fmt
        assert not image_export.ACTIVE
        if fmt == "svg":
            text = path.read_text(encoding="utf-8")
            sizes = [int(float(w)) for w in re.findall(r'<image[^>]*width="([\d.]+)"', text)]
            assert "<path" in text and all(w <= 64 for w in sizes), "the views are drawn as vectors (only style icons are images)"
            # Text at its screen size inside the 2x group (font sizes in screen pixels), not resized again for 300 DPI.
            fonts = [float(x) for x in re.findall(r'font-size="([\d.]+)"', text)]
            screen = win.font().pointSizeF() * root.logicalDpiY() / 72
            assert fonts and max(fonts) < 2 * screen, (max(fonts), screen)
            width_mm = float(re.search(r'<svg\b[^>]*width="([\d.]+)mm"', text).group(1))
            assert abs(width_mm - root.width() * 2 / 300 * 25.4) < 0.1, "prints at its pixels / DPI"
            texts = re.findall(r"<text [^>]*>", text)
            assert texts and all('fill-rule="nonzero"' in t for t in texts), "glyphs whose outlines overlap stay solid"
        else:
            img = QImage(str(path))
            assert (img.width(), img.height()) == (root.width() * 2, root.height() * 2), fmt
            assert round(img.dotsPerMeterX() * 0.0254) == 300, fmt
            assert len({img.pixel(x, y) for x in range(0, img.width(), 37) for y in range(0, img.height(), 37)}) > 3, fmt


def test_a_project_carries_the_genetic_code_and_the_codon_table_in_use(win, tmp_path):
    """A .odp holds the code and the codon usage table in use: on another machine (here: after deleting them) they are added and
    used; a table or code of the same name with other numbers is kept beside, as "(from project)"."""
    import codon_usage
    import genetic_code
    import tool_settings
    saved = (tool_settings.value(genetic_code.TOOL, "table"), codon_usage.active_name())
    aas, starts = genetic_code.recoded(11, ["AGA"])
    code_id = genetic_code.save_code("Shared recoded host", aas, starts)
    table = codon_usage._table("Shared usage", {c: 10.0 + i for i, c in enumerate(codon_usage.CODONS)})
    codon_usage.save_user(table)
    try:
        tool_settings.set_value(genetic_code.TOOL, "table", f"{code_id} – Shared recoded host")
        codon_usage.set_active("Shared usage")
        project = tmp_path / "shared.odp"
        assert win._write(project)
        with zipfile.ZipFile(project) as z:
            meta = json.loads(z.read("project.json"))
        assert meta["genetic_code"]["name"] == "Shared recoded host" and meta["genetic_code"]["aas"] == aas
        assert meta["codon_usage"]["name"] == "Shared usage" and meta["codon_usage"]["per_thousand"]["TTT"] == 10.0
        # another machine: neither is there, and the standard code is in use
        genetic_code.delete_code(code_id)
        codon_usage.delete_user("Shared usage")
        assert genetic_code.table_id() == 1 and "Shared usage" not in codon_usage.all_tables()
        other = MainWindow()
        assert other.load(project)
        assert genetic_code.TABLES[genetic_code.table_id()].name == "Shared recoded host" and genetic_code.removed_codons(genetic_code.table_id()) == ["AGA"]
        assert codon_usage.active_name() == "Shared usage" and codon_usage.active().per_thousand["TTT"] == 10.0
        assert "Genetic code set to" in other.statusBar().currentMessage() and "Codon usage table set to" in other.statusBar().currentMessage()
        other.close()
        # the same names with other numbers: both are kept
        mine = genetic_code.table_id()
        other_aas, other_starts = genetic_code.recoded(11, ["AGG"])
        genetic_code.delete_code(mine)
        genetic_code.save_code("Shared recoded host", other_aas, other_starts)
        codon_usage.save_user(codon_usage._table("Shared usage", {c: 1.0 for c in codon_usage.CODONS}))
        tool_settings.set_value(genetic_code.TOOL, "table", genetic_code.DEFAULT_LABEL)
        again = MainWindow()
        assert again.load(project)
        names = {t.name for t in genetic_code.TABLES.values()}
        assert {"Shared recoded host", "Shared recoded host (from project)"} <= names
        assert codon_usage.active_name() == "Shared usage (from project)" and codon_usage.active().per_thousand["TTT"] == 10.0
        assert codon_usage.all_tables()["Shared usage"].per_thousand["TTT"] == 1.0   # yours is untouched
        again.close()
    finally:
        for t in [t for t in genetic_code.TABLES.values() if t.id >= genetic_code.FIRST_USER and t.name.startswith("Shared recoded host")]:
            genetic_code.delete_code(t.id)
        for name in [n for n in codon_usage.all_tables() if n.startswith("Shared usage")]:
            codon_usage.delete_user(name)
        tool_settings.set_value(genetic_code.TOOL, "table", saved[0])
        codon_usage.set_active(saved[1])
