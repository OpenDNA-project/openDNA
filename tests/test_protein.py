"""The Protein view: mass, isoelectric point, numbering, following the selection, selecting residues selects codons."""
import os
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ["OPENDNA_DATA_DIR"] = tempfile.mkdtemp(prefix="opendna_test_")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "app"))

import pytest
from PySide6.QtWidgets import QApplication

from opendna_core import Document

import protein
from assembly_panel import FragmentEditor

app = QApplication.instance() or QApplication([])
# 5' flank, M K W H stop on the forward strand (CDS 11..22), a feature inside it, a 3' flank
SEQ = "CCCCCCCCCC" + "ATGAAATGGCATTAA" + "GGGGGGGGGG"


def make():
    ed = FragmentEditor(Document.new_sequence("t", SEQ, False), Path("t.gb"))
    ed.add_features([("CDS", "11..25", [("label", "tiny"), ("codon_start", "1"), ("db_xref", "PDB:1abc")]),
                     ("misc_feature", "17..22", [("label", "motif")])])
    return ed


def test_mass_pi_charge_and_layout():
    """Mass, isoelectric point, charge, numbered layout and the index <-> text offset round trip."""
    # Gly-Gly: 2 x 57.0519 + water; a lysine-only chain is basic, an aspartate-only chain acidic
    assert protein.mass("GG") == pytest.approx(132.1191, abs=1e-3)
    assert protein.isoelectric_point("KKKKKK") > 10 and protein.isoelectric_point("DDDDDD") < 4
    assert protein.charge("KKKKKK", 7.0) > 5 > 0 > protein.charge("DDDDDD", 7.0)
    assert protein.until_stop("MK*AA") == "MK"
    # --- the layout and the index round trip
    text = protein.layout("A" * 120)
    lines = text.split("\n")
    assert len(lines) == 3 and lines[0].startswith("     1  AAAAAAAAAA AAAAAAAAAA") and lines[1].lstrip().startswith("51")
    for i in (0, 9, 10, 49, 50, 119):
        line, off = protein.offset_of(i)
        assert protein.index_at(off, line) == i, i


def test_follows_the_selection_and_selecting_residues_selects_codons():
    """The view follows the selection and pins a CDS; selecting residues selects their codons."""
    # --- follows the selection and pins a CDS
    ed = make()
    ed.set_cursor(0)
    view = protein.ProteinView(ed)
    assert "Select a CDS" in view.title.text()
    cds = next(f for f in ed.features if f.kind == "CDS")
    ed.select_feature(cds.index)
    view.refresh()
    assert "MKWH*" in view.text.toPlainText().replace(" ", "") and "4 aa" in view.props.text() and "tiny" in view.title.text()
    assert not view.frame.isEnabled()
    assert view.features.count() == 1 and "aa 3–4" in view.features.item(0).text()
    # a feature inside the CDS keeps the protein shown and highlights its residues
    motif = next(f for f in ed.features if f.kind == "misc_feature")
    ed.select_feature(motif.index)
    view.refresh()
    assert "tiny" in view.title.text() and view.text.marked == {2, 3} and view.text.extraSelections() == [], \
        "W and H are on the Selection track; the letters themselves stay plain"
    # a stretch of DNA outside it is translated in the chosen frame
    ed.select(0, 9)
    view.refresh()
    assert "Translation of 1–9" in view.title.text() and view.frame.isEnabled() and "PPP" in view.text.toPlainText()
    view.frame.setCurrentIndex(3)                       # −1: the reverse complement GGGGGGGGG -> GGG
    assert "GGG" in view.text.toPlainText().replace(" ", "")
    # an edit that changes the protein is shown (the CDS is selected again: the old one is no longer current)
    ed.select_feature(cds.index)
    view.refresh()
    ed.set_cursor(13)
    ed.type_bases("T")
    view.refresh(True)
    assert view.src is None or view.src["version"] == ed.version
    # --- selecting residues selects codons
    ed = make()
    view = protein.ProteinView(ed)
    ed.select_feature(next(f for f in ed.features if f.kind == "CDS").index)
    view.refresh()
    cursor = view.text.textCursor()
    block = view.text.document().firstBlock().position()
    line, off = protein.offset_of(1)                      # K
    cursor.setPosition(block + off)
    cursor.setPosition(block + protein.offset_of(2)[1] + 1, cursor.MoveMode.KeepAnchor)   # through W
    view.text.setTextCursor(cursor)
    assert (ed.start, ed.end) == (13, 19)                 # codons AAA TGG
    assert "tiny" in view.title.text()                    # the CDS stays shown


def test_structure_in_3d_follows_the_selection_both_ways():
    """A structure got for the CDS shows under its letters, the selection is highlighted on it, a click on a residue selects its
    codons (Shift+click: up to it), and another protein hides it."""
    from PySide6.QtCore import QPoint, Qt
    from PySide6.QtTest import QTest
    import structure
    ed = make()
    view = protein.ProteinView(ed)
    view.resize(500, 700)
    view.show()
    ed.select_feature(next(f for f in ed.features if f.kind == "CDS").index)
    view.refresh()
    assert not view.view3d.isVisible()
    assert view.src["links"] == ["1ABC"], "the feature's PDB links go to the 3D structure window"
    atoms =[("CA", 3.8 * i, i % 2, 0.0, 90.0 - 20 * i, res) for i, res in enumerate(("MET", "LYS", "TRP", "HIS"))] + [("CB", 7.6, 1.5, 0.0, 50.0, "TRP")]
    model = "\n".join(f"ATOM  {n + 1:5d}  {name:<3s} {res} A{['MET', 'LYS', 'TRP', 'HIS'].index(res) + 1:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00{b:6.2f}"
                      for n, (name, x, y, z, b, res) in enumerate(atoms))
    view._got_structure("MKWH", structure.Structure(model, "ESMFold", "test model", True))
    canvas = view.view3d.canvas
    assert view.view3d.isVisible() and view.view3d.mode.currentText() == "Confidence"
    assert canvas.highlight == set(), "the whole CDS is selected: nothing to point out"
    ed.select_feature(next(f for f in ed.features if f.kind == "misc_feature").index)   # covers W and H
    view.refresh()
    assert canvas.highlight == {2, 3}
    view.grab()   # draws the tube, the highlighted residues and their atoms
    # --- a click on a residue selects its codon, Shift+click selects up to another
    def click(k, modifier=Qt.NoModifier):
        x, y, _ = canvas.project(canvas.cas[k])
        QTest.mouseClick(canvas, Qt.LeftButton, modifier, QPoint(round(x), round(y)))
    click(1)                                              # K
    assert (ed.start, ed.end) == (13, 16) and canvas.highlight == {1}
    click(3, Qt.ShiftModifier)                            # through H
    assert (ed.start, ed.end) == (13, 22) and canvas.highlight == {1, 2, 3}
    # --- another protein: no structure for it, the view hides
    ed.select(0, 9)
    view.refresh()
    assert not view.view3d.isVisible()


def test_tracks_under_the_letters():
    """Tracks under each line of letters: the selection always, the others switched with the buttons (and remembered); the model's
    tracks only when the protein has a 3D model; each cell says its value on hover."""
    from PySide6.QtCore import QPoint
    import protein_tracks
    import structure
    from structure_view import aa_colour
    from protein_tracks import GAP, PAD, ROW
    assert aa_colour("Hydrophobicity", "I") == "#e8544e" and aa_colour("Hydrophobicity", "R") == "#5b9be0"
    assert aa_colour("Size", "G") == "#fff6d5" and aa_colour("Size", "W") == "#a0581c"
    assert aa_colour("Charge", "A") is None and aa_colour("Property", "*") is None
    ed = make()
    view = protein.ProteinView(ed)
    view.resize(700, 500)
    view.show()
    ed.select_feature(next(f for f in ed.features if f.kind == "CDS").index)
    view.refresh()
    text, buttons = view.text, view.tracks.buttons
    assert text.tracks == ["Selection", "Hydrophobicity", "Size"], "Structure is on, but there is no model yet"
    assert not buttons["Structure"].isEnabled() and "Get a 3D structure first" in buttons["Structure"].toolTip()
    room = text.document().firstBlock().blockFormat().bottomMargin()
    assert room == PAD + 3 * (ROW + GAP), "room under each line for the three tracks"
    # --- switching a track on makes room for it and is remembered
    buttons["Charge"].setChecked(True)
    assert text.tracks == ["Selection", "Hydrophobicity", "Size", "Charge"] and "Charge" in protein_tracks.saved_tracks()
    assert text.document().firstBlock().blockFormat().bottomMargin() == PAD + 4 * (ROW + GAP)
    buttons["Size"].setChecked(False)
    assert "Size" not in text.tracks and "Size" not in protein_tracks.saved_tracks()
    # --- hovering a cell: the track and the amino acid under it
    xs, w, y0 = text._rows(text.document().firstBlock())

    def at(track, i):
        return text.track_at(QPoint(round(xs[i] + w / 2), round(y0 + text.tracks.index(track) * (ROW + GAP) + ROW / 2)))
    assert at("Hydrophobicity", 2) == ("Hydrophobicity", 2) and at("Charge", 1) == ("Charge", 1)
    assert text.track_at(QPoint(round(xs[0]), round(y0 - 8))) is None, "over the letters: no track"
    assert text.describe("Hydrophobicity", 2) == "Hydrophobicity · W3: -0.9 (Kyte-Doolittle)"
    assert text.describe("Charge", 1) == "Charge · K2: positive" and text.describe("Size", 2) == "Size · W3: 228 Å³"
    view.grab()   # draws every track
    # --- a model: its structure track, and its confidence when it is a prediction (M and K only are in it)
    model = "\n".join(f"ATOM  {i + 1:5d}  CA  {res} A{i + 1:4d}    {3.8 * i:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00{b:6.2f}"
                      for i, (res, b) in enumerate((("MET", 95.0), ("LYS", 40.0))))
    view._got_structure("MKWH", structure.Structure(model, "ESMFold", "half a model", True))
    assert buttons["Structure"].isEnabled() and "Structure" in text.tracks and text.model == {0: ("-", 95.0), 1: ("-", 40.0)}
    buttons["Confidence"].setChecked(True)
    assert text.tracks[:3] == ["Selection", "Structure", "Confidence"]
    assert text.describe("Confidence", 1) == "Confidence · K2: pLDDT 40" and text.describe("Structure", 3) == "Structure · H4: not in the model"
    view.grab()
    # --- another protein, without a model: its tracks go
    ed.select(0, 9)
    view.refresh()
    assert "Structure" not in text.tracks and not buttons["Structure"].isEnabled() and text.model == {}


def test_a_pane_in_the_window(monkeypatch):
    from main import MainWindow
    errors = []
    monkeypatch.setattr(sys, "excepthook", lambda *exc: errors.append(exc))
    win = MainWindow()
    win.show()
    win.show_document(Document.new_sequence("t", SEQ, False), Path("t.gb"))
    view = win.workspace.ensure(protein.ProteinView)
    win.editor.select(10, 25)
    view.refresh()
    assert "4 aa + stop" in view.title.text(), "a stop codon is not an amino acid: the header agrees with the properties line"
    win.grab()
    win.editor.doc.modified = False
    win.close()
    assert not errors, errors
