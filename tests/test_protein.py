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
    ed.add_features([("CDS", "11..25", [("label", "tiny"), ("codon_start", "1")]), ("misc_feature", "17..22", [("label", "motif")])])
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
    assert "tiny" in view.title.text() and [s.cursor.selectedText() for s in view.text.extraSelections()] == ["WH"]
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
    atoms = [("CA", 3.8 * i, i % 2, 0.0, 90.0 - 20 * i, res) for i, res in enumerate(("MET", "LYS", "TRP", "HIS"))] + [("CB", 7.6, 1.5, 0.0, 50.0, "TRP")]
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


def test_letter_colours():
    """The letters are coloured by property, hydrophobicity, size or charge, or like the 3D view; the 3D view has the same schemes."""
    import structure
    from structure_view import PROPERTY, aa_colour
    assert aa_colour("Hydrophobicity", "I") == "#e8544e" and aa_colour("Hydrophobicity", "R") == "#5b9be0"
    assert aa_colour("Size", "G") == "#fff6d5" and aa_colour("Size", "W") == "#a0581c"
    assert aa_colour("Charge", "A") is None and aa_colour("Property", "*") is None
    ed = make()
    view = protein.ProteinView(ed)
    view.show()
    ed.select_feature(next(f for f in ed.features if f.kind == "CDS").index)
    view.refresh()

    def backgrounds():
        """The background colour of each letter of the first line."""
        block = view.text.document().firstBlock()
        return {block.text()[r.start]: r.format.background().color().name() for r in block.layout().formats()}
    assert backgrounds() == {}, "plain letters to begin with"
    view.colour_by.setCurrentText("Property")
    assert backgrounds() == {"M": PROPERTY["M"].lower(), "K": PROPERTY["K"].lower(), "W": PROPERTY["W"].lower(), "H": PROPERTY["H"].lower()}
    view.colour_by.setCurrentText("Charge")
    assert set(backgrounds()) == {"K", "H"}
    # --- like the 3D view: only the residues in the model, in its colours
    view.colour_by.setCurrentText("Same as 3D")
    assert backgrounds() == {}, "no structure: nothing to copy"
    model = "\n".join(f"ATOM  {i + 1:5d}  CA  {res} A{i + 1:4d}    {3.8 * i:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00{b:6.2f}"
                      for i, (res, b) in enumerate((("MET", 95.0), ("LYS", 40.0))))   # M and K only
    view._got_structure("MKWH", structure.Structure(model, "ESMFold", "half a model", True))
    assert backgrounds() == {"M": "#0053d6", "K": "#ff7d45"}, "pLDDT colours, as on the model"
    view.view3d.mode.setCurrentText("Size")
    assert backgrounds() == {"M": aa_colour("Size", "M"), "K": aa_colour("Size", "K")}
    assert view.view3d.canvas.colour(0) == aa_colour("Size", "M")


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
