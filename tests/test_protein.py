"""The Protein view: mass, isoelectric point, numbering, following the selection, selecting residues selects codons."""
import os
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ["GENOMEEDIT_DATA_DIR"] = tempfile.mkdtemp(prefix="genomeedit_test_")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "app"))

import pytest
from PySide6.QtWidgets import QApplication

from genomeedit_core import Document

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
