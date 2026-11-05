"""Bad input fed to the app's input fields (docs/INPUT_AUDIT.md): each one gives a friendly message or is refused, never a
traceback, a nonsense value or a frozen window."""
import os
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ["OPENDNA_DATA_DIR"] = tempfile.mkdtemp(prefix="opendna_test_")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "app"))

import pytest
from PySide6.QtCore import QEvent, Qt
from PySide6.QtGui import QGuiApplication, QKeyEvent
from PySide6.QtWidgets import QApplication, QDialog, QDialogButtonBox, QLabel, QLineEdit, QMessageBox, QPlainTextEdit

import primer_panel
from editor import clean_pasted_protein, dna_field, protein_field
from feature_dialog import FeatureDialog
from opendna_core import Document, tm
from primer_db import db
from test_cas9 import close, dna, open_doc

app = QApplication.instance() or QApplication([])


@pytest.fixture
def win(monkeypatch):
    w = open_doc(Document.new_sequence("t", dna(500), False))
    w.messages, w.errors = [], []
    monkeypatch.setattr(QMessageBox, "warning", staticmethod(lambda parent, title, text, *a, **k: w.messages.append(text)))
    w.editor.error.connect(w.errors.append)
    yield w
    close(w)


def paste(field, text):
    QGuiApplication.clipboard().setText(text)
    field.setFocus()
    app.sendEvent(field, QKeyEvent(QEvent.KeyPress, Qt.Key_V, Qt.ControlModifier))


# ---- sequence fields --------------------------------------------------------------------------------------------------------

def test_sequence_field_paste():
    """A pasted DNA sequence is cleaned (header, numbers, line breaks) or refused whole; a pasted protein likewise."""
    # --- DNA field ---
    f = dna_field(QLineEdit())
    paste(f, ">seq1 some record\n1 acgt acgt\n9 ACGT\n")
    assert f.text() == "acgtacgtACGT"   # header, numbers and line breaks dropped
    f.clear()
    paste(f, "ACGT is not DNA: hello")
    assert f.text() == ""              # refused whole, with a tooltip
    f.insert("ACGU xyz")               # inserted through the validator: refused
    assert "x" not in f.text().lower()
    # --- protein field ---
    assert clean_pasted_protein(">p\nMKT AYI\n10 AKQ") == "MKTAYIAKQ"
    with pytest.raises(ValueError):
        clean_pasted_protein("MKT123?!")
    f = protein_field(QLineEdit())
    paste(f, ">p\nMKTAYI\nAKQ")
    assert f.text() == "MKTAYIAKQ"


# ---- primer dialog ----------------------------------------------------------------------------------------------------------

def test_primer_dialogs(win, tmp_path):
    """Primer name / sequence problems are named; the primer dialog limits its fields and refuses bad input; the PCR conditions
    never give a nonsense Tm; the options dialog refuses a missing genome file."""
    # --- problem messages ---
    for name, seq, expect in [
        ("", "ACGTACGTAC", "name"),
        ("a\tb", "ACGTACGTAC", "control character"),
        ("ok_short", "ACGTAC", "at least 8"),
        ("ok_long", "A" * 201, "at most"),
        ("ok_bad", "ACGTACGTXX", "not DNA"),
    ]:
        assert expect in primer_panel.primer_problem(name, seq), (name, expect)
    # --- primer dialog ---
    win.editor.select(10, 40)
    d = primer_panel.PrimerDialog(win, win.editor)
    assert d.name.maxLength() == primer_panel.MAX_NAME and d.sequence.maxLength() == primer_panel.MAX_PRIMER
    d.name.setText("x" * 500)
    assert len(d.name.text()) == primer_panel.MAX_NAME
    d.name.setText("   ")
    d.accept()                                     # a blank name: a message, nothing saved, still open
    assert win.messages and "name" in win.messages[-1].lower() and d.result() != QDialog.Accepted
    d.name.setText("good_primer")
    d.sequence.setText("ACGT")
    d.accept()
    assert "at least 8" in win.messages[-1] and db().find("good_primer") is None
    d.sequence.setText("ACGT ACGT\nACGT")          # spaces and a line break are fine
    d.accept()
    assert db().find("good_primer")["sequence"] == "ACGTACGTACGT"
    # a selection too long to be a primer leaves the box empty instead of filling it with a genome
    win.editor.select(0, 400)
    assert primer_panel.PrimerDialog(win, win.editor).sequence.text() == ""
    # --- conditions ---
    for key, minimum in primer_panel.CONDITION_MINIMUM.items():
        s = primer_panel.settings()
        s.setValue(f"pcr/{key}", 0)
        try:
            assert primer_panel.conditions()[key] >= minimum
            assert tm("ATGCATGCATGCATGCATGC", **primer_panel.conditions()) > 0
        finally:
            s.remove(f"pcr/{key}")
    box = primer_panel.OptionsDialog(None).boxes
    assert box["na_mm"].minimum() >= 1 and box["primer_nm"].minimum() >= 1
    # --- options dialog ---
    d = primer_panel.OptionsDialog(win)
    d.genome_check.setChecked(True)
    d.genome.setText(str(tmp_path / "nope.fa"))
    d.accept()
    assert "no file" in win.messages[-1].lower() and d.result() != QDialog.Accepted


# ---- feature dialog -------------------------------------------------------------------------------------------------------

def feature_attempt(win, loc="10..20", qualifiers="", name="x"):
    win.editor.select(10, 50)
    d = FeatureDialog(win, win.editor)
    d.location.setText(loc)
    d.qualifiers.setPlainText(qualifiers)
    d.name.setText(name)
    before = len(win.editor.features)
    d.accept()
    return len(win.editor.features) - before, d


def test_feature_and_new_sequence_dialogs(win, monkeypatch):
    """The feature dialog refuses an empty / unreadable / outside location and bad qualifiers with a message, and accepts a
    good one; the new-sequence dialog stays open on bad bases."""
    # --- feature dialog ---
    added, d = feature_attempt(win, loc="  ")
    assert added == 0 and "location" in win.messages[-1].lower()
    added, _ = feature_attempt(win, loc="banana")
    assert added == 0 and "Cannot read the location" in win.messages[-1] and "complement" in win.messages[-1]
    added, _ = feature_attempt(win, loc="5..99999")
    assert added == 0 and "outside" in win.messages[-1]
    added, _ = feature_attempt(win, qualifiers="=orphan value")
    assert added == 0 and "key=value" in win.messages[-1]
    added, _ = feature_attempt(win, qualifiers="two words=x")
    assert added == 0
    added, d = feature_attempt(win, loc=" 10..20 ", qualifiers="/note=\"fine\"\n\n  gene = g1 ", name="é ü")
    assert added == 1
    assert d.name.maxLength() >= 200
    assert dict(win.editor.features[-1].qualifiers)["note"] == "fine"
    # --- new sequence dialog ---
    seen = {}

    def fake_exec(dialog):
        dialog.findChild(QPlainTextEdit).setPlainText("ACGT hello world")
        dialog.findChild(QLineEdit).setText("n" * 500)
        dialog.findChild(QDialogButtonBox).button(QDialogButtonBox.Ok).click()   # OK with bad bases: stays open, says why
        seen["open"] = dialog.result() != QDialog.Accepted
        seen["problem"] = [lab.text() for lab in dialog.findChildren(QLabel) if "Not DNA" in lab.text()]
        seen["name_max"] = dialog.findChild(QLineEdit).maxLength()
        return QDialog.Rejected

    monkeypatch.setattr(QDialog, "exec", fake_exec)
    assert win.new_file() is None
    assert seen["open"] and seen["problem"] and seen["name_max"] <= 100


# ---- go-to / locus box, new sequence, search ---------------------------------------------------------------------------------

def locus_view(win):
    pane = win.workspace.panes(win.workspace.currentWidget())[0]
    pane.set_view("Alignment")
    return pane.build()


def test_locus_box(win):
    """Bad numbers in the locus box give a message instead of an exception; good input is accepted."""
    view = locus_view(win)
    for text in [",", "-", ":", "0", "99999999999999999999", "5000000", "–"]:
        win.errors.clear()
        view.locus.setText(text)
        view.go()                               # no exception; a message instead
        assert win.errors and win.errors[-1], text
    for text in ("12", "  100 - 200 "):
        win.errors.clear()
        view.locus.setText(text)
        view.go()
        assert not win.errors, text
    assert (win.editor.start, win.editor.end) == (99, 200)


def test_small_limits_and_links(win):
    """The search box is bounded; linked min / max spin boxes stay ordered; a tab name is one short line; the image export
    refuses a giant raster; a saved primer pair name is checked."""
    # --- search box ---
    assert win.search_bar.box.maxLength() == 5000
    win.search_bar.box.setText(">seq\nACGT" + "N" * 6000)
    win.search_bar.run(1)                         # no exception, no hang
    # --- linked min / max spin boxes ---
    from primer_design import DesignDialog
    win.editor.select(100, 200)
    d = DesignDialog(win, win.editor)
    d.min_len.setValue(30)
    assert d.max_len.value() >= 30
    d.max_len.setValue(15)
    assert d.min_len.value() <= 15

    from assembly_panel import AssemblyPanel
    p = AssemblyPanel(win)
    p.min_ov.setValue(50)
    assert p.max_ov.value() >= 50
    p.max_ov.setValue(20)
    assert p.min_ov.value() <= 20
    # --- tab name ---
    ws = win.workspace
    ws.set_tab_name(0, "a\nvery   long " + "x" * 300)
    assert "\n" not in ws.tabText(0) and len(ws.tabText(0)) <= 60
    # --- giant raster ---
    import image_export
    win.resize(1400, 850)
    app.processEvents()
    d = image_export.ImageExportDialog(win, win.workspace)
    d.dpi.setValue(2400)
    d.scale.setValue(1000)
    assert not d.ok.isEnabled() and "Too large" in d.warning.text()
    d.format.setCurrentText(next(k for k, v in image_export.FORMATS.items() if v == "svg"))
    assert d.ok.isEnabled()
    d.format.setCurrentIndex(0)
    d.dpi.setValue(300)
    assert d.ok.isEnabled()
    # --- pair name ---
    from primer_design import DesignDialog
    win.editor.select(200, 260)
    d = DesignDialog(win, win.editor)
    d.pairs = [{"fwd": {"seq": "ACGTACGTACGTACGTAC"}, "rev": {"seq": "TTGCAAGCTTGCAAGCTT"}}]   # one row to save
    d.table.setRowCount(1)
    d.table.setCurrentCell(0, 0)
    assert d.save("bad\tname") is None and win.messages and "control character" in win.messages[-1]
    assert d.save("x" * 100) is None and "too long" in win.messages[-1]


# ---- small spin box links and limits ----------------------------------------------------------------------------------------


