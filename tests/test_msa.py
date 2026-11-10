"""Multiple sequence alignment: the built-in aligner on a planted set (substitutions and one indel), the consensus and identity
figures, the aligned-FASTA reader, and the Multiple Alignment view opened on a multi-record file (a column click moves the cursor)."""
import os
import random
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ["OPENDNA_DATA_DIR"] = tempfile.mkdtemp(prefix="opendna_test_")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "app"))

import pytest
from PySide6.QtCore import Qt
from PySide6.QtTest import QTest
from PySide6.QtWidgets import QApplication

import msa
from main import MainWindow
from msa_view import MsaView

app = QApplication.instance() or QApplication([])
DEMO = Path(__file__).resolve().parent.parent / "demo" / "msa_demo.fasta"
random.seed(3)
BASE = "".join(random.choice("ACGT") for _ in range(200))


def planted():
    """Four copies of BASE: one identical, one with the base at 50 changed, one with the base at 120 changed, one missing 3 bases at 80."""
    swap = lambda s, i: s[:i] + ("A" if s[i] != "A" else "C") + s[i + 1:]
    return ["ref", "sub50", "sub120", "del80"], [BASE, swap(BASE, 50), swap(BASE, 120), BASE[:80] + BASE[83:]]


def test_planted_substitutions_and_one_indel_are_aligned_with_the_right_consensus_and_identity():
    """Rows keep their sequences, the deletion becomes a 3-column gap, the consensus is the original, identity counts gaps as differences."""
    names, seqs = planted()
    a = msa.align(names, seqs)
    assert [r.replace("-", "") for r in a.rows] == seqs
    assert a.width == 200 and a.rows[3][80:83] == "---" and a.rows[0] == BASE
    assert a.consensus() == BASE            # every variant is outvoted 3 to 1 or 2 to 2 (gaps count only when most rows have one)
    assert a.identity(0) == 100.0
    assert a.identity(1) == pytest.approx(100 * 199 / 200)
    assert a.identity(3) == pytest.approx(100 * 197 / 200)   # 3 gap columns where the reference has bases
    a.reference = None                                         # against the consensus, which here is the reference
    assert a.identities() == [a.identity(i, 0) for i in range(4)]
    cons = a.conservation()
    assert cons[0] == 1.0 and cons[50] == 0.75 and cons[81] == 0.75   # the substituted base, and a base the deletion lacks


def test_columns_and_positions_map_both_ways_through_gaps():
    a = msa.Alignment(["x", "y"], ["AC--GT", "ACTTGT"])
    assert [a.position(0, c) for c in range(6)] == [0, 1, 2, 2, 2, 3]
    assert a.column_of(0, 2) == 4 and a.column_of(1, 2) == 2


def test_aligned_fasta_text_from_mafft_or_clustal_is_read_by_name():
    names, rows = msa.parse_aligned_fasta(">a\nAC-\nGT\n>b desc\nACTG\nT-\n")
    assert names == ["a", "b desc"] and rows == ["AC-GT", "ACTGT-"]


def test_the_aligner_refuses_one_sequence_and_sequences_too_long_for_it():
    with pytest.raises(ValueError):
        msa.align(["a"], ["ACGT"])
    with pytest.raises(ValueError):
        msa.align(["a", "b"], ["A" * (msa.MAX_BUILTIN + 1), "A"])


def test_the_view_opens_on_a_multi_record_file_and_a_column_click_moves_the_cursor():
    """The demo FASTA is aligned when the view opens; clicking a column puts the open record's cursor at that column's base."""
    win = MainWindow()
    win.show()
    assert win.load(str(DEMO))
    view = MsaView(win.editor)
    view.show()
    a = view.alignment
    assert a is not None and len(a.rows) == 5 and a.names[0] == "reference"
    assert [r.replace("-", "") for r in a.rows][0] == win.editor.doc.sequence()
    assert a.identity(0) == 100.0 and all(i > 95 for i in a.identities())
    open_row = view.open_row()
    col = 300   # right of the 6-base gap that clone_3 lacks: the positions differ between rows there
    QTest.mouseClick(view.canvas.viewport(), Qt.LeftButton, pos=view.canvas.viewport().rect().center())  # (any click must not crash)
    view._column_clicked(col)
    assert win.editor.start == a.position(open_row, col)
    assert f"Column {col + 1}" in view.status.text() and "clone_3" in view.status.text()
    win.close()
