"""Folder Search: a sequence taken from a plasmid in demo/ is found in the right file, at the right place and strand; random DNA
finds nothing; a hit in a multi-record file names its record; an unreadable file is skipped and counted, not fatal."""
import os
import random
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ["OPENDNA_DATA_DIR"] = tempfile.mkdtemp(prefix="opendna_test_")
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "app"))

from PySide6.QtWidgets import QApplication, QWidget

from opendna_core import Document

from editor import reverse_complement as rc
from folder_search import FolderSearchDialog, hits_text, search_folder

app = QApplication.instance() or QApplication([])
DEMO = HERE.parent / "demo"


def test_a_piece_of_a_demo_plasmid_is_found_in_its_file_at_its_place_and_strand():
    text = Document(str(DEMO / "pFG051_Chlor.gb")).sequence()
    hits, skipped, files = search_folder(DEMO, text[1000:1060])
    top = [h for h in hits if h.path.name == "pFG051_Chlor.gb"][0]
    assert (top.start, top.end, top.strand, top.identity) == (1000, 1060, 1, 1.0)
    assert not skipped and files >= 5
    flipped, _, _ = search_folder(DEMO, rc(text[1000:1060]))
    top = [h for h in flipped if h.path.name == "pFG051_Chlor.gb"][0]
    assert (top.start, top.end, top.strand) == (1000, 1060, -1)


def test_random_dna_finds_nothing_in_the_demo_folder():
    random.seed(7)
    query = "".join(random.choice("ACGT") for _ in range(60))
    assert search_folder(DEMO, query)[0] == []


def test_a_hit_inside_a_multi_record_file_names_its_record_and_an_unreadable_file_is_skipped_and_counted(tmp_path):
    random.seed(3)
    a, b = ("".join(random.choice("ACGT") for _ in range(2000)) for _ in range(2))
    (tmp_path / "many.fa").write_text(f">first\n{a}\n>second\n{b}\n")
    (tmp_path / "broken.gb").write_text("this is not a GenBank file\n")
    hits, skipped, files = search_folder(tmp_path, b[500:560])
    assert (files, skipped) == (2, 1)
    assert [(h.record, h.start, h.end, h.strand) for h in hits] == [("second", 500, 560, 1)]
    assert "many.fa\tsecond\t501-560\tforward\t100.0%\t60" in hits_text(hits)


def test_cancelling_stops_the_search_before_the_next_file(tmp_path):
    (tmp_path / "a.fa").write_text(">a\nACGT\n")
    assert search_folder(tmp_path, "ACGTACGTACGT", cancelled=lambda: True)[0] == []


def test_the_dialog_lists_hits_remembers_the_folder_and_opens_the_hit(tmp_path, monkeypatch):
    random.seed(5)
    seq = "".join(random.choice("ACGT") for _ in range(1500))
    (tmp_path / "p.fa").write_text(f">p\n{seq}\n")
    opened = []
    win = QWidget()
    win.editor, win.open_path = None, opened.append
    dlg = FolderSearchDialog(win)
    dlg.query.setPlainText(seq[200:260])
    dlg.folder.setText(str(tmp_path))
    dlg.start()
    dlg.thread.join()
    dlg._poll()
    assert dlg.table.rowCount() == 1 and dlg.table.item(0, 2).text() == "201–260"
    dlg.copy_results()
    assert "p.fa" in QApplication.clipboard().text()
    dlg.open_hit(0)
    assert opened == [str(tmp_path / "p.fa")]
    assert FolderSearchDialog(win).folder.text() == str(tmp_path)
