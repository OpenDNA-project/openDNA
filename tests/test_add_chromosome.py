"""Import > DNA > Sequence File as New Chromosome: the sequences of another file become new records of the open file. A one-record
GenBank file becomes a multi-record file, a FASTA file gets more sequences; the records already in the file are copied unchanged,
and the open record keeps its unsaved edits."""
import os
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ["OPENDNA_DATA_DIR"] = tempfile.mkdtemp(prefix="opendna_test_")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "app"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import pytest
from PySide6.QtCore import Qt
from PySide6.QtWidgets import QApplication, QDialog, QFileDialog, QMessageBox

from opendna_core import Document, Fasta

import add_record
import multigb
from main import MainWindow
from sample_files import DNA_FILE, PLASMID, PLASMID2

app = QApplication.instance() or QApplication([])


def needs(*paths):
    if not all(p.exists() for p in paths):
        pytest.skip("demo files missing")


def copy(source, folder, name=None):
    target = folder / (name or source.name)
    target.write_bytes(source.read_bytes())
    return target


def accept_dialog(monkeypatch, untick=(), rename=None):
    """The Import as New Chromosome dialog answers OK at once, with rows `untick` unticked and {row: name} typed as names."""
    seen = {}

    def fake_exec(dialog):
        seen["dialog"] = dialog
        for row in untick:
            dialog.table.item(row, 0).setCheckState(Qt.Unchecked)
        for row, name in (rename or {}).items():
            dialog.table.item(row, 0).setText(name)
        problem = dialog.check()
        seen["problem"] = problem
        return QDialog.Rejected if problem else QDialog.Accepted
    monkeypatch.setattr(add_record.AddRecordsDialog, "exec", fake_exec)
    return seen


def opened(path):
    win = MainWindow()
    win.show()
    assert win.load(str(path))
    return win


def test_one_record_genbank_becomes_multi_record(tmp_path, monkeypatch):
    """A plasmid file gets a second record: the first is copied byte for byte, the window shows the new one, the map lists both,
    and the first comes back with its unsaved edit, which Save writes in place."""
    needs(PLASMID, PLASMID2)
    path = copy(PLASMID, tmp_path)
    original = path.read_bytes()
    win = opened(path)
    first = win.editor.doc.name
    before = win.editor.doc.sequence()
    win.editor.set_cursor(10)
    win.editor.type_bases("ACGT")   # an unsaved edit: kept, not written by the import
    accept_dialog(monkeypatch)
    assert win.import_as_chromosome(str(PLASMID2))
    added = Document(str(PLASMID2)).name
    assert multigb.is_multi_genbank(path)
    assert [n for n, *_ in multigb.records(path)] == [first, added]
    assert path.read_bytes().startswith(original), "the record already in the file is unchanged"
    other = multigb.read_record(path, added)
    plasmid2 = Document(str(PLASMID2))
    assert other.sequence() == plasmid2.sequence() and other.circular == plasmid2.circular
    assert len(other.features()) >= len([f for f in plasmid2.features() if f.kind != "primer_bind"])
    ed = win.editor
    assert ed.multi_genbank and ed.fasta_record == added and {n for n, _ in ed.fasta_records} == {first, added}
    assert win.picker.count() == 2
    assert win.record_docs[first].modified, "the edited record stays in memory with its edit"
    ed.record_requested.emit(first)
    app.processEvents()
    assert win.editor.fasta_record == first and win.editor.doc.sequence() == before[:10] + "ACGT" + before[10:]
    assert win.save()
    assert multigb.read_record(path, first).sequence() == before[:10] + "ACGT" + before[10:]
    assert multigb.read_record(path, added).sequence() == plasmid2.sequence(), "the added record is untouched by Save"
    assert not list(tmp_path.glob("*.saving"))


def test_names_are_unique_and_can_be_changed(tmp_path, monkeypatch):
    """The same plasmid added to itself gets "_2"; a name typed in the dialog is used (spaces made "_"); a name the file already
    has is refused, and nothing is written."""
    needs(PLASMID)
    path = copy(PLASMID, tmp_path)
    win = opened(path)
    name = win.editor.doc.name
    accept_dialog(monkeypatch)
    assert win.import_as_chromosome(str(PLASMID))
    assert [n for n, *_ in multigb.records(path)] == [name, f"{name}_2"]
    accept_dialog(monkeypatch, rename={0: "my copy"})
    assert win.import_as_chromosome(str(PLASMID))
    assert [n for n, *_ in multigb.records(path)] == [name, f"{name}_2", "my_copy"]
    size = path.stat().st_size
    seen = accept_dialog(monkeypatch, rename={0: name})
    assert not win.import_as_chromosome(str(PLASMID))
    assert "already has" in seen["problem"] and path.stat().st_size == size


def test_multi_fasta_source_and_target(tmp_path, monkeypatch):
    """Sequences of a multi-FASTA file are added to a FASTA file (only those ticked); then a GenBank plasmid is added as bases
    only. The .fai index is rebuilt with the new sequences, and the line endings of a Windows file are kept."""
    needs(PLASMID)
    target = tmp_path / "genome.fa"
    target.write_bytes(b">chrA first\r\n" + b"ACGT" * 30 + b"\r\n>chrB\r\n" + b"GGCC" * 10)   # no line ending at the end
    Fasta(str(target))   # writes genome.fa.fai, which must not be used afterwards
    source = tmp_path / "more.fasta"
    source.write_text(">p1 plasmid one\n" + "A" * 70 + "\n>p2\n" + "C" * 5 + "\n>chrA\n" + "T" * 8 + "\n")
    win = opened(target)
    assert win.editor.fasta_record and not win.editor.multi_genbank
    accept_dialog(monkeypatch, untick=(1,))   # p1 and chrA (renamed chrA_2: the file has a chrA), not p2
    assert win.import_as_chromosome(str(source))
    names = dict(Fasta(str(target)).sequences())
    assert names == {"chrA": 120, "chrB": 40, "p1": 70, "chrA_2": 8}
    data = target.read_bytes()
    assert b"\n" not in data.replace(b"\r\n", b""), "only Windows line endings"
    assert b">p1\r\n" + b"A" * 60 + b"\r\n" + b"A" * 10 + b"\r\n" in data   # 60 bases a line (the core's FASTA reader keeps no description)
    assert Document.from_fasta(str(target), "chrB").sequence() == "GGCC" * 10, "the old last line was ended, not joined"
    assert win.editor.fasta_record == "p1" and not win.editor.multi_genbank and win.picker.count() == 4
    accept_dialog(monkeypatch)
    assert win.import_as_chromosome(str(PLASMID))
    plasmid = Document(str(PLASMID))
    assert Document.from_fasta(str(target), plasmid.name).sequence() == plasmid.sequence()
    assert win.picker.count() == 5


def test_dna_source_keeps_primers_and_dna_target_is_saved_first(tmp_path, monkeypatch):
    """A .dna file's primers are written as primer_bind features of the new record (as Save As GenBank does). A .dna file open in
    the window is never changed: it is saved as GenBank first, and the chromosome added to that file."""
    needs(PLASMID, DNA_FILE)
    path = copy(PLASMID, tmp_path)
    win = opened(path)
    accept_dialog(monkeypatch)
    assert win.import_as_chromosome(str(DNA_FILE))
    added = multigb.read_record(path, Document(str(DNA_FILE)).name)
    assert any(f.kind == "primer_bind" for f in added.features())
    # a .dna file open: Save As GenBank first
    dna = copy(DNA_FILE, tmp_path)
    original = dna.read_bytes()
    win = opened(dna)
    saved = tmp_path / "from_dna.gb"
    monkeypatch.setattr(QMessageBox, "question", lambda *a, **k: QMessageBox.Save)
    monkeypatch.setattr(QFileDialog, "getSaveFileName", lambda *a, **k: (str(saved), ""))
    accept_dialog(monkeypatch)
    assert win.import_as_chromosome(str(PLASMID))
    assert dna.read_bytes() == original
    assert [n for n, *_ in multigb.records(saved)] == [Document(str(DNA_FILE)).name, Document(str(PLASMID)).name]
    assert win.editor.path == saved and win.editor.multi_genbank


def test_cancelled_dialog_writes_nothing(tmp_path, monkeypatch):
    needs(PLASMID, PLASMID2)
    path = copy(PLASMID, tmp_path)
    original = path.read_bytes()
    win = opened(path)
    monkeypatch.setattr(add_record.AddRecordsDialog, "exec", lambda dialog: QDialog.Rejected)
    assert not win.import_as_chromosome(str(PLASMID2))
    assert path.read_bytes() == original and not win.editor.multi_genbank and win.editor.fasta_record is None


def test_record_name():
    assert multigb.record_name(" chr 1 ") == "chr_1"
    assert multigb.record_name("", ["sequence"]) == "sequence_2"
    assert multigb.record_name("a", ["a", "a_2"]) == "a_3"
