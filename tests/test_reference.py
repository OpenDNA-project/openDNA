"""The reference genome of the off-target checks and the other records of a multi-record FASTA: every record is searched
(Primers and Spacer Finder), NCBI "use as reference" keeps them all, tracks follow the record, saving says what it leaves."""
import io
import os
import random
import sys
import tempfile
import time
import zipfile
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ["GENOMEEDIT_DATA_DIR"] = tempfile.mkdtemp(prefix="genomeedit_test_")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "app"))

from PySide6.QtCore import QSettings
from PySide6.QtTest import QTest
from PySide6.QtWidgets import QApplication, QMessageBox

from genomeedit_core import Document

import cas9
import ncbi
import primer_panel
import reference
from editor import Editor, reverse_complement as rc
from reference import Reference

app = QApplication.instance() or QApplication([])
rng = random.Random(5)
ACGT = "ACGT"


def dna(n):
    return "".join(rng.choice(ACGT) for _ in range(n))


def write_fasta(path, records, width=60):
    with open(path, "w", newline="") as fh:
        for name, seq in records:
            fh.write(f">{name} test\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + "\n")


def other(c):
    return ACGT[(ACGT.index(c) + 1) % 4]


PRIMER = dna(24)


@pytest.fixture
def genome(tmp_path, monkeypatch):
    """A 3-record FASTA: chr1 has nothing, chr2 the whole primer (top strand), chr3 its reverse complement; and the options
    (a throw-away settings file) that point the checks at it."""
    chr1 = dna(900)
    chr2 = dna(500) + PRIMER + dna(400)
    chr3 = dna(300) + rc(PRIMER) + dna(700)
    path = tmp_path / "genome.fa"
    write_fasta(path, [("chr1", chr1), ("chr2", chr2), ("chr3", chr3)])
    s = QSettings(str(tmp_path / "options.ini"), QSettings.IniFormat)
    monkeypatch.setattr(primer_panel, "settings", lambda: s)
    s.setValue("primers/genome_check", "true")
    s.setValue("primers/genome", str(path))
    s.setValue("primers/homology_seed", 12)
    s.setValue("primers/homology_mismatches", 0)
    return path, (chr1, chr2, chr3)


def plasmid_editor(seq):
    return Editor(Document.new_sequence("plasmid", seq, True), Path("p.gb"))


def genome_check(ed, primer):
    return {n: (lv, t) for n, lv, t in primer_panel.primer_checks(ed, primer)}["Non-specific (genome)"]


def test_genome_checks(genome, tmp_path):
    """The off-target check of a primer: it sees copies on the second and third record; never compares the open record with
    itself; says a search that ran out of time or was cancelled is not conclusive; runs in the primer dialog on a worker
    thread; and a 3' copy on a later record is yellow, a clean primer green."""
    # --- copies on the second and third record ---
    path, _ = genome
    ed = plasmid_editor(dna(300) + PRIMER + dna(300))
    level, text = genome_check(ed, PRIMER)
    assert level == "bad", "the whole primer is on chr2 and chr3, not in the first record"
    assert "chr2:" in text and "chr3:" in text and "chr1:" not in text and "501–524" in text
    ref = primer_panel.genome()
    assert ref.count == 3 and ref.label().startswith("3 sequences, ")
    # --- the open record is never compared with itself ---
    path, (chr1, chr2, chr3) = genome
    open_chr2 = Editor(Document.from_fasta(str(path), "chr2"), path)
    # the primer is on chr2 (open: left out) and chr3 (reverse strand): reported on chr3 only
    level, text = genome_check(open_chr2, PRIMER)
    assert level == "bad" and "chr3:" in text and "chr2:" not in text
    # a single-record file that is the open sequence
    one = path.with_name("one.fa")
    write_fasta(one, [("chr2", chr2)])
    primer_panel.settings().setValue("primers/genome", str(one))
    assert genome_check(open_chr2, PRIMER) == ("info", "this sequence is the genome")
    primer_panel.settings().setValue("primers/genome", str(path))   # back to the 3-record genome for what follows
    # --- out of time or cancelled ---
    path, _ = genome
    primer_panel._genomes.clear()   # (a merged test: start from a Reference with no cached searches)
    ref = primer_panel.genome()
    slow = ref.search(PRIMER, PRIMER, 0, 12, seconds=-1)   # no time at all
    assert slow.stopped == "time" and slow.searched == 0 and not slow.complete
    seen = []
    stopped = Reference(path).search(PRIMER, PRIMER, 0, 12, progress=lambda *a: seen.append(a), cancel=lambda: len(seen) >= 2)
    assert stopped.stopped == "cancelled" and stopped.searched < 3 and seen[0][2] == "chr1" and seen[0][4] == 3
    # in the check: a partial search with no hit is yellow, never the green "no match"
    ed = plasmid_editor(dna(300) + PRIMER + dna(300))
    ref._cache.clear()
    ref._cache[(PRIMER, PRIMER, 0, 12, ("plasmid", ed.length))] = slow
    level, text = genome_check(ed, PRIMER)
    assert level == "warn" and "0 of 3" in text and "not conclusive" in text
    primer_panel._genomes.clear()
    # --- the primer dialog searches on a worker thread ---
    ed = plasmid_editor(dna(300) + PRIMER + dna(300))
    d = primer_panel.PrimerDialog(None, ed)
    d.sequence.setText(PRIMER)
    d._update_report()
    t = time.perf_counter()
    while "chr3:" not in d.report.text() and time.perf_counter() - t < 10:
        QTest.qWait(30)
    assert "chr2:" in d.report.text() and "chr3:" in d.report.text()
    d.close()
    # --- 3' copy yellow, clean primer green (rewrites the genome file) ---
    path, (chr1, chr2, chr3) = genome
    copy3 = "".join(other(c) for c in "TCA") + PRIMER[-14:]            # three other bases, then the 3' end of the primer
    write_fasta(path, [("chr1", chr1), ("chr2", dna(500)), ("chr3", dna(200) + copy3 + dna(300))])
    primer_panel._genomes.clear()
    ed = plasmid_editor(dna(300) + PRIMER + dna(300))
    level, text = genome_check(ed, PRIMER)
    assert level == "warn" and "chr3:" in text and "chr2:" not in text
    write_fasta(path, [("chr1", chr1), ("chr2", dna(500)), ("chr3", dna(500))])
    primer_panel._genomes.clear()
    level, text = genome_check(ed, PRIMER)
    assert level == "ok" and "3 sequences" in text


def test_reference_searches_every_record(tmp_path, monkeypatch):
    """A multi-record GenBank reference is read one record at a time; the guide finder counts the sites of every record; a
    long record is searched in overlapping windows without losing or doubling a site."""
    # --- multi-record GenBank reference ---
    a = Document.new_sequence("recA", dna(400), True)
    b = Document.new_sequence("recB", dna(300) + PRIMER + dna(300), False)
    (tmp_path / "a.gb").parent.mkdir(exist_ok=True)
    a.save_genbank(tmp_path / "a.gb")
    b.save_genbank(tmp_path / "b.gb")
    both = tmp_path / "both.gb"
    both.write_text((tmp_path / "a.gb").read_text() + (tmp_path / "b.gb").read_text())
    with pytest.raises(OSError, match="reference genome"):
        Document(str(both))   # still not openable as one document, with a message that says what to do
    ref = Reference(both)
    assert ref.count == 2 and [r[0] for r in ref.records] == ["recA", "recB"] and ref.records[0][2] is True
    r = ref.search(PRIMER, PRIMER, 0, 12)
    assert [h[0] for h in r.whole] == ["recB"] and r.searched == 2 and r.complete
    assert [n for n, _, _ in ref.sequences()] == ["recA", "recB"]
    # --- guide finder ---
    guide = "GACTCGTACGTAGCATCGAT"
    sp = cas9.NUCLEASES["SpCas9 (NGG, 20 nt)"]
    filler = lambda n: "".join(rng.choice("AT") for _ in range(n))
    seq = filler(30) + guide + "TGG" + filler(30)
    (g,) = [x for x in cas9.candidates(seq, False, sp) if x.guide == guide]
    path = Path(tempfile.mkdtemp()) / "g.fa"
    write_fasta(path, [("c1", filler(100)), ("c2", filler(40) + guide + "AGG" + filler(40)), ("c3", filler(40) + guide + "CGG" + filler(40))])
    (spec,) = cas9.scan(seq, False, "plasmid", [g], sp, 2, Reference(path))
    assert spec.genome[0] == 2 and spec.genome_name == "3 sequences", "one perfect copy on c2 and one on c3"
    assert cas9.verdict(spec, g)[0] == "bad"
    # --- overlapping windows ---
    seq = dna(1000)
    seq = seq[:395] + PRIMER + seq[419:]             # spans the window edge at 400
    seq = seq[:700] + rc(PRIMER) + seq[724:]
    path = tmp_path / "w.fa"
    write_fasta(path, [("chrW", seq)])
    whole = Reference(path).search(PRIMER, PRIMER, 0, 12)
    monkeypatch.setattr(reference, "WINDOW", 100)    # ten windows
    windowed = Reference(path).search(PRIMER, PRIMER, 0, 12)
    assert sorted(windowed.whole) == sorted(whole.whole) and len(windowed.whole) == 2
    assert sorted(windowed.hits) == sorted(whole.hits)
    assert [h[1] for h in sorted(windowed.whole)] == [395, 700]


# ---- NCBI: use as reference keeps every sequence -----------------------------------------------------------------------------

def test_ncbi_members_are_joined_with_their_newlines_and_the_reference_is_the_fasta_itself(monkeypatch, tmp_path):
    import json
    from ncbi_dialog import GenomeDialog
    import main as M
    report = {"reports": [{"accession": "GCF_000000002.1", "organism": {"organism_name": "Madeup multi"},
                           "assembly_info": {"assembly_name": "ASM2", "assembly_level": "Chromosome", "release_date": "2026-01-01"},
                           "assembly_stats": {"total_sequence_length": "24", "number_of_component_sequences": 2}}], "total_count": 1}
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as z:   # the first member has no final newline
        z.writestr("ncbi_dataset/data/GCF_000000002.1/chr1.fna", ">chr1\nACGTACGTACGT")
        z.writestr("ncbi_dataset/data/GCF_000000002.1/chr2.fna", ">chr2\nTTTTGGGGCCCC\n")
    monkeypatch.setattr(ncbi, "fetch", lambda url, timeout=60: json.dumps(report))
    monkeypatch.setattr(ncbi, "open_url", lambda url, headers, timeout=120: io.BytesIO(buf.getvalue()))
    monkeypatch.setattr(ncbi, "genome_folder", lambda a: tmp_path / a.accession)
    monkeypatch.setattr(QMessageBox, "question", lambda *a, **k: QMessageBox.Yes)
    s = QSettings(str(tmp_path / "options.ini"), QSettings.IniFormat)
    import primer_db
    monkeypatch.setattr(primer_db, "settings", lambda: s)
    w = M.MainWindow()
    dialog = GenomeDialog(w)
    dialog.query.setText("Madeup multi")
    dialog.search()
    t = time.perf_counter()
    while not dialog.assemblies and time.perf_counter() - t < 10:
        QTest.qWait(20)
    assert dialog.reference_button.isEnabled(), "two sequences: usable as a reference too"
    dialog.use_as_reference()
    t = time.perf_counter()
    while "reference genome" not in dialog.status.text() and time.perf_counter() - t < 10:
        QTest.qWait(20)
    fasta = tmp_path / "GCF_000000002.1" / "GCF_000000002.1.fna"
    assert s.value("primers/genome") == str(fasta), "the FASTA itself, every record"
    assert "2 sequences" in dialog.status.text()
    assert [n for n, _ in __import__("genomeedit_core").Fasta(str(fasta)).sequences()] == ["chr1", "chr2"]
    w.close()


# ---- the other records elsewhere ---------------------------------------------------------------------------------------------

def test_other_records(tmp_path):
    """A contig name matches ignoring chr and case; tracks follow the record and a mismatching single contig is said; saving
    one record as GenBank says the others stay behind."""
    # --- contig names ---
    import main as M
    assert M.match_contig(["chr1", "chr2"], "chr2") == "chr2"
    assert M.match_contig(["1", "2"], "chr2") == "2" and M.match_contig(["Chr1", "Chr2"], "chr1") == "Chr1"
    assert M.match_contig(["1", "chr1"], "chrX") is None and M.match_contig(["chr1", "1"], "chr1") == "chr1"
    assert M.match_contig(["chr7"], "chr1") is None
    # --- tracks follow the record ---
    import main as M
    fa = tmp_path / "g.fa"
    write_fasta(fa, [("chr1", dna(300)), ("chr2", dna(300))])
    gff = tmp_path / "t.gff3"
    gff.write_text("##gff-version 3\nchr1\tx\tgene\t5\t50\t.\t+\t.\tName=a1\nchr2\tx\tgene\t10\t60\t.\t+\t.\tName=b2\n")
    lonely = tmp_path / "only7.gff3"
    lonely.write_text("##gff-version 3\nchr7\tx\tgene\t5\t50\t.\t+\t.\tName=z\n")
    from PySide6.QtWidgets import QMessageBox
    warned = []
    real_warning = QMessageBox.warning
    QMessageBox.warning = staticmethod(lambda parent, title, text, *a, **k: warned.append(text))
    try:
        w = M.MainWindow()
        assert w.load(str(fa))
        entry = w.load_track(str(gff))
        assert entry.seqid == "chr1" and not warned, "a matching track: no warning"
        assert w.load_track(str(lonely)).seqid == "chr7", "a track of another contig is still loaded (the user may know better)..."
        assert warned and "'chr7'" in warned[-1] and "'chr1'" in warned[-1], "... after a warning that names both"
    finally:
        QMessageBox.warning = real_warning
    w.switch_sequence("chr2")
    assert w.editor.fasta_record == "chr2"
    assert [(t.name, t.seqid) for t in w.editor.tracks] == [("t.gff3", "chr2")], "the track follows; the chr7 one has no chr2"
    assert "only7.gff3" in w.statusBar().currentMessage()
    w.editor.doc.modified = False
    w.close()
    # --- saving one record ---
    import main as M
    fa = tmp_path / "g.fa"
    write_fasta(fa, [("chr1", dna(300)), ("chr2", dna(300)), ("chr3", dna(300))])
    w = M.MainWindow()
    assert w.load(str(fa))
    assert M.other_records_note(w.editor, w.picker.count()).startswith(" (only chr1: the other 2 sequences of g.fa")
    assert w._write(tmp_path / "chr1.gb")
    assert "the other 2" in w.statusBar().currentMessage()
    w.editor.doc.modified = False
    w.close()


