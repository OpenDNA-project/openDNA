"""The reference genome of the off-target checks: every record of a multi-record FASTA or GenBank file is searched (primer checks
and the guide finder), and reference genomes come from NCBI Datasets (network faked)."""
import io
import json
import os
import random
import sys
import tempfile
import zipfile
from pathlib import Path

import pytest

os.environ["OPENDNA_DATA_DIR"] = tempfile.mkdtemp(prefix="opendna_test_")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "app"))

from PySide6.QtCore import QSettings

from opendna_core import Document, Fasta

import cas9
import ncbi
import primer_panel
import reference
from editor import Editor, reverse_complement as rc
from reference import Reference

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
    primer_panel._genomes.clear()
    yield path, (chr1, chr2, chr3)
    primer_panel._genomes.clear()


def plasmid_editor(seq):
    return Editor(Document.new_sequence("plasmid", seq, True), Path("p.gb"))


def genome_check(ed, primer):
    return {n: (lv, t) for n, lv, t in primer_panel.primer_checks(ed, primer)}["Non-specific (genome)"]


def test_off_target_search_covers_every_record(genome, tmp_path, monkeypatch):
    """The off-target check of a primer sees copies on the second and third record and never compares the open record with
    itself; a search that ran out of time or was cancelled is not conclusive (yellow, never green); a 3' copy is yellow and a
    clean primer green. A multi-record GenBank reference is read one record at a time; the guide finder counts the sites of every
    record; a long record searched in overlapping windows neither loses nor doubles a site."""
    # --- copies on the second and third record
    path, (chr1, chr2, chr3) = genome
    ed = plasmid_editor(dna(300) + PRIMER + dna(300))
    level, text = genome_check(ed, PRIMER)
    assert level == "bad", "the whole primer is on chr2 and chr3, not in the first record"
    assert "chr2:" in text and "chr3:" in text and "chr1:" not in text and "501–524" in text
    ref = primer_panel.genome()
    assert ref.count == 3 and ref.label().startswith("3 sequences, ")

    # --- the open record is never compared with itself
    open_chr2 = Editor(Document.from_fasta(str(path), "chr2"), path)
    level, text = genome_check(open_chr2, PRIMER)
    assert level == "bad" and "chr3:" in text and "chr2:" not in text
    one = path.with_name("one.fa")
    write_fasta(one, [("chr2", chr2)])
    primer_panel.settings().setValue("primers/genome", str(one))
    assert genome_check(open_chr2, PRIMER) == ("info", "this sequence is the genome")
    primer_panel.settings().setValue("primers/genome", str(path))

    # --- out of time or cancelled
    primer_panel._genomes.clear()
    ref = primer_panel.genome()
    slow = ref.search(PRIMER, PRIMER, 0, 12, seconds=-1)   # no time at all
    assert slow.stopped == "time" and slow.searched == 0 and not slow.complete
    seen = []
    stopped = Reference(path).search(PRIMER, PRIMER, 0, 12, progress=lambda *a: seen.append(a), cancel=lambda: len(seen) >= 2)
    assert stopped.stopped == "cancelled" and stopped.searched < 3 and seen[0][2] == "chr1" and seen[0][4] == 3
    ed = plasmid_editor(dna(300) + PRIMER + dna(300))
    ref._cache.clear()
    ref._cache[(PRIMER, PRIMER, 0, 12, ("plasmid", ed.length))] = slow
    level, text = genome_check(ed, PRIMER)
    assert level == "warn" and "0 of 3" in text and "not conclusive" in text

    # --- a 3' copy is yellow, a clean primer green (rewrites the genome file)
    copy3 = "".join(other(c) for c in "TCA") + PRIMER[-14:]            # three other bases, then the 3' end of the primer
    write_fasta(path, [("chr1", chr1), ("chr2", dna(500)), ("chr3", dna(200) + copy3 + dna(300))])
    primer_panel._genomes.clear()
    level, text = genome_check(ed, PRIMER)
    assert level == "warn" and "chr3:" in text and "chr2:" not in text
    write_fasta(path, [("chr1", chr1), ("chr2", dna(500)), ("chr3", dna(500))])
    primer_panel._genomes.clear()
    level, text = genome_check(ed, PRIMER)
    assert level == "ok" and "3 sequences" in text

    # --- a multi-record GenBank reference
    a = Document.new_sequence("recA", dna(400), True)
    b = Document.new_sequence("recB", dna(300) + PRIMER + dna(300), False)
    a.save_genbank(tmp_path / "a.gb")
    b.save_genbank(tmp_path / "b.gb")
    both = tmp_path / "both.gb"
    both.write_text((tmp_path / "a.gb").read_text() + (tmp_path / "b.gb").read_text())
    ref = Reference(both)
    assert ref.count == 2 and [r[0] for r in ref.records] == ["recA", "recB"] and ref.records[0][2] is True
    r = ref.search(PRIMER, PRIMER, 0, 12)
    assert [h[0] for h in r.whole] == ["recB"] and r.searched == 2 and r.complete
    assert [n for n, _, _ in ref.sequences()] == ["recA", "recB"]

    # --- the guide finder counts the sites of every record
    guide = "GACTCGTACGTAGCATCGAT"
    sp = cas9.NUCLEASES["SpCas9 (NGG, 20 nt)"]
    filler = lambda n: "".join(rng.choice("AT") for _ in range(n))
    seq = filler(30) + guide + "TGG" + filler(30)
    (g,) = [x for x in cas9.candidates(seq, False, sp) if x.guide == guide]
    gpath = tmp_path / "g.fa"
    write_fasta(gpath, [("c1", filler(100)), ("c2", filler(40) + guide + "AGG" + filler(40)), ("c3", filler(40) + guide + "CGG" + filler(40))])
    (spec,) = cas9.scan(seq, False, "plasmid", [g], sp, 2, Reference(gpath))
    assert spec.genome[0] == 2 and spec.genome_name == "3 sequences", "one perfect copy on c2 and one on c3"
    assert cas9.verdict(spec, g)[0] == "bad"

    # --- overlapping windows
    seq = dna(1000)
    seq = seq[:395] + PRIMER + seq[419:]             # spans the window edge at 400
    seq = seq[:700] + rc(PRIMER) + seq[724:]
    wpath = tmp_path / "w.fa"
    write_fasta(wpath, [("chrW", seq)])
    whole = Reference(wpath).search(PRIMER, PRIMER, 0, 12)
    monkeypatch.setattr(reference, "WINDOW", 100)    # ten windows
    windowed = Reference(wpath).search(PRIMER, PRIMER, 0, 12)
    assert sorted(windowed.whole) == sorted(whole.whole) and len(windowed.whole) == 2
    assert sorted(windowed.hits) == sorted(whole.hits)
    assert [h[1] for h in sorted(windowed.whole)] == [395, 700]


REPORT = {"reports": [{"accession": "GCF_000000001.1", "organism": {"organism_name": "Madeup coli"},
                       "assembly_info": {"assembly_name": "ASM1v1", "assembly_level": "Complete Genome", "release_date": "2026-01-01",
                                         "refseq_category": "reference genome"},
                       "assembly_stats": {"total_sequence_length": "84", "number_of_component_sequences": 2}}], "total_count": 1}


def fake_zip():
    """An NCBI Datasets download: two chromosome files (the first without a final newline) and a GFF."""
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as z:
        z.writestr("README.md", "x")
        z.writestr("ncbi_dataset/data/GCF_000000001.1/chr1.fna", ">chr1 Madeup\n" + "ACGT" * 15)
        z.writestr("ncbi_dataset/data/GCF_000000001.1/chr2.fna", ">chr2\nTTTTGGGGCCCCAAAATTTTGGGG\n")
        z.writestr("ncbi_dataset/data/GCF_000000001.1/genomic.gff", "chr1\tRefSeq\tgene\t5\t20\t.\t+\t.\tName=abcA\n")
    return buf.getvalue()


def test_ncbi_genome_search_and_download(monkeypatch, tmp_path):
    """Searching NCBI Datasets by name or accession (an empty search is refused), and downloading a genome once: its chromosome
    files joined into one FASTA (with their newlines, every record readable), its GFF beside it, no partial file left, and the
    second request served from disk (network faked)."""
    # --- search ---
    urls = []
    monkeypatch.setattr(ncbi, "fetch", lambda url, timeout=60: urls.append(url) or json.dumps(REPORT))
    found = ncbi.search_genomes("Madeup coli")
    assert "taxon/Madeup%20coli/dataset_report" in urls[0] and "reference_only=true" in urls[0]
    a = found[0]
    assert (a.accession, a.length, a.sequences, a.level) == ("GCF_000000001.1", 84, 2, "Complete Genome")
    ncbi.search_genomes("gcf_000000001.1")
    assert "accession/GCF_000000001.1/dataset_report" in urls[1]
    with pytest.raises(ncbi.NcbiError):
        ncbi.search_genomes("  ")

    # --- download, once ---
    monkeypatch.setattr(ncbi, "genome_folder", lambda assembly: tmp_path / assembly.accession)
    calls, sizes = [], []
    monkeypatch.setattr(ncbi, "open_url", lambda url, headers, timeout=120: calls.append(url) or io.BytesIO(fake_zip()))
    fasta, gff, cached = ncbi.download_genome(a, progress=sizes.append)
    assert not cached and sizes and "abcA" in gff.read_text()
    assert "GENOME_FASTA" in calls[0] and "GENOME_GFF" in calls[0]
    assert fasta == tmp_path / "GCF_000000001.1" / "GCF_000000001.1.fna"
    assert Fasta(str(fasta)).sequences() == [("chr1", 60), ("chr2", 24)], "the members are joined with their newlines"
    assert not list(fasta.parent.glob("*.part"))
    assert ncbi.download_genome(a)[2] and len(calls) == 1, "kept: not downloaded twice"
