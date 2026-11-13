"""Files of several records (a GenBank file of several LOCUS records, a multi-FASTA): one record is read and saved back among the
others unchanged; records are added, reordered and removed with every other record copied byte for byte."""
import os
import sys
import tempfile
from pathlib import Path

import pytest

os.environ["OPENDNA_DATA_DIR"] = tempfile.mkdtemp(prefix="opendna_test_")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "app"))

from opendna_core import Document, Fasta

import add_record
import fileio
import multigb
from editor import Editor
from main import match_contig, other_records_note
from sample_files import DNA_FILE, PLASMID, PLASMID2, VIBRIO


def three_records(folder):
    """A GenBank file of three records: the two demo plasmids and a copy of the first named "copy"."""
    if not (PLASMID.exists() and PLASMID2.exists()):
        pytest.skip("demo plasmids missing")
    first = PLASMID.read_text(encoding="utf-8").rstrip("\n") + "\n"
    copy = first.replace(Document(str(PLASMID)).name, "copy", 1)
    path = folder / "three.gb"
    path.write_bytes((first + PLASMID2.read_text(encoding="utf-8").rstrip("\n") + "\n" + copy).encode("utf-8"))   # "\n" endings
    return path


def chunks(path):
    """{name: the bytes of the record} of a GenBank file."""
    data = path.read_bytes()
    _, spans = multigb.record_spans(path)
    return {n: data[s:e] for n, s, e in spans}


def names(path):
    return [n for n, *_ in multigb.records(path)]


def test_one_record_is_read_and_saved_in_place(tmp_path, monkeypatch):
    """A GenBank file of several records is recognised and listed in file order; a record reads as its own document; saving one
    record (with an edit) changes only that record; a save whose final rename fails leaves the file exactly as it was; a track's
    contig name matches ignoring chr and case; saving one record of a FASTA as GenBank says the others stay behind; the real
    two-chromosome Vibrio file lists both."""
    path = three_records(tmp_path)
    a, b, c = names(path)
    # --- recognised and listed
    assert multigb.is_multi_genbank(path) and not multigb.is_multi_genbank(PLASMID)
    assert [(n, length) for n, length, _ in multigb.records(path)] == [
        (a, Document(str(PLASMID)).length), (b, Document(str(PLASMID2)).length), ("copy", Document(str(PLASMID)).length)]
    assert multigb.read_record(path, b).sequence() == Document(str(PLASMID2)).sequence()
    with pytest.raises(OSError, match="no record named"):
        multigb.read_record(path, "nope")

    # --- saving one record changes only that one
    before = chunks(path)
    doc = multigb.read_record(path, a)
    seq = doc.sequence()
    doc.insert(10, "ACGTACGT")
    multigb.save_document(path, a, doc)
    assert multigb.read_record(path, a).sequence() == seq[:10] + "ACGTACGT" + seq[10:]
    after = chunks(path)
    assert after[b] == before[b] and after[c] == before[c], "the other records are untouched"
    assert names(path) == [a, b, c] and not list(tmp_path.glob("*.saving"))

    # --- a save that cannot replace the file leaves it as it was
    saved = path.read_bytes()

    def locked(tmp, target, *args, **kwargs):
        raise PermissionError(13, "The process cannot access the file because it is being used by another process")
    with monkeypatch.context() as m:
        m.setattr(fileio, "replace", locked)
        doc.insert(0, "TTTT")
        with pytest.raises(PermissionError):
            multigb.save_document(path, a, doc)
    assert path.read_bytes() == saved and not list(tmp_path.glob("*.saving")), "nothing half written, no stray file"

    # --- contig names of tracks; the note when one record of a FASTA is saved alone
    assert match_contig(["chr1", "chr2"], "chr2") == "chr2"
    assert match_contig(["1", "2"], "chr2") == "2" and match_contig(["Chr1", "Chr2"], "chr1") == "Chr1"
    assert match_contig(["1", "chr1"], "chrX") is None and match_contig(["chr1", "1"], "chr1") == "chr1"
    assert match_contig(["chr7"], "chr1") is None
    fa = tmp_path / "g.fa"
    fa.write_text(">chr1\nACGT\n>chr2\nGGCC\n>chr3\nTTAA\n")
    ed = Editor(Document.from_fasta(str(fa), "chr1"), fa)
    ed.fasta_record = "chr1"
    assert other_records_note(ed, 3).startswith(" (only chr1: the other 2 sequences of g.fa")
    assert other_records_note(ed, 1) == ""

    # --- the real two-chromosome file
    if VIBRIO.exists():
        assert multigb.is_multi_genbank(VIBRIO)
        assert [(n, length) for n, length, _ in multigb.records(VIBRIO)] == [("NC_002505", 2961149), ("NC_002506", 1072315)]


def test_records_are_added_reordered_and_removed_byte_for_byte(tmp_path):
    """Records are reordered and removed with each copied byte for byte (GenBank, and FASTA with a byte order mark, CRLF, no final
    line ending and a stale index), and bad orders are refused; a new record gets a free name; a plasmid added to a one-record
    GenBank file keeps the first byte for byte (a .dna source brings its primers as primer_bind features); FASTA sequences are
    added with the file's own line endings, the old last line ended and the index rebuilt."""
    # --- reorder and remove (GenBank)
    path = three_records(tmp_path)
    before = chunks(path)
    a, b, c = names(path)
    multigb.rewrite_records(path, [c, a, b])
    assert names(path) == [c, a, b] and chunks(path) == before
    multigb.rewrite_records(path, [c, b])
    assert names(path) == [c, b] and chunks(path) == {c: before[c], b: before[b]}
    for bad in ([], [c, c], [c, "nope"]):
        with pytest.raises(OSError):
            multigb.rewrite_records(path, bad)
    assert names(path) == [c, b] and not list(tmp_path.glob("*.saving"))

    # --- reorder and remove (FASTA)
    fa = tmp_path / "g.fa"
    fa.write_bytes(b"\xef\xbb\xbf>chr1 one\r\nACGT\r\nAC\r\n>chr2\r\nGGGG\r\n>chr3\r\nTT")
    Fasta(str(fa))
    multigb.rewrite_records(fa, ["chr3", "chr1", "chr2"])
    assert fa.read_bytes() == b"\xef\xbb\xbf>chr3\r\nTT\r\n>chr1 one\r\nACGT\r\nAC\r\n>chr2\r\nGGGG\r\n"
    assert [n for n, _ in Fasta(str(fa)).sequences()] == ["chr3", "chr1", "chr2"]
    assert Document.from_fasta(str(fa), "chr1").sequence() == "ACGTAC"
    multigb.rewrite_records(fa, ["chr2"])
    assert dict(Fasta(str(fa)).sequences()) == {"chr2": 4}

    # --- free record names
    assert multigb.record_name(" chr 1 ") == "chr_1"
    assert multigb.record_name("", ["sequence"]) == "sequence_2"
    assert multigb.record_name("a", ["a", "a_2"]) == "a_3"

    # --- a GenBank file gets a second record (the same plasmid: "_2"), and a .dna one with its primers
    one = tmp_path / "one.gb"
    one.write_bytes(PLASMID.read_bytes())
    original = one.read_bytes()
    first = Document(str(one)).name
    (src,) = add_record.source_records(PLASMID)
    name = multigb.record_name(src.name, add_record.current_names(Editor(Document(str(one)), one)))
    assert name == f"{first}_2"
    multigb.append_genbank(one, [multigb.renamed_genbank(src.load(), name)])
    assert names(one) == [first, name] and one.read_bytes().startswith(original), "the record already there is unchanged"
    assert multigb.read_record(one, name).sequence() == Document(str(PLASMID)).sequence()
    assert Document(str(PLASMID)).name == first, "renaming the copy leaves the source document alone"
    if DNA_FILE.exists():
        dna = Document(str(DNA_FILE))
        primers = [("primer_bind", f"{s + 1}..{s + len(seq)}", [("label", n), ("note", f"sequence: {seq}")])
                   for n, seq in dna.primers() for s in [dna.sequence().upper().find(seq.upper())] if s >= 0]
        assert primers
        multigb.append_genbank(one, [multigb.renamed_genbank(dna, "from_dna", primers)])
        assert names(one) == [first, name, "from_dna"]
        assert any(f.kind == "primer_bind" for f in multigb.read_record(one, "from_dna").features())

    # --- FASTA sequences added: the file's line endings, the old last line ended, the index rebuilt
    target = tmp_path / "genome.fa"
    target.write_bytes(b">chrA first\r\n" + b"ACGT" * 30 + b"\r\n>chrB\r\n" + b"GGCC" * 10)   # no line ending at the end
    Fasta(str(target))   # writes genome.fa.fai, which must not be used afterwards
    source = tmp_path / "more.fasta"
    source.write_text(">p1 plasmid one\n" + "A" * 70 + "\n>p2\n" + "C" * 5 + "\n>chrA\n" + "T" * 8 + "\n")
    assert [(r.name, r.length) for r in add_record.source_records(source)] == [("p1", 70), ("p2", 5), ("chrA", 8)]
    taken = [n for n, _ in Fasta(str(target)).sequences()]
    picked = [r for r in add_record.source_records(source) if r.name != "p2"]   # p2 unticked
    added = []
    for r in picked:
        added.append((multigb.record_name(r.name, taken + [n for n, _ in added]), r.load()))
    multigb.append_fasta(target, added)
    assert dict(Fasta(str(target)).sequences()) == {"chrA": 120, "chrB": 40, "p1": 70, "chrA_2": 8}
    data = target.read_bytes()
    assert b"\n" not in data.replace(b"\r\n", b""), "only Windows line endings"
    assert b">p1\r\n" + b"A" * 60 + b"\r\n" + b"A" * 10 + b"\r\n" in data   # 60 bases a line
    assert Document.from_fasta(str(target), "chrB").sequence() == "GGCC" * 10, "the old last line was ended, not joined"
