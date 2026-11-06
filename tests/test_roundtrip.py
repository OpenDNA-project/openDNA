"""File I/O round trips: GenBank / .dna / FASTA files written by us read back the same (by us and by Biopython),
odd content survives, short and corrupt files are handled, saving and autosave never damage the original.

Real files in the repo root are used when present (skipped otherwise); everything is written in tmp_path."""
import os
import random
import sys
import tempfile
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ["OPENDNA_DATA_DIR"] = tempfile.mkdtemp(prefix="opendna_test_")
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "app"))

from Bio import SeqIO
from PySide6.QtWidgets import QApplication, QMessageBox

from opendna_core import Document, Fasta

import main as M
import primer_db
from editor import Editor, primer_sequence

app = QApplication.instance() or QApplication([])

# Biopython warns about the width of the LOCUS line we write for short names; the parse itself is right.
pytestmark = pytest.mark.filterwarnings("ignore::Bio.BiopythonParserWarning")

from sample_files import GENOME as GENOME_DNA, PLASMID as PLASMID_GB, DNA_FILE as PLASMID_DNA   # copies (see sample_files.py)
DEMOS = sorted((ROOT / "demo").glob("demo_*.gb"))   # the made-up demo files, not the example plasmids and genome


@pytest.fixture(autouse=True)
def private_primer_db(tmp_path, monkeypatch):
    monkeypatch.setattr(primer_db, "_db", primer_db.PrimerDB(tmp_path / "primers.tsv"))


# ---- helpers ---------------------------------------------------------------------------------------------------------------

def snap(doc):
    """Everything a save must keep: name, definition, topology, bases (as typed) and every feature as written."""
    return {
        "name": doc.name, "definition": doc.definition, "circular": doc.circular, "bases": doc.sequence(),
        "features": [(f.kind, f.location, f.qualifiers) for f in doc.features()],
    }


def bio_read(path):
    with open(path, encoding="utf-8") as fh:          # our files are UTF-8; Biopython would use the console's code page
        return SeqIO.read(fh, "genbank")


def our_parts(f):
    return sorted((s, e, f.strand) for s, e in f.spans)


def bio_parts(f):
    return sorted((int(p.start), int(p.end), p.strand) for p in f.location.parts)


def bio_qualifiers(f):
    return {k: list(v) for k, v in f.qualifiers.items()}


def our_qualifiers(f):
    out = {}
    for k, v in f.qualifiers:
        out.setdefault(k, []).append("" if v is None else v)
    return out


def assert_biopython_agrees(path, doc):
    """Our reading of the file in memory = Biopython's reading of the file on disk."""
    rec = bio_read(path)
    assert str(rec.seq) == doc.sequence().upper()                    # Biopython upper-cases; we keep the case as typed
    assert rec.annotations["topology"] == ("circular" if doc.circular else "linear")
    assert rec.id == doc.name or rec.name == doc.name
    feats = doc.features()
    assert len(rec.features) == len(feats)
    for ours, theirs in zip(feats, rec.features):
        assert theirs.type == ours.kind
        assert bio_parts(theirs) == our_parts(ours), ours.location
        assert bio_qualifiers(theirs) == our_qualifiers(ours), ours.location


def rich_document(circular=True, seed=3):
    """A made-up record with joins, complements, fuzzy ends, order() and awkward qualifiers."""
    rng = random.Random(seed)
    doc = Document.new_sequence("rich", "".join(rng.choice("ACGT") for _ in range(400)), circular)
    doc.add_features([
        ("source", "1..400", [("organism", "synthetic construct"), ("mol_type", "other DNA")]),
        ("CDS", "join(10..20,30..40)", [("label", "joined"), ("codon_start", "1")]),
        ("misc_feature", "complement(join(60..70,80..90))", [("label", "complement of a join")]),
        ("misc_feature", "<100..>120", [("label", "fuzzy both ends")]),
        ("misc_feature", "<150..180", [("label", "fuzzy start")]),
        ("misc_feature", "200..>230", [("label", "fuzzy end")]),
        ("misc_feature", "order(250..260,300..310)", [("label", "ordered")]),
        ("gene", "complement(320..390)", [("gene", "abcD"), ("locus_tag", "T_0001"), ("pseudo", None)]),
    ])
    return doc


def save_and_reload(doc, path):
    doc.save_genbank(path)
    return Document(path)


def genbank_sources(extra=()):
    """(id, path) of every real or demo file that is present (plus the `extra` ones when present)."""
    out = [("plasmid_dna", PLASMID_DNA), ("plasmid_gb", PLASMID_GB)] + [(p.stem, p) for p in DEMOS] + list(extra)
    return [(i, p) for i, p in out if p.exists()]


# ---- .dna / GenBank -> GenBank -> reload -------------------------------------------------------------------------------

def test_real_file_round_trips(tmp_path):
    """Real and demo files: load, save as GenBank and reload is identical (bases, feature locations, qualifiers with their order,
    topology, name, definition); saving a reloaded file again is byte-identical (twice over); and Biopython reads our file the same."""
    sources = genbank_sources(extra=[("genome_dna", GENOME_DNA)])
    if not sources:
        pytest.skip("no real or demo files")
    # --- load, save as GenBank, reload is identical
    for pid, path in sources:
        doc = Document(path)
        before = snap(doc)
        assert before["features"], f"{pid}: the file should have annotations to make this meaningful"
        if pid == "genome_dna":   # the 4.6 Mb E. coli genome, the largest file
            assert doc.length > 4_500_000 and doc.circular and len(doc.features()) > 4000
        out = tmp_path / "out.gb"
        doc.save_genbank(out)
        back = Document(out)
        after = snap(back)
        for key in before:
            assert after[key] == before[key], f"{pid}: {key}"
        assert [f.spans for f in back.features()] == [f.spans for f in doc.features()], pid
        assert [f.strand for f in back.features()] == [f.strand for f in doc.features()], pid
        assert not back.modified, pid

    # --- save, reload, save is byte-identical (no drift between saves), the genome too
    for pid, path in genbank_sources(extra=[("genome_dna", GENOME_DNA)]):
        first, second, third = tmp_path / "1.gb", tmp_path / "2.gb", tmp_path / "3.gb"
        Document(path).save_genbank(first)
        Document(first).save_genbank(second)
        Document(second).save_genbank(third)
        assert first.read_bytes() == second.read_bytes() == third.read_bytes(), pid

    # --- a file saved by us is read the same by Biopython
    for pid, path in genbank_sources():
        doc = Document(path)
        out = tmp_path / "out.gb"
        doc.save_genbank(out)
        assert_biopython_agrees(out, doc)


def test_rich_document_round_trip(tmp_path):
    """A made-up record with joins, complements, fuzzy ends and awkward qualifiers: save, reload, save is byte-identical, and
    Biopython reads it the same, circular and linear."""
    # --- save, reload, save is byte-identical
    first, second = tmp_path / "1.gb", tmp_path / "2.gb"
    doc = rich_document()
    back = save_and_reload(doc, first)
    assert snap(back) == snap(doc)
    back.save_genbank(second)
    assert first.read_bytes() == second.read_bytes()

    # --- read the same by Biopython, circular and linear
    for circular in (True, False):
        doc = rich_document(circular)
        out = tmp_path / "rich.gb"
        doc.save_genbank(out)
        assert_biopython_agrees(out, doc)
        rec = bio_read(out)
        by_label = {f.qualifiers.get("label", [""])[0]: f for f in rec.features}
        assert by_label["fuzzy both ends"].location.start.__class__.__name__ == "BeforePosition", circular
        assert by_label["fuzzy both ends"].location.end.__class__.__name__ == "AfterPosition", circular
        assert by_label["complement of a join"].location.strand == -1, circular


# ---- topology --------------------------------------------------------------------------------------------------------------

def test_topology_survives_genbank_save_and_load(tmp_path):
    for circular in (True, False):
        doc = Document.new_sequence("topo", "ACGTACGTTTGACCA" * 7, circular)
        out = tmp_path / "t.gb"
        doc.save_genbank(out)
        assert (" circular " in out.read_text()) == circular and (" linear " in out.read_text()) != circular
        back = Document(out)
        assert back.circular is circular and back.sequence() == doc.sequence()
        # Toggling it after loading is saved too (and is one undo step that gives the old topology back).
        back.set_topology(not circular)
        back.save_genbank(out)
        assert Document(out).circular is (not circular)
        assert back.undo() and back.circular is circular


# ---- locations -------------------------------------------------------------------------------------------------------------

LOCATIONS = [
    "10..20", "complement(10..20)", "join(10..20,30..40)", "complement(join(10..20,30..40))",
    "join(complement(30..40),complement(10..20))", "join(10..20,complement(30..40))",
    "<10..20", "10..>20", "<10..>20", "order(10..20,30..40)", "join(10..20,30..40,50..60)",
    "join(350..400,1..20)",                                           # across the origin of the circle
]


def test_location_forms_survive_a_genbank_round_trip(tmp_path):
    """Every location form survives, written as given and read the same by Biopython; the gap of a between location has the
    same span as Biopython's."""
    # --- every location form
    for location in LOCATIONS:
        doc = Document.new_sequence("loc", "ACGT" * 100, True)
        doc.add_features([("misc_feature", location, [("label", "x")])])
        (before,) = doc.features()
        out = tmp_path / "loc.gb"
        doc.save_genbank(out)
        (after,) = Document(out).features()
        assert (after.location, after.spans, after.strand) == (before.location, before.spans, before.strand), location
        assert after.location == location, location                       # written in the form it was given
        rec = bio_read(out)                                               # (a part's own strand is Biopython's only: a join can mix them)
        assert sorted((a, b) for a, b, _ in bio_parts(rec.features[0])) == sorted(after.spans), location

    # --- `20^21` is the gap between bases 20 and 21 (1-based): Biopython puts it at 0-based position 20
    doc = Document.new_sequence("gap", "ACGT" * 20, False)
    doc.add_features([("misc_feature", "20^21", [("label", "gap")])])
    out = tmp_path / "gap.gb"
    doc.save_genbank(out)
    (f,) = Document(out).features()
    assert f.location == "20^21"                                      # the text itself survives
    rec = bio_read(out)
    assert (int(rec.features[0].location.start), int(rec.features[0].location.end)) == (20, 20)
    assert f.spans == [(20, 20)]


# ---- qualifiers ------------------------------------------------------------------------------------------------------------

def words(n, prefix="word"):
    return " ".join(f"{prefix}{i}" for i in range(n))


def test_qualifier_text_survives_and_is_written_as_documented(tmp_path):
    """Awkward qualifier text (quotes doubled and read back, whitespace runs become one space as documented in document.rs
    `unwrap_text`, long text wraps at spaces inside 80 columns, unicode kept, a flag has no value); one very long word is never
    cut; a value that fills the last line is read by Biopython without a trailing space."""
    # --- awkward text
    long_text = words(70)                                             # about 500 characters, wraps at spaces
    unicode_text = "é ü ñ 日本語 µ → ✓ " * 12
    given = [
        ("label", "special"),
        ("note", 'say "hi" and ""twice"" here'),
        ("note", '"starts and ends with quotes"'),
        ("note", '"'),
        ("note", "line one\nline two\n\nline   three\tand a tab"),
        ("note", long_text),
        ("note", unicode_text),
        ("note", "abc " * 14 + 'x"y "z" ""' + " end"),               # quotes near the wrapping column
        ("product", "x=y; /z=(1..2) [brackets] {braces}"),
        ("translation", "MKV LLA\nGGT*"),                              # whitespace removed in a translation
        ("pseudo", None),
    ]
    doc = Document.new_sequence("quals", "ACGT" * 50, False)
    doc.add_features([("misc_feature", "10..90", given)])
    out = tmp_path / "q.gb"
    doc.save_genbank(out)
    expected = [(k, v if v is None else ("".join(v.split()) if k == "translation" else " ".join(v.split())))
                for k, v in given]
    (back,) = Document(out).features()
    assert back.qualifiers == expected
    assert [q for q in doc.features()[0].qualifiers] == expected      # canonical form is what the document holds too
    text = out.read_text(encoding="utf-8")
    assert max(len(line.encode()) for line in text.splitlines()) <= 80   # bytes: the writer counts bytes, not characters
    assert "日本語" in text and "/pseudo\n" in text
    # Saving the reloaded file again changes nothing.
    again = tmp_path / "q2.gb"
    Document(out).save_genbank(again)
    assert again.read_bytes() == out.read_bytes()
    # Biopython reads the same values (it joins wrapped lines with a space, exactly the canonical form).
    rec = bio_read(out)
    got = bio_qualifiers(rec.features[0])
    assert got["note"] == [v for k, v in expected if k == "note"]
    assert got["translation"] == ["MKVLLAGGT*"]
    assert got["pseudo"] == [""] and got["product"] == ["x=y; /z=(1..2) [brackets] {braces}"]

    # --- a word longer than a line is never cut (cutting would add a space when read back): the line runs long
    word = "x" * 300
    doc = Document.new_sequence("w", "ACGT" * 20, False)
    doc.add_features([("misc_feature", "1..10", [("note", f"before {word} after"), ("note", word)])])
    out = tmp_path / "w.gb"
    doc.save_genbank(out)
    (back,) = Document(out).features()
    assert back.qualifiers == [("note", f"before {word} after"), ("note", word)]

    # --- a value that fills the last line is read without a trailing space by Biopython
    # every length around the page width, so the last line is full for some of them (the closing quote must fit on it)
    for n in range(40, 90):
        note = "sequence: " + dna(n, n)
        doc = Document.new_sequence("w", "ACGT" * 20, False)
        doc.add_features([("primer_bind", "1..10", [("label", "p"), ("note", note)])])
        out = tmp_path / "w.gb"
        doc.save_genbank(out)
        assert Document(out).features()[0].qualifiers[1] == ("note", note)       # our own reader
        assert bio_read(out).features[0].qualifiers["note"] == [note], n         # and Biopython: no trailing space
        # a line is at most 79 wide, but a single unbreakable word (and the closing quote after it) may be longer
        assert all(len(line) <= 79 or " " not in line.strip() for line in out.read_text().splitlines()), n


# ---- primers ---------------------------------------------------------------------------------------------------------------

def primer_features(doc):
    return [(f.label, primer_sequence(f), f.location) for f in doc.features() if f.kind == "primer_bind"]


def test_primers_survive_save(tmp_path):
    """primer_bind features keep name and sequence; the primers of a .dna file survive a save as GenBank; the primers of the
    database are written by the Editor (only those imported to the file) and read back."""
    # --- primer_bind features keep name and sequence
    doc = Document.new_sequence("p", "ACGTTGCATGCCGATTAGCCAAGTTGCCATG" * 4, False)
    given = [
        ("primer_bind", "5..24", [("label", "fwd"), ("note", "sequence: TTTTTACGTTGCATGCCGATTAG")]),
        ("primer_bind", "complement(60..80)", [("label", "rev"), ("note", "sequence: cgatccggCATGGCAACTTGGCTAATC")]),
    ]
    out = tmp_path / "p.gb"
    doc.save_genbank(out, given)                                      # the app passes the current primer sites like this
    back = Document(out)
    assert primer_features(back) == [("fwd", "TTTTTACGTTGCATGCCGATTAG", "5..24"),
                                     ("rev", "cgatccggCATGGCAACTTGGCTAATC", "complement(60..80)")]
    assert len(doc.features()) == 0, "the document itself does not get the primer sites"

    # --- the primers of a .dna file survive the save as GenBank
    if PLASMID_DNA.exists():
        doc = Document(PLASMID_DNA)
        named = dict(doc.primers())
        assert len(named) >= 5
        out = tmp_path / "pdna.gb"
        doc.save_genbank(out)
        back = primer_features(Document(out))
        assert back == primer_features(doc) and back
        assert {n: s for n, s, _ in back} == named                        # every primer of the .dna, name and full sequence

    # --- Editor.save writes a primer_bind for each exact site of a primer that belongs to the file (imported to it); a primer only
    # in the database is not written. Reopened, the sites have name and sequence.
    doc = Document.new_sequence("pe", "GATTACAGGCTTAACCGGTACCATGGCCAATTGGCCTTAAGG" * 3 + "TTTT", True)
    primer = doc.sequence(5, 27)
    primer_db.db().add("myprimer", primer)
    ed = Editor(doc, tmp_path / "pe.gb")
    ed.save(tmp_path / "pe.gb")
    assert not primer_features(Document(tmp_path / "pe.gb")), "only in the database: not part of the file"
    assert ed.import_primers([("myprimer", primer)]) == 1
    ed.save(tmp_path / "pe.gb")
    back = Document(tmp_path / "pe.gb")
    found = primer_features(back)
    assert found and all(n == "myprimer" and s == primer for n, s, _ in found)
    assert any(loc == "6..27" for _, _, loc in found)
    assert not doc.modified and len(doc.history()[0]) == 1, "writing the sites adds no undo step (the import is the one)"


# ---- very short sequences --------------------------------------------------------------------------------------------------

def test_tiny_and_empty_sequences_save_and_reload(tmp_path):
    """One base and other tiny sequences (both topologies) save and reload; an empty sequence too; a sequence edited down to
    nothing saves without error."""
    # --- one base and other tiny sequences
    for bases in ("A", "g", "N", "AC", "acgtNRYKMSWBDHVU"):
        for circular in (False, True):
            doc = Document.new_sequence("tiny", bases, circular)
            doc.add_features([("misc_feature", f"1..{len(bases)}", [("label", "all of it")])])
            out = tmp_path / "tiny.gb"
            doc.save_genbank(out)
            back = Document(out)
            assert (back.sequence(), back.length, back.circular, back.name) == (bases, len(bases), circular, "tiny"), (bases, circular)
            assert snap(back) == snap(doc), (bases, circular)
            rec = bio_read(out)
            assert str(rec.seq) == bases.upper() and len(rec.features) == 1, (bases, circular)

    # --- an empty sequence
    for circular in (False, True):
        doc = Document.new_sequence("empty", "", circular)
        out = tmp_path / "empty.gb"
        doc.save_genbank(out)
        back = Document(out)
        assert (back.length, back.sequence(), back.circular, back.name) == (0, "", circular, "empty"), circular

    # --- a sequence edited down to nothing
    doc = Document.new_sequence("shrunk", "ACGT", False)
    doc.delete(0, 4)
    out = tmp_path / "s.gb"
    doc.save_genbank(out)                                              # saving works ...
    assert out.exists() and not doc.modified and doc.length == 0


# ---- letter case -----------------------------------------------------------------------------------------------------------

def test_gap_letter_and_letter_case_survive(tmp_path):
    """A gap letter survives a GenBank round trip; lower and mixed case bases are kept as typed and written consistently."""
    # --- a gap letter
    doc = Document.new_sequence("gap", "ACGT-ACGT", False)
    out = tmp_path / "gap.gb"
    doc.save_genbank(out)
    assert Document(out).sequence() == "ACGT-ACGT"

    # --- lower and mixed case bases
    src = tmp_path / "case.gb"
    src.write_text(
        "LOCUS       case                      60 bp    DNA     linear UNK 01-JAN-1970\n"
        "FEATURES             Location/Qualifiers\n"
        "     misc_feature    5..20\n"
        '                     /label="x"\n'
        "ORIGIN      \n"
        "        1 acgtacgtac GTACGTACGT nnnnRYKMac gtACGTacgt ACGTACGTac gtacgtacgt\n"
        "//\n")
    doc = Document(src)
    assert doc.sequence() == "acgtacgtacGTACGTACGTnnnnRYKMacgtACGTacgtACGTACGTacgtacgtacgt"
    assert doc.sequence(8, 14) == "acGTAC"                              # slices keep the case too
    first, second = tmp_path / "1.gb", tmp_path / "2.gb"
    doc.save_genbank(first)
    assert Document(first).sequence() == doc.sequence()
    # Written as typed, in the usual 10-base blocks, and the same on every later save.
    origin = first.read_text().split("ORIGIN")[1].splitlines()[1]
    assert origin.split()[1:] == ["acgtacgtac", "GTACGTACGT", "nnnnRYKMac", "gtACGTacgt", "ACGTACGTac", "gtacgtacgt"]
    Document(first).save_genbank(second)
    assert first.read_bytes() == second.read_bytes()
    assert str(bio_read(first).seq) == doc.sequence().upper()
    # An edit keeps the case of what was typed.
    doc.insert(0, "TTT")
    assert doc.sequence()[:6] == "TTTacg"


# ---- FASTA -----------------------------------------------------------------------------------------------------------------

def dna(n, seed):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def fasta_text(records, width=60, eol="\n", final=True, blank_after=0):
    """records: [(header, bases)] -> text wrapped at `width`."""
    lines = []
    for header, bases in records:
        lines.append(">" + header)
        lines += [bases[i:i + width] for i in range(0, len(bases), width)]
    text = eol.join(lines) + (eol if final else "") + eol * blank_after
    return text


def write(path, text):
    path.write_bytes(text.encode("utf-8"))
    return path


RECORDS = [("chr1 first one", dna(250, 1)), ("chr2", dna(130, 2)), ("chr3 last", dna(61, 3)), ("empty_one", "")]


def test_fasta_reading_and_export(tmp_path):
    """A document exported to FASTA is read by Biopython with the same id and sequence; the multi-record index reads each record;
    line-ending, blank-line and case variants all give the same bases."""
    # --- the FASTA of a document is read by Biopython with the same id and sequence
    src = write(tmp_path / "src.fa", fasta_text(RECORDS))
    doc = Document.from_fasta(src, "chr2")
    assert (doc.name, doc.length, doc.circular, doc.sequence()) == ("chr2", 130, False, RECORDS[1][1])
    dst = tmp_path / "export.fa"
    doc.write_fasta(src, dst)
    recs = {r.id: r for r in SeqIO.parse(dst, "fasta")}
    assert list(recs) == ["chr1", "chr2", "chr3", "empty_one"]
    assert str(recs["chr2"].seq) == doc.sequence()
    for (header, bases) in (RECORDS[0], RECORDS[2]):
        assert str(recs[header.split()[0]].seq) == bases
    assert recs["chr1"].description == "chr1 first one"
    # An edited document is exported the same way.
    doc.insert(10, "GATTACA")
    doc.write_fasta(src, dst)
    assert str({r.id: r for r in SeqIO.parse(dst, "fasta")}["chr2"].seq) == doc.sequence()
    assert doc.modified, "exporting does not mark the document saved (the caller does)"

    # --- the multi-record index reads each record
    path = write(tmp_path / "multi.fa", fasta_text(RECORDS))
    fa = Fasta(path)
    assert fa.sequences() == [("chr1", 250), ("chr2", 130), ("chr3", 61), ("empty_one", 0)]
    bio = {r.id: str(r.seq) for r in SeqIO.parse(path, "fasta")}
    for name, length in fa.sequences():
        assert fa.fetch(name, 0, length) == bio[name]
        assert fa.fetch(name, 0, 10_000) == bio[name], "an end past the sequence is clipped"
    # Pieces across line breaks, and the first and last base.
    chr1 = bio["chr1"]
    for start, end in [(0, 1), (59, 61), (55, 130), (249, 250), (119, 121), (0, 250)]:
        assert fa.fetch("chr1", start, end) == chr1[start:end]
    # Each record opens as its own document.
    for name, length in fa.sequences():
        assert Document.from_fasta(path, name).sequence() == bio[name]
    assert (tmp_path / "multi.fa.fai").exists(), "the index is written beside the file"
    assert Fasta(path).sequences() == fa.sequences(), "and reused on the next open"

    # --- line ending and case variants give the same bases
    bases = {"chr1": dna(150, 7), "chr2": dna(77, 8)}
    records = list(bases.items())
    variants = {
        "crlf": fasta_text(records, eol="\r\n"),
        "no final newline": fasta_text(records, final=False),
        "trailing blank lines": fasta_text(records, blank_after=3),
        "lower case": fasta_text([(h, b.lower()) for h, b in records]),
        "crlf no final newline": fasta_text(records, eol="\r\n", final=False),
        "crlf blank lines": fasta_text(records, eol="\r\n", blank_after=2),
        "mixed case": fasta_text([(h, b[:50] + b[50:].lower()) for h, b in records]),
        "one long line": fasta_text(records, width=10_000),
    }
    for i, (variant, text) in enumerate(variants.items()):
        path = write(tmp_path / f"v{i}.fa", text)                       # one file per variant: no stale index from the last one
        fa = Fasta(path)
        bio = {r.id: str(r.seq) for r in SeqIO.parse(path, "fasta")}
        assert fa.sequences() == [("chr1", 150), ("chr2", 77)], variant
        for name, length in fa.sequences():
            got = fa.fetch(name, 0, length)
            assert got.upper() == bases[name] and got == bio[name], f"{variant}: no \\r, no blank line in the bases; case as in the file"
            assert Document.from_fasta(path, name).sequence() == got, variant
        assert fa.fetch("chr1", 148, 152) == bio["chr1"][148:150], variant


def test_fasta_writing_after_an_edit(tmp_path):
    """Writing a FASTA after an edit equals the string edit of the file (to another file and in place), and leaves the other
    records untouched even with CRLF."""
    # --- equals the string edit of the file
    records = [(h, b) for h, b in RECORDS if b]
    text = fasta_text(records)
    src = write(tmp_path / "genome.fa", text)
    doc = Document.from_fasta(src, "chr2")
    doc.insert(5, "GGGGGGG")
    doc.delete(100, 110)
    doc.replace(0, 3, "tt")
    new_chr2 = RECORDS[1][1]                                            # the same edits on a plain string, in the same order
    new_chr2 = new_chr2[:5] + "GGGGGGG" + new_chr2[5:]
    new_chr2 = new_chr2[:100] + new_chr2[110:]
    new_chr2 = "tt" + new_chr2[3:]
    assert doc.sequence() == new_chr2
    expected = fasta_text([(h, new_chr2 if h.startswith("chr2") else b) for h, b in records])
    # To another file: the source stays as it was, the copy is exactly the string edit.
    dst = tmp_path / "edited.fa"
    doc.write_fasta(src, dst)
    assert dst.read_bytes() == expected.encode() and src.read_bytes() == text.encode()
    # In place, as File > Save does: same result; the old index is replaced by a good one.
    Fasta(src)                                                          # make the index first, so there is an old one to replace
    doc.write_fasta(src, src)
    assert src.read_bytes() == expected.encode()
    fa = Fasta(src)
    assert fa.fetch("chr2", 0, 10_000) == new_chr2 and fa.fetch("chr1", 0, 10_000) == RECORDS[0][1]
    assert fa.fetch("chr3", 0, 10_000) == RECORDS[2][1]
    assert not list(tmp_path.glob("*.tmp"))

    # --- the other records are untouched even with CRLF
    text = fasta_text(RECORDS[:3], eol="\r\n")
    src = write(tmp_path / "crlf.fa", text)
    doc = Document.from_fasta(src, "chr2")
    doc.insert(0, "AAAA")
    doc.write_fasta(src, src)
    out = src.read_bytes()
    raw = text.encode()
    before, after = raw[:raw.index(b">chr2")], raw[raw.index(b">chr3"):]
    assert out.startswith(before + b">chr2") and out.endswith(after), "chr1 and chr3 are copied byte for byte"
    # (the edited record itself is written with plain "\n" line ends; the index and Biopython cope with the mix)
    assert Fasta(src).fetch("chr2", 0, 999) == "AAAA" + RECORDS[1][1]
    recs = {r.id: str(r.seq) for r in SeqIO.parse(src, "fasta")}
    assert recs == {"chr1": RECORDS[0][1], "chr2": "AAAA" + RECORDS[1][1], "chr3": RECORDS[2][1]}


# ---- saving never damages a file -------------------------------------------------------------------------------------------

def test_a_failed_save_leaves_the_original_alone(tmp_path):
    """A failed GenBank save leaves the original file and the document alone (bad primer site, missing folder, a folder in the
    way) and the real save afterwards works; a failed FASTA write leaves the source alone; a failed save in the editor keeps the
    path and the modified flag."""
    # --- a failed GenBank save
    out = tmp_path / "keep.gb"
    doc = rich_document()
    doc.save_genbank(out)
    original = out.read_bytes()
    doc.insert(5, "GATTACA")                                            # unsaved edit
    assert doc.modified
    # A primer feature that does not fit the sequence is refused before anything is written.
    with pytest.raises(OSError):
        doc.save_genbank(out, [("primer_bind", "1..99999", [("label", "bad")])])
    assert out.read_bytes() == original and not list(tmp_path.glob("*.tmp"))
    assert doc.modified, "a failed save does not mark the document saved"
    # A folder that does not exist: an error, nothing created.
    with pytest.raises(OSError):
        doc.save_genbank(tmp_path / "no_such_folder" / "x.gb")
    assert not (tmp_path / "no_such_folder").exists() and doc.modified
    # A name that is taken by a folder: an error, the folder is untouched.
    folder = tmp_path / "taken.gb"
    folder.mkdir()
    (folder / "inside.txt").write_text("x")
    with pytest.raises(OSError):
        doc.save_genbank(folder)
    assert folder.is_dir() and (folder / "inside.txt").read_text() == "x" and doc.modified
    # And the real save afterwards works and replaces the file in one go.
    doc.save_genbank(out)
    assert out.read_bytes() != original and not doc.modified and Document(out).sequence() == doc.sequence()

    # --- a failed FASTA write leaves the source alone
    text = fasta_text(RECORDS[:2])
    src = write(tmp_path / "g.fa", text)
    doc = Document.new_sequence("not_in_the_file", "ACGT", False)
    with pytest.raises(OSError, match="no sequence named"):
        doc.write_fasta(src, src)
    assert src.read_bytes() == text.encode() and not list(tmp_path.glob("*.tmp"))
    ok = Document.from_fasta(src, "chr1")
    with pytest.raises(OSError):
        ok.write_fasta(src, tmp_path / "no_such_folder" / "x.fa")
    assert src.read_bytes() == text.encode()

    # --- a failed save in the editor keeps the path and the modified flag
    doc = Document.new_sequence("ed", "ACGT" * 10, False)
    ed = Editor(doc, tmp_path / "ed.gb")
    ed.set_cursor(0)
    assert ed.paste("GATTACA") and doc.modified
    with pytest.raises(OSError):
        ed.save(tmp_path / "no_such_folder" / "ed.gb")
    assert ed.path == tmp_path / "ed.gb" and doc.modified
    ed.save(tmp_path / "ed.gb")
    assert not doc.modified and Document(tmp_path / "ed.gb").sequence() == "GATTACA" + "ACGT" * 10


# ---- autosave and recovery -------------------------------------------------------------------------------------------------

@pytest.fixture
def work_file(tmp_path):
    """A copy of a demo vector, with an old modification time (the recovery copy is judged by being newer)."""
    if not DEMOS:
        pytest.skip("no demo files")
    path = tmp_path / "vector.gb"
    path.write_bytes((ROOT / "demo" / "demo_vector.gb").read_bytes())
    old = path.stat().st_mtime - 100
    os.utime(path, (old, old))
    return path


def open_in_window(path, monkeypatch, answer=None):
    asked = []

    def ask(window, path, recovery):   # the "Recover Unsaved Changes" prompt: Recover Changes (Yes) or Discard Changes
        asked.append(f"{Path(path).name} has unsaved changes ({recovery.name})")
        return answer == QMessageBox.Yes

    monkeypatch.setattr(M.MainWindow, "_ask_recover", ask)
    win = M.MainWindow()
    assert win.load(str(path))
    return win, asked


def test_autosave_and_recovery(work_file, tmp_path, monkeypatch):
    """Autosave (every minute) writes a recovery copy and never touches the original (also for a .dna file, which gets a
    .autosave.gb beside it); declined or stale recovery copies are dropped and the original opens; unsaved edits, features
    included, are offered again after a crash and recovered; saving then removes the copy."""
    original = work_file.read_bytes()
    # --- autosave writes a recovery copy and never touches the original
    win, _ = open_in_window(work_file, monkeypatch)
    assert win.autosave_timer.isActive() and win.autosave_timer.interval() == M.AUTOSAVE_MS == 60_000
    ed = win.editor
    win.autosave()
    assert not M.autosave_path(work_file).exists(), "nothing changed: nothing written"
    ed.set_cursor(10)
    assert ed.paste("GATTACA")
    win.autosave()
    recovery = M.autosave_path(work_file)
    assert recovery == work_file.with_name("vector.gb.autosave.gb") and recovery.exists()
    assert work_file.read_bytes() == original, "the user's file is only written by Save"
    assert ed.doc.modified, "autosave does not mark the document saved"
    copy = Document(recovery)
    assert copy.sequence() == ed.doc.sequence() and len(copy.sequence()) == len(Document(work_file).sequence()) + 7
    assert not list(work_file.parent.glob("*.tmp"))
    stamp = recovery.stat().st_mtime_ns
    win.autosave()                                                     # no new change since: not rewritten
    assert recovery.stat().st_mtime_ns == stamp
    # A second edit, a second autosave: the copy follows.
    assert ed.paste("CCC") and (win.autosave(), Document(recovery).sequence() == ed.doc.sequence())[1]
    assert work_file.read_bytes() == original
    ed.doc.modified = False
    win.close()
    assert not recovery.exists(), "closing after discarding removes the recovery copy"

    # --- declined or stale recovery copies are dropped and the original opens
    win, _ = open_in_window(work_file, monkeypatch)
    win.editor.set_cursor(0)
    win.editor.paste("AAAA")
    win.autosave()
    recovery = M.autosave_path(work_file)
    win.autosave_file = None
    win.editor.doc.modified = False
    win.close()
    # Declined.
    win2, asked = open_in_window(work_file, monkeypatch, answer=QMessageBox.No)
    assert len(asked) == 1 and not recovery.exists()
    assert win2.editor.doc.sequence() == Document(work_file).sequence() and not win2.editor.doc.modified
    win2.close()
    # Stale: older than the file itself (the file was saved after): dropped without asking.
    recovery.write_bytes(original)
    old = work_file.stat().st_mtime - 50
    os.utime(recovery, (old, old))
    win3, asked = open_in_window(work_file, monkeypatch, answer=QMessageBox.Yes)
    assert asked == [] and not recovery.exists() and not win3.editor.doc.modified
    win3.close()

    # --- a .dna file's autosave goes beside it and the .dna is never written
    if PLASMID_DNA.exists():
        work = tmp_path / "plasmid.dna"
        work.write_bytes(PLASMID_DNA.read_bytes())
        win, _ = open_in_window(work, monkeypatch)
        win.editor.set_cursor(0)
        win.editor.paste("GATTACA")
        win.autosave()
        recovery_dna = M.autosave_path(work)
        assert recovery_dna.name == "plasmid.dna.autosave.gb" and recovery_dna.exists()
        assert work.read_bytes() == PLASMID_DNA.read_bytes()
        assert Document(recovery_dna).sequence() == win.editor.doc.sequence()
        win.editor.doc.modified = False
        win.close()

    # --- unsaved edits are offered again after a crash and recovered (last: it saves, changing the original)
    win, _ = open_in_window(work_file, monkeypatch)
    win.editor.set_cursor(0)
    win.editor.paste("TTTTTTTT")
    win.editor.add_features([("misc_feature", "30..50", [("label", "mine")])])   # a feature added by the user comes back too
    edited = win.editor.doc.sequence()
    win.autosave()
    recovery = M.autosave_path(work_file)
    # The crash: the window is gone, the recovery file stays.
    win.autosave_file = None
    win.editor.doc.modified = False
    win.close()
    assert recovery.exists() and work_file.read_bytes() == original
    # Next start: asked, answer yes -> the edited document, still unsaved, features intact.
    win2, asked = open_in_window(work_file, monkeypatch, answer=QMessageBox.Yes)
    assert len(asked) == 1 and "unsaved changes" in asked[0]
    assert win2.editor.doc.sequence() == edited and win2.editor.doc.modified and win2.windowTitle().startswith("*")
    assert len(win2.editor.features) == len(Document(work_file).features()) + 1
    assert [f.location for f in win2.editor.features if f.label == "mine"] == ["30..50"]
    assert work_file.read_bytes() == original, "recovering does not write the original either"
    # Saving it makes it permanent and removes the recovery copy.
    assert win2.save()
    assert Document(work_file).sequence() == edited and not recovery.exists()
    win2.close()


# ---- corrupt and truncated files -------------------------------------------------------------------------------------------

def cut_points(n):
    return sorted({*range(1, 40), *range(40, n, max(1, n // 120))})


def test_a_truncated_file_is_a_clear_error_never_a_crash_or_a_silent_partial_document(tmp_path):
    for label, source in (("dna", PLASMID_DNA), ("genbank", DEMOS[0] if DEMOS else PLASMID_GB)):
        if not source.exists():
            continue
        data = source.read_bytes()
        full = Document(source)
        path = tmp_path / ("cut" + source.suffix)
        errors = 0
        for cut in cut_points(len(data)):
            path.write_bytes(data[:len(data) - cut])
            try:
                doc = Document(path)
            except OSError as e:                                            # the only kind of failure: a normal Python error
                errors += 1
                assert str(e), "with a message"
            else:                                                           # opened: then it must be the whole record
                assert doc.sequence() == full.sequence() and len(doc.features()) == len(full.features()), (label, cut)
        assert errors > 100, label


def test_corrupt_and_foreign_files_give_an_error_message(tmp_path, monkeypatch):
    """Corrupt and foreign files (also a damaged .dna) give an error message; a GenBank file with two records is refused instead
    of silently dropping one; opening a corrupt file in the window shows an error and keeps the open document."""
    # --- corrupt and foreign files
    cases = {
        "empty.gb": b"",
        "empty.dna": b"\x09",
        "text.gb": b"hello world\nthis is not a sequence file\n",
        "binary.gb": bytes(range(256)) * 4,
        "locus_only.gb": b"LOCUS       x                 10 bp    DNA     linear UNK 01-JAN-1970\n",
        "fasta_as.gb": b">a\nACGT\n",
        "bad_packet.dna": b"\x09\x00\x00\x00\x10" + b"\xff" * 5,
    }
    for name, content in cases.items():
        path = tmp_path / name
        path.write_bytes(content)
        with pytest.raises(OSError) as err:
            Document(path)
        assert str(err.value).strip(), name
    with pytest.raises(FileNotFoundError):
        Document(tmp_path / "missing.gb")
    # Damage in the middle of a real .dna file.
    if PLASMID_DNA.exists():
        data = bytearray(PLASMID_DNA.read_bytes())
        for i in range(100, len(data), 53):
            data[i] = 0xFF
        bad = tmp_path / "scrambled.dna"
        bad.write_bytes(bytes(data))
        with pytest.raises(OSError):
            Document(bad)

    if not DEMOS:
        return
    # --- a GenBank file with two records is refused instead of silently dropping one
    one = DEMOS[0].read_text()
    path = tmp_path / "two.gb"
    path.write_text(one + one)
    with pytest.raises(OSError, match="2 GenBank records"):
        Document(path)

    # --- opening a corrupt file in the window shows an error and keeps the open document
    bad = tmp_path / "bad.gb"
    bad.write_bytes(DEMOS[0].read_bytes()[:300])
    shown = []
    monkeypatch.setattr(QMessageBox, "critical", lambda *a, **k: shown.append(a[1:]))
    win = M.MainWindow()
    assert win.load(str(DEMOS[0]))
    name, length = win.editor.doc.name, win.editor.doc.length
    assert win.load(str(bad)) is False
    assert len(shown) == 1 and "Cannot Open File" in shown[0][0] and "bad.gb" in shown[0][1]
    assert (win.editor.doc.name, win.editor.doc.length) == (name, length)
    win.close()
