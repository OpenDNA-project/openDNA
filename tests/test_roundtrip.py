"""File round trips: GenBank / .dna / FASTA files written by us read back the same (by us and by Biopython), odd content
survives, and the real demo files keep everything through a save.

Real files (copies of demo/, see sample_files.py) are used when present (skipped otherwise); everything is written in tmp_path."""
import os
import random
import sys
import tempfile
from pathlib import Path

import pytest

os.environ["OPENDNA_DATA_DIR"] = tempfile.mkdtemp(prefix="opendna_test_")
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "app"))

from Bio import SeqIO

from opendna_core import Document, Fasta

import primer_db
from editor import Editor, arrow, is_directional, primer_sequence
from sample_files import GENOME as GENOME_DNA, PLASMID as PLASMID_GB, PLASMID2, DNA_FILE as PLASMID_DNA   # copies

# Biopython warns about the width of the LOCUS line we write for short names; the parse itself is right.
pytestmark = pytest.mark.filterwarnings("ignore::Bio.BiopythonParserWarning")

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


def genbank_sources(extra=()):
    """(id, path) of every real or demo file that is present (plus the `extra` ones when present)."""
    out = [("plasmid_dna", PLASMID_DNA), ("plasmid_gb", PLASMID_GB)] + [(p.stem, p) for p in DEMOS] + list(extra)
    return [(i, p) for i, p in out if p.exists()]


def dna(n, seed):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def fasta_text(records, width=60, eol="\n", final=True, blank_after=0):
    """records: [(header, bases)] -> text wrapped at `width`."""
    lines = []
    for header, bases in records:
        lines.append(">" + header)
        lines += [bases[i:i + width] for i in range(0, len(bases), width)]
    return eol.join(lines) + (eol if final else "") + eol * blank_after


def write(path, text):
    path.write_bytes(text.encode("utf-8"))
    return path


RECORDS = [("chr1 first one", dna(250, 1)), ("chr2", dna(130, 2)), ("chr3 last", dna(61, 3)), ("empty_one", "")]
FA = ">chr1 first\nACGTACGTAC\nGTAC\n>chr2\nTTTTTTTTTT\nTT\n"


# ---- real files --------------------------------------------------------------------------------------------------------------

def test_real_files_round_trip_and_agree_with_biopython(tmp_path):
    """The real and demo files (the 4.6 Mb genome too): load, save as GenBank and reload is identical; saving again is
    byte-identical; Biopython reads our file the same. The .dna plasmid matches its GenBank export (bases, primers, feature
    direction), and an edit inside a feature saved as GenBank moves the feature, Biopython agrees, and undo puts it back."""
    sources = genbank_sources(extra=[("genome_dna", GENOME_DNA)])
    if not sources:
        pytest.skip("no real or demo files")
    # --- load, save as GenBank, reload is identical; save, reload, save is byte-identical
    for pid, path in sources:
        doc = Document(path)
        before = snap(doc)
        assert before["features"], f"{pid}: the file should have annotations to make this meaningful"
        # --- what we read from the ORIGINAL file = Biopython's reading of it (bases, every feature, part and qualifier; the
        # genome too: about a second). Otherwise a reader that drops a qualifier on open saves, reloads and agrees with itself.
        if path.suffix == ".gb":
            assert_biopython_agrees(path, doc)
        if pid == "genome_dna":   # the 4.6 Mb E. coli genome, the largest file
            assert doc.length > 4_500_000 and doc.circular and len(doc.features()) > 4000
        first, second, third = tmp_path / "1.gb", tmp_path / "2.gb", tmp_path / "3.gb"
        doc.save_genbank(first)
        back = Document(first)
        after = snap(back)
        for key in before:
            assert after[key] == before[key], f"{pid}: {key}"
        assert [f.spans for f in back.features()] == [f.spans for f in doc.features()], pid
        assert [f.strand for f in back.features()] == [f.strand for f in doc.features()], pid
        assert not back.modified, pid
        back.save_genbank(second)
        Document(second).save_genbank(third)
        assert first.read_bytes() == second.read_bytes() == third.read_bytes(), pid
        # --- a file saved by us is read the same by Biopython (not the genome: Biopython is slow on it)
        if pid != "genome_dna":
            assert_biopython_agrees(first, doc)

    # --- the bundled example plasmids open with their sizes
    for path, length in ((PLASMID_GB, 5936), (PLASMID2, 3780)):
        if path.exists():
            doc = Document(path)
            assert doc.features() and doc.circular and doc.length == length, path.name
            assert_biopython_agrees(path, doc)   # pFG051_Chlor.gb is not among the sources above

    if not (PLASMID_DNA.exists() and PLASMID2.exists()):
        return
    # --- the .dna file matches its GenBank export: same bases, its primer list, directions read from the file
    dna_doc, gb = Document(PLASMID_DNA), Document(PLASMID2)
    assert dna_doc.length == gb.length == 3780 and dna_doc.circular and gb.circular
    assert dna_doc.sequence().upper() == gb.sequence().upper()
    assert len(dna_doc.primers()) == 15, "the .dna file's primer list"
    feats = dna_doc.features()
    cds = [f for f in feats if f.kind == "CDS"]
    assert cds and all(arrow(f) != 0 for f in cds), "CDSs keep their arrow"
    assert any(arrow(f) == 0 for f in feats if f.kind == "misc_feature"), "a .dna file marks features without direction"
    assert not is_directional("misc_feature", []) and is_directional("CDS", [])
    assert is_directional("misc_feature", [("directional", "yes")]) and not is_directional("promoter", [("directional", "no")])

    # --- edit the .dna inside a feature, save as GenBank, Biopython agrees, undo puts it back
    doc = Document(PLASMID_DNA)
    target = next(f for f in doc.features() if f.kind != "source" and len(f.spans) == 1
                  and f.spans[0][1] - f.spans[0][0] > 20)
    start, end = target.spans[0]
    doc.insert(start + 10, "GATTACA")
    assert doc.features()[target.index].spans == [(start, end + 7)] and doc.modified
    out = tmp_path / "edited.gb"
    doc.save_genbank(out)
    assert not doc.modified
    rec = SeqIO.read(out, "genbank")
    assert str(rec.seq) == doc.sequence().upper() and len(rec.features) == len(doc.features())
    bio = rec.features[target.index]
    assert (int(bio.location.start), int(bio.location.end)) == (start, end + 7)
    assert doc.undo() and doc.features()[target.index].spans == [(start, end)]


# ---- made-up records ----------------------------------------------------------------------------------------------------------

LOCATIONS = [
    "10..20", "complement(10..20)", "join(10..20,30..40)", "complement(join(10..20,30..40))",
    "join(complement(30..40),complement(10..20))", "join(10..20,complement(30..40))",
    "<10..20", "10..>20", "<10..>20", "order(10..20,30..40)", "join(10..20,30..40,50..60)",
    "join(350..400,1..20)",                                           # across the origin of the circle
]


def test_unusual_records_round_trip(tmp_path):
    """Records at the edges of the format survive a GenBank save and reload, and Biopython reads them the same: joins,
    complements, fuzzy ends and order() (byte-identical on a second save); every location form, written as given; a between
    location; both topologies (and a topology toggle, one undo step); one-base, tiny and empty sequences; a gap letter; lower
    and mixed case kept as typed."""
    # --- a rich record: save, reload, save is byte-identical; Biopython agrees, circular and linear
    first, second = tmp_path / "1.gb", tmp_path / "2.gb"
    doc = rich_document()
    doc.save_genbank(first)
    back = Document(first)
    assert snap(back) == snap(doc)
    back.save_genbank(second)
    assert first.read_bytes() == second.read_bytes()
    for circular in (True, False):
        doc = rich_document(circular)
        out = tmp_path / "rich.gb"
        doc.save_genbank(out)
        assert_biopython_agrees(out, doc)
        by_label = {f.qualifiers.get("label", [""])[0]: f for f in bio_read(out).features}
        assert by_label["fuzzy both ends"].location.start.__class__.__name__ == "BeforePosition", circular
        assert by_label["fuzzy both ends"].location.end.__class__.__name__ == "AfterPosition", circular
        assert by_label["complement of a join"].location.strand == -1, circular

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
    # `20^21` is the gap between bases 20 and 21 (1-based): Biopython puts it at 0-based position 20
    doc = Document.new_sequence("gap", "ACGT" * 20, False)
    doc.add_features([("misc_feature", "20^21", [("label", "gap")])])
    out = tmp_path / "gap.gb"
    doc.save_genbank(out)
    (f,) = Document(out).features()
    assert f.location == "20^21" and f.spans == [(20, 20)]
    rec = bio_read(out)
    assert (int(rec.features[0].location.start), int(rec.features[0].location.end)) == (20, 20)

    # --- topology, and a toggle after loading (one undo step)
    for circular in (True, False):
        doc = Document.new_sequence("topo", "ACGTACGTTTGACCA" * 7, circular)
        out = tmp_path / "t.gb"
        doc.save_genbank(out)
        assert (" circular " in out.read_text()) == circular and (" linear " in out.read_text()) != circular
        back = Document(out)
        assert back.circular is circular and back.sequence() == doc.sequence()
        back.set_topology(not circular)
        back.save_genbank(out)
        assert Document(out).circular is (not circular)
        assert back.undo() and back.circular is circular

    # --- one base and other tiny sequences, both topologies
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
    # an empty sequence, and one edited down to nothing
    for circular in (False, True):
        doc = Document.new_sequence("empty", "", circular)
        out = tmp_path / "empty.gb"
        doc.save_genbank(out)
        back = Document(out)
        assert (back.length, back.sequence(), back.circular, back.name) == (0, "", circular, "empty"), circular
    # edited down to nothing: the file written holds no bases (not the 4 it had), reloads empty, Biopython reads 0 bases
    for circular in (False, True):
        doc = Document.new_sequence("shrunk", "ACGT", circular)
        doc.delete(0, 4)
        out = tmp_path / "s.gb"
        doc.save_genbank(out)
        assert out.exists() and not doc.modified and doc.length == 0
        back = Document(out)
        assert (back.length, back.sequence(), back.circular, back.name) == (0, "", circular, "shrunk"), circular
        text = out.read_text(encoding="utf-8")
        assert "ORIGIN" in text and not any(c in "ACGTacgt" for c in text.split("ORIGIN", 1)[1]), circular
        assert len(bio_read(out).seq) == 0, circular

    # --- a gap letter
    doc = Document.new_sequence("gap", "ACGT-ACGT", False)
    out = tmp_path / "gapletter.gb"
    doc.save_genbank(out)
    assert Document(out).sequence() == "ACGT-ACGT"

    # --- lower and mixed case bases kept as typed, written in 10-base blocks, the same on every save
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
    doc.save_genbank(first)
    assert Document(first).sequence() == doc.sequence()
    origin = first.read_text().split("ORIGIN")[1].splitlines()[1]
    assert origin.split()[1:] == ["acgtacgtac", "GTACGTACGT", "nnnnRYKMac", "gtACGTacgt", "ACGTACGTac", "gtacgtacgt"]
    Document(first).save_genbank(second)
    assert first.read_bytes() == second.read_bytes()
    assert str(bio_read(first).seq) == doc.sequence().upper()
    doc.insert(0, "TTT")
    assert doc.sequence()[:6] == "TTTacg", "an edit keeps the case of what was typed"


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
    again = tmp_path / "q2.gb"
    Document(out).save_genbank(again)
    assert again.read_bytes() == out.read_bytes(), "saving the reloaded file again changes nothing"
    # Biopython reads the same values (it joins wrapped lines with a space, exactly the canonical form).
    got = bio_qualifiers(bio_read(out).features[0])
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
    database are written by the Editor (only those imported to the file) and read back, adding no undo step."""
    # --- primer_bind features keep name and sequence
    doc = Document.new_sequence("p", "ACGTTGCATGCCGATTAGCCAAGTTGCCATG" * 4, False)
    given = [
        ("primer_bind", "5..24", [("label", "fwd"), ("note", "sequence: TTTTTACGTTGCATGCCGATTAG")]),
        ("primer_bind", "complement(60..80)", [("label", "rev"), ("note", "sequence: cgatccggCATGGCAACTTGGCTAATC")]),
    ]
    out = tmp_path / "p.gb"
    doc.save_genbank(out, given)                                      # the app passes the current primer sites like this
    assert primer_features(Document(out)) == [("fwd", "TTTTTACGTTGCATGCCGATTAG", "5..24"),
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

    # --- Editor.save writes a primer_bind for each exact site of a primer imported to the file; one only in the database is not
    doc = Document.new_sequence("pe", "GATTACAGGCTTAACCGGTACCATGGCCAATTGGCCTTAAGG" * 3 + "TTTT", True)
    primer = doc.sequence(5, 27)
    primer_db.db().add("myprimer", primer)
    ed = Editor(doc, tmp_path / "pe.gb")
    ed.save(tmp_path / "pe.gb")
    assert not primer_features(Document(tmp_path / "pe.gb")), "only in the database: not part of the file"
    assert ed.import_primers([("myprimer", primer)]) == 1
    ed.save(tmp_path / "pe.gb")
    found = primer_features(Document(tmp_path / "pe.gb"))
    assert found and all(n == "myprimer" and s == primer for n, s, _ in found)
    assert any(loc == "6..27" for _, _, loc in found)
    assert not doc.modified and len(doc.history()[0]) == 1, "writing the sites adds no undo step (the import is the one)"


# ---- FASTA -----------------------------------------------------------------------------------------------------------------

def test_fasta_files_are_read_like_biopython(tmp_path):
    """The multi-record FASTA index reads each record (pieces across line breaks, ends clipped, index written and reused);
    line-ending, blank-line, case and byte-order-mark variants all give the same bases; each record opens as a document."""
    # --- the multi-record index reads each record
    path = write(tmp_path / "multi.fa", fasta_text(RECORDS))
    fa = Fasta(path)
    assert fa.sequences() == [("chr1", 250), ("chr2", 130), ("chr3", 61), ("empty_one", 0)]
    bio = {r.id: str(r.seq) for r in SeqIO.parse(path, "fasta")}
    for name, length in fa.sequences():
        assert fa.fetch(name, 0, length) == bio[name]
        assert fa.fetch(name, 0, 10_000) == bio[name], "an end past the sequence is clipped"
        assert Document.from_fasta(path, name).sequence() == bio[name], "each record opens as its own document"
    for start, end in [(0, 1), (59, 61), (55, 130), (249, 250), (119, 121), (0, 250)]:   # across line breaks, first and last base
        assert fa.fetch("chr1", start, end) == bio["chr1"][start:end]
    doc = Document.from_fasta(path, "chr2")
    assert (doc.name, doc.length, doc.circular) == ("chr2", 130, False)
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

    # --- a byte order mark before the first header
    path = tmp_path / "bom.fa"
    path.write_bytes(b"\xef\xbb\xbf" + FA.encode())
    fasta = Fasta(str(path))
    assert fasta.sequences() == [("chr1", 14), ("chr2", 12)]
    assert fasta.fetch("chr1", 0, 14) == "ACGTACGTACGTAC" and fasta.fetch("chr2", 0, 4) == "TTTT"
    assert Document.from_fasta(str(path), "chr1").sequence() == "ACGTACGTACGTAC"


def test_fasta_writing_is_the_string_edit_of_the_file(tmp_path):
    """Writing a FASTA after an edit equals the string edit of the file (to another file and in place, replacing the index),
    Biopython reads the export with the same ids; the other records are copied byte for byte, with CRLF or a byte order mark
    kept; exporting does not mark the document saved."""
    # --- the export is read by Biopython with the same ids and sequences
    src = write(tmp_path / "src.fa", fasta_text(RECORDS))
    doc = Document.from_fasta(src, "chr2")
    dst = tmp_path / "export.fa"
    doc.insert(10, "GATTACA")
    doc.write_fasta(src, dst)
    recs = {r.id: r for r in SeqIO.parse(dst, "fasta")}
    assert list(recs) == ["chr1", "chr2", "chr3", "empty_one"] and recs["chr1"].description == "chr1 first one"
    assert str(recs["chr2"].seq) == doc.sequence()
    for header, bases in (RECORDS[0], RECORDS[2]):
        assert str(recs[header.split()[0]].seq) == bases
    assert doc.modified, "exporting does not mark the document saved (the caller does)"

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
    dst = tmp_path / "edited.fa"
    doc.write_fasta(src, dst)
    assert dst.read_bytes() == expected.encode() and src.read_bytes() == text.encode(), "to another file: the source is untouched"
    Fasta(src)                                                          # make the index first, so there is an old one to replace
    doc.write_fasta(src, src)                                           # in place, as File > Save does
    assert src.read_bytes() == expected.encode()
    fa = Fasta(src)
    assert fa.fetch("chr2", 0, 10_000) == new_chr2 and fa.fetch("chr1", 0, 10_000) == RECORDS[0][1]
    assert fa.fetch("chr3", 0, 10_000) == RECORDS[2][1]
    assert not list(tmp_path.glob("*.tmp"))

    # --- the other records are untouched with CRLF, and the edited one is written with CRLF too
    text = fasta_text(RECORDS[:3], eol="\r\n")
    src = write(tmp_path / "crlf.fa", text)
    doc = Document.from_fasta(src, "chr2")
    doc.insert(0, "AAAA")
    doc.write_fasta(src, src)
    out, raw = src.read_bytes(), text.encode()
    before, after = raw[:raw.index(b">chr2")], raw[raw.index(b">chr3"):]
    assert out.startswith(before + b">chr2") and out.endswith(after), "chr1 and chr3 are copied byte for byte"
    assert Fasta(src).fetch("chr2", 0, 999) == "AAAA" + RECORDS[1][1]
    assert {r.id: str(r.seq) for r in SeqIO.parse(src, "fasta")} == {
        "chr1": RECORDS[0][1], "chr2": "AAAA" + RECORDS[1][1], "chr3": RECORDS[2][1]}
    path = tmp_path / "crlf2.fa"
    path.write_bytes(FA.replace("\n", "\r\n").encode())
    out = tmp_path / "out.fa"
    Document.new_sequence("chr2", "G" * 70, False).write_fasta(str(path), str(out))
    data = out.read_bytes()
    assert b"\n" not in data.replace(b"\r\n", b""), "every line ends with CRLF"
    assert data.startswith(b">chr1 first\r\nACGTACGTAC\r\nGTAC\r\n>chr2\r\n" + b"G" * 60 + b"\r\n" + b"G" * 10 + b"\r\n")

    # --- a byte order mark on the first header is kept when that record is replaced
    path = tmp_path / "bom.fa"
    path.write_bytes(b"\xef\xbb\xbf" + FA.encode())
    Document.new_sequence("chr1", "GGGG", False).write_fasta(str(path), str(out))
    assert out.read_bytes() == b"\xef\xbb\xbf>chr1 first\nGGGG\n>chr2\nTTTTTTTTTT\nTT\n"
