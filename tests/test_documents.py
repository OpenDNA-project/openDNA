"""Documents and the files behind them, without a window: saving never damages a file, corrupt files are refused, editing is
exact through undo / redo (bases, features, edited-base ranges, primer steps of the History), copy keeps annotations, the text
exports and the Open DNA Project file, and the gene groups of the features."""
import json
import os
import random
import sys
import tempfile
import threading
import time
import zipfile
from pathlib import Path
from types import SimpleNamespace

import pytest

os.environ["OPENDNA_DATA_DIR"] = tempfile.mkdtemp(prefix="opendna_test_")
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "app"))

from opendna_core import Document

import exports
import fileio
import primer_db
import selftest
from editor import Editor, reverse_complement
from file_history import FileHistory
from history_panel import group
from main import autosave_path
from primer_db import add_primer, add_primers, db, remove_primers, update_primer
from sample_files import DNA_FILE as PLASMID_DNA, PLASMID as PLASMID_GB

DEMOS = sorted((ROOT / "demo").glob("demo_*.gb"))   # the made-up demo files
DEMO = ROOT / "demo" / "demo_vector.gb"
rng = random.Random(8)
SEQ = "".join(rng.choice("ACGT") for _ in range(600))


@pytest.fixture(autouse=True)
def private_primer_db(tmp_path, monkeypatch):
    monkeypatch.setattr(primer_db, "_db", primer_db.PrimerDB(tmp_path / "primers.tsv"))


def editor(seq=SEQ, circular=False, features=()):
    doc = Document.new_sequence("t", seq, circular)
    if features:
        doc.add_features(list(features))
    return Editor(doc, Path("t.gb"))


def rich_document():
    doc = Document.new_sequence("rich", "".join(random.Random(3).choice("ACGT") for _ in range(400)), True)
    doc.add_features([("CDS", "join(10..20,30..40)", [("label", "joined")]),
                      ("misc_feature", "<100..>120", [("label", "fuzzy")])])
    return doc


# ---- saving never damages a file --------------------------------------------------------------------------------------------

def test_saving_never_damages_the_original(tmp_path):
    """A failed GenBank save (bad primer site, missing folder, a folder in the way) leaves the original file, the document's
    modified flag and no temporary file behind, and the real save afterwards works; a failed FASTA write leaves the source alone;
    a failed save in the editor keeps the path and the flag; the final rename waits for a reader to let go (Windows) but not when
    waiting cannot help; the autosave copy goes beside the file and does not count as a save."""
    # --- a failed GenBank save
    out = tmp_path / "keep.gb"
    doc = rich_document()
    doc.save_genbank(out)
    original = out.read_bytes()
    doc.insert(5, "GATTACA")                                            # unsaved edit
    with pytest.raises(OSError):                                        # a primer site that does not fit: refused before writing
        doc.save_genbank(out, [("primer_bind", "1..99999", [("label", "bad")])])
    assert out.read_bytes() == original and not list(tmp_path.glob("*.tmp"))
    assert doc.modified, "a failed save does not mark the document saved"
    with pytest.raises(OSError):                                        # a folder that does not exist: nothing created
        doc.save_genbank(tmp_path / "no_such_folder" / "x.gb")
    assert not (tmp_path / "no_such_folder").exists() and doc.modified
    folder = tmp_path / "taken.gb"                                      # a name taken by a folder: the folder is untouched
    folder.mkdir()
    (folder / "inside.txt").write_text("x")
    with pytest.raises(OSError):
        doc.save_genbank(folder)
    assert folder.is_dir() and (folder / "inside.txt").read_text() == "x" and doc.modified
    assert not list(tmp_path.glob("*.tmp")), "no half-finished file is left beside it"
    doc.save_genbank(out)                                               # the real save afterwards works
    assert out.read_bytes() != original and not doc.modified and Document(out).sequence() == doc.sequence()

    # --- a failed FASTA write leaves the source alone
    text = ">chr1 first\nACGTACGTAC\nGTAC\n>chr2\nTTTTTTTTTT\nTT\n"
    src = tmp_path / "g.fa"
    src.write_text(text)
    with pytest.raises(OSError, match="no sequence named"):
        Document.new_sequence("not_in_the_file", "ACGT", False).write_fasta(src, src)
    assert src.read_text() == text and not list(tmp_path.glob("*.tmp"))
    ok = Document.from_fasta(src, "chr1")
    with pytest.raises(OSError):
        ok.write_fasta(src, tmp_path / "no_such_folder" / "x.fa")
    blocked = tmp_path / "blocked.fa"
    blocked.mkdir()
    with pytest.raises(OSError):
        ok.write_fasta(src, blocked)
    assert src.read_text() == text and not list(tmp_path.glob("*.tmp"))

    # --- a failed save in the editor keeps the path and the modified flag
    ed = Editor(Document.new_sequence("ed", "ACGT" * 10, False), tmp_path / "ed.gb")
    ed.set_cursor(0)
    assert ed.paste("GATTACA") and ed.doc.modified
    with pytest.raises(OSError):
        ed.save(tmp_path / "no_such_folder" / "ed.gb")
    assert ed.path == tmp_path / "ed.gb" and ed.doc.modified
    ed.save(tmp_path / "ed.gb")
    assert not ed.doc.modified and Document(tmp_path / "ed.gb").sequence() == "GATTACA" + "ACGT" * 10

    # --- the rename over a file another reader holds waits for it (Windows refuses it at once)
    target, tmp = tmp_path / "a.gb", tmp_path / "a.gb.tmp"
    target.write_text("old")
    tmp.write_text("new")
    reader = open(target, "rb")   # on Windows this blocks the rename (opened without delete sharing)
    threading.Timer(0.5, reader.close).start()
    start = time.time()
    fileio.replace(tmp, target)
    assert target.read_text() == "new" and not tmp.exists()
    assert os.name != "nt" or time.time() - start >= 0.4, "it did wait on Windows"
    # ... but a folder in the way is an error at once, not after 5 s
    tmp.write_text("new")
    start = time.time()
    with pytest.raises(OSError):
        fileio.replace(tmp, folder)
    assert time.time() - start < 2

    # --- the autosave copy: beside the file (one per FASTA sequence), never the original, not a save
    assert autosave_path(tmp_path / "plasmid.dna").name == "plasmid.dna.autosave.gb"
    assert autosave_path(tmp_path / "genome.fa", "chr2").name == "genome.fa.chr2.autosave.gb"
    saved = out.read_bytes()
    doc.insert(0, "CCC")
    doc.write_genbank(str(autosave_path(out)))                          # what autosave writes
    assert out.read_bytes() == saved and doc.modified, "the user's file is only written by Save"
    assert Document(autosave_path(out)).sequence() == doc.sequence()


# ---- corrupt and truncated files -------------------------------------------------------------------------------------------

def test_corrupt_and_truncated_files_are_a_clear_error(tmp_path):
    """A truncated .dna or GenBank file is an OSError with a message or the whole record, never a crash or a silent partial
    document; corrupt and foreign files (a damaged .dna too) are refused with a message; a GenBank file of two records is refused
    as one document instead of silently dropping one, saying what to do."""
    # --- truncated at every point near the end and across the file
    for label, source in (("dna", PLASMID_DNA), ("genbank", DEMOS[0] if DEMOS else PLASMID_GB)):
        if not source.exists():
            continue
        data = source.read_bytes()
        full = Document(source)
        path = tmp_path / ("cut" + source.suffix)
        errors = 0
        for cut in sorted({*range(1, 40), *range(40, len(data), max(1, len(data) // 120))}):
            path.write_bytes(data[:len(data) - cut])
            try:
                doc = Document(path)
            except OSError as e:                                            # the only kind of failure: a normal Python error
                errors += 1
                assert str(e), "with a message"
            else:                                                           # opened: then it must be the whole record
                assert doc.sequence() == full.sequence() and len(doc.features()) == len(full.features()), (label, cut)
        assert errors > 100, label

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
    if PLASMID_DNA.exists():                                            # damage in the middle of a real .dna file
        data = bytearray(PLASMID_DNA.read_bytes())
        for i in range(100, len(data), 53):
            data[i] = 0xFF
        bad = tmp_path / "scrambled.dna"
        bad.write_bytes(bytes(data))
        with pytest.raises(OSError):
            Document(bad)

    # --- two GenBank records in one file are not one document
    if DEMOS:
        one = DEMOS[0].read_text()
        path = tmp_path / "two.gb"
        path.write_text(one + one)
        with pytest.raises(OSError, match="2 GenBank records") as err:
            Document(path)
        assert "reference genome" in str(err.value), "the message says what such a file can be used for"


# ---- editing -------------------------------------------------------------------------------------------------------------

def test_random_edits_keep_bases_and_features_then_undo_redo_exactly():
    """500 random edits (insert, delete, replace, overwrite, paste with features, feature add / delete, reverse complement; rotate
    and edits across the origin on plasmids): after each one the bases equal a plain string and every feature covers exactly the
    bases the per-base identity model says. Then undo everything and redo everything, snapshot by snapshot. Three seeds, linear
    and circular (app/selftest.py holds the machinery; the app runs short versions)."""
    for circular in (False, True):
        for seed in range(3):
            assert selftest.run_random_edits(seed, circular, 500) == 500, (seed, circular)


def test_edited_base_ranges_follow_edits_undo_redo_and_saves(tmp_path):
    """The bases edited since the last save: nothing when opened; a typed run is one stretch, replace and paste mark their bases,
    a deletion marks nothing; later edits shift or shrink the ranges; undo / redo keep them consistent; a save clears them, an
    autosave does not; an edit across the circular origin survives undo, redo and rotation."""
    # --- opened, typed, replaced, pasted
    base = "".join(random.Random(5).choice("ACGT") for _ in range(200))
    ed = editor(base)
    assert ed.doc.edited_ranges() == []
    ed.set_cursor(10)
    for b in "ACG":
        ed.type_bases(b)
    assert ed.doc.edited_ranges() == [(10, 13)], "typing a run is one stretch"
    ed.select(50, 54)
    ed.type_bases("TT")           # replace four bases by two
    assert ed.doc.edited_ranges() == [(10, 13), (50, 52)]
    ed.select(100, 100)
    ed.paste("GGGGG")
    assert ed.doc.edited_ranges()[-1] == (100, 105)

    # --- later edits shift and shrink the ranges
    ed = editor(base)
    ed.set_cursor(100)
    ed.type_bases("AAAA")
    ed.select(0, 10)
    ed.backspace()                # a deletion before: the stretch moves left, nothing new is marked
    assert ed.doc.edited_ranges() == [(90, 94)]
    ed.set_cursor(0)
    ed.type_bases("CC")           # an insertion before: it moves right
    assert ed.doc.edited_ranges() == [(0, 2), (92, 96)]
    ed.select(93, 95)
    ed.backspace()                # a deletion inside shrinks it
    assert ed.doc.edited_ranges() == [(0, 2), (92, 94)]

    # --- undo and redo
    ed = editor(base)
    ed.set_cursor(20)
    ed.type_bases("GG")
    ed.select(5, 8)
    ed.backspace()
    assert ed.doc.edited_ranges() == [(17, 19)]
    ed.undo()                     # the deleted bases come back: the stretch is where it was
    assert ed.doc.edited_ranges() == [(20, 22)]
    ed.undo()                     # undoing an insertion removes its highlight
    assert ed.doc.edited_ranges() == []
    ed.redo()
    assert ed.doc.edited_ranges() == [(20, 22)]
    ed.redo()
    assert ed.doc.edited_ranges() == [(17, 19)]

    # --- a save clears them, an autosave does not
    ed = editor(base)
    ed.set_cursor(30)
    ed.type_bases("TTT")
    ed.doc.write_genbank(str(tmp_path / "auto.gb"))        # what autosave and the project writer use
    assert ed.doc.edited_ranges() == [(30, 33)], "autosave must not count as a save"
    ed.save(tmp_path / "t.gb")
    assert ed.doc.edited_ranges() == []
    ed.set_cursor(0)
    ed.type_bases("A")
    assert ed.doc.edited_ranges() == [(0, 1)]
    ed.undo()
    ed.undo()                     # back before the save: the removed bases are gone, nothing to mark
    assert ed.doc.edited_ranges() == []
    ed.redo()
    assert ed.doc.edited_ranges() == [], "redo back onto the saved state: nothing differs from the file"
    ed.doc.modified = False       # FASTA and project saves mark the document saved this way
    assert ed.doc.edited_ranges() == []

    # --- across the circular origin
    ed = editor(base, circular=True)
    ed.select(190, 205)           # 190..199 and 0..4: passes the origin
    ed.type_bases("GGGGGG")
    ranges = ed.doc.edited_ranges()
    assert sum(b - a for a, b in ranges) == 6
    ed.undo()
    assert ed.doc.edited_ranges() == []
    ed.redo()
    assert ed.doc.edited_ranges() == ranges
    ed.set_cursor(0)
    ed.type_bases("A")
    ed.doc.rotate(50)             # rotation keeps them too
    assert sum(b - a for a, b in ed.doc.edited_ranges()) == 7
    ed.doc.undo()
    ed.doc.undo()
    ed.doc.undo()
    assert ed.doc.edited_ranges() == []

    # --- the same across the origin, with the bases as the oracle for the positions. (`base` above is all G: it builds a
    # new Random(5) for every base. Here no G, so the typed G's can be found in the string.)
    pick = random.Random(5)
    circ = "".join(pick.choice("ACT") for _ in range(200))
    ed = editor(circ, circular=True)
    ed.select(190, 205)           # removes 190..199 and 0..4; the G's go where the selection ended, after base[5:190]
    ed.type_bases("GGGGGG")
    model = circ[5:190] + "GGGGGG"
    assert ed.doc.sequence() == model and model.index("G") == 185
    assert ed.doc.edited_ranges() == [(185, 191)]
    assert "".join(ed.doc.sequence()[a:b] for a, b in ed.doc.edited_ranges()) == "GGGGGG"
    ed.set_cursor(0)
    ed.type_bases("A")
    marked = {0} | set(range(186, 192))
    model = "A" + model
    ed.doc.rotate(50)             # every marked base moves 50 back, around the circle
    model = model[50:] + model[:50]
    assert ed.doc.sequence() == model
    now = {p for a, b in ed.doc.edited_ranges() for p in range(a, b)}
    assert now == {(p - 50) % len(model) for p in marked}
    assert sorted(model[p] for p in now) == sorted("AGGGGGG")


def texts(ed):
    return [text for _, _, text in group(ed.doc.history()[0])]


def names():
    return [r["name"] for r in db().rows]


def test_undo_history_steps(monkeypatch, tmp_path):
    """The undo history: annotation changes are named steps; primer database changes are steps that undo and redo in order with
    the base edits without making the document modified (edit, rename, delete, import; a primer changed or taken meanwhile is
    left alone, with a message); the DNA lock refuses base edits and their undo but not feature edits; a file of several records
    keeps one log across them, in the order of the edits, following undo and redo."""
    FWD, REV = SEQ[100:124], SEQ[300:324]

    def fresh_db():
        monkeypatch.setattr(primer_db, "_db", primer_db.PrimerDB(tmp_path / f"primers_{time.perf_counter_ns()}.tsv"))

    # --- annotations
    ed = editor()
    ed.add_features([("misc_feature", "11..40", [("label", "alpha")])])
    index = ed.features[0].index
    ed.set_feature(index, "misc_feature", "11..60", [("label", "alpha")])
    ed.set_feature(index, "misc_feature", "11..60", [("label", "beta")])
    ed.delete_feature(index)
    assert texts(ed) == ["Added feature alpha", "Changed feature alpha", "Added feature beta; Deleted feature alpha", "Deleted feature beta"]

    # --- primer changes undo and redo with the edits
    fresh_db()
    ed = editor()
    ed.doc.save_genbank(tmp_path / "t.gb")   # saved: not modified
    assert add_primer(ed, "pFwd", FWD) == "added" and names() == ["pFwd"]
    assert not ed.doc.modified, "a primer database change does not make the document modified"
    assert [s.label for s in ed.primer_sites] == ["pFwd"]
    ed.select(10, 10)
    ed.type_bases("A")   # a base edit between the two primer steps
    assert add_primer(ed, "pRev", REV) == "added"
    assert texts(ed) == ["Added primer pFwd", "Inserted 1 bp at 11: A", "Added primer pRev"]
    assert ed.undo() and names() == ["pFwd"] and ed.doc.modified
    assert ed.undo() and ed.doc.sequence() == SEQ and not ed.doc.modified and names() == ["pFwd"]
    assert ed.undo() and names() == [] and not ed.primer_sites, "the primer is no longer matched on the sequence"
    assert not ed.undo()
    assert ed.redo() and names() == ["pFwd"] and ed.redo() and ed.doc.modified and ed.redo() and names() == ["pFwd", "pRev"]
    assert {s.label for s in ed.primer_sites} == {"pFwd", "pRev"}
    assert ed.undo() and ed.undo()
    ed.select(0, 0)
    ed.type_bases("C")
    assert not ed.redo() and names() == ["pFwd"], "a new edit drops what could have been redone, primer steps included"

    # --- editing and deleting a primer
    fresh_db()
    ed = editor()
    add_primer(ed, "p1", FWD, "first")
    add_primer(ed, "p2", REV)
    row_before = dict(db().find("p1"))
    assert update_primer(ed, "p1", "p1b", FWD, "renamed") and names() == ["p1b", "p2"]
    assert texts(ed)[-1] == "Renamed primer p1 to p1b"
    assert update_primer(ed, "p1b", "p1b", FWD, "edited") and texts(ed)[-1] == "Edited primer p1b"
    assert not update_primer(ed, "p1b", "p2", FWD, ""), "a name that belongs to another primer is refused"
    assert len(texts(ed)) == 4, "and is not a history step"
    remove_primers(ed, ["p1b", "p2"])
    assert names() == [] and texts(ed)[-1] == "Deleted primers p1b, p2"
    assert ed.undo() and names() == ["p1b", "p2"] and db().find("p1b")["note"] == "edited"
    assert ed.undo() and db().find("p1b")["note"] == "renamed"
    assert ed.undo() and names() == ["p1", "p2"] and db().find("p1") == row_before, "back to the original, date included"
    assert ed.redo() and ed.redo() and ed.redo() and names() == []
    assert ed.undo_steps(3) == 3 and names() == ["p1", "p2"] and ed.undo_steps(-3) == 3 and names() == [], "several steps at once"

    # --- importing primers; a primer changed or taken meanwhile is left alone
    fresh_db()
    ed = editor()
    report = add_primers(ed, [("a1", FWD), ("a2", REV), ("a3", SEQ[400:424])], "Imported")
    assert report["added"] == 3 and texts(ed) == ["Imported 3 primers: a1, a2, a3"]
    assert add_primers(ed, [("a1", FWD)])["exists"] == 1 and len(texts(ed)) == 1, "nothing added, nothing recorded"
    errors = []
    ed.error.connect(errors.append)
    update_primer(None, "a2", "a2", REV, "changed elsewhere")   # another document changes a2 meanwhile
    assert ed.undo() and names() == ["a2"] and db().find("a2")["note"] == "changed elsewhere", "a1 and a3 go, a2 stays as it is"
    assert errors and "a2" in errors[-1] and "changed since" in errors[-1]
    assert ed.redo() and names() == ["a1", "a2", "a3"] and "a2" in errors[-1], "redo puts a1 and a3 back; a2 is there already"
    remove_primers(ed, ["a1"])
    add_primer(None, "a1", SEQ[500:524])
    assert ed.undo() and db().find("a1")["sequence"] == SEQ[500:524] and "a1" in errors[-1] and "used now" in errors[-1]

    # --- the DNA lock
    ed = editor()
    messages = []
    ed.error.connect(messages.append)
    ed.set_cursor(5)
    assert ed.type_bases("ACGT")       # unlocked: a base edit
    ed.dna_locked = True
    assert not ed.type_bases("A") and not ed.backspace() and not ed.paste("ACGT")
    assert messages and "locked" in messages[0]
    ed.add_features([("misc_feature", "1..20", [("label", "ok")])])   # features still change
    assert ed.undo(), "the feature: allowed"
    assert not ed.undo() and ed.doc.sequence() != SEQ, "the typed bases: refused while locked"
    ed.dna_locked = False
    assert ed.undo() and ed.doc.sequence() == SEQ

    # --- one log across the records of a file
    h = FileHistory()
    h.observe("A", 0, 0)
    h.observe("B", 0, 0)                 # starting points
    h.observe("A", 1, 0)
    h.observe("B", 1, 0)
    h.observe("A", 2, 0)
    assert h.undo_log == ["A", "B", "A"] and h.next_undo() == "A"
    h.observe("A", 1, 1)                 # undo in A
    assert h.undo_log == ["A", "B"] and h.redo_log == ["A"] and h.next_redo() == "A" and h.next_undo() == "B"
    h.observe("B", 0, 1)                 # undo in B
    assert h.undo_log == ["A"] and h.redo_log == ["A", "B"] and h.next_redo() == "B"
    h.observe("B", 1, 0)                 # redo in B
    assert h.undo_log == ["A", "B"] and h.redo_log == ["A"]
    h.observe("A", 2, 0)                 # a new edit in A: A's undone step cannot come back
    assert h.undo_log == ["A", "B", "A"] and h.redo_log == []
    h.reset()
    assert h.undo_log == h.redo_log == [] and h.next_undo() is None
    h.undo_log, h.redo_log = ["A", "B", "A"], ["B", "A"]   # each record's own summaries, in the order of the edits
    undo, redo = h.steps({"A": (["a1", "a2"], ["a3"]), "B": (["b1"], ["b2"])})
    assert undo == [("A", "a1"), ("B", "b1"), ("A", "a2")]
    assert redo == [("A", "a3"), ("B", "b2")], "redo order: the next one first"


# ---- copy ----------------------------------------------------------------------------------------------------------------

def features_by_label(ed):
    return {next(v for k, v in f.qualifiers if k == "label"): f for f in ed.doc.copy_features(0, ed.length) if f.kind != "source"}


def test_copy_reverse_complement_keeps_annotations_flipped():
    """Copy Reverse Complement keeps the annotations of the selection flipped (a join keeps both parts), a plain copy keeps them
    unflipped; it works across the origin, and the reverse complement of the reverse complement is the plain copy."""
    seq = "".join(random.Random(9).choice("ACGT") for _ in range(300))
    feats = [("gene", "61..100", [("label", "top")]), ("gene", "complement(111..140)", [("label", "bottom")]),
             ("CDS", "join(51..60,121..130)", [("label", "split")])]
    # --- reverse complement: flipped, pasted with the bases
    ed = editor(seq, features=feats)
    ed.select(40, 160)
    mime = ed.copy_reverse_complement_mime()
    assert mime.text() == reverse_complement(seq[40:160])
    target = editor("A" * 50)
    target.set_cursor(10)
    assert target.paste_mime(mime)
    assert target.doc.sequence()[10:130] == reverse_complement(seq[40:160])
    pasted = features_by_label(target)
    assert set(pasted) == {"top", "bottom", "split"}
    for label, strand in (("top", -1), ("bottom", 1), ("split", -1)):
        assert pasted[label].strand == strand, f"{label}: the strand is flipped"
    # "top" covered 61..100 of the original: in the reversed copy of 41..160, pasted at 10
    assert pasted["top"].location == "complement(%d..%d)" % (10 + 160 - 100 + 1, 10 + 160 - 61 + 1)
    assert pasted["split"].location.startswith("complement(join("), "a join keeps both parts, in the new reading order"
    # Exact locations by hand: a base at 0-based p of the copied 40..159 lands at 10 + 159 - p. bottom [110,140) -> [30,60);
    # split's parts [50,60) -> [110,120) and [120,130) -> [40,50), listed in the new reading order.
    assert pasted["bottom"].location == "31..60"
    assert pasted["split"].location == "complement(join(41..50,111..120))"

    # --- the plain copy keeps them unflipped
    ed = editor(seq, features=feats)
    ed.select(40, 160)
    target = editor("A" * 50)
    target.set_cursor(0)
    assert target.paste_mime(ed.copy_mime())
    pasted = features_by_label(target)
    assert pasted["top"].strand == 1 and pasted["bottom"].strand == -1

    # --- across the origin, and twice is the original
    ed = editor(seq, circular=True, features=feats)
    ed.select(280, 320)   # 280..299 and 0..19
    mime = ed.copy_reverse_complement_mime()
    assert mime.text() == reverse_complement(seq[280:] + seq[:20])
    middle = editor("A" * 10)
    middle.set_cursor(0)
    middle.paste_mime(mime)
    middle.select(0, 40)
    assert middle.copy_reverse_complement_mime().text() == seq[280:] + seq[:20]

    # --- a feature across the origin goes along, at the place the bases say (`seq` above is one letter repeated: it builds a
    # new Random(9) for every base; this sequence is varied so the bases can tell positions apart)
    pick = random.Random(9)
    varied = "".join(pick.choice("ACGT") for _ in range(300))
    ed = editor(varied, circular=True, features=[("misc_feature", "join(291..300,1..10)", [("label", "ori")])])
    ed.select(280, 320)   # copied: varied[280:300] + varied[0:20]; ori is copy[10:30], which the reversal keeps at [10,30)
    target = editor("A" * 10)
    target.set_cursor(0)
    assert target.paste_mime(ed.copy_reverse_complement_mime())
    ori = features_by_label(target)["ori"]
    assert ori.location == "complement(11..30)"
    assert target.doc.sequence()[10:30] == reverse_complement(varied[290:300] + varied[:10])


# ---- exports and the project file ------------------------------------------------------------------------------------------

def test_text_exports_and_the_project_file(tmp_path):
    """The text exports of a document (GFF3, BED, NCBI feature table with minus-strand parts as start > stop, proteins of a CDS
    on each strand, FASTA wrapped at 70, the submission pair); an Open DNA Project (.odp) round-trips the sequence, the features,
    the tracks (also by their path beside the project), the display settings and the DNA lock."""
    ed = Editor(Document(str(DEMO)), DEMO)
    # --- GFF3, BED, feature table, proteins
    gff = exports.gff3_text(ed)
    assert gff.startswith("##gff-version 3") and gff.count("\n") >= len(ed.features) + 2
    cds = exports.cds_features(ed)
    assert cds, "the demo vector has a CDS"
    f = cds[0]
    line = next(l for l in gff.splitlines() if l.split("\t")[2:3] == ["CDS"])
    assert line.split("\t")[6] == ("+" if f.strand > 0 else "-")
    assert len(exports.bed_text(ed).splitlines()) == sum(len(f.spans) for f in ed.features)
    tbl = exports.tbl_text(ed).splitlines()
    assert tbl[0] == f">Feature {exports.seq_id(ed)}"
    first = next(l for l in tbl if l.endswith("\tCDS"))
    a, b = (int(x.strip("<>")) for x in first.split("\t")[:2])
    assert (a < b) == (f.strand > 0), "the table lists minus-strand parts as start > stop"
    assert not any("ApEinfo" in l or "\tlabel\t" in l for l in tbl)
    assert exports.proteins_fasta(ed).startswith(">")
    gene = "ATGAAACCCGGGTAA"   # a real CDS on each strand (the demo files are random sequences)
    built = Editor(Document.new_sequence("t", "TT" + gene + "TT" + reverse_complement(gene) + "TT", False), "t.gb")
    built.add_features([("CDS", "3..17", [("label", "fwd")]), ("CDS", "complement(20..34)", [("label", "rev")])])
    assert [exports.protein_of(built, f) for f in exports.cds_features(built)] == ["MKPG", "MKPG"]
    lines = exports.tbl_text(built).splitlines()
    assert "3\t17\tCDS" in lines and "34\t20\tCDS" in lines and "\t\t\tnote\trev" in lines
    # GFF3 (1-based, end inclusive) and BED (0-based start, end exclusive) coordinates by hand for CDS 3..17 and
    # complement(20..34), and the bases at each GFF3 line are the CDS on its strand
    sid = exports.seq_id(built)
    body = [l for l in exports.gff3_text(built).splitlines() if not l.startswith("#")]
    assert body == [f"{sid}\topenDNA\tCDS\t3\t17\t.\t+\t.\tID=f0;Name=fwd", f"{sid}\topenDNA\tCDS\t20\t34\t.\t-\t.\tID=f1;Name=rev"]
    assert exports.bed_text(built).splitlines() == [f"{sid}\t2\t17\tfwd\t0\t+", f"{sid}\t19\t34\trev\t0\t-"]
    whole = built.doc.sequence()
    for l in body:
        _, _, _, start, end, _, strand, _, _ = l.split("\t")
        got = whole[int(start) - 1:int(end)]
        assert (got if strand == "+" else reverse_complement(got)) == gene, l
    assert exports.translate("ATGAAATAA") == "MK*"
    assert exports.fasta_text("x", "A" * 75).splitlines()[1:] == ["A" * 70, "A" * 5]
    # --- the submission pair
    fsa, tbl_file = exports.write_submission(ed, tmp_path / "sub.tbl")
    assert fsa.read_text().startswith(f">{exports.seq_id(ed)}") and tbl_file.exists()

    # --- the project file
    track = tmp_path / "extra.gff3"
    track.write_text(f"{ed.doc.name}\tx\tgene\t10\t50\t.\t+\t.\tName=g1\n")
    ed.tracks.append(SimpleNamespace(name=track.name, path=str(track), seqid=ed.doc.name))   # as the window's load_track does
    ed.set_display(enzymes=True, hidden_kinds={"misc_feature"})
    ed.dna_locked = True
    site = ed.doc.sequence()[100:120]
    assert ed.doc.sequence().upper().count(site.upper()) == 1 and reverse_complement(site).upper() not in ed.doc.sequence().upper()
    assert ed.import_primers([("pX", site)]) == 1     # a primer of this file: written as a primer_bind feature
    project = tmp_path / "work.odp"
    exports.save_project(ed, project)
    with zipfile.ZipFile(project) as z:
        assert set(z.namelist()) == {"sequence.gb", "project.json"}
        meta = json.loads(z.read("project.json"))
    assert meta["version"] == exports.PROJECT_VERSION and meta["dna_locked"]
    assert meta["tracks"][0]["relative"] == "extra.gff3" and meta["tracks"][0]["path"] == str(track)
    assert meta["display"]["enzymes"] and meta["display"]["hidden_kinds"] == ["misc_feature"]
    assert not list(tmp_path.glob("*.part")), "written through a temporary file"
    doc, back = exports.open_project(project)
    assert doc.sequence() == ed.doc.sequence() and back == meta
    # every feature comes back whole (kind, location, qualifiers), plus the file's primer at the site its bases say
    want = [(f.kind, f.location, list(f.qualifiers)) for f in ed.doc.features()]
    want.append(("primer_bind", "101..120", [("label", "pX"), ("note", f"sequence: {site}")]))
    assert sorted((f.kind, f.location, list(f.qualifiers)) for f in doc.features()) == sorted(want)
    moved = tmp_path / "moved"
    moved.mkdir()
    (moved / "extra.gff3").write_text(track.read_text())
    entry = dict(meta["tracks"][0], path=str(tmp_path / "gone.gff3"))
    assert exports.project_file(entry, moved / "work.odp") == str(moved / "extra.gff3"), "found beside the moved project"
    bad = tmp_path / "bad.odp"
    bad.write_bytes(b"not a zip")
    with pytest.raises(OSError, match="not an Open DNA Project"):
        exports.open_project(bad)


# ---- gene groups ---------------------------------------------------------------------------------------------------------

def grouped_genbank(path):
    """1,200 bp with six genes (each a gene + CDS sharing a /locus_tag), a misc_feature in gene F and one in no gene."""
    feats = []
    for i, name in enumerate("ABCDEF"):
        s = 1 + i * 180
        loc = f"{s}..{s + 149}" if i % 2 == 0 else f"complement({s}..{s + 149})"
        feats.append(f"     gene            {loc}\n                     /locus_tag=\"T_{name}\"\n                     /gene=\"gen{name}\"")
        feats.append(f"     CDS             {loc}\n                     /locus_tag=\"T_{name}\"\n                     /label=cds{name}")
    feats.append("     misc_feature    901..920\n                     /locus_tag=\"T_F\"\n                     /label=tagF")
    feats.append("     misc_feature    1101..1150\n                     /label=lonely")
    seq = ("ATGAAACGCATTAGCACCACC" * 60)[:1200]
    origin = "\n".join(f"{i + 1:>9} " + " ".join(seq[j:j + 10] for j in range(i, min(i + 60, 1200), 10))
                       for i in range(0, 1200, 60))
    path.write_text("LOCUS       groups      1200 bp    DNA     linear   UNK 01-JAN-1980\nFEATURES             Location/Qualifiers\n"
                    + "\n".join(feats) + "\nORIGIN\n" + origin + "\n//\n")
    return path


def test_gene_groups_and_display_switches(tmp_path):
    """Features sharing a /locus_tag form one gene group; a display switch (collapse genes, hide a type) lays the features out
    again without reading them all from the core again (it also dropped the Local Align hits); the count of each type follows
    the edits."""
    ed = Editor(Document(str(grouped_genbank(tmp_path / "g.gb"))), tmp_path / "g.gb")
    assert len(ed.groups) == 6 and len(ed.groups["T_F"].parts) == 2, "gene F: its CDS and the tagged misc_feature"
    assert ed.kind_counts == {"gene": 6, "CDS": 6, "misc_feature": 2}
    ed.query_hits = ["a hit"]
    reads = []
    real = ed.doc.features
    ed.doc = type("Spy", (), {"__getattr__": lambda s, k: getattr(real.__self__, k),
                              "features": lambda s: reads.append(1) or real()})()
    ed.set_display(hidden_kinds={"misc_feature"})
    assert not reads and ed.query_hits == ["a hit"]
    assert not any(f.kind == "misc_feature" for f in ed.shown)
    ed.set_display(collapse_genes=not ed.collapse_genes)
    ed.set_display(hidden_kinds=())
    assert not reads and any(f.kind == "misc_feature" for f in ed.shown)
    ed.doc = real.__self__
    ed.add_features([("repeat_region", "1001..1010", [("label", "r")])])
    assert ed.kind_counts["repeat_region"] == 1 and ed.kind_counts["gene"] == 6
