"""Two pieces of app logic the editor owes the rest of the program: "CDS only" (a gene at exactly a CDS's bases is hidden, and
features still cover the same bases after edits) and the History texts written for every kind of edit."""
import os
import random
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ["OPENDNA_DATA_DIR"] = tempfile.mkdtemp(prefix="opendna_test_")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "app"))

from opendna_core import Document

from editor import Editor
from history_panel import group


def random_dna(n, seed):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def editor_of(seq, circular=False, features=(), name="t"):
    doc = Document.new_sequence(name, seq, circular)
    if features:
        doc.add_features(list(features))
    return Editor(doc, Path(f"{name}.gb"))


def covered(ed, f):
    """The bases a feature covers, in order."""
    return "".join(ed.doc.sequence(s, e) for s, e in f.spans)


def by_label(ed, label):
    return next(f for f in ed.features if f.label == label)


GENE_SEQ = random_dna(1500, 11)


def gene_features():
    q = lambda label: [("label", label)]
    return [
        ("gene", "101..400", q("gA")),                 # longer than its CDS: kept
        ("CDS", "151..390", q("cA")),
        ("promoter", "61..100", q("pA")),
        ("gene", "900..1000", q("abc")),               # exactly its CDS's bases and strand: hidden
        ("CDS", "900..1000", q("abcCDS")),
        ("gene", "complement(1100..1200)", q("g2")),   # same bases as a CDS on the other strand: kept
        ("CDS", "1100..1200", q("c2")),
    ]


def test_cds_only_hides_a_gene_at_exactly_a_cds_and_edits_keep_the_bases():
    """A gene feature at exactly the bases and strand of a CDS is hidden by "CDS only" (every view), the others are kept, and
    after inserts, deletes inside the gene and an undo every feature still covers exactly the bases it did."""
    # --- which genes are hidden ---
    ed = editor_of(GENE_SEQ, features=gene_features())
    shown = lambda: {f.label for f in ed.shown}
    assert ed.cds_only and shown() == {"gA", "cA", "pA", "abcCDS", "g2", "c2"}
    assert "abc" not in {f.label for f in ed.visible_features(0, ed.length)}
    ed.set_display(cds_only=False)
    assert shown() == {f.label for f in ed.features}

    # --- still the same bases after edits ---
    ed = editor_of(GENE_SEQ, features=gene_features())
    before = {f.label: covered(ed, f) for f in ed.features}
    ed.set_cursor(0)
    ed.type_bases("A" * 40)                                          # before everything: all shift
    ed.select(300, 310)
    ed.backspace()                                                   # inside the gene (and its CDS): both shrink
    after = {f.label: covered(ed, f) for f in ed.features}
    assert after["cA"] == before["cA"][:110] + before["cA"][120:]     # bases 260..270 of the sequence were the CDS's 110th..120th
    assert after["gA"] == before["gA"][:160] + before["gA"][170:]
    for label in ("pA", "abcCDS", "abc", "g2", "c2"):                # not touched by the delete: same bases
        assert after[label] == before[label], label
    # An insert inside the gene but outside its CDS grows the gene and moves the CDS: still the same bases.
    cds_text, gene_len = covered(ed, by_label(ed, "cA")), len(covered(ed, by_label(ed, "gA")))
    ed.set_cursor(40 + 125)
    ed.type_bases("CCCCC")
    assert covered(ed, by_label(ed, "cA")) == cds_text
    assert len(covered(ed, by_label(ed, "gA"))) == gene_len + 5
    # Everything undone: the original bases.
    assert ed.undo_steps(3) == 3                                  # the three edits (the features were there from the start)
    assert {f.label: covered(ed, f) for f in ed.features} == before


def test_history_texts():
    """The text History shows for each kind of edit (insert, delete, short and long replace, paste with features,
    feature added / changed / renamed / deleted, reverse complement) and for a circular sequence (rotation, linearising,
    edits across the origin)."""
    # --- each kind of edit ---
    seq = random_dna(300, 4)
    ed = editor_of(seq)
    expected = []

    def step(text):
        expected.append(text)
        assert [t for _, _, t in group(ed.doc.history()[0])] == expected

    ed.set_cursor(10)
    ed.type_bases("GGG")
    step("Inserted 3 bp at 11: GGG")
    ed.select(10, 13)
    ed.backspace()
    step("Deleted 3 bp at 11–13: " + "GGG")
    old = seq[20:23]
    new = "".join("A" if c != "A" else "C" for c in old)
    ed.select(20, 23)
    ed.type_bases(new)
    step(f"{old} > {new} at 21")                                      # short replace: both sides spelt out
    ed.select(30, 60)
    ed.type_bases("AAAAA")
    step("Replaced 30 bp with 5 bp at 31")                           # longer than 12: counts only
    ed.set_cursor(0)
    ed.type_bases("C" * 40)
    step("Inserted 40 bp at 1: " + "C" * 30 + "…")                   # the first 30 bases, then an ellipsis
    gone = ed.doc.sequence(100, 140)     # the bases about to go, read from the sequence itself (not the undo record)
    assert len(gone) == 40
    ed.select(100, 140)
    ed.backspace()
    step("Deleted 40 bp at 101–140: " + gone[:30] + "…")                   # the first 30 of the 40 bases
    # Paste: the bases and the features that came with them.
    ed.set_cursor(200)
    ed.paste("GGGCCC", clip=("GGGCCC", [("misc_feature", "1..3", [("label", "pf")])]))
    step("Pasted 6 bp at 201: GGGCCC; Added feature pf")
    ed.set_cursor(250)
    ed.paste("TTTAAA", clip=("TTTAAA", [("misc_feature", "1..3", [("label", "q1")]), ("misc_feature", "4..6", [("label", "q2")])]))
    step("Pasted 6 bp at 251: TTTAAA; Added features q1, q2")
    # Features: add, change, rename (a new name is a new feature and the old one gone), delete.
    ed.add_features([("misc_feature", "10..20", [("label", "mine")])])
    step("Added feature mine")
    index = by_label(ed, "mine").index
    ed.set_feature(index, "misc_feature", "12..22", [("label", "mine")])
    step("Changed feature mine")
    ed.set_feature(index, "misc_feature", "12..22", [("label", "renamed")])
    step("Added feature renamed; Deleted feature mine")
    ed.delete_feature(index)
    step("Deleted feature renamed")
    ed.add_features([("misc_feature", f"{i * 10 + 1}..{i * 10 + 5}", [("label", f"n{i}")]) for i in range(5)])
    step("Added features n0, n1, n2, n3…")                            # at most four names, then an ellipsis
    ed.reverse_complement()
    step("Reverse-complemented the whole sequence")
    # Undone, the row keeps its text and waits to be redone.
    ed.undo()
    assert [t for _, _, t in group(ed.doc.history()[1])] == ["Reverse-complemented the whole sequence"]

    # --- a circular sequence ---
    ed = editor_of(random_dna(200, 5), circular=True)
    ed.doc.rotate(7)
    ed.linearize(0)
    ed.set_circular(True)
    ed.linearize(50)
    ed.set_circular(True)
    ed.select(190, 205)                  # across the origin
    ed.type_bases("GG")
    texts = [t for _, _, t in group(ed.doc.history()[0])]
    assert texts[0] == "Moved base 1 (broke the circle elsewhere)"
    assert texts[1] == "Made the sequence linear"                    # broken at the origin: nothing cut
    assert texts[2] == "Made the sequence circular"
    assert texts[3] == "Made the sequence linear (broke the circle elsewhere)"
    assert texts[4] == "Made the sequence circular"
    assert texts[5] == "Replaced 15 bp across the origin from 191 with GG"
    n = ed.length
    ed.select(n - 5, n + 10)
    ed.backspace()                       # deleting across the origin: from 0-based n - 5, so base n - 4
    assert group(ed.doc.history()[0])[-1][2] == f"Replaced 15 bp across the origin from {n - 4} (deleted)"

    # Set Origin to Selection Start: the circle turns so base 10 (0-based) becomes base 1; nothing happens on a line or at base 1.
    seq = ed.doc.sequence()
    assert ed.set_origin(10) and ed.doc.circular and ed.doc.sequence() == seq[10:] + seq[:10] and ed.start == 0
    assert not ed.set_origin(0)
    ed.set_circular(False)
    assert not ed.set_origin(5)
