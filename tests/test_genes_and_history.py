"""Two pieces of app logic the editor owes the rest of the program: gene groups (the features that share a tag folded
under their gene, and still covering the same bases after edits) and the History texts written for every kind of edit."""
import os
import random
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ["OPENDNA_DATA_DIR"] = tempfile.mkdtemp(prefix="opendna_test_")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "app"))

from opendna_core import Document

import genes
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
    q = lambda label, **kw: [("label", label)] + list(kw.items())
    return [
        ("gene", "101..400", q("gA", gene="gA", locus_tag="L1")),
        ("mRNA", "121..395", q("mA", locus_tag="L1")),
        ("CDS", "151..390", q("cA", locus_tag="L1")),
        ("promoter", "61..100", q("pA")),                  # near the gene but with no tag: alone
        ("misc_feature", "600..650", q("solo")),
        ("CDS", "700..760", q("orphan", locus_tag="L9")),  # a tag with no gene feature: alone
        ("gene", "900..1000", q("abc", gene="abc")),       # only /gene: grouped by it
        ("CDS", "910..990", q("abcCDS", gene="abc")),
        ("gene", "1100..1200", q("g2", locus_tag="L2")),   # a gene with no parts
    ]


def test_gene_groups_fold_their_parts_and_survive_edits():
    """Features sharing a tag fold under their gene and the others are left alone (which of them the views are given,
    and the switch that unfolds them); after inserts, deletes inside the gene and an undo, every feature still covers
    exactly the bases it did."""
    # --- fold features sharing a tag and leave the others alone ---
    ed = editor_of(GENE_SEQ, features=gene_features())
    assert set(ed.groups) == {"L1", "abc", "L2"}
    g = ed.groups["L1"]
    assert g.master.label == "gA" and [p.label for p in g.parts] == ["mA", "cA"]   # mRNA starts before the CDS: same rank, by position
    assert ed.groups["abc"].master.label == "abc" and [p.label for p in ed.groups["abc"].parts] == ["abcCDS"]
    assert ed.groups["L2"].parts == []
    grouped = {by_label(ed, x).index for x in ("mA", "cA", "abcCDS")}
    assert set(ed.part_of) == grouped and ed.part_of[by_label(ed, "cA").index] == "L1"
    for alone in ("pA", "solo", "orphan"):     # in no group: not a part, not a master
        f = by_label(ed, alone)
        assert f.index not in ed.part_of and all(x.master is not f for x in ed.groups.values())
    # The map draws the gene, not its parts; features in no group are drawn as they are. A gene whose only part is one CDS
    # ("abc") is not drawn: the CDS is, in every view (View > Hide Gene Arrow when its CDS is Alone).
    assert ed.cds_only and {f.label for f in ed.map_features} == {"gA", "pA", "solo", "orphan", "abcCDS", "g2"}
    assert "abc" not in {f.label for f in ed.shown} and "abcCDS" in {f.label for f in ed.shown}
    assert "abc" not in {f.label for f in ed.visible_features(0, ed.length)}
    ed.set_display(cds_only=False)
    assert {f.label for f in ed.map_features} == {"gA", "pA", "solo", "orphan", "abc", "g2"} and "abc" in {f.label for f in ed.shown}
    ed.set_display(cds_only=False, collapse_genes=False)
    assert {f.label for f in ed.map_features} == {f.label for f in ed.features}

    # --- still the same bases after edits ---
    ed = editor_of(GENE_SEQ, features=gene_features())
    before = {f.label: covered(ed, f) for f in ed.features}
    ed.set_cursor(0)
    ed.type_bases("A" * 40)                                          # before everything: all shift
    ed.select(300, 310)
    ed.backspace()                                                   # inside the gene (and its CDS): all shrink
    assert set(ed.groups) == {"L1", "abc", "L2"} and [p.label for p in ed.groups["L1"].parts] == ["mA", "cA"]
    after = {f.label: covered(ed, f) for f in ed.features}
    assert after["cA"] == before["cA"][:110] + before["cA"][120:]     # bases 260..270 of the sequence were the CDS's 110th..120th
    assert after["gA"] == before["gA"][:160] + before["gA"][170:] and len(after["mA"]) == len(before["mA"]) - 10
    for label in ("pA", "solo", "orphan", "abcCDS", "abc", "g2"):   # not touched by the delete: same bases
        assert after[label] == before[label], label
    # The gene still spans all its parts.
    lo, hi = genes.span_of([ed.groups["L1"].master])
    for p in ed.groups["L1"].parts:
        plo, phi = genes.span_of([p])
        assert lo <= plo and phi <= hi
    # An insert inside the gene but outside its CDS grows the gene and the mRNA, and moves the CDS: still the same bases.
    cds_text, gene_len = covered(ed, by_label(ed, "cA")), len(covered(ed, by_label(ed, "gA")))
    ed.set_cursor(40 + 125)
    ed.type_bases("CCCCC")
    assert covered(ed, by_label(ed, "cA")) == cds_text
    assert len(covered(ed, by_label(ed, "gA"))) == gene_len + 5 and ed.groups["L1"].master.label == "gA"
    # Everything undone: the original bases and groups.
    assert ed.undo_steps(3) == 3                                  # the three edits (the features were there from the start)
    assert {f.label: covered(ed, f) for f in ed.features} == before and set(ed.groups) == {"L1", "abc", "L2"}


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
    ed.select(100, 140)
    ed.backspace()
    step("Deleted 40 bp at 101–140: " + ed.doc.history()[0][-1][1] + "…")      # the first 30 of the 40 bases
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
    ed.select(ed.length - 5, ed.length + 10)
    ed.backspace()                       # deleting across the origin
    last = group(ed.doc.history()[0])[-1][2]
    assert last.startswith("Replaced 15 bp across the origin from ") and last.endswith(" (deleted)")
