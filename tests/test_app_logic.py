"""App logic that sits between the Rust core and the screen (docs/TESTING.md section 6): gene groups, the display filter,
History texts, Assembly checks (Sidewinder, bad input), Local Align results, the primer checks' thresholds, the feature table,
settings, autosave and recovery, the self-test panel. Everything is built in memory (random sequences, a private primer
database); no real file of the repo is needed."""
import os
import random
import sys
import tempfile
from collections import Counter
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ["OPENDNA_DATA_DIR"] = tempfile.mkdtemp(prefix="opendna_test_")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "app"))

from PySide6.QtCore import QSettings, Qt
from PySide6.QtWidgets import QAbstractItemView, QApplication, QMessageBox, QPushButton, QWidget

from opendna_core import Document

import assembly as A
import assembly_panel as P
import genes
import primer_db
import primer_panel as PP
import selftest
import tool_settings
from display_panel import DisplayPanel
from editor import Editor, reverse_complement
from feature_dialog import FeatureDialog
from feature_table import FeatureTable
from history_panel import NOW, HistoryPanel, group
from local_align import FORWARD_TEXT, INVERTED_TEXT, LocalAlignDialog
from main import MainWindow

app = QApplication.instance() or QApplication([])


@pytest.fixture(autouse=True)
def private_primer_db(tmp_path, monkeypatch):
    """A primer database of this test's own (the shared one may be big, and tests here add primers)."""
    monkeypatch.setattr(primer_db, "_db", primer_db.PrimerDB(tmp_path / "primers.tsv"))


@pytest.fixture
def errors(monkeypatch):
    """Qt swallows exceptions raised in paint/event handlers; collect them instead."""
    caught = []
    monkeypatch.setattr(sys, "excepthook", lambda *exc: caught.append(exc))
    yield caught
    assert not caught, caught


def random_dna(n, seed):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def editor_of(seq, circular=False, features=(), name="t"):
    doc = Document.new_sequence(name, seq, circular)
    if features:
        doc.add_features(list(features))
    return Editor(doc, Path(f"{name}.gb"))


def window(doc, path="t.gb"):
    win = MainWindow()
    win.resize(1200, 800)
    win.show()
    win.show_document(doc, Path(path))
    app.processEvents()
    return win


def close(win):
    win.editor.doc.modified = False   # no "save changes?" question
    win.close()


def covered(ed, f):
    """The bases a feature covers, in order."""
    return "".join(ed.doc.sequence(s, e) for s, e in f.spans)


def by_label(ed, label):
    return next(f for f in ed.features if f.label == label)


# ---- gene groups ----------------------------------------------------------------------------------------------------------

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


def test_gene_groups():
    """Gene groups fold the features sharing a tag and leave the others alone (table, map, display), and still cover the
    same bases after edits (shift, shrink, grow, undo)."""
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
    # The table lists a gene as one row (its parts under it once opened); features in no group are rows of their own.
    table = FeatureTable(ed)
    top = [r[3] for r in table.model.shown]
    assert sorted(top) == sorted(f.index for f in ed.features if f.index not in grouped)
    assert [r[4] for r in table.model.rows if r[3] == g.master.index] == ["L1"]
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


# ---- display filter ------------------------------------------------------------------------------------------------------


def test_display_filter_hides_exactly_the_types_switched_off(errors):
    win = window(Document.new_sequence("t", GENE_SEQ, False))
    ed = win.editor
    ed.set_display(cds_only=False)   # (the default hides a gene whose only part is one CDS)
    ed.add_features(gene_features())
    panel = DisplayPanel(ed)
    kinds = Counter(f.kind for f in ed.features)
    assert {k: int(box.text().rsplit("(", 1)[1].rstrip(")")) for k, box in panel.kinds.items()} == dict(kinds)   # the counts shown
    assert all(box.isChecked() for box in panel.kinds.values())

    panel.kinds["promoter"].setChecked(False)
    assert ed.hidden_kinds == {"promoter"}
    assert {f.index for f in ed.shown} == {f.index for f in ed.features if f.kind != "promoter"}
    assert {f.index for f in ed.visible_features(0, ed.length)} == {f.index for f in ed.shown}   # what the views draw
    assert not any(f.kind == "promoter" for f in ed.map_features) and set(ed.map_lanes) == {f.index for f in ed.map_features}
    assert len(ed.features) == sum(kinds.values())   # still in the file ...
    model = FeatureTable(ed).model
    assert len(model.rows) + sum(map(len, model.children.values())) == len(ed.features)   # ... and in the table
    win.map_view.grab()
    win.sequence_view.grab()

    # Hiding the gene type hides the genes only: their parts are then drawn on the map, not folded into a gene that is not there.
    panel.kinds["gene"].setChecked(False)
    assert ed.hidden_kinds == {"promoter", "gene"}
    drawn = {f.label for f in ed.map_features}
    assert {"mA", "cA", "abcCDS"} <= drawn and not {"gA", "abc", "g2", "pA"} & drawn
    assert drawn == {f.label for f in ed.shown}
    win.map_view.grab()
    panel.kinds["gene"].setChecked(True)
    panel.kinds["promoter"].setChecked(True)
    assert not ed.hidden_kinds and {f.label for f in ed.map_features} == {"gA", "pA", "solo", "orphan", "abc", "g2"}

    # None / All, and the filter follows an edit (a type that appears or vanishes is added to or removed from the list).
    none = next(b for b in panel.findChildren(QPushButton) if b.text() == "None")
    none.click()
    assert ed.hidden_kinds == set(kinds) and ed.shown == [] and ed.map_features == []
    next(b for b in panel.findChildren(QPushButton) if b.text() == "All").click()
    assert not ed.hidden_kinds
    panel.kinds["CDS"].setChecked(False)
    ed.set_cursor(5)
    ed.type_bases("A")           # the filter lives with the document: an edit keeps it
    ed.undo()
    assert ed.hidden_kinds == {"CDS"} and not any(f.kind == "CDS" for f in ed.shown)
    ed.add_features([("terminator", "1400..1420", [("label", "t1")])])
    assert "terminator" in panel.kinds and panel.kinds["terminator"].isChecked() and panel.kinds["CDS"].text().startswith("CDS")
    assert not panel.kinds["CDS"].isChecked()
    # Each document has its own filter (per file, while it is open).
    other = window(Document.new_sequence("o", GENE_SEQ, False), "o.gb")
    other.editor.add_features(gene_features())
    assert not other.editor.hidden_kinds and other.editor.shown
    close(other)
    close(win)


# ---- History texts ---------------------------------------------------------------------------------------------------------


def test_history_texts():
    """History texts for each kind of edit (insert, delete, replace, paste, features, reverse complement, redone rows
    after the marker) and for a circular sequence (rotation, linearising, edits across the origin)."""
    # --- each kind of edit ---
    seq = random_dna(300, 4)
    ed = editor_of(seq)
    panel = HistoryPanel(ed)
    expected = []

    def step(text):
        expected.append(text)
        assert [t for _, _, t in group(ed.doc.history()[0])] == expected
        assert [panel.list.item(i).text() for i in range(panel.list.count())] == ["Opened t.gb"] + expected + [NOW]

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
    # Redone rows are listed after the current marker, greyed; the texts are the same.
    ed.undo()
    texts = [panel.list.item(i).text() for i in range(panel.list.count())]
    assert texts[-2:] == [NOW, "Reverse-complemented the whole sequence"]

    # --- a circular sequence ---
    seq = random_dna(200, 5)
    ed = editor_of(seq, circular=True)
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


# ---- Assembly: Sidewinder and bad input --------------------------------------------------------------------------------------


def piece(name, seq, features=()):
    doc = Document.new_sequence(name, seq)
    if features:
        doc.add_features([(k, loc, [("label", lab)]) for k, loc, lab in features])
    return A.Frag(name, doc)


def sidewinder_parts(seed):
    return [piece("a", random_dna(100, seed), [("promoter", "1..50", "pro")]),
            piece("b", random_dna(120, seed + 1000), [("CDS", "10..60", "cds1")]),
            piece("c", random_dna(90, seed + 2000))]


def test_sidewinder():
    """Sidewinder: toeholds, oligos and the product with its features; the product is linear even when a circle is asked
    for (Gibson closes the loop; fewer than two fragments is an error)."""
    # --- toeholds, oligos and product ---
    for seed in range(60):          # a random draw can hold a poor toehold (the checks are right to say so): take one that is fine
        parts = sidewinder_parts(seed)
        d = A.analyse(parts, 2, False)
        if d.ok:
            break
    assert d.ok, "no seed gave a clean Sidewinder design"
    product = "".join(f.seq for f in parts)
    assert d.product() == product and len(product) == 310
    half = A.TOEHOLD // 2
    bounds = [100, 220]
    assert [j.kind for j in d.joints] == ["toehold", "toehold"] and [(j.left, j.right) for j in d.joints] == [("a", "b"), ("b", "c")]
    for j, b in zip(d.joints, bounds):
        assert len(j.seq) == A.TOEHOLD
        assert j.seq == product[j.pos:j.pos + A.TOEHOLD]             # taken from the product itself: scarless
        assert abs(j.pos + half - b) <= A.MAX_SHIFT                  # straddles the boundary, or slid a little to dodge a repeat
        assert product.count(j.seq) == 1
    # The coding strands tile the product: each ends with its toehold, the next starts right after it.
    coding = {r[0]: r[1] for r in d.rows if r[0].endswith("-coding")}
    ends = [j.pos + A.TOEHOLD for j in d.joints]
    assert [coding[n] for n in ("a-coding", "b-coding", "c-coding")] == [product[:ends[0]], product[ends[0]:ends[1]], product[ends[1]:]]
    assert "".join(coding.values()) == product
    assert coding["a-coding"].endswith(d.joints[0].seq) and coding["b-coding"].endswith(d.joints[1].seq)
    # The barcode strand is the reverse complement plus an 18 nt arm (none on the last fragment); the arms differ.
    barcode = {r[0]: r[1] for r in d.rows if r[0].endswith("-barcode")}
    arms = []
    for name, c in coding.items():
        full = barcode[name.replace("coding", "barcode")]
        assert full.startswith(reverse_complement(c))
        arms.append(full[len(c):])
    assert [len(a) for a in arms] == [A.BARCODE, A.BARCODE, 0] and len(set(arms[:2])) == 2
    assert all(len(r[1]) <= A.MAX_OLIGO for r in d.rows if r[3] == "Sidewinder oligo")
    # The product carries each fragment's features, moved to where the fragment sits, and no marker of its own...
    doc = A.simulate(parts, d)
    feats = {f.label: f for f in doc.features()}
    assert feats["pro"].spans == [(0, 50)] and feats["cds1"].spans == [(100 + 9, 100 + 60)]
    assert not {"a", "b", "c"} & set(feats)
    assert doc.sequence(0, None).upper() == product
    # ...the preview of the panel names the fragments.
    marks = {f.label: f for f in A.simulate(parts, d, preview=True).features()}
    assert [(marks[n].spans[0]) for n in ("a", "b", "c")] == [(0, 100), (100, 220), (220, 310)]

    # --- the product is linear even when a circle is asked for ---
    parts =sidewinder_parts(3)
    linear, asked_circular = A.analyse(parts, 2, False), A.analyse(sidewinder_parts(3), 2, True)
    for d in (linear, asked_circular):
        assert not d.circular and len(d.joints) == 2          # two junctions, not three: nothing joins the last to the first
        assert any(c[0] == "Linear product" for c in d.checks)
        assert not A.simulate(parts, d).circular
    assert linear.product() == asked_circular.product()
    gibson = [piece("a", random_dna(300, 1)), piece("b", random_dna(300, 2))]
    A.auto_primers(gibson, 0, True)
    assert A.analyse(gibson, 0, True).circular and len(A.analyse(gibson, 0, True).joints) == 2   # Gibson closes the loop
    with pytest.raises(ValueError, match="at least two"):
        A.analyse(parts[:1], 2, False)


def bad_checks(d):
    return [c for c in d.checks if c[1] == "bad"]


def find_text(d, name):
    return next(t for n, _, t in d.checks if n == name)


def sidewinder_arms(d):
    """The barcode arms: what each barcode strand carries past the reverse complement of its coding strand."""
    coding = {r[0]: r[1] for r in d.rows if r[0].endswith("-coding")}
    arms = [r[1][len(coding[r[0].replace("-barcode", "-coding")]):] for r in d.rows if r[0].endswith("-barcode")]
    return [a for a in arms if a]


def test_sidewinder_toeholds_and_barcodes():
    """A toehold that occurs twice in the product slides along the junction until it is unique (still scarless);
    the barcodes share no 8-mer with each other or with a toehold; the same fragments give the same oligos, and
    two products of the same length do not share barcodes."""
    # --- a repeat at the first junction: the toehold slides off it ---
    a, b = random_dna(100, 7), random_dna(120, 8)
    planted = (a + b)[95:105]                   # the toehold the boundary itself would give
    parts = [piece("a", a), piece("b", b), piece("c", random_dna(40, 9) + planted + random_dna(40, 10))]
    d = A.analyse(parts, 2, False)
    product = d.product()
    assert product.count(planted) == 2, "the repeat must really be in the product"
    assert d.joints[0].seq != planted and product.count(d.joints[0].seq) == 1
    assert abs(d.joints[0].pos + A.TOEHOLD // 2 - 100) <= A.MAX_SHIFT     # slid, not moved somewhere else entirely
    assert find_text(d, "Toehold 1 unique") == "occurs once in the product"
    coding = [r[1] for r in d.rows if r[0].endswith("-coding")]
    assert "".join(coding) == product, "sliding only moves which oligo carries the bases: the product is unchanged"

    # --- the barcodes are orthogonal to each other and to every toehold ---
    kmers = lambda s: {x[i:i + 8] for x in (s, reverse_complement(s)) for i in range(len(x) - 7)}
    arms = sidewinder_arms(d)
    assert len(arms) == 2 and all(len(x) == A.BARCODE for x in arms)
    assert not kmers(arms[0]) & kmers(arms[1])
    for arm in arms:
        assert 8 <= sum(c in "GC" for c in arm) <= 10 and not any(c * 4 in arm for c in "ACGT")
        assert all(not kmers(arm) & kmers(j.seq) for j in d.joints)

    # --- same fragments, same oligos; two products of the same length, different barcodes ---
    again = A.analyse([piece(f.name, f.seq) for f in parts], 2, False)
    assert [r[1] for r in again.rows] == [r[1] for r in d.rows], "an order must be reproducible"
    one, two = A.analyse(sidewinder_parts(11), 2, False), A.analyse(sidewinder_parts(12), 2, False)
    assert len(one.product()) == len(two.product())        # same length, different sequence
    assert sidewinder_arms(one) != sidewinder_arms(two), "barcodes are drawn from the product, not from its length"


def test_assembly_bad_input():
    """A tiny or empty fragment gives a bad check, never an exception (every method, circular or linear, 0 / 5 / 18 bp);
    an empty stretch of a source is a fragment with no bases and the panel shows the failed check; a repeated overhang is
    a bad check and blocks the assembly."""
    # --- a tiny or empty fragment gives a bad check, never an exception ---
    for method, method_name in ((0, "gibson"), (1, "golden gate"), (2, "sidewinder")):
        for circular in (True, False):
            for size in (0, 5, 18):
                case = f"{method_name}, {'circular' if circular else 'linear'}, {size} bp"
                seqs = [random_dna(300, 1), random_dna(size, 2), random_dna(250, 3)]
                parts = [piece(n, s) for n, s in zip("ABC", seqs)]
                d = A.analyse(parts, method, circular)                  # as it is: its ends cannot meet its neighbours'
                assert not d.ok and bad_checks(d) and all(text for _, _, text in bad_checks(d)), case
                A.auto_primers(parts, method, circular)                 # designing primers for it must not fail either
                d = A.analyse(parts, method, circular)
                if method == 2 or size == 0:
                    assert not d.ok and bad_checks(d), case             # Sidewinder needs two toeholds in it; nothing joins an empty one
                else:                                                   # tails can give a short one its overlaps; what is left is the length warning
                    assert dict((n, lv) for n, lv, _ in d.checks)["Fragment length"] == "warn" and "B" in find_text(d, "Fragment length"), case
                if size == 0:
                    assert parts[1].doc.length == 0 and parts[1].seq == "", case

    # --- an empty stretch of a source is a fragment with no bases ---
    src =Document.new_sequence("s", random_dna(100, 6))
    for start, end in ((5, 5), (0, 0), (50, 20)):
        f = A.Frag("e", src, start, end)
        assert f.seq == "" and f.core_seq == "" and f.tails == (0, 0)
    # ... and the Assembly panel shows the failed check instead of raising.
    class Win(QWidget):
        def open_documents(self):
            return []
    panel = P.AssemblyPanel(Win())
    for i, n in enumerate((300, 0, 250)):
        panel._set_source(i, Document.new_sequence(f"s{i}", random_dna(n, i)), f"s{i}.gb", None)
    assert panel.design is not None and not panel.design.ok and bad_checks(panel.design)
    assert panel.status.text().startswith("✗") and not panel.buttons["perform"].isEnabled()

    # --- a repeated overhang is a bad check and blocks the assembly ---
    def gg(name, left, core, right):
        return piece(name, "tt" + "GGTCTC" + "a" + left + core + right + "t" + "GAGACC" + "tt")
    def cores():
        rng = random.Random(8)
        out = []
        while len(out) < 3:
            s = "".join(rng.choice("ACGT") for _ in range(250))
            if not any(x in s for x in ("GGTCTC", "GAGACC", "CGTCTC", "GAGACG")):
                out.append(s)
        return out
    def ring(ov):
        c = cores()
        return [gg(f"g{i}", ov[i], c[i], ov[(i + 1) % 3]) for i in range(3)]
    d = A.analyse(ring(["GGAG", "GGAG", "TCGC"]), 1, True, "BsaI")
    flagged = [c for c in d.checks if c[0].endswith("unique") and c[1] == "bad"]
    assert flagged and "same overhang" in flagged[0][2] and not d.ok
    assert A.analyse(ring(["GGAG", "TCGC", "CAGT"]), 1, True, "BsaI").ok   # the same ring with three different ones is fine


# ---- Local Align ---------------------------------------------------------------------------------------------------------

ALIGN_SEQ = random_dna(3000, 21)


def align_hits(win, query):
    dlg = LocalAlignDialog(win, win.editor)
    dlg.query.setPlainText(query)
    dlg.run()
    app.processEvents()
    return dlg


def test_each_local_align_hit_is_its_own_row_in_the_alignment_view(errors):
    """Overlapping / repeated hits sit on rows of their own, labelled in the names column; the dialog's table and the view
    agree on order and numbering and picking either picks the other."""
    from genome_view import HITS_HEAD, LEFT
    unit = ALIGN_SEQ[200:260]
    seq = ALIGN_SEQ[:600] + unit + ALIGN_SEQ[600:900] + reverse_complement(unit) + ALIGN_SEQ[900:1200] + unit + ALIGN_SEQ[1200:]
    win = window(Document.new_sequence("rows", seq, False))
    ed = win.editor
    dlg = align_hits(win, unit)
    assert len(dlg.hits) == 4, "the original, two direct copies and one inverted one"
    assert ed.query_hits is dlg.hits, "the view draws the table's list: same order and numbers"
    canvas = win.genome_view.canvas
    win.grab()
    name, top, height = next(r for r in canvas.layout_rows()[0] if r[0] == "Query hits")
    lane = canvas.lane_h("Query hits")
    assert height == HITS_HEAD + 4 * lane, "one row per alignment"
    bars = {a[1]: r for r, a in canvas.hits if a[0] == "qhit"}
    assert sorted(bars) == [0, 1, 2, 3]
    assert len({round(r.top()) for r in bars.values()}) == 4, "no two hits on the same row"
    assert [round(bars[i].top()) for i in range(4)] == [round(top + HITS_HEAD + i * lane) for i in range(4)]
    assert canvas.hit_label(1).startswith("#1 ") and any(canvas.hit_label(i + 1).split()[1] == "↺" for i in range(4)),         "the inverted hit's label carries the inversion mark"
    assert sum(h.strand < 0 for h in dlg.hits) == 1
    # Click the label of row 3 (names column): the hit is selected and its row in the table.
    from PySide6.QtCore import QPointF, Qt
    from PySide6.QtTest import QTest
    y = top + HITS_HEAD + 2 * lane + lane / 2
    assert canvas.hit_row_at(y) == 2 and canvas.hit_row_at(top + 2) == -1
    QTest.mouseClick(canvas, Qt.LeftButton, pos=QPointF(LEFT / 2, y).toPoint())
    h3 = dlg.hits[2]
    assert (ed.start, ed.end) == (h3.start, h3.end) and ed.current_query_hit == 2 and dlg.table.currentRow() == 2
    # Click the table's row 1: the view's picked row follows.
    dlg.table.selectRow(0)
    assert ed.current_query_hit == 0 and (ed.start, ed.end) == (dlg.hits[0].start, dlg.hits[0].end)
    # Painting with a picked row, a condensed track, and a track made shorter than its rows does not fail.
    win.grab()
    canvas.set_style("Query hits", mode="condensed")
    win.grab()
    assert canvas.layout_rows()[0][next(i for i, r in enumerate(canvas.layout_rows()[0]) if r[0] == "Query hits")][2]         == HITS_HEAD + 4 * canvas.lane_h("Query hits")
    canvas.set_style("Query hits", mode="expanded", height=HITS_HEAD + lane)
    win.grab()
    canvas.set_style("Query hits", height=None)
    dlg.clear()
    assert not ed.query_hits and ed.current_query_hit == -1
    close(win)


def test_local_align_hits(errors):
    """Local Align: exact copies on either strand, one SNP, nothing, bad queries, the hits as a track in the Alignment
    view; a copy across the origin of a circle."""
    win =window(Document.new_sequence("t", ALIGN_SEQ, False))
    ed = win.editor
    piece_ = ALIGN_SEQ[1000:1100]
    # A pasted exact copy: one hit, at the right place, on the strand it is on.
    dlg = align_hits(win, piece_)
    assert [(h.start, h.end, h.strand, h.identity) for h in dlg.hits] == [(1000, 1100, 1, 1.0)]
    assert dlg.table.rowCount() == 1 and dlg.table.item(0, 2).text() == "1,001–1,100" and dlg.table.item(0, 3).text() == FORWARD_TEXT
    assert dlg.table.item(0, 6).text() == "0 mismatches, 0 ins, 0 del" and dlg.status.text() == "1 hit"
    dlg = align_hits(win, reverse_complement(piece_))
    assert [(h.start, h.end, h.strand) for h in dlg.hits] == [(1000, 1100, -1)] and dlg.table.item(0, 3).text() == INVERTED_TEXT
    # One copy of each strand in the same sequence: two hits, one per strand.
    both = Document.new_sequence("b", ALIGN_SEQ[:1500] + reverse_complement(piece_) + ALIGN_SEQ[1600:], False)
    win2 = window(both, "b.gb")
    dlg2 = align_hits(win2, piece_)
    assert sorted((h.start, h.end, h.strand) for h in dlg2.hits) == [(1000, 1100, 1), (1500, 1600, -1)]
    close(win2)
    # One SNP: one mismatch, at the right base, nothing else.
    snp = list(piece_)
    snp[40] = "A" if snp[40] != "A" else "C"
    dlg = align_hits(win, "".join(snp))
    top = dlg.hits[0]
    assert len(dlg.hits) == 1 and (top.start, top.end) == (1000, 1100) and top.insertions == [] and top.deletions == []
    assert [m[0] for m in top.mismatches] == [1040] and top.mismatches[0][1] == snp[40] and abs(top.identity - 0.99) < 1e-6
    assert dlg.table.item(0, 6).text() == "1 mismatch, 0 ins, 0 del"
    # The hits become a track in the Alignment view.
    canvas = win.genome_view.canvas
    win.grab()
    assert len(ed.query_hits) == 1 and "Query hits" in canvas.kinds
    assert any(a[0] == "qhit" and (ed.query_hits[a[1]].start, ed.query_hits[a[1]].end) == (1000, 1100) for _, a in canvas.hits)
    # A sequence that is nowhere in it: no hit, no track.
    for seed in (31, 32, 33):
        dlg = align_hits(win, random_dna(80, seed))
        assert dlg.hits == [] and dlg.status.text() == "No hit" and dlg.table.rowCount() == 0
    win.grab()
    assert ed.query_hits == [] and "Query hits" not in canvas.kinds
    # Text that is not a query says why.
    dlg = align_hits(win, "AC")
    assert dlg.status.text() == "Query too short" and dlg.hits == []
    dlg = align_hits(win, "ACGT" + "Z" * 10)
    assert dlg.hits == [] and dlg.status.text() and dlg.status.text() != "No hit"
    close(win)

    # --- a copy across the origin of a circle ---
    ed =editor_of(ALIGN_SEQ, circular=True)
    query = ALIGN_SEQ[-30:] + ALIGN_SEQ[:30]
    parent = QWidget()
    dlg = LocalAlignDialog(parent, ed)
    dlg.query.setPlainText(query)
    dlg.win = type("W", (), {"set_genome_visible": lambda self, on: None,
                             "genome_view": type("G", (), {"set_view": lambda self, a, b: None})()})()
    dlg.run()
    assert [(h.start, h.end, h.strand) for h in dlg.hits] == [(2970, 3030, 1)]
    assert dlg.table.item(0, 2).text() == "2,971–30"            # the end is written on past the origin


# ---- Primer checks -------------------------------------------------------------------------------------------------------

DEFAULT_LIMITS = {"tm_low": 52.0, "tm_high": 65.0, "tm_margin": 7.0, "gc_low": 40, "gc_high": 60, "gc_margin": 10,
                  "clamp_min": 1, "clamp_max": 3, "hairpin_bad": -6.0, "dimer3_bad": -8.0, "dimer_bad": -12.0}


@pytest.fixture
def checks(monkeypatch):
    """primer_checks with the documented default limits, the default salt conditions, no genome and a Tm we choose."""
    limits = dict(DEFAULT_LIMITS)
    tm = {"value": 60.0}
    monkeypatch.setattr(PP, "limit", lambda key: limits[key])
    monkeypatch.setattr(PP, "conditions", lambda: {key: default for key, _, default, _ in PP.CONDITIONS})
    monkeypatch.setattr(PP, "melting_temp", lambda seq: tm["value"])
    monkeypatch.setattr(PP, "genome", lambda: None)
    ed = editor_of(random_dna(400, 9))

    def run(seq, name=None):
        result = PP.primer_checks(ed, seq)
        return {n: lv for n, lv, _ in result} if name is None else next(lv for n, lv, _ in result if n == name)
    run.limits, run.tm = limits, tm
    return run


def test_primer_settings_and_check_defaults():
    """The primer settings live in the data folder, not the registry; the primer check defaults are the documented
    limits and the band function flips at them."""
    # --- the primer settings live in the data folder, not the registry ---
    s = primer_db.settings()
    assert s.format() == QSettings.IniFormat and Path(s.fileName()) == primer_db.data_dir() / "primers.ini"
    assert PP.settings is primer_db.settings
    s.setValue("pcr/mg_mm", 3.0)
    s.sync()
    try:
        assert PP.conditions()["mg_mm"] == 3.0
    finally:
        s.remove("pcr/mg_mm")
        s.sync()

    # --- the defaults are the documented limits ---
    declared = {p.key: p.default for p in tool_settings.TOOLS[PP.CHECKS]}
    assert declared == DEFAULT_LIMITS
    assert [PP.band(v, 52, 65, 7) for v in (51.9, 52, 65, 65.1, 72, 72.1, 45, 44.9)] == \
        ["warn", "ok", "ok", "warn", "warn", "bad", "warn", "bad"]


def gc_primer(n_gc, n=20):
    body = list("G" * n_gc + "A" * (n - n_gc))
    random.Random(3).shuffle(body)
    return "".join(body)


def test_primer_check_levels(checks):
    """The Tm, GC, 3' clamp and hairpin / self-dimer levels change at their thresholds (and follow the limit settings);
    the verdicts roll up and short primers have none."""
    # --- Tm ---
    seq = "ACGTAGCTAGCTAGGATCCA"
    for tm, level in [(51.9, "warn"), (52.0, "ok"), (58.0, "ok"), (65.0, "ok"), (65.1, "warn"), (72.0, "warn"), (72.1, "bad"),
                      (45.0, "warn"), (44.9, "bad"), (20.0, "bad")]:
        checks.tm["value"] = tm
        assert checks(seq, "Tm") == level, tm
    checks.limits["tm_low"], checks.limits["tm_margin"] = 55.0, 2.0      # the limits are the settings: moving them moves the verdict
    checks.tm["value"] = 53.0
    assert checks(seq, "Tm") == "warn"
    checks.tm["value"] = 52.9
    assert checks(seq, "Tm") == "bad"

    # --- GC ---
    checks.limits.update(DEFAULT_LIMITS)
    got ={k: checks(gc_primer(k), "GC") for k in (4, 5, 6, 7, 8, 12, 13, 14, 15, 20)}   # of 20 bases: 5 % per base
    assert got == {4: "bad", 5: "bad", 6: "warn", 7: "warn", 8: "ok", 12: "ok", 13: "warn", 14: "warn", 15: "bad", 20: "bad"}
    checks.limits["gc_low"], checks.limits["gc_margin"] = 45, 5
    assert [checks(gc_primer(k), "GC") for k in (8, 9, 10)] == ["warn", "ok", "ok"]
    assert checks("gcgcgcgcATATATAT", "GC") == "ok"   # lower case counts as the same bases (50 %)

    # --- 3' clamp ---
    checks.limits.update(DEFAULT_LIMITS)
    ends ={"ATATA": 0, "ATATG": 1, "AGCTA": 2, "AGCTG": 3, "CGCTG": 4, "GCGCG": 5}
    got = {n: checks("ACTGATGCATGACAT" + end, "3′ clamp") for end, n in ends.items()}
    assert got == {0: "warn", 1: "ok", 2: "ok", 3: "ok", 4: "warn", 5: "warn"}
    checks.limits["clamp_max"] = 4
    assert checks("ACTGATGCATGACAT" + "CGCTG", "3′ clamp") == "ok"

    # --- hairpin and self-dimer ---
    checks.limits.update(DEFAULT_LIMITS)
    clean,loop, strong = "ATGCCGATTAGCCTAGCAAGTCCGA", "GGATCCAAAAAGGATCC", "GCGCGCGCAAAAGCGCGCGC"
    assert checks(clean, "Hairpin") == "ok" and checks(clean, "Self-dimer") == "ok"
    assert checks(loop, "Hairpin") == "warn" and checks(strong, "Hairpin") == "bad"     # usual limit (3' end < -2 or any < -3), then < -6
    assert checks("CAGTCGGAGTCCGACTG", "Self-dimer") == "warn" and checks("ACGCGTACGCGTACGCGT", "Self-dimer") == "bad"
    # The exact edge: the verdict flips where the limit crosses the primer's own ΔG.
    s = PP.structure(strong)
    checks.limits["hairpin_bad"] = s["hairpin"] - 0.01
    assert checks(strong, "Hairpin") == "warn"
    checks.limits["hairpin_bad"] = s["hairpin"] + 0.01
    assert checks(strong, "Hairpin") == "bad"
    d = PP.structure("ACGCGTACGCGTACGCGT")
    checks.limits.update(dimer3_bad=d["self3"] - 0.01, dimer_bad=d["self"] - 0.01)
    assert checks("ACGCGTACGCGTACGCGT", "Self-dimer") == "warn"
    checks.limits["dimer3_bad"] = d["self3"] + 0.01
    assert checks("ACGCGTACGCGTACGCGT", "Self-dimer") == "bad"
    checks.limits.update(dimer3_bad=-50.0, dimer_bad=d["self"] + 0.01)     # only the "anywhere" limit
    assert checks("ACGCGTACGCGTACGCGT", "Self-dimer") == "bad"

    # --- the verdicts roll up and short primers have none ---
    assert PP.primer_checks(editor_of("ACGT" * 50), "ACGT") == [("Primer", "info", "Type at least 8 bases.")]
    assert PP.worst([("a", "ok", ""), ("b", "warn", ""), ("c", "info", "")]) == "warn"
    assert PP.worst([("a", "ok", ""), ("b", "bad", ""), ("c", "warn", "")]) == "bad"
    assert PP.worst([("a", "info", "")]) == "info" and PP.worst([]) == "info"
    html = PP.checks_html([("Tm", "ok", "60.0 °C"), ("GC", "bad", "10 %")])
    assert PP.MARKS["ok"][0] in html and PP.MARKS["bad"][0] in html and "60.0 °C" in html


# ---- Feature table ---------------------------------------------------------------------------------------------------------


def table_features():
    return [("promoter", "300..350", [("label", "zeta")]), ("CDS", "100..250", [("label", "Alpha")]),
            ("terminator", "complement(500..520)", [("label", "mid")]), ("CDS", "10..40", [("label", "beta")])]


def column(table, c):
    return [table.model.data(table.model.index(r, c)) for r in range(table.model.rowCount())]


def test_feature_table(monkeypatch):
    """The feature table sorts by each column and filters by text; editing a feature updates the document and the table
    in one undo step (and a location the sequence cannot hold changes nothing)."""
    # --- sorts by each column and filters by text ---
    ed = editor_of(random_dna(800, 2), features=table_features())
    table = FeatureTable(ed)
    assert table.view.editTriggers() == QAbstractItemView.NoEditTriggers   # cells are not edited in place (see the next test)
    names = lambda: column(table, 0)
    table.view.sortByColumn(0, Qt.AscendingOrder)
    assert names() == ["Alpha", "beta", "mid", "zeta"]                         # case-insensitive
    table.view.sortByColumn(0, Qt.DescendingOrder)
    assert names() == ["zeta", "mid", "beta", "Alpha"]
    table.view.sortByColumn(1, Qt.AscendingOrder)
    assert column(table, 1) == ["CDS", "CDS", "promoter", "terminator"]
    table.view.sortByColumn(2, Qt.AscendingOrder)                             # by position
    assert names() == ["beta", "Alpha", "zeta", "mid"]
    assert column(table, 2) == ["10..40", "100..250", "300..350", "complement(500..520)"]
    table.view.sortByColumn(3, Qt.AscendingOrder)                             # by length
    assert column(table, 3) == ["21", "31", "51", "151"] and names() == ["mid", "beta", "zeta", "Alpha"]
    table.view.sortByColumn(3, Qt.DescendingOrder)
    assert column(table, 3) == ["151", "51", "31", "21"]
    table.view.sortByColumn(4, Qt.AscendingOrder)
    assert column(table, 4)[0] == "←" and column(table, 4)[1:] == ["→"] * 3   # the reverse-strand feature first
    table.view.sortByColumn(2, Qt.AscendingOrder)

    for text, expected in [("cds", ["beta", "Alpha"]), ("ALPHA", ["Alpha"]), ("complement", ["mid"]), ("term", ["mid"]),
                           ("a", ["beta", "Alpha", "zeta", "mid"]), ("nothing like this", []), ("", ["beta", "Alpha", "zeta", "mid"])]:
        table.filter.setText(text)       # the box under the table: case does not matter, name / type / location are searched
        assert names() == expected, text
    assert len(ed.features) == 4          # filtering shows less, it deletes nothing

    # --- editing a feature updates the document and the table in one undo step ---
    ed = editor_of(random_dna(800, 2), features=table_features())
    table = FeatureTable(ed)
    asked = []
    ed.edit_feature_requested.connect(asked.append)
    row = table.model.row_of(by_label(ed, "Alpha").index)
    table.view.doubleClicked.emit(table.model.index(row, 0))       # a double click asks for the editor of that feature
    assert asked == [by_label(ed, "Alpha").index]
    steps = lambda: len(ed.doc.history()[0])
    before = steps()
    dialog = FeatureDialog(None, ed, asked[0])
    dialog.name.setText("Omega")
    dialog.kind.setCurrentText("gene")
    dialog.location.setText("120..260")
    dialog.accept()
    assert steps() == before + 1                                   # one step for name, type and location together
    f = by_label(ed, "Omega")
    assert (f.kind, f.location) == ("gene", "120..260") and not any(x.label == "Alpha" for x in ed.features)
    texts = [table.model.data(table.model.index(r, c)) for r in range(table.model.rowCount()) for c in (0, 1, 2)]
    assert "Omega" in texts and "gene" in texts and "120..260" in texts and "Alpha" not in texts   # the table follows
    assert ed.doc.history()[0][-1][5:7] == (["Omega"], ["Alpha"])
    assert ed.undo()
    assert by_label(ed, "Alpha").location == "100..250" and steps() == before
    assert "Alpha" in column(table, 0) and "Omega" not in column(table, 0)
    ed.redo()
    assert by_label(ed, "Omega").location == "120..260"
    # A location the sequence cannot hold changes nothing and adds no step.
    monkeypatch.setattr(QMessageBox, "warning", lambda *a, **k: None)
    bad = FeatureDialog(None, ed, by_label(ed, "Omega").index)
    bad.location.setText("5..999999")
    bad.accept()
    assert by_label(ed, "Omega").location == "120..260" and steps() == before + 1


# ---- Settings --------------------------------------------------------------------------------------------------------------

TOOL = "Zz test tool"


@pytest.fixture
def settings_sandbox():
    """A tool of our own in the settings, and the shared settings.ini put back as it was (other tests read it)."""
    store = tool_settings.store()
    snapshot = {k: store.value(k) for k in store.allKeys()}
    tool_settings.register(TOOL, [
        tool_settings.Param("ratio", "A decimal", 1.5, 0, 10, step=0.5), tool_settings.Param("count", "A count", 3, 0, 10),
        tool_settings.Param("flag", "A switch", True), tool_settings.Param("mode", "A choice", "fast", choices=["fast", "slow"])])
    tool_settings._cache.clear()
    heard = []
    tool_settings.notifier.changed.connect(heard.append)
    yield heard
    tool_settings.notifier.changed.disconnect(heard.append)
    tool_settings.TOOLS.pop(TOOL, None)
    store = tool_settings.store()
    for key in store.allKeys():
        if key not in snapshot:
            store.remove(key)
        elif store.value(key) != snapshot[key]:
            store.setValue(key, snapshot[key])
    store.sync()
    tool_settings._cache.clear()


def tab_button(dialog, tool, text):
    page = dialog.page_for(tool)
    return next(b for b in page.findChildren(QPushButton) if b.text() == text)


def test_settings(settings_sandbox):
    """Settings persist in settings.ini and fall back to defaults; the dialog saves, notifies and restores defaults."""
    # --- persist in settings.ini and fall back to defaults ---
    assert [tool_settings.value(TOOL, k) for k in ("ratio", "count", "flag", "mode")] == [1.5, 3, True, "fast"]
    for key, v in (("ratio", 2.5), ("count", 7), ("flag", False), ("mode", "slow")):
        tool_settings.set_value(TOOL, key, v)
    ini = tool_settings.data_dir() / "settings.ini"
    assert ini.exists() and "ratio=2.5" in ini.read_text()
    tool_settings._cache.clear()                                       # as after a restart: read again from the file
    assert [tool_settings.value(TOOL, k) for k in ("ratio", "count", "flag", "mode")] == [2.5, 7, False, "slow"]   # typed back from text
    fresh = QSettings(str(ini), QSettings.IniFormat)
    assert fresh.value(f"{TOOL}/count") == "7" or int(fresh.value(f"{TOOL}/count")) == 7
    # A value the file cannot give back as the declared type, or a choice that no longer exists: the default.
    tool_settings.set_value(TOOL, "count", "abc")
    tool_settings.set_value(TOOL, "mode", "gone")
    assert tool_settings.value(TOOL, "count") == 3 and tool_settings.value(TOOL, "mode") == "fast"

    # --- the dialog saves, notifies and restores defaults ---
    heard = settings_sandbox
    heard.clear()
    store = tool_settings.store()
    for key in ("ratio", "count", "flag", "mode"):   # back to the declared defaults, as for a fresh test
        store.remove(f"{TOOL}/{key}")
    tool_settings._cache.clear()
    dialog =tool_settings.SettingsDialog(None, TOOL)
    widgets = {k: w for (name, k), (_, w) in dialog.widgets.items() if name == TOOL}
    assert dialog.tabs.currentPage() == TOOL
    assert widgets["ratio"].value() == 1.5 and widgets["count"].value() == 3 and widgets["flag"].isChecked()
    dialog.accept()                                                    # nothing changed: nothing saved, nobody told
    assert heard == []
    dialog = tool_settings.SettingsDialog(None, TOOL)
    widgets = {k: w for (name, k), (_, w) in dialog.widgets.items() if name == TOOL}
    widgets["ratio"].setValue(4.5)
    widgets["count"].setValue(9)
    widgets["flag"].setChecked(False)
    widgets["mode"].setCurrentText("slow")
    dialog.accept()
    assert heard == [TOOL]                                             # one notice per tool that changed
    assert [tool_settings.value(TOOL, k) for k in ("ratio", "count", "flag", "mode")] == [4.5, 9, False, "slow"]
    tool_settings._cache.clear()
    assert tool_settings.value(TOOL, "ratio") == 4.5                   # and it is in the file
    # "Restore Defaults" puts the declared values back in the form; they are saved only by OK.
    heard.clear()
    dialog = tool_settings.SettingsDialog(None, TOOL)
    widgets = {k: w for (name, k), (_, w) in dialog.widgets.items() if name == TOOL}
    tab_button(dialog, TOOL, "Restore Defaults").click()
    assert [widgets["ratio"].value(), widgets["count"].value(), widgets["flag"].isChecked(), widgets["mode"].currentText()] == \
        [1.5, 9 - 6, True, "fast"]
    assert tool_settings.value(TOOL, "ratio") == 4.5 and heard == []
    dialog.accept()
    assert heard == [TOOL] and [tool_settings.value(TOOL, k) for k in ("ratio", "count", "flag", "mode")] == [1.5, 3, True, "fast"]
    # A change that reaches a tool: the primer checks read their limits from here.
    dialog = tool_settings.SettingsDialog(None, PP.CHECKS)
    w = dialog.widgets[(PP.CHECKS, "tm_low")][1]
    old = tool_settings.value(PP.CHECKS, "tm_low")
    heard.clear()
    w.setValue(40.0)
    dialog.accept()
    assert heard == [PP.CHECKS] and PP.limit("tm_low") == 40.0 != old


# ---- Self-test panel -------------------------------------------------------------------------------------------------------


def test_self_test_panel_counts_passes_and_failures_and_shows_the_message(monkeypatch):
    def passes():
        pass

    def fails():
        raise AssertionError("one plus one is not three")

    def fails_silently():
        raise AssertionError()      # no message of its own

    def crashes():
        raise ValueError("boom")

    monkeypatch.setattr(selftest, "CHECKS", [("Area", "passes", passes), ("Area", "fails", fails),
                                             ("Other", "fails silently", fails_silently), ("Other", "crashes", crashes)])
    panel = selftest.SelfTestPanel()
    assert panel.status.text() == "Not run" and panel.table.rowCount() == 4
    panel.run_all()
    assert panel.status.text() == "1 passed, 3 failed"
    cell = lambda r, c: panel.table.item(r, c).text()
    assert [cell(r, 2) for r in range(4)] == [PP.MARKS["ok"][0]] + [PP.MARKS["bad"][0]] * 3
    assert cell(0, 4) == "" and cell(1, 4) == "one plus one is not three"
    assert cell(2, 4) == "assertion failed" and cell(3, 4) == "ValueError: boom"        # a crash is a failure with its type
    assert [cell(r, 0) for r in range(4)] == ["Area", "Area", "Other", "Other"] and cell(1, 1) == "fails"
    assert all(cell(r, 3).isdigit() for r in range(4))                                # a time for each
    # All passing: "N passed, 0 failed" (the real checks run under pytest in test_selftest.py).
    monkeypatch.setattr(selftest, "CHECKS", [("Area", "a", passes), ("Area", "b", passes)])
    panel = selftest.SelfTestPanel()
    panel.run_all()
    assert panel.status.text() == "2 passed, 0 failed"
