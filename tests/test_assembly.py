"""Assembly: the engine (junctions found at the real ends of the fragments, primers that add what is
missing, the exact product), the overhang checks, and the Assembly tab."""
import os
import random
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ["OPENDNA_DATA_DIR"] = tempfile.mkdtemp(prefix="opendna_test_")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "app"))

import pytest
from PySide6.QtCore import QPoint, QPointF, Qt
from PySide6.QtGui import QMouseEvent, QPixmap
from PySide6.QtWidgets import QApplication, QWidget

from opendna_core import Document

import assembly as A
import assembly_panel as P
import overhangs as OH
from citations import citations_html

app = QApplication.instance() or QApplication([])


@pytest.fixture(autouse=True)
def private_primer_db(tmp_path, monkeypatch):
    """These tests save primers: into a database of their own, not the one the other test files share (a bigger one makes
    the genome test of test_gui slow)."""
    import primer_db
    monkeypatch.setattr(primer_db, "_db", primer_db.PrimerDB(tmp_path / "primers.tsv"))
rng = random.Random(1)
rc = A.rc


def dna(n):
    """Random bases without a BsaI / BsmBI site on either strand (a test adds the ones it wants)."""
    while True:
        s = "".join(rng.choice("ACGT") for _ in range(n))
        if not any(site in s for site in ("GGTCTC", "GAGACC", "CGTCTC", "GAGACG")):
            return s


def frag(name, seq, features=(), **kw):
    doc = Document.new_sequence(name, seq)
    if features:
        doc.add_features([(k, loc, [("label", lab)]) for k, loc, lab in features])
    return A.Frag(name, doc, **kw)


def frags():
    return [frag("A", dna(300)), frag("B", dna(400)), frag("C", dna(250))]


def vector_seq(tail=1396):
    """A destination vector: the two sites face away from the backbone, and its own end bases (what the first and last
    insert must overhang with) are two good overhangs, GGAG and TCGC."""
    return dna(296) + "GGAG" + "T" + "GAGACC" + dna(120) + "GGTCTC" + "A" + "TCGC" + dna(tail)


def vector_frag():
    return A.Frag("vec", Document.new_sequence("vec", vector_seq(), True), vector=True)


def find_check(d, name):
    return next(c for c in d.checks if c[0] == name)


def rotation_of(a, b):
    return len(a) == len(b) and a in b + b


def joined_bases(d, j):
    """What the product holds at the junction's position (reading on round a circle)."""
    p, n = d.product(), len(d.product())
    return "".join(p[q % n] for q in range(j.pos, j.pos + len(j.seq)))


# ---- Fragment ----------------------------------------------------------------------------------------


# ---- Gibson ---------------------------------------------------------------------------------------------


def frags_like(done):
    """Plain copies of fragments that already have primers (their cores)."""
    return [A.Frag(f.name, f.core, vector=f.vector) for f in done]


def test_fragment_and_gibson_engine():
    """A fragment is a PCR product and flips with its primers; Gibson: ends that do not fit are reported then fixed by new primers, existing overlap and uniqueness; every junction sits where it says in the product."""
    # --- a fragment is a PCR product and flips with its primers ---
    template = dna(600)
    doc = Document.new_sequence("t", template)
    fwd = ("F", "GGGGCC" + template[100:122], 22)
    rev = ("R", "TTTAAA" + rc(template[300:322]), 22)
    f = A.Frag("x", doc, 100, 322, fwd, rev)
    assert f.seq == "GGGGCC" + template[100:322] + rc("TTTAAA") and f.tails == (6, 6)
    assert f.core_seq == template[100:322]
    before, tails = f.seq, f.tails
    f.flip()
    assert f.seq == rc(before) and f.tails == tails and f.flipped
    assert f.fwd[0] == "R" and f.rev[0] == "F", "the primers swap roles when the fragment is reverse-complemented"
    f.flip()
    assert f.seq == before

    # --- Gibson: ends that do not fit, then new primers ---
    f = frags()
    d = A.analyse(f, 0, True)
    assert [j.kind for j in d.joints] == ["none"] * 3 and not d.ok
    assert find_check(d, "Circle closes")[1] == "bad" and "no homology" in d.problems[0]
    assert A.auto_primers(f, 0, True) == 3
    d = A.analyse(f, 0, True)
    assert d.ok and [(j.left, j.right) for j in d.joints] == [("A", "B"), ("B", "C"), ("C", "A")]
    assert all(20 <= len(j.seq) <= 40 and A.joint_tm(j.seq) >= 52 for j in d.joints)
    assert len(d.rows) == 6 and len(d.product()) == 950
    # The very product the older designer built from the plain fragments (same bases, same loop).
    assert rotation_of(d.product(), A.design_gibson(frags_like(f), True).product())
    assert find_check(d, "Circle closes")[1] == "ok"

    # --- Gibson: existing overlap and uniqueness ---
    shared = dna(30)
    d = A.analyse([frag("a", dna(200) + shared), frag("b", shared + dna(300))], 0, False)
    assert d.ok and d.joints[0].kind == "existing overlap" and len(d.joints[0].seq) >= 30 and len(d.product()) == 530
    # The same bases twice in the product: the ends could join in the wrong place, a red error.
    twice = A.analyse([frag("a", dna(200) + shared), frag("b", shared + dna(100) + shared + dna(100))], 0, False)
    assert find_check(twice, "Homology unique (a → b)")[1] == "bad" and not twice.ok
    # Short homology still fits but is flagged.
    short = A.analyse([frag("a", dna(200) + shared[:16]), frag("b", shared[:16] + dna(300))], 0, False)
    assert short.joints[0].kind == "existing overlap" and find_check(short, "Homology (a → b)")[1] == "warn"

    # --- every junction sits where it says in the product ---
    f = frags()
    A.auto_primers(f, 0, True)
    designs = {"gibson circular": A.analyse(f, 0, True), "gibson linear": A.analyse(f, 0, False),
               "sidewinder": A.analyse(frags(), 2, False)}
    g = [vector_frag(), frag("B", dna(500)), frag("C", dna(300))]
    A.auto_primers(g, 1, True)
    designs["golden gate"] = A.analyse(g, 1, True)
    for name, d in designs.items():
        for j in d.joints:
            if j.kind != "none":
                assert joined_bases(d, j) == j.seq, f"{name}: {j.left} → {j.right}"
    for j in designs["golden gate"].joints:
        assert (j.a[0], j.a[1], j.b[0], j.b[1]) == (j.pos, j.pos + 4, j.pos, j.pos + 4)
    for j in designs["sidewinder"].joints:
        assert j.a[0] == j.b[0] == j.pos + A.TOEHOLD


# ---- Golden Gate --------------------------------------------------------------------------------------------


def gg_fragment(name, left, core, right, spacer="a", pad="tt"):
    """A PCR product the way Golden Gate wants it: site + spacer + overhang, the core, overhang + spacer + site."""
    return frag(name, pad + "GGTCTC" + spacer + left + core + right + spacer.replace("a", "t") + "GAGACC" + pad)


def ring(overhangs, cores=None):
    n = len(overhangs)
    cores = cores or [dna(250) for _ in range(n)]
    return [gg_fragment(f"g{i}", overhangs[i], cores[i], overhangs[(i + 1) % n]) for i in range(n)]


def test_golden_gate_engine_and_overhang_checks():
    """Golden Gate: vector and inserts become one circle, the overhang each end really carries is read, poor overhang choices are warned about, and the paper's overhang sets."""
    # --- vector and inserts become one circle ---
    f = [vector_frag(), frag("ins", dna(700)), frag("term", dna(220))]
    d = A.analyse(f, 1, True, "BsaI")
    assert [j.kind for j in d.joints] == ["none"] * 3 and not d.ok
    assert A.auto_primers(f, 1, True, "BsaI") == 2, "the vector is cut, not amplified"
    d = A.analyse(f, 1, True, "BsaI")
    assert d.ok and len(d.rows) == 4 and all(len(j.seq) == 4 for j in d.joints)
    assert len({j.seq for j in d.joints}) == 3
    assert rotation_of(d.product(), A.design_golden_gate(frags_like(f), True, "BsaI").product())
    assert all(dict((s, lv) for s, lv, _ in d.jchecks[j])["unique"] == "ok" for j in range(3))
    assert find_check(d, "BsaI sites left in product")[1] == "ok"
    # Reverse-complement the insert: its ends no longer fit until it gets new primers.
    f[1].flip()
    assert not A.analyse(f, 1, True, "BsaI").ok
    A.auto_primers(f, 1, True, "BsaI")
    assert A.analyse(f, 1, True, "BsaI").ok
    # A destination vector needs exactly two sites, facing away from the backbone.
    with pytest.raises(ValueError, match="needs exactly 2"):
        A.analyse([A.Frag("v", Document.new_sequence("v", dna(500), True), vector=True), frag("x", dna(100))], 1, True)
    same_way = "GGTCTCA" + dna(300) + "GGTCTCA" + dna(500)
    with pytest.raises(ValueError, match="point the same way"):
        A.analyse([A.Frag("v", Document.new_sequence("v", same_way, True), vector=True), frag("x", dna(100))], 1, True)

    # --- the overhang each end really carries ---
    ov = ["GGAG", "TCGC", "CAGT"]            # three of the paper's validated overhangs
    d = A.analyse(ring(ov), 1, True, "BsaI")
    assert d.ok and [j.seq for j in d.joints] == ov[1:] + ov[:1]
    assert len(d.product()) == 3 * (4 + 250), "each piece is its left overhang + core; the right overhang is the next piece's"
    assert all(dict((s, lv) for s, lv, _ in d.jchecks[j])["validated"] == "ok" for j in range(3))
    assert not [c for c in d.checks if c[1] in ("bad", "warn")]
    # Two ends that disagree: not joined, and it says which bases.
    bad = ring(ov)
    bad[1] = gg_fragment("g1", "TCGA", dna(250), "CAGT")
    d = A.analyse(bad, 1, True, "BsaI")
    assert d.joints[0].kind == "none" and "TCGC" in d.problems[0] and "TCGA" in d.problems[0] and not d.ok
    # An end with no site.
    nosite = ring(ov)
    nosite[2] = frag("g2", "CAGT" + dna(250) + "GGAG")
    d = A.analyse(nosite, 1, True, "BsaI")
    assert "no BsaI site" in d.problems[1] and not d.ok
    # A site inside the fragment is a red error naming it.
    inner = ring(ov, [dna(120) + "GGTCTC" + dna(120), dna(250), dna(250)])
    d = A.analyse(inner, 1, True, "BsaI")
    assert find_check(d, "BsaI sites in g0")[1] == "bad" and not d.ok

    # --- poor overhang choices ---
    def chips(ov):
        d = A.analyse(ring(ov), 1, True, "BsaI")
        return {(short, level) for cs in d.jchecks.values() for short, level, _ in cs}
    assert ("unique", "bad") in chips(["GGAG", "GGAG", "TCGC"])
    assert ("not palindromic", "bad") in chips(["GATC", "TCGC", "CAGT"])
    assert ("not palindromic", "bad") not in chips(["GGTT", "TCGC", "CAGT"])
    close = chips(["GGTT", "GGTA", "TCGC"])
    assert ("orthogonal", "warn") in close and ("validated", "info") in close
    assert ("GC mix", "warn") in chips(["AAAA", "TCGC", "CAGT"])
    assert ("orthogonal", "bad") in chips(["CAGT", "ACTG", "TCGC"]), "a reverse complement is certain to cross-ligate"
    assert not [c for c in chips(["GGAG", "TCGC", "CAGT"]) if c[1] in ("bad", "warn")]

    # --- the paper's overhang sets ---
    assert len(OH.VALIDATED) == 25 == len(set(OH.VALIDATED)) and all(len(o) == 4 for o in OH.VALIDATED)
    assert not [c for c in OH.assess(OH.VALIDATED) if c[1] in ("bad", "warn")]
    assert [c for c in OH.assess(OH.POOR_EXAMPLE) if c[1] == "warn"], "the paper's low-fidelity set is flagged"
    assert OH.pair_risk("GGAG", "CTCC") == "bad" and OH.pair_risk("GGAG", "GGAA") == "warn"
    assert {s: lv for s, lv, _ in OH.junction_checks("GGAG", ["TCGC", "CAGT"])} == {
        "unique": "ok", "not palindromic": "ok", "orthogonal": "ok", "GC mix": "ok", "validated": "ok"}
    assert {s: lv for s, lv, _ in OH.junction_checks("GGAA", ["GGAG"])}["orthogonal"] == "warn"
    # The designer prefers the validated overhangs, and never repeats one.
    assert A.pick_overhangs([[(0, "ACCA"), (1, "GGAG")]], [None], 4) == [(1, "GGAG")]
    assert A.pick_overhangs([[(0, "GGAG")], [(0, "GGAG"), (1, "TCGC")]], [None, None], 4) == [(0, "GGAG"), (1, "TCGC")]


# ---- Sidewinder, product, views ---------------------------------------------------------------------------------


def test_sidewinder_product_views_and_citations():
    """Sidewinder and the simulated product, the product keeps the features of every method, the junction view shows the sketch, reason and tails, and the citations name the sources."""
    # --- Sidewinder and the simulated product ---
    f = [frag("a", dna(100), [("promoter", "1..50", "p")]), frag("b", dna(120)), frag("c", dna(90))]
    d = A.analyse(f, 2, False)
    assert len(d.joints) == 2 and d.problems == [None, None] and d.ok
    doc = A.simulate(f, d)
    assert doc.length == 310 and not [x for x in doc.features() if x.kind == "misc_feature"], "the product carries no marks of its own"
    assert {x.label for x in A.simulate(f, d, preview=True).features() if x.kind == "misc_feature"} >= {"a", "b", "c"}, "the preview names the fragments"
    long = A.analyse(frags(), 2, False)
    assert not long.ok and find_check(long, f"Oligos at most {A.MAX_OLIGO} nt")[1] == "bad", "Sidewinder wants short fragments"
    with pytest.raises(ValueError, match="at least two"):
        A.analyse(f[:1], 2, False)

    # --- the product keeps the features of every method ---
    for method in (0, 1):
        f = ([vector_frag()] if method else [frag("v", dna(500), [("CDS", "50..200", "gene")])]) + [
            frag("i", dna(600), [("promoter", "10..80", "pro")])]
        A.auto_primers(f, method, True)
        d = A.analyse(f, method, True)
        assert d.ok, [(c, d.problems) for c in d.checks if c[1] == "bad"]
        doc = A.simulate(f, d)
        assert [x.label for x in doc.features() if x.label == "pro"] == ["pro"]
        assert doc.circular and doc.length == len(d.product())
    piece = A.piece_doc(f, d, 0)
    assert piece.length == len(d.seqs[0])

    # --- the junction view ---
    f = frags()
    plain = A.junction_html(A.analyse(f, 0, True))
    assert "no homology" in plain and plain.count("✗") >= 3
    A.auto_primers(f, 0, True)
    f[1].flip()
    A.auto_primers(f, 0, True)
    d = A.analyse(f, 0, True)
    html = A.junction_html(d)
    assert "reverse-complemented" in html and "as given" in html and d.note[:20] in html
    assert html.count("<table") == len(d.joints) and "[A]" in html and "bottom" in html
    for j in d.joints:   # the sketch: the joining bases and, under them, their complement
        assert j.seq in html and j.seq.translate(A.COMP) in html
    assert html.count("background-color:" + A.TAIL_BG) >= 2 * 2 * 3 * 5, "tails show lowercase on cream"
    g = [vector_frag(), frag("i", dna(500)), frag("t", dna(200))]
    A.auto_primers(g, 1, True)
    gh = A.junction_html(A.analyse(g, 1, True))
    gd = A.analyse(g, 1, True)
    assert "overhang" in gh and all(c and len(c) == 2 and j.seq in c[0] and j.seq in c[1] for c, j in zip(gd.cuts, gd.joints)),         "the two ends with the cut marked are shown at every junction, vector included"

    # --- the citations name the sources ---
    html = citations_html()
    for needle in ("10.1021/acssynbio.8b00333", "10.1038/nmeth.1318", "10.1371/journal.pone.0003647", "Sidewinder",
                   "SantaLucia", "PySide6", "PyO3", "10.1126/science.1225829", "10.1038/nbt.2647", "10.1038/nature14299",
                   "10.1038/s41597-024-03571-y", "10.1093/nar/15.5.2343", "10.1016/S0022-2836(05)80005-0",
                   "10.1073/pnas.71.4.1342"):
        assert needle in html


# ---- the product keeps every ORF on its own bases ---------------------------------------------------------------
# A scientist cuts ORFs out of genes and assembles them. Wherever a neighbour supplies the bases of a junction (a Gibson
# overlap, a Golden Gate overhang), the product must still annotate each ORF over all of its bases: a CDS that starts
# or ends a few bases late reads another frame, and the protein is wrong.

STOPS = {"TAA", "TAG", "TGA"}
SENSE = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT" if a + b + c not in STOPS]


def orf_frag(name, codons, seed, strand=1):
    """A fragment that is one ORF (ATG, sense codons, TAA) with its CDS annotated over all of it, the way an ORF is cut out
    of a gene. Fixed by `seed`; no BsaI / BsmBI site on either strand, so every method accepts it. Returns (fragment, the
    ORF's bases as designed)."""
    r = random.Random(seed)
    while True:
        orf = "ATG" + "".join(r.choice(SENSE) for _ in range(codons)) + "TAA"
        if not any(site in orf for site in ("GGTCTC", "GAGACC", "CGTCTC", "GAGACG")):
            break
    doc = Document.new_sequence(name, orf if strand > 0 else rc(orf))
    where = f"1..{len(orf)}" if strand > 0 else f"complement(1..{len(orf)})"
    doc.add_features([("CDS", where, [("label", f"{name} ORF"), ("codon_start", "1")])])
    return A.Frag(name, doc), orf


def three_orfs(seed=0, codons=(90, 110, 100), strands=(1, 1, 1)):
    """Fragments A, B, C, each one ORF, and the ORFs they must keep: {CDS label: bases}."""
    made = [orf_frag(name, n, seed * 10 + i, s) for i, (name, n, s) in enumerate(zip("ABC", codons, strands))]
    return [f for f, _ in made], {f.name + " ORF": orf for f, orf in made}


def read_cds(doc, feature):
    """The bases a feature reads: its parts in order, reverse-complemented on the minus strand. Taken from the sequence and
    the base ranges, not from the location text, so it does not depend on how the location is written."""
    seq = doc.sequence(0, None).upper()
    bases = "".join(seq[s:e] for s, e in feature.spans)
    return rc(bases) if feature.strand < 0 else bases


def assert_orfs_intact(product, orfs):
    """Every ORF is annotated in the product over exactly the bases it was designed with, so it still translates."""
    from Bio.Seq import Seq
    cds = {f.label: f for f in product.features() if f.kind == "CDS"}
    for label, bases in orfs.items():
        assert label in cds, f"{label} is gone from the product"
        reads = read_cds(product, cds[label])
        assert reads == bases, f"{label} reads {len(reads)} bases in the product instead of {len(bases)}: {reads[:30]}... {reads[-30:]}"
        assert dict(cds[label].qualifiers).get("codon_start") == "1", f"{label} lost its reading frame"
        protein = str(Seq(reads).translate())
        assert protein[0] == "M" and protein.endswith("*") and "*" not in protein[:-1], f"{label} does not translate as an ORF"


@pytest.mark.parametrize("overlap", range(20, 26))   # half of the overlap sits on the neighbour: 10 to 13 bases, so every frame offset
def test_gibson_with_new_primers_keeps_every_orf_on_its_bases(overlap):
    f, orfs = three_orfs()
    A.auto_primers(f, 0, False, min_ov=overlap, max_ov=overlap)
    d = A.analyse(f, 0, False)
    assert [j.kind for j in d.joints] == ["existing overlap"] * 2 and len(d.product()) == sum(map(len, d.seqs)) - 2 * overlap
    assert_orfs_intact(A.simulate(f, d), orfs)


@pytest.mark.parametrize("overlap", range(20, 26))
def test_a_circular_gibson_keeps_every_orf_including_the_one_that_closes_the_circle(overlap):
    f, orfs = three_orfs()
    A.auto_primers(f, 0, True, min_ov=overlap, max_ov=overlap)
    d = A.analyse(f, 0, True)
    assert [j.kind for j in d.joints] == ["existing overlap"] * 3 and len(d.product()) == sum(map(len, d.seqs)) - 3 * overlap
    product = A.simulate(f, d)
    assert product.circular
    assert_orfs_intact(product, orfs)


@pytest.mark.parametrize("seed", range(6))   # the overhang is picked where it scores best: it takes 0 to 4 bases off the end of an ORF
def test_golden_gate_keeps_every_orf_including_the_bases_in_an_overhang(seed):
    f, orfs = three_orfs(seed)
    g = [vector_frag()] + f
    A.auto_primers(g, 1, True)
    d = A.analyse(g, 1, True)
    assert [j.kind for j in d.joints] == ["overhang"] * 4
    product = A.simulate(g, d)
    assert product.circular
    assert_orfs_intact(product, orfs)


def test_gibson_keeps_an_orf_that_starts_inside_an_overlap_the_fragments_already_share():
    """Each fragment ends with the first 30 bases of the next ORF, unannotated (the next ORF's start is in the overlap)."""
    (a, orf_a), (b, orf_b), (c, orf_c) = [orf_frag(name, n, 40 + n) for name, n in zip("ABC", (60, 70, 80))]
    for frag_, orf, nxt in ((a, orf_a, orf_b), (b, orf_b, orf_c)):
        doc = Document.new_sequence(frag_.name, orf + nxt[:30])
        doc.add_features([("CDS", f"1..{len(orf)}", [("label", f"{frag_.name} ORF"), ("codon_start", "1")])])
        frag_.src, frag_.end = doc, doc.length
        frag_.rebuild()
    d = A.analyse([a, b, c], 0, False)
    assert [j.kind for j in d.joints] == ["existing overlap"] * 2 and all(len(j.seq) == 30 for j in d.joints)
    assert_orfs_intact(A.simulate([a, b, c], d), {"A ORF": orf_a, "B ORF": orf_b, "C ORF": orf_c})


@pytest.mark.parametrize("method, circular", [(0, False), (0, True), (1, True)])
def test_an_orf_on_the_minus_strand_keeps_its_bases_too(method, circular):
    """Its start sits at the fragment's right end and its stop at the left end: the other side of every junction."""
    f, orfs = three_orfs(strands=(1, -1, 1))
    g = ([vector_frag()] if method else []) + f
    A.auto_primers(g, method, circular)
    d = A.analyse(g, method, circular)
    assert all(j.kind != "none" for j in d.joints)
    assert_orfs_intact(A.simulate(g, d), orfs)


def test_sidewinder_keeps_every_orf():
    f, orfs = three_orfs(codons=(30, 40, 35))
    d = A.analyse(f, 2, False)
    assert len(d.joints) == 2
    assert_orfs_intact(A.simulate(f, d), orfs)


def test_a_feature_inside_an_existing_overlap_is_annotated_once():
    """Two stretches of one plasmid that overlap both carry the annotation of the overlap: the product shows it once."""
    src = Document.new_sequence("plasmid", dna(1000))
    src.add_features([("misc_feature", "476..495", [("label", "in overlap")]), ("promoter", "51..120", [("label", "left")]),
                      ("terminator", "801..900", [("label", "right")])])
    f = [A.Frag("left half", src, start=0, end=500), A.Frag("right half", src, start=470, end=1000)]
    d = A.analyse(f, 0, False)
    assert d.joints[0].kind == "existing overlap" and len(d.joints[0].seq) == 30 and len(d.product()) == 1000
    product = A.simulate(f, d)
    found = sorted((x.label, x.spans) for x in product.features() if x.label in ("in overlap", "left", "right"))
    assert found == [("in overlap", [(475, 495)]), ("left", [(50, 120)]), ("right", [(800, 900)])]


# ---- the panel -------------------------------------------------------------------------------------------------------


class Win(QWidget):
    """Just what the panel asks of the main window."""
    def __init__(self, docs=()):
        super().__init__()
        self.docs, self.opened = list(docs), None

    def open_documents(self):
        return [(d.name, d) for d in self.docs]

    def open_product(self, doc):
        self.opened = doc


def source_doc(name, n, features=(), circular=False):
    doc = Document.new_sequence(name, dna(n), circular)
    if features:
        doc.add_features([(k, loc, [("label", lab)]) for k, loc, lab in features])
    return doc


def filled_panel(method=0):
    docs = [source_doc("vector", 900, [("CDS", "100..400", "gene")]), source_doc("insert", 600), source_doc("term", 300)]
    panel = P.AssemblyPanel(Win(docs))
    panel.method.setCurrentIndex(method)
    for i, d in enumerate(docs):
        panel._set_source(i, A.snapshot(d), f"{d.name}.gb", None)
    return panel


def click(strip, x, kind=QMouseEvent.Type.MouseButtonPress, button=Qt.LeftButton, y=None):
    pos = QPointF(x, strip.boxes[0].center().y() if y is None else y)
    strip.event(QMouseEvent(kind, pos, QPointF(strip.mapToGlobal(QPoint(int(x), 0))), button, button, Qt.NoModifier))


def click_tile(tile):
    pos = QPointF(tile.width() / 2, tile.height() / 2)
    tile.event(QMouseEvent(QMouseEvent.Type.MouseButtonPress, pos, QPointF(tile.mapToGlobal(QPoint(int(pos.x()), int(pos.y())))),
                           Qt.LeftButton, Qt.LeftButton, Qt.NoModifier))


def test_panel_boxes_check_design_perform_and_golden_gate_vector():
    """The panel starts with gray boxes and follows the fragment count; clicking a gray box asks for a source; it checks, designs and performs automatically; Golden Gate vector and enzyme."""
    # --- gray boxes and the fragment count ---
    panel = P.AssemblyPanel(Win())
    assert panel.method.count() == 3 and len(panel.slots) == 3 and panel.slots == [None] * 3
    assert "still empty" in panel.status.text() and not panel.buttons["perform"].isEnabled()
    panel.count.setValue(5)
    assert len(panel.slots) == 5 and len(panel.strip.slots) == 5
    panel.count.setValue(1)
    assert len(panel.slots) == 1
    assert panel.count.minimum() == 1 and panel.count.maximum() == P.MAX_FRAGMENTS
    panel.strip.resize(700, panel.strip.minimumHeight())
    panel.strip.render(QPixmap(panel.strip.size()))      # the gray boxes paint

    # --- clicking a gray box ---
    panel = P.AssemblyPanel(Win())
    asked = []
    panel.strip.source_asked.disconnect(panel._ask_source)   # the real one opens a menu and waits for a click
    panel.strip.source_asked.connect(lambda i, pos: asked.append(i))
    panel.strip.resize(700, panel.strip.minimumHeight())
    panel.strip.render(QPixmap(panel.strip.size()))
    click(panel.strip, panel.strip.boxes[1].center().x())
    assert asked == [1] and panel.strip.current == 1
    panel = filled_panel()
    asked.clear()
    panel.strip.source_asked.disconnect(panel._ask_source)
    panel.strip.source_asked.connect(lambda i, pos: asked.append(i))
    panel.strip.resize(900, panel.strip.minimumHeight())
    panel.strip.render(QPixmap(panel.strip.size()))
    box = panel.strip.boxes[2]
    click(panel.strip, box.center().x())
    assert asked == [] and panel.strip.current == 2 and "term" in panel.which.text(), "selects it and shows its file below"
    click(panel.strip, panel.strip._source_rect(box).center().x(), y=panel.strip._source_rect(box).center().y())
    assert asked == [2], "the ▾ of a filled box asks for another source"

    # --- checks automatically, designs and performs ---
    panel = filled_panel(0)
    assert "do not fit" in panel.status.text() and not panel.buttons["perform"].isEnabled()
    assert [c.level for c in panel.strip.cards] == ["bad"] * 3
    panel.add_homology()
    assert panel.design.ok and "checks out" in panel.status.text() and panel.buttons["perform"].isEnabled()
    assert "bad" not in [c.level for c in panel.strip.cards] and len(panel.design.rows) == 6
    panel.perform()
    assert panel.win.opened.length == 1800 and panel.win.opened.circular
    panel.strip.render(QPixmap(panel.strip.size()))
    # Reordering or flipping breaks it again, and the panel says so at once.
    panel._flip_at(1)
    assert not panel.design.ok and not panel.buttons["perform"].isEnabled()
    panel.add_homology()
    assert panel.design.ok
    panel._move_to(2, 0)   # a rotation of the circle: every junction still fits
    assert [f.name for f in panel.slots] == ["term", "vector", "insert"] and panel.design.ok
    panel._move_to(1, 0)   # a real reordering does not
    assert [f.name for f in panel.slots] == ["vector", "term", "insert"] and not panel.design.ok

    # --- Golden Gate vector and enzyme ---
    panel = filled_panel(1)
    panel.slots[0].src.circular = True
    panel.strip.current = 0
    panel.toggle_vector()
    assert panel.slots[0].vector and "needs exactly 2" in panel.status.text(), "the plain vector has no sites"
    v = Document.new_sequence("vec", vector_seq(500), True)
    panel._set_source(0, A.snapshot(v), "vec.gb", None)
    panel.strip.current = 0
    assert panel.slots[0].vector, "a new source keeps the fragment's role"
    panel.add_homology()
    assert panel.design.ok and panel.buttons["perform"].isEnabled()
    panel._show_product()
    ov = [f for f in panel.product_view.editor.doc.features() if "overhang" in f.label]
    assert len(ov) == 3 and panel.product_view.editor.doc.sequence(0, None).upper() == panel.design.product()
    panel.enzyme.setCurrentText("BsmBI / Esp3I")
    assert panel.design is None and "0 BsmBI site(s)" in panel.status.text(), "the vector has no BsmBI sites"
    assert not panel.buttons["perform"].isEnabled()


def test_panel_sequence_view_selection_pcr_and_saving(monkeypatch, tmp_path):
    """The panel uses the program's sequence view read-only, selects a stretch in it and uses it, makes a PCR fragment from two primers (adding the tail), and saves fragments and primers."""
    # --- the sequence view is read-only ---
    from PySide6.QtGui import QKeyEvent
    panel = filled_panel(0)
    panel.strip.current = 1
    panel._show_source(1)
    view, ed = panel.view, panel.sources[1].editor
    assert isinstance(view, P.SequenceView) and view.editor is ed and panel.seqstack.currentWidget() is view
    before, version = ed.doc.sequence(), ed.version
    for key, mods, text in ((Qt.Key_A, Qt.NoModifier, "a"), (Qt.Key_Backspace, Qt.NoModifier, ""), (Qt.Key_Delete, Qt.NoModifier, ""),
                            (Qt.Key_X, Qt.ControlModifier, "x"), (Qt.Key_V, Qt.ControlModifier, "v")):
        view.keyPressEvent(QKeyEvent(QKeyEvent.Type.KeyPress, key, mods, text))
    ed.select(5, 9)
    view.keyPressEvent(QKeyEvent(QKeyEvent.Type.KeyPress, Qt.Key_Right, Qt.NoModifier, ""))   # moving is allowed
    view.mousePressEvent(QMouseEvent(QMouseEvent.Type.MouseButtonPress, QPointF(40, 40), QPointF(40, 40), Qt.RightButton,
                                     Qt.RightButton, Qt.NoModifier))   # no editing menu
    assert ed.doc.sequence() == before and ed.version == version, "nothing can be typed, cut, pasted or deleted"
    assert ed.start == ed.end == 9
    view.resize(900, 500)
    view.show()
    view.viewport().repaint()

    # --- a stretch selected in the view ---
    panel = filled_panel(0)
    panel.strip.current = 1
    panel._show_source(1)
    ed = panel.sources[1].editor
    ed.select(10, 30)
    assert "Selected 11–30 (20 bp)" in panel.chosen.text() and panel.by_region.isChecked()
    panel.use_selection()
    f = panel.slots[1]
    assert (f.start, f.end) == (10, 30) and len(f.seq) == 20 and f.name == "insert 11-30"
    assert panel.strip.current == 1 and "insert.gb 11" in panel.chosen.text()
    panel.use.click()
    assert panel.slots[1].name == "insert 11-30", "the same stretch keeps its name"
    panel._show_source(1)
    assert (ed.start, ed.end) == (10, 30), "coming back to a fragment shows its stretch selected"
    ed.select(5, 5)
    panel.use_selection()
    assert "Select a stretch" in panel.chosen.text()

    # --- PCR with two primers adds the tail ---
    from primer_db import db
    template = dna(1000)
    doc = Document.new_sequence("tmpl", template)
    db().add("F1", "GGTCTCAAATG" + template[200:224])
    db().add("R1", "TTT" + rc(template[600:624]))
    panel = P.AssemblyPanel(Win([doc]))
    panel._set_source(0, A.snapshot(doc), "tmpl.gb", None)
    ed = panel.sources[0].editor
    fw, rv = [s for s in ed.primer_sites if s.label == "F1"][0], [s for s in ed.primer_sites if s.label == "R1"][0]
    assert (fw.strand, rv.strand) == (1, -1) and fw.tail > 0
    assert panel.forward.count() == 1 and panel.reverse.count() == 1 and panel.by_primers.isEnabled()
    ed.select_primer(fw.key)     # click a primer in the view...
    assert "shift-click" in panel.chosen.text()
    ed.select_primer(rv.key, extend=True)   # ...and shift-click the other: the stretch they amplify
    assert panel.by_primers.isChecked() and "F1 + R1" in panel.chosen.text()
    panel.use_selection()
    f = panel.slots[0]
    assert (f.start, f.end) == (fw.start, rv.end) and f.fwd[0] == "F1" and f.rev[0] == "R1"
    tf = ("GGTCTCAAATG" + template[200:224])[:-(fw.end - fw.start)]
    assert f.seq.startswith(tf + template[fw.start:fw.start + 5]) and f.tails[0] == len(tf)
    assert f.seq.endswith(rc(f.rev[1][:len(f.rev[1]) - (rv.end - rv.start)])) and f.tails[1] == len(f.rev[1]) - (rv.end - rv.start)
    assert "PCR with F1 + R1" in panel.chosen.text()
    # The primers of the file itself are found too, without being in the database.
    d2 = Document.new_sequence("withprimer", template)
    d2.add_features([("primer_bind", "300..323", [("label", "mine"), ("note", "sequence: " + template[299:323])])])
    p2 = P.AssemblyPanel(Win())
    p2._set_source(0, A.snapshot(d2), "withprimer.gb", None)
    assert [s.label for s in p2.sources[0].editor.primer_sites if s.label == "mine"] == ["mine"]
    assert p2.forward.count() >= 1

    # --- saves fragments and primers ---
    panel = filled_panel(0)
    panel.add_homology()
    out = tmp_path / "frags"      # apart from the primer database file in tmp_path
    out.mkdir()
    monkeypatch.setattr(P.QFileDialog, "getExistingDirectory", lambda *a, **k: str(out))
    assert panel.save_fragments() == 3
    saved = sorted(p.name for p in out.iterdir())
    assert saved == ["insert.gb", "term.gb", "vector.gb"]
    back = Document(str(out / "insert.gb"))
    assert back.length == len(panel.slots[1].seq) and back.sequence(0, None).upper() == panel.slots[1].seq
    panel.save_primers()
    assert "added" in panel.status.text()


def test_assembly_is_a_tab_opened_from_the_menu():
    import main as M
    import workspace
    assert "Assembly" not in workspace.VIEWS, "not in the drop-down of a pane"
    win = M.MainWindow()
    win.show()
    assert win.assembly is None
    win.show_assembly()                                  # needs no open file: a tab on its own
    assert win.tool_host is not None and win.tool_host.tabText(0) == "Assembly" and win.tool_host.currentWidget() is win.assembly
    win.tool_host.close_tab(0)
    assert win.tool_host is None and win.assembly.parent() is None, "closing the tab keeps the panel for next time"
    win.assembly.count.setValue(2)
    win.show_assembly()
    # Opening a file moves the tab into the workspace, next to the others.
    path = Path(tempfile.mkdtemp()) / "p.gb"
    Document.new_sequence("p", dna(500)).save_genbank(path)
    win.load(str(path))
    ws = win.workspace
    assert win.tool_host is None and ws.count() == 2 and ws.tabText(1) == "Assembly" and ws.currentWidget() is win.assembly
    assert [ws.is_layout(i) for i in range(2)] == [True, False]
    assert len(ws.dump()) == 1, "the tool tab is not part of the saved layouts"
    assert win.assembly.count.value() == 2
    ws.setCurrentIndex(0)
    assert win.sequence_view is not None and ws.currentIndex() == 0, "the other tabs are all still there"
    ws.setCurrentIndex(1)
    assert win.sequence_view is not None, "the editor's views are reachable while the tab is shown"
    ws.duplicate_tab(1)
    assert ws.count() == 2 and ws.can_close(1) and not ws.can_close(0), "the only tab of panes stays"
    assert [n for n, _ in win.open_documents()] == ["p.gb"]
    from sourcefiles import recent_files
    assert str(path.resolve()) in recent_files()
    win.assembly._open_path(0, str(path))
    assert win.assembly.slots[0].name == "p" and win.assembly.sources[0].doc.length == 500
    win.assembly.count.setValue(1)
    win.assembly.add_homology()
    assert win.assembly.design.ok, "one fragment closes on itself"
    win.assembly.perform()
    # Another file: the tab and its content survive the new workspace.
    win.load(str(path))
    assert win.workspace.indexOf(win.assembly) >= 0 and win.assembly.slots[0].name == "p"
    ws = win.workspace
    ws.close_tab(ws.indexOf(win.assembly))
    assert ws.count() == 1 and win.assembly.slots[0].name == "p"
    win.show_assembly()
    assert ws.count() == 2 and ws.currentWidget() is win.assembly and win.assembly.slots[0].name == "p"
    win.close()


# ---- the cards, the primers on the sequence ---------------------------------------------------------------------------


def test_every_method_puts_its_junction_info_on_a_card_and_the_map_shows_them():
    """Every method puts its junction info on a card between the fragments; a vector with a palindromic overhang is flagged; the map has no tabs and shows cards in place of diamonds."""
    # --- cards of every method ---
    g = frags()
    A.auto_primers(g, 0, True)
    gib = A.junction_cards(A.analyse(g, 0, True))
    assert len(gib) == 3 and [c.closing for c in gib] == [False, False, True]
    for c, j in zip(gib, A.analyse(g, 0, True).joints):
        assert j.seq in c.lines and j.seq.translate(A.COMP) in c.lines, "the joining bases on both strands"
        assert "bp homology" in c.title and "°C" in c.title and "GC" in c.title
        assert [n for n, _, _ in c.chips] == ["length", "Tm", "GC", "unique", "hairpin", "distinct"]
        assert c.brief[0] == f"{len(j.seq)} bp" and "°C" in c.brief[1] and "%" in c.brief[1], "the two lines of the compact pill"
    h = [vector_frag(), frag("i", dna(500)), frag("t", dna(300))]
    A.auto_primers(h, 1, True)
    gg = A.junction_cards(A.analyse(h, 1, True))
    assert all("overhang" in c.title and c.lines[0].startswith("…") and "|" in c.lines[0] for c in gg)
    assert all(c.brief[0].endswith(" bp overhang") for c in gg), "a pill names what joins when it is not homology"
    assert [n for n, _, _ in gg[0].chips] == ["ends fit", "unique", "not palindromic", "orthogonal", "GC mix", "validated"]
    sw = A.junction_cards(A.analyse([frag("a", dna(100)), frag("b", dna(100)), frag("c", dna(100))], 2, False))
    assert len(sw) == 2 and all("toehold" in c.title and c.lines[-1].startswith("barcode ") for c in sw)
    assert all(c.brief[0].endswith(" nt toehold") for c in sw)
    assert [n for n, _, _ in sw[0].chips] == ["GC", "unique", "distinct"]
    # An end that does not fit is a red card that says why.
    red = A.junction_cards(A.analyse(frags(), 0, True))
    assert all(c.level == "bad" and "do not fit" in c.title and c.chips[0][1] == "bad" and c.text for c in red)
    assert all(c.brief == ("ends do not fit", "") for c in red)

    # --- a palindromic vector overhang is flagged ---
    seq = dna(296) + "AATT" + "T" + "GAGACC" + dna(120) + "GGTCTC" + "A" + "TCGC" + dna(1396)
    f = [A.Frag("vec", Document.new_sequence("vec", seq, True), vector=True), frag("i", dna(700)), frag("t", dna(250))]
    A.auto_primers(f, 1, True)
    d = A.analyse(f, 1, True)
    assert not d.ok and A.junction_cards(d)[0].level == "bad"
    assert ("not palindromic", "bad") in {(s, lv) for s, lv, _ in d.jchecks[0]}

    # --- the map shows a pill for each junction in place of diamonds ---
    panel = filled_panel(0)
    assert not hasattr(panel, "tabs") and not hasattr(panel, "junction_view") and not hasattr(panel, "table")
    assert not [w for w in panel.findChildren(QWidget) if w.__class__.__name__ in ("QTabWidget", "QTableWidget")]
    assert "Gibson et al" not in " ".join(w.text() for w in panel.findChildren(P.QLabel)), "references live in Help only"
    panel.add_homology()
    m = panel.strip
    m.resize(m.minimumWidth(), m.minimumHeight())
    m.render(QPixmap(m.size()))
    assert len(m.card_rects) == 3 and all(r.width() >= m.PILL_MIN and r.height() < m.BOX_H for r in m.card_rects), "compact pills"
    assert m.boxes[0].width() == m.BOX_W <= 130, "fragments are small"
    for i in range(3):   # boxes and junctions alternate, left to right
        assert m.boxes[i].right() < m.card_rects[i].left() and (i == 2 or m.card_rects[i].right() < m.boxes[i + 1].left())
    assert panel.map_area.height() >= m.minimumHeight()
    # Clicking a junction opens its detail.
    r = m.card_rects[1]
    m.event(QMouseEvent(QMouseEvent.Type.MouseButtonPress, r.center(), QPointF(m.mapToGlobal(QPoint(int(r.center().x()), 0))),
                        Qt.LeftButton, Qt.LeftButton, Qt.NoModifier))
    assert len(panel.dialogs) == 1 and "Junction 2" in panel.dialogs[0].windowTitle()
    panel.show_report()
    assert panel.dialogs[-1].windowTitle() == "Assembly Report"


def test_a_circular_assembly_loops_back_to_fragment_1_and_the_product_is_a_ring():
    """The map says a circle is a circle: an arrow runs from the last junction back into Fragment 1 (even before the boxes are
    filled), each junction is a compact pill whose verdict names what is wrong, and the product beside the map is a ring built from
    the fragments' own lengths and colours (a bar when linear). Ticking Circular moves nothing below the map."""
    # --- the arrow, from the first gray boxes on ---
    panel = P.AssemblyPanel(Win())
    m, tile = panel.strip, panel.product_tile
    assert panel.circular.isChecked() and m.loop is not None and m.boxes[0].left() >= m.LOOP_X
    leave, y, right, top, left = m.loop
    assert leave == m.boxes[-1].right() and top < m.boxes[0].top() and left < m.boxes[0].left(), "leaves the last box, runs above, enters Fragment 1"
    assert tile.circular and not tile.ready() and "once every junction fits" in tile.toolTip()
    m.render(QPixmap(m.size()))
    tile.render(QPixmap(tile.size()))      # a dashed gray ring
    height = (m.height(), panel.map_area.height(), tile.height())
    # --- junctions that do not fit: red pills that say so ---
    panel = filled_panel(0)
    m, tile = panel.strip, panel.product_tile
    assert [c.level for c in m.cards] == ["bad"] * 3 and m.loop is not None
    main, second, level = m._pill_lines(m.cards[0])
    assert "do not fit" in main and second == "click for why" and level == "bad"
    assert m.cards[2].closing and m._pill_lines(m.cards[2])[0].startswith("↻"), "the last junction closes the circle"
    assert (m.height(), panel.map_area.height(), tile.height()) == height, "one height whatever the junctions"
    # --- fine: pills with Tm and GC (a check that is only a warning says which one), a ring with a dot on every junction ---
    panel.add_homology()
    levels = [c.level for c in m.cards]
    assert panel.design.ok and "bad" not in levels and set(levels) <= {"ok", "warn"}
    for card in m.cards:
        main, second, level = m._pill_lines(card)
        assert main.endswith(" bp") and level == card.level, (main, second)
        if card.level == "ok":
            assert main.startswith(("✓", "↻ ✓")) and "°C" in second and "%" in second, (main, second)
        else:
            failed = m._failed(card)
            assert failed and second == " · ".join(failed[:3]) and main.startswith(("?", "↻ ?")), (main, second)
    assert (m.height(), panel.map_area.height(), tile.height()) == height
    m.mouseMoveEvent(QMouseEvent(QMouseEvent.Type.MouseMove, m.card_rects[0].center(), QPointF(m.mapToGlobal(QPoint(0, 0))),
                                 Qt.NoButton, Qt.NoButton, Qt.NoModifier))
    assert "homology" in m.toolTip() and "Tm" in m.toolTip(), "the tooltip is the rest of the card in words"
    assert tile.circular and tile.ready() and tile.text.endswith("circular")
    d = panel.design
    assert [c for _, c in tile.parts] == [A.FRAG_COLORS[i] for i in range(3)] and tile.levels == levels
    assert sum(b for b, _ in tile.parts) == len(d.product()), "the ring is the whole product, each fragment its own share of it"
    ends = tile._ends()
    assert ends == sorted(ends) and abs(ends[-1] - 1.0) < 1e-9 and ends[0] > 0
    for selected in (False, True):
        tile.selected = selected
        tile.render(QPixmap(tile.size()))
    m.resize(m.minimumWidth(), m.minimumHeight())
    m.render(QPixmap(m.size()))
    # --- not a circle: no arrow, a bar, and still the same height ---
    panel.circular.setChecked(False)
    assert m.loop is None and not tile.circular and tile.ready() and tile.text.endswith("linear") and len(tile.levels) == 2
    assert (m.height(), panel.map_area.height(), tile.height()) == height
    m.render(QPixmap(m.size()))
    tile.render(QPixmap(tile.size()))
    # --- a very short fragment is still seen on the ring ---
    tile.refresh(True, [(1, "#4fc3f7"), (5000, "#9ccc65")], ["ok", "ok"], "5,001 bp circular")
    assert tile._ends()[0] >= tile.MIN_SHARE * 0.9
    tile.render(QPixmap(tile.size()))
    # --- the map takes all the room up to the product; a scroll bar only when the fragments do not fit ---
    docs = [source_doc("vector", 900), source_doc("insert", 600), source_doc("term", 300)]
    panel = P.AssemblyPanel(Win(docs))
    panel.resize(1500, 900)
    panel.show()
    app.processEvents()
    for i, d in enumerate(docs):   # shown first and filled after, as in the program (the map grows as the junctions come)
        panel._set_source(i, A.snapshot(d), f"{d.name}.gb", None)
    panel.add_homology()
    app.processEvents()
    at = lambda w: w.mapTo(panel, QPoint(0, 0))
    bar = panel.map_area.horizontalScrollBar()
    assert bar.maximum() == 0, "three fragments fit: no scroll bar"
    gap = at(panel.product_tile).x() - (at(panel.map_area).x() + panel.map_area.width())
    assert 0 <= gap <= 8 and panel.map_area.width() > panel.strip.minimumWidth() + 200, "the map stretches all the way to the product"
    panel.count.setValue(12)
    app.processEvents()
    assert bar.maximum() > 0, "twelve do not fit: it scrolls"
    panel.count.setValue(3)
    app.processEvents()
    assert bar.maximum() == 0
    panel.close()


def test_the_product_carries_the_primers_that_made_the_fragments_and_marks_what_they_added():
    """Every primer of the fragments (designed or given) is a primer_bind feature of the product, at the site where it binds
    (a primer across the origin included), so the preview, the file Perform assembly opens and the primer lanes show them; the
    preview marks what each primer's tail added (orange) apart from the fragments (their colours) and the joint (slate)."""
    f, orfs = three_orfs()
    A.auto_primers(f, 0, True, min_ov=20, max_ov=20)
    d = A.analyse(f, 0, True)
    wanted = {(p[0], p[1].upper()) for x in f for p in (x.fwd, x.rev)}
    assert len(wanted) == 6
    for preview in (False, True):
        doc = A.simulate(f, d, preview=preview)
        sites = [x for x in doc.features() if x.kind == "primer_bind"]
        assert {(x.label, P.primer_sequence(x)) for x in sites} == wanted, "each primer once, with its whole sequence in the note"
        assert all(sum(e - s for s, e in x.spans) == len(P.primer_sequence(x)) for x in sites), "the whole primer binds: its tail is the neighbour's own end"
        assert {x.label: x.strand for x in sites} == {name: (1 if name.endswith("-F") else -1) for name, _ in wanted}
    assert any(len(x.spans) == 2 for x in sites), "one of them runs across the origin of the circle"
    # The program draws them as primers: a view of the product looks for them by their bases.
    from editor import Editor
    ed = Editor(doc, "product")
    assert {(n, q.upper()) for n, q in ed.file_primer_set.values()} == wanted
    assert {x.label for x in ed.primer_sites} >= {name for name, _ in wanted}
    # What the primers added, and the colours that keep it apart from the rest.
    assert A.TAIL_COLOR not in A.FRAG_COLORS and A.JOINT_COLOR not in A.FRAG_COLORS and A.TAIL_COLOR != A.JOINT_COLOR
    colour = lambda x: dict(x.qualifiers).get("ApEinfo_fwdcolor")
    tails = [x for x in doc.features() if colour(x) == A.TAIL_COLOR]
    joints_ = [x for x in doc.features() if colour(x) == A.JOINT_COLOR]
    assert len(joints_) == 3 and all(x.label.startswith("homology") for x in joints_)
    assert sum(e - s for x in tails for s, e in x.spans) == 60, "three overlaps of 20 bp, all of them made by tails (10 nt from each primer)"
    # Only a preview marks anything: the product itself is the fragments' own annotations and the primers, and the product view of
    # the panel leaves the tails to the junction inspector.
    plain = A.simulate(f, d)
    assert not [x for x in plain.features() if x.label in d.names or x.label.startswith(("tail", "homology"))], "no fragment, joint or tail marks"
    no_tails = A.simulate(f, d, preview=True, tails=False)
    assert [x for x in no_tails.features() if x.label in d.names] and not [x for x in no_tails.features() if x.label.startswith("tail")]
    assert len([x for x in no_tails.features() if colour(x) == A.JOINT_COLOR]) == 3
    for j in range(3):
        facts = A.junction_facts(d, j)
        assert (facts.length, facts.shared, len(facts.added)) == (20, 0, 2) and sorted(c for *_, c in facts.added) == [10, 10]
        assert {side for _, _, side, _ in facts.added} == {"start", "end"}, "one primer on each side of the junction"
        text = A.junction_summary(d, j)
        assert "20 bp homology" in text and text.count("adds 10 nt") == 2 and "start of Fragment" in text and "end of Fragment" in text
    # Perform assembly opens a product that has them.
    panel = filled_panel(0)
    panel.add_homology()
    panel.perform()
    assert {x.label for x in panel.win.opened.features() if x.kind == "primer_bind"} == {f"{n}-{k}" for n in ("vector", "insert", "term") for k in "FR"}
    assert not [x for x in panel.win.opened.features() if x.label in panel.design.names or x.label.startswith(("tail", "homology"))], \
        "Perform assembly adds no fragment-name, joint or tail annotation"
    # Golden Gate: the primers are found by the bases that are in the product (site and spacer are cut off: a tail hangs off).
    g = [vector_frag()] + three_orfs()[0]
    A.auto_primers(g, 1, True)
    dg = A.analyse(g, 1, True)
    sites = [x for x in A.simulate(g, dg).features() if x.kind == "primer_bind"]
    assert len(sites) == 2 * len(g) - 2 and "overhang" in A.junction_summary(dg, 1)
    # Sidewinder: oligos, not primers (they bind nowhere).
    sw = three_orfs(codons=(30, 40, 35))[0]
    assert not [x for x in A.simulate(sw, A.analyse(sw, 2, False)).features() if x.kind == "primer_bind"]


def test_what_the_fragments_already_share_is_not_credited_to_a_primer():
    src = Document.new_sequence("plasmid", dna(1000))
    f = [A.Frag("left", src, start=0, end=500), A.Frag("right", src, start=470, end=1000)]
    d = A.analyse(f, 0, False)
    facts = A.junction_facts(d, 0)
    assert (facts.length, facts.added, facts.shared) == (30, [], 30)
    assert "comes from the fragments' own bases" in A.junction_summary(d, 0) and "no primer adds anything" in A.junction_summary(d, 0)
    assert not [x for x in A.simulate(f, d, preview=True).features() if x.kind == "primer_bind" or x.label.startswith("tail")]
    # Ten bases that both have and ten more that a primer's tail brings: half is credited to the primer, half is shared.
    bases = src.sequence(0, None)
    f = [A.Frag("left", src, start=0, end=500), A.Frag("right", src, start=490, end=1000, fwd=("right-F", bases[480:512], 22))]
    d = A.analyse(f, 0, False)
    facts = A.junction_facts(d, 0)
    assert (facts.length, facts.added, facts.shared) == (20, [("right-F", 1, "start", 10)], 10)
    assert A.junction_summary(d, 0).endswith("primer right-F adds 10 nt to the start of Fragment 2; the other 10 bp were already shared by the two fragments.")
    tails = [x for x in A.simulate(f, d, preview=True).features() if x.label.startswith("tail")]
    assert [x.spans for x in tails] == [[(480, 490)]], "the tail, where the left fragment's own bases are"


def test_the_junction_inspector_is_the_programs_sequence_view_on_the_top_strand_at_the_joint():
    """Clicking a junction opens the product around the joint in the sequence view with the second strand off: the joint is
    selected, each fragment is a bar, what a primer added is orange, the primers are drawn and answer to the pointer, a sentence says
    what each primer adds, and the checks are listed (the dialog for the junction that closes the circle wraps the origin)."""
    from hover import describe
    panel = filled_panel(0)
    panel.add_homology()
    d = panel.design
    for j in range(3):
        dlg = P.JunctionDialog(panel, panel.slots, d, j)
        view, ed = dlg.view, dlg.view.editor
        n = len(d.names)
        assert isinstance(view, P.ReadOnlyView) and not ed.show_complement and view.strand_h == view.lh, "the top strand alone"
        joint = len(d.joints[j].seq)
        assert (ed.start, ed.end) == (50, 50 + joint) and ed.doc.length == joint + 100, "the joint, selected, with 50 bases each side"
        assert ed.doc.sequence(50, 50 + joint) == d.joints[j].seq, "the bases of the junction"
        labels = {x.label for x in ed.doc.features()}
        assert d.names[j % n] in labels and d.names[(j + 1) % n] in labels and any(x.startswith("homology") for x in labels)
        assert sum(1 for x in ed.doc.features() if x.label.startswith("tail")) == 2, "a tail from each primer"
        assert {panel.slots[j % n].rev[0], panel.slots[(j + 1) % n].fwd[0]} <= {s.label for s in ed.primer_sites}, "the two primers at the joint"
        text = " ".join(w.text() for w in dlg.findChildren(P.QLabel))
        assert f"Junction {j + 1}" in text and "adds 10 nt" in text and "top strand" in text and "of the product" in text
        assert "length" in dlg.findChildren(P.QTextBrowser)[0].toPlainText()
        # The pointer: a primer and a tail answer as they do in every view.
        view.resize(900, 400)
        view.show()
        app.processEvents()
        bpl = view.bases_per_line()
        site = next(s for s in ed.primer_sites if s.start < 50 + joint and s.end > 50 and s.strand > 0)
        row = site.start // bpl
        y = view.offset(row) - view.verticalScrollBar().value() + view._primer_lane_top(row, site.strand, ed.primer_lanes[site.key]) + view.primer_h / 2
        key = view.hover_target(view.left + (site.start + 2 - row * bpl) * view.cw, y)
        assert key == site.key and site.label in describe(ed, key) and "Tm" in describe(ed, key)
        tail = next(x for x in ed.features if x.label.startswith("tail"))
        assert "primer" in describe(ed, tail.index) and "not in its own sequence" in describe(ed, tail.index)
        dlg.close()
    # A junction whose ends do not fit has no product to show: its two ends and the reason.
    bad = filled_panel(0)
    dlg = P.JunctionDialog(bad, bad.slots, bad.design, 0)
    assert dlg.view is None and "do not fit" in " ".join(w.text() for w in dlg.findChildren(P.QLabel))
    assert "ends" in dlg.findChildren(P.QTextBrowser)[0].toPlainText()
    # The panel opens this dialog for a click on the pill of the map.
    panel._show_junction(1)
    assert isinstance(panel.dialogs[-1], P.JunctionDialog) and panel.dialogs[-1].view is not None
    # The window is cut at the ends of a linear product, and wraps the origin for the junction that closes a circle.
    lin = filled_panel(0)
    lin.circular.setChecked(False)
    lin.add_homology()
    doc = A.simulate(lin.slots, lin.design, preview=True)
    window, at, first = A.junction_window(doc, lin.design, 0)
    assert first >= 0 and at == lin.design.joints[0].pos - first and window.length <= len(lin.design.joints[0].seq) + 100
    n = len(d.product())
    wrap, at, first = A.junction_window(A.simulate(panel.slots, d, preview=True), d, 2)
    assert at == 50 and first == (d.joints[2].pos % n - 50) % n and wrap.length == len(d.joints[2].seq) + 100


def test_designed_primers_are_temporary_until_saved_and_explain_themselves_on_hover():
    from hover import describe
    panel = filled_panel(0)
    panel.strip.current = 1
    panel._show_source(1)
    ed = panel.sources[1].editor
    mine = lambda f: {s.label: s for s in ed.primer_sites if s.label in (f.fwd[0], f.rev[0]) and s.sequence in (f.fwd[1], f.rev[1])}
    assert not [s for s in ed.primer_sites if s.label.startswith("insert-")]
    panel.add_homology()
    f = panel.slots[1]
    sites = mine(f)
    assert set(sites) == {f.fwd[0], f.rev[0]} and all(s.temporary for s in sites.values()), "designed, not in the database: temporary"
    assert (sites[f.fwd[0]].strand, sites[f.rev[0]].strand) == (1, -1)
    assert sites[f.fwd[0]].tail == len(f.fwd[1]) - f.fwd[2] > 0, "the view draws the tail beside the primer"
    assert panel.primer_info.isHidden(), "no box listing the designed primers (only Sidewinder's oligos, which bind nowhere, are listed)"
    # Hovering the primer in the sequence view says what it is.
    view = panel.view
    view.resize(900, 500)
    view.show()
    app.processEvents()
    site = sites[f.fwd[0]]
    row_y = view.offset(0) - view.verticalScrollBar().value()
    y = row_y + view._primer_lane_top(0, site.strand, ed.primer_lanes[site.key]) + view.primer_h / 2
    x = view.left + (site.start + 3) * view.cw
    assert view.hover_target(x, y) == site.key
    tip = describe(ed, site.key)
    assert f.fwd[0] in tip and "temporary primer" in tip and "Tm" in tip and "hairpin" in tip and "Save primers" in tip
    assert f.fwd[1][:site.tail].lower() in tip, "its tail, in lowercase"
    # Saving the primers makes them regular primers: no longer temporary, same place.
    panel.save_primers()
    sites = mine(f)
    assert set(sites) == {f.fwd[0], f.rev[0]} and not any(s.temporary for s in sites.values())
    assert "temporary" not in describe(ed, sites[f.fwd[0]].key) and "primer" in describe(ed, sites[f.fwd[0]].key)
    # Reverse-complemented, the same oligos sit on the other strands (and are new, temporary ones: the old are in the database).
    panel._flip_at(1)
    panel.add_homology()
    f = panel.slots[1]
    sites = mine(f)
    assert (sites[f.fwd[0]].strand, sites[f.rev[0]].strand) == (-1, 1)
    view.viewport().repaint()   # paints the primers, their tails and the selection


def test_actions_menu_pointer_product_box_and_sidewinder_panel():
    """Actions are columns on the right with Perform last and highlighted; the right-click menu of a fragment; a hand pointer over boxes and cards; the product box; the Sidewinder panel shows cards and the oligos."""
    # --- actions are columns on the right ---
    from PySide6.QtWidgets import QPushButton
    panel = filled_panel(0)
    panel.resize(1500, 950)
    panel.show()
    app.processEvents()
    at = lambda w: w.mapTo(panel, QPoint(0, 0))
    b = panel.buttons
    assert set(b) == {"homology", "save", "primers", "perform"}, "no move, reverse-complement or clear buttons: the right-click menu has them"
    assert not [w for w in panel.findChildren(QPushButton) if w.text() in ("◀", "▶", "Reverse complement", "Clear", "Vector")]
    assert b["homology"].text() == "Design primers"
    # The verdict (every junction checks out) and the buttons are a column on the right of the map.
    map_right = at(panel.map_area).x() + panel.map_area.width()
    for w in (panel.status, b["homology"], b["save"], b["primers"], b["perform"]):
        assert at(w).x() >= map_right - 2, w
    ys = [at(b[k]).y() for k in ("homology", "save", "primers", "perform")]
    assert ys == sorted(ys) and at(panel.status).y() < ys[0], "the verdict on top, Perform assembly at the bottom"
    # Perform assembly and Use as fragment are drawn in the highlight colour of the theme.
    assert b["perform"].objectName() == "primary" and panel.use.objectName() == "primary"
    assert "palette(highlight)" in panel.styleSheet() and "font-weight: 700" in panel.styleSheet()
    # What the fragment is made of: a column on the right of the sequence view.
    view_right = at(panel.seqstack).x() + panel.seqstack.width()
    for w in (panel.by_region, panel.by_primers, panel.forward, panel.reverse, panel.use, panel.chosen):
        assert at(w).x() >= view_right - 2, w
    panel.close()

    # --- right-click menu of a fragment ---
    panel = filled_panel(1)
    strip = panel.strip
    asked = []
    strip.menu_asked.disconnect(panel._context_menu)   # the real one opens a menu and waits for a click
    strip.menu_asked.connect(lambda i, pos: asked.append(i))
    strip.resize(1400, strip.minimumHeight())
    strip.render(QPixmap(strip.size()))
    click(strip, strip.boxes[1].center().x(), button=Qt.RightButton)
    assert asked == [1], "a right-click on a box asks for its menu"
    menu = panel._build_context_menu(1)
    texts = [a.text() for a in menu.actions() if a.text()]
    assert texts == ["Source", "Move left", "Move right", "Reverse complement", "Destination vector (cut at its two sites)",
                     "Clear this fragment"]
    act = lambda m, name: next(a for a in m.actions() if a.text() == name)
    assert [a.text() for a in act(menu, "Source").menu().actions()] == ["Currently open files", "Recently opened files", "Browse…"]
    names = lambda: [f.name if f else None for f in panel.slots]
    act(menu, "Move left").trigger()
    assert names() == ["insert", "vector", "term"]
    menu = panel._build_context_menu(0)
    assert not act(menu, "Move left").isEnabled() and act(menu, "Move right").isEnabled(), "the first fragment cannot move left"
    before = panel.slots[0].seq
    act(menu, "Reverse complement").trigger()
    assert panel.slots[0].flipped and panel.slots[0].seq == rc(before)
    vector = act(panel._build_context_menu(0), "Destination vector (cut at its two sites)")
    assert vector.isCheckable() and not vector.isChecked()
    vector.trigger()
    assert panel.slots[0].vector and act(panel._build_context_menu(0), "Destination vector (cut at its two sites)").isChecked()
    act(panel._build_context_menu(0), "Clear this fragment").trigger()
    assert names()[0] is None and panel.sources[0] is None
    # A gray box offers only its source; Gibson has no destination vector.
    assert [a.text() for a in panel._build_context_menu(0).actions() if a.text()] == ["Source"]
    panel.method.setCurrentIndex(0)
    assert "Destination vector (cut at its two sites)" not in [a.text() for a in panel._build_context_menu(1).actions()]

    # --- the pointer is a hand over boxes and cards ---
    panel = filled_panel(0)
    panel.add_homology()
    strip = panel.strip
    strip.resize(1400, strip.minimumHeight())
    strip.render(QPixmap(strip.size()))

    def hover(point):
        strip.mouseMoveEvent(QMouseEvent(QMouseEvent.Type.MouseMove, point, QPointF(strip.mapToGlobal(QPoint(int(point.x()), int(point.y())))),
                                         Qt.NoButton, Qt.NoButton, Qt.NoModifier))
        return strip.cursor().shape()

    assert hover(strip.boxes[0].center()) == Qt.PointingHandCursor
    assert hover(strip.card_rects[0].center()) == Qt.PointingHandCursor
    assert hover(QPointF(strip.width() - 3, 2)) == Qt.ArrowCursor

    # --- the product tile ---
    panel = filled_panel(0)
    strip, tile = panel.strip, panel.product_tile
    strip.resize(1400, strip.minimumHeight())
    strip.render(QPixmap(strip.size()))
    assert tile.text is None and not tile.ready(), "gray until every junction fits"
    click_tile(tile)
    assert not panel.product_shown
    panel.add_homology()
    assert panel.design.ok and tile.text.endswith("circular") and tile.ready()
    click_tile(tile)
    assert panel.product_shown and tile.selected and panel.view is panel.product_view
    doc = panel.product_view.editor.doc
    product = panel.design.product()
    assert doc.sequence(0, None).upper() == product and "Product" in panel.which.text()
    assert not (panel.use.isEnabled() or panel.by_region.isEnabled() or panel.forward.isEnabled()), "the fragment controls rest"
    labels = {f.label: f for f in doc.features()}
    for f in panel.slots:   # one coloured bar per fragment
        assert f.name in labels and ("ApEinfo_fwdcolor", A.FRAG_COLORS[panel.slots.index(f)]) in labels[f.name].qualifiers
    assert not [f for f in doc.features() if f.label.startswith("tail")], "the temporary product view has no 'tail N nt' bars (the junction inspector has)"
    homologies = [f for f in doc.features() if "homology" in f.label]
    assert len(homologies) == 3, "three fragments in a circle: three junctions, the last wraps round the origin"
    assert all(sum(e - s for s, e in f.spans) == len(panel.design.joints[k].seq) for k, f in enumerate(homologies))
    panel.view.viewport().repaint()
    # Changing the design redraws the preview; clicking a fragment goes back to its source.
    panel._flip_at(1)
    assert panel.product_shown or not tile.text, "redrawn, or gone when it can no longer be made"
    strip.render(QPixmap(strip.size()))
    panel._pick(0)
    assert not panel.product_shown and panel.product_view is None and not tile.selected and panel.view is panel.sources[0].view
    assert panel.by_region.isEnabled() and panel.forward.isEnabled() and panel.use.isEnabled(), "the fragment controls are back"

    # --- Sidewinder panel ---
    for seed in range(100):
        r = random.Random(seed)
        panel = P.AssemblyPanel(Win())
        panel.method.setCurrentIndex(2)
        for i in range(3):
            panel._set_source(i, Document.new_sequence(f"s{i}", "".join(r.choice("ACGT") for _ in range(100))), f"s{i}.gb", None)
        if panel.design and panel.design.ok:
            break
    assert panel.design and panel.design.ok and len(panel.strip.cards) == 2 and not panel.strip.circular
    assert not panel.buttons["homology"].isEnabled()
    panel.strip.current = 1
    panel._show_source(1)
    assert not panel.primer_info.isHidden()
    info = panel.primer_info.toPlainText()
    assert "s1-coding" in info and "s1-barcode" in info
    panel._show_product()
    assert panel.product_shown and [f for f in panel.product_view.editor.doc.features() if "toehold" in f.label]
    panel.perform()
    assert panel.win.opened.length == 300
