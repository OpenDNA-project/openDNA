"""Assembly engines: Gibson, Golden Gate, Sidewinder and restriction / ligation. The products are checked base for base,
with the features and primers they must carry; no window, panel or drawing is involved."""
import os
import random
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ["OPENDNA_DATA_DIR"] = tempfile.mkdtemp(prefix="opendna_test_")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "app"))

import pytest

from opendna_core import Document

import assembly as A
import assembly_panel as P
import overhangs as OH
from editor import reverse_complement

rng = random.Random(1)
rc = A.rc


# ---- fragments made of random bases -----------------------------------------------------------------------------------


def dna(n):
    """Random bases without a BsaI / BsmBI site on either strand (a test adds the ones it wants)."""
    while True:
        s = "".join(rng.choice("ACGT") for _ in range(n))
        if not any(site in s for site in ("GGTCTC", "GAGACC", "CGTCTC", "GAGACG")):
            return s


def random_dna(n, seed):
    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))


def frag(name, seq, features=(), **kw):
    doc = Document.new_sequence(name, seq)
    if features:
        doc.add_features([(k, loc, [("label", lab)]) for k, loc, lab in features])
    return A.Frag(name, doc, **kw)


def frags():
    return [frag("A", dna(300)), frag("B", dna(400)), frag("C", dna(250))]


def frags_like(done):
    """Plain copies of fragments that already have primers (their cores)."""
    return [A.Frag(f.name, f.core, vector=f.vector) for f in done]


def vector_seq(tail=1396):
    """A destination vector: the two sites face away from the backbone, and its own end bases (what the first and last
    insert must overhang with) are two good overhangs, GGAG and TCGC."""
    return dna(296) + "GGAG" + "T" + "GAGACC" + dna(120) + "GGTCTC" + "A" + "TCGC" + dna(tail)


def vector_frag():
    return A.Frag("vec", Document.new_sequence("vec", vector_seq(), True), vector=True)


def find_check(d, name):
    return next(c for c in d.checks if c[0] == name)


def find_text(d, name):
    return next(t for n, _, t in d.checks if n == name)


def rotation_of(a, b):
    return len(a) == len(b) and a in b + b


def joined_bases(d, j):
    """What the product holds at the junction's position (reading on round a circle)."""
    p, n = d.product(), len(d.product())
    return "".join(p[q % n] for q in range(j.pos, j.pos + len(j.seq)))


def gg_fragment(name, left, core, right, spacer="a", pad="tt"):
    """A PCR product the way Golden Gate wants it: site + spacer + overhang, the core, overhang + spacer + site."""
    return frag(name, pad + "GGTCTC" + spacer + left + core + right + spacer.replace("a", "t") + "GAGACC" + pad)


def ring(overhangs, cores=None):
    n = len(overhangs)
    cores = cores or [dna(250) for _ in range(n)]
    return [gg_fragment(f"g{i}", overhangs[i], cores[i], overhangs[(i + 1) % n]) for i in range(n)]


# ---- Gibson -------------------------------------------------------------------------------------------------------------


def test_fragment_and_gibson_products():
    """A fragment is a PCR product that flips with its primers; Gibson reports ends that do not fit, fixes them with new
    primers, uses homology the fragments already share, flags a repeated or short one, and every junction sits in the
    product where it says it does."""
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

    # --- ends that do not fit, then new primers ---
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

    # --- homology the fragments already share, and its uniqueness ---
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


# ---- Golden Gate ---------------------------------------------------------------------------------------------------------


def test_golden_gate_products_and_overhang_checks():
    """Golden Gate: a vector and its inserts become one circle, the overhang each end really carries is read, poor overhang
    choices are judged, and the paper's validated set passes every check."""
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


# ---- Sidewinder ----------------------------------------------------------------------------------------------------------


def piece(name, seq, features=()):
    doc = Document.new_sequence(name, seq)
    if features:
        doc.add_features([(k, loc, [("label", lab)]) for k, loc, lab in features])
    return A.Frag(name, doc)


def sidewinder_parts(seed):
    return [piece("a", random_dna(100, seed), [("promoter", "1..50", "pro")]),
            piece("b", random_dna(120, seed + 1000), [("CDS", "10..60", "cds1")]),
            piece("c", random_dna(90, seed + 2000))]


def sidewinder_arms(d):
    """The barcode arms: what each barcode strand carries past the reverse complement of its coding strand."""
    coding = {r[0]: r[1] for r in d.rows if r[0].endswith("-coding")}
    arms = [r[1][len(coding[r[0].replace("-barcode", "-coding")]):] for r in d.rows if r[0].endswith("-barcode")]
    return [a for a in arms if a]


def test_sidewinder_oligos_toeholds_and_barcodes():
    """Sidewinder: the oligos tile the scarless product, a toehold that occurs twice slides until it is unique, the
    barcodes are orthogonal, the product is linear and carries the fragments' features, and the same fragments give the
    same oligos."""
    # --- toeholds, oligos and the product ---
    for seed in range(60):          # a random draw can hold a poor toehold (the checks are right to say so): take one that is fine
        parts = sidewinder_parts(seed)
        d = A.analyse(parts, 2, False)
        if d.ok:
            break
    assert d.ok, "no seed gave a clean Sidewinder design"
    product = "".join(f.seq for f in parts)
    assert d.product() == product and len(product) == 310
    half = A.TOEHOLD // 2
    assert [j.kind for j in d.joints] == ["toehold"] * 2 and [(j.left, j.right) for j in d.joints] == [("a", "b"), ("b", "c")]
    for j, b in zip(d.joints, [100, 220]):
        assert len(j.seq) == A.TOEHOLD
        assert j.seq == product[j.pos:j.pos + A.TOEHOLD]             # taken from the product itself: scarless
        assert abs(j.pos + half - b) <= A.MAX_SHIFT                  # straddles the boundary, or slid a little to dodge a repeat
        assert product.count(j.seq) == 1
    # The coding strands tile the product: each ends with its toehold, the next starts right after it.
    coding = {r[0]: r[1] for r in d.rows if r[0].endswith("-coding")}
    ends = [j.pos + A.TOEHOLD for j in d.joints]
    assert [coding[n] for n in ("a-coding", "b-coding", "c-coding")] == [product[:ends[0]], product[ends[0]:ends[1]], product[ends[1]:]]
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
    # ...while the panel's preview names them.
    marks = {f.label: f for f in A.simulate(parts, d, preview=True).features()}
    assert [(marks[n].spans[0]) for n in ("a", "b", "c")] == [(0, 100), (100, 220), (220, 310)]

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
    assert "".join(r[1] for r in d.rows if r[0].endswith("-coding")) == product, "sliding only moves which oligo carries the bases"

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

    # --- the product is linear even when a circle is asked for, and needs at least two fragments ---
    parts = sidewinder_parts(3)
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
    long = A.analyse(frags(), 2, False)
    assert not long.ok and find_check(long, f"Oligos at most {A.MAX_OLIGO} nt")[1] == "bad", "Sidewinder wants short fragments"


# ---- every ORF keeps its own bases ---------------------------------------------------------------------------------------
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
        assert reads == bases, f"{label} reads {len(reads)} bases in the product instead of {len(bases)}"
        assert dict(cds[label].qualifiers).get("codon_start") == "1", f"{label} lost its reading frame"
        protein = str(Seq(reads).translate())
        assert protein[0] == "M" and protein.endswith("*") and "*" not in protein[:-1], f"{label} does not translate as an ORF"


def test_products_keep_every_orf_on_its_own_bases():
    """Whatever supplies the bases of a junction, every ORF is still annotated over all of its bases: Gibson (linear and
    circular, every frame offset), Golden Gate (the overhang eats into an end), Sidewinder, an ORF on the minus strand and
    one that starts inside an overlap the fragments already share."""
    # --- Gibson with new primers: half of the overlap sits on the neighbour (10 to 13 bases: every frame offset) ---
    for overlap in range(20, 26):
        f, orfs = three_orfs()
        A.auto_primers(f, 0, False, min_ov=overlap, max_ov=overlap)
        d = A.analyse(f, 0, False)
        assert [j.kind for j in d.joints] == ["existing overlap"] * 2 and len(d.product()) == sum(map(len, d.seqs)) - 2 * overlap
        assert_orfs_intact(A.simulate(f, d), orfs)
        # the same, closing the circle: the junction over the origin too
        f, orfs = three_orfs()
        A.auto_primers(f, 0, True, min_ov=overlap, max_ov=overlap)
        d = A.analyse(f, 0, True)
        assert [j.kind for j in d.joints] == ["existing overlap"] * 3 and len(d.product()) == sum(map(len, d.seqs)) - 3 * overlap
        product = A.simulate(f, d)
        assert product.circular
        assert_orfs_intact(product, orfs)

    # --- Golden Gate: the overhang is picked where it scores best, taking 0 to 4 bases off the end of an ORF ---
    for seed in range(6):
        f, orfs = three_orfs(seed)
        g = [vector_frag()] + f
        A.auto_primers(g, 1, True)
        d = A.analyse(g, 1, True)
        assert [j.kind for j in d.joints] == ["overhang"] * 4
        product = A.simulate(g, d)
        assert product.circular
        assert_orfs_intact(product, orfs)

    # --- Sidewinder ---
    f, orfs = three_orfs(codons=(30, 40, 35))
    d = A.analyse(f, 2, False)
    assert len(d.joints) == 2
    assert_orfs_intact(A.simulate(f, d), orfs)

    # --- an ORF on the minus strand: its start sits at the fragment's right end, its stop at the left end ---
    for method, circular in ((0, False), (0, True), (1, True)):
        f, orfs = three_orfs(strands=(1, -1, 1))
        g = ([vector_frag()] if method else []) + f
        A.auto_primers(g, method, circular)
        d = A.analyse(g, method, circular)
        assert all(j.kind != "none" for j in d.joints)
        assert_orfs_intact(A.simulate(g, d), orfs)

    # --- an ORF that starts inside an overlap the fragments already share (unannotated on the left fragment) ---
    (a, orf_a), (b, orf_b), (c, orf_c) = [orf_frag(name, n, 40 + n) for name, n in zip("ABC", (60, 70, 80))]
    for frag_, orf, nxt in ((a, orf_a, orf_b), (b, orf_b, orf_c)):
        doc = Document.new_sequence(frag_.name, orf + nxt[:30])
        doc.add_features([("CDS", f"1..{len(orf)}", [("label", f"{frag_.name} ORF"), ("codon_start", "1")])])
        frag_.src, frag_.end = doc, doc.length
        frag_.rebuild()
    d = A.analyse([a, b, c], 0, False)
    assert [j.kind for j in d.joints] == ["existing overlap"] * 2 and all(len(j.seq) == 30 for j in d.joints)
    assert_orfs_intact(A.simulate([a, b, c], d), {"A ORF": orf_a, "B ORF": orf_b, "C ORF": orf_c})


def test_the_product_carries_the_features_and_the_primers_of_its_fragments():
    """The simulated product keeps every feature once (an overlap's annotation included), carries each primer as a
    primer_bind feature at the site where it binds, and the junction facts say what each primer's tail added and what the
    fragments already shared."""
    # --- a feature inside an existing overlap is annotated once ---
    src = Document.new_sequence("plasmid", dna(1000))
    src.add_features([("misc_feature", "476..495", [("label", "in overlap")]), ("promoter", "51..120", [("label", "left")]),
                      ("terminator", "801..900", [("label", "right")])])
    f = [A.Frag("left half", src, start=0, end=500), A.Frag("right half", src, start=470, end=1000)]
    d = A.analyse(f, 0, False)
    assert d.joints[0].kind == "existing overlap" and len(d.joints[0].seq) == 30 and len(d.product()) == 1000
    product = A.simulate(f, d)
    found = sorted((x.label, x.spans) for x in product.features() if x.label in ("in overlap", "left", "right"))
    assert found == [("in overlap", [(475, 495)]), ("left", [(50, 120)]), ("right", [(800, 900)])]

    # --- the features of every method survive, in a circular product ---
    for method in (0, 1):
        f = ([vector_frag()] if method else [frag("v", dna(500), [("CDS", "50..200", "gene")])]) + [
            frag("i", dna(600), [("promoter", "10..80", "pro")])]
        A.auto_primers(f, method, True)
        d = A.analyse(f, method, True)
        assert d.ok, [(c, d.problems) for c in d.checks if c[1] == "bad"]
        doc = A.simulate(f, d)
        assert [x.label for x in doc.features() if x.label == "pro"] == ["pro"]
        assert doc.circular and doc.length == len(d.product())
    assert A.piece_doc(f, d, 0).length == len(d.seqs[0])

    # --- every primer of the fragments is a primer_bind feature of the product, where it binds ---
    f, orfs = three_orfs()
    A.auto_primers(f, 0, True, min_ov=20, max_ov=20)
    d = A.analyse(f, 0, True)
    wanted = {(p[0], p[1].upper()) for x in f for p in (x.fwd, x.rev)}
    assert len(wanted) == 6
    for preview in (False, True):
        doc = A.simulate(f, d, preview=preview)
        sites = [x for x in doc.features() if x.kind == "primer_bind"]
        assert {(x.label, P.primer_sequence(x)) for x in sites} == wanted, "each primer once, with its whole sequence in the note"
        assert all(sum(e - s for s, e in x.spans) == len(P.primer_sequence(x)) for x in sites), "the whole primer binds"
        assert {x.label: x.strand for x in sites} == {name: (1 if name.endswith("-F") else -1) for name, _ in wanted}
    assert any(len(x.spans) == 2 for x in sites), "one of them runs across the origin of the circle"
    # The program draws them as primers: a view of the product looks for them by their bases.
    from editor import Editor
    ed = Editor(doc, "product")
    assert {(n, q.upper()) for n, q in ed.file_primer_set.values()} == wanted
    assert {x.label for x in ed.primer_sites} >= {name for name, _ in wanted}

    # --- what the primers added, kept apart from the fragments and the joints ---
    assert A.TAIL_COLOR not in A.FRAG_COLORS and A.JOINT_COLOR not in A.FRAG_COLORS and A.TAIL_COLOR != A.JOINT_COLOR
    colour = lambda x: dict(x.qualifiers).get("ApEinfo_fwdcolor")
    tails = [x for x in doc.features() if colour(x) == A.TAIL_COLOR]
    joints_ = [x for x in doc.features() if colour(x) == A.JOINT_COLOR]
    assert len(joints_) == 3 and all(x.label.startswith("homology") for x in joints_)
    assert sum(e - s for x in tails for s, e in x.spans) == 60, "three overlaps of 20 bp, all made by tails (10 nt from each primer)"
    # Only a preview marks anything: the product itself is the fragments' annotations and the primers.
    plain = A.simulate(f, d)
    assert not [x for x in plain.features() if x.label in d.names or x.label.startswith(("tail", "homology"))]
    no_tails = A.simulate(f, d, preview=True, tails=False)
    assert [x for x in no_tails.features() if x.label in d.names] and not [x for x in no_tails.features() if x.label.startswith("tail")]
    assert len([x for x in no_tails.features() if colour(x) == A.JOINT_COLOR]) == 3
    for j in range(3):
        facts = A.junction_facts(d, j)
        assert (facts.length, facts.shared, len(facts.added)) == (20, 0, 2) and sorted(c for *_, c in facts.added) == [10, 10]
        assert {side for _, _, side, _ in facts.added} == {"start", "end"}, "one primer on each side of the junction"
        text = A.junction_summary(d, j)
        assert "20 bp homology" in text and text.count("adds 10 nt") == 2 and "start of Fragment" in text and "end of Fragment" in text
    # Golden Gate: the primers are found by the bases that are in the product; Sidewinder's oligos bind nowhere.
    g = [vector_frag()] + three_orfs()[0]
    A.auto_primers(g, 1, True)
    dg = A.analyse(g, 1, True)
    sites = [x for x in A.simulate(g, dg).features() if x.kind == "primer_bind"]
    assert len(sites) == 2 * len(g) - 2 and "overhang" in A.junction_summary(dg, 1)
    sw = three_orfs(codons=(30, 40, 35))[0]
    assert not [x for x in A.simulate(sw, A.analyse(sw, 2, False)).features() if x.kind == "primer_bind"]

    # --- what the fragments already share is not credited to a primer ---
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
    assert A.junction_summary(d, 0).endswith(
        "primer right-F adds 10 nt to the start of Fragment 2; the other 10 bp were already shared by the two fragments.")
    tails = [x for x in A.simulate(f, d, preview=True).features() if x.label.startswith("tail")]
    assert [x.spans for x in tails] == [[(480, 490)]], "the tail, where the left fragment's own bases are"


# ---- restriction / ligation ------------------------------------------------------------------------------------------

SITES = {"EcoRI": "GAATTC", "BamHI": "GGATCC", "BglII": "AGATCT", "XhoI": "CTCGAG", "EcoRV": "GATATC", "SmaI": "CCCGGG", "PstI": "CTGCAG",
         "NsiI": "ATGCAT", "HindIII": "AAGCTT"}


def clean(n):
    """Random bases with none of the sites these tests cut with, on either strand."""
    while True:
        s = dna(n)
        if not any(site in s or rc(site) in s for site in SITES.values()):
            return s


def ligation(vector, insert, enzymes, features=(), insert_features=()):
    """The fragments (a circular vector cut as a ring, a linear insert) and their design with the two enzymes."""
    v = frag("vec", vector, features, vector=True)
    v.doc.circular = True
    i = frag("ins", insert, insert_features)
    return [v, i], A.analyse([v, i], 3, True, enzymes=enzymes)


def directional():
    """A vector with an EcoRI and a BamHI site 20 bases apart, and an insert cut from EcoRI to BamHI."""
    vec = clean(300) + "GAATTC" + clean(20) + "GGATCC" + clean(400)
    mid = clean(200)
    return vec, clean(5) + "GAATTC" + mid + "GGATCC" + clean(5), mid


def test_restriction_ligation_products():
    """The ligase's exact sequence for a directional cut, for compatible but different enzymes, and for a blunt one; the
    features of the vector and the insert travel into the product, and a real cut-and-ligate of the demo files."""
    # --- the exact product of a directional EcoRI / BamHI ligation ---
    vec, ins, mid = directional()
    (v, i), d = ligation(vec, ins, ("EcoRI", "BamHI"))
    assert d.ok and [j.kind for j in d.joints] == ["overhang", "overhang"] and [j.seq for j in d.joints] == ["AATT", "GATC"]
    # EcoRI G^AATTC and BamHI G^GATCC: the backbone is the long arc, from the BamHI cut round to the EcoRI cut, then the insert as cut.
    expected = vec[327:] + vec[:301] + "AATTC" + mid + "G"
    assert d.product() == expected and d.circular
    assert [joined_bases(d, j) for j in d.joints] == ["AATT", "GATC"], "the overhangs sit where the joints say"
    assert "GAATTC" in expected and "GGATCC" in expected + expected[:5], "both sites come back, one at each end of the insert"
    assert A.simulate([v, i], d).sequence(0, None).upper() == expected
    assert all(c.title.startswith("4 bp overhang") for c in A.junction_cards(d)) and d.cuts[0][0].count("|") == 2

    # --- compatible but different enzymes join, and the hybrid site is cut by neither ---
    vector = clean(300) + "GGATCC" + clean(400)
    insert = clean(5) + "AGATCT" + clean(200) + "AGATCT" + clean(5)   # BglII at both ends, BamHI on the vector
    _, d = ligation(vector, insert, ("BamHI", "BglII"))
    assert d.ok and [j.seq for j in d.joints] == ["GATC", "GATC"]
    assert all({"BamHI", "BglII"} <= set(t.split()) and "5′ GATC overhang" in t for j in range(2) for s, _, t in d.jchecks[j] if s == "ends fit")
    assert all("hybrid junction" in t and "both sites destroyed" in t for j in range(2) for s, _, t in d.jchecks[j] if s == "sites")
    product = d.product() + d.product()[:5]   # a circle: a join can sit over the origin
    assert "GGATCT" in product and "AGATCC" in product and "GGATCC" not in product and "AGATCT" not in product

    # --- a blunt ligation has no overhang and still gives the exact product ---
    blunt_vec = clean(300) + "GATATC" + clean(400)
    mid2 = clean(150)
    (v, i), d = ligation(blunt_vec, clean(5) + "CCCGGG" + mid2 + "CCCGGG" + clean(5), ("EcoRV", "SmaI"))
    assert d.ok and [j.kind for j in d.joints] == ["blunt end", "blunt end"] and A.junction_cards(d)[0].title == "blunt end"
    assert d.product() == blunt_vec[303:] + blunt_vec[:303] + "GGG" + mid2 + "CCC"
    assert "blunt end" in A.junction_html(d) and A.simulate([v, i], d, preview=True).length == len(d.product())

    # --- the features of the vector and the insert travel into the product ---
    # a vector gene that runs over the origin, one in the discarded stretch, and an insert gene on its bottom strand
    feats = [("CDS", "join(700..720,1..30)", "wraps"), ("misc_feature", "310..320", "dropped"), ("CDS", "40..100", "plain")]
    (v, i), d = ligation(vec, ins, ("EcoRI", "BamHI"), feats, [("CDS", "complement(20..120)", "gene")])
    doc = A.simulate([v, i], d)
    found = {x.label: x for x in doc.features()}
    assert "dropped" not in found, "a feature only in the stretch that is thrown away goes with it"
    wraps = sum(e - s for x in doc.features() if x.label == "wraps" for s, e in x.spans)
    for label, src_seq in (("plain", vec[39:100]), ("gene", ins[19:120])):
        s, e = found[label].spans[0]
        assert doc.sequence(s, e).upper() == src_seq.upper(), label
    assert found["gene"].strand < 0 and found["plain"].strand > 0
    assert wraps == 51 and d.product().upper() == doc.sequence(0, None).upper()

    # --- a real cut and ligate from the demo files ---
    demo = Path(__file__).resolve().parent.parent / "demo"
    whole = Document(str(demo / "pBAD30_eilR.gb"))
    seq = whole.sequence(0, None).upper()
    bam, xho = seq.index("GGATCC"), seq.index("CTCGAG")
    v = A.Frag("pBAD30", A.snapshot(Document(str(demo / "pBAD30_RC.gb"))), vector=True)
    cassette = A.Frag("eilR cassette", A.snapshot(whole), bam - 4, xho + 10)   # BamHI to XhoI, cut out as a linear fragment
    d = A.analyse([v, cassette], 3, True, enzymes=("BamHI", "XhoI"))
    assert d.ok, [c for c in d.checks if c[1] == "bad"]
    # The vector is cut at BamHI and XhoI and its long arc (XhoI round the origin to BamHI) kept; the cassette closes it.
    assert d.product().upper() == seq[xho + 1:] + seq[:bam + 1] + seq[bam + 1:xho + 1] and len(d.product()) == len(seq)
    doc = A.simulate([v, cassette], d)
    assert any(x.label == "eilR" and doc.sequence(*x.spans[0]).upper().startswith("ATG") for x in doc.features())
    assert find_check(d, "Vector self-ligation (pBAD30)")[1] == "ok" and find_check(d, "Orientation")[1] == "ok"


def test_restriction_ligation_end_rules_and_warnings():
    """Which ends ligate (bases, kind and strand), the orientation and self-ligation warnings, what a mismatch says, the
    errors for a vector or insert that is not cut twice, and the pieces a digest gives."""
    # --- sticky ends join only when their bases and kind agree; blunt joins blunt ---
    cuts = {"EcoRI": (11, 15), "BamHI": (11, 15), "XhoI": (11, 15), "EcoRV": (13, 13), "SmaI": (13, 13), "PstI": (15, 11), "NsiI": (15, 11)}
    end = lambda name: A.cut_end("A" * 10 + SITES[name] + "C" * 10, 0, name, *cuts[name])
    assert A.ends_fit(end("EcoRI"), end("EcoRI")), "the same sticky end"
    assert not A.ends_fit(end("EcoRI"), end("BamHI")), "AATT is not GATC"
    assert not A.ends_fit(end("EcoRI"), end("EcoRV")), "sticky never joins blunt"
    assert A.ends_fit(end("EcoRV"), end("SmaI")), "blunt joins blunt, whichever enzyme"
    assert A.ends_fit(end("PstI"), end("NsiI")) and end("PstI").kind == "3′", "a 3′ overhang pairs with the same 3′ overhang"
    assert not A.ends_fit(end("PstI"), A.End("x", "5′", end("PstI").seq, 0, 0, "")), "the same bases on the other strand do not pair"

    # --- two identical ends: either orientation, and the user is told ---
    vec = clean(300) + "GAATTC" + clean(400)
    insert = clean(5) + "GAATTC" + clean(200) + "GAATTC" + clean(5)
    _, d = ligation(vec, insert, ("EcoRI", "EcoRI"))
    orientation = find_check(d, "Orientation")
    assert d.ok and orientation[1] == "warn" and "ins" in orientation[2] and "2 different products" in orientation[2]
    assert orientation in d.global_checks, "it is shown on the panel, not hidden on a card"
    directional_vec, directional_ins, _ = directional()
    assert find_check(ligation(directional_vec, directional_ins, ("EcoRI", "BamHI"))[1], "Orientation")[1] == "ok"
    blunt = clean(300) + "GATATC" + clean(400)
    assert find_check(ligation(blunt, clean(5) + "GATATC" + clean(200) + "GATATC" + clean(5), ("EcoRV", "EcoRV"))[1], "Orientation")[1] == "warn"

    # --- a vector whose two ends are compatible can close on itself ---
    check = find_check(d, "Vector self-ligation (vec)")
    assert check[1] == "warn" and "close on itself" in check[2] and "Dephosphorylate" in check[2] and check in d.global_checks
    both = clean(300) + "GGATCC" + clean(20) + "AGATCT" + clean(400)   # BamHI and BglII leave the same ends
    assert find_check(ligation(both, clean(5) + "GGATCC" + clean(100) + "AGATCT" + clean(5), ("BamHI", "BglII"))[1],
                      "Vector self-ligation (vec)")[1] == "warn"
    assert find_check(ligation(directional_vec, directional_ins, ("EcoRI", "BamHI"))[1], "Vector self-ligation (vec)")[1] == "ok"

    # --- ends that do not match say so, and flipping the insert fixes it ---
    (v, i), d = ligation(directional_vec, A.rc(directional_ins), ("EcoRI", "BamHI"))
    assert not d.ok and d.joints[0].kind == "none" and "ends do not match" in d.problems[0] and "EcoRI" in d.problems[0] and "BamHI" in d.problems[0]
    assert A.junction_cards(d)[0].level == "bad" and A.junction_cards(d)[0].brief == ("ends do not fit", "")
    i.flip()
    assert A.analyse([v, i], 3, True, enzymes=("EcoRI", "BamHI")).ok
    _, blunt_design = ligation(clean(300) + "GATATC" + clean(400), clean(5) + "GAATTC" + clean(100) + "GAATTC" + clean(5), ("EcoRV", "EcoRI"))
    assert not blunt_design.ok, "a blunt vector end does not take a sticky insert"
    with pytest.raises(ValueError, match="needs 1 or 2"):
        ligation(clean(500), clean(100), ("EcoRI", "BamHI"))
    with pytest.raises(ValueError, match="cut once"):
        ligation(directional_vec, clean(50) + "GAATTC" + clean(50), ("EcoRI", "BamHI"))

    # --- a digest gives the fragments between the cuts, each holding both whole sites ---
    seq = "CCAAGG" + clean(40) + "GAATTC" + clean(60) + "GGATCC" + clean(30) + "GAATTC" + clean(20) + "GGATCC" + clean(10)
    for circular in (False, True):
        doc = Document.new_sequence("x", seq, circular)
        pieces = A.digest_pieces(doc, ("EcoRI", "BamHI"))
        assert [(p.left.enzyme, p.right.enzyme) for p in pieces] == \
            [("EcoRI", "BamHI"), ("BamHI", "EcoRI"), ("EcoRI", "BamHI")] + [("BamHI", "EcoRI")] * circular
        for p in pieces:   # cut again as a fragment, it gives the same two ends and no further site
            f = A.Frag("f", doc, p.start, p.end)
            kept, _, left, right, extra, _ = A.cut_fragment(f, ("EcoRI", "BamHI"))
            assert (left.enzyme, right.enzyme) == (p.left.enzyme, p.right.enzyme) and not extra and len(kept) == (p.kept[1] - p.kept[0]) % len(seq)
    assert [len(A.cut_fragment(A.Frag("f", Document.new_sequence("x", seq), p.start, p.end), ("EcoRI", "EcoRI"))[0])
            for p in A.digest_pieces(Document.new_sequence("x", seq), ("EcoRI", "EcoRI"))] == [102]


def test_assembly_reports_bad_input_instead_of_failing():
    """A tiny or empty fragment is a bad check, never an exception (every method, circular or linear); an empty stretch of
    a source is a fragment with no bases; a repeated Golden Gate overhang blocks the assembly."""
    # --- a tiny or empty fragment gives a bad check, never an exception ---
    bad_checks = lambda d: [c for c in d.checks if c[1] == "bad"]
    for method, method_name in ((0, "gibson"), (1, "golden gate"), (2, "sidewinder")):
        for circular in (True, False):
            for size in (0, 5, 18):
                case = f"{method_name}, {'circular' if circular else 'linear'}, {size} bp"
                parts = [piece(n, s) for n, s in zip("ABC", [random_dna(300, 1), random_dna(size, 2), random_dna(250, 3)])]
                d = A.analyse(parts, method, circular)                  # as it is: its ends cannot meet its neighbours'
                assert not d.ok and bad_checks(d) and all(text for _, _, text in bad_checks(d)), case
                A.auto_primers(parts, method, circular)                 # designing primers for it must not fail either
                d = A.analyse(parts, method, circular)
                if method == 2 or size == 0:
                    assert not d.ok and bad_checks(d), case             # Sidewinder needs two toeholds; nothing joins an empty fragment
                else:                                                   # tails can give a short one its overlaps: the length warning is left
                    assert dict((n, lv) for n, lv, _ in d.checks)["Fragment length"] == "warn" and "B" in find_text(d, "Fragment length"), case
                if size == 0:
                    assert parts[1].doc.length == 0 and parts[1].seq == "", case

    # --- an empty stretch of a source is a fragment with no bases ---
    src = Document.new_sequence("s", random_dna(100, 6))
    for start, end in ((5, 5), (0, 0), (50, 20)):
        f = A.Frag("e", src, start, end)
        assert f.seq == "" and f.core_seq == "" and f.tails == (0, 0)

    # --- a repeated overhang is a bad check and blocks the assembly ---
    d = A.analyse(ring(["GGAG", "GGAG", "TCGC"]), 1, True, "BsaI")
    flagged = [c for c in d.checks if c[0].endswith("unique") and c[1] == "bad"]
    assert flagged and "same overhang" in flagged[0][2] and not d.ok
    assert A.analyse(ring(["GGAG", "TCGC", "CAGT"]), 1, True, "BsaI").ok   # the same ring with three different ones is fine
