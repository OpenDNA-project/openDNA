"""The quick local aligner (`opendna_core.align`), the Local Align query text and junction descriptions, and Folder Search.

The aligner is checked on planted queries with exact expected coordinates: strands, indels, max_hits, overhangs, partial
matches, the origin of a circle, and structural-variant queries ([A] + 200 unmatched bases + [B] + [C inverted] + [D]).
Nothing here needs a window.
"""
import os
import random
import sys
import tempfile
from pathlib import Path
from types import SimpleNamespace

import pytest

os.environ["OPENDNA_DATA_DIR"] = tempfile.mkdtemp(prefix="opendna_test_")
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "app"))

import opendna_core as core
from opendna_core import Document

from folder_search import hits_text, search_folder
from local_align import MAX_QUERY, QueryError, clean_query, junction_text, summarize

DEMO = HERE.parent / "demo"
COMP = str.maketrans("ACGTacgt", "TGCAtgca")
COMP_BASE = dict(zip("ACGT", "TGCA"))


def rc(s):
    return s.translate(COMP)[::-1]


def rand_dna(rng, n, alphabet="ACGT"):
    return "".join(rng.choice(alphabet) for _ in range(n))


def pick_unslidable(rng, template, query_len, gap_len=0, insert=""):
    """A position `a` in the template where deleting template[a:a+gap_len] (or inserting `insert` before a) has exactly one
    alignment: the neighbouring bases differ, so the gap cannot slide left or right."""
    while True:
        a = rng.randrange(40, len(template) - 40 - gap_len)
        if gap_len:
            b = a + gap_len
            if template[a - 1] != template[b - 1] and template[a] != template[b]:
                return a
        elif insert[-1] != template[a - 1] and insert[0] != template[a]:
            return a


def test_align_single_queries_coordinates_indels_and_overhangs():
    """Local alignment of one query: reverse strand and coordinates; indels at their exact positions; max_hits, no hit and
    degenerate inputs; a query longer than the template and across the origin; queries hanging off a template end (linear and
    circular), whose ends match but whose middle does not, and an end overlap below the reporting threshold."""
    # --- align query on the reverse strand reports strand and coordinates ---
    rng = random.Random(51)
    t = rand_dna(rng, 3000)
    # a flank whose last base cannot, by chance, extend the alignment on either side
    flank = rand_dna(rng, 49) + next(b for b in "ACGT" if b != t[999] and b != COMP_BASE[t[1300]])
    query = flank + t[1000:1300]
    h = core.align(t, query)[0]
    assert (h[0], h[1], h[2], h[5], h[6]) == (1000, 1300, 1, 50, 350) and h[4] == 1.0
    # the reverse complement of the query: the template coordinates are the same, strand -1, the query part is still given in the
    # query's own orientation, i.e. in the reverse complement of the query that was typed
    h = core.align(t, rc(query))[0]
    assert (h[0], h[1], h[2], h[4]) == (1000, 1300, -1, 1.0)
    assert (h[5], h[6]) == (0, 300), "the 50 base flank is at the end of the reverse-complemented query: not aligned"
    h = core.align(t, flank + rc(t[1000:1300]))[0]
    assert (h[0], h[1], h[2], h[5], h[6]) == (1000, 1300, -1, 50, 350)
    # score: +2 per matching base
    assert core.align(t, t[1000:1300])[0][3] == 600
    # --- align indels are found at their exact positions ---
    rng = random.Random(52)
    t = rand_dna(rng, 3000)
    start, end = 1000, 1400
    window = t[start:end]
    # A 5 base deletion from the query: reported as the template span [a, a + 5) missing from the query
    a = start + pick_unslidable(rng, window, 0, gap_len=5)
    deleted = t[start:a] + t[a + 5:end]
    for query, strand in ((deleted, 1), (rc(deleted), -1)):
        h = core.align(t, query)[0]
        assert (h[0], h[1], h[2]) == (start, end, strand), h
        assert h[9] == [(a, a + 5)] and h[8] == [] and h[7] == [], h
    # A 3 base insertion in the query, before template position p
    ins = "GGG"
    for _ in range(100):
        p = start + pick_unslidable(rng, window, 0, insert=ins)
        break
    inserted = t[start:p] + ins + t[p:end]
    h = core.align(t, inserted)[0]
    assert (h[0], h[1]) == (start, end) and h[8] == [(p, 3)] and h[9] == [] and h[7] == [], h
    # One mismatch: its template position and the base the query has
    q = list(t[start:end])
    pos = start + 150
    wrong = "A" if t[pos] != "A" else "C"
    q[150] = wrong
    h = core.align(t, "".join(q))[0]
    assert h[7] == [(pos, wrong)] and h[8] == [] and h[9] == [] and h[4] == pytest.approx(399 / 400)
    # --- align max hits no hit and degenerate inputs ---
    rng = random.Random(53)
    unit = rand_dna(rng, 80)
    t = rand_dna(rng, 300) + unit + rand_dna(rng, 300) + unit + rand_dna(rng, 300) + rc(unit) + rand_dna(rng, 300)
    hits = core.align(t, unit)
    assert sorted((h[0], h[2]) for h in hits) == [(300, 1), (680, 1), (1060, -1)]
    assert len(core.align(t, unit, max_hits=2)) == 2
    assert len(core.align(t, unit, max_hits=1)) == 1 and core.align(t, unit, max_hits=0) == []
    assert [h[0] for h in core.align(t, unit, max_hits=10)] == [h[0] for h in hits], "same ones every time"
    # the best first: scores never increase
    assert all(a[3] >= b[3] for a, b in zip(hits, hits[1:]))
    # random query: no hit; a query that is too short (fewer than 4 bases) and an empty template: nothing
    assert core.align(t, rand_dna(random.Random(54), 300)) == []
    assert core.align(t, "ACG") == [] and core.align("", unit) == []
    # --- align query longer than the template and across the origin ---
    rng = random.Random(55)
    t = rand_dna(rng, 150)
    longer = rand_dna(rng, 60) + t + rand_dna(rng, 60)
    h = core.align(t, longer)
    assert len(h) == 1 and (h[0][0], h[0][1], h[0][2]) == (0, 150, 1)
    assert (h[0][5], h[0][6]) == (60, 210), "the part of the query that aligned"
    h = core.align(t, rc(longer))
    assert (h[0][0], h[0][1], h[0][2]) == (0, 150, -1)
    # A circle: the query runs through the origin, so end passes the length. On a line, each half is found where it is.
    t = rand_dna(rng, 1000)
    q = t[900:] + t[:100]
    h = core.align(t, q, circular=True)[0]
    assert (h[0], h[1], h[2], h[4]) == (900, 1100, 1, 1.0)
    h = core.align(t, rc(q), circular=True)[0]
    assert (h[0], h[1], h[2]) == (900, 1100, -1)
    # On a line, two hits: q[100:200] is t[0:100] and q[0:100] is t[900:1000] (both score 200, then by position); typed as
    # the reverse complement, the same template places with the query parts mirrored
    lin = core.align(t, q, circular=False)
    assert [(x[0], x[1], x[2], x[5], x[6]) for x in lin] == [(0, 100, 1, 100, 200), (900, 1000, 1, 0, 100)]
    lin = core.align(t, rc(q), circular=False)
    assert [(x[0], x[1], x[2], x[5], x[6]) for x in lin] == [(0, 100, -1, 0, 100), (900, 1000, -1, 100, 200)]
    # a hit at the very start of a circle is reported once (not again as a copy past the end)
    assert len(core.align(t, t[:80], circular=True)) == 1
    # Audit (see docs/TEST_CATALOGUE.md section 8): the tests above cover a query LONGER than the template (flanks on
    # both sides) and a query that runs through the origin. These cover the in-between case: the query starts or ends
    # inside the template and the rest of it hangs off the template's end, on a line (no wrapping possible).
    rng = random.Random(61)
    t = rand_dna(rng, 1000)
    overhang = rand_dna(rng, 60)            # 60 bases that exist nowhere in the template
    # The query's first 40 bases are the template's LAST 40; its last 60 bases hang off the end.
    q = t[-40:] + overhang
    h = core.align(t, q)
    assert len(h) == 1, "one hit, and nothing invented for the overhang"
    assert (h[0][0], h[0][1], h[0][2], h[0][4]) == (960, 1000, 1, 1.0), "ends exactly at the template end, never past it"
    assert (h[0][5], h[0][6]) == (0, 40), "only the query's first 40 bases aligned"
    # The same query reverse-complemented: same template place, strand -1, and the aligned part is reported in the
    # query as typed, i.e. now at its END.
    h = core.align(t, rc(q))
    assert len(h) == 1 and (h[0][0], h[0][1], h[0][2], h[0][4]) == (960, 1000, -1, 1.0)
    assert (h[0][5], h[0][6]) == (60, 100)
    # Mirror: the query's LAST 40 bases are the template's FIRST 40; its first 60 bases hang off the start.
    q = overhang + t[:40]
    h = core.align(t, q)
    assert len(h) == 1 and (h[0][0], h[0][1], h[0][2], h[0][4]) == (0, 40, 1, 1.0), "starts at 0, never before it"
    assert (h[0][5], h[0][6]) == (60, 100)
    h = core.align(t, rc(q))
    assert len(h) == 1 and (h[0][0], h[0][1], h[0][2], h[0][4]) == (0, 40, -1, 1.0)
    assert (h[0][5], h[0][6]) == (0, 40)
    # A query that really continues over the origin (end of the template followed by its start) on a LINE: both
    # halves align, each where it belongs (the longer one first), and they are two hits, not one wrapped hit.
    wrap = t[-40:] + t[:60]
    hits = core.align(t, wrap)
    assert [(x[0], x[1], x[5], x[6]) for x in hits] == [(0, 60, 40, 100), (960, 1000, 0, 40)]
    # --- align overhanging query on a circular template wraps across the origin ---
    rng = random.Random(61)                  # same template and overhang as the linear test above
    t = rand_dna(rng, 1000)
    overhang = rand_dna(rng, 60)
    # On a circle the template continues past its last base with its first one: a query made of the template's end
    # followed by its start aligns completely across the origin, `end` passing the length (1060 > 1000).
    wrap = t[-40:] + t[:60]
    for query, strand in ((wrap, 1), (rc(wrap), -1)):
        h = core.align(t, query, circular=True)
        assert len(h) == 1, "one wrapped hit, not two halves"
        assert (h[0][0], h[0][1], h[0][2], h[0][4]) == (960, 1060, strand, 1.0)
        assert (h[0][5], h[0][6]) == (0, 100), "the whole query aligned"
    # The mirror: the same hit is found from a query that starts in the template's start and the other way round does not matter.
    # An overhang that is NOT the template's start still aligns only its overlap: it starts at the right base, covers the
    # 40 base overlap, and (the template being circular) may only extend by a few chance matches into the wrapped start.
    h = core.align(t, t[-40:] + overhang, circular=True)
    assert len(h) == 1 and h[0][0] == 960 and h[0][2] == 1 and h[0][5] == 0
    assert h[0][1] >= 1000 and h[0][6] >= 40 and h[0][6] < 60, "the random overhang does not align as a whole"
    assert h[0][1] - 1000 <= 12, "past the origin: only chance matches, no long invented stretch"
    # A hit at the start of the circle that is also reachable as a copy after the end is reported once.
    assert len(core.align(t, overhang + t[:40], circular=True)) == 1
    # --- align query whose ends match but whose middle does not ---
    # What the aligner returns when only the extremities of a query match the template. Before the structural-variant work a
    # group of nearby diagonals gave only its best local piece; now every piece above the threshold is its own hit.
    rng = random.Random(62)
    t = rand_dna(rng, 1000)
    # (1) Same diagonal (the middle is a SUBSTITUTION of equal length) and a middle too long to bridge: BOTH ends are reported,
    # the longer one (60 bases at 490..550) first, then the 40 base left end (300..340, a chance base may extend it).
    q = t[300:340] + rand_dna(rng, 150) + t[490:550]
    h = core.align(t, q)
    assert len(h) == 2 and (h[0][0], h[0][1], h[0][2], h[0][4]) == (490, 550, 1, 1.0)
    assert (h[0][5], h[0][6]) == (190, 250)
    assert h[1][0] == 300 and 340 <= h[1][1] <= 342 and h[1][5] == 0 and h[1][2] == 1
    h = core.align(t, rc(q))
    assert len(h) == 2 and (h[0][0], h[0][1], h[0][2]) == (490, 550, -1) and (h[0][5], h[0][6]) == (0, 60)
    assert h[1][0] == 300 and h[1][2] == -1
    # (2) Different diagonals (the middle is an insertion relative to the template: the template jumps ahead): the two ends
    # are two separate clusters, so BOTH are reported, the better one first.
    q = t[300:340] + rand_dna(rng, 150) + t[600:660]
    hits = core.align(t, q)
    assert len(hits) == 2 and hits[0][3] >= hits[1][3]
    assert (hits[0][0], hits[0][1], hits[0][5], hits[0][6]) == (600, 660, 190, 250)
    assert hits[1][0] == 300 and 340 <= hits[1][1] <= 342 and hits[1][5] == 0, "the left end (may extend by a chance base)"
    # (3) A short middle (10 mismatching bases) is cheaper than the ends are worth: one bridged hit over the whole query,
    # identity below 1, with the middle shown as mismatches (and maybe some gaps). 60 mismatching bases are not bridged any
    # more (the aligner never crosses a stretch scoring 40 points below its best, see MAX_DRAWDOWN in align.rs): two pieces.
    q = t[300:340] + rand_dna(rng, 10) + t[350:410]
    h = core.align(t, q)
    assert len(h) == 1 and (h[0][0], h[0][1], h[0][5], h[0][6]) == (300, 410, 0, 110)
    assert h[0][4] < 0.97 and h[0][7] and (h[0][7] + h[0][8] + h[0][9]), "the middle shows as mismatches / gaps"
    q = t[300:340] + rand_dna(rng, 60) + t[400:460]
    h = core.align(t, q)
    got = sorted((x[0], x[5]) for x in h)   # (a chance base may extend the right piece to its left by a base or two)
    assert len(got) == 2 and got[0] == (300, 0) and 396 <= got[1][0] <= 400 and got[1][1] == got[1][0] - 300
    # --- align end overlap below the reporting threshold gives no hit ---
    # A 60 base random query that shares only a few bases with the template's end or start. A hit needs min(1.2 * length, 40)
    # points, i.e. 20 matching bases (+2 each) for a query of 34 bases or more; below that nothing is reported.
    rng = random.Random(61)
    t = rand_dna(rng, 1000)
    overhang = rand_dna(rng, 60)
    for k in (5, 10, 15):
        assert core.align(t, overhang + t[:k]) == [], f"{k} base overlap at the template start"
        assert core.align(t, t[-k:] + overhang) == [], f"{k} base overlap at the template end"
        assert core.align(t, rc(t[-k:] + overhang)) == []
        assert core.align(t, t[-k:] + overhang, circular=True) == []
    # Exactly at the threshold the overlap is a hit (score 2 * 20 = 40), so the limit is where the docstring says.
    h = core.align(t, overhang + t[:20])
    assert [(x[0], x[1], x[3], x[5], x[6]) for x in h] == [(0, 20, 40, 60, 80)]
    h = core.align(t, t[-20:] + overhang)
    assert [(x[0], x[1], x[3], x[5], x[6]) for x in h] == [(980, 1000, 40, 0, 20)]


PIECE, GAP = 100, 200
# Where the four pieces start in the 2000 bp template: (A, B, C, D). C is the one that is inverted in the query.
LAYOUTS = {
    "far apart, same order": (100, 600, 1100, 1600),
    "far apart, shuffled": (1400, 900, 500, 100),
    "close together": (500, 650, 800, 950),
    "A B D collinear (same diagonal, 200 bp gap)": (500, 800, 1500, 1000),
    "A and B touch in the template (a 200 bp insertion)": (500, 600, 1200, 1600),
    "B and D touch in the template": (200, 800, 1300, 900),
}


def scenario(layout, circular=False, seed=1, wrap_d=False):
    """(template, query, expected) for the A-gap-B-C(inverted)-D query.

    `expected` is the list of (qstart, qend, strand, start, end) in query order. The random draw is repeated until no
    piece could be extended by a chance base into its neighbour (otherwise a hit may legitimately be one base longer,
    which would make exact coordinates a coin flip)."""
    rng = random.Random(seed)
    a, b, c, d = layout
    if wrap_d:
        d = 1950                      # D = template[1950:2000] + template[0:50] on a circle
    while True:
        t = rand_dna(rng, 2000)
        tt = t + t[:PIECE]            # lets a slice run over the end
        pieces = [tt[a:a + PIECE], tt[b:b + PIECE], rc(tt[c:c + PIECE]), tt[d:d + PIECE]]
        spacer = rand_dna(rng, GAP)
        query = pieces[0] + spacer + pieces[1] + pieces[2] + pieces[3]
        q0 = [0, PIECE + GAP, PIECE + GAP + PIECE, PIECE + GAP + 2 * PIECE]
        starts, strands = [a, b, c, d], [1, 1, -1, 1]
        ok = True
        for q, s, strand in zip(q0, starts, strands):
            lq, rq = query[q - 1] if q else None, query[q + PIECE] if q + PIECE < len(query) else None
            if strand == 1:   # query base left of the piece pairs with the template base left of it, etc.
                lt, rt = tt[(s - 1) % 2000], tt[(s + PIECE) % 2000]
            else:             # reversed: the left query base pairs with the complement of the template base at the end
                lt, rt = tt[(s + PIECE) % 2000].translate(COMP), tt[(s - 1) % 2000].translate(COMP)
            ok &= (lq is None or lq != lt) and (rq is None or rq != rt)
        if ok:
            expected = [(q, q + PIECE, strand, s, s + PIECE) for q, s, strand in zip(q0, starts, strands)]
            return t, query, expected


TOL = 2   # bases


def summary(hits):
    """Hits as (qstart, qend, strand, start, end), in query order."""
    return sorted((h[5], h[6], h[2], h[0], h[1]) for h in hits)


def assert_pieces(hits, expected, tol=TOL, case=""):
    """Exactly one hit per expected piece (qstart, qend, strand, start, end), none left over. Each hit must COVER its piece
    and may stick out by at most `tol` bases per side: a local alignment may legitimately extend a piece into the unmatched
    neighbour by a few chance matches (even with a small gap), which is not a missing or wrong piece."""
    got = summary(hits)
    assert len(got) == len(expected), f"{case}: {len(got)} hits, wanted {len(expected)}: {got}"
    for g, e in zip(got, sorted(expected)):
        assert g[2] == e[2], f"{case}: strand of {e}: {g}"
        assert e[0] - tol <= g[0] <= e[0] and e[1] <= g[1] <= e[1] + tol, f"{case}: query range of {e}: {g}"
        assert e[3] - tol <= g[3] <= e[3] and e[4] <= g[4] <= e[4] + tol, f"{case}: template range of {e}: {g}"


def test_align_finds_every_piece_of_a_structural_variant_query():
    """Every piece of an A-gap-B-C(inverted)-D query is reported with exact coordinates and strand, for six layouts on lines
    and circles and typed as its reverse complement; a piece over the origin; max_hits and no merging over a gap; a duplicated
    piece; a piece with a few differences next to an exact one."""
    for name in LAYOUTS:
        for circular in (False, True):
            case = f"{name}, circular={circular}"
            t, query, expected = scenario(LAYOUTS[name], circular=circular, seed=len(name))
            hits = core.align(t, query, circular=circular)
            assert_pieces(hits, expected, case=case)
            assert all(h[4] >= 0.95 and h[3] >= 2 * PIECE - 10 for h in hits), f"{case}: each piece is (almost) exact"
            # the same query typed as its own reverse complement: the same pieces, every strand flipped, query coordinates mirrored
            mirrored = core.align(t, rc(query), circular=circular)
            m = len(query)
            assert_pieces(mirrored, [(m - qe, m - qs, -strand, s, e) for qs, qe, strand, s, e in expected], case=case + ", mirrored")
    # --- piece over the origin of a circle ---
    t, query, expected = scenario(LAYOUTS["far apart, same order"], seed=7, wrap_d=True)
    hits = core.align(t, query, circular=True)
    # D is template[1950:2050] = the end followed by the start: `end` passes the template length
    assert expected[3][3:] == (1950, 2050)
    assert_pieces(hits, expected)
    assert max(h[1] for h in hits) >= 2050
    # on a line the same piece is found as two halves (each end of the template), nothing else is lost
    line = core.align(t, query, circular=False)
    assert_pieces(line, expected[:3] + [(500, 550, 1, 1950, 2000), (550, 600, 1, 0, 50)])
    # --- max_hits and no merging over a gap ---
    t, query, expected = scenario(LAYOUTS["A B D collinear (same diagonal, 200 bp gap)"], seed=3)
    assert len(core.align(t, query, max_hits=2)) == 2 and len(core.align(t, query, max_hits=1)) == 1
    assert core.align(t, query, max_hits=0) == []
    hits = core.align(t, query)
    # no hit spans the 200 bp gap: A, B and D (all on one diagonal) are three pieces of 100 bases, not one of 400+
    assert max(h[1] - h[0] for h in hits) <= PIECE + TOL and len(hits) == 4
    assert all(a[3] >= b[3] for a, b in zip(hits, hits[1:])), "best first"
    # --- duplicated piece in the query ---
    # A tandem duplication: the query holds the same 100 bp twice. Both copies align to the one template place.
    rng = random.Random(5)
    t = rand_dna(rng, 2000)
    q = rand_dna(rng, 60) + t[700:800] + rand_dna(rng, 60) + t[700:800] + rand_dna(rng, 60)
    hits = core.align(t, q)
    assert_pieces(hits, [(60, 160, 1, 700, 800), (220, 320, 1, 700, 800)])
    # --- piece with a few differences ---
    rng = random.Random(6)
    t = rand_dna(rng, 2000)
    messy = list(t[300:420])
    messy[30] = "A" if messy[30] != "A" else "C"
    del messy[80:83]
    q = t[900:1000] + rand_dna(rng, 150) + "".join(messy)
    hits = sorted(core.align(t, q), key=lambda h: h[5])
    assert len(hits) == 2
    assert (hits[0][0], hits[0][1], hits[0][2]) == (900, 1000, 1) and hits[0][4] == 1.0
    assert 300 - TOL <= hits[1][0] <= 300 and 420 <= hits[1][1] <= 420 + TOL and hits[1][7] and hits[1][9] and hits[1][4] < 1.0


def hit(qstart, qend, start, end, strand=1):
    return SimpleNamespace(qstart=qstart, qend=qend, start=start, end=end, strand=strand)


def test_local_align_query_cleaning_and_junction_texts():
    """clean_query accepts FASTA, lowercase, RNA, digits and gaps and says what is wrong otherwise; the junction and summary
    texts name insertions, deletions, inversions and duplications."""
    # --- clean_query accepts ---
    assert clean_query(">seq1 description\nacgt acgt\n  ggcc\n") == "ACGTACGTGGCC"
    assert clean_query("1 acgtacgt 9\n11 ggcctt") == "ACGTACGTGGCCTT"            # a GenBank ORIGIN block
    assert clean_query("ACGU UGCA") == "ACGTTGCA"                                  # RNA: U becomes T
    assert clean_query("AC-GT--AC") == "ACGTAC"                                    # alignment gap dashes
    assert clean_query("ACGTNNRYACGT") == "ACGTNNRYACGT"                           # IUPAC codes are DNA
    # several FASTA records: the first one only
    assert clean_query(">a\nAAAACCCC\n>b\nGGGGTTTT\n") == "AAAACCCC"
    assert clean_query("\n\n  >a\r\nACGT\r\nACGT\r\n") == "ACGTACGT"
    # --- clean_query errors ---
    for text, words in (("", "empty"), ("  \n ", "empty"), (">only a header", "only a FASTA header"), ("123 456", "No DNA")):
        with pytest.raises(QueryError, match=words):
            clean_query(text)
    with pytest.raises(QueryError, match="protein"):
        clean_query("MKVLAAGIVGLLLAQPSEEF")                                          # amino acids
    with pytest.raises(QueryError, match="Not DNA"):
        clean_query("ACGT hello 12")                                                 # H, O: not protein-only letters
    with pytest.raises(QueryError, match=r"too long.*200,000"):
        clean_query("A" * (MAX_QUERY + 1))
    assert len(clean_query("ACGT" * (MAX_QUERY // 4))) == MAX_QUERY                 # exactly the limit is fine
    # --- junction and summary texts ---
    n = 5000
    a = hit(0, 100, 1000, 1100)
    assert junction_text(a, hit(100, 200, 1100, 1200), n, False) == "contiguous"
    # whole texts, the distances worked out by hand (query gap = next qstart - qend; template gap = next start - end)
    assert junction_text(a, hit(300, 400, 1100, 1200), n, False) == "insertion in the query: 200 bp unmatched in the query"
    assert junction_text(a, hit(100, 200, 1150, 1250), n, False) == "deletion from the query: 50 bp of the template skipped"
    # 900 - 1100 = -200: the template goes back 200 bp
    assert junction_text(a, hit(100, 200, 900, 1000), n, False) == (
        "template jumps back (duplication or pieces out of order): template back 200 bp")
    # strand change: the closest pair of ends, min(|1500-1100|, |1000-1600|, |1500-1000|, |1600-1100|) = 400
    assert junction_text(a, hit(100, 200, 1500, 1600, -1), n, False) == (
        "inversion breakpoint (strand changes): 400 bp of the template skipped")
    # 50 bp of query between the pieces where 100 bp of template are skipped
    assert junction_text(a, hit(150, 250, 1200, 1300), n, False) == (
        "unmatched query in place of template: 50 bp unmatched in the query; 100 bp of the template skipped")
    # the pieces share 20 query bases, the template runs straight on
    assert junction_text(a, hit(80, 180, 1100, 1200), n, False) == "contiguous: queries overlap by 20 bp"
    # two reverse pieces: the template runs backwards along the query
    r1, r2 = hit(0, 100, 1000, 1100, -1), hit(100, 200, 900, 1000, -1)
    assert junction_text(r1, r2, n, False) == "contiguous"
    # a circle: the short way round the origin is a small gap, not a jump back of nearly the whole template
    assert junction_text(hit(0, 100, 4900, 5000), hit(100, 200, 50, 150), n, True) == (
        "deletion from the query: 50 bp of the template skipped")
    assert summarize([a], n, False) == "1 piece: 1 forward"
    s = summarize([a, hit(300, 400, 1100, 1200), hit(400, 500, 3000, 3100, -1)], n, False)
    assert s.startswith("3 pieces: 2 forward, 1 inverted (reverse complement); gaps in the query between pieces: 200 bp unmatched (#1→#2)")
    assert "insertion in the query (#1→#2)" in s and "inversion breakpoint (#2→#3)" in s
    assert summarize([], n, False) == ""


def test_folder_search_finds_a_query_in_its_file_record_place_and_strand(tmp_path):
    """Folder Search: a piece of a demo plasmid is found in its file at its place on both strands; random DNA finds nothing; a
    hit in a multi-record file names its record; an unreadable file is skipped and counted; cancelling stops the search."""
    text = Document(str(DEMO / "pFG051_Chlor.gb")).sequence()
    hits, skipped, files = search_folder(DEMO, text[1000:1060])
    top = [h for h in hits if h.path.name == "pFG051_Chlor.gb"][0]
    assert (top.start, top.end, top.strand, top.identity) == (1000, 1060, 1, 1.0)
    assert not skipped and files >= 5
    flipped, _, _ = search_folder(DEMO, rc(text[1000:1060]))
    top = [h for h in flipped if h.path.name == "pFG051_Chlor.gb"][0]
    assert (top.start, top.end, top.strand) == (1000, 1060, -1)
    # --- random DNA finds nothing ---
    random.seed(4171)   # not a seed any generated demo file uses: seed(7) is demo/msa_demo.fasta's own first 60 bases
    query = "".join(random.choice("ACGT") for _ in range(60))
    assert search_folder(DEMO, query)[0] == []
    # --- a multi-record file and an unreadable one ---
    random.seed(3)
    a, b = ("".join(random.choice("ACGT") for _ in range(2000)) for _ in range(2))
    (tmp_path / "many.fa").write_text(f">first\n{a}\n>second\n{b}\n")
    (tmp_path / "broken.gb").write_text("this is not a GenBank file\n")
    hits, skipped, files = search_folder(tmp_path, b[500:560])
    assert (files, skipped) == (2, 1)
    assert [(h.record, h.start, h.end, h.strand) for h in hits] == [("second", 500, 560, 1)]
    assert "many.fa\tsecond\t501-560\tforward\t100.0%\t60" in hits_text(hits)
    # --- cancelling ---
    (tmp_path / "a.fa").write_text(">a\nACGT\n")
    assert search_folder(tmp_path, "ACGTACGTACGT", cancelled=lambda: True)[0] == []
    # a query that IS found (('second', 500, 560, +1) above) finds nothing when cancelled before the first file
    assert search_folder(tmp_path, b[500:560], cancelled=lambda: True) == ([], 0, 3)
    # cancelled once the first file is done: no second file is read
    calls = []
    search_folder(tmp_path, b[500:560], progress=lambda done, total: calls.append(done), cancelled=lambda: len(calls) >= 1)
    assert calls == [1]
