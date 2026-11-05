"""Structural-variant style queries for the quick aligner (`genomeedit_core.align`).

The scenario: a query made of [A: 100 bp of the template] + [200 bp that match nothing] + [B: 100 bp of the template]
+ [C: 100 bp of the template as its REVERSE COMPLEMENT, what a genome inversion gives] + [D: 100 bp of the template].
Every matching piece must be reported on its own, with exact template coordinates, query coordinates and strand.
The template pieces come from different non-overlapping places (far apart, close together, collinear on one diagonal,
touching, wrapped over the origin of a circle). Nothing here needs Qt.
"""
import random
import sys
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "app"))
import genomeedit_core as core

COMP = str.maketrans("ACGT", "TGCA")


def rc(s):
    return s.translate(COMP)[::-1]


def rand_dna(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))


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


def test_every_piece_of_a_structural_variant_query_is_found():
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


def test_alignment_edge_cases():
    """A piece running over the origin of a circle is found with the others (two halves on a line); max_hits keeps the best
    pieces and pieces are not merged over a mismatching gap; a duplicated query piece is reported twice on the same template
    place; a piece with a few differences next to an exact piece is still two pieces."""
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


