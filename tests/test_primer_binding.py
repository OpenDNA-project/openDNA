"""Primer binding (docs/TESTING.md section 3): where a primer anneals, with mismatches, tails, both strands, circular
templates, the annealing extension, primer pairs and their product, the non-specific priming report, a FASTA genome.

The rules, from core/src/primers.rs (what these tests pin down):
- the primer's 3' SEED (15 bases by default, or the whole primer when it is shorter) must match the template with at most
  `max_mismatches` differences, and its last 3 bases must match exactly (the polymerase starts there);
- from the seed, the annealed part is EXTENDED toward the primer's 5' end (+1 per matching base, -2 per mismatch, ending at the
  best score and stopping once 3 below it). Mismatches in that extension are counted but NOT limited by `max_mismatches`;
  what the extension does not reach is a 5' tail that does not bind and is not counted;
- primers shorter than 8 bases, or a seed shorter than 4, never bind; case and U/T say nothing; any other letter (N, IUPAC)
  is simply a base that differs;
- sites are (start, end, strand, mismatches) on the top strand, strand +1 = the primer has the top-strand sequence, `end` may
  pass the length on a circular template; the list follows the scan along the template (by the position of the seed).
"""
import os
import random
import sys
import tempfile
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ["OPENDNA_DATA_DIR"] = tempfile.mkdtemp(prefix="opendna_test_")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "app"))

from PySide6.QtCore import QSettings
from PySide6.QtWidgets import QApplication

from opendna_core import Document, Fasta, find_primers

import editor as editor_module
import primer_db
import primer_panel
import tools
from editor import Editor
from primer_db import add_primer

app = QApplication.instance() or QApplication([])
ROOT = Path(__file__).resolve().parent.parent
rng = random.Random(11)
COMP = str.maketrans("ACGTacgt", "TGCAtgca")
SEED = 15


def rc(s):
    return s[::-1].translate(COMP)


def dna(n, letters="ACGT"):
    return "".join(rng.choice(letters) for _ in range(n))


def other(base):
    """A base that is not `base`."""
    return "ACGT"[("ACGT".index(base) + 1) % 4]


def mutate(seq, positions):
    out = list(seq)
    for k in positions:
        out[k] = other(out[k])
    return "".join(out)


def primer_for(core, strand, positions=()):
    """The primer that anneals to the top-strand bases `core` (strand +1: the same bases; -1: their reverse complement) with the
    bases at `positions` changed. Positions count from the primer's 5' end, as everywhere in the app."""
    return mutate(core if strand > 0 else rc(core), positions)


def sites(template, primer, circular=False, max_mismatches=1, seed=SEED):
    return find_primers(template, circular, [primer], max_mismatches, seed)[0]


def sites_at(template, primer, lo, hi, circular=False, max_mismatches=1, seed=SEED):
    """Only the sites whose 3' end is where the primer was cut from, [lo, hi): a short seed or a high limit also hits other
    places of a random template by chance, and those are not what is being tested."""
    return [x for x in sites(template, primer, circular, max_mismatches, seed)
            if (x[2] > 0 and x[1] == hi) or (x[2] < 0 and x[0] == lo)]


@pytest.fixture(autouse=True)
def private_state(tmp_path, monkeypatch):
    """Never the user's primer database or settings: a private database, and the primer mismatch limit and the primer
    panel's settings read from a file in the temp folder."""
    monkeypatch.setattr(primer_db, "_db", primer_db.PrimerDB(tmp_path / "primers.tsv"))
    monkeypatch.setattr(editor_module, "binding_mismatches", lambda: 2)
    monkeypatch.setattr(primer_panel, "binding_mismatches", lambda: 1)
    ini = QSettings(str(tmp_path / "options.ini"), QSettings.IniFormat)
    monkeypatch.setattr(primer_panel, "settings", lambda: ini)
    monkeypatch.setattr(primer_panel, "genome", lambda: None)


# ---- the rules written again, independently of the Rust code ---------------------------------------------------------------

def reference_sites(template, circular, primer, max_mismatches, seed):
    """The documented rules, in plain Python: same list, same order (a scan along the template)."""
    t = template.upper()
    p = primer.upper().replace("U", "T")
    n = len(t)
    if len(p) < 8 or n == 0 or seed < 4:
        return []
    seed = min(len(p), seed, n)

    def base(i):
        if circular:
            return t[i % n]
        return t[i] if 0 <= i < n else None

    def anneal(i, query, exact_at_the_end):
        """Mismatches of the seed `query` placed at template position i, or None when it does not bind there."""
        bad = 0
        for k, q in enumerate(query):
            b = base(i + k)
            if b is None:
                return None
            if b != q:
                if (k + 3 >= seed if exact_at_the_end else k < 3) or bad == max_mismatches:
                    return None
                bad += 1
        return bad

    def extend(pairs):
        score = best = best_len = bad = best_bad = 0
        for j, pair in enumerate(pairs):
            if pair is None:
                break
            score, bad = (score + 1, bad) if pair[0] == pair[1] else (score - 2, bad + 1)
            if score > best:
                best, best_len, best_bad = score, j + 1, bad
            if best - score > 3:
                break
        return best_len, best_bad

    reach = len(p) - seed
    back = rc(p)
    found = []
    for i in range(n if circular else n + 1 - seed):
        bad = anneal(i, p[len(p) - seed:], True)
        if bad is not None:
            ext, ext_bad = extend(
                [None if base(i - 1 - j) is None else (base(i - 1 - j), p[len(p) - seed - 1 - j]) for j in range(reach)])
            found.append(((i - ext) % n, min(seed + ext, n), 1, bad + ext_bad))
        bad = anneal(i, back[:seed], False)
        if bad is not None:
            ext, ext_bad = extend(
                [None if base(i + seed + j) is None else (base(i + seed + j), back[seed + j]) for j in range(reach)])
            found.append((i, min(seed + ext, n), -1, bad + ext_bad))
    return [(start, start + length, strand, bad) for start, length, strand, bad in found]


def test_the_rust_search_agrees_with_the_rules_written_in_python():
    """Random templates (with repeats), primers cut from them with changes, both strands and topologies, several seeds and limits."""
    r = random.Random(5)
    total = 0
    for _ in range(1500):
        unit = "".join(r.choice("ACGT") for _ in range(r.randint(2, 5)))
        template = "".join(r.choice("ACGT") for _ in range(r.randint(12, 90)))
        template = template[:8] + unit * r.randint(0, 8) + template[8:]
        circular = r.random() < 0.5
        start = r.randrange(len(template))
        primer = (template + template)[start:start + r.randint(6, 30)] if circular else template[start:start + r.randint(6, 30)]
        primer = "".join(r.choice("ACGT") if r.random() < 0.06 else c for c in primer)
        if r.random() < 0.5:
            primer = rc(primer)
        if r.random() < 0.3:
            primer = "".join(r.choice("ACGT") for _ in range(r.randint(0, 6))) + primer
        limit, seed = r.randint(0, 3), r.choice([4, 8, 12, 15, 20])
        expected = reference_sites(template, circular, primer, limit, seed)
        got = sites(template, primer, circular, limit, seed)
        assert got == expected, (template, circular, primer, limit, seed)
        total += len(got)
    assert total > 500, "the random cases must actually contain sites"


def test_document_binding_sites_and_find_primers_agree():
    template = dna(300)
    primer = template[100:124]
    for circular in (False, True):
        doc = Document.new_sequence("t", template, circular)
        for p in (primer, rc(primer), mutate(primer, [12]), "GAATTC" + primer):
            for limit in (0, 1, 2):
                assert doc.binding_sites(p, limit) == sites(template, p, circular, limit)
                assert doc.binding_sites(p, limit, 12) == sites(template, p, circular, limit, 12)


# ---- single mismatches ---------------------------------------------------------------------------------------------------------

def expected_single_mismatch(k, length, start, end, strand, seed=SEED, limit=1):
    """The site of a primer of `length` with its base k (from the 5' end) changed, as the rules say; `strand` +1: the primer
    is bound at [start, end) on its 3' side; -1: the same coordinates, the 5' end of the primer on the right."""
    if k >= length - 3:
        return []                                 # the last 3 bases at the 3' end must match: it never binds
    if k >= length - min(seed, length):           # inside the seed: counted, and limited by `limit`
        return [(start, end, strand, 1)] if limit >= 1 else []
    if k >= 3:                                    # in the extension, crossed when 3 matching bases follow it
        return [(start, end, strand, 1)]
    cut = k + 1                                   # 0, 1, 2: too close to the 5' end to be worth crossing: a tail
    return [(start + cut, end, 1, 0)] if strand > 0 else [(start, end - cut, -1, 0)]


def test_single_mismatch_positions_on_both_strands():
    """Every single-mismatch position is reported with count 1 on both strands (lengths 16-30), and the reverse strand has the same coordinates as the forward case."""
    # --- every single-mismatch position ---
    for strand in [1, -1]:
        for length in [16, 20, 24, 30]:
            template = dna(400)
            core = template[100:100 + length]
            for k in range(length):
                got = sites(template, primer_for(core, strand, [k]))
                assert got == expected_single_mismatch(k, length, 100, 100 + length, strand), (strand, length, k)

    # --- reverse strand sites have the forward coordinates ---
    template = dna(400)
    core = template[100:124]
    for k in (4, 9, 15, 20):
        fwd = sites(template, mutate(core, [k]))
        # the reverse primer that pairs with the same bases, changed at the same template base
        rev = sites(template, rc(mutate(core, [k])))
        assert [(s, e, m) for s, e, _, m in fwd] == [(s, e, m) for s, e, _, m in rev] == [(100, 124, 1)]
        assert (fwd[0][2], rev[0][2]) == (1, -1)


# ---- the limit, the seed, the 3' end ---------------------------------------------------------------------------------------------------


def test_the_limit_counts_mismatches_and_decides_what_binds():
    """Two and three mismatches are counted and the limit decides; a site with exactly the limit in the seed is found and one more is not (both strands)."""
    # --- two and three mismatches, the limit decides ---
    for strand in [1, -1]:
        template = dna(400)
        core = template[100:124]
        for positions in ([12], [10, 14], [9, 15, 19], [10, 12, 14, 16]):
            primer = primer_for(core, strand, positions)
            for limit in range(0, 5):
                got = sites_at(template, primer, 100, 124, False, limit)
                if len(positions) <= limit:
                    assert got == [(100, 124, strand, len(positions))], (strand, positions, limit)
                else:
                    assert got == [], (strand, positions, limit)

    # --- exactly the limit in the seed is found, one more is not ---
    for strand in [1, -1]:
        template = dna(500)
        core = template[200:230]
        for limit in (1, 2, 3):
            for seed in (8, 12, 15, 20):
                inside = list(range(len(core) - seed, len(core) - 3))     # seed bases that may differ
                first = inside[:: max(1, len(inside) // (limit + 1))][:limit]
                assert len(first) == limit
                primer = primer_for(core, strand, first)
                assert sites_at(template, primer, 200, 230, False, limit, seed) == [(200, 230, strand, limit)], (strand, limit, seed)
                too_many = sorted(set(first) | {inside[-1] if inside[-1] not in first else inside[-2]})
                assert len(too_many) == limit + 1
                primer = primer_for(core, strand, too_many)
                assert sites_at(template, primer, 200, 230, False, limit, seed) == [], (strand, limit, seed)


def test_the_seed_and_the_3_prime_end_and_too_short_primers():
    """A mismatch in the seed is rejected for any seed and just outside it is allowed (both strands); a seed shorter than 4 or a primer shorter than 8 never binds."""
    # --- a mismatch in the seed, for any seed ---
    for strand in [1, -1]:
        template = dna(500)
        core = template[200:230]
        length = len(core)
        for seed in range(8, 27):      # (the base just outside the seed must still have 3 more beyond it)
            def primer_with(k):
                return primer_for(core, strand, [k])
            # limit 0: the first base of the seed, and every base of it up to the last 3, rejects the site
            for k in range(length - seed, length - 3):
                assert sites_at(template, primer_with(k), 200, 230, False, 0, seed) == [], (strand, seed, k)
            # the last 3 bases never differ, whatever the limit
            for k in range(length - 3, length):
                assert sites_at(template, primer_with(k), 200, 230, False, 5, seed) == [], (strand, seed, k)
            # just outside the seed the extension crosses it (3 matching bases follow it): found, counted, limit 0
            k = length - seed - 1
            assert sites_at(template, primer_with(k), 200, 230, False, 0, seed) == [(200, 230, strand, 1)], (strand, seed, k)

    # --- too short a seed or primer never binds ---
    template = dna(60)
    assert sites(template, template[10:30], False, 1, 3) == []
    assert sites(template, template[10:30], False, 1, 0) == []
    for n in range(0, 8):
        assert sites(template, template[10:10 + n]) == []
    assert sites("", template[10:30]) == []
    assert sites(template, "") == []


# ---- tails -----------------------------------------------------------------------------------------------------------------------------


def test_tails_are_not_counted_and_a_matching_tail_extends_the_site():
    """A tail of any length 0 to 30 is not counted and does not move the site (both strands, linear and circular); a tail that matches the template extends the site and the count stays right."""
    # --- a non-matching tail of any length ---
    for strand in [1, -1]:
        for circular in [False, True]:
            template = "T" * 40 + dna(24, "ACG") + "T" * 40
            core = template[40:64]
            for length in range(0, 31):
                tail = dna(length, "ACG" if strand > 0 else "CGT")
                primer = tail + (core if strand > 0 else rc(core))
                assert sites(template, primer, circular, 0) == [(40, 64, strand, 0)], (length, strand, circular)
                assert sites(template, primer, circular, 2) == [(40, 64, strand, 0)], (length, strand, circular)

    # --- a matching tail extends the site ---
    for strand in [1, -1]:
        template = dna(300)
        for extra in range(0, 31):
            if strand > 0:
                primer = template[100 - extra:124]               # the "tail" is the bases just upstream: it anneals too
                expected = (100 - extra, 124, 1, 0)
            else:
                primer = rc(template[100:124 + extra])
                expected = (100, 124 + extra, -1, 0)
            assert sites(template, primer, False, 0) == [expected], (strand, extra)
        # the same with one mismatch in the matching tail, crossed because 3+ matching bases follow it
        primer = mutate(template[60:124], [10]) if strand > 0 else rc(mutate(template[100:164], [53]))
        expected = (60, 124, 1, 1) if strand > 0 else (100, 164, -1, 1)
        assert sites(template, primer, False, 0) == [expected], strand


# ---- repeats, palindromes, order ------------------------------------------------------------------------------------------------


def test_repeats_tandem_repeats_and_palindromes_are_returned_once_in_order():
    """Repeated sites are each returned once in order; tandem repeat sites overlap but are never duplicated; a palindromic primer binds both strands at the same place once each."""
    # --- repeated sites ---
    core = dna(24)
    filler = [dna(60, "AT") for _ in range(5)]
    template = filler[0] + core + filler[1] + rc(core) + filler[2] + core + filler[3] + core + filler[4]
    found = sites(template, core, False, 0)
    # positions of the copies, worked out from the pieces
    top1 = len(filler[0])
    bottom = top1 + 24 + len(filler[1])
    top2 = bottom + 24 + len(filler[2])
    top3 = top2 + 24 + len(filler[3])
    assert found == [(top1, top1 + 24, 1, 0), (bottom, bottom + 24, -1, 0), (top2, top2 + 24, 1, 0), (top3, top3 + 24, 1, 0)]
    assert len(set(found)) == len(found)
    assert [s for s, *_ in found] == sorted(s for s, *_ in found)
    # the reverse primer finds the same four places with the strands swapped
    assert sites(template, rc(core), False, 0) == [(s, e, -t, m) for s, e, t, m in found]

    # --- tandem repeats ---
    unit = "GATTACA"
    template = dna(30, "CG") + unit * 8 + dna(30, "CG")
    primer = unit * 3
    found = sites(template, primer, False, 0)
    # the 3' seed fits every 7 bases along the repeat, and each one anneals all 21 bases
    assert found == [(30 + 7 * i, 51 + 7 * i, 1, 0) for i in range(6)]
    assert len(set(found)) == len(found)
    assert found == reference_sites(template, False, primer, 0, SEED)

    # --- palindromic primer ---
    palindrome = "GAATTCAGCTGCAGCTGAATTC"
    assert palindrome == rc(palindrome)
    template = "T" * 30 + palindrome + "T" * 30
    found = sites(template, palindrome, False, 0)
    assert sorted(found) == [(30, 52, -1, 0), (30, 52, 1, 0)]
    # no duplicates, and a random template gives no duplicate or out-of-order scan either
    template = dna(2000)
    primer = template[500:524]
    for circular in (False, True):
        found = sites(template, primer, circular, 3, 8)
        assert len(set(found)) == len(found)
        anchors = [e - 8 if strand > 0 else s for s, e, strand, _ in found]
        assert anchors == sorted(anchors)


# ---- the origin and the ends --------------------------------------------------------------------------------------------------------


def test_sites_across_the_origin_and_mismatches_on_each_side_of_it():
    """A site across the origin on both strands (and after a rotation); a mismatch on each side of the origin."""
    # --- site across the origin ---
    for strand in [1, -1]:
        template = dna(200)
        n = len(template)
        across = template[n - 12:] + template[:12]                 # 24 bases, the origin between the two halves
        primer = across if strand > 0 else rc(across)
        assert sites(template, primer, True, 0) == [(n - 12, n + 12, strand, 0)], strand
        assert sites(template, primer, False, 0) == [], f"a linear template has no such site (strand {strand})"
        # the same site after a rotation: it moves, and no longer crosses the origin
        rotated = template[40:] + template[:40]
        assert sites(rotated, primer, True, 0) == [(n - 12 - 40, n + 12 - 40, strand, 0)], strand

    # --- a mismatch on each side of the origin ---
    for strand in [1, -1]:
        template = dna(200)
        n = len(template)
        across = template[n - 12:] + template[:12]               # template bases n-12..n-1, then 0..11
        # top-strand positions in `across` (0..11 left of the origin, 12..23 right of it), all inside the 15-base seed of both strands
        for left, right in ((9, 12), (11, 12), (10, 14), (11, 13)):
            changed = mutate(across, [left, right])
            primer = changed if strand > 0 else rc(changed)
            assert sites(template, primer, True, 2, SEED) == [(n - 12, n + 12, strand, 2)], (left, right, strand)
            assert sites(template, primer, True, 1, SEED) == [], f"two mismatches in the seed, one allowed (strand {strand})"
            # one side at a time
            for k in (left, right):
                one = mutate(across, [k])
                assert sites(template, one if strand > 0 else rc(one), True, 1) == [(n - 12, n + 12, strand, 1)], (k, strand)


def test_template_ends_linear_versus_circular_and_primers_longer_than_the_template():
    """A linear template does not find a site across its end and a circular one does; a primer longer than the template and very short primers do not panic."""
    # --- linear versus circular ends ---
    template = dna(200)
    n = len(template)
    primer = template[n - 10:] + template[:14]
    assert sites(template, primer, False, 0) == []
    assert sites(template, primer, True, 0) == [(n - 10, n + 14, 1, 0)]
    assert sites(template, rc(primer), False, 0) == []
    assert sites(template, rc(primer), True, 0) == [(n - 10, n + 14, -1, 0)]
    # a tail that runs off the end of a linear template is just a tail ...
    tail = "".join(other(c) for c in template[n - 5:])       # does not match what is at the other end of a circle
    primer = tail + template[:24]
    assert sites(template, primer, False, 0) == [(0, 24, 1, 0)]
    assert sites(template, primer, True, 0) == [(0, 24, 1, 0)]
    # ... but when it does match, on a circle it is part of the site (through the origin)
    primer = template[n - 5:] + template[:24]
    assert sites(template, primer, False, 0) == [(0, 24, 1, 0)]
    assert sites(template, primer, True, 0) == [(n - 5, n + 24, 1, 0)]
    # the same at the right end, for a reverse primer
    primer = rc(template[n - 24:] + template[:5])
    assert sites(template, primer, False, 0) == [(n - 24, n, -1, 0)]
    assert sites(template, primer, True, 0) == [(n - 24, n + 5, -1, 0)]

    # --- primers longer than the template, very short primers ---
    for n in (1, 2, 5, 8, 12, 20):
        template = dna(n)
        for length in range(0, 41):
            primer = (template * 41)[:length]
            for circular in (False, True):
                for seed in (4, 8, 15, 30):
                    got = sites(template, primer, circular, 1, seed)
                    assert got == reference_sites(template, circular, primer, 1, seed)
                    for s, e, strand, mm in got:
                        assert 0 <= s < n and 0 < e - s <= n and strand in (1, -1) and mm >= 0
                    if length < 8:
                        assert got == []
    # a primer longer than a circular template winds round it: the site is the whole template at most
    template = "ACGGTTCAGGCT"
    got = sites(template, template + "ACGG", True, 0)
    assert got and all(e - s <= 12 for s, e, *_ in got)
    # a primer between 8 and the seed: the whole primer is the seed, one mismatch allowed except in its last 3 bases
    template = dna(100)
    for length in range(8, SEED + 1):
        core = template[40:40 + length]
        assert sites(template, core, False, 0) == [(40, 40 + length, 1, 0)]
        assert sites(template, mutate(core, [0]), False, 0) == []
        assert sites(template, mutate(core, [0]), False, 1) == [(40, 40 + length, 1, 1)]
        assert sites_at(template, mutate(core, [length - 1]), 40, 40 + length, False, 3) == [], "the last base never differs"


# ---- case, U, N and IUPAC ------------------------------------------------------------------------------------------------------------

def test_lower_case_u_and_other_letters_in_primer_or_template():
    flank = "T" * 20
    core = "ACGTTGCAGGCTAGCTAGGA"
    template = flank + core + flank
    expected = [(20, 40, 1, 0)]
    assert sites(template, core, False, 0) == expected
    assert sites(template, core.lower(), False, 0) == expected, "case says nothing"
    assert sites(template.lower(), core, False, 0) == expected
    assert sites(template, core.replace("T", "U"), False, 0) == expected, "U reads as T"
    doc = Document.new_sequence("t", template.lower(), False)        # a document keeps the case it is given
    assert doc.binding_sites(core, 0) == expected
    # a letter that is not a base or an IUPAC code is a base that differs, allowed only where a mismatch is
    assert sites(template, core[:15] + "X" + core[16:], False, 0) == []
    assert sites(template, core[:15] + "X" + core[16:], False, 1) == [(20, 40, 1, 1)]
    # an IUPAC code in the primer matches every base it stands for (a match, not a mismatch), even at the 3' end
    assert sites(template, core[:15] + "N" + core[16:], False, 0) == expected
    assert sites(template, core[:15] + "Y" + core[16:], False, 0) == [(20, 40, 1, 0)], "Y is C or T: the template has T"
    assert sites(template, core[:15] + "R" + core[16:], False, 0) == [], "R is A or G: the template has T"
    assert sites(template, core[:15] + "R" + core[16:], False, 1) == [(20, 40, 1, 1)]
    assert sites(template, core[:-1] + "N", False, 0) == expected
    assert sites(template, core[:-1] + "Y", False, 3) == [], "Y excludes the template's A"
    # an N (or a letter that is not a base) in the template never matches
    with_n = flank + core[:10] + "N" + core[11:] + flank
    assert sites(with_n, core, False, 0) == []
    assert sites(with_n, core, False, 1) == [(20, 40, 1, 1)]


# ---- the extension -----------------------------------------------------------------------------------------------------------------


def test_the_extension_stops_three_below_the_best_and_follows_the_ends_and_the_mirror():
    """The extension stops three below the best (m = matches the template, x = does not, from the seed toward the 5' end; both strands); it does not run past the end of a linear template but wraps on a circle; a site and its mirror image are the same site on the other strand; a lone mismatch gives the mirrored result."""
    # --- the extension stops three below the best ---
    for pattern, extension, mismatches in [
        ("mmm", 3, 0),
        ("mmmxmmm", 7, 1),           # a lone mismatch followed by 3 good bases is crossed
        ("mmxmm", 2, 0),             # only 2 follow it: not worth it, the rest is a tail
        ("mmxxmmmm", 2, 0),          # two in a row end it
        ("xmmm", 4, 1),              # even right next to the seed
        ("xmm", 0, 0),
        ("mxmxmmmmm", 9, 2),         # 3 below the best is not "more than 3": it goes on and wins later
        ("mxmxmxm", 1, 0),           # 4 below the best stops it
        ("mmmmmxxmmmmmmmm", 5, 0),   # the 2nd mismatch ends it: the matches beyond are a tail
        ("", 0, 0),
    ]:
        for strand in [1, -1]:
            # pattern = from the seed toward the primer's 5' end: m = matches the template, x = does not.
            seed_bases = dna(15)
            outward = dna(len(pattern))                              # the template next to the seed, going outward from it
            pairs = [other(b) if c == "x" else b for b, c in zip(outward, pattern)]   # the primer's bases there
            if strand > 0:
                template = dna(60) + outward[::-1] + seed_bases + dna(60)
                at = 60 + len(pattern)                               # where the seed starts
                primer = "".join(reversed(pairs)) + seed_bases
                expected = (at - extension, at + 15, 1, mismatches)
            else:
                template = dna(60) + seed_bases + outward + dna(60)
                at = 60
                primer = rc(seed_bases + "".join(pairs))
                expected = (at, at + 15 + extension, -1, mismatches)
            assert sites(template, primer, False, 0) == [expected], (pattern, strand)

    # --- the extension at the ends of a template ---
    for strand in [1, -1]:
        template = dna(120)
        n = len(template)
        core = template[:20]
        # the primer reaches 10 bases beyond the start of the template (forward) or its end (reverse)
        if strand > 0:
            overhang = "".join(other(c) for c in template[n - 10:])
            primer_linear = overhang + core
            primer_circle = template[n - 10:] + core
            assert sites(template, primer_linear, False, 0) == [(0, 20, 1, 0)], strand
            assert sites(template, primer_circle, False, 0) == [(0, 20, 1, 0)], "linear: nothing before the first base"
            assert sites(template, primer_circle, True, 0) == [(n - 10, n + 20, 1, 0)], "circular: through the origin"
        else:
            core = template[n - 20:]
            primer_circle = rc(core + template[:10])
            assert sites(template, primer_circle, False, 0) == [(n - 20, n, -1, 0)], strand
            assert sites(template, primer_circle, True, 0) == [(n - 20, n + 10, -1, 0)], strand
        # no site ever starts before 0 or ends after the template on a linear one
        for _ in range(200):
            start = rng.randrange(n - 10)
            p = template[start:start + rng.randint(8, 40)]
            p = p if rng.random() < 0.5 else rc(p)
            for s, e, *_ in sites(template, p, False, 1, 8):
                assert 0 <= s < e <= n, strand

    # --- a site and its mirror image ---
    for circular in [False, True]:
        # Reverse complementing the template turns every site into the mirrored site on the other strand (coordinates, count,
        # extension all mirrored): the reverse strand follows the same rules as the forward one. Includes lone mismatches.
        r = random.Random(23)
        checked = 0
        for _ in range(600):
            template = "".join(r.choice("ACGT") for _ in range(r.randint(60, 120)))
            n = len(template)
            start = r.randrange(n - 30)
            primer = template[start:start + r.randint(14, 36)]
            primer = "".join(r.choice("ACGT") if r.random() < 0.08 else c for c in primer)
            limit = r.randint(0, 2)
            here = sites(template, primer, circular, limit, SEED)
            mirrored = sites(rc(template), primer, circular, limit, SEED)
            want = sorted(((n - e) % n if circular else n - e, strand * -1, m, e - s) for s, e, strand, m in here)
            got = sorted((s, strand, m, e - s) for s, e, strand, m in mirrored)
            assert got == want, (template, primer, circular, limit)
            checked += len(here)
        assert checked > 300, circular

    # --- a lone mismatch mirrored on the reverse strand ---
    template = dna(300)
    core = template[100:124]
    # 4 bases from the 5' end: crossed (4 matching bases follow it); 2: too close, the end becomes a tail
    assert sites(template, primer_for(core, 1, [4]), False, 0) == [(100, 124, 1, 1)]
    assert sites(template, primer_for(core, -1, [4]), False, 0) == [(100, 124, -1, 1)]
    assert sites(template, primer_for(core, 1, [2]), False, 0) == [(103, 124, 1, 0)]
    assert sites(template, primer_for(core, -1, [2]), False, 0) == [(100, 121, -1, 0)]


# ---- mismatch positions (the Editor draws them) -------------------------------------------------------------------------------------------

def editor_with(template, circular, primers):
    ed = Editor(Document.new_sequence("t", template, circular), Path("t.gb"))
    ed.set_extra(primers)
    return ed


def test_mismatch_positions_are_template_coordinates():
    """Mismatch positions (what the Editor draws) are template coordinates on both strands, with two mismatches and across the origin."""
    # --- on both strands ---
    template = dna(400)
    core = template[100:130]
    for tail in ("", "GGATCC", "CCATGGAAGCTT"):
        # a tail that matches nothing next to the site (every base differs from the template there)
        fwd_tail = mutate(template[100 - len(tail):100], range(len(tail)))
        rev_tail = rc(mutate(template[130:130 + len(tail)], range(len(tail))))
        for k in (3, 5, 12, 20, 25):
            # the template base 100+k is wrong in both primers; the tail is on each primer's 5' end
            fwd = fwd_tail + mutate(core, [k])
            rev = rev_tail + rc(mutate(core, [k]))
            ed = editor_with(template, False, [("f", fwd), ("r", rev)])
            by_name = {s.label: s for s in ed.primer_sites}
            assert by_name["f"].mismatch_at == [100 + k] and by_name["f"].mismatches == 1, (tail, k)
            assert (by_name["f"].start, by_name["f"].tail) == (100, len(tail))
            assert by_name["r"].mismatch_at == [100 + k] and by_name["r"].mismatches == 1, (tail, k)
            assert (by_name["r"].end, by_name["r"].tail) == (130, len(tail))
    # an exact primer has none
    ed = editor_with(template, False, [("f", core), ("r", rc(core))])
    assert [s.mismatch_at for s in ed.primer_sites] == [[], []]

    # --- two mismatches and across the origin ---
    template = dna(300)
    n = len(template)
    across = template[n - 15:] + template[:15]
    changed = mutate(across, [4, 20])                  # one on each side of the origin: template n-11 and 5
    for strand in (1, -1):
        (site,) = editor_with(template, True, [("p", changed if strand > 0 else rc(changed))]).primer_sites
        assert site.strand == strand and site.spans == [(n - 15, n), (0, 15)]
        assert site.mismatch_at == [n - 11, 5] and site.mismatches == 2
        (site,) = editor_with(template, True, [("p", across if strand > 0 else rc(across))]).primer_sites
        assert site.mismatch_at == [] and site.mismatches == 0


# ---- primer pairs and their product ------------------------------------------------------------------------------------------------------

def pair_editor(template, circular, f_site, r_site, f_tail="", r_tail=""):
    fwd = f_tail + template[f_site[0]:f_site[1]]
    rev = r_tail + rc(template[r_site[0]:r_site[1]])
    return editor_with(template, circular, [("F", fwd), ("R", rev)])


def site_named(ed, label):
    return next(s for s in ed.primer_sites if s.label == label)


def test_primer_pair_products_linear_circular_and_across_the_origin():
    """A primer pair's product: linear with tails, a reverse primer upstream (none when linear, through the origin on a circle), across the origin, a mismatched primer, and the pair with the shortest product is picked."""
    # --- linear with tails ---
    template = dna(500)
    ed = pair_editor(template, False, (100, 124), (300, 324), "GGAATTCC", "CCAAGCTTGG")
    f, r = site_named(ed, "F"), site_named(ed, "R")
    assert (f.strand, r.strand) == (1, -1) and (f.start, f.end, r.start, r.end) == (100, 124, 300, 324)
    assert (f.tail, r.tail) == (8, 10)
    assert tools.product_range(ed, f, r) == (100, 324)
    doc = tools.pcr_product(ed, f, r)
    assert doc.sequence() == "GGAATTCC" + template[100:324] + rc("CCAAGCTTGG")
    assert doc.length == 224 + 8 + 10 and not doc.circular
    # the dialog (built, not run) says the same size and makes the same product
    dialog = tools.PcrDialog(None, ed)
    assert f"{224 + 8 + 10} bp" in dialog.info.text() and dialog.pair()[0].label == "F"
    assert dialog.product().sequence() == doc.sequence()

    # --- reverse primer upstream of the forward one ---
    template = dna(500)
    # the forward primer at 350, the reverse one at 50..74
    ed = pair_editor(template, False, (350, 374), (50, 74))
    f, r = site_named(ed, "F"), site_named(ed, "R")
    assert tools.product_range(ed, f, r) is None
    assert "not downstream" in tools.PcrDialog(None, ed).info.text()
    ed = pair_editor(template, True, (350, 374), (50, 74), "GG", "CC")
    f, r = site_named(ed, "F"), site_named(ed, "R")
    assert tools.product_range(ed, f, r) == (350, 574)
    doc = tools.pcr_product(ed, f, r)
    assert doc.sequence() == "GG" + template[350:] + template[:74] + "GG"
    assert doc.length == 150 + 74 + 4
    assert "228 bp" in tools.PcrDialog(None, ed).info.text()
    # two forward primers, or two reverse primers, make no pair: the dialog needs one of each
    ed = editor_with(template, False, [("F1", template[50:74]), ("F2", template[100:124])])
    assert tools.PcrDialog(None, ed).pair() is None

    # --- across the origin ---
    template = dna(400)
    n = len(template)
    # the forward primer itself crosses the origin; the reverse one is further on
    fwd = template[n - 12:] + template[:12]
    tail = rc(mutate(template[124:130], range(6)))            # a 5' tail that matches nothing next to the site
    rev = tail + rc(template[100:124])
    ed = editor_with(template, True, [("F", fwd), ("R", rev)])
    f, r = site_named(ed, "F"), site_named(ed, "R")
    assert (f.start, f.end, f.strand) == (n - 12, n + 12, 1)
    assert (r.start, r.end, r.strand) == (100, 124, -1)
    assert tools.product_range(ed, f, r) == (n - 12, n + 124)
    doc = tools.pcr_product(ed, f, r)
    assert doc.sequence() == template[n - 12:] + template[:124] + rc(tail)
    assert doc.length == 12 + 124 + 6
    # a mismatched primer changes the product the way a real PCR does
    mutated = mutate(template[100:124], [8])
    ed = editor_with(template, False, [("F", template[50:74]), ("R", rc(mutated))])
    f, r = site_named(ed, "F"), site_named(ed, "R")
    assert r.mismatches == 1 and (f.strand, r.strand) == (1, -1)
    assert tools.pcr_product(ed, f, r).sequence() == template[50:100] + mutated

    # --- the shortest product is picked ---
    template = dna(600)
    ed = editor_with(template, False, [("F1", template[50:74]), ("F2", template[200:224]),
                                       ("R1", rc(template[300:324])), ("R2", rc(template[500:524]))])
    dialog = tools.PcrDialog(None, ed)
    f, r = dialog.pair()
    assert (f.label, r.label) == ("F2", "R1")
    assert f"{324 - 200} bp" in dialog.info.text()


# ---- a genome in a FASTA file ---------------------------------------------------------------------------------------------------------

def write_fasta(path, records, width=60, crlf=False):
    nl = "\r\n" if crlf else "\n"
    with open(path, "w", newline="") as fh:
        for name, seq in records:
            fh.write(f">{name} test{nl}")
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + nl)


def test_genome_search_in_a_fasta_file_finds_the_same_sites_as_the_search_in_memory(tmp_path):
    g = random.Random(2)
    chr1 = "".join(g.choice("ACGT") for _ in range(900))
    chr2 = "".join(g.choice("ACGT") for _ in range(1400))
    primer = "TTAC" + chr2[300:324]                       # a tail, then 24 bases that bind chr2 at 300
    chr2 = chr2[:700] + rc(chr2[300:324]) + chr2[724:]   # and a copy on the other strand
    chr1 = chr1[:100] + chr2[300:312] + chr1[112:]       # a partial copy of the 3' end of the site in chr1 is not a binding site
    path = tmp_path / "genome.fa"
    write_fasta(path, [("chr1", chr1), ("chr2", chr2)])
    fa = Fasta(str(path))
    assert [name for name, _ in fa.sequences()] == ["chr1", "chr2"]
    for name, seq in (("chr1", chr1), ("chr2", chr2)):
        doc = Document.from_fasta(str(path), name)
        assert doc.sequence() == seq
        for limit in (0, 1, 2):
            for p in (primer, rc(primer), mutate(primer, [16]), primer[4:]):
                assert doc.binding_sites(p, limit) == sites(seq, p, False, limit) == reference_sites(seq, False, p, limit, SEED)
        # the exact sites of the whole annealing part, from the file's own motif search (reads the file by line, both strands)
        exact = [(s, e, strand, 0) for s, e, strand, mm in fa.search_motif(name, primer[4:], 0, False, True)]
        assert sorted(exact) == sorted(h for h in doc.binding_sites(primer[4:], 0) if h[1] - h[0] == len(primer) - 4)
    assert len(Document.from_fasta(str(path), "chr2").binding_sites(primer, 0)) == 2
    assert Document.from_fasta(str(path), "chr1").binding_sites(primer, 0) == []
    # lines of any width and CRLF endings give the same sites
    path2 = tmp_path / "genome_crlf.fa"
    write_fasta(path2, [("chr1", chr1), ("chr2", chr2)], width=70, crlf=True)
    assert Document.from_fasta(str(path2), "chr2").binding_sites(primer, 1) == Document.from_fasta(str(path), "chr2").binding_sites(primer, 1)


# ---- non-specific priming: yellow when the 3' end matches elsewhere, red when the whole primer does -----------------------------------

def checks_by_name(ed, primer):
    return {name: (level, text) for name, level, text in primer_panel.primer_checks(ed, primer)}


def test_non_specific_priming_in_the_sequence_and_in_the_genome(monkeypatch):
    """Non-specific priming: in the sequence yellow for a 3' end matching elsewhere and red for a full match; in the genome check."""
    # --- in the sequence ---
    monkeypatch.setattr(primer_panel, "homology_settings", lambda: (12, 0))
    template = dna(700)
    primer = template[200:224]
    # unique: green
    level, _ = checks_by_name(editor_with(template, False, []), primer)["Non-specific (sequence)"]
    assert level == "ok"
    # the last 12 bases of the primer elsewhere, after 3 different bases, so that it is not a binding site itself
    partial = template[:400] + "".join(other(c) for c in primer[-15:-12]) + primer[-12:] + template[415:]
    ed = editor_with(partial, False, [])
    checks = primer_panel.primer_checks(ed, primer)
    by_name = {name: (level, text) for name, level, text in checks}
    level, text = by_name["Non-specific (sequence)"]
    assert level == "warn" and "404–415 →" in text, "a warning, and where its 3' end matches"
    import theme
    assert primer_panel.MARKS["warn"][1] == theme.hex("status_warn"), "the theme's warning colour"
    assert by_name["Binds here"][0] == "ok", "its own site is not counted as mis-priming"
    assert primer_panel.MARKS["warn"][1] in primer_panel.checks_html(checks)
    # the same from a reverse primer: its 3' end is the start of `primer`, matched elsewhere on the top strand
    partial_rev = template[:400] + primer[:12] + "".join(other(c) for c in primer[12:15]) + template[415:]
    assert checks_by_name(editor_with(partial_rev, False, []), rc(primer))["Non-specific (sequence)"][0] == "warn"
    # the whole primer a second time: red
    ed = editor_with(template[:400] + primer + template[424:], False, [])
    checks = primer_panel.primer_checks(ed, primer)
    level, text = {name: (level, text) for name, level, text in checks}["Non-specific (sequence)"]
    assert level == "bad" and "2 places" in text
    assert primer_panel.MARKS["bad"][1] == "#c62828" and primer_panel.MARKS["bad"][1] in primer_panel.checks_html(checks)
    assert primer_panel.worst(checks) == "bad"
    # a lower-case 5' tail is not part of what must match elsewhere; a reverse primer is judged the same way
    ed = editor_with(template, False, [])
    assert checks_by_name(ed, "gaattcgg" + primer)["Non-specific (sequence)"][0] == "ok"
    assert checks_by_name(ed, rc(primer))["Non-specific (sequence)"][0] == "ok"
    assert checks_by_name(editor_with(template[:400] + primer + template[424:], False, []), rc(primer))["Non-specific (sequence)"][0] == "bad"

    # --- in the genome ---
    monkeypatch.setattr(primer_panel, "homology_settings", lambda: (12, 0))
    template = dna(700)
    primer = template[200:224]
    ed = editor_with(template, False, [])
    assert checks_by_name(ed, primer)["Non-specific (genome)"][0] == "info", "no genome set"
    other_genome = dna(5000)
    for genome_seq, level in ((other_genome, "ok"),
                              (other_genome[:1000] + "".join(other(c) for c in primer[-15:-12]) + primer[-12:] + other_genome[1015:], "warn"),
                              (other_genome[:1000] + rc(primer) + other_genome[1024:], "bad")):
        genome = Document.new_sequence("genome", genome_seq, True)
        monkeypatch.setattr(primer_panel, "genome", lambda genome=genome: genome)
        got, text = checks_by_name(ed, primer)["Non-specific (genome)"]
        assert got == level, text
    # a genome that is the sequence itself is not compared with itself
    same = Document.new_sequence("t", template, False)
    monkeypatch.setattr(primer_panel, "genome", lambda: same)
    assert checks_by_name(ed, primer)["Non-specific (genome)"][0] == "info"


# ---- real files --------------------------------------------------------------------------------------------------------------------------

def test_our_sites_are_the_binding_sites_stored_in_each_real_dna_file():
    from sample_files import DNA_FILE as path   # a copy (see sample_files.py)
    name = path.name
    if not path.exists():
        pytest.skip(f"missing {name}")
    doc = Document(str(path))
    if not doc.primers():
        pytest.skip(f"{name} has no primers, so there is nothing stored in the file to compare with")
    n = doc.length
    stored = {(f.label, tuple(f.spans), f.strand) for f in doc.features() if f.kind == "primer_bind"}
    def bases(spans):
        return {p % n for a, b in spans for p in range(a, b)}

    compared = 0
    for primer_name, sequence in doc.primers():
        for start, end, strand, _ in doc.binding_sites(sequence, 0):
            spans = ((start, end),) if end <= n else ((start, n), (0, end - n))
            # The same site, or (a primer with a 5' tail the file also matched) inside the stored one, ending at the same 3' base.
            three = (end - 1) % n if strand > 0 else start % n
            assert any(name == primer_name and s == strand and bases(spans) <= bases(sg) and three in bases(sg)
                       for name, sg, s in stored), primer_name
            compared += 1
    assert compared >= len(stored) >= 15, (compared, len(stored))   # the comparison is never empty: one site per stored primer at least


# ---- saving primer sites ------------------------------------------------------------------------------------------------------------------

def test_primer_sites_written_on_save_do_not_change_the_document_or_its_undo_stack(tmp_path):
    template = dna(600)
    path = tmp_path / "t.gb"
    ed = Editor(Document.new_sequence("t", template, True), path)
    ed.add_features([("misc_feature", "11..40", [("label", "alpha")])])
    ed.select(300, 300)
    ed.type_bases("ACGTACGTAC")
    ed.undo()
    ed.redo()
    template = ed.doc.sequence()
    assert add_primer(ed, "pFwd", template[100:124]) == "added"
    assert add_primer(ed, "pRev", rc(template[400:424])) == "added"
    assert len(ed.primer_sites) == 2

    def snapshot():
        d = ed.doc
        feats = [(f.kind, f.label, f.location, tuple(f.spans), f.strand, tuple(f.qualifiers)) for f in d.features()]
        return d.sequence(), feats, d.history(), d.circular, d.length, d.modified

    before = snapshot()
    ed.save(path)
    assert snapshot()[:5] == before[:5], "saving changes nothing in the document, its features or its history"
    assert ed.doc.modified is False
    assert not [f for f in ed.doc.features() if f.kind == "primer_bind"], "the sites are written, never added to the document"
    # the file has the sites, found on the current sequence
    saved = Document(str(path))
    written = sorted((f.label, tuple(f.spans), f.strand) for f in saved.features() if f.kind == "primer_bind")
    assert written == [("pFwd", ((100, 124),), 1), ("pRev", ((400, 424),), -1)]
    assert saved.sequence() == template
    # saving again writes the same bytes and still changes nothing
    first = path.read_bytes()
    ed.save(path)
    assert path.read_bytes() == first and snapshot()[:5] == before[:5]
    # Undo goes back through the primers and the edits as before, the save is not a step
    steps = len(ed.doc.history()[0])
    assert steps == 4        # feature, insert, two primers
    while ed.doc.history()[0]:
        assert ed.undo()
    assert ed.doc.sequence() == template[:300] + template[310:] and not ed.doc.features()
    assert primer_db.db().rows == [] and ed.primer_sites == []
    # the sites follow an edit made after the save
    ed2 = Editor(Document.new_sequence("t", template, True), tmp_path / "t2.gb")
    add_primer(ed2, "pFwd", template[100:124])
    ed2.select(0, 0)
    ed2.type_bases("AAAAAAAAAA")
    ed2.save(tmp_path / "t2.gb")
    moved = Document(str(tmp_path / "t2.gb"))
    assert [(f.label, tuple(f.spans)) for f in moved.features() if f.kind == "primer_bind"] == [("pFwd", ((110, 134),))]


# ---- degenerate (IUPAC) primers ------------------------------------------------------------------------------------------------------

def test_a_primer_with_three_degenerate_positions_binds_every_site_it_matches_on_both_strands():
    core = dna(30)
    # the same 30 bases with the three degenerate positions (4, 17, 25) taking each allowed base in turn
    def variant(r, y, n):
        v = list(core)
        v[4], v[17], v[25] = r, y, n
        return "".join(v)
    primer = variant("R", "Y", "N")
    template = "".join(dna(60, "T") + v for v in (variant("A", "C", "G"), variant("G", "T", "T"), variant("A", "T", "A"))) + dna(60, "T")
    forward = [(60, 90, 1, 0), (150, 180, 1, 0), (240, 270, 1, 0)]
    assert sites(template, primer, False, 0) == forward
    # on the reverse-complemented template the same primer binds the other strand (R pairs with Y there), still no mismatch
    found = sites(rc(template), primer, False, 0)
    assert len(found) == 3 and all(strand == -1 and mm == 0 for _, _, strand, mm in found)
    # a base the code does not allow, inside the 3' seed, is still a mismatch: Y (C or T) against an A
    off = template[:60] + variant("A", "A", "G") + template[90:]
    assert [x[0] for x in sites(off, primer, False, 0)] == [150, 240]


def test_a_degenerate_base_in_the_template_still_matches_only_the_same_letter():
    core = dna(30)
    template = "T" * 40 + core[:20] + "N" + core[21:] + "T" * 40
    assert sites(template, core, False, 0) == []
    assert sites(template, core[:20] + "N" + core[21:], False, 0) == [(40, 70, 1, 0)], "N against N is the same letter"


def test_the_tm_range_brackets_the_tm_of_every_expanded_oligo():
    from itertools import product
    from opendna_core import tm, tm_range
    degenerate = "AGCGGRTAACAATYTCACACNGGA"
    codes = dict(R="AG", Y="CT", N="ACGT")
    oligos = ["".join(p) for p in product(*[codes.get(c, c) for c in degenerate])]
    tms = [tm(o) for o in oligos]
    lo, hi = tm_range(degenerate)
    assert len(oligos) == 16
    assert lo == pytest.approx(min(tms), abs=1e-6) and hi == pytest.approx(max(tms), abs=1e-6)
    assert lo < hi
    assert tm(degenerate) is None, "tm itself stays for one plain oligo"
    plain = "AGCGGATAACAATTTCACACAGGA"
    assert tm_range(plain) == pytest.approx((tm(plain), tm(plain)))
    assert tm_range("ACGTX") is None


def test_the_primer_panel_shows_a_tm_range_for_a_degenerate_primer_and_a_plain_tm_for_a_plain_one():
    plain = "AGCGGATAACAATTTCACACAGGA"
    assert primer_panel.tm_text(plain) == f"{primer_panel.melting_temp(plain):.1f}"
    shown = primer_panel.tm_text("AGCGGRTAACAATYTCACACNGGA")
    low, high = (float(x) for x in shown.split("–"))
    assert low < high
    assert primer_panel.tm_text("ACGTX") == ""
    assert "not calculated" in primer_panel.structure_text("AGCGGRTAACAATYTCACACNGGA")
    assert primer_panel.tm_sort_key("AGCGGRTAACAATYTCACACNGGA") == pytest.approx(low, abs=0.06)


def test_a_degenerate_primer_can_be_saved_to_the_database():
    assert primer_db.primer_problem("fam", "AGCGGRTAACAATYTCACACNGGA") == ""
