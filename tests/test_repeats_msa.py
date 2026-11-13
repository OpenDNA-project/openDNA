"""Repeats and multiple sequence alignment: the repeat finder on planted repeats (and nothing in random DNA), the shared
stretches of two sequences, and the built-in multiple aligner with its consensus, identity, column mapping and aligned-FASTA
reader."""
import os
import random
import sys
import tempfile
from pathlib import Path

os.environ["OPENDNA_DATA_DIR"] = tempfile.mkdtemp(prefix="opendna_test_")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "app"))

import pytest

import msa
import repeats

COMP = str.maketrans("ACGTacgt", "TGCAtgca")
rng = random.Random(11)
random.seed(3)   # (the module-level generator, as before: the planted 3-base deletion must not be able to slide)
BASE = "".join(random.choice("ACGT") for _ in range(200))


def rc(s):
    return s.translate(COMP)[::-1]


def dna(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def kinds(result, kind):
    return [r for r in result.repeats if r.kind == kind]


def test_repeats_are_found_where_planted_and_not_in_random_dna():
    """Planted direct (3 copies), inverted and tandem repeats come back at their places; random DNA has none at the default
    thresholds; thresholds and mismatches decide what is listed; a repeat over the origin of a circle; a long input is
    capped and says so; the stretches two sequences share, in both orientations, and a sequence against itself."""
    d, h, unit = dna(30), dna(16), "GATTACAGCT"
    seq = dna(100) + d + dna(150) + d + dna(80) + d + dna(60) + h + dna(90) + rc(h) + dna(70) + unit * 5 + dna(100)
    found = repeats.find_repeats(seq)
    (direct,) = kinds(found, "direct")
    # (a random neighbouring base may match by chance and lengthen a repeat by a base: positions are checked within that)
    near = lambda got, want: len(got) == len(want) and all(abs(g - w) <= 2 for g, w in zip(got, want))
    assert 30 <= direct.length <= 32 and near([c[0] for c in direct.copies], [seq.find(d), seq.find(d, 140), seq.rfind(d)])
    assert all(d in seq[s - 2:e + 2] for s, e, _ in direct.copies)
    (inv,) = kinds(found, "inverted")
    assert 16 <= inv.length <= 18 and near([c[0] for c in inv.copies], [seq.find(h), seq.find(rc(h))])
    assert [c[2] for c in inv.copies] == [1, -1]
    (tandem,) = kinds(found, "tandem")
    assert tandem.unit == 10 and len(tandem.copies) == 5 and abs(tandem.copies[0][0] - seq.find(unit)) <= 2
    assert not found.capped
    # --- random DNA ---
    assert not repeats.find_repeats(dna(20_000)).repeats
    assert not repeats.find_repeats(dna(5_000), circular=True).repeats
    assert len(repeats.find_repeats(dna(20_000), max_loop=None).repeats) > 0, "a 12-base stem anywhere in 20 kb is chance, hence the loop limit"
    # --- thresholds and mismatches ---
    d = dna(14)
    seq = dna(80) + d + dna(80) + d + dna(80)
    assert not kinds(repeats.find_repeats(seq), "direct")
    assert len(kinds(repeats.find_repeats(seq, min_direct=12), "direct")) == 1
    long_one = dna(30)
    bent = long_one[:15] + ("A" if long_one[15] != "A" else "C") + long_one[16:]
    seq = dna(80) + long_one + dna(80) + bent + dna(80)
    assert not [r for r in kinds(repeats.find_repeats(seq), "direct") if r.length >= 30]
    (r,) = kinds(repeats.find_repeats(seq, mismatches=1), "direct")
    assert r.length >= 30 and r.mismatches >= 1
    # --- over the origin of a circle ---
    d = dna(25)
    near = lambda got, want: len(got) == len(want) and all(abs(g - w) <= 2 for g, w in zip(got, want))
    seq = d[10:] + dna(100) + d + dna(100) + d[:10]
    (r,) = kinds(repeats.find_repeats(seq, circular=True), "direct")
    assert 25 <= r.length <= 28 and near(sorted(c[0] for c in r.copies), [115, 240])
    # --- a long input is capped ---
    result = repeats.find_repeats("ACGT" * 200_000)   # nothing but a microsatellite: every k-mer is skipped
    assert result.capped and "skipped" in result.capped
    # --- two sequences ---
    d, e = dna(40), dna(30)
    a = dna(50) + d + dna(60) + e + dna(50)
    b = dna(30) + rc(e) + dna(40) + d + dna(20)
    found = {r[0]: r for r in repeats.pair_repeats(a, b, 20).repeats}
    near = lambda got, want: all(abs(g - w) <= 2 for g, w in zip(got, want))   # (a random flank base may extend a stretch)
    assert near(found["direct"][1:], (50, 100, 40)) and near(found["inverted"][1:], (150, 30, 30))
    assert not repeats.pair_repeats(dna(3000), dna(3000), 20).repeats
    seq = dna(100) + d + dna(100) + d + dna(100)
    (r,) = repeats.pair_repeats(seq, seq, 30).repeats   # itself: the pair once, not the main diagonal
    assert r[0] == "direct" and near(r[1:4], (100, 240, 40))
    direct, inverted, _ = repeats.pair_dots(a, b)
    assert direct and inverted


def planted():
    """Four copies of BASE: one identical, one with the base at 50 changed, one with the base at 120 changed, one missing 3 bases at 80."""
    swap = lambda s, i: s[:i] + ("A" if s[i] != "A" else "C") + s[i + 1:]
    return ["ref", "sub50", "sub120", "del80"], [BASE, swap(BASE, 50), swap(BASE, 120), BASE[:80] + BASE[83:]]


def test_multiple_alignment_consensus_identity_and_parsing():
    """The built-in aligner on planted substitutions and one indel (rows keep their sequences, the deletion is a 3-column gap,
    consensus and identity); columns and positions map both ways through gaps; aligned FASTA is read by name; one sequence
    or too long a sequence is refused."""
    names, seqs = planted()
    a = msa.align(names, seqs)
    assert [r.replace("-", "") for r in a.rows] == seqs
    assert a.width == 200 and a.rows[3][80:83] == "---" and a.rows[0] == BASE
    assert a.consensus() == BASE            # every variant is outvoted 3 to 1 or 2 to 2 (gaps count only when most rows have one)
    assert a.identity(0) == 100.0
    assert a.identity(1) == pytest.approx(100 * 199 / 200)
    assert a.identity(3) == pytest.approx(100 * 197 / 200)   # 3 gap columns where the reference has bases
    a.reference = None                                         # against the consensus, which here is the reference
    assert a.identities() == [a.identity(i, 0) for i in range(4)]
    cons = a.conservation()
    assert cons[0] == 1.0 and cons[50] == 0.75 and cons[81] == 0.75   # the substituted base, and a base the deletion lacks
    # --- columns and positions through gaps ---
    a = msa.Alignment(["x", "y"], ["AC--GT", "ACTTGT"])
    assert [a.position(0, c) for c in range(6)] == [0, 1, 2, 2, 2, 3]
    assert a.column_of(0, 2) == 4 and a.column_of(1, 2) == 2
    # --- aligned FASTA ---
    names, rows = msa.parse_aligned_fasta(">a\nAC-\nGT\n>b desc\nACTG\nT-\n")
    assert names == ["a", "b desc"] and rows == ["AC-GT", "ACTGT-"]
    # --- refusals ---
    with pytest.raises(ValueError):
        msa.align(["a"], ["ACGT"])
    with pytest.raises(ValueError):
        msa.align(["a", "b"], ["A" * (msa.MAX_BUILTIN + 1), "A"])
