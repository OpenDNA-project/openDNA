"""The sequence algorithms of the Rust core (`opendna_core`): melting temperature, hairpins and dimers, GC bins, ORFs,
restriction cuts and digests, find_parts.

Every test has an independent reference: Biopython (Tm, Restriction, translation), a hand count, a small model written here
(ORFs), or a value the docs / code comments promise. Nothing here needs Qt.
"""
import math
import random
import sys
from itertools import product
from pathlib import Path

import pytest
from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "app"))

import opendna_core as core
from opendna_core import Document, tm, tm_range

import enzymes
import tools
from enzymes import Enzyme

from sample_files import PLASMID   # a copy (see sample_files.py)

COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def rc(s):
    return s.translate(COMP)[::-1]


def rand_dna(rng, n, alphabet="ACGT"):
    return "".join(rng.choice(alphabet) for _ in range(n))


# ---- Tm ------------------------------------------------------------------------------------------------------------------

# SantaLucia 1998 "unified" nearest-neighbour parameters (Table 1 of the paper): dH kcal/mol, dS cal/(K mol), copied here from
# the paper so the hand calculation below shares no code with the Rust one.
SL98 = {"AA": (-7.9, -22.2), "TT": (-7.9, -22.2), "AT": (-7.2, -20.4), "TA": (-7.2, -21.3), "CA": (-8.5, -22.7), "TG": (-8.5, -22.7),
        "GT": (-8.4, -22.4), "AC": (-8.4, -22.4), "CT": (-7.8, -21.0), "AG": (-7.8, -21.0), "GA": (-8.2, -22.2), "TC": (-8.2, -22.2),
        "CG": (-10.6, -27.2), "GC": (-9.8, -24.4), "GG": (-8.0, -19.9), "CC": (-8.0, -19.9)}
SL98_INIT = {"GC": (0.1, -2.8), "AT": (2.3, 4.1)}   # initiation, per duplex end


def tm_by_hand_1m_nacl(seq, primer_nm):
    """Tm (C) at 1 M NaCl straight from the SantaLucia 1998 formula: dH / (dS + R ln(CT / 4)), CT / 4 for a non-self-complementary
    strand (no salt correction at 1 M: that is the reference condition of the parameters)."""
    dh = ds = 0.0
    for end in (seq[0], seq[-1]):
        h, s = SL98_INIT["AT" if end in "AT" else "GC"]
        dh, ds = dh + h, ds + s
    for pair in zip(seq, seq[1:]):
        h, s = SL98["".join(pair)]
        dh, ds = dh + h, ds + s
    return dh * 1000 / (ds + 1.987 * math.log(primer_nm * 1e-9 / 4)) - 273.15


# Three well-known primers (M13 forward, lac promoter primer, T7 promoter primer).
LITERATURE_PRIMERS = ["GTAAAACGACGGCCAGT", "AGCGGATAACAATTTCACACAGGA", "TAATACGACTCACTATAGGG"]


def biopython_tm(seq, na, mg, dntp, primer_nm):
    # dnac1 = dnac2 = C/2 makes Biopython use C/4, Primer3's convention for our formula.
    return mt.Tm_NN(seq, nn_table=mt.DNA_NN3, saltcorr=5, Na=na, Mg=mg, dNTPs=dntp,
                    dnac1=primer_nm / 2, dnac2=primer_nm / 2)


def test_melting_temperature_matches_santalucia_and_biopython():
    """Tm: three known primers against the SantaLucia 1998 formula and Biopython; random primers under several salt
    conditions and self-complementary primers against Biopython; the direction of length, GC, salt, Mg, dNTP and primer
    concentration; the Tm range of a degenerate primer brackets every expanded oligo."""
    # --- tm three known primers against the santalucia 1998 formula and biopython ---
    for seq in LITERATURE_PRIMERS:
        hand = tm_by_hand_1m_nacl(seq, 250)
        ours = core.tm(seq, na_mm=1000, mg_mm=0, dntp_mm=0, primer_nm=250)
        assert ours == pytest.approx(hand, abs=0.5), seq
        assert ours == pytest.approx(hand, abs=1e-6), seq   # same formula: in fact equal to rounding
        # Biopython with no salt correction and the same concentrations (dnac1 = dnac2 = C / 2 gives C / 4)
        bio = mt.Tm_NN(seq, nn_table=mt.DNA_NN3, saltcorr=0, Na=1000, dnac1=125, dnac2=125)
        assert ours == pytest.approx(bio, abs=0.5), seq
        # The same primer at the default PCR conditions is a few degrees lower and in the usual 45-70 C range.
        assert 45 < core.tm(seq) < 70
    # --- tm rises with length and gc ---
    rng = random.Random(3)
    base = rand_dna(rng, 40)
    tms = [core.tm(base[:n]) for n in range(10, 41)]
    assert all(b > a for a, b in zip(tms, tms[1:])), "each added base raises Tm (random seq, nearest-neighbour terms all negative dH)"
    # Same length, more G/C: AT-only < mixed < GC-only
    at, mixed, gc = "ATTATAATTATAATATTAAT" * 1, "ATGCATGCATGCATGCATGC", "GCGCCGGCGCCGGCGCCGGC"
    assert core.tm(at) < core.tm(mixed) < core.tm(gc)
    # Case does not matter, U reads as T
    assert core.tm("acgtacgtacgtacgt") == core.tm("ACGTACGTACGTACGT")
    assert core.tm("ACGUACGUACGUACGU") == core.tm("ACGTACGTACGTACGT")
    # --- tm salt magnesium dntp and primer concentration directions ---
    seq = "AGCGGATAACAATTTCACACAGGA"
    base = core.tm(seq, na_mm=50, mg_mm=1.5, dntp_mm=0.6, primer_nm=50)
    assert core.tm(seq, na_mm=100, mg_mm=1.5, dntp_mm=0.6, primer_nm=50) > base       # more Na+ stabilises
    assert core.tm(seq, na_mm=20, mg_mm=1.5, dntp_mm=0.6, primer_nm=50) < base
    assert core.tm(seq, na_mm=50, mg_mm=3.0, dntp_mm=0.6, primer_nm=50) > base         # more Mg2+ stabilises
    assert core.tm(seq, na_mm=50, mg_mm=1.5, dntp_mm=1.0, primer_nm=50) < base         # dNTPs bind Mg2+: less free Mg2+, lower Tm
    assert core.tm(seq, na_mm=50, mg_mm=1.5, dntp_mm=0.6, primer_nm=500) > base        # more primer: higher Tm
    # With dNTPs >= Mg2+ there is no free Mg2+ to correct for, so it is the same as none at all.
    assert core.tm(seq, na_mm=50, mg_mm=0.5, dntp_mm=0.6) == core.tm(seq, na_mm=50, mg_mm=0, dntp_mm=0)
    # --- random primers under several conditions against Biopython; non-bases give None ---
    for conditions in [(50, 1.5, 0.6, 50), (50, 0, 0, 250), (100, 2, 0.8, 500)]:
        rng = random.Random(1)
        for _ in range(200):
            seq = "".join(rng.choice("ACGT") for _ in range(rng.randint(12, 40)))
            if seq == seq[::-1].translate(str.maketrans("ACGT", "TGCA")):
                continue  # palindromes: Biopython and Primer3 use different concentrations
            assert tm(seq, *conditions) == pytest.approx(biopython_tm(seq, *conditions), abs=0.01), (conditions, seq)
    assert tm("ACGTNACGT") is None
    # --- self-complementary primers: whole concentration and the symmetry term, as Biopython's selfcomp mode ---
    palindromes = ["GAATTCGAATTC", "ACGTACGTACGT", "GGATCCGGATCC", "CGCGAATTCGCG", "AAGCTTAAGCTT", "GCGGCCGCGCGGCCGC"]
    for seq in palindromes:
        assert seq == seq[::-1].translate(str.maketrans("ACGT", "TGCA")), seq
        for na, mg, dntp, primer_nm in [(50, 1.5, 0.6, 50), (50, 0, 0, 250), (100, 2, 0.8, 500)]:
            expected = mt.Tm_NN(seq, nn_table=mt.DNA_NN3, saltcorr=5, Na=na, Mg=mg, dNTPs=dntp, dnac1=primer_nm, dnac2=0,
                                selfcomp=True)
            assert tm(seq, na, mg, dntp, primer_nm) == pytest.approx(expected, abs=0.01), (seq, na, mg, dntp, primer_nm)
            # ... and visibly different from the non-self-complementary formula, so the branch matters
            assert abs(expected - biopython_tm(seq, na, mg, dntp, primer_nm)) > 1, seq
    # --- a degenerate primer: tm_range brackets the Tm of every expanded oligo ---
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


# ---- hairpins and dimers -------------------------------------------------------------------------------------------------


def test_hairpins_and_dimers():
    """Hairpin and self-dimer of palindromes, none for poly-A or too-short primers; energy signs, symmetry and ordering
    (stem length, GC, loop size, 3' end, salt)."""
    # --- hairpin and self dimer of palindromes and none for poly a ---
    # GGATCC is its own reverse complement: stem GGATCC, loop TTTT, stem GGATCC.
    hp_any, hp_3p = core.hairpin("GGATCCTTTTGGATCC")
    assert hp_any < -1 and hp_3p == hp_any, "the stem includes the 3' end here"
    # A palindromic primer pairs with a copy of itself.
    pal = "GAATTCGAATTC"
    any_dg, three_dg = core.dimer(pal, pal)
    assert any_dg < -5 and three_dg == any_dg
    # Poly-A pairs with nothing; neither does a primer that is too short for a stem.
    assert core.hairpin("A" * 20) == (0.0, 0.0)
    assert core.dimer("A" * 20, "A" * 20) == (0.0, 0.0)
    assert core.hairpin("AC") == (0.0, 0.0) and core.hairpin("") == (0.0, 0.0)
    assert core.dimer("", "ACGT") == (0.0, 0.0)
    # --- dimer and hairpin energy signs and ordering ---
    rng = random.Random(5)
    for _ in range(300):
        a, b = rand_dna(rng, rng.randint(8, 30)), rand_dna(rng, rng.randint(8, 30))
        for any_dg, three_dg in (core.dimer(a, b), core.hairpin(a)):
            assert any_dg <= 0.0 and three_dg <= 0.0          # stable pairing has a negative dG, none is exactly 0
            assert any_dg <= three_dg                         # the best of all dimers is at least as good as the best 3'-end one
        assert core.dimer(a, b) == pytest.approx(core.dimer(b, a))   # the pair is symmetric
    # A longer perfectly paired stretch is more stable than a shorter one; GC pairs more than AT.
    p = "AGCGGATAACAATTTCACACAGGA"
    assert core.dimer(p, rc(p))[0] < core.dimer(p[:12], rc(p[:12]))[0] < core.dimer(p[:8], rc(p[:8]))[0] < 0
    assert core.dimer("GCGCGCGC", rc("GCGCGCGC"))[0] < core.dimer("ATATATAT", rc("ATATATAT"))[0]
    # Same stem, longer loop: less stable (loop term grows with its size); more GC in the stem: more stable
    short_loop, long_loop = core.hairpin("GCGCGCTTTTGCGCGC")[0], core.hairpin("GCGCGC" + "T" * 12 + "GCGCGC")[0]
    assert short_loop < long_loop < 0
    assert core.hairpin("GCGCGCTTTTGCGCGC")[0] < core.hairpin("ATATATTTTTATATAT")[0]
    # A hairpin away from the 3' end does not count as a 3' hairpin
    inner = core.hairpin("GCGCGCTTTTGCGCGC" + "A" * 10)
    assert inner[0] < -3 and inner[1] == 0.0
    # More salt: more stable (more negative)
    assert core.hairpin("GCGCGCTTTTGCGCGC", na_mm=200)[0] < core.hairpin("GCGCGCTTTTGCGCGC", na_mm=20)[0]
    assert core.dimer(pal := "GAATTCGAATTC", pal, na_mm=200)[0] < core.dimer(pal, pal, na_mm=20)[0]
    # --- exact energies by hand from the SL98 table above ---
    # dG37 = dH - 310.15 * (dS + 0.368 * (N - 1) * ln(Na / 1 M)) / 1000; dimers add the initiation of both ends, hairpins
    # the SantaLucia & Hicks 2004 loop term (3: 3.5, 4: 3.5, 5: 3.3, 6: 4.0 kcal/mol). No Mg, no dNTP: Na_eq = Na.
    def stem_by_hand(stem, na_mm, init):
        dh = sum(SL98[a + b][0] for a, b in zip(stem, stem[1:]))
        ds = sum(SL98[a + b][1] for a, b in zip(stem, stem[1:]))
        if init:
            for end in (stem[0], stem[-1]):
                h, s = SL98_INIT["AT" if end in "AT" else "GC"]
                dh, ds = dh + h, ds + s
        return dh - 310.15 * (ds + 0.368 * (len(stem) - 1) * math.log(na_mm / 1000)) / 1000
    salt = dict(mg_mm=0, dntp_mm=0)
    # GCGCGCGC with itself: 4 GC + 3 CG steps and two GC ends: dH -70.8, dS -184.8, dG -70.8 + 310.15 * 0.1848 = -13.48428
    assert stem_by_hand("GCGCGCGC", 1000, True) == pytest.approx(-13.48428, abs=1e-9)
    assert core.dimer("GCGCGCGC", "GCGCGCGC", na_mm=1000, **salt) == pytest.approx((-13.48428, -13.48428), abs=1e-6)
    # Hairpin: 6 bp GCGCGC stem (dH -50.6, dS -127.6: -11.02486) + loop of 4 (3.5) = -7.52486, at 1 M and at 50 mM Na
    stem6 = stem_by_hand("GCGCGC", 1000, False)
    assert stem6 == pytest.approx(-11.02486, abs=1e-9)
    assert core.hairpin("GCGCGCTTTTGCGCGC", na_mm=1000, **salt) == pytest.approx((stem6 + 3.5,) * 2, abs=1e-6)
    assert stem_by_hand("GCGCGC", 50, False) + 3.5 == pytest.approx(-5.815267, abs=1e-6)   # dS += 0.368 * 5 * ln(0.05)
    assert core.hairpin("GCGCGCTTTTGCGCGC", na_mm=50, **salt) == pytest.approx((-5.815267,) * 2, abs=1e-5)
    # 3' flag: one or two unpairable bases after the stem leave no 3' hairpin; + GC makes a new one ending at the 3' base
    # (G0..C5 pair the last six bases G12..C17, loop TTTTGC of 6: -11.02486 + 4.0)
    for tail in ("A", "CC", "A" * 10):
        assert core.hairpin("GCGCGCTTTTGCGCGC" + tail, na_mm=1000, **salt) == pytest.approx((stem6 + 3.5, 0.0), abs=1e-6), tail
    assert core.hairpin("GCGCGCTTTTGCGCGCGC", na_mm=1000, **salt) == pytest.approx((stem6 + 3.5, stem6 + 4.0), abs=1e-6)


def test_gc_bins_match_an_exact_count():
    """GC bins: all G / all A / case / N, also zoomed out on blocks; block counts match the exact count up to block rounding;
    windows past the end of a circle, empty bins and more bins than bases."""
    # --- gc bins all g all a and case ---
    assert Document.new_sequence("g", "G" * 100).gc_bins(0, 100, 4) == [1.0] * 4
    assert Document.new_sequence("a", "A" * 100).gc_bins(0, 100, 4) == [0.0] * 4
    assert Document.new_sequence("c", "gcgcgcgc").gc_bins(0, 8, 2) == [1.0, 1.0], "lower case counts"
    assert Document.new_sequence("n", "NNGC").gc_bins(0, 4, 1) == [0.5], "N is not G or C, but it is in the denominator"
    # Zoomed out the blocks are used (bins wider than 4 * 1024 bases): same exact answers for the uniform cases, also when the
    # last block is not full.
    assert Document.new_sequence("g", "G" * 10_000).gc_bins(0, 10_000, 1) == [1.0]
    assert Document.new_sequence("a", "A" * 10_001).gc_bins(0, 10_001, 2) == [0.0, 0.0]
    # Mixed: half-and-half halves
    half = Document.new_sequence("h", "A" * 5120 + "G" * 5120)       # the halves meet on a block edge (5 * 1024)
    assert half.gc_bins(0, 10_240, 2) == [0.0, 1.0]
    # --- gc bins from blocks match the exact count up to block rounding ---
    rng = random.Random(9)
    seq = rand_dna(rng, 100_000)
    doc = Document.new_sequence("r", seq)
    exact = [sum(c in "GC" for c in seq[i * 25_000:(i + 1) * 25_000]) / 25_000 for i in range(4)]
    from_blocks = doc.gc_bins(0, 100_000, 4)             # bins of 25 000 > 4 blocks: cached block sums, edges rounded to 1024
    assert from_blocks == pytest.approx(exact, abs=0.02)
    narrow = doc.gc_bins(0, 100_000, 100)               # bins of 1000 < 4 blocks: read base by base, exact
    assert narrow == pytest.approx([sum(c in "GC" for c in seq[i * 1000:(i + 1) * 1000]) / 1000 for i in range(100)], abs=1e-6)
    # --- edges: an end past the sequence, empty bins, more bins than bases ---
    circ = Document.new_sequence("c", "G" * 10 + "A" * 10, True)
    # The bins never wrap round the origin of a circular sequence: an end past the length is cut to the length.
    assert circ.gc_bins(0, 1000, 2) == circ.gc_bins(0, 20, 2) == [1.0, 0.0]
    assert circ.gc_bins(5, 15, 1) == [0.5]
    # More bins than bases: the bins that hold no base are NaN, the others are the base itself.
    tiny = Document.new_sequence("t", "GA").gc_bins(0, 2, 6)
    assert len(tiny) == 6
    assert [x for x in tiny if not math.isnan(x)] == [1.0, 0.0] and sum(math.isnan(x) for x in tiny) == 4
    # An empty or inverted window, or no bins at all
    assert all(math.isnan(x) for x in circ.gc_bins(7, 7, 3)) and len(circ.gc_bins(7, 7, 3)) == 3
    assert all(math.isnan(x) for x in circ.gc_bins(15, 5, 2))
    assert circ.gc_bins(0, 20, 0) == []


# ---- ORFs ----------------------------------------------------------------------------------------------------------------

STOPS = {"TAA", "TAG", "TGA"}


def orfs_by_hand(seq, min_codons):
    """ORFs the way the docs describe them: ATG to stop, six frames, the first ATG after the previous stop; (start, end, strand)."""
    out = []
    n = len(seq)
    for strand, s in ((1, seq), (-1, rc(seq))):
        for frame in range(3):
            start = None
            for i in range(frame, n - 2, 3):
                codon = s[i:i + 3]
                if start is None and codon == "ATG":
                    start = i
                elif codon in STOPS:
                    if start is not None and (i + 3 - start) // 3 >= min_codons:
                        out.append((start, i + 3, 1) if strand == 1 else (n - i - 3, n - start, -1))
                    start = None
    return sorted(out, key=lambda o: o[0])


def circ_orfs_by_hand(seq, min_codons):
    """ORFs of a ring: from every stop, walk back codon by codon round the ring to the furthest ATG before the previous
    in-frame stop (an ORF is at most one lap). Start in [0, n), end past n when it crosses the origin; (start, end, strand)."""
    n, out = len(seq), set()
    for strand, s in ((1, seq), (-1, rc(seq))):
        ring = s + s
        codon = lambda i: ring[i % n:i % n + 3]
        for p in range(n):
            if codon(p) not in STOPS:
                continue
            best = None
            for k in range(1, n // 3 + 1):
                c = codon(p - 3 * k)
                if c in STOPS:
                    break
                if c == "ATG" and 3 * k + 3 <= n:
                    best = k
            if best is not None and best + 1 >= min_codons:
                a, size = (p - 3 * best) % n, 3 * best + 3
                b = (-(a + size)) % n   # bottom strand: the same bases counted from the top-strand start
                out.add((a, a + size, 1) if strand == 1 else (b, b + size, -1))
    return sorted(out)


def test_orfs_match_a_brute_force_scan():
    """ORFs: minimum length, stop coordinates, both strands, every stop codon; nested, back-to-back, open, missing start or
    stop; agreement with a brute-force scan on random sequences; ORFs across the origin of a circle, on both strands and
    on random rings."""
    # --- orfs minimum length stop coordinates and strands ---
    #            0  3   6   9   12
    seq = "CC" + "ATG" + "AAA" + "AAA" + "TAA" + "CC"   # 4 codons (ATG AAA AAA TAA) at 2..14
    doc = Document.new_sequence("o", seq)
    assert doc.orfs(4) == [(2, 14, 1)]
    assert doc.orfs(5) == [], "minimum length counts the stop codon: 4 codons are not enough for 5"
    assert doc.orfs(1) == [(2, 14, 1)]
    # The end includes the stop codon, the start is the ATG
    s, e, strand = doc.orfs(4)[0]
    assert seq[s:s + 3] == "ATG" and seq[e - 3:e] in STOPS
    # The same ORF on the bottom strand has the same coordinates on the top strand
    assert Document.new_sequence("o", rc(seq)).orfs(4) == [(2, 14, -1)]
    # Each of the three stop codons ends an ORF; lower case is read
    for stop in STOPS:
        d = Document.new_sequence("o", ("ccatgaaa" + stop.lower() + "cc"))
        assert d.orfs(3) == [(2, 11, 1)], stop
    # --- orfs nested back to back open and missing start or stop ---
    # ATG ATG AAA TAA: the longest ORF (first ATG) is reported once, not the inner one too
    assert Document.new_sequence("n", "ATGATGAAATAA").orfs(1) == [(0, 12, 1)]
    # Two ORFs in the same frame, separated by a stop
    two = "ATGAAATAA" + "ATGCCCTGA"
    assert Document.new_sequence("n", two).orfs(1) == [(0, 9, 1), (9, 18, 1)]
    # No ATG: nothing, even with stops and many codons
    assert Document.new_sequence("n", "AAA" * 30 + "TAA").orfs(1) == []
    # ATG but no stop to the end of the sequence: not an ORF
    assert Document.new_sequence("n", "ATG" + "AAA" * 30).orfs(1) == []
    # A stop before any ATG does not make an ORF
    assert Document.new_sequence("n", "TAATAGTGA").orfs(1) == []
    # --- orfs agree with a brute force scan on random sequences ---
    rng = random.Random(21)
    for n in (500, 3001):
        seq = rand_dna(rng, n)
        for min_codons in (1, 10, 30):
            assert Document.new_sequence("r", seq).orfs(min_codons) == orfs_by_hand(seq, min_codons)
    # And each ORF is a real one: codons from ATG to the first in-frame stop, on its own strand
    seq = rand_dna(rng, 6000)
    found = Document.new_sequence("r", seq).orfs(40)
    assert found
    for s, e, strand in found:
        orf = seq[s:e] if strand == 1 else rc(seq[s:e])
        codons = [orf[i:i + 3] for i in range(0, len(orf), 3)]
        assert len(orf) % 3 == 0 and codons[0] == "ATG" and codons[-1] in STOPS and not STOPS & set(codons[:-1])
    # --- orf across the origin of a circular sequence ---
    # ATG AAA AAA TAA read through the origin (the last 6 bases of the ORF are at the start of the file):
    # the filler is all C, which holds no ATG and no stop
    seq = "AAATAA" + "C" * 30 + "ATGAAA"       # circular: ...ATG AAA | AAA TAA...
    # ATG at 36, n = 42: the ORF ends 6 bases past the origin, at 42 + 6 = 48 (end past the length, like a wrapping feature)
    assert Document.new_sequence("c", seq, True).orfs(4) == [(36, 48, 1)]
    # The same ring written the other way round: the same bases, now on the bottom strand
    assert Document.new_sequence("c", rc(seq), True).orfs(4) == [(36, 48, -1)]
    # An ATG just after the origin, inside an ORF that started before it, is not a second ORF: one stop, one ORF, the longest
    # (this ring used to give [(0, 9, 1), (39, 54, 1)])
    nested = "ATGAAATAA" + "C" * 30 + "ATGAAA"
    assert Document.new_sequence("c", nested, True).orfs(3) == [(39, 54, 1)]
    assert Document.new_sequence("c", rc(nested), True).orfs(3) == [(36, 51, -1)]
    # Random rings (n not a multiple of 3, so the frame changes at the origin): exactly the ORFs of circ_orfs_by_hand, and
    # each is a real ORF read round the ring.
    for _ in range(3):
        seq = rand_dna(rng, 3001)
        ring = seq + seq
        for min_codons in (1, 10, 30):
            found = Document.new_sequence("r", seq, True).orfs(min_codons)
            assert sorted(circ_orfs_by_hand(seq, min_codons)) == sorted(found), min_codons
            for s, e, strand in found:
                assert 0 <= s < 3001 and e - s < 3001 and (e - s) // 3 >= min_codons
                orf = ring[s:e] if strand == 1 else rc(ring[s:e])
                codons = [orf[i:i + 3] for i in range(0, len(orf), 3)]
                assert len(orf) % 3 == 0 and codons[0] == "ATG" and codons[-1] in STOPS and not STOPS & set(codons[:-1])


def _digest_fragments(cut_positions, n, circular):
    """Fragment sizes from the top-strand cut positions of one digest."""
    cuts = sorted(set(cut_positions))
    if not cuts:
        return [n]
    if circular:
        return sorted([b - a for a, b in zip(cuts, cuts[1:])] + [cuts[0] + n - cuts[-1]])
    edges = [0] + [c for c in cuts if 0 < c < n] + [n]
    return sorted(b - a for a, b in zip(edges, edges[1:]))


def _biopython_cuts(name, seq, circular):
    from Bio import Restriction
    n = len(seq)
    found = getattr(Restriction, name).search(Seq(seq), linear=not circular)
    # Biopython: 1-based position of the first base after the cut; ours: the 0-based gap, the same number minus one
    return sorted({(p - 1) % n if circular else p - 1 for p in found})


def test_restriction_cuts_and_digests_match_biopython():
    """Restriction cuts: overhangs on a hand-made sequence and every common enzyme planted in both orientations; both
    strands, sites across the origin, degenerate (IUPAC) sites; digests of the real plasmid and of random circles with every
    enzyme planted equal Biopython's Restriction module."""
    # --- cuts and overhangs on a hand made sequence ---
    common = {e.name: e for e in enzymes.COMMON}
    # EcoRI G^AATTC: 5' overhang AATT; KpnI GGTAC^C: 3' overhang GTAC; EcoRV GAT^ATC: blunt
    s = "A" * 10 + "GAATTC" + "A" * 10 + "GGTACC" + "A" * 10 + "GATATC" + "A" * 10
    cuts = enzymes.find_cuts(s, False)
    assert cuts["EcoRI"] == [(11, 15)] and s[11:15] == "AATT"      # top cut after G, bottom cut 4 further: the overhang between
    assert cuts["KpnI"] == [(31, 27)] and s[27:31] == "GTAC"       # top cut after GGTAC, bottom 4 earlier: a 3' overhang
    assert cuts["EcoRV"] == [(45, 45)] and common["EcoRV"].blunt
    # Every common enzyme on a planted site, in both orientations: the cut offsets declared for it, mirrored for the other strand
    pad = "A" * 30
    for e in enzymes.COMMON:
        top, bottom = e.cuts[0]
        L = len(e.site)
        assert enzymes.find_cuts(pad + e.site + pad, False, [e])[e.name][0] == (30 + top, 30 + bottom), e.name
        if enzymes.revcomp_site(e.site) != e.site:
            assert enzymes.find_cuts(pad + rc(e.site) + pad, False, [e])[e.name] == [(30 + L - bottom, 30 + L - top)], e.name
    # --- both strands and sites across the origin ---
    bsai = next(e for e in enzymes.COMMON if e.name == "BsaI")
    pad = "A" * 20
    # BsaI GGTCTC(1/5): cuts 1 base after the site on the top strand, 5 on the bottom
    fwd = pad + "GGTCTC" + pad
    rev = pad + rc("GGTCTC") + pad               # GAGACC: the site on the bottom strand, cut upstream, mirrored
    assert enzymes.find_cuts(fwd, False, [bsai])["BsaI"] == [(27, 31)]
    assert enzymes.find_cuts(rev, False, [bsai])["BsaI"] == [(15, 19)]
    assert enzymes.find_cuts(fwd, False, [bsai])["BsaI"][0][1] - 27 == 4, "4 base 5' overhang"
    # BsaI: the cuts of the bottom-strand site are mirrored copies of the top-strand ones
    # Across the origin of a circle (GAATTC split) the site is found, on a line it is not
    circ = "ATTC" + "A" * 20 + "GA"
    assert enzymes.find_cuts(circ, True)["EcoRI"] == [(25, 3)]   # the top cut is just after the G at the end of the file
    assert "EcoRI" not in enzymes.find_cuts(circ, False)
    # site_spans: both strands, once for a palindrome, across the origin the end passes the length
    ecori = next(e for e in enzymes.COMMON if e.name == "EcoRI")
    assert enzymes.site_spans(circ, True, ecori) == [(24, 30, 1)]
    spans = enzymes.site_spans(fwd + rev, False, bsai)
    assert [(a, b - a, st) for a, b, st in spans] == [(20, 6, 1), (66, 6, -1)]
    # a type IIS cut that falls beyond the end of a linear sequence is not a cut; on a circle it wraps round
    assert "BsaI" not in enzymes.find_cuts("GGTCTC", False, [bsai])
    assert enzymes.find_cuts("GGTCTCA", True, [bsai])["BsaI"] == [(0, 4)], "on a circle it wraps (7 mod 7 = 0, 11 mod 7 = 4)"
    # --- degenerate sites and both sided cutters ---
    # BsmBI CGTCTC(1/5) like BsaI; AccI GT^MKAC (M = A/C, K = G/T); a degenerate N spacer
    bsmbi = next(e for e in enzymes.COMMON if e.name == "BsmBI")
    assert bsmbi.site == "CGTCTC" and bsmbi.cuts == ((7, 11),)
    s = "T" * 10 + "CGTCTC" + "T" * 10
    assert enzymes.find_cuts(s, False, [bsmbi])["BsmBI"] == [(17, 21)]
    acci = Enzyme("AccI", "GTMKAC", ((2, 4),), False)
    ok = ["GTAGAC", "GTCGAC", "GTATAC", "GTCTAC"]
    notok = ["GTGGAC", "GTTGAC", "GTAGCC"]
    for site in ok:
        assert enzymes.find_cuts("AAA" + site + "AAA", False, [acci]) == {"AccI": [(5, 7)]}, site
    for site in notok:
        assert enzymes.find_cuts("AAA" + site + "AAA", False, [acci]) == {}, site
    # The degenerate site is read on the other strand too
    assert enzymes.find_cuts("AAA" + rc("GTAGAC") + "AAA", False, [acci]) != {}
    # N in the middle of a site: BglI GCCNNNN^NGGC (3' overhang)
    bgli = Enzyme("BglI", "GCCNNNNNGGC", ((7, 4),), False)
    assert enzymes.find_cuts("TT" + "GCCATCGTGGC" + "TT", False, [bgli]) == {"BglI": [(9, 6)]}
    # --- the fragment sizes the app shows (tools.fragments) ---
    # One cut in a circle opens it: one fragment of the whole length, from the cut back round to it
    assert tools.fragments(100, True, {10: ["A"]}) == [(100, 10, 10, ["A"], ["A"])]
    assert tools.fragments(100, True, {10: ["A"], 60: ["B"]}) == [(50, 10, 60, ["A"], ["B"]), (50, 60, 10, ["B"], ["A"])]
    assert tools.fragments(100, True, {90: ["A"], 20: ["B"]}) == [(70, 20, 90, ["B"], ["A"]), (30, 90, 20, ["A"], ["B"])]
    # Two cuts in a line: three pieces, the ends of the line bounding the outer two
    assert tools.fragments(100, False, {10: ["A"], 60: ["B"]}) == [(10, 0, 10, [], ["A"]), (50, 10, 60, ["A"], ["B"]),
                                                                   (40, 60, 100, ["B"], [])]
    assert tools.fragments(100, False, {}) == [(100, 0, 100, [], [])]

    def digest_sizes(seq, circular, e):
        """Fragment sizes of a single digest, as the Digest tool builds them: {top-strand cut: [enzyme]} into tools.fragments."""
        cuts = {t: [e.name] for t, _ in enzymes.find_cuts(seq, circular, [e]).get(e.name, [])}
        return sorted(size for size, *_ in tools.fragments(len(seq), circular, cuts))

    def biopython_digest(name, seq, circular):
        from Bio import Restriction
        return sorted(len(f) for f in getattr(Restriction, name).catalyse(Seq(seq), linear=not circular))

    # --- the real plasmid, every common enzyme, against Biopython ---
    if PLASMID.exists():   # (a missing sample file only leaves this part out)
        doc = Document(str(PLASMID))
        seq, n = doc.sequence(0, None), doc.length
        assert doc.circular
        checked = 0
        for e in enzymes.COMMON:        # every common enzyme, 20-odd of which cut this plasmid
            ours = sorted({t for t, _ in enzymes.find_cuts(seq, True, [e]).get(e.name, [])})
            theirs = _biopython_cuts(e.name, seq, True)
            assert ours == theirs, e.name
            assert _digest_fragments(ours, n, True) == _digest_fragments(theirs, n, True)
            assert sum(_digest_fragments(ours, n, True)) == n
            checked += bool(ours)
            assert digest_sizes(seq, True, e) == biopython_digest(e.name, seq, True), e.name
        assert checked >= 10, "at least 10 enzymes actually cut this plasmid"
        # Opened as a line (cut positions only inside the sequence)
        lin = seq
        for e in enzymes.COMMON[:20]:
            ours = sorted({t for t, _ in enzymes.find_cuts(lin, False, [e]).get(e.name, [])})
            assert ours == _biopython_cuts(e.name, lin, False), e.name
            assert digest_sizes(lin, False, e) == biopython_digest(e.name, lin, False), e.name
    # --- every enzyme planted in random circles, one across the origin, against Biopython ---
    rng = random.Random(33)
    n = 4000
    for circle in range(5):
        seq = rand_dna(rng, n)
        # plant every enzyme's site in a random orientation at a random place
        for e in enzymes.COMMON:
            site = e.site if rng.random() < .5 else rc(e.site)
            p = rng.randrange(n - len(site))
            seq = seq[:p] + site + seq[p + len(site):]
        # and one across the origin: the end of the file holds the start of the site, the start of the file its end. EcoRI and
        # BsaI take turns: both at the same origin, the second would overwrite the first.
        e, k = ((enzymes.COMMON[14], 3), (enzymes.COMMON[9], 2))[circle % 2]   # EcoRI, BsaI
        assert e.name == ("EcoRI", "BsaI")[circle % 2]
        seq = e.site[k:] + seq[len(e.site) - k:n - k] + e.site[:k]
        assert len(seq) == n and seq[-k:] + seq[:len(e.site) - k] == e.site
        across = (n - k + e.cuts[0][0]) % n   # the top-strand cut of the site that starts k bases before the end
        assert across in {t for t, _ in enzymes.find_cuts(seq, True, [e])[e.name]}, (e.name, "cut across the origin")
        cutting = 0
        for e in enzymes.COMMON:
            for circular in (True, False):
                ours = sorted({t for t, _ in enzymes.find_cuts(seq, circular, [e]).get(e.name, [])})
                assert ours == _biopython_cuts(e.name, seq, circular), (e.name, circular)
                # the sizes the app shows, against the hand formula on Biopython's cuts (not catalyse: on a ring whose
                # site crosses the origin Biopython's catalyse returns nonsense, e.g. [0, 4517] for a 4000 bp BsaI digest)
                assert digest_sizes(seq, circular, e) == _digest_fragments(_biopython_cuts(e.name, seq, circular), n, circular)
                if circular:
                    assert sum(_digest_fragments(ours, n, True)) == n
                    cutting += bool(ours)
        assert cutting >= 40, "nearly every planted enzyme cuts (a later planting can cover part of an earlier one)"


# ---- find_parts ----------------------------------------------------------------------------------------------------------


def mutate(seq, positions):
    out = list(seq)
    for p in positions:
        out[p] = {"A": "C", "C": "G", "G": "T", "T": "A"}[out[p]]
    return "".join(out)



def test_find_parts_identity_coding_and_repeated_parts():
    """find_parts: forward and reverse hits and the identity thresholds; short parts need an exact match; a circle wraps;
    the coding flag and min_protein; overlapping, repeated and tandem parts."""
    # --- find parts forward reverse and the identity thresholds ---
    rng = random.Random(41)
    part = rand_dna(rng, 100)
    left, right = rand_dna(rng, 300), rand_dna(rng, 300)
    # mismatches 15 bases apart so that exact 12-mers remain for the seeds
    spaced = [8, 24, 40, 56, 72, 88]
    for strand, make in ((1, lambda s: s), (-1, rc)):
        assert core.find_parts(left + make(part) + right, False, [part], [False]) == [(0, 300, 400, strand, 1.0, False)]
        # 4 mismatches: identity 0.96 >= 0.95 found; 6: 0.94 is below 0.95 so no hit...
        four = core.find_parts(left + make(mutate(part, spaced[:4])) + right, False, [part], [False])
        assert len(four) == 1 and four[0][1:4] == (300, 400, strand) and four[0][4] == pytest.approx(0.96)
        six_seq = left + make(mutate(part, spaced)) + right
        assert core.find_parts(six_seq, False, [part], [False]) == []
        # ... until the threshold is lowered
        six = core.find_parts(six_seq, False, [part], [False], min_dna=0.9)
        assert len(six) == 1 and six[0][4] == pytest.approx(0.94) and six[0][1:4] == (300, 400, strand)
        # and a very strict threshold drops the imperfect one
        assert core.find_parts(left + make(mutate(part, spaced[:1])) + right, False, [part], [False], min_dna=0.995) == []
    # --- find parts short parts need an exact match and circular wraps ---
    rng = random.Random(42)
    short = "GAATTCAC"          # shorter than the 12-letter index words: exact search only
    host = rand_dna(rng, 200, "AC") + short + rand_dna(rng, 200, "AC")
    assert [h[1:4] for h in core.find_parts(host, False, [short], [False])] == [(200, 208, 1)]
    assert core.find_parts(host.replace(short, mutate(short, [3])), False, [short], [False], min_dna=0.5) == []
    # Across the origin of a circle: found with a start near the end and an end past the length; on a line, not
    part = rand_dna(rng, 60)
    circ = part[20:] + rand_dna(rng, 500) + part[:20]
    n = len(circ)
    assert core.find_parts(circ, True, [part], [False]) == [(0, n - 20, n + 40, 1, 1.0, False)]
    assert core.find_parts(circ, False, [part], [False]) == []
    # --- find parts coding flag and min protein ---
    table = CodonTable.unambiguous_dna_by_id[1].forward_table
    by_aa = {}
    for codon, aa in table.items():
        by_aa.setdefault(aa, []).append(codon)
    rng = random.Random(43)
    protein = "M" + "".join(rng.choice([a for a in by_aa if len(by_aa[a]) > 1 and a != "M"]) for _ in range(59))
    original = "".join(rng.choice(by_aa[a]) for a in protein) + "TAA"
    # Recoded: every codon replaced by a different synonymous one
    recoded = "".join(rng.choice([c for c in by_aa[a] if c != original[3 * i:3 * i + 3]] or by_aa[a]) for i, a in enumerate(protein)) + "TAG"
    identity_dna = sum(a == b for a, b in zip(original, recoded)) / len(original)
    assert identity_dna < 0.9, "the recoded copy must be far apart as DNA for this test to mean anything"
    assert str(Seq(recoded).translate()) == str(Seq(original).translate())
    host = rand_dna(rng, 400, "AT") + recoded + rand_dna(rng, 400, "AT")
    # coding: found as a protein, identity 1.0, the whole part length including the stop
    hits = core.find_parts(host, False, [original], [True])
    assert [(h[1], h[2], h[3], h[4], h[5]) for h in hits] == [(400, 400 + len(original), 1, 1.0, True)]
    # not coding: it is a DNA search only and the recoded copy is not found
    assert core.find_parts(host, False, [original], [False]) == []
    # on the other strand
    rev = core.find_parts(rc(host), False, [original], [True])
    assert [(h[3], h[4], h[5]) for h in rev] == [(-1, 1.0, True)] and rev[0][1] == len(host) - 400 - len(original)
    # min_protein: change 3 of 59 amino acids (5 %) -> found at 0.9 and with the identity reported, not at 0.97
    codons = [recoded[i:i + 3] for i in range(0, len(recoded), 3)]
    for i in (10, 30, 50):     # spread out so that 5-letter words survive
        other = next(a for a in by_aa if a != protein[i] and a != "M")
        codons[i] = by_aa[other][0]
    changed = rand_dna(rng, 400, "AT") + "".join(codons) + rand_dna(rng, 400, "AT")
    found = core.find_parts(changed, False, [original], [True])
    assert len(found) == 1 and found[0][5] is True
    prot_len = len(protein)           # the stop is dropped from the protein before comparing
    assert found[0][4] == pytest.approx(1 - 3 / prot_len)
    assert core.find_parts(changed, False, [original], [True], min_protein=0.97) == []
    assert len(core.find_parts(changed, False, [original], [True], min_protein=0.9)) == 1
    # --- find parts overlapping and repeated parts ---
    rng = random.Random(44)
    a = rand_dna(rng, 120)
    b = a[60:] + rand_dna(rng, 60)               # b starts half way through a: they overlap on the template
    host = rand_dna(rng, 100) + a + rand_dna(rng, 60) + rand_dna(rng, 100)
    hits = core.find_parts(host, False, [a, b[:60]], [False, False])
    # two different parts that overlap are both reported (only the same part on the same strand competes with itself)
    assert sorted((h[0], h[1], h[2]) for h in hits) == [(0, 100, 220), (1, 160, 220)]
    # The same part twice, apart: two hits, sorted by position; and once on each strand
    host = rand_dna(rng, 50) + a + rand_dna(rng, 50) + rc(a) + rand_dna(rng, 50) + a + rand_dna(rng, 50)
    hits = core.find_parts(host, False, [a], [False])
    assert [(h[1], h[3]) for h in hits] == [(50, 1), (220, -1), (390, 1)]
    # Tandem copies with an offset overlap a good alignment: only the best of overlapping ones of the same strand is kept
    unit = rand_dna(rng, 40)
    tandem = rand_dna(rng, 30) + unit * 3 + rand_dna(rng, 30)
    spans = [(h[1], h[2]) for h in core.find_parts(tandem, False, [unit], [False])]
    assert spans == [(30, 70), (70, 110), (110, 150)]
