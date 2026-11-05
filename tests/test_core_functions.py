"""The Rust core through its Python API (`genomeedit_core`): melting temperature, hairpins and dimers, GC bins, ORFs, restriction
sites, find_parts, translation, local alignment, BAM reads / pileup / variants, GFF / BED tracks, feature locations and
indexed genome FASTA files.

Every test has an independent reference: Biopython (Tm, Restriction, translation), a hand count, a small model written
here (BAM pileup, ORFs, density bins), or a value the docs / code comments promise. Nothing here needs Qt.
"""
import itertools
import math
import random
import struct
import sys
import zlib
from pathlib import Path

import pytest
from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "app"))

import genomeedit_core as core
from genomeedit_core import Document

import enzymes
from enzymes import Enzyme

from sample_files import DNA_FILE as PLASMID   # a copy (see sample_files.py)
COMP = str.maketrans("ACGTacgt", "TGCAtgca")


COMP_BASE = dict(zip("ACGT", "TGCA"))


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


def test_melting_temperature():
    """Tm of three known primers against the SantaLucia 1998 formula and Biopython; Tm rises with length and GC, ignores case,
    reads U as T; salt, magnesium, dNTP and primer concentration move it the right way."""
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


# ---- GC bins -------------------------------------------------------------------------------------------------------------


def test_gc_bins_uniform_mixed_and_block_counts():
    """GC bins: all G / all A / case / N, also zoomed out on blocks; block counts match the exact count up to block rounding."""
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


def test_gc_bins_edges_end_past_the_sequence_empty_bins_and_more_bins_than_bases():
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


def test_orfs_lengths_strands_and_nesting():
    """ORFs: minimum length, stop coordinates, both strands, every stop codon; nested, back-to-back, open, missing start or stop."""
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


def test_orfs_agree_with_a_brute_force_scan_and_wrap_the_origin():
    """ORFs agree with a brute-force scan on random sequences (and each is a real one); an ORF across the origin of a circle."""
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
    orfs = Document.new_sequence("c", seq, True).orfs(4)
    assert len(orfs) == 1 and orfs[0][2] == 1


# ---- restriction sites ---------------------------------------------------------------------------------------------------


def test_cuts_overhangs_strands_origin_and_degenerate_sites():
    """Restriction cuts and overhangs on a hand-made sequence and every common enzyme planted; both strands and sites across
    the origin; degenerate sites (IUPAC) and cutters on both sides."""
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


def test_digest_on_the_real_plasmid_equals_biopython_restriction():
    if not PLASMID.exists():
        pytest.skip(f"missing {PLASMID.name}")
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
    assert checked >= 10, "at least 10 enzymes actually cut this plasmid"
    # Opened as a line (cut positions only inside the sequence)
    lin = seq
    for e in enzymes.COMMON[:20]:
        ours = sorted({t for t, _ in enzymes.find_cuts(lin, False, [e]).get(e.name, [])})
        assert ours == _biopython_cuts(e.name, lin, False), e.name


def test_digest_with_every_enzyme_planted_in_a_random_circle_equals_biopython():
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


def test_find_parts_identity_thresholds_short_parts_and_the_origin():
    """find_parts: forward and reverse hits and the identity thresholds; short parts need an exact match; a circle wraps."""
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


def test_find_parts_coding_and_repeated_parts():
    """find_parts: the coding flag and min_protein; overlapping, repeated and tandem parts."""
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


# ---- translation (through the protein search, which translates the three frames of both strands) -------------------------

ALL_CODONS = ["".join(c) for c in itertools.product("TCAG", repeat=3)]


def test_translation_of_all_codons_incomplete_codons_and_n():
    """Translation through the protein search: all 64 codons against the standard table; incomplete last codon; codons with N."""
    # --- translation of all 64 codons against the standard table ---
    table = CodonTable.unambiguous_dna_by_id[1]
    for codon in ALL_CODONS:
        aa = "*" if codon in table.stop_codons else table.forward_table[codon]
        assert str(Seq(codon).translate()) == aa
        # forward strand, frame 1: found with the exact amino acid, no mismatches, in the right place
        assert (0, 3, 1, 0, 1) in core.find_protein(codon, aa), codon
        # and the reverse strand: the codon's reverse complement on the top strand reads as the same amino acid
        assert (0, 3, -1, 0, 1) in core.find_protein(rc(codon), aa), codon
        # a different amino acid is no exact match there
        other = "A" if aa != "A" else "G"
        assert (0, 3, 1, 0, 1) not in core.find_protein(codon, other), codon
    # all 64 together as one 64-residue protein, in order
    dna = "".join(ALL_CODONS)
    protein = str(Seq(dna).translate())
    assert len(protein) == 64
    assert core.find_protein(dna, protein) == [(0, 192, 1, 0, 1)]
    # the same, in lower case
    assert core.find_protein(dna.lower(), protein) == [(0, 192, 1, 0, 1)]
    # --- translation incomplete last codon and codons with n ---
    dna = "ATGGCCAAAGG"          # M A K and a last incomplete codon GG
    forward = lambda hits: [h for h in hits if h[2] == 1]
    assert forward(core.find_protein(dna, "MAK")) == [(0, 9, 1, 0, 1)]
    assert forward(core.find_protein(dna, "MAKG")) == [], "the two leftover bases are not a codon"
    # Other frames read the same sequence shifted: TGG CCA AAG -> W P K
    assert (1, 10, 1, 0, 2) in core.find_protein(dna, "WPK")
    # A codon with an N is X: only an X in the query (or a mismatch allowed) matches it
    n_dna = "ATGNNNGCC"                                         # M X A
    assert forward(core.find_protein(n_dna, "MXA")) == [(0, 9, 1, 0, 1)]
    assert forward(core.find_protein(n_dna, "MAA")) == []
    assert forward(core.find_protein(n_dna, "MAA", max_mismatches=1)) == [(0, 9, 1, 1, 1)]


# ---- local alignment -----------------------------------------------------------------------------------------------------


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


def test_align_strands_indels_hits_and_query_longer_than_the_template():
    """Local alignment: reverse strand and coordinates; indels at their exact positions; max_hits, no hit and degenerate
    inputs; a query longer than the template and across the origin."""
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
    # A circle: the query runs through the origin, so end passes the length. On a line, only the longer half is found.
    t = rand_dna(rng, 1000)
    q = t[900:] + t[:100]
    h = core.align(t, q, circular=True)[0]
    assert (h[0], h[1], h[2], h[4]) == (900, 1100, 1, 1.0)
    h = core.align(t, rc(q), circular=True)[0]
    assert (h[0], h[1], h[2]) == (900, 1100, -1)
    lin = core.align(t, q, circular=False)
    assert lin and lin[0][1] <= 1000 and lin[0][1] - lin[0][0] == 100
    # a hit at the very start of a circle is reported once (not again as a copy past the end)
    assert len(core.align(t, t[:80], circular=True)) == 1


def test_align_overhangs_and_partial_matches():
    """Local alignment of partial matches: a query hanging off the end of a linear template, wrapping on a circle, whose ends
    match but whose middle does not, and an end overlap below the reporting threshold."""
    # --- align query hanging off the end of a linear template aligns only the overlap ---
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


# ---- BAM ---------------------------------------------------------------------------------------------------------------


def _reg2bin(beg, end):
    """The BAM index bin of [beg, end) (SAM spec section 5.3)."""
    end -= 1
    for shift, offset in ((14, 4681), (17, 585), (20, 73), (23, 9), (26, 1)):
        if beg >> shift == end >> shift:
            return offset + (beg >> shift)
    return 0


def _bgzf(data):
    """BGZF (blocked gzip) of `data`, with the empty end-of-file block the format asks for."""
    out = b""
    for i in range(0, len(data), 0xff00):
        chunk = data[i:i + 0xff00]
        comp = zlib.compressobj(6, zlib.DEFLATED, -15)
        body = comp.compress(chunk) + comp.flush()
        size = len(body) + 25            # block size minus 1: 18 header + body + 8 trailer - 1
        out += struct.pack("<BBBBIBBHBBHH", 31, 139, 8, 4, 0, 0, 255, 6, 66, 67, 2, size)
        out += body + struct.pack("<II", zlib.crc32(chunk), len(chunk))
    return out + bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")


CIGAR_CODES = {"M": 0, "I": 1, "D": 2, "N": 3, "S": 4, "H": 5, "P": 6, "=": 7, "X": 8}
SEQ_CODES = {c: i for i, c in enumerate("=ACMGRSVTWYHKDBN")}


def read(name, pos, cigar, seq, flag=0, mapq=60, qual=None, ref_id=0, mate_ref=-1, mate_pos=-1, tlen=0):
    """One read: `cigar` like [(5, "S"), (45, "M")]; `qual` a list of base qualities, or None for 'not stored'."""
    return dict(name=name, pos=pos, cigar=cigar, seq=seq, flag=flag, mapq=mapq, qual=qual, ref_id=ref_id,
                mate_ref=mate_ref, mate_pos=mate_pos, tlen=tlen)


def write_bam(path, refs, reads):
    """A coordinate-sorted BAM written with nothing but the standard library: `refs` = [(name, length)]; reads sorted here."""
    text = b"@HD\tVN:1.6\tSO:coordinate\n"
    head = b"BAM\1" + struct.pack("<i", len(text)) + text + struct.pack("<i", len(refs))
    for name, length in refs:
        head += struct.pack("<i", len(name) + 1) + name.encode() + b"\0" + struct.pack("<i", length)
    body = b""
    for r in sorted(reads, key=lambda r: (r["ref_id"] if r["ref_id"] >= 0 else 1 << 30, r["pos"])):
        cigar = b"".join(struct.pack("<I", n << 4 | CIGAR_CODES[op]) for n, op in r["cigar"])
        ref_len = sum(n for n, op in r["cigar"] if op in "MDN=X")
        seq = r["seq"]
        packed = bytearray((len(seq) + 1) // 2)
        for i, base in enumerate(seq):
            packed[i // 2] |= SEQ_CODES[base] << (4 if i % 2 == 0 else 0)
        qual = bytes(r["qual"]) if r["qual"] is not None else b"\xff" * len(seq)
        name = r["name"].encode() + b"\0"
        end = r["pos"] + max(ref_len, 1)
        rec = struct.pack("<iiBBHHHiiii", r["ref_id"], r["pos"], len(name), r["mapq"], _reg2bin(r["pos"], end), len(r["cigar"]),
                          r["flag"], len(seq), r["mate_ref"], r["mate_pos"], r["tlen"])
        rec += name + cigar + bytes(packed) + qual
        body += struct.pack("<i", len(rec)) + rec
    Path(path).write_bytes(_bgzf(head + body))


def build_read(name, pos, ops, ref, subs=None, ins="", clip="T", **kw):
    """A read made from the reference along `ops` ("5S20M3D10M2I15M" as a list of (n, op)): substituted bases `subs`
    {reference position: base}, inserted bases `ins` (split over the I operations), soft-clipped bases `clip`."""
    subs, seq, r, ins_at = subs or {}, "", pos, 0
    for n, op in ops:
        if op == "M":
            seq += "".join(subs.get(r + k, ref[r + k]) for k in range(n))
            r += n
        elif op == "D":
            r += n
        elif op == "I":
            seq += ins[ins_at:ins_at + n]
            ins_at += n
        elif op == "S":
            seq += clip * n
    return read(name, pos, ops, seq, **kw)


def model_walk(r):
    """('M', ref position, base, quality) / ('I', ref position it goes before, bases) / ('D', ref position) for a read."""
    pos, q, events = r["pos"], 0, []
    for n, op in r["cigar"]:
        if op in "M=X":
            for k in range(n):
                qual = 255 if r["qual"] is None else r["qual"][q + k]
                events.append(("M", pos + k, r["seq"][q + k], qual))
            pos, q = pos + n, q + n
        elif op == "I":
            events.append(("I", pos, r["seq"][q:q + n]))
            q += n
        elif op in "DN":
            events += [("D", pos + k) for k in range(n)]
            pos += n
        elif op == "S":
            q += n
    return events


def usable(r):
    return not r["flag"] & (0x4 | 0x400 | 0x200)           # unmapped, duplicate, failed QC: ignored


def model_pileup(reads, ref_id, start, end, min_baseq=0, min_mapq=0):
    """Hand-counted pileup: per position [A, C, G, T, N, deleted, insertion before]."""
    counts = [[0] * 7 for _ in range(end - start)]
    for r in reads:
        if r["ref_id"] != ref_id or not usable(r) or r["mapq"] < min_mapq:
            continue
        for ev in model_walk(r):
            if not start <= ev[1] < end:
                continue
            slot = counts[ev[1] - start]
            if ev[0] == "M" and (ev[3] == 255 or ev[3] >= min_baseq):
                slot["ACGT".index(ev[2]) if ev[2] in "ACGT" else 4] += 1
            elif ev[0] == "D":
                slot[5] += 1
            elif ev[0] == "I":
                slot[6] += 1
    return counts


def model_depth(reads, ref_id, length):
    depth = [0] * length
    for r in reads:
        if r["ref_id"] == ref_id and usable(r):
            for ev in model_walk(r):
                if ev[0] == "M" and ev[1] < length:
                    depth[ev[1]] += 1
    return depth


@pytest.fixture(scope="module")
def ref_a():
    return rand_dna(random.Random(61), 2000)


@pytest.fixture(scope="module")
def hand_bam(tmp_path_factory, ref_a):
    """One read of every kind on chrA (2000 bp) and chrB (500 bp)."""
    ref = ref_a
    swap = lambda p, avoid=(): next(b for b in "ACGT" if b != ref[p] and b not in avoid)
    reads = [
        build_read("perfect", 10, [(50, "M")], ref, qual=[40] * 50),
        build_read("reverse", 30, [(50, "M")], ref, subs={45: swap(45)}, flag=0x10, qual=[40] * 50),
        build_read("softclip", 100, [(5, "S"), (45, "M")], ref, clip="T", qual=[30] * 50),
        build_read("deletion", 200, [(20, "M"), (5, "D"), (25, "M")], ref, qual=[30] * 45),
        build_read("insertion", 300, [(20, "M"), (3, "I"), (27, "M")], ref, ins="GGG", qual=[30] * 50),
        read("unmapped", 310, [], "ACGTACGTAC", flag=0x4, mapq=0, qual=[30] * 10),
        build_read("duplicate", 100, [(50, "M")], ref, flag=0x400, qual=[40] * 50),
        build_read("qcfail", 100, [(50, "M")], ref, flag=0x200, qual=[40] * 50),
        build_read("lowmapq", 500, [(50, "M")], ref, subs={510: swap(510)}, mapq=5, qual=[40] * 50),
        build_read("lowbq", 500, [(50, "M")], ref, subs={520: swap(520)}, qual=[40] * 20 + [10] + [40] * 29),
        build_read("hardclip", 700, [(10, "H"), (40, "M"), (10, "H")], ref, qual=[40] * 40),
        build_read("pair1", 900, [(50, "M")], ref, flag=0x1 | 0x40 | 0x20, mate_ref=1, mate_pos=100, tlen=0, qual=[40] * 50),
        build_read("pair2", 1000, [(50, "M")], ref, flag=0x1 | 0x40 | 0x20, mate_ref=0, mate_pos=1150, tlen=200, qual=[40] * 50),
        read("chrB_read", 10, [(30, "M")], "ACGT" * 7 + "AC", ref_id=1, qual=[40] * 30),
    ]
    path = tmp_path_factory.mktemp("bam") / "hand.bam"
    write_bam(path, [("chrA", 2000), ("chrB", 500)], reads)
    return core.Bam(str(path)), reads, path


def test_bam_hand_made_reads_pileup_and_coverage(hand_bam, ref_a):
    """BAM with one read of every kind: sequences and the .bai; reads (soft clips, deletions, insertions, strand, ignored
    reads, pairing, windows); pileup and coverage against a hand count."""
    # --- bam sequences and index written next to the file ---
    bam, _, path = hand_bam
    assert bam.sequences() == [("chrA", 2000), ("chrB", 500)]
    assert path.with_name(path.name + ".bai").exists(), "an unindexed sorted BAM gets its .bai written"
    with pytest.raises(OSError):
        bam.reads("chrZ", 0, 10, "", 10)
    # --- bam reads soft clips deletions insertions strand and ignored reads ---
    bam, _, _ = hand_bam
    got, more = bam.reads("chrA", 0, 2000, ref_a, 1000)
    assert not more
    by_name = {r[4]: r for r in got}
    # unmapped, duplicate and QC-failed reads are not listed; hard-clipped one is (it reads the aligned part)
    assert set(by_name) == {"perfect", "reverse", "softclip", "deletion", "insertion", "lowmapq", "lowbq", "hardclip", "pair1", "pair2"}
    p = by_name["perfect"]
    assert (p[0], p[1], p[2], p[3], p[5], p[6], p[7]) == (10, 60, 1, 60, [], [], [])
    r = by_name["reverse"]
    assert (r[0], r[1], r[2]) == (30, 80, -1) and r[8] & 0x10
    assert len(r[5]) == 1 and r[5][0][0] == 45 and r[5][0][1] != ref_a[45]
    s = by_name["softclip"]
    assert (s[0], s[1], s[5]) == (100, 145, []), "the 5 clipped T are not aligned: the read starts at 100 and has no mismatch"
    d = by_name["deletion"]
    assert (d[0], d[1], d[7]) == (200, 250, [(220, 225)]), "20 + 5 deleted + 25 = 50 reference bases"
    i = by_name["insertion"]
    assert (i[0], i[1], i[6], i[5]) == (300, 347, [(320, 3)], []), "inserted bases are not mismatches"
    h = by_name["hardclip"]
    assert (h[0], h[1]) == (700, 740)
    # pairing: flags, the mate position (0-based), the mate on another sequence, the insert size
    a, b = by_name["pair1"], by_name["pair2"]
    assert (a[9], a[10], a[11]) == (100, True, 0)
    assert (b[9], b[10], b[11]) == (1150, False, 200)
    assert by_name["perfect"][9:] == (None, False, 0)
    # a window gets only the reads that overlap it; limit sets the flag
    assert [r[4] for r in bam.reads("chrA", 100, 130, ref_a[100:130], 100)[0]] == ["softclip"]
    assert [r[4] for r in bam.reads("chrA", 0, 200, ref_a[:200], 2)[0]] == ["perfect", "reverse"]
    assert bam.reads("chrA", 0, 200, ref_a[:200], 2)[1] is True
    assert bam.reads("chrA", 0, 200, ref_a[:200], 4)[1] is False, "exactly `limit` reads: nothing more"
    assert [r[4] for r in bam.reads("chrB", 0, 100, "A" * 100, 10)[0]] == ["chrB_read"]
    # the reference is compared case-insensitively
    assert bam.reads("chrA", 0, 200, ref_a[:200].lower(), 100)[0][1][5] == bam.reads("chrA", 0, 200, ref_a[:200], 100)[0][1][5]
    # --- bam pileup hand count on the designed reads ---
    bam, reads, _ = hand_bam
    # 100..150: the soft-clipped read (45 aligned bases from 100) only; duplicate and QC-failed reads at the same place are ignored
    pile = bam.pileup("chrA", 100, 150)
    assert [sum(c[:4]) for c in pile] == [1] * 45 + [0] * 5
    assert all(c[4:] == [0, 0, 0] for c in pile)
    for i in range(45):
        assert pile[i]["ACGT".index(ref_a[100 + i])] == 1
    # deletion: 220..225 counted as deleted (and in the depth), insertion counted before 320
    d = bam.pileup("chrA", 215, 230)
    assert [c[5] for c in d] == [0] * 5 + [1] * 5 + [0] * 5
    assert [sum(c[:6]) for c in d] == [1] * 15
    ins = bam.pileup("chrA", 315, 325)
    assert [c[6] for c in ins] == [0, 0, 0, 0, 0, 1, 0, 0, 0, 0]
    # everything, against the hand count
    assert bam.pileup("chrA", 0, 2000) == model_pileup(reads, 0, 0, 2000)
    # a window that runs past the end of the sequence is cut to it
    assert len(bam.pileup("chrA", 1990, 3000)) == 10
    # filters: the base with quality 10 is dropped at min_baseq 20, a mapq 5 read at min_mapq 20; the reads with no
    # filter-relevant values are untouched
    p0, pq, pm = bam.pileup("chrA", 500, 550), bam.pileup("chrA", 500, 550, min_baseq=20), bam.pileup("chrA", 500, 550, min_mapq=20)
    assert sum(sum(c[:4]) for c in p0) == 100
    assert sum(sum(c[:4]) for c in pq) == 99 and sum(p0[20][:4]) - sum(pq[20][:4]) == 1
    assert sum(sum(c[:4]) for c in pm) == 50 and pm == model_pileup(reads, 0, 500, 550, min_mapq=20)
    assert pq == model_pileup(reads, 0, 500, 550, min_baseq=20)
    assert sum(sum(c[:4]) for c in bam.pileup("chrA", 500, 550, min_baseq=41)) == 0
    assert sum(sum(c[:4]) for c in bam.pileup("chrA", 500, 550, min_mapq=61)) == 0
    # --- bam coverage bins equal a hand count ---
    bam, reads, _ = hand_bam
    depth = model_depth(reads, 0, 2000)
    # narrow bins (< 512 bases): counted exactly from the reads; mean depth of each bin
    for start, end, bins in ((0, 200, 20), (95, 155, 6), (200, 260, 60)):
        width = (end - start) // bins
        expect = [sum(depth[start + b * width:start + (b + 1) * width]) / width for b in range(bins)]
        assert bam.coverage("chrA", start, end, bins) == pytest.approx(expect, abs=1e-6), (start, end, bins)
    # wide bins (>= 4 blocks of 128): from the per-block totals computed once for the whole sequence
    assert not bam.has_depth("chrA")
    wide = bam.coverage("chrA", 0, 1024, 2)
    assert bam.has_depth("chrA")
    assert wide == pytest.approx([sum(depth[0:512]) / 512, sum(depth[512:1024]) / 512], abs=1e-6)
    assert bam.coverage("chrA", 0, 2048, 4) == pytest.approx([sum(depth[i * 512:(i + 1) * 512]) / 512 for i in range(4)], abs=1e-6)
    # an empty window or no bins
    assert bam.coverage("chrA", 50, 50, 3) == [0.0, 0.0, 0.0] and bam.coverage("chrA", 0, 100, 0) == []


@pytest.fixture(scope="module")
def random_bam(tmp_path_factory, ref_a):
    """300 reads of random position, length, strand, mapping quality, base qualities and CIGAR (with soft clips, deletions
    and insertions), with sequencing errors; plus a few that are ignored."""
    rng = random.Random(62)
    reads = []
    for i in range(300):
        pos = rng.randrange(0, 1900)
        ops = [(rng.randint(15, 60), "M")]
        if rng.random() < .15:
            ops += [(rng.randint(1, 6), "D"), (rng.randint(10, 30), "M")]
        if rng.random() < .15:
            ops += [(rng.randint(1, 4), "I"), (rng.randint(10, 30), "M")]
        if rng.random() < .3:
            ops = [(rng.randint(1, 8), "S")] + ops
        if rng.random() < .3:
            ops += [(rng.randint(1, 8), "S")]
        subs = {pos + rng.randrange(0, 60): rng.choice("ACGTN") for _ in range(rng.randint(0, 3))}
        ref_len = sum(n for n, op in ops if op in "MD")
        if pos + ref_len > 2000:
            continue
        n_bases = sum(n for n, op in ops if op in "MIS")
        flag = rng.choice([0, 0, 0x10, 0x10, 0x4, 0x400, 0x200]) if i % 10 == 0 else rng.choice([0, 0x10])
        if flag == 0x4:
            reads.append(read(f"r{i}", pos, [], rand_dna(rng, 20), flag=flag, mapq=0, qual=[30] * 20))
            continue
        r = build_read(f"r{i}", pos, ops, ref_a, subs={p: b for p, b in subs.items()}, ins=rand_dna(rng, 10), clip="C",
                       flag=flag, mapq=rng.choice([0, 10, 30, 60, 60, 60]),
                       qual=[rng.choice([2, 10, 20, 30, 40]) for _ in range(n_bases)])
        reads.append(r)
    path = tmp_path_factory.mktemp("bam") / "random.bam"
    write_bam(path, [("chrA", 2000)], reads)
    return core.Bam(str(path)), reads


def test_bam_random_reads_and_variants(random_bam, snp_bam, ref_a, hand_bam):
    """BAM pileup and coverage equal the model on random reads; variant calling: a planted SNP and the thresholds, the
    reference base, deletions and insertions, random reads against the model."""
    # --- bam pileup and coverage equal the model on random reads ---
    bam, reads = random_bam
    assert len(reads) > 250
    for min_baseq, min_mapq in ((0, 0), (20, 0), (0, 30), (30, 60), (41, 0)):
        assert bam.pileup("chrA", 0, 2000, min_baseq, min_mapq) == model_pileup(reads, 0, 0, 2000, min_baseq, min_mapq), (min_baseq, min_mapq)
        assert bam.pileup("chrA", 700, 913, min_baseq, min_mapq) == model_pileup(reads, 0, 700, 913, min_baseq, min_mapq)
    # keyword form too
    assert bam.pileup("chrA", 0, 100, min_baseq=20, min_mapq=10) == model_pileup(reads, 0, 0, 100, 20, 10)
    depth = model_depth(reads, 0, 2000)
    assert bam.coverage("chrA", 0, 2000, 100) == pytest.approx([sum(depth[i * 20:(i + 1) * 20]) / 20 for i in range(100)], abs=1e-6)
    assert bam.coverage("chrA", 0, 1024, 2) == pytest.approx([sum(depth[0:512]) / 512, sum(depth[512:1024]) / 512], abs=1e-6)
    # the depth the pileup reports (A..T, N, deleted) is not the coverage: deleted bases are in the first, not in the second
    assert any(c[5] for c in bam.pileup("chrA", 0, 2000))
    # --- variants planted snp at half is called and thresholds exclude it ---
    bam, reads, alt = snp_bam
    found = bam.variants("chrA", ref_a, min_depth=5, min_fraction=0.2)
    assert [v for v in found if v[0] in (120, 125, 130)] == [(120, "snp", ref_a[120], alt(120), 5, 10), (130, "snp", ref_a[130], alt(130), 3, 10)]
    assert not [v for v in found if v[0] == 125], "1 of 10 reads is 10 %: below the default 20 %"
    # the fraction: 50 % passes at 0.5, not at 0.6; 30 % passes at 0.3 and 0.2, not at 0.4
    at = lambda **kw: {v[0] for v in bam.variants("chrA", ref_a, **kw)}
    assert at(min_depth=5, min_fraction=0.5) == {120}
    assert at(min_depth=5, min_fraction=0.6) == set()
    assert at(min_depth=5, min_fraction=0.4) == {120}
    assert at(min_depth=5, min_fraction=0.05) == {120, 125, 130}
    # the depth: ten reads cover these bases, so 10 passes and 11 does not
    assert at(min_depth=10, min_fraction=0.5) == {120}
    assert at(min_depth=11, min_fraction=0.5) == set()
    # filters: with every base quality 40, min_baseq 41 leaves no depth at all
    assert bam.variants("chrA", ref_a, 1, 0.1, min_baseq=41) == []
    assert bam.variants("chrA", ref_a, 1, 0.1, min_mapq=61) == []
    # all of it equals the model, every position
    brief = [(p, r, n, d) for p, k, r, a, n, d in bam.variants("chrA", ref_a, min_depth=1, min_fraction=0.1) if k == "snp"]
    assert brief == model_snps(reads, ref_a, 1, 0.1)
    # --- variants reference base is compared and reported correctly ---
    bam, _, alt = snp_bam
    # Hand a reference that has the other base at 120: now the five reads that carry `alt` agree with it and the other five disagree
    swapped = ref_a[:120] + alt(120) + ref_a[121:]
    v = [x for x in bam.variants("chrA", swapped, min_depth=5, min_fraction=0.5) if x[0] == 120]
    assert v == [(120, "snp", alt(120), ref_a[120], 5, 10)]
    # the reference in lower case is the same reference
    assert bam.variants("chrA", ref_a.lower(), 5, 0.2) == bam.variants("chrA", ref_a, 5, 0.2)
    # a reference N at the position is no SNP at all, and a read N is not an alternative base
    with_n = ref_a[:120] + "N" + ref_a[121:]
    assert not [x for x in bam.variants("chrA", with_n, 5, 0.5) if x[0] == 120]
    with pytest.raises(OSError):
        bam.variants("chrZ", ref_a)
    # --- variants deletions and insertions and the random reads ---
    bam, _, _ = hand_bam
    found = bam.variants("chrA", ref_a, min_depth=1, min_fraction=0.5)
    got = {(v[0], v[1]): v for v in found}
    # the deletion read: 5 deleted bases merged into one variant (220..225), the reference bases listed
    assert got[(220, "del")] == (220, "del", ref_a[220:225], "", 1, 1)
    # the insertion read: GGG before 320
    assert got[(320, "ins")] == (320, "ins", "", "GGG", 1, 1)
    # the reverse read's substitution at 45: one of the two reads covering it (50 %)
    snp45 = [v for v in found if v[0] == 45]
    assert len(snp45) == 1 and snp45[0][1:3] == ("snp", ref_a[45]) and snp45[0][4:] == (1, 2)
    # the soft-clipped T bases are not variants: no call before 100 from that read
    assert not [v for v in found if 95 <= v[0] < 100]
    # Random reads against the model (SNPs only; the order of equal alternatives is not specified so the alt is checked by count)
    bam2, reads = random_bam
    snps = [v for v in bam2.variants("chrA", ref_a, 5, 0.3, 10, 20) if v[1] == "snp"]
    assert [(p, r, n, d) for p, k, r, a, n, d in snps] == model_snps(reads, ref_a, 5, 0.3, 10, 20)
    counts = model_pileup(reads, 0, 0, 2000, 10, 20)
    assert all(counts[p]["ACGT".index(a)] == n for p, k, r, a, n, d in snps)


def model_snps(reads, ref, min_depth, min_fraction, min_baseq=0, min_mapq=0):
    counts = model_pileup(reads, 0, 0, len(ref), min_baseq, min_mapq)
    out = []
    for pos, c in enumerate(counts):
        depth = sum(c[:6])
        if depth < max(min_depth, 1):
            continue
        want = "ACGT".index(ref[pos])
        n = max(c[s] for s in range(4) if s != want)
        if n > 0 and n >= min_fraction * depth:
            out.append((pos, ref[pos], n, depth))
    return out


@pytest.fixture(scope="module")
def snp_bam(tmp_path_factory, ref_a):
    """Ten reads over 100..175; at 120 five of them carry a SNP (50 %), at 125 one (10 %), at 130 three (30 %)."""
    reads = []
    alt = lambda p: next(b for b in "ACGT" if b != ref_a[p])
    for i in range(10):
        subs = {}
        if i < 5:
            subs[120] = alt(120)
        if i == 5:
            subs[125] = alt(125)
        if i in (6, 7, 8):
            subs[130] = alt(130)
        reads.append(build_read(f"s{i}", 100 + i // 2, [(70, "M")], ref_a, subs=subs, qual=[40] * 70, flag=0x10 if i % 2 else 0))
    path = tmp_path_factory.mktemp("bam") / "snp.bam"
    write_bam(path, [("chrA", 2000)], reads)
    return core.Bam(str(path)), reads, alt


# ---- tracks: GFF3, GTF, BED ----------------------------------------------------------------------------------------------


def test_gff3_attributes_and_malformed_lines(tmp_path):
    """GFF3: escaped attributes, names and coordinates; malformed lines are errors, not crashes."""
    # --- gff3 escaped attributes names and coordinates ---
    path = tmp_path / "t.gff3"
    path.write_text("##gff-version 3\n"
                    "chr1\tsrc\tregion\t1\t1000\t.\t+\t.\tID=chr1\n"                              # whole-sequence entries are skipped
                    "chr1\tsrc\tgene\t10\t100\t.\t+\t.\tID=g1;Name=big%2Cgene\n"                   # %2C = comma
                    "chr1\tsrc\tgene\t120\t180\t.\t-\t.\tID=g2;Name=a%3Bb%3Dc\n"                   # ; and = escaped
                    "chr1\tsrc\tgene\t200\t300\t.\t.\t.\tName=100%25 %C3%A9\n"                      # %25 = %, %C3%A9 = e acute
                    "chr1\tsrc\tCDS\t310\t320\t.\t+\t0\tID=c1;gene=araC\n"
                    "chr1\tsrc\texon\t330\t340\t.\t+\t.\tID=e1;gene_name=bad%ZZ;note=x\n"           # a bad escape stays as written
                    "chr1\tsrc\tmRNA\t350\t360\t.\t+\t.\tgene_id \"x1\"; transcript_id \"t1\"; Name \"gtfname\"\n"   # GTF quoting, Name first
                    "chr1\tsrc\ttRNA\t370\t380\t.\t+\t.\tcolor=red\n"                              # no name at all: the type
                    "chr1\tsrc\tgene\t400\t400\t.\t+\t.\tName=onebase\n"                            # one base long
                    "##FASTA\n>chr1\nACGT\n")
    t = core.Track(str(path))
    items, more = t.query("chr1", 0, 1000, 100)
    assert not more
    assert [i[3] for i in items] == ["big,gene", "a;b=c", "100% é", "araC", "bad%ZZ", "gtfname", "tRNA", "onebase"], "in start order"
    by_name = {i[3]: i for i in items}
    assert by_name["big,gene"][:3] == (9, 100, 1), "1-based inclusive 10..100 is 0-based 9..100"
    assert by_name["a;b=c"][:3] == (119, 180, -1)
    assert by_name["100% é"][2] == 0, "'.' strand is unknown"
    assert by_name["araC"][4] == "CDS" and by_name["araC"][:3] == (309, 320, 1), "no Name: the gene attribute names it"
    assert by_name["gtfname"][:3] == (349, 360, 1)
    assert by_name["tRNA"][4] == "tRNA"
    assert by_name["onebase"][:2] == (399, 400)
    assert "chr1" in [s for s, _ in t.sequences()] and t.sequences() == [("chr1", len(items))]
    assert not [i for i in items if i[4] in ("region", "chromosome", "source")]
    # Name beats gene_name, gene, locus_tag, gene_id, ID, transcript_id: the order of the preference list
    order = ["Name", "gene_name", "gene", "locus_tag", "gene_id", "ID", "transcript_id"]
    for k in range(len(order)):
        attrs = ";".join(f"{key}=v_{key}" for key in reversed(order[k:]))
        p = tmp_path / f"o{k}.gff3"
        p.write_text(f"c\ts\tgene\t1\t5\t.\t+\t.\t{attrs}\n")
        assert core.Track(str(p)).query("c", 0, 10)[0][0][3] == f"v_{order[k]}"
    # --- gff malformed lines are errors not crashes ---
    for i, line in enumerate(["chr1\ts\tgene\t10\t5\t.\t+\t.\tName=x",       # end before start
                              "chr1\ts\tgene\t0\t5\t.\t+\t.\tName=x",        # GFF is 1-based: 0 is no position
                              "chr1\ts\tgene\tten\t20\t.\t+\t.\tName=x",     # not a number
                              "chr1\ts\tgene\t1\t20\t.\t+",                  # too few columns
                              "just some text"]):
        p = tmp_path / f"bad{i}.gff3"
        p.write_text("##gff-version 3\n" + line + "\n")
        with pytest.raises(OSError):
            core.Track(str(p))
    with pytest.raises(OSError):
        core.Track(str(tmp_path / "no_such_file.gff3"))


def test_bed_columns_headers_and_errors(tmp_path):
    """BED with 3 to 12 columns; headers, comments, blank lines, zero-length features and errors."""
    # --- bed with 3 to 12 columns ---
    full = ["chr1", "100", "200", "geneX", "500", "-", "110", "190", "255,0,0", "2", "10,20", "0,80"]
    for columns in range(3, 13):
        p = tmp_path / f"t{columns}.bed"
        p.write_text("\t".join(full[:columns]) + "\n" + " ".join(["chr2", "5", "15"] + full[3:columns]) + "\n")   # tabs, then spaces
        t = core.Track(str(p))
        (a,) = t.query("chr1", 0, 1000)[0]
        assert (a[0], a[1], a[4]) == (100, 200, "bed"), f"{columns} columns: BED is 0-based, end-exclusive: used as written"
        assert a[3] == ("geneX" if columns >= 4 else ""), f"{columns} columns: the name is column 4"
        assert a[2] == (-1 if columns >= 6 else 0), f"{columns} columns: the strand is column 6"
        (b,) = t.query("chr2", 0, 1000)[0]
        assert (b[0], b[1], b[2], b[3]) == (5, 15, a[2], a[3]), f"{columns} columns: space separated gives the same"
        assert t.sequences() == [("chr1", 1), ("chr2", 1)], columns
    # --- bed headers comments blank lines zero length and errors ---
    p = tmp_path / "t.bed"
    p.write_text("track name=x\nbrowser position chr1:1-100\n# comment\n\nchr1\t10\t10\tempty\nchr1\t20\t30\tn\t0\t.\nchr1\t40\t50\tp\t0\t+\n")
    t = core.Track(str(p))
    items = t.query("chr1", 0, 100)[0]
    assert [(i[0], i[1], i[2], i[3]) for i in items] == [(10, 10, 0, "empty"), (20, 30, 0, "n"), (40, 50, 1, "p")]
    # a zero-length feature is found in a window that contains its position, not in one that starts at it
    assert [i[3] for i in t.query("chr1", 5, 15)[0]] == ["empty"]
    assert t.query("chr1", 10, 20)[0] == []
    assert [i[3] for i in t.query("chr1", 0, 11)[0]] == ["empty"]
    for bad in ("chr1\tnot-a-number\t10\n", "chr1\t30\t10\n", "chr1\t30\n"):
        p.write_text(bad)
        with pytest.raises(OSError):
            core.Track(str(p))


@pytest.fixture(scope="module")
def track_of_many(tmp_path_factory):
    rng = random.Random(71)
    path = tmp_path_factory.mktemp("track") / "many.bed"
    rows = []
    for i in range(400):
        s = rng.randrange(0, 9_000)
        rows.append((s, s + rng.randint(1, 600), f"feat{i}"))
    path.write_text("".join(f"chr1\t{s}\t{e}\t{n}\n" for s, e, n in rows))
    return core.Track(str(path)), rows


def test_track_query_density_and_search(track_of_many, tmp_path):
    """Tracks: query equals a brute-force overlap and the limit flag; density bins equal the overlap count; search is
    case-insensitive and limited."""
    # --- track query equals a brute force overlap and the limit flag ---
    t, rows = track_of_many
    for start, end in ((0, 10_000), (1000, 1100), (4321, 4322), (8900, 9700), (20_000, 30_000)):
        expect = sorted(n for s, e, n in rows if s < end and e > start)
        items, more = t.query("chr1", start, end, 100_000)
        assert sorted(i[3] for i in items) == expect and not more
        assert [i[0] for i in items] == sorted(i[0] for i in items), "start order"
    everything = t.query("chr1", 0, 10_000, 100_000)[0]
    assert len(everything) == 400
    # limit: exactly `limit` items, and the flag says whether there were more
    items, more = t.query("chr1", 0, 10_000, 399)
    assert len(items) == 399 and more is True
    items, more = t.query("chr1", 0, 10_000, 400)
    assert len(items) == 400 and more is False
    assert t.query("chr1", 0, 10_000, 1)[1] is True
    assert t.query("chr2", 0, 10_000) == ([], False), "unknown sequence: nothing, no error"
    # default limit is 5000: the 400 come in one piece
    assert len(t.query("chr1", 0, 10_000)[0]) == 400
    # --- track density bins equal the overlap count and sum to the item count ---
    t, rows = track_of_many
    # Each bin counts the items overlapping it, so an item over several bins is counted in each of them
    for start, end, bins in ((0, 10_000, 100), (0, 10_000, 7), (2000, 6000, 40)):
        width = (end - start) / bins
        expect = [sum(1 for s, e, _ in rows if s < start + (b + 1) * width and e > start + b * width) for b in range(bins)]
        got = t.density("chr1", start, end, bins)
        assert got == expect, (start, end, bins)
        assert len(got) == bins
    # Items that each sit inside one bin: the bins add up to the number of items
    p = tmp_path / "tidy.bed"
    p.write_text("".join(f"c\t{s}\t{s + 5}\tf{s}\n" for s in range(2, 1000, 10)))      # 100 items of 5 bases in bins of 10
    tidy = core.Track(str(p))
    assert sum(tidy.density("c", 0, 1000, 100)) == 100
    assert tidy.density("c", 0, 1000, 100) == [1] * 100
    # No items in a window, no bins, an empty window
    assert tidy.density("c", 2000, 3000, 4) == [0, 0, 0, 0] and tidy.density("c", 0, 100, 0) == []
    assert tidy.density("c", 50, 50, 3) == [0, 0, 0]
    assert tidy.density("nope", 0, 100, 2) == [0, 0]
    # --- track search is case insensitive and limited ---
    p = tmp_path / "s.bed"
    p.write_text("chr1\t0\t10\tLacZ\nchr1\t20\t30\tlacY\nchr2\t0\t5\tÉCOLI-gene\nchr2\t10\t15\tother\n" + "".join(f"chr3\t{i}\t{i + 1}\tx{i}lac\n" for i in range(10)),
                 encoding="utf-8")
    t = core.Track(str(p))
    for needle in ("lac", "LAC", "LaC"):
        got = t.search(needle, 100)
        assert sorted(g[4] for g in got) == sorted(["LacZ", "lacY"] + [f"x{i}lac" for i in range(10)]), needle
    assert [g[4] for g in t.search("écoli", 10)] == ["ÉCOLI-gene"], "unicode case folding too"
    assert t.search("ecoli", 10) == [], "no accent folding"
    assert len(t.search("lac", 3)) == 3 and len(t.search("lac")) == 12
    row = t.search("lacz")[0]
    assert row == ("chr1", 0, 10, 0, "LacZ", "bed")
    assert t.search("zzz") == []


# ---- feature locations ---------------------------------------------------------------------------------------------------

# GenBank location text -> how it is written back (and the spans, 0-based end-exclusive, and the strand)
LOCATION_FORMS = [
    ("10..20", "10..20", [(9, 20)], 1),
    ("5", "5", [(4, 5)], 1),
    ("5..5", "5", [(4, 5)], 1),
    ("1..200", "1..200", [(0, 200)], 1),
    ("complement(10..20)", "complement(10..20)", [(9, 20)], -1),
    ("join(1..5,10..20)", "join(1..5,10..20)", [(0, 5), (9, 20)], 1),
    ("order(1..5,10..20)", "order(1..5,10..20)", [(0, 5), (9, 20)], 1),
    ("complement(join(1..5,10..20))", "complement(join(1..5,10..20))", [(0, 5), (9, 20)], -1),
    ("complement(order(1..5,10..20))", "complement(order(1..5,10..20))", [(0, 5), (9, 20)], -1),
    ("join(complement(1..5),complement(10..20))", "join(complement(1..5),complement(10..20))", [(0, 5), (9, 20)], -1),
    ("<10..20", "<10..20", [(9, 20)], 1),
    ("10..>20", "10..>20", [(9, 20)], 1),
    ("<1..>20", "<1..>20", [(0, 20)], 1),
    ("complement(<10..20)", "complement(<10..20)", [(9, 20)], -1),
    ("10^11", "10^11", [(10, 10)], 1),                                                     # the gap is where base 11 starts
    ("200^1", "200^1", [(0, 0)], 1),                                                       # across the origin of the 200 bp circle
    ("one-of(5,8)", "one-of(5,8)", [(4, 5), (7, 8)], 1),
    ("bond(1..5,10..20)", "bond(1..5,10..20)", [(0, 5), (9, 20)], 1),
    ("join(1..5,J00194.1:100..202)", "join(1..5,J00194.1:100..202)", [(0, 5)], 1),          # a part on another record
    ("join(1..5,gap(10),20..30)", "join(1..5,gap(10),20..30)", [(0, 5), (19, 30)], 1),      # a gap of 10 bases
]


def test_locations_round_trip_and_refusals():
    """Feature locations: every GenBank form parses and is written back as a fixed point; bad ones are refused."""
    # --- location parse and format round trip ---
    for text, written, spans, strand in LOCATION_FORMS:
        doc = Document.new_sequence("l", "ACGT" * 50, True)
        doc.add_features([("misc_feature", text, [("label", "x")])])
        f = doc.features()[-1]
        assert f.location == written and f.spans == spans and f.strand == strand, text
        # the written form is a fixed point: it parses to itself
        doc.add_features([("misc_feature", f.location, [("label", "y")])])
        g = doc.features()[-1]
        assert (g.location, g.spans, g.strand) == (f.location, f.spans, f.strand), text
        # set_feature takes the same text
        index = doc.set_feature(None, "gene", text, [])
        assert doc.features()[index].location == written, text
    # --- locations outside the sequence or not locations are refused ---
    doc = Document.new_sequence("l", "ACGT" * 50, True)
    for bad in ("10..300", "0..5", "20..10", "bad", ""):
        with pytest.raises((OSError, ValueError)):
            doc.add_features([("misc_feature", bad, [])])
    assert len(doc.features()) == 0, "a refused feature leaves nothing behind"
    with pytest.raises(OSError):
        doc.add_features([("two words", "1..5", [])])


def test_locations_in_genbank_files(tmp_path):
    """Locations in GenBank files: a fuzzy end range is kept word for word; the forms round trip through save."""
    # --- a location with a fuzzy end range is kept word for word ---
    # `10..(20.30)` is not understood as a range, but it must never be cut to its first base: a file keeps it as written
    # (the feature has no extent, so it is not drawn), and typing one in is refused with a message.
    gb = ("LOCUS       t                        200 bp    DNA     circular UNK 01-JAN-1980\n"
          "FEATURES             Location/Qualifiers\n"
          "     misc_feature    10..(20.30)\n                     /label=\"f\"\n"
          "     misc_feature    join(1..5,(10.20)..(30.40))\n                     /label=\"g\"\n"
          "ORIGIN\n        1 " + "acgt" * 50 + "\n//\n")
    src = tmp_path / "fz.gb"
    src.write_text(gb)
    doc = Document(str(src))
    assert [f.location for f in doc.features()] == ["10..(20.30)", "join(1..5,(10.20)..(30.40))"]
    out = tmp_path / "fz2.gb"
    doc.write_genbank(str(out))
    assert [f.location for f in Document(str(out)).features()] == ["10..(20.30)", "join(1..5,(10.20)..(30.40))"]
    with pytest.raises(OSError, match="outside"):
        doc.add_features([("misc_feature", "10..(20.30)", [])])
    # --- locations from a genbank file round trip through save ---
    text = "".join(f"     misc_feature    {loc}\n                     /label=\"f{i}\"\n" for i, (loc, *_rest) in enumerate(LOCATION_FORMS)
                    if "gap" not in loc)
    gb = ("LOCUS       t                        200 bp    DNA     circular UNK 01-JAN-1980\nDEFINITION  .\nFEATURES             Location/Qualifiers\n"
          + text + "ORIGIN\n        1 " + "acgt" * 50 + "\n//\n")
    src = tmp_path / "t.gb"
    src.write_text(gb)
    doc = Document(str(src))
    before = [(f.location, f.spans, f.strand) for f in doc.features()]
    assert len(before) == len(LOCATION_FORMS) - 1
    out = tmp_path / "out.gb"
    doc.write_genbank(str(out))
    again = Document(str(out))
    assert [(f.location, f.spans, f.strand) for f in again.features()] == before
    # A location that spans two sides of the origin on a circle: join(190..200,1..10), with its strand
    doc.add_features([("misc_feature", "complement(join(190..200,1..10))", [])])
    f = doc.features()[-1]
    assert f.spans == [(189, 200), (0, 10)] and f.strand == -1


# ---- genome FASTA --------------------------------------------------------------------------------------------------------


@pytest.fixture(scope="module")
def genome(tmp_path_factory):
    """Two chromosomes wrapped at 60 columns: one of 250 bases, one that is exactly two full lines (120)."""
    rng = random.Random(81)
    seqs = {"chr1": rand_dna(rng, 250), "chr2": rand_dna(rng, 120), "chr3": rand_dna(rng, 7)}
    path = tmp_path_factory.mktemp("genome") / "g.fa"
    path.write_text("".join(f">{n} some description\n" + "\n".join(s[i:i + 60] for i in range(0, len(s), 60)) + "\n" for n, s in seqs.items()))
    return core.Fasta(str(path)), seqs, path


def test_genome_fasta_fetch_and_motif_search(genome):
    """Indexed genome FASTA: sequences and fetch (first / last base, line breaks, the end), random ranges equal slicing, a
    search across a line break and the IUPAC search on the file agrees with find_motif."""
    # --- genome fasta sequences and fetch first last base line breaks and the end ---
    fa, seqs, path = genome
    assert fa.sequences() == [("chr1", 250), ("chr2", 120), ("chr3", 7)]
    assert path.with_name(path.name + ".fai").exists()
    for name, s in seqs.items():
        n = len(s)
        assert fa.fetch(name, 0, 1) == s[0], "first base"
        assert fa.fetch(name, n - 1, n) == s[-1], "last base"
        assert fa.fetch(name, 0, n) == s, "everything"
        assert fa.fetch(name, 0, n + 100) == s, "an end past the sequence is cut to it"
        assert fa.fetch(name, n - 3, n + 100) == s[-3:]
    s = seqs["chr1"]
    # across line breaks: 59|60, 119|120, 179|180, a range of several lines
    for a, b in ((58, 62), (59, 61), (60, 61), (59, 60), (118, 122), (179, 181), (10, 245), (0, 60), (0, 61), (60, 120)):
        assert fa.fetch("chr1", a, b) == s[a:b], (a, b)
    # the last base of a sequence whose length is a multiple of the line length (no short last line)
    assert fa.fetch("chr2", 119, 120) == seqs["chr2"][119] and fa.fetch("chr2", 60, 120) == seqs["chr2"][60:]
    # empty or inverted ranges: nothing
    assert fa.fetch("chr1", 5, 5) == "" and fa.fetch("chr1", 9, 5) == "" and fa.fetch("chr1", 250, 260) == ""
    with pytest.raises(OSError):
        fa.fetch("chrZ", 0, 5)
    # --- genome fasta random ranges equal slicing ---
    fa, seqs, _ = genome
    rng = random.Random(82)
    for _ in range(300):
        name = rng.choice(list(seqs))
        a = rng.randrange(0, len(seqs[name]))
        b = rng.randrange(a, len(seqs[name]) + 5)
        assert fa.fetch(name, a, b) == seqs[name][a:b]
    # --- genome fasta search across a line break and motif search agrees with find motif ---
    fa, seqs, _ = genome
    s = seqs["chr1"]
    across = s[57:63]                                # six bases straddling the first line break
    assert any(p == 57 for p, st in fa.search("chr1", across) if st == 1)
    # the IUPAC search on the file equals the one on a string
    for motif, mism in (("GARTTC", 0), ("ACGT", 1), ("NNGG", 0)):
        on_file = fa.search_motif("chr1", motif, max_mismatches=mism)
        on_string = core.find_motif(s, motif, max_mismatches=mism)
        assert on_file == on_string


def test_genome_fasta_search_both_strands_and_overlapping_hits(tmp_path):
    rng = random.Random(83)
    filler = lambda n: rand_dna(rng, n, "AC")     # no G or T, so no accidental hit for a motif that has both
    motif = "GATTACA"
    seq = filler(30) + motif + filler(40) + rc(motif) + filler(25) + "AAAAAAAAA" + filler(5) + "GAATTC" + filler(10)
    path = tmp_path / "s.fa"
    path.write_text(">s\n" + "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60)) + "\n")
    fa = core.Fasta(str(path))
    assert fa.search("s", motif) == [(30, 1), (77, -1)]
    assert fa.search("s", rc(motif)) == [(30, -1), (77, 1)], "the reverse complement of the motif finds the same places, strands swapped"
    assert fa.search("s", motif.lower()) == fa.search("s", motif), "case does not matter"
    # overlapping hits: AAAA in nine A is found at six positions, one after the other (and as TTTT on the other strand)
    a_start = seq.index("AAAAAAAAA")
    assert [p for p, st in fa.search("s", "AAAA") if a_start <= p <= a_start + 5] == list(range(a_start, a_start + 6))
    # a palindrome is reported once, on the top strand
    e = seq.index("GAATTC")
    assert fa.search("s", "GAATTC") == [(e, 1)]
    # against a brute-force scan, for several motifs, including ones across a line break (60)
    for m in ("AC", "ACA", "CA" * 5, "GAT", "TACAC"):
        expect = []
        for i in range(len(seq) - len(m) + 1):
            w = seq[i:i + len(m)]
            if w == m:
                expect.append((i, 1))
            elif w == rc(m):
                expect.append((i, -1))
        assert fa.search("s", m) == expect, m
    assert fa.search("s", "") == [] and fa.search("s", "G" * 500) == []


