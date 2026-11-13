"""Restriction enzymes (where each strand is cut, the common list, the filter, the cut tooltips) and CRISPR guides
(candidates, scores, off-target counts, and the way a saved guide site is stored in a GenBank file)."""
import os
import random
import re
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ["OPENDNA_DATA_DIR"] = tempfile.mkdtemp(prefix="opendna_test_")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "app"))

import pytest
from Bio import SeqIO

import opendna_core as core
from opendna_core import Document

import cas9
import cut_sites as cs
import enzymes
import exports
from editor import Editor, reverse_complement as rc
from enzymes import Enzyme

pytestmark = pytest.mark.filterwarnings("ignore::Bio.BiopythonParserWarning")

ENZ = {e.name: e for e in enzymes.COMMON}
SP = cs.NUCLEASES["SpCas9 (NGG, 20 nt)"]
SA = cas9.NUCLEASES["SaCas9 (NNGRRT, 21 nt)"]
A12 = cs.NUCLEASES["Cas12a (TTTV, 23 nt)"]
GUIDE = "GACTCGTACGTAGCATCGAT"   # 20 nt, 50 % GC, no GG or CC in it: no PAM of its own on either strand
GUIDE23 = GUIDE + "GCA"         # 23 nt for Cas12a
FILL = "ATATTATAATTATATATTAA"
rng = random.Random(7)


def at_dna(n):
    """Filler of A and T only: no G or C, so no PAM on either strand (the sites are the ones a test puts in)."""
    return "".join(rng.choice("AT") for _ in range(n))


def plain(html):
    """A tooltip as the reader sees it."""
    return re.sub(r"<[^>]+>", "", html.replace("<br>", "\n"))


def tip(name, seq, circular=False, which=0):
    top, bottom = enzymes.find_cuts(seq, circular, [ENZ[name]])[name][which]
    return plain(cs.enzyme_tip(ENZ[name], seq, circular, top, bottom))


# ---- restriction enzymes -------------------------------------------------------------------------------------------


def test_cut_finder_the_common_list_and_the_filter():
    """Degenerate sites, enzymes that cut on both sides, both strands and the origin of a circle; the common list cuts
    where it always did and describes itself the way suppliers write it; and the filter keeps an enzyme only when every
    criterion holds."""
    # Enzymes written the way a longer list would hold them: site (IUPAC) and the cuts as offsets from the site's start.
    acci = Enzyme("AccI", "GTMKAC", ((2, 4),), False)               # GT^MKAC, M = A/C, K = G/T
    bgli = Enzyme("BglI", "GCCNNNNNGGC", ((7, 4),), False)          # GCCNNNN^NGGC: a 3' overhang
    bcgi = Enzyme("BcgI", "CGANNNNNNTGC", ((-10, -12), (24, 22)), False)   # cuts on both sides of its site
    bsai = Enzyme("BsaI", "GGTCTC", ((7, 11),), False)
    ecori = Enzyme("EcoRI", "GAATTC", ((1, 5),), False)
    flank = "AT" * 20   # 40 bases that hold no site of these enzymes on either strand

    # --- AccI GT^MKAC: M is A or C, K is G or T.
    s = "AAA" + "GTAGAC" + "TTT" + "GTCTAC" + "AAA" + "GTGGAC"
    assert enzymes.find_cuts(s, False, [acci])["AccI"] == [(5, 7), (14, 16)], "GTAGAC and GTCTAC match; GTGGAC does not"
    # --- BglI leaves a 3' overhang: the top cut is after the bottom one on the top strand's numbering.
    assert enzymes.find_cuts("TT" + "GCCAAGTTGGC" + "TT", False, [bgli])["BglI"] == [(9, 6)]
    # --- BcgI cuts on both sides of CGANNNNNNTGC: two cut pairs per site, 2 nt 3' overhangs.
    assert enzymes.find_cuts(flank + "CGAGGGGGGTGC" + flank, False, [bcgi])["BcgI"] == [(30, 28), (64, 62)]
    # --- a site that is not its own reverse complement is found on the other strand, with the cuts mirrored.
    fwd = enzymes.find_cuts("AA" + "GGTCTC" + "A" * 20, False, [bsai])["BsaI"]
    rev = enzymes.find_cuts("A" * 20 + "GAGACC" + "AA", False, [bsai])["BsaI"]
    assert fwd == [(9, 13)] and rev == [(15, 19)] and rev[0][1] - rev[0][0] == 4
    # --- across the origin of a circle.
    assert enzymes.find_cuts("ATTCAAAAAAGA", True, [ecori])["EcoRI"] == [(11, 3)]

    # --- the common list ---
    assert ENZ["EcoRI"].cuts == ((1, 5),) and ENZ["KpnI"].cuts == ((5, 1),) and ENZ["EcoRV"].blunt
    assert ENZ["BsaI"].cuts == ((7, 11),) and ENZ["SapI"].cuts == ((8, 11),)
    assert enzymes.find_cuts("TTGAATTCAA", False)["EcoRI"] == [(3, 7)]
    assert enzymes.find_cuts("AAGATATCAA", False)["EcoRV"] == [(5, 5)], "a blunt cutter cuts both strands in the same place"
    assert enzymes.find_sites("TTGAATTCAA" + "GAGACCTTTTT" + "TCAA", False) == {"EcoRI": [3], "BsaI": [5]}
    assert enzymes.describe(ENZ["EcoRI"]) == "GAATTC: G^AATTC, 5′ overhang 4 nt"
    assert enzymes.describe(ENZ["KpnI"]).endswith("3′ overhang 4 nt") and "blunt" in enzymes.describe(ENZ["EcoRV"])
    assert "more place" in enzymes.describe(bcgi)
    # A Type IIS enzyme cuts outside its site: written the way suppliers do, GGTCTC(1/5), so it is clear that the
    # overhang's sequence is whatever the user puts between the site and the cut (the point of Golden Gate).
    assert enzymes.describe(ENZ["BsaI"]) == "GGTCTC: GGTCTC(1/5), 5′ overhang 4 nt"
    assert enzymes.describe(ENZ["BbsI"]).startswith("GAAGAC: GAAGAC(2/6),")
    assert enzymes.describe(ENZ["SapI"]) == "GCTCTTC: GCTCTTC(1/4), 5′ overhang 3 nt"

    # --- the filter: how often it cuts, site size, palindromic, Type IIS, kind of end; all of them must hold ---
    from enzymes import EnzymeFilter, wanted
    anything = dict(min_sites=1, max_sites=None)
    assert EnzymeFilter() == EnzymeFilter(1, 1, 6), "the default is the unique cutters of 6+ bp"
    assert not wanted(EnzymeFilter(), ENZ["AluI"], 1) and wanted(EnzymeFilter(), ENZ["EcoRI"], 1)
    assert wanted(EnzymeFilter(palindromic=True, **anything), ENZ["EcoRI"], 1)
    assert not wanted(EnzymeFilter(palindromic=True, **anything), ENZ["BsaI"], 1), "GGTCTC is not its own reverse complement"
    assert wanted(EnzymeFilter(type_iis=True, **anything), ENZ["BsaI"], 1)
    assert not wanted(EnzymeFilter(type_iis=True, **anything), ENZ["EcoRI"], 1)
    assert wanted(EnzymeFilter(type_iis=False, **anything), ENZ["EcoRI"], 1)
    assert wanted(EnzymeFilter(ends="5", **anything), ENZ["EcoRI"], 1) and wanted(EnzymeFilter(ends="3", **anything), ENZ["KpnI"], 1)
    assert wanted(EnzymeFilter(ends="blunt", **anything), ENZ["EcoRV"], 1)
    assert not wanted(EnzymeFilter(ends="blunt", **anything), ENZ["EcoRI"], 1)
    assert wanted(EnzymeFilter(min_size=8, **anything), ENZ["NotI"], 1) and not wanted(EnzymeFilter(min_size=8, **anything), ENZ["EcoRI"], 1)
    both = EnzymeFilter(palindromic=True, ends="5", min_size=6, **anything)
    assert wanted(both, ENZ["EcoRI"], 1) and not wanted(both, ENZ["KpnI"], 1) and not wanted(both, ENZ["EcoRV"], 1), "all criteria must hold"
    # How often it cuts is a criterion too: 2 to 3 sites of at most 4 bp.
    narrow = EnzymeFilter(2, 3, 1, 4)
    assert wanted(narrow, ENZ["AluI"], 3) and not wanted(narrow, ENZ["AluI"], 1) and not wanted(narrow, ENZ["EcoRI"], 2)


def test_what_a_cut_says_about_itself():
    """The tooltip of a 5' overhang, a 3' overhang, a blunt and a Type IIS enzyme, of a site on the bottom strand and of
    one across the origin of a circle; every common enzyme has one at each of its cuts."""
    # --- 5' overhang ---
    t = tip("EcoRI", "TTGAATTCAA")
    assert "EcoRI" in t and "G^AATTC" in t and "bases 3–8" in t and "Type IIS" not in t
    assert "Cuts after base 3 on the top strand and after base 7 on the bottom strand." in t
    assert "5′ overhang of 4 nt: AATT" in t
    assert "top strand)" in t, "a palindromic enzyme has one site, listed as the top strand"

    # --- 3' overhang and blunt ---
    t = tip("KpnI", "TTGGTACCAA")
    assert "GGTAC^C" in t and "Cuts after base 7 on the top strand and after base 3 on the bottom strand." in t
    assert "3′ overhang of 4 nt: GTAC" in t
    t = tip("EcoRV", "TTGATATCAA")
    assert "GAT^ATC" in t and "Cuts after base 5 on the top strand and after base 5 on the bottom strand." in t
    assert "Blunt ends" in t and "overhang of" not in t

    # --- a Type IIS enzyme says how far from its site it cuts ---
    t = tip("BsaI", "TTGGTCTCAAGGAATTTT")   # GGTCTC at bases 3-8; cuts 1 nt after on the top strand, 5 after on the bottom
    assert "GGTCTC(1/5)" in t and "5′ overhang of 4 nt: AGGA" in t
    assert "Cuts after base 9 on the top strand and after base 13 on the bottom strand." in t
    assert "Type IIS: it cuts outside its site, 1 nt after the site on the top strand and 5 nt after the site on the bottom strand." in t

    # --- a site on the bottom strand ---
    t = tip("BsaI", "TTTTAAAAGGAAGAGACCTT")   # GAGACC = GGTCTC on the bottom strand, at bases 13-18; the cuts are before it
    assert "bases 13–18, on the bottom strand, where it reads GGTCTC; GAGACC on the top strand" in t
    assert "Type IIS" in t and "before the site" in t
    top, bottom = enzymes.find_cuts("TTTTAAAAGGAAGAGACCTT", False, [ENZ["BsaI"]])["BsaI"][0]
    assert (top, bottom) == (7, 11) and f"after base {top} on the top strand and after base {bottom} on the bottom" in t

    # --- a site across the origin of a circle ---
    t = tip("EcoRI", "ATTCAAAAAAGA", circular=True)   # GAATTC = GA at the end + ATTC at the start
    assert "bases 11–4" in t and "Cuts after base 11 on the top strand and after base 3 on the bottom strand." in t
    assert "5′ overhang of 4 nt: AATT" in t

    # --- every common enzyme has a tooltip at each cut ---
    seq = "".join(e.site.replace("N", "A") for e in enzymes.COMMON) * 2
    for name, cuts in enzymes.find_cuts(seq, True, enzymes.COMMON).items():
        for top, bottom in cuts:
            t = plain(cs.enzyme_tip(ENZ[name], seq, True, top, bottom))
            assert name in t and "Cuts after base" in t and "Recognition site" in t


# ---- Cas9 / Cas12a guides ---------------------------------------------------------------------------------------------


def test_guide_candidates_cuts_and_scores():
    """Guides on both strands with their cut, in a stretch and across the origin of a circle; SaCas9's longer PAM;
    Cas12a's 5' PAM and staggered cut; the score is a checklist and the best guides come first."""
    # --- SpCas9 on both strands ---
    top = at_dna(40) + GUIDE + "TGG" + at_dna(40)
    (g,) = cas9.candidates(top, False, SP)
    assert (g.start, g.end, g.strand, g.guide, g.pam) == (40, 63, 1, GUIDE, "TGG")
    assert g.cut == 40 + 20 - 3, "Cas9 cuts 3 bp from the PAM, inside the protospacer"
    bottom = at_dna(40) + rc(GUIDE + "AGG") + at_dna(40)
    (r,) = cas9.candidates(bottom, False, SP)
    assert (r.start, r.end, r.strand, r.guide, r.pam) == (40, 63, -1, GUIDE, "AGG")
    assert r.cut == 40 + 3 + 3, "3 bp inside the protospacer, counted from the PAM on the left"
    assert bottom[r.start:r.end] == rc(r.guide + r.pam)

    # --- an offset stretch, a circle, and SaCas9 ---
    seq = at_dna(60) + GUIDE + "CGG" + at_dna(30)
    (g,) = cas9.candidates(seq[50:], False, SP, offset=50)
    assert g.start == 60 and g.guide == GUIDE, "a stretch is read at its place on the whole sequence"
    rotated = (GUIDE + "TGG")[12:] + at_dna(80) + (GUIDE + "TGG")[:12]
    (c,) = cas9.candidates(rotated, True, SP, wrap=True)
    assert c.guide == GUIDE and c.start == len(rotated) - 12 and c.end > len(rotated)
    assert not cas9.candidates(rotated, False, SP), "a line has no site across its ends"
    (s,) = [c for c in cas9.candidates(at_dna(30) + GUIDE + "A" + "ATGAAT" + at_dna(30), False, SA) if c.strand > 0]
    assert len(s.guide) == 21 and s.pam == "ATGAAT"
    assert not [c for c in cas9.candidates(at_dna(30) + GUIDE + "A" + "ATGCCT" + at_dna(30), False, SA) if c.strand > 0], "ATGCCT is no NNGRRT"

    # --- Cas12a: the PAM is on the 5' side and the two strands are cut in different places ---
    pad = "CA" * 25             # no T run, no AAA: no PAM of its own on either strand
    top = pad + "TTTA" + GUIDE23 + pad
    (g12,) = cas9.candidates(top, False, A12)
    assert (g12.start, g12.end, g12.strand, g12.guide, g12.pam) == (50, 77, 1, GUIDE23, "TTTA")
    assert (g12.cut, g12.bottom) == (50 + 4 + 18, 50 + 4 + 23), "the PAM's strand is cut after 18 bases, the other after 23"
    bottom12 = pad + rc("TTTC" + GUIDE23) + pad
    (r12,) = cas9.candidates(bottom12, False, A12)
    assert (r12.start, r12.end, r12.strand, r12.guide, r12.pam) == (50, 77, -1, GUIDE23, "TTTC")
    assert (r12.cut, r12.bottom) == (50, 55) and bottom12[r12.start:r12.end] == rc(r12.pam + r12.guide)
    assert [x.start for x in cas9.candidates(pad + "TTTT" + GUIDE23 + pad, False, A12)] == [51], "TTTT is no PAM; the TTTG after its first T is"
    assert "5′ overhang" in cas9.cut_text(A12, g12, len(top))
    assert "3 bp" in cas9.cut_text(SP, cas9.candidates(GUIDE + "TGG" + pad, False, SP)[0], 99)

    # --- the score is a plain checklist, and the best guides come first ---
    assert cas9.score_guide(GUIDE) == (100, [])
    s, notes = cas9.score_guide("ACGTACGTTTTTACGTACGA")
    assert s < 70 and any("TTTT" in n for n in notes)
    s, notes = cas9.score_guide("GGGGGGGGGGGGGGGGGGGG")
    assert s <= 45 and any("GC" in n for n in notes) and any("run" in n for n in notes)
    s, notes = cas9.score_guide("ATATATATATTATAATATAT")
    assert s <= 60 and any("GC" in n for n in notes)
    assert cas9.score_guide("ACGATCGTAGCTAGCTAGCA")[0] == 100
    found = cas9.candidates(at_dna(30) + GUIDE + "TGG" + at_dna(20) + "ACGTACGTTTTTACGTACGA" + "TGG" + at_dna(30), False, SP)
    assert [g.score for g in found] == sorted((g.score for g in found), reverse=True) and found[0].guide == GUIDE


def test_other_guide_sites_here_and_in_a_genome():
    """The other sites a guide has here and in a reference genome, counted by mismatches, with the verdict they give; a
    mismatch next to the PAM counts for more; the Rust search finds exactly what a plain Python one does."""
    bend = lambda guide, positions: "".join("ACGT"["ACGT".index(c) - 1] if i in positions else c for i, c in enumerate(guide))
    near = bend(GUIDE, [2])               # one mismatch, far from the PAM
    far = bend(GUIDE, [0, 1, 5, 19])      # four mismatches, one next to the PAM
    seq = at_dna(30) + GUIDE + "TGG" + at_dna(30) + near + "AGG" + at_dna(30) + far + "CGG" + at_dna(30)
    (g,) = [x for x in cas9.candidates(seq, False, SP) if x.guide == GUIDE]
    (spec,) = cas9.scan(seq, False, "plasmid", [g], SP, 3, None)
    assert spec.here == [0, 1, 0, 0] and spec.genome is None, "one near site; the four-mismatch one is beyond 3"
    assert cas9.verdict(spec, g)[0] == "warn"
    # A genome holding the guide twice (a perfect copy and one with two mismatches).
    two = bend(GUIDE, [3, 9])
    ref = Document.new_sequence("ref", at_dna(50) + GUIDE + "AGG" + at_dna(40) + two + "TGG" + at_dna(50), True)
    (spec,) = cas9.scan(seq, False, "plasmid", [g], SP, 3, ref)
    assert spec.genome == [1, 0, 1, 0] and spec.genome_name == "ref" and cas9.verdict(spec, g)[0] == "warn"
    # Two perfect copies in the genome, or another perfect site here, is bad.
    ref2 = Document.new_sequence("ref2", GUIDE + "AGG" + at_dna(30) + GUIDE + "CGG" + at_dna(30), False)
    (spec,) = cas9.scan(seq, False, "plasmid", [g], SP, 3, ref2)
    assert spec.genome[0] == 2 and cas9.verdict(spec, g)[0] == "bad"
    twice = at_dna(20) + GUIDE + "TGG" + at_dna(20) + GUIDE + "AGG" + at_dna(20)
    (g1, g2) = sorted(cas9.candidates(twice, False, SP), key=lambda x: x.start)
    (s1, s2) = cas9.scan(twice, False, "p", [g1, g2], SP, 3, None)
    assert s1.here[0] == 1 and s2.here[0] == 1 and cas9.verdict(s1, g1)[0] == "bad"
    # The genome is the sequence itself: it is not scanned twice.
    (spec,) = cas9.scan(seq, False, "plasmid", [g], SP, 3, Document.new_sequence("plasmid", seq, False))
    assert spec.genome is None and cas9.verdict(None, g)[0] == "info"
    # A Cas12a site is counted with its PAM before the protospacer.
    pad = "CA" * 25
    text = pad + "TTTA" + GUIDE23 + pad + "TTTC" + GUIDE23 + pad
    guides = cas9.candidates(text, False, A12)
    assert len(guides) == 2 and cas9.scan(text, False, "t", guides[:1], A12, 1)[0].here == [1, 0]

    # --- seed mismatches (the bases next to the PAM) ---
    assert cas9.seed_mismatches(1, 20) == 0, "base 0 is at the far end"
    assert cas9.seed_mismatches(1 << 19, 20) == 1, "base 19 touches the PAM"
    assert cas9.seed_mismatches((1 << 19) | (1 << 7) | 1, 20) == 1 and cas9.seed_mismatches((1 << 8) | (1 << 12), 20) == 2

    # --- the Rust search matches a plain Python one ---
    seq = "".join(rng.choice("ACGT") for _ in range(4000))
    guides = [seq[i:i + 20] for i in (100, 1000, 2500)]
    for guide, found in zip(guides, core.find_guide_sites(seq, False, guides, "NGG", 2)):
        expect = set()
        for strand, s in ((1, seq), (-1, rc(seq))):
            for i in range(len(s) - 22):
                if s[i + 21:i + 23] == "GG" and sum(a != b for a, b in zip(s[i:i + 20], guide)) <= 2:
                    expect.add((i if strand > 0 else len(seq) - i - 23, strand))
        assert {(h[0], h[2]) for h in found} == expect


def guide_doc(strand, nuclease=SP, circular=False, text=None):
    """A document with one guide site, as the Spacer Finder saves it."""
    site = text if text is not None else GUIDE + "TGG"
    if strand < 0:
        site = rc(site)
    seq = FILL * 3 + site + FILL * 3
    lead = len(FILL) * 3
    doc = Document.new_sequence("g", seq, circular)
    loc = f"{lead + 1}..{lead + len(site)}"
    note = f"{nuclease.name.split(' (')[0]} guide {site}; cuts 3 bp from the PAM"   # the readable note the Spacer Finder writes
    doc.add_features([("misc_feature", f"complement({loc})" if strand < 0 else loc,
                       [("label", "gRNA GACTCGTA"), ("note", note), (cs.CAS9_KEY, cs.cas9_qualifier(nuclease))])])
    return doc, lead, len(seq)


def test_a_saved_guide_site_is_read_back_and_round_trips_in_the_file():
    """What the feature of a saved guide says (protospacer, PAM, cut, on both strands, Cas12a's two cuts, across the
    origin), the qualifier value read back, the file it makes (byte for byte, a plain misc_feature to other readers, the
    old key of the program's former name) and the exports that leave the program's own qualifier out."""
    # --- Cas9 on both strands: protospacer, PAM and the blunt cut ---
    for strand in (1, -1):
        doc, lead, n = guide_doc(strand)
        (f,) = [f for f in doc.features() if f.kind == "misc_feature"]
        layout = cs.feature_layout(f, n)
        t = plain(cs.cas9_tip(layout, doc.sequence, n, "gRNA GACTCGTA"))
        assert f"5′-{GUIDE}-3′" in t and "PAM TGG (NGG)" in t and "SpCas9" in t, f"strand {strand}"
        assert "blunt cut" in t and "3 bp from the PAM" in t, f"strand {strand}"
        cut = lead + 20 - 3 if strand > 0 else lead + 3 + 3   # the gap, as the Spacer Finder computes it
        assert f"after base {cut} on both strands" in t and cs.cut_places(SP, lead, strand) == (cut, cut)
        assert layout.gaps == [cut]
        if strand > 0:   # PAM on the right for the top strand, on the left for the bottom
            assert layout.proto == [(lead, lead + 20)] and layout.pam == [(lead + 20, lead + 23)]
        else:
            assert layout.pam == [(lead, lead + 3)] and layout.proto == [(lead + 3, lead + 23)]

    # --- Cas12a has two cuts ---
    doc, lead, n = guide_doc(1, A12, text="TTTA" + "GACTCGTACGTAGCATCGATAGC")
    (f,) = [f for f in doc.features() if f.kind == "misc_feature"]
    t = plain(cs.cas9_tip(cs.feature_layout(f, n), doc.sequence, n, "g"))
    assert "Cas12a" in t and "PAM TTTA (TTTV)" in t and "5 nt 5′ overhang" in t
    assert f"after base {lead + 4 + 18} on the top strand and after base {lead + 4 + 23} on the bottom strand" in t

    # --- a guide site across the origin of a circle ---
    site = GUIDE + "TGG"
    seq = site[10:] + "ATATTATAATTATATATTAA" * 2 + site[:10]   # starts 13 bases before the end and goes on at the start
    n = len(seq)
    doc = Document.new_sequence("c", seq, True)
    doc.add_features([("misc_feature", f"join({n - 10 + 1}..{n},1..13)", [("label", "x"), (cs.CAS9_KEY, cs.cas9_qualifier(SP))])])
    (f,) = [f for f in doc.features() if f.kind == "misc_feature"]
    layout = cs.feature_layout(f, n)
    assert layout.proto_cells[0] == n - 10 and layout.proto_cells[-1] == 9, "protospacer runs from the end across the origin"
    assert plain(cs.cas9_tip(layout, doc.sequence, n)).count(GUIDE) == 1

    # --- the qualifier is read back, a damaged one is ignored, and the finder and the view agree on the cut ---
    for nuc in cs.NUCLEASES.values():
        assert cs.parse_cas9(cs.cas9_qualifier(nuc)) == nuc
    assert cs.parse_cas9("nonsense") is None and cs.parse_cas9("name=X;pam=NGG") is None and cs.parse_cas9("") is None
    seq = FILL * 3 + GUIDE + "TGG" + FILL
    (g,) = [g for g in cas9.candidates(seq, False, SP) if g.strand > 0]
    assert cs.cas9_layout(SP, [(g.start + i) % len(seq) for i in range(23)], 1, len(seq)).gaps == [g.cut]
    (r,) = [g for g in cas9.candidates(rc(seq), False, SP) if g.strand < 0]
    assert cs.cas9_layout(SP, [r.start + i for i in range(23)], -1, len(seq)).gaps == [r.cut]

    # --- the file: byte for byte, and a plain misc_feature to another reader ---
    tmp = Path(tempfile.mkdtemp(prefix="opendna_guide_"))
    doc, lead, n = guide_doc(1)
    p1, p2 = tmp / "a.gb", tmp / "b.gb"
    doc.save_genbank(str(p1))
    back = Document(str(p1))
    back.save_genbank(str(p2))
    assert p1.read_bytes() == p2.read_bytes(), "load and save again must not change the file"
    (f,) = [f for f in back.features() if f.kind == "misc_feature"]
    layout = cs.feature_layout(f, back.length)
    assert layout is not None and layout.gaps == [lead + 20 - 3]
    assert "".join(back.sequence(a, b) for a, b in layout.proto) == GUIDE
    assert "".join(back.sequence(a, b) for a, b in layout.pam) == "TGG"
    with open(p1, encoding="utf-8") as fh:
        rec = SeqIO.read(fh, "genbank")
    (feat,) = [x for x in rec.features if x.type == "misc_feature"]
    assert feat.qualifiers["label"] == ["gRNA GACTCGTA"] and feat.qualifiers[cs.CAS9_KEY] == [cs.cas9_qualifier(SP)]
    assert GUIDE in feat.qualifiers["note"][0], "the note another reader sees spells the guide out"
    assert int(feat.location.start) == lead and feat.location.strand == 1
    text = p1.read_text(encoding="utf-8")
    assert f'/{cs.CAS9_KEY}="{cs.cas9_qualifier(SP)}"' in text and "OpenGenEdit" not in text, "new files carry the new key"
    assert all(len(line) <= 80 for line in text.splitlines())
    # a reverse-strand guide keeps its strand in the file
    rpath = tmp / "r.gb"
    rdoc, rlead, _ = guide_doc(-1)
    rdoc.save_genbank(str(rpath))
    with open(rpath, encoding="utf-8") as fh:
        (feat,) = [x for x in SeqIO.read(fh, "genbank").features if x.type == "misc_feature"]
    assert feat.location.strand == -1 and int(feat.location.start) == rlead
    assert cs.feature_layout([f for f in Document(str(rpath)).features() if f.kind == "misc_feature"][0],
                             len(rdoc.sequence(0, None))).pam == [(rlead, rlead + 3)]

    # --- a guide saved before the rename to openDNA still works, and exports leave the key out ---
    old = tmp / "old.gb"
    old.write_text(text.replace("/OpenDNA_cas9=", "/OpenGenEdit_cas9="), encoding="utf-8")
    old_ed = Editor(Document(str(old)), old)
    (f,) = [f for f in old_ed.features if f.kind == "misc_feature"]
    assert [k for k, _ in f.qualifiers if k.endswith("_cas9")] == ["OpenGenEdit_cas9"], "the old file is read as it was written"
    layout = cs.feature_layout(f, old_ed.length)
    assert layout is not None and layout.gaps == [lead + 20 - 3]
    assert "OpenGenEdit_cas9" not in exports.tbl_text(old_ed) and "OpenGenEdit_cas9" not in exports.gff3_text(old_ed)
    new_ed = Editor(back, p1)
    assert cs.CAS9_KEY not in exports.tbl_text(new_ed) and cs.CAS9_KEY not in exports.gff3_text(new_ed)
    assert GUIDE in exports.tbl_text(new_ed), "the plain note stays"
