"""Plasmid Checkup (app/checkup.py): promoter / terminator / RBS prediction, transcription units and the verdict of each
CDS, on constructs built here from known parts (filler and genes are random, with fixed seeds), plus the transcription
unit features (app/units.py) those checkups add to a document."""
import os
import random
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ["OPENDNA_DATA_DIR"] = tempfile.mkdtemp(prefix="opendna_test_")
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "app"))

import pytest

from opendna_core import Document

import checkup as C
import exports
import genetic_code
import tool_settings
import units
from editor import Editor, reverse_complement

# A consensus sigma-70 promoter (-35 TTGACA, 17 bp spacer, -10 TATAAT), the rrnB T1 terminator, a strong RBS.
PROMOTER = "TTGACAATTAATCATCCGGCTCGTATAATGTGTGGAA"
T1 = "AAGGCTCAGTCGAAAGACTGGGCCTTTCGTTTTATCTGTTGTTTGTCGG"
RBS = "AAGGAGGTAACAT"
KOZAK = "GCCACC"   # gccRcc before ATG; gene() puts a G at +4 only by chance, so the tests set it


def filler(n, seed):
    """No T at all: no -10 box, no T tract, so nothing is predicted in it."""
    rnd = random.Random(seed)
    return "".join(rnd.choice("ACG") for _ in range(n))


def gene(aa=60, seed=1):
    rnd = random.Random(seed)
    standard = genetic_code.TABLES[1].codons
    sense = [c for c in standard if standard[c] != "*"]
    return "ATG" + "".join(rnd.choice(sense) for _ in range(aa)) + "TAA"


POLYA = "AATAAA" + filler(40, 10)


def construct(parts, circular=True):
    """Editor on the concatenated parts; `parts`: [(name or None, kind, seq, strand)] (strand -1: written reversed)."""
    seq, features = "", []
    for name, kind, s, strand in parts:
        s = s if strand > 0 else reverse_complement(s)
        if name:
            loc = f"{len(seq) + 1}..{len(seq) + len(s)}"
            features.append((kind, loc if strand > 0 else f"complement({loc})", [("label", name)]))
        seq += s
    ed = Editor(Document.new_sequence("t", seq, circular), Path("t.gb"))
    if features:
        ed.add_features(features)
    return ed


@pytest.fixture(autouse=True)
def no_library(request, monkeypatch):
    """Only what each test builds is used: the library parts are tested once, by the test named so."""
    if "library" not in request.node.name:
        real = tool_settings.value
        monkeypatch.setattr(tool_settings, "value",
                            lambda tool, key: False if (tool, key) == (C.TOOL, "library") else real(tool, key))


def verdicts(ed):
    return {c.feature.label: c for c in C.run(ed).cds}


def test_bacterial_cds_verdicts():
    """The verdict of each CDS in bacteria: an annotated unit expresses it; with only the CDS annotated, the promoter, RBS
    and terminator are predicted from their motifs (on either strand) and the verdict relies on them; a terminator before
    or inside the CDS blocks / truncates it; the wrong strand and a missing RBS; coding problems and regulation notes."""
    # --- an annotated unit expresses its CDS ---
    ed = construct([(None, "", filler(200, 1), 1), ("Pcons", "promoter", PROMOTER, 1), ("rbs", "RBS", RBS, 1),
                    ("gfp", "CDS", gene(), 1), ("T1", "terminator", T1, 1), (None, "", filler(300, 2), 1)])
    r = C.run(ed)
    c = {c.feature.label: c for c in r.cds}["gfp"]
    assert c.verdict == C.YES and c.unit.promoter.label == "Pcons" and c.rbs.label == "rbs", c.reasons
    unit = next(u for u in r.units if u.promoter.label == "Pcons")
    assert unit.terminator.label == "T1" and unit.cds == ["gfp"]

    # --- predicted elements: nothing annotated but the CDS, on the plus and the minus strand ---
    for strand in (1, -1):
        g = gene()
        whole = PROMOTER + "AAAC" + RBS + g + "AAAC" + T1
        ed = construct([(None, "", filler(150, 3), 1), (None, "", whole, strand), (None, "", filler(150, 4), 1)])
        s = 150 + (len(PROMOTER + "AAAC" + RBS) if strand > 0 else len(T1) + 4)
        loc = f"{s + 1}..{s + len(g)}"
        ed.add_features([("CDS", loc if strand > 0 else f"complement({loc})", [("label", "g")])])
        r = C.run(ed)
        kinds = {(e.kind, e.origin, e.strand) for e in r.elements}
        assert {("promoter", "predicted", strand), ("terminator", "predicted", strand),
                ("RBS", "predicted", strand)} <= kinds, f"strand {strand}"
        c = r.cds[0]
        assert c.verdict == C.LIKELY and "relies on predicted elements" in c.reasons, (strand, c.reasons)

    # --- a terminator before the CDS blocks it ---
    ed = construct([("P", "promoter", PROMOTER, 1), ("T1", "terminator", T1, 1), (None, "", filler(100, 5), 1),
                    ("rbs", "RBS", RBS, 1), ("lost", "CDS", gene(), 1), (None, "", filler(100, 6), 1)])
    c = verdicts(ed)["lost"]
    assert c.verdict == C.NO and "stopped by T1" in c.reasons[0], c.reasons

    # --- a terminator inside the CDS truncates it ---
    g = gene(80)
    ed = construct([("P", "promoter", PROMOTER, 1), ("rbs", "RBS", RBS, 1), (None, "", g[:120], 1),
                    ("Tin", "terminator", T1, 1), (None, "", g[120:], 1), (None, "", filler(100, 7), 1)])
    start = len(PROMOTER + RBS)
    ed.add_features([("CDS", f"{start + 1}..{start + len(g) + len(T1)}", [("label", "cut")])])
    c = verdicts(ed)["cut"]
    assert c.verdict == C.NO and "inside the CDS" in c.reasons[0], c.reasons

    # --- wrong strand and no RBS ---
    ed = construct([("P", "promoter", PROMOTER, 1), (None, "", "AAACAAACAAAC", 1), ("anti", "CDS", gene(), -1),
                    ("T1", "terminator", T1, 1), (None, "", filler(100, 8), 1),
                    ("P2", "promoter", PROMOTER, 1), (None, "", "CCCCCCCCCCCC", 1), ("bare", "CDS", gene(seed=2), 1),
                    ("T1b", "terminator", T1, 1)], circular=False)
    v = verdicts(ed)
    assert v["anti"].verdict == C.NO and "no bacterial promoter upstream" in v["anti"].reasons[0], v["anti"].reasons
    assert v["bare"].verdict == C.UNLIKELY and "no RBS" in v["bare"].reasons[0], v["bare"].reasons

    # --- coding problems and regulation notes ---
    broken = "ATG" + "GCT" * 10 + "TAA" + "GCT" * 10 + "TAA"
    ed = construct([("Ptac", "promoter", PROMOTER, 1), ("rbs", "RBS", RBS, 1), ("broken", "CDS", broken, 1),
                    ("rbs2", "RBS", RBS, 1), ("lacI", "CDS", gene(seed=3), 1), ("T1", "terminator", T1, 1)])
    v = verdicts(ed)
    assert v["broken"].verdict == C.UNLIKELY and any("truncated protein" in r for r in v["broken"].reasons)
    assert any("IPTG-inducible (lacI on the plasmid)" in r for r in v["lacI"].reasons), v["lacI"].reasons
    assert v["lacI"].verdict == C.YES   # inducible is a note, not a doubt


def test_eukaryotic_hosts_and_unit_ends():
    """A CMV unit ends at its polyA signal and needs a Kozak context, not an RBS (a weak Kozak and a missing polyA are
    doubts); a U6 unit (Pol III) is not translated; eukaryotic cells get a column of their own because the plasmid has a
    eukaryotic promoter; on a shuttle vector each unit ends at its own kind of terminator and each CDS gets the verdict of
    each host."""
    # --- eukaryotic and Pol III promoters ---
    egfp = "ATGG" + gene()[4:]
    ed = construct([("CMV promoter", "promoter", filler(200, 9), 1), (None, "", KOZAK, 1), ("egfp", "CDS", egfp, 1),
                    ("pA", "polyA_signal", POLYA, 1),
                    ("U6 promoter", "promoter", filler(200, 11), 1), ("odd", "CDS", gene(seed=4), 1),
                    ("pA2", "polyA_signal", POLYA, 1)])
    r = C.run(ed)
    assert r.hosts == (C.BACTERIA, C.EUKARYOTES)
    v = {c.feature.label: c for c in r.cds}
    egfp_euk = v["egfp"].hosts[C.EUKARYOTES]
    assert egfp_euk.verdict == C.YES and "strong Kozak" in egfp_euk.rbs.label, egfp_euk.reasons
    assert v["egfp"].hosts[C.BACTERIA].verdict == C.NO and v["egfp"].host == C.EUKARYOTES
    assert v["odd"].hosts[C.EUKARYOTES].verdict == C.NO and "Pol III" in v["odd"].hosts[C.EUKARYOTES].reasons[0]

    # --- weak Kozak and no polyA ---
    ed = construct([("CMV promoter", "promoter", filler(200, 14), 1), (None, "", "TTTTTT", 1),
                    ("weak", "CDS", "ATGA" + gene()[4:], 1), (None, "", filler(300, 15), 1)], circular=False)
    euk = verdicts(ed)["weak"].hosts[C.EUKARYOTES]
    assert any("weak Kozak" in x for x in euk.reasons) and any("no polyA" in x for x in euk.reasons)
    assert euk.verdict == C.UNLIKELY

    # --- bacterial and eukaryotic units stop at their own terminators ---
    egfp = "ATGG" + gene(seed=5)[4:]
    ed = construct([("Pcons", "promoter", PROMOTER, 1), ("rbs", "RBS", RBS, 1), ("kanR", "CDS", gene(seed=6), 1),
                    ("pA", "polyA_signal", POLYA, 1), ("rbs2", "RBS", RBS, 1), ("after", "CDS", gene(seed=7), 1),
                    ("T1", "terminator", T1, 1), (None, "", filler(100, 16), 1),
                    ("CMV promoter", "promoter", filler(200, 17), 1), (None, "", KOZAK, 1), ("egfp", "CDS", egfp, 1),
                    ("T1b", "terminator", T1, 1), ("bGH polyA", "polyA_signal", POLYA, 1),
                    (None, "", filler(100, 18), 1)])
    r = C.run(ed)
    by = {u.promoter.label: u for u in r.units}
    assert by["Pcons"].terminator.label == "T1", "polyA does not stop the bacterial polymerase"
    assert by["CMV promoter"].terminator.label == "bGH polyA", "rrnB T1 does not stop the eukaryotic one"
    v = {c.feature.label: c for c in r.cds}
    assert v["after"].hosts[C.BACTERIA].verdict == C.YES
    assert v["kanR"].hosts[C.EUKARYOTES].verdict == C.NO and v["kanR"].hosts[C.BACTERIA].verdict == C.YES
    assert v["egfp"].hosts[C.EUKARYOTES].verdict == C.YES and v["egfp"].hosts[C.BACTERIA].verdict == C.NO


def test_how_far_a_unit_runs_and_how_it_ends():
    """With no terminator a unit runs round a circle or to the end of a line (the far CDS gets a note on its distance);
    a unit can pass the origin and be added as an mRNA feature; +1 is exact only where the start base is known; units end
    in their own ways (intrinsic, rho-dependent?, polyA, runs on)."""
    # --- a single promoter without a terminator ---
    for circular in (True, False):
        ed = construct([("P", "promoter", PROMOTER, 1), ("rbs", "RBS", RBS, 1), ("near", "CDS", gene(), 1),
                        (None, "", filler(3500, 19), 1), ("rbs2", "RBS", RBS, 1), ("far", "CDS", gene(seed=8), 1),
                        (None, "", filler(100, 20), 1)], circular)
        r = C.run(ed)
        u = r.units[0]
        assert u.terminator is None and u.b - u.a == ed.length - (0 if circular else u.a)
        v = {c.feature.label: c for c in r.cds}
        assert v["near"].verdict == v["far"].verdict == C.YES
        assert any("no terminator" in x for x in v["near"].reasons)
        assert any("kb downstream of P" in x for x in v["far"].reasons), v["far"].reasons

    # --- a unit around the origin, added as an mRNA feature ---
    ed = construct([("rbs", "RBS", RBS, 1), ("g", "CDS", gene(), 1), ("T1", "terminator", T1, 1),
                    (None, "", filler(300, 12), 1), ("P", "promoter", PROMOTER, 1)])
    r = C.run(ed)
    assert r.cds[0].verdict == C.YES, r.cds[0].reasons
    u = r.units[0]
    assert C.unit_extent(u, ed.length)[1] > ed.length, "the unit passes the origin"
    ed.add_features([C.unit_feature(u, ed.length)])
    mrna = [f for f in ed.features if f.kind == "mRNA"]
    assert len(mrna) == 1 and len(mrna[0].spans) == 2

    # --- exact +1 only: an annotated TSS or a T7 promoter; a promoter's 3' end or a predicted -10 box is an estimate ---
    t7 = "TAATACGACTCACTATAGGG"
    ed = construct([("T7 promoter", "promoter", t7[:19], 1), ("g", "CDS", gene(), 1), ("T1", "terminator", T1, 1),
                    (None, "", filler(60, 37), 1), ("Pcons", "promoter", PROMOTER, 1), ("TSS", "misc_feature", "A", 1),
                    ("g2", "CDS", gene(seed=10), 1), ("T1b", "terminator", T1, 1), (None, "", filler(60, 38), 1),
                    ("Pplain", "promoter", PROMOTER, 1), ("g3", "CDS", gene(seed=11), 1), ("T1c", "terminator", T1, 1)],
                   False)
    by = {u.promoter.label: u for u in C.run(ed).units}
    assert by["T7 promoter"].tss_exact and by["T7 promoter"].tss == 17
    tss_at = len(t7[:19] + gene() + T1) + 60 + len(PROMOTER)
    assert by["Pcons"].tss_exact and by["Pcons"].tss == tss_at, "the annotated TSS base"
    assert not by["Pplain"].tss_exact
    ed.add_features([C.unit_feature(u, ed.length) for u in by.values()])
    exact = {f.label.split(" →")[0]: units.geometry(f, ed.length, False).exact for f in ed.features if units.is_tu(f)}
    assert exact == {"T7 promoter": True, "Pcons": True, "Pplain": False}

    # --- how units end ---
    ed = construct([("Pa", "promoter", PROMOTER, 1), ("T1", "terminator", T1, 1), (None, "", filler(50, 33), 1),
                    ("Pb", "promoter", PROMOTER, 1), ("tR?", "terminator", filler(60, 32), 1), (None, "", filler(50, 34), 1),
                    ("CMV promoter", "promoter", filler(150, 35), 1), ("bGH polyA", "polyA_signal", POLYA, 1)], False)
    by = {u.promoter.label: u for u in C.run(ed).units}
    assert by["Pa"].ending == "intrinsic" and "hairpin" in by["Pa"].ending_detail
    assert by["Pb"].ending == "rho-dependent?"
    assert by["CMV promoter"].ending == "polyA"
    assert C.unit_label(by["Pa"]) == "Pa → T1 (intrinsic)" and C.unit_label(by["CMV promoter"]) == "CMV promoter → bGH polyA"
    lone = construct([("P", "promoter", PROMOTER, 1), (None, "", filler(80, 36), 1)])
    u = C.run(lone).units[0]
    assert u.ending == "runs on" and C.unit_label(u) == "P → runs on"


def test_origins_and_prediction_thresholds():
    """Replication origins are named (with the pir note and the ColE1-group / no-origin remarks), the thresholds keep
    chance hits rare in a 50 kb random sequence, and motifs near the start of a linear sequence are placed where they are."""
    # --- replication origins ---
    ed = construct([("pUC ori", "rep_origin", filler(600, 21), 1), ("p15A ori", "rep_origin", filler(500, 22), 1),
                    ("ColE1 origin", "rep_origin", filler(600, 23), 1), ("R6K gamma ori", "rep_origin", filler(400, 24), 1),
                    ("f1 ori", "rep_origin", filler(450, 25), 1)])
    r = C.run(ed)
    names = {o.element.label: o for o in r.origins}
    assert names["pUC ori"].name == "pUC" and names["p15A ori"].name == "p15A"
    assert names["ColE1 origin"].name == "pMB1 / ColE1" and names["f1 ori"].host == "Phage"
    assert "pir" in names["R6K gamma ori"].note
    with_pir = construct([("R6K gamma ori", "rep_origin", filler(400, 27), 1), ("pir", "CDS", gene(seed=12), 1)])
    assert C.run(with_pir).origins[0].note == "pir is on the plasmid"
    assert any("2 origins of the ColE1 group" in x for x in r.origin_notes)
    lonely = construct([("f1 ori", "rep_origin", filler(450, 26), 1), ("g", "CDS", gene(), 1)])
    assert any("no bacterial replication origin" in x for x in C.run(lonely).origin_notes)

    # --- a random sequence has few predictions ---
    rnd = random.Random(7)
    st = C.Strand("".join(rnd.choice("ACGT") for _ in range(50_000)), 1, True)
    assert len(C.predict_promoters(st, 10)) <= 8
    assert len(C.predict_terminators(st, 18)) <= 8

    # --- motifs near the start of a linear sequence: a line's slices start at its first base ---
    assert C.sd_match(C.Strand("AGGAGGTTTTATG" + "GCT" * 20, 1, False), 10, 20) == (6, 0, 6)


def test_library_parts_on_a_real_plasmid():
    """pBAD30 with eilR behind the araBAD promoter (araC on the plasmid). The file annotates its own
    terminator and origins, so the library adds nothing here — it is there for what a file leaves out,
    and the same sequence without annotations gets its terminator from the library instead."""
    path = ROOT / "demo" / "pBAD30_eilR.gb"
    ed = Editor(Document(str(path)), path)
    r = C.run(ed)
    assert not any(e.origin == "library" for e in r.elements), "nothing to add to an annotated plasmid"
    bare = Editor(Document.new_sequence("bare", ed.doc.sequence(), True), path)
    assert any(e.origin == "library" and e.kind == "terminator" for e in C.run(bare).elements)
    v = {c.feature.label: c for c in r.cds}
    assert v["eilR"].verdict == C.YES and v["eilR"].unit.promoter.label == "araBAD promoter", v["eilR"].reasons
    assert any("arabinose" in x and "araC on the plasmid" in x for x in v["eilR"].reasons)


def test_unit_features_keep_their_geometry_and_survive_a_file_round_trip(tmp_path):
    """A transcription unit written as a feature: its +1, its thin promoter piece and thick transcript, on both strands
    and across the origin (a plain mRNA is not a unit); its own display switch; and a file saved before the rename to
    openDNA (/OpenGenEdit_tu) still reads as the same units, while exports leave the program's qualifier out."""
    ed = construct([("Pcons", "promoter", PROMOTER, 1), ("rbs", "RBS", RBS, 1), ("gfp", "CDS", gene(), 1),
                    ("T1", "terminator", T1, 1), (None, "", filler(150, 30), 1), ("Prev", "promoter", PROMOTER, -1),
                    (None, "", filler(100, 31), 1)])
    ed.add_features([C.unit_feature(u, ed.length) for u in C.run(ed).units])

    # --- geometry on both strands and across the origin ---
    fwd = next(f for f in ed.features if units.is_tu(f) and f.strand > 0)
    g = units.geometry(fwd, ed.length, True)
    tss = len(PROMOTER)
    assert g.tss == tss and g.stopped and g.end == len(PROMOTER + RBS + gene() + T1)
    assert g.pieces == [(0, tss, False), (tss, g.end, True)], "thin over the promoter, thick from +1"
    rev = next(f for f in ed.features if units.is_tu(f) and f.strand < 0)
    g = units.geometry(rev, ed.length, True)
    assert not g.stopped, "nothing stops it: it runs on (all round the circle)"
    p_end = len(PROMOTER + RBS + gene() + T1) + 150   # Prev's left end, where +1 is on the minus strand
    assert g.tss == p_end
    thick = [(a, b) for a, b, heavy in g.pieces if heavy]
    thin = [(a, b) for a, b, heavy in g.pieces if not heavy]
    assert thin == [(p_end, p_end + len(PROMOTER))] and sum(b - a for a, b in thick) == ed.length - len(PROMOTER)
    assert len(rev.spans) == 2, "a full turn starting inside the sequence crosses the origin"
    plain = next(f for f in ed.features if f.label == "gfp")
    assert units.geometry(plain, ed.length, True) is None and units.tu_of(plain) is None, "a plain feature is not a unit"

    # --- units have a display switch of their own, apart from the mRNA type ---
    assert sum(units.is_tu(f) for f in ed.shown) == 2
    ed.set_display(units=False)
    assert not any(units.is_tu(f) for f in ed.shown) and not any(units.is_tu(f) for f in ed.map_features)
    assert exports.project_meta(ed, "x.odp")["display"]["units"] is False
    ed.set_display(hidden_kinds={"mRNA"}, units=True)
    assert sum(units.is_tu(f) for f in ed.shown) == 2, "hiding the mRNA type leaves units to their own switch"

    # --- the file: the new key, and units saved under the old one still read back the same ---
    new, old = tmp_path / "new.gb", tmp_path / "old.gb"
    ed.doc.save_genbank(str(new))
    text = new.read_text(encoding="utf-8")
    assert "/OpenDNA_tu=" in text and "OpenGenEdit" not in text, "new files carry the new key"
    old.write_text(text.replace("/OpenDNA_tu=", "/OpenGenEdit_tu="), encoding="utf-8")
    back = Editor(Document(str(old)), old)
    before = [units.geometry(f, ed.length, True) for f in ed.features if units.is_tu(f)]
    after = [units.geometry(f, back.length, True) for f in back.features if units.is_tu(f)]
    assert len(after) == 2 and after == before, "the old units read back exactly as the new ones"
    assert sum(units.is_tu(f) for f in back.shown) == 2
    back.set_display(units=False)
    assert not any(units.is_tu(f) for f in back.shown)
    assert units.TU_KEY not in exports.gff3_text(ed), "the program's own qualifier stays out of exports"
    assert "OpenGenEdit_tu" not in exports.gff3_text(back) and "OpenGenEdit_tu" not in exports.tbl_text(back)
