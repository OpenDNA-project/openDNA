"""Codon work: recoding a CDS for a codon usage table without ever changing the protein, reading a CDS under two genetic
codes, and the usage tables themselves (CSV in and out, tables of your own, the adaptation index). openDNA ships no codon
usage table (no licence allows it), so every table here is made up."""
import os
import random
import re
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ["OPENDNA_DATA_DIR"] = tempfile.mkdtemp(prefix="opendna_test_")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "app"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import pytest

import opendna_core
from opendna_core import Document

import checkup
import codon_tools as ct
import codon_usage as cu
import genetic_code
from assembly_panel import FragmentEditor
from editor import reverse_complement
from enzymes import COMMON, Enzyme, find_cuts
from sample_files import PLASMID

BY_NAME = {e.name: e for e in COMMON}
AVOID = [BY_NAME[n] for n in ("EcoRI", "BsaI", "BsmBI", "HindIII", "NcoI", "XbaI", "BamHI", "KpnI")]


def uneven_table(seed=7):
    """A made-up host: every codon between 1 and 40 per thousand, so there is always a clear best codon."""
    rng = random.Random(seed)
    return cu._table("uneven", {c: rng.uniform(1, 40) for c in cu.CODONS})


def flat(name="flat", **changes):
    """A table where every codon is used 15.625 per thousand, except the codons in `changes`."""
    freq = {c: 15.625 for c in cu.CODONS}
    freq.update(changes)
    return cu._table(name, freq)


def make_editor():
    return FragmentEditor(Document(PLASMID), Path("t.gb"))


def cds(ed, label=None, strand=None):
    """A real CDS of the demo plasmid (the one called `label`, or the first one on `strand`), selected."""
    f = next(f for f in ed.features if f.kind == "CDS" and (label is None or f.label == label) and (strand is None or f.strand == strand))
    ed.select_feature(f.index)
    return f


def forward(ed, src):
    """The CDS as the coding strand now reads in the document."""
    bases = ed.doc.sequence(src.start, src.end)
    return bases if src.strand > 0 else reverse_complement(bases)


@pytest.fixture
def clean():
    """No table of yours; the table in use left as found."""
    for p in cu.folder().glob("*.json"):
        p.unlink()
    saved = cu.active_name()
    yield
    cu.set_active(saved)
    for p in cu.folder().glob("*.json"):
        p.unlink()


# ---- recoding -------------------------------------------------------------------------------------------------------


def test_recoding_preserves_the_protein_and_obeys_the_constraints():
    """The one thing that must never be wrong: the same protein, on either strand and with either method. On top of
    that the recoded CDS carries none of the forbidden sites (its two edges included), stays inside the GC window, keeps
    its runs short and holds no RBS — and a rule that cannot be met is reported with where, not hidden."""
    # --- same protein, either strand, either method ---
    for strand in (1, -1):
        for method in (ct.MATCHED, ct.HIGHEST):
            ed = make_editor()
            cds(ed, strand=strand)
            src = ct.cds_source(ed)
            result = ct.recode(src.seq, genetic_code.codon_map(src.table), uneven_table(), method=method)
            assert result.seq != src.seq and len(result.seq) == len(src.seq), (strand, method)
            assert opendna_core.translate(result.seq, src.table) == opendna_core.translate(src.seq, src.table), (strand, method)

    # --- forbidden sites, the GC window, the run limit and the RBS ---
    ed = make_editor()
    cds(ed, strand=-1)
    src = ct.cds_source(ed)
    result = ct.recode(src.seq, genetic_code.codon_map(src.table), uneven_table(), avoid=AVOID, gc_range=(40, 60), window=50,
                       max_run=5, rbs=True, left=src.left, right=src.right)
    assert result.unmet == []
    assert sorted(set(find_cuts(src.left + result.seq + src.right, False, AVOID))) == [], "also across the CDS's two edges"
    assert ct.RBS not in result.seq
    assert max(len(run) for run in re.findall(r"A+|C+|G+|T+", result.seq)) <= 5
    windows = [ct.gc(result.seq[i:i + 50]) for i in range(len(result.seq) - 49)]
    assert min(windows) >= 40 and max(windows) <= 60

    # --- a rule that cannot be met: Met-Trp has one codon each, so ATGTGG cannot be silently removed ---
    site = Enzyme("Impossible", "ATGTGG", ((1, 5),), False)
    result = ct.recode("GCT" + "ATGTGG" + "GCT", genetic_code.codon_map(1), uneven_table(), avoid=[site], rbs=False)
    assert result.unmet == ["Impossible site at bases 4–9"]
    assert opendna_core.translate(result.seq, 1) == "AMWA"


def test_codon_choice_follows_the_table_and_the_edit_is_one_undo_step():
    """Highest frequency always takes the most used codon, frequency-matched follows the table's shares, a codon the
    table never uses is avoided unless it is the only one; the same seed gives the same gene and another seed another
    one; a selection that is not whole codons is refused; applying a recoding is one edit that undo puts back, and a
    locked DNA refuses it."""
    code = genetic_code.codon_map(1)
    table = uneven_table()
    # --- highest frequency ---
    seq = "".join(random.Random(1).choices(sorted(code), k=300))
    result = ct.recode(seq, code, table, method=ct.HIGHEST, rbs=False, max_run=99)
    best = {aa: max((c for c in code if code[c] == aa), key=lambda c: table.per_thousand[c]) for aa in set(code.values())}
    assert all(best[code[result.seq[i:i + 3]]] == result.seq[i:i + 3] for i in range(0, len(seq), 3))

    # --- frequency-matched: Leucine's six codons come out close to their share ---
    leu = cu._table("t", {"CTG": 50, "CTA": 10, "CTC": 10, "CTT": 10, "TTA": 10, "TTG": 10})
    result = ct.recode("CTG" * 3000, code, leu, rbs=False, max_run=99)
    share = {c: sum(result.seq[i:i + 3] == c for i in range(0, 9000, 3)) / 3000 for c in ("CTG", "CTA")}
    assert share["CTG"] == pytest.approx(0.5, abs=0.04) and share["CTA"] == pytest.approx(0.1, abs=0.03)

    # --- a codon the table never uses is not chosen unless it is the only choice ---
    never = cu._table("t", {c: 10 for c in cu.CODONS} | {"CTG": 0})
    out = ct.recode("CTG" * 200, code, never, rbs=False, max_run=99)
    assert "CTG" not in {out.seq[i:i + 3] for i in range(0, 600, 3)}
    only = cu._table("t", {c: 10 for c in cu.CODONS} | {"TGG": 0})   # Trp has no other codon
    assert ct.recode("TGG" * 10, code, only, rbs=False, max_run=99).seq == "TGG" * 10

    # --- the same seed gives the same gene; a sequence that is not whole codons is refused ---
    ed = make_editor()
    cds(ed, strand=1)
    src = ct.cds_source(ed)
    scode = genetic_code.codon_map(src.table)
    assert ct.recode(src.seq, scode, table, seed=3).seq == ct.recode(src.seq, scode, table, seed=3).seq
    assert ct.recode(src.seq, scode, table, seed=3).seq != ct.recode(src.seq, scode, table, seed=4).seq
    with pytest.raises(ValueError, match="whole number of codons"):
        ct.recode("ATGAA", code, table)

    # --- what the selection is read as ---
    ed.select(0, 10)
    assert "whole number of codons" in ct.cds_source(ed)
    ed.select(5, 5)
    assert "Select a CDS" in ct.cds_source(ed)
    ed.select(0, 30)
    picked = ct.cds_source(ed)
    assert picked.seq == ed.doc.sequence(0, 30) and picked.strand == 1 and picked.start == 0

    # --- applying is one edit, undo puts the bases back, and the features stay where they were ---
    ed = make_editor()
    cds(ed, strand=-1)
    before, features = ed.doc.sequence(), [(g.kind, g.location, g.label) for g in ed.features]
    src = ct.cds_source(ed)
    result = ct.recode(src.seq, genetic_code.codon_map(src.table), table, avoid=AVOID)
    assert ed.overwrite_bases(src.start, reverse_complement(result.seq))
    assert ed.doc.sequence() != before and len(ed.doc.sequence()) == len(before)
    assert forward(ed, src) == result.seq
    assert [(g.kind, g.location, g.label) for g in ed.features] == features
    assert ed.undo() and ed.doc.sequence() == before

    # --- a locked DNA refuses it ---
    ed.dna_locked = True
    assert not ed.overwrite_bases(ed.start, ed.doc.sequence(ed.start, ed.end))
    assert ed.doc.sequence() == before


def test_comparing_a_cds_under_two_genetic_codes():
    """TGA is a stop in the standard code and Trp in Mycoplasma (table 4): a gene written for Mycoplasma is cut short in
    E. coli. The final stop is not counted as internal, equal codes differ nowhere, and the drawn comparison numbers its
    rows and marks the differences."""
    seq = "ATG" + "TGA" + "GCT" + "TAA"
    standard, mycoplasma = genetic_code.codon_map(1), genetic_code.codon_map(4)
    moved = ct.compare(seq, mycoplasma, standard)
    assert (moved.a, moved.b) == ("MWA*", "M*A*")
    assert moved.gained == [2] and moved.lost == [] and moved.differ == [2]
    back = ct.compare(seq, standard, mycoplasma)
    assert back.lost == [2] and back.gained == []
    # the final stop is not an internal one, and two codes that say the same thing differ nowhere
    same = ct.compare("ATGGCTTGA", genetic_code.codon_map(1), genetic_code.codon_map(11))
    assert same.differ == [] and same.gained == [] and same.lost == []
    final = ct.compare("ATGGCTTGA", genetic_code.codon_map(1), genetic_code.codon_map(4))
    assert final.differ == [3] and final.gained == []   # the last codon: a stop became Trp, not an internal stop
    # the drawing numbers the rows, marks the differences and shows the stops
    html = ct.comparison_html(ct.compare("ATG" + "TGA" * 12 + "TAA", genetic_code.codon_map(4), genetic_code.codon_map(1)),
                              "Code 4", "Code 1")
    assert "Code 4" in html and "Code 1" in html and "▲" in html and "*" in html and "14" in html


# ---- the usage tables and the genetic codes --------------------------------------------------------------------------


def test_usage_tables_their_files_and_what_they_say_about_a_gene(clean, tmp_path, monkeypatch):
    """No table ships with the app; a table (with the genetic code beside it) survives a CSV round trip and bad files are
    refused; tables of yours live in the data folder and one that cannot be written is reported, not lost; the adaptation
    index and the rare codons follow the table, and so do the checkup's notes; a recoded code reads its removed codons
    as X."""
    import fileio
    # --- nothing ships, and nothing breaks without a table ---
    assert cu.all_tables() == {} and cu.active() is None
    assert checkup.codon_notes("ATG" + "CTA" * 20 + "TAA") == []   # no table: no note, no error

    # --- CSV in and out, and the files that are refused ---
    t = flat(AAA=40.0)
    t.name, t.source = "My organism", "a paper"
    cu.write_csv(t, tmp_path / "t.csv")
    back = cu.read_csv(tmp_path / "t.csv")
    assert back.name == "My organism" and back.source == "a paper" and back.per_thousand["AAA"] == 40.0 and len(back.per_thousand) == 64
    pairs = tmp_path / "pairs.txt"   # a plain list of codon, frequency (U accepted for T, no header)
    pairs.write_text("\n".join(f"{c.replace('T', 'U')}\t{i}" for i, c in enumerate(cu.CODONS)), encoding="utf-8")
    assert cu.read_csv(pairs).per_thousand["TTT"] == 0 and cu.read_csv(pairs).per_thousand["GGG"] == 63
    short = tmp_path / "short.csv"
    short.write_text("codon,aa,per_thousand\nTTT,F,10\n", encoding="utf-8")
    with pytest.raises(ValueError, match="64 codons"):
        cu.read_csv(short)
    bad = tmp_path / "bad.csv"
    bad.write_text("codon,aa,per_thousand\nTTT,F,10\nXYZ,F,10\n", encoding="utf-8")
    with pytest.raises(ValueError, match="not a codon"):
        cu.read_csv(bad)
    rows = "\n".join(f"{c},{'B' if c == 'TTT' else 'F'},15.6" for c in cu.CODONS)
    (tmp_path / "aa.csv").write_text("codon,aa,per_thousand\n" + rows, encoding="utf-8")
    with pytest.raises(ValueError, match="not an amino acid"):
        cu.read_csv(tmp_path / "aa.csv")
    plain = tmp_path / "plain.csv"
    plain.write_text("codon,per_thousand\n" + "\n".join(f"{c},{15.625}" for c in cu.CODONS), encoding="utf-8")
    assert cu.read_csv(plain).aas is None, "a file with no amino-acid column is still just a usage table"

    # --- one CSV carries the usage table and the genetic code shown beside it ---
    aas, starts = genetic_code.recoded(11, ["AGA", "AGG"])
    code_id = genetic_code.save_code("Export test code", aas, starts)
    try:
        both = flat("Export test host", CTG=50.0, CTA=0.5)
        path = tmp_path / "both.csv"
        cu.write_csv(both, path, genetic_code.TABLES[code_id])
        text = path.read_text(encoding="utf-8")
        assert "# code: Export test code" in text and "codon,aa,per_thousand,start" in text
        assert "AGA,X," in text and "ATG,M," in text   # the removed codon and the start codon travel too
        read = cu.read_csv(path)
        assert read.aas == aas and read.starts == starts and read.code_name == "Export test code"
        assert read.per_thousand["CTG"] == 50.0 and read.per_thousand["CTA"] == 0.5
    finally:
        genetic_code.delete_code(code_id)

    # --- tables of yours live in the data folder ---
    mine = cu.duplicate(flat("Original"), "Recoded host")
    cu.save_user(mine)
    assert (cu.folder() / "Recoded_host.json").is_file() and cu.folder().parent == Path(os.environ["OPENDNA_DATA_DIR"])
    assert list(cu.all_tables()) == ["Recoded host"] and cu.all_tables()["Recoded host"].note == "Copy of Original."
    with pytest.raises(ValueError, match="needs a name"):
        cu.save_user(cu.duplicate(mine, "  "))
    (cu.folder() / "broken.json").write_text("{not json", encoding="utf-8")
    assert list(cu.user_tables()) == ["Recoded host"]   # a stray file is ignored
    cu.set_active("Recoded host")
    assert cu.active().name == "Recoded host"
    cu.delete_user("Recoded host")
    assert cu.all_tables() == {} and cu.active_name() == "" and cu.active() is None   # the table in use went

    # --- a table that cannot be written is reported, not lost ---
    def locked(tmp, target, *a, **k):
        raise PermissionError(13, "being used by another process")
    monkeypatch.setattr(fileio, "replace", locked)
    with pytest.raises(ValueError, match="could not be saved.*being used by another process"):
        cu.save_user(flat("Mine"))
    assert cu.all_tables() == {}
    aas, starts = genetic_code.recoded(11, ["AGA"])
    with pytest.raises(ValueError, match="could not be saved"):
        genetic_code.save_code("A code of mine", aas, starts)
    assert not any(t.name == "A code of mine" for t in genetic_code.TABLES.values())
    monkeypatch.undo()

    # --- the adaptation index and the rare codons follow the table ---
    code = genetic_code.codon_map(1)
    table = flat(CTG=50.0, CTA=1.0)   # Leu: CTG is the best codon, CTA is used 1/50 as much
    assert cu.cai("CTGCTGCTG", table, code) == pytest.approx(1.0)
    assert cu.cai("CTACTA", table, code) == pytest.approx(1 / 50)
    mixed = cu.cai("CTGCTA", table, code)
    assert 1 / 50 < mixed < 1 and mixed == pytest.approx((1 / 50) ** 0.5)   # the geometric mean
    assert cu.cai("ATGTGGTAA", table, code) is None   # Met, Trp and the stop have no synonym: nothing to score
    assert cu.rare_codons("ATGCTGCTAGGT", table, code) == [(3, "CTA")]
    assert cu.cai("CTA", flat(CTA=0.0, CTG=10.0), code) == pytest.approx(0.01)   # never used counts as 0.01, not 0
    removed_free = cu.adaptiveness(flat(CTA=0.0), code)
    assert removed_free["CTA"] == 0 and "TAA" not in removed_free

    # --- the checkup's note about a poor adaptation ---
    cu.save_user(flat("Host", CTG=50.0, CTA=0.5))
    cu.set_active("Host")
    assert checkup.codon_notes("ATG" + "CTG" * 20 + "TAA") == []   # all best codons: nothing to say
    notes = checkup.codon_notes("ATG" + "CTA" * 20 + "TAA")
    assert len(notes) == 1 and "adaptation index 0.01 for Host" in notes[0] and "20 rare codons (CTA)" in notes[0]
    cu.set_active("Missing table")   # a table that is gone: another one is used, never an error
    assert isinstance(checkup.codon_notes("ATG" + "CTA" * 20 + "TAA"), list)

    # --- recoded genetic codes: the removed codons read as X, and codes of yours live in the data folder ---
    assert genetic_code.removed_codons(101) == ["TCA", "TCG", "TAG"] and len(genetic_code.removed_codons(102)) == 7
    assert opendna_core.translate("ATGTCGTAGTCATAA", 101) == "MXXX*"   # TCG, TAG, TCA are not used; TAA still stops
    assert genetic_code.table_id() in genetic_code.TABLES
    aas, starts = genetic_code.recoded(11, ["AGA"])
    code_id = genetic_code.save_code("Test recoded host", aas, starts, note="made in a test")
    try:
        assert code_id >= genetic_code.FIRST_USER and f"{code_id} – Test recoded host" in genetic_code.LABELS
        assert (genetic_code.folder() / f"{code_id}.json").is_file()
        assert opendna_core.translate("AGAAGG", code_id) == "XR"   # AGA removed, AGG still arginine in table 11
        with pytest.raises(ValueError, match="already a genetic code"):
            genetic_code.save_code("test RECODED host", aas, starts)
        with pytest.raises(ValueError, match="built-in"):
            genetic_code.save_code("Edit a built-in", aas, starts, id=101)
        with pytest.raises(ValueError, match="64"):
            genetic_code.save_code("Short", "FFLL", "----")
    finally:
        genetic_code.delete_code(code_id)
    assert code_id not in genetic_code.TABLES and not (genetic_code.folder() / f"{code_id}.json").exists()

    # --- and the checkup flags the codons the chosen code does not use ---
    import tool_settings
    saved = tool_settings.value(genetic_code.TOOL, "table")
    try:
        tool_settings.set_value(genetic_code.TOOL, "table", next(x for x in genetic_code.LABELS if x.startswith("101 ")))
        assert "uses codons this genetic code does not use: TCA, TCG" in checkup.protein_problems("ATGTCGAAATCATAA")
        tool_settings.set_value(genetic_code.TOOL, "table", genetic_code.DEFAULT_LABEL)
        assert checkup.protein_problems("ATGTCGAAATCATAA") == []
    finally:
        tool_settings.set_value(genetic_code.TOOL, "table", saved)
