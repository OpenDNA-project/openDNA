"""The codon tools (app/codon_tools.py): recoding a CDS for a codon usage table without ever changing the protein, and reading a CDS
under two genetic codes. The CDSs are the real ones of the demo plasmid; the usage tables are made up (openDNA ships none)."""
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
from PySide6.QtWidgets import QApplication

import opendna_core
from opendna_core import Document

import codon_tools as ct
import codon_usage
import genetic_code
from assembly_panel import FragmentEditor
from editor import reverse_complement
from enzymes import COMMON, Enzyme, find_cuts
from sample_files import PLASMID

app = QApplication.instance() or QApplication([])
BY_NAME = {e.name: e for e in COMMON}
AVOID = [BY_NAME[n] for n in ("EcoRI", "BsaI", "BsmBI", "HindIII", "NcoI", "XbaI", "BamHI", "KpnI")]


def uneven_table(seed=7):
    """A made-up host: every codon between 1 and 40 per thousand, so there is always a clear best codon."""
    rng = random.Random(seed)
    return codon_usage._table("uneven", {c: rng.uniform(1, 40) for c in codon_usage.CODONS})


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


def has_sites(seq, enzymes):
    return sorted(set(find_cuts(seq, False, enzymes)))


@pytest.mark.parametrize("strand", [1, -1])
@pytest.mark.parametrize("method", [ct.MATCHED, ct.HIGHEST])
def test_a_recoded_real_cds_translates_to_the_identical_protein_on_either_strand(strand, method):
    """The one thing that must never be wrong: same protein, read with the opendna core's own translate."""
    ed = make_editor()
    cds(ed, strand=strand)
    src = ct.cds_source(ed)
    result = ct.recode(src.seq, genetic_code.codon_map(src.table), uneven_table(), method=method)
    assert result.seq != src.seq and len(result.seq) == len(src.seq)
    assert opendna_core.translate(result.seq, src.table) == opendna_core.translate(src.seq, src.table)


def test_the_recoded_cds_has_none_of_the_forbidden_sites_obeys_the_gc_window_and_run_limit_and_has_no_rbs():
    ed = make_editor()
    cds(ed, strand=-1)
    src = ct.cds_source(ed)
    result = ct.recode(src.seq, genetic_code.codon_map(src.table), uneven_table(), avoid=AVOID, gc_range=(40, 60), window=50,
                       max_run=5, rbs=True, left=src.left, right=src.right)
    assert result.unmet == []
    assert has_sites(src.left + result.seq + src.right, AVOID) == []   # also across the CDS's two edges
    assert ct.RBS not in result.seq
    assert max(len(run) for run in re.findall(r"A+|C+|G+|T+", result.seq)) <= 5
    windows = [ct.gc(result.seq[i:i + 50]) for i in range(len(result.seq) - 49)]
    assert min(windows) >= 40 and max(windows) <= 60


def test_highest_frequency_always_takes_the_most_used_synonymous_codon():
    table, code = uneven_table(), genetic_code.codon_map(1)
    seq = "".join(random.Random(1).choices(sorted(code), k=300))
    result = ct.recode(seq, code, table, method=ct.HIGHEST, rbs=False, max_run=99)
    best = {aa: max((c for c in code if code[c] == aa), key=lambda c: table.per_thousand[c]) for aa in set(code.values())}
    assert all(best[code[result.seq[i:i + 3]]] == result.seq[i:i + 3] for i in range(0, len(seq), 3))


def test_frequency_matched_follows_the_table_instead_of_always_taking_the_best_codon():
    """Leucine has six codons: sampled in proportion, each should come out close to its share, not the best one every time."""
    table = codon_usage._table("t", {"CTG": 50, "CTA": 10, "CTC": 10, "CTT": 10, "TTA": 10, "TTG": 10})
    code = genetic_code.codon_map(1)
    result = ct.recode("CTG" * 3000, code, table, rbs=False, max_run=99)
    share = {c: sum(result.seq[i:i + 3] == c for i in range(0, 9000, 3)) / 3000 for c in ("CTG", "CTA", "TTA")}
    assert share["CTG"] == pytest.approx(0.5, abs=0.04) and share["CTA"] == pytest.approx(0.1, abs=0.03)


def test_a_codon_the_table_never_uses_is_not_chosen_unless_it_is_the_only_choice():
    table = codon_usage._table("t", {c: 10 for c in codon_usage.CODONS} | {"CTG": 0})
    result = ct.recode("CTG" * 200, genetic_code.codon_map(1), table, rbs=False, max_run=99)
    assert "CTG" not in {result.seq[i:i + 3] for i in range(0, 600, 3)}
    only = codon_usage._table("t", {c: 10 for c in codon_usage.CODONS} | {"TGG": 0})   # Trp has no other codon
    assert ct.recode("TGG" * 10, genetic_code.codon_map(1), only, rbs=False, max_run=99).seq == "TGG" * 10


def test_a_rule_that_cannot_be_met_is_reported_with_where_instead_of_being_hidden():
    """Met-Trp has one codon each, so a 'site' ATGTGG cannot be removed by a silent change."""
    site = Enzyme("Impossible", "ATGTGG", ((1, 5),), False)
    result = ct.recode("GCT" + "ATGTGG" + "GCT", genetic_code.codon_map(1), uneven_table(), avoid=[site], rbs=False)
    assert result.unmet == ["Impossible site at bases 4–9"]
    assert opendna_core.translate(result.seq, 1) == "AMWA"


def test_two_runs_with_the_same_seed_agree_and_another_seed_gives_another_gene():
    ed = make_editor()
    cds(ed, strand=1)
    src = ct.cds_source(ed)
    code, table = genetic_code.codon_map(src.table), uneven_table()
    assert ct.recode(src.seq, code, table, seed=3).seq == ct.recode(src.seq, code, table, seed=3).seq
    assert ct.recode(src.seq, code, table, seed=3).seq != ct.recode(src.seq, code, table, seed=4).seq


def test_a_coding_sequence_that_is_not_whole_codons_is_refused_with_a_sentence():
    with pytest.raises(ValueError, match="whole number of codons"):
        ct.recode("ATGAA", genetic_code.codon_map(1), uneven_table())


def test_applying_a_recoding_is_one_edit_that_undo_restores_exactly_and_features_stay_where_they_were():
    ed = make_editor()
    cds(ed, strand=-1)
    before, features = ed.doc.sequence(), [(g.kind, g.location, g.label) for g in ed.features]
    src = ct.cds_source(ed)
    result = ct.recode(src.seq, genetic_code.codon_map(src.table), uneven_table(), avoid=AVOID)
    assert ed.overwrite_bases(src.start, reverse_complement(result.seq))
    assert ed.doc.sequence() != before and len(ed.doc.sequence()) == len(before)
    assert forward(ed, src) == result.seq
    assert [(g.kind, g.location, g.label) for g in ed.features] == features
    assert ed.undo()
    assert ed.doc.sequence() == before


def test_a_locked_dna_refuses_the_recoding():
    ed = make_editor()
    cds(ed, strand=1)
    ed.dna_locked = True
    before = ed.doc.sequence()
    assert not ed.overwrite_bases(ed.start, ed.doc.sequence(ed.start, ed.end))
    assert ed.doc.sequence() == before


def test_the_selection_is_read_as_a_cds_only_when_it_is_whole_codons_on_one_piece():
    ed = make_editor()
    ed.select(0, 10)
    assert "whole number of codons" in ct.cds_source(ed)
    ed.select(5, 5)
    assert "Select a CDS" in ct.cds_source(ed)
    ed.select(0, 30)
    src = ct.cds_source(ed)
    assert src.seq == ed.doc.sequence(0, 30) and src.strand == 1 and src.start == 0


def test_the_recode_dialog_previews_then_applies_and_one_undo_puts_the_bases_back():
    codon_usage.save_user(uneven_table())
    ed = make_editor()
    cds(ed, strand=1)
    before = ed.doc.sequence()
    dlg = ct.RecodeDialog(None, ed)
    assert dlg.method.currentText() == ct.MATCHED   # the better default
    assert not dlg.apply_button.isEnabled()
    dlg.run_recode()
    assert dlg.apply_button.isEnabled() and "protein is identical" in dlg.report.text() and "GC" in dlg.report.text()
    dlg.rbs.toggle()   # an option changed: the preview no longer describes what Apply would do
    assert not dlg.apply_button.isEnabled()
    dlg.run_recode()
    dlg.apply()
    assert ed.doc.sequence() != before
    assert ed.undo() and ed.doc.sequence() == before


def test_the_recode_dialog_says_why_when_nothing_is_selected_or_there_is_no_usage_table():
    ed = make_editor()
    ed.select(3, 3)
    assert "Select a CDS" in ct.RecodeDialog(None, ed).head.text()
    for p in codon_usage.folder().glob("*.json"):
        p.unlink()
    cds(ed, strand=1)
    assert "no codon usage table" in ct.RecodeDialog(None, ed).head.text()


# ---- the same CDS under two genetic codes ----------------------------------------------------------------------------------
def test_comparing_flags_an_internal_stop_that_appears_under_the_second_code_and_one_that_disappears():
    """TGA is a stop in the standard code and Trp in Mycoplasma (table 4): a gene written for Mycoplasma is cut short in E. coli."""
    seq = "ATG" + "TGA" + "GCT" + "TAA"
    standard, mycoplasma = genetic_code.codon_map(1), genetic_code.codon_map(4)
    moved = ct.compare(seq, mycoplasma, standard)
    assert (moved.a, moved.b) == ("MWA*", "M*A*")
    assert moved.gained == [2] and moved.lost == [] and moved.differ == [2]
    back = ct.compare(seq, standard, mycoplasma)
    assert back.lost == [2] and back.gained == []


def test_comparing_does_not_count_the_final_stop_as_internal_and_sees_no_difference_between_equal_codes():
    seq = "ATGGCTTGA"
    res = ct.compare(seq, genetic_code.codon_map(1), genetic_code.codon_map(11))
    assert res.differ == [] and res.gained == [] and res.lost == []
    final = ct.compare(seq, genetic_code.codon_map(1), genetic_code.codon_map(4))
    assert final.differ == [3] and final.gained == []   # the last codon: a stop became Trp, not an internal stop


def test_the_comparison_drawn_numbers_the_rows_marks_the_differences_and_shows_the_stops():
    res = ct.compare("ATG" + "TGA" * 12 + "TAA", genetic_code.codon_map(4), genetic_code.codon_map(1))
    html = ct.comparison_html(res, "Code 4", "Code 1")
    assert "Code 4" in html and "Code 1" in html and "▲" in html and "*" in html and "14" in html


def test_the_compare_dialog_names_the_internal_stop_and_lets_you_choose_both_codes():
    ed = make_editor()
    cds(ed, strand=1)
    seq = ct.cds_source(ed).seq
    stops = [i for i in range(0, len(seq) - 3, 3) if seq[i:i + 3] == "TGG"]
    assert stops, "the demo CDS has at least one Trp codon to turn into TGA"
    i = stops[0]
    ed.overwrite_bases(ed.start + i, "TGA")   # a Trp that Mycoplasma reads the same as the standard code reads a stop
    ed.select_feature(next(f for f in ed.features if f.kind == "CDS" and f.strand == 1).index)
    dlg = ct.CompareDialog(None, ed)
    dlg.code_a.setCurrentIndex(dlg.code_a.findData(4))
    dlg.code_b.setCurrentIndex(dlg.code_b.findData(1))
    assert "gains an internal stop" in dlg.summary.text() and str(i // 3 + 1) in dlg.summary.text()
    dlg.code_b.setCurrentIndex(dlg.code_b.findData(4))
    assert "same under both codes" in dlg.summary.text()
