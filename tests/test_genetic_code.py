"""Genetic codes: every embedded NCBI table against Biopython (used here only, never by the app), translation through the
protein search, alternative start codons and table-dependent stops in the ORF finder, a CDS read with its own /transl_table,
/codon_start and /transl_except, exports following the setting, and the /translation refresh."""
import itertools
import os
import sys
import tempfile
from pathlib import Path

os.environ["OPENDNA_DATA_DIR"] = tempfile.mkdtemp(prefix="opendna_test_")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "app"))

import pytest
from Bio.Data import CodonTable
from Bio.Seq import Seq

import opendna_core as core
import genetic_code
import tool_settings
from editor import Editor, reverse_complement
from opendna_core import Document, find_protein, translate
from sequence_view import cds_codons, codons

COMP = str.maketrans("ACGTacgt", "TGCAtgca")
ALL_CODONS = ["".join(c) for c in itertools.product("TCAG", repeat=3)]


def rc(s):
    return s.translate(COMP)[::-1]


# Tables left out on purpose: the meaning of a codon there depends on its place in the gene (see core/src/code.rs).
LEFT_OUT = {27, 28, 31}
IDS = sorted(CodonTable.unambiguous_dna_by_id.keys() - LEFT_OUT)


@pytest.fixture
def code_setting():
    """Leave the genetic code settings as found (the settings folder is shared by all test modules)."""
    saved = {k: tool_settings.value(genetic_code.TOOL, k) for k in ("table", "alt_starts")}
    yield
    for k, v in saved.items():
        tool_settings.set_value(genetic_code.TOOL, k, v)



def choose(table, alt=False):
    label = next(x for x in genetic_code.LABELS if x.split()[0] == str(table))
    tool_settings.set_value(genetic_code.TOOL, "table", label)
    tool_settings.set_value(genetic_code.TOOL, "alt_starts", alt)


def cds_doc(seq, *features):
    doc = Document.new_sequence("t", seq, False)
    doc.add_features(list(features))
    return doc


def cds(location, label="g", **quals):
    return ("CDS", location, [("label", label)] + [(k, v) for k, v in quals.items()])



def test_translation_tables_match_biopython():
    """Every codon of every embedded table on both strands, and the start / stop sets, match Biopython; tables differ where
    they should; Rust and Python agree; the protein search translates all 64 codons, incomplete codons and N; ORFs follow
    the table's starts (alternative ones only when allowed) and stops."""
    # --- the embedded tables are the expected ones
    assert [i for i in sorted(genetic_code.TABLES) if not genetic_code.is_custom(i)] == IDS   # (the recoded codes, 100+, are not NCBI's)
    assert IDS[0] == 1 and 11 in IDS and 4 in IDS and 2 in IDS
    # --- every codon of every table matches Biopython on both strands, and the start and stop sets
    for table in IDS:
        ref = CodonTable.unambiguous_dna_by_id[table]
        bases = "ACGT"
        all_codons = [a + b + c for a in bases for b in bases for c in bases]
        expected = {c: "*" if c in ref.stop_codons else ref.forward_table[c] for c in all_codons}
        assert len(expected) == 64
        for c in all_codons:
            assert translate(c, table) == expected[c], (table, c)
            assert genetic_code.codon_map(table)[c] == expected[c], (table, c)
        # The reverse strand: a CDS on the bottom strand reads the reverse complement of the top strand.
        top = "".join(reverse_complement(c) for c in all_codons)
        doc = Document.new_sequence("rev", top, False)
        got = codons(doc, [(0, len(top))], -1, 0, table)
        # Reading backwards, the codons come out in the reverse order of how they were laid down.
        assert [aa for _, aa, _, _ in got] == [expected[c] for c in reversed(all_codons)], table
        # And the forward strand in one go.
        assert translate("".join(all_codons), table) == "".join(expected[c] for c in all_codons), table
        t = genetic_code.TABLES[table]
        assert sorted(t.stop_codons) == sorted(ref.stop_codons), table
        assert sorted(t.start_codons) == sorted(ref.start_codons), table
    # --- tables differ where they should
    assert translate("TGA", 1) == "*" and translate("TGA", 4) == "W"        # mycoplasma: TGA is Trp
    assert translate("AGA", 1) == "R" and translate("AGA", 2) == "*"        # vertebrate mitochondria: AGA/AGG stop
    assert translate("ATA", 1) == "I" and translate("ATA", 2) == "M"
    assert translate("CTG", 12) == "S" and translate("CTG", 1) == "L"       # alternative yeast nuclear
    assert translate("TGA", 25) == "G"                                       # Gracilibacteria
    # Table 11 has the standard amino acids and stops, but more start codons.
    assert genetic_code.TABLES[11].codons == genetic_code.TABLES[1].codons
    assert set(genetic_code.TABLES[1].start_codons) < set(genetic_code.TABLES[11].start_codons)
    # An unknown table falls back to the standard code rather than failing.
    assert translate("TGA", 99) == "*"
    # --- Rust and Python agree on the data
    from opendna_core import genetic_codes
    assert [t[0] for t in genetic_codes()] == IDS
    assert all(len(t[2]) == len(t[3]) == 64 for t in genetic_codes())
    # --- translation through the protein search (standard table) ---
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
    # --- ORFs: alternative start codons and stops follow the table ---
    gtg = Document.new_sequence("g", "CCGTGAAAAAATAACC", False)         # GTG AAA AAA TAA
    assert gtg.orfs(4) == []
    assert gtg.orfs(4, 11) == [], "alternative starts are off by default"
    assert gtg.orfs(4, 11, True) == [(2, 14, 1)]
    assert gtg.orfs(4, 1, True) == [], "table 1 has no GTG start"
    ttg = Document.new_sequence("t", "CCTTGAAAAAATAACC", False)
    assert ttg.orfs(4, 1, True) == [(2, 14, 1)] and ttg.orfs(4, 1) == []
    # The reverse strand too.
    assert Document.new_sequence("r", reverse_complement("CCGTGAAAAAATAACC"), False).orfs(4, 11, True) == [(2, 14, -1)]
    # --- stops follow the table
    seq = Document.new_sequence("m", "ATGTGAAAATAA", False)             # ATG TGA AAA TAA
    assert seq.orfs(4) == [] and seq.orfs(4, 4) == [(0, 12, 1)]          # table 4: TGA is not a stop
    assert Document.new_sequence("v", "ATGAGAAAATAA", False).orfs(1, 2) == [(0, 6, 1)]   # table 2: AGA stops
    # --- protein motif search uses the table
    seq = "AAATGAAAAAA"   # ... TGA AAA: W K in table 4, stop K in table 1 (frame +3 starts at base 2)
    assert find_protein(seq, "WK", 0, 100, 4)
    assert not find_protein(seq, "WK", 0, 100, 1)


STALE = "MQ"


def refresh_doc():
    """#0 stale, #1 codon_start 2 stale, #2 pseudo, #3 transl_except, #4 no /translation, #5 table 4, #6 partial without frame."""
    seq = "ATGAAATAA" + "GATGAAATAA" + "ATGAAATAA" + "ATGTGATAA" + "ATGAAATAA" + "ATGTGAAAATAA" + "ATGAAATAA"
    return cds_doc(seq,
                   cds("1..9", translation=STALE),
                   cds("10..19", label="cs2", codon_start="2", translation=STALE),
                   cds("20..28", label="pseudo", pseudo="", translation=STALE),
                   cds("29..37", label="exc", transl_except="(pos:32..34,aa:Sec)", translation="MQ"),
                   cds("38..46", label="none"),
                   cds("47..58", label="t4", transl_table="4", translation=STALE),
                   cds("<59..67", label="partial", translation=STALE))


def trans(ed, label):
    return dict(next(f for f in ed.features if f.label == label).qualifiers).get("translation")


def test_cds_translation_and_translation_refresh(code_setting):
    """A CDS is read with its own /transl_table (beating the setting), /codon_start and /transl_except (an unreproducible one
    is flagged); exports follow the setting; the /translation refresh updates stale ones, leaves the rest alone, is one undo
    step, uses the setting's table for a CDS without its own, and matches what the views draw."""
    # --- /transl_table regardless of the setting
    doc = cds_doc("ATGTGAGTGTAA", cds("1..12", transl_table="4"), cds("1..12", label="plain"))
    ed = Editor(doc, Path("t.gb"))
    own, plain = ed.features
    for setting in (1, 11):
        choose(setting)
        assert "".join(aa for _, aa, _, _ in cds_codons(doc, own)) == "MWV*"
    choose(4)
    assert "".join(aa for _, aa, _, _ in cds_codons(doc, plain)) == "MWV*", "no /transl_table: the setting applies"
    choose(1)
    assert "".join(aa for _, aa, _, _ in cds_codons(doc, plain)) == "M*V*"
    # --- /codon_start and the start codon as methionine
    doc = cds_doc("ATGTGAGTGTAA" + "TTGAAATAA",
                  cds("1..12", codon_start="2"), cds("13..21", label="ttg"), cds("<13..21", label="partial", codon_start="1"))
    two, ttg, partial = Editor(doc, Path("t.gb")).features
    assert "".join(aa for _, aa, _, _ in cds_codons(doc, two)) == "CEC"
    # A start codon of the table at the start of a complete CDS is shown (and exported) as M, like GenBank's /translation.
    assert "".join(aa for _, aa, _, _ in cds_codons(doc, ttg)) == "MK*"
    assert "".join(aa for _, aa, _, _ in cds_codons(doc, partial)) == "LK*"
    # --- /transl_except is applied and an unreproducible one is flagged
    seq = "ATGTGATAA"
    sec = cds_doc(seq, cds("1..9", transl_except="(pos:4..6,aa:Sec)"))
    assert "".join(aa for _, aa, _, _ in cds_codons(sec, Editor(sec, Path("t.gb")).features[0])) == "MU*"
    odd = cds_doc(seq, cds("1..9", transl_except="(pos:4..5,aa:TERM)"))
    f = Editor(odd, Path("t.gb")).features[0]
    assert "not reproducible" in genetic_code.cds_note(odd, f)
    assert "".join(aa for _, aa, _, _ in cds_codons(odd, f)) == "M**", "drawn plainly, and flagged"
    # --- exports use the chosen table ---
    import exports
    doc = cds_doc("CCATGTGAAAATAACC", cds("3..14"))
    ed = Editor(doc, Path("t.gb"))
    f = ed.features[0]
    choose(1)
    assert exports.protein_of(ed, f) == "M*K", "an internal stop is kept: only the final one is stripped"
    assert exports.translate("ATGTGA") == "M*"
    choose(4)
    assert exports.protein_of(ed, f) == "MWK"
    assert exports.translate("ATGTGA") == "MW"
    # --- /translation refresh ---
    # --- updates stale translations and leaves the rest alone
    choose(1)
    ed = Editor(refresh_doc(), Path("t.gb"))
    before = {f.label: f.qualifiers for f in ed.features}
    depth = len(ed.doc.history()[0])
    report = ed.refresh_translations()
    status = {ed.features[i].label: s for i, s, _ in report}
    assert status == {"g": "updated", "cs2": "updated", "pseudo": "pseudo", "exc": "updated", "t4": "updated", "partial": "partial"}
    assert trans(ed, "g") == "MK"
    assert trans(ed, "cs2") == "MK", "/codon_start=2 is honoured"
    assert trans(ed, "t4") == "MWK", "/transl_table=4 wins over the setting"
    assert trans(ed, "exc") == "MU", "the selenocysteine /transl_except is applied"
    assert trans(ed, "pseudo") == STALE and trans(ed, "partial") == STALE
    assert trans(ed, "none") is None, "a CDS without /translation never gets one"
    # Only /translation values changed.
    for f in ed.features:
        assert [kv for kv in f.qualifiers if kv[0] != "translation"] == [kv for kv in before[f.label] if kv[0] != "translation"]
    # One undo step; a second run does nothing and adds no step.
    assert len(ed.doc.history()[0]) == depth + 1
    assert all(s == "unchanged" or s in ("pseudo", "partial") for _, s, _ in ed.refresh_translations())
    assert len(ed.doc.history()[0]) == depth + 1
    ed.undo()
    assert {f.label: f.qualifiers for f in ed.features} == before
    # --- uses the settings table for a CDS without its own
    choose(4)
    doc = cds_doc("ATGTGAAAATAA", cds("1..12", translation="M"))
    ed = Editor(doc, Path("t.gb"))
    assert [s for _, s, _ in ed.refresh_translations()] == ["updated"]
    assert trans(ed, "g") == "MWK"
    choose(1)
    assert [s for _, s, why in ed.refresh_translations()] == ["other"], "the standard code: a stop inside the CDS, left alone"
    assert trans(ed, "g") == "MWK"
    # --- matches the translation drawn in the views: Rust (refresh) and Python (drawing, FASTA export) follow the same
    # plan: selenocysteine, codon_start, tables, strands
    import exports
    choose(11)
    seq = "ATGTGATAA"
    top = reverse_complement(seq)
    doc = cds_doc(top + "GATGAAATAA" + "TTGAAATAA",
                  cds("complement(1..9)", label="rev", transl_except="(pos:complement(4..6),aa:Sec)", translation="X"),
                  cds("10..19", label="cs", codon_start="2", translation="X"),
                  cds("20..28", label="ttg", translation="X"))
    ed = Editor(doc, Path("t.gb"))
    ed.refresh_translations()
    for f in ed.features:
        assert trans(ed, f.label) == exports.protein_of(ed, f), f.label
    assert trans(ed, "rev") == "MU" and trans(ed, "ttg") == "MK"
    # --- the sort of report text
    report = [(0, "updated", ""), (1, "updated", ""), (2, "updated", ""), (3, "pseudo", ""), (4, "unchanged", "")]
    assert genetic_code.summarize(report) == "Updated /translation of 3 CDSs, skipped 1 pseudo"
    assert genetic_code.summarize([(0, "updated", "")]) == "Updated /translation of 1 CDS"
