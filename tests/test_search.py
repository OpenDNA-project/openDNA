"""The search engine: DNA with IUPAC / mismatches, proteins, names, primers, enzymes, scope and FASTA files."""
import os
import random
import sys
import tempfile
import threading
from pathlib import Path

import pytest

os.environ["OPENDNA_DATA_DIR"] = tempfile.mkdtemp(prefix="opendna_test_")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "app"))


from opendna_core import Document

import primer_db
import search
from assembly_panel import FragmentEditor
from editor import reverse_complement
from search import Options

rng = random.Random(5)
SEQ = "".join(rng.choice("ACGT") for _ in range(2000))
# Clean the random part of the sites used below so they occur exactly where put.
SEQ = SEQ.replace("GAATTC", "GAATTG").replace("GGATCC", "GGATCG")
SEQ = SEQ[:300] + "GAATTC" + SEQ[306:700] + "ACGTACGTAC" + SEQ[710:1200] + reverse_complement("GGATCCAAAT") + SEQ[1210:]


@pytest.fixture(autouse=True)
def private_primer_db(tmp_path, monkeypatch):
    monkeypatch.setattr(primer_db, "_db", primer_db.PrimerDB(tmp_path / "primers.tsv"))


def editor(seq=SEQ, circular=False):
    return FragmentEditor(Document.new_sequence("t", seq, circular), Path("t.gb"))


def spans(hits):
    return [(h.start, h.end) for h in hits]



def test_sequence_search_engine(tmp_path):
    """The search engine: DNA on both strands with IUPAC, mismatches and their limit, a site across a circle's origin, protein
    motifs on both strands; names, details, primers and enzymes found by name only with their option; the selection scope and
    the search in a FASTA file."""
    # --- DNA: both strands and IUPAC ---
    ed = editor()
    hits, _ = search.search_sequence(ed, "GGATCCAAAT", Options())
    assert [(h.start, h.end, h.strand) for h in hits] == [(1200, 1210, -1)]
    hits, _ = search.search_sequence(ed, "GGATCCAAAT", Options(both=False))
    assert hits == []
    # R = A or G: "GAATTC" is "GRATTC"; with IUPAC off the R is not a base, so nothing matches (and a name search finds nothing)
    assert (300, 306) in spans(search.search_sequence(ed, "GRATTC", Options())[0])
    assert search.search_sequence(ed, "GRATTC", Options(iupac=False))[0] == []
    # N in the query matches anything
    assert (700, 710) in spans(search.search_sequence(ed, "ACGTNCGTAC", Options())[0])
    # --- fuzzy matching and its limit ---
    ed = editor()
    probe = SEQ[700:710]
    mutated = probe[:4] + ("A" if probe[4] != "A" else "C") + probe[5:]
    assert (700, 710) not in spans(search.search_sequence(ed, mutated, Options())[0])
    hits, _ = search.search_sequence(ed, mutated, Options(fuzzy=1))
    assert any((h.start, h.end, h.mismatches) == (700, 710, 1) for h in hits)
    # never more than a third of the query, whatever the setting: 1 mismatch in 4 bases, 3 in 30
    assert search.allowed_mismatches(Options(fuzzy=3), 4) == 1 and search.allowed_mismatches(Options(fuzzy=3), 30) == 3
    # --- circular site across the origin ---
    ed = editor("TTCAAAAAAAAAAAAAAGA", circular=True)   # GATTC across the origin: GA at the end, TTC at the start
    hits, _ = search.search_sequence(ed, "GATTC", Options(both=False))
    assert spans(hits) == [(17, 22)]   # ends past the length: it wraps
    ed.set_search_hits([(17, 22)])
    assert ed.search_spans == [(0, 3), (17, 19)]   # split at the origin for the views
    assert ed.hits_between(0, 5) == [(0, 3)] and ed.hits_between(10, 16) == []
    # --- protein motif in both strands ---
    ed = editor("ATG" + "GCTAAATGGCAT" + "TAA" + "AAAA")    # M A K W H *
    hits, _ = search.search_sequence(ed, "AKW", Options(translations=True))
    assert [(h.start, h.end, h.kind, h.strand) for h in hits] == [(3, 12, "protein", 1)]
    assert search.search_sequence(ed, "AKW", Options())[0] == []     # translations off: AKW is not DNA, so a name search: nothing
    rc = editor(reverse_complement(ed.doc.sequence()))
    hits, _ = search.search_sequence(rc, "AXW", Options(translations=True))      # X = any
    # --- names, details, primers, enzymes ---
    ed = editor()
    ed.add_features([("CDS", "101..400", [("gene", "abcD"), ("product", "kinase"), ("db_xref", "needle")])])
    assert spans(search.search_sequence(ed, "abcd", Options())[0]) == [(100, 400)]
    assert search.search_sequence(ed, "needle", Options())[0] == []                      # a note is not a name...
    assert spans(search.search_sequence(ed, "needle", Options(details=True))[0]) == [(100, 400)]  # ...unless details are on
    # primers: the name finds where it binds (on by default)
    primer_db.add_primer(ed, "fwdprimer", SEQ[500:520])
    assert spans(search.search_sequence(ed, "fwdprim", Options())[0]) == [(500, 520)]
    assert search.search_sequence(ed, "fwdprim", Options(primers=False))[0] == []
    # enzymes: EcoRI finds its one site, nothing found without the option (the word is not DNA either)
    hits, _ = search.search_sequence(ed, "EcoRI", Options(enzymes=True))
    assert [(h.start, h.end, h.kind) for h in hits] == [(300, 306, "enzyme")]
    assert search.search_sequence(ed, "EcoRI", Options())[0] == []
    # --- selection scope ---
    ed = editor()
    ed.select(250, 350)
    hits, note = search.search_sequence(ed, "GAATTC", Options(scope="selection"))
    assert spans(hits) == [(300, 306)] and not note
    ed.select(400, 500)
    assert search.search_sequence(ed, "GAATTC", Options(scope="selection"))[0] == []
    ed.set_cursor(10)
    hits, note = search.search_sequence(ed, "GAATTC", Options(scope="selection"))
    assert spans(hits) == [(300, 306)] and "Nothing is selected" in note
    # --- file scope ---
    a = "ACGT" * 50 + "GGTACC" + "TTTT" * 20
    b = "TTTT" * 30 + "GGTACC" + "GGTACC"[::-1] + "ACAC" * 10
    path = tmp_path / "g.fa"
    path.write_text(f">chrA\n{a}\n>chrB\n{b}\n")
    got = {}
    motif = search.search_file(path, ["chrA", "chrB"], "chrA", a, "GGTACC", Options(), threading.Event(),
                               lambda name, hits: got.__setitem__(name, spans(hits)))
    assert motif == "GGTACC"
    assert got["chrA"] == [(200, 206)] and got["chrB"] == [(120, 126)]   # GGTACC is a palindrome: one hit per site
    assert search.search_file(path, ["chrA"], "chrA", a, "kinase", Options(), threading.Event(), lambda *_: None) is None
    # --- the editor keeps the strand each hit matched (0: no strand given, or a palindrome found on both) ---
    ed = editor("ACGT" * 30)
    ed.set_search_hits([(10, 20, -1), (40, 50)])
    assert ed.search_strand[(10, 20)] == -1 and ed.search_strand[(40, 50)] == 0
    ed.set_search_hits([(5, 15, 1), (5, 15, -1)])
    assert ed.search_strand[(5, 15)] == 0
    ed.set_search_hits([])
    assert ed.search_strand == {} and ed.search_spans == []
