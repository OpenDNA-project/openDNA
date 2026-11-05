"""Tm against Biopython's nearest-neighbour implementation (same model), and primer
binding against the sites stored in the real .dna file (skipped if absent)."""
import random
from pathlib import Path

import pytest
from Bio.SeqUtils import MeltingTemp as mt

from opendna_core import Document, tm

ROOT = Path(__file__).resolve().parent.parent
from sample_files import DNA_FILE   # a copy of the real .dna file (see sample_files.py)


def biopython_tm(seq, na, mg, dntp, primer_nm):
    # dnac1 = dnac2 = C/2 makes Biopython use C/4, Primer3's convention for our formula.
    return mt.Tm_NN(seq, nn_table=mt.DNA_NN3, saltcorr=5, Na=na, Mg=mg, dNTPs=dntp,
                    dnac1=primer_nm / 2, dnac2=primer_nm / 2)


def test_tm_matches_biopython():
    for conditions in [(50, 1.5, 0.6, 50), (50, 0, 0, 250), (100, 2, 0.8, 500)]:
        rng = random.Random(1)
        for _ in range(200):
            seq = "".join(rng.choice("ACGT") for _ in range(rng.randint(12, 40)))
            if seq == seq[::-1].translate(str.maketrans("ACGT", "TGCA")):
                continue  # palindromes: Biopython and Primer3 use different concentrations
            assert tm(seq, *conditions) == pytest.approx(biopython_tm(seq, *conditions), abs=0.01), (conditions, seq)
    assert tm("ACGTNACGT") is None


def test_tm_of_self_complementary_primers_matches_biopython():
    """The branch the test above leaves out: a primer that is its own reverse complement pairs with itself, so the duplex
    uses the whole primer concentration (not C/4) and a symmetry term of -1.4 e.u. (SantaLucia 1998). Biopython's selfcomp
    mode is the same model; treating these primers like the others would put Tm several degrees off."""
    palindromes = ["GAATTCGAATTC", "ACGTACGTACGT", "GGATCCGGATCC", "CGCGAATTCGCG", "AAGCTTAAGCTT", "GCGGCCGCGCGGCCGC"]
    for seq in palindromes:
        assert seq == seq[::-1].translate(str.maketrans("ACGT", "TGCA")), seq
        for na, mg, dntp, primer_nm in [(50, 1.5, 0.6, 50), (50, 0, 0, 250), (100, 2, 0.8, 500)]:
            expected = mt.Tm_NN(seq, nn_table=mt.DNA_NN3, saltcorr=5, Na=na, Mg=mg, dNTPs=dntp, dnac1=primer_nm, dnac2=0,
                                selfcomp=True)
            assert tm(seq, na, mg, dntp, primer_nm) == pytest.approx(expected, abs=0.01), (seq, na, mg, dntp, primer_nm)
            # ... and visibly different from the non-self-complementary formula, so the branch matters
            assert abs(expected - biopython_tm(seq, na, mg, dntp, primer_nm)) > 1, seq


def test_every_primer_the_dna_file_annotates_is_found_on_its_strand_with_the_same_3_prime_end():
    """The other direction of test_primer_binding.py::test_our_sites_are_the_binding_sites_stored_in_each_real_dna_file (every
    site we find is one the file stores): every primer_bind the real .dna file stores is found by us, on the same strand, inside
    the stored site and ending at the same 3' base. The stored site can be longer: SnapGene also marks a 5' tail that matches
    by chance (oAC812, oAC813), which does not anneal from the 3' end and which we leave out."""
    if not DNA_FILE.exists():
        pytest.skip(f"missing {DNA_FILE.name}")
    doc = Document(str(DNA_FILE))
    n = doc.length
    primers = dict(doc.primers())
    stored = [(f.label, f.spans, f.strand) for f in doc.features() if f.kind == "primer_bind"]
    assert len(stored) == 15 and len(primers) == 15, "the file's 15 primers, each with a stored site"

    def bases(spans):
        return {p % n for a, b in spans for p in range(a, b)}

    shorter = []
    for name, spans, strand in stored:
        found = [(a, b, s) for a, b, s, _ in doc.binding_sites(primers[name], 0) if s == strand]
        three = max(bases(spans), key=lambda p: (p - spans[0][0]) % n) if strand > 0 else spans[0][0]
        ours = [(a, b) for a, b, _ in found if bases([(a, b)]) <= bases(spans)
                and three == ((b - 1) % n if strand > 0 else a % n)]
        assert ours, f"{name}: stored at {spans} on strand {strand}, found {found}"
        if len(bases([ours[0]])) < len(bases(spans)):
            shorter.append(name)
    assert sorted(shorter) == ["oAC812", "oAC813"], shorter
