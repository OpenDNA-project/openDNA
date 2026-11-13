"""Protein properties (app/protein.py): mass, isoelectric point, charge, and the numbered layout of the letters."""
import os
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ["OPENDNA_DATA_DIR"] = tempfile.mkdtemp(prefix="opendna_test_")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "app"))

import pytest

import protein


def test_mass_pi_charge_and_numbered_layout():
    """Mass, isoelectric point, charge, the protein up to its stop, the numbered layout and the residue index <-> text offset
    round trip."""
    # --- Gly-Gly: 2 x 57.0519 + water; a lysine-only chain is basic, an aspartate-only chain acidic
    assert protein.mass("GG") == pytest.approx(132.1191, abs=1e-3)
    assert protein.isoelectric_point("KKKKKK") > 10 and protein.isoelectric_point("DDDDDD") < 4
    assert protein.charge("KKKKKK", 7.0) > 5 > 0 > protein.charge("DDDDDD", 7.0)
    assert protein.until_stop("MK*AA") == "MK"

    # --- the layout and the index round trip
    text = protein.layout("A" * 120)
    lines = text.split("\n")
    assert len(lines) == 3 and lines[0].startswith("     1  AAAAAAAAAA AAAAAAAAAA") and lines[1].lstrip().startswith("51")
    for i in (0, 9, 10, 49, 50, 119):
        line, off = protein.offset_of(i)
        assert protein.index_at(off, line) == i, i

    # --- every residue's mass against Biopython's ProtParam (independent table), alone and in one chain
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    all20 = "ACDEFGHIKLMNPQRSTVWY"
    assert protein.mass(all20) == pytest.approx(ProteinAnalysis(all20).molecular_weight(), abs=0.05)
    for a in all20:
        assert protein.mass(a) - 18.0153 == pytest.approx(ProteinAnalysis(a).molecular_weight() - 18.015, abs=0.01), a

    # --- charge at each side chain's own pKa, by hand from the EMBOSS iep table typed here (not read from the module):
    # the side chain is half charged there; the termini (N 8.6, C 3.6) add their own Henderson-Hasselbalch terms
    def n_term(ph):
        return 1 / (1 + 10 ** (ph - 8.6))

    def c_term(ph):
        return 1 / (1 + 10 ** (3.6 - ph))
    for aa, pka in (("K", 10.8), ("R", 12.5), ("H", 6.5)):
        assert protein.charge(aa, pka) == pytest.approx(n_term(pka) + 0.5 - c_term(pka), abs=1e-6), aa
    for aa, pka in (("D", 3.9), ("E", 4.1), ("C", 8.5), ("Y", 10.1)):
        assert protein.charge(aa, pka) == pytest.approx(n_term(pka) - 0.5 - c_term(pka), abs=1e-6), aa
    assert protein.isoelectric_point("G") == pytest.approx((8.6 + 3.6) / 2, abs=0.01)

    # --- offset_of / index_at against the layout text itself: the i-th letter of the text sits where offset_of(i) says,
    # and index_at on that letter gives i back (both wrong the same way no longer passes)
    assert [protein.offset_of(i) for i in (0, 10, 50)] == [(0, 8), (0, 19), (1, 8)]
    seen = 0
    for n, line in enumerate(lines):
        for off, ch in enumerate(line):
            if off >= 8 and ch == "A":
                assert protein.offset_of(seen) == (n, off) and protein.index_at(off, n) == seen, (n, off)
                seen += 1
    assert seen == 120
