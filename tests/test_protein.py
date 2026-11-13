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
