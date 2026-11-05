"""The real sequence files the tests use, from demo/: copied once into a temporary folder, so a test that saves (or a bug that
saves where it should not) can never change the files shipped with the program.

- PLASMID: pBAD30_eilR.gb, the main plasmid of the tests (araC / eilR, 6 primer sites, one across the origin).
- PLASMID2: pFG051_Chlor.gb, a second plasmid (15 primer sites, many misc_features).
- DNA_FILE: pFG051_Chlor.dna, the same plasmid as the original .dna file (not cleaned up): only for tests of reading .dna files.
- GENOME: BW25113_2026.gb, the E. coli genome (4.6 Mb).
"""
import shutil
import tempfile
from pathlib import Path

DEMO = Path(__file__).resolve().parent.parent / "demo"
_folder = Path(tempfile.mkdtemp(prefix="opendna_samples_"))


def _copy(name):
    source = DEMO / name
    target = _folder / name
    if source.exists():
        shutil.copyfile(source, target)
    return target   # missing source: a path that does not exist, so the tests that need it skip


PLASMID = _copy("pBAD30_eilR.gb")
PLASMID2 = _copy("pFG051_Chlor.gb")
DNA_FILE = _copy("pFG051_Chlor.dna")
GENOME = _copy("BW25113_2026.gb")
