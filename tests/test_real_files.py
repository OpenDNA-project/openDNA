"""Checks against the real plasmid files (copies of demo/, see sample_files.py): the .dna reader against the GenBank export of
the same plasmid, and an edit saved as GenBank read back by Biopython. The genome round trip is in test_roundtrip.py."""
import time
from pathlib import Path

import pytest
from Bio import SeqIO

from genomeedit_core import Document

ROOT = Path(__file__).resolve().parent.parent
from sample_files import PLASMID as PLASMID_GB, PLASMID2, DNA_FILE as PLASMID_DNA   # copies (see sample_files.py)


def need(*paths):
    missing = [p.name for p in paths if not p.exists()]
    if missing:
        pytest.skip(f"missing {missing}")


def timed_open(path):
    t = time.perf_counter()
    doc = Document(path)
    print(f"\n  {path.name}: {doc!r} opened in {time.perf_counter() - t:.2f} s")
    return doc


def test_plasmid_files(tmp_path):
    """The .dna file and the GenBank file of the same plasmid hold the same sequence (the GenBank one was cleaned up, so its
    features differ); the other GenBank plasmid opens; phase 1 milestone: open .dna, insert bases inside a feature, save .gb,
    Biopython agrees."""
    need(PLASMID_DNA, PLASMID2, PLASMID_GB)
    # --- the .dna file matches its GenBank export
    dna, gb = timed_open(PLASMID_DNA), timed_open(PLASMID2)
    assert dna.length == gb.length == 3780 and dna.circular and gb.circular
    assert dna.sequence().upper() == gb.sequence().upper()
    assert len(dna.primers()) == 15, "the .dna file's primer list"

    # --- the other GenBank plasmid opens
    plasmid = timed_open(PLASMID_GB)
    assert plasmid.circular and plasmid.length == 5936

    # --- edit the .dna, save as GenBank, Biopython agrees
    doc = timed_open(PLASMID_DNA)
    print("  primers:", doc.primers())
    target = next(f for f in doc.features() if f.kind != "source" and len(f.spans) == 1
                  and f.spans[0][1] - f.spans[0][0] > 20)
    start, end = target.spans[0]
    doc.insert(start + 10, "GATTACA")
    moved = doc.features()[target.index]
    assert moved.spans == [(start, end + 7)]
    assert doc.modified

    out = tmp_path / "edited.gb"
    doc.save_genbank(out)
    assert not doc.modified

    rec = SeqIO.read(out, "genbank")
    assert str(rec.seq) == doc.sequence().upper()  # Biopython uppercases on read; we keep case
    assert len(rec.features) == len(doc.features())
    bio = rec.features[target.index]
    assert (int(bio.location.start), int(bio.location.end)) == (start, end + 7)

    assert doc.undo()
    assert doc.features()[target.index].spans == [(start, end)]
