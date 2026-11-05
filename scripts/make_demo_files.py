"""Write the demo fragments used by docs/TUTORIAL.md into demo/ (random but seeded, so the same
files every time): a destination vector with a dropout between two outward-facing BsaI sites, a
promoter + reporter insert and a terminator. They are made-up sequences, not real parts.

    python scripts/make_demo_files.py [folder]
"""
import random
import sys
from pathlib import Path

from opendna_core import Document

rng = random.Random(5)


def dna(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def save(folder, name, bases, features, circular=False):
    doc = Document.new_sequence(name, bases, circular)
    doc.add_features([(kind, loc, [("label", label)]) for kind, loc, label in features])
    path = folder / f"{name}.gb"
    doc.save_genbank(path)
    print(f"{path}: {doc.length} bp, {len(features)} features")


def main(folder):
    folder.mkdir(parents=True, exist_ok=True)
    # Vector: backbone, a BsaI site facing left, the dropout, a BsaI site facing right, backbone.
    # The sites point outwards, so BsaI releases the dropout together with both sites and the
    # rest of the plasmid (1,700 bp) is the backbone the inserts are ligated into.
    save(folder, "demo_vector", dna(300) + "TGAGACC" + dna(120) + "GGTCTCA" + dna(1400), [
        ("rep_origin", "1..250", "ori"), ("CDS", "complement(700..1500)", "AmpR"),
        ("misc_feature", "308..427", "dropout")], circular=True)
    save(folder, "demo_insert", dna(700), [
        ("promoter", "1..80", "pDemo"), ("CDS", "101..700", "reporter")])
    save(folder, "demo_terminator", dna(220), [("terminator", "20..200", "T-demo")])


if __name__ == "__main__":
    main(Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent.parent / "demo")
