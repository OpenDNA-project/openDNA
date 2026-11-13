"""Variants from reads written into the sequence or annotated (app/variants.py): the code behind the Alignment view's
right-click "Apply" and "Annotate". Applying changes the user's bases, so each kind of change is checked against a plain
string edit of the same sequence, together with the order of the edits, the DNA lock, and the `variation` features."""
import os
import sys
import tempfile
from pathlib import Path
from types import SimpleNamespace

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ["OPENDNA_DATA_DIR"] = tempfile.mkdtemp(prefix="opendna_test_")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "app"))

from opendna_core import Document

from assembly_panel import FragmentEditor
from variants import annotate_variants, apply_variants, can_apply, change_text, variant_feature

SEQ = "ACGTACGTAAGGCCTTACGTACGTAAGGCCTT"   # 32 bases


def variant(kind, pos, ref="", alt="", count=8, depth=10):
    """A variant as alignments.py hands it over (pos 0-based: a SNP's base, a deletion's first base, the base an insertion goes before)."""
    return SimpleNamespace(kind=kind, pos=pos, ref=ref, alt=alt, count=count, depth=depth)


def editor(seq=SEQ):
    ed = FragmentEditor(Document.new_sequence("v", seq, False), Path("v.gb"))
    messages = []
    ed.error.connect(messages.append)
    return ed, messages


def entry():
    """The BAM entry whose variant list and cache apply_variants resets."""
    return SimpleNamespace(name="sample.bam", variants=["old"], cache="old")


def test_variants_are_written_into_the_bases_or_annotated_beside_them():
    """Each kind of change is named, only a settled insertion can be applied, and the bases that come out equal the same
    edit made on a plain string — one undo step each, written last position first so every position still means a base of
    the original, with the coordinates-moved message only after an indel; the DNA lock refuses everything. Annotating
    instead leaves the bases alone and writes one `variation` feature per variant, over the bases it is about."""
    # --- the label of each change, and what can be applied ---
    assert change_text(variant("snp", 4, "A", "G")) == "A > G"
    assert change_text(variant("del", 4, "ACGTA")) == "deletion of 5 bp (ACGTA)"
    assert change_text(variant("del", 4, "A" * 11)) == "deletion of 11 bp"   # more than 10 bases: not spelled out
    assert change_text(variant("ins", 4, alt="TT")) == "insertion of TT"
    assert change_text(variant("ins", 4, alt="")) == "insertion of ?"
    assert can_apply(variant("snp", 4, "A", "G")) and can_apply(variant("del", 4, "AC")) and can_apply(variant("ins", 4, alt="TT"))
    assert not can_apply(variant("ins", 4, alt="")) and not can_apply(variant("ins", 4, alt="TNT"))

    # --- each kind of change equals the string edit ---
    cases = [
        (variant("snp", 5, "C", "T"), SEQ[:5] + "T" + SEQ[6:], False),
        (variant("del", 8, "AAGG"), SEQ[:8] + SEQ[12:], True),
        (variant("ins", 8, alt="CCC"), SEQ[:8] + "CCC" + SEQ[8:], True),
        (variant("snp", 0, "A", "G"), "G" + SEQ[1:], False),                       # first base
        (variant("del", 31, "T"), SEQ[:31], True),                                 # last base
        (variant("ins", 32, alt="GG"), SEQ + "GG", True),                          # after the last base
    ]
    for v, expected, moved in cases:
        ed, messages = editor()
        e = entry()
        assert apply_variants(ed, e, [v]) == 1, v
        assert ed.doc.sequence() == expected, v
        assert len(ed.doc.history()[0]) == 1, v   # one undo step
        # a SNP keeps the coordinates (the scan runs again: variants None); an indel moves them (no flags until re-mapped)
        assert e.variants == ([] if moved else None) and e.cache is None, v
        assert any("map the reads again" in m for m in messages) == moved, (v, messages)

    # --- several changes: written right to left, whatever order they are chosen in; each one undoes ---
    chosen = [variant("ins", 20, alt="TTT"), variant("snp", 2, "G", "C"), variant("del", 10, "GGCC"),
              variant("ins", 25, alt="NN")]
    expected = SEQ[:2] + "C" + SEQ[3:10] + SEQ[14:20] + "TTT" + SEQ[20:]
    for order in (chosen, chosen[::-1]):
        ed, _ = editor()
        assert apply_variants(ed, entry(), list(order)) == 3, "the unsettled insertion is left out of the count"
        assert ed.doc.sequence() == expected
        undone = 0
        while ed.doc.sequence() != SEQ:
            assert ed.doc.undo(), "nothing left to undo before the original came back"
            undone += 1
        assert undone == 3 and len(ed.doc.history()[1]) == 3   # and three redo steps

    # --- nothing is applied while the DNA is locked, or when nothing can be ---
    ed, messages = editor()
    ed.dna_locked = True
    e = entry()
    assert apply_variants(ed, e, [variant("snp", 5, "C", "T")]) == 0
    assert ed.doc.sequence() == SEQ and ed.doc.history() == ([], [])   # no undo step either
    assert e.variants == ["old"] and e.cache == "old"
    assert any("locked" in m for m in messages)
    ed.dna_locked = False
    assert apply_variants(ed, e, []) == 0
    assert apply_variants(ed, e, [variant("ins", 3, alt="")]) == 0
    assert ed.doc.sequence() == SEQ and e.variants == ["old"]

    # --- annotated instead: what each `variation` feature covers ---
    kind, loc, quals = variant_feature(variant("snp", 5, "C", "T", count=7, depth=10), 32, "sample.bam")
    assert (kind, loc) == ("variation", "6")
    assert dict(quals) == {"label": "C > T", "replace": "T", "note": "7 of 10 reads (70 %) in sample.bam"}
    assert variant_feature(variant("del", 8, "AAGG"), 32)[1] == "9..12"
    assert dict(variant_feature(variant("del", 8, "AAGG"), 32)[2])["replace"] == ""
    assert variant_feature(variant("ins", 8, alt="CCC"), 32)[1] == "8..9"   # between base 8 and base 9 (1-based)
    assert variant_feature(variant("ins", 0, alt="C"), 32)[1] == "1..1"     # before the first base
    assert variant_feature(variant("ins", 32, alt="C"), 32)[1] == "32..32"  # after the last base
    assert dict(variant_feature(variant("snp", 1, "C", "A", count=0, depth=0), 32)[2])["note"] == "0 of 0 reads (0 %)"

    # --- annotating: one feature per variant, one undo step, the bases untouched ---
    ed, _ = editor()
    assert annotate_variants(ed, entry(), [variant("snp", 5, "C", "T"), variant("del", 8, "AAGG")]) == 2
    features = [f for f in ed.doc.features() if f.kind == "variation"]
    assert sorted(f.location for f in features) == ["6", "9..12"]
    assert ed.doc.sequence(5, 6) == "C" and ed.doc.sequence(8, 12) == "AAGG", "each feature covers its variant's reference bases"
    assert ed.doc.sequence() == SEQ
    assert ed.doc.undo()
    assert not [f for f in ed.doc.features() if f.kind == "variation"]
    assert annotate_variants(ed, entry(), []) == 0
