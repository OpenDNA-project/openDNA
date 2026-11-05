"""Variants from reads written into the sequence or annotated (app/variants.py): the code behind the Alignment view's
right-click "Apply" and "Annotate". Applying changes the user's bases, so each kind of change is checked against a plain
string edit of the same sequence, together with the order of the edits, the DNA lock, and the `variation` features."""
import os
import sys
import tempfile
from pathlib import Path
from types import SimpleNamespace

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ["GENOMEEDIT_DATA_DIR"] = tempfile.mkdtemp(prefix="genomeedit_test_")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "app"))

from PySide6.QtWidgets import QApplication

from genomeedit_core import Document

from assembly_panel import FragmentEditor
from variants import annotate_variants, apply_variants, can_apply, change_text, variant_feature

app = QApplication.instance() or QApplication([])

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


def test_change_text_and_which_variants_can_be_applied():
    """The label of each kind of change, and an insertion whose bases the reads did not settle (empty or N) is not applied."""
    assert change_text(variant("snp", 4, "A", "G")) == "A > G"
    assert change_text(variant("del", 4, "ACGTA")) == "deletion of 5 bp (ACGTA)"
    assert change_text(variant("del", 4, "A" * 11)) == "deletion of 11 bp"   # more than 10 bases: not spelled out
    assert change_text(variant("ins", 4, alt="TT")) == "insertion of TT"
    assert change_text(variant("ins", 4, alt="")) == "insertion of ?"
    assert can_apply(variant("snp", 4, "A", "G")) and can_apply(variant("del", 4, "AC"))
    assert can_apply(variant("ins", 4, alt="TT"))
    assert not can_apply(variant("ins", 4, alt="")) and not can_apply(variant("ins", 4, alt="TNT"))


def test_each_kind_of_change_equals_the_string_edit():
    """A SNP replaces one base, a deletion removes its bases, an insertion goes before the base at its position; the bases
    equal the same edit made on a plain string, and the coordinates-moved message comes only after an indel."""
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


def test_several_changes_are_applied_last_position_first_and_each_undoes():
    """Changes are written right to left, so every position still means the base of the original sequence; the result equals
    the string edit in that order whatever order they are chosen in. Each change is one undo step; undoing them all gives the
    original back. Insertions with unsettled bases are left out of the count."""
    chosen = [variant("ins", 20, alt="TTT"), variant("snp", 2, "G", "C"), variant("del", 10, "GGCC"),
              variant("ins", 25, alt="NN")]
    expected = SEQ[:2] + "C" + SEQ[3:10] + SEQ[14:20] + "TTT" + SEQ[20:]
    for order in (chosen, chosen[::-1]):
        ed, _ = editor()
        assert apply_variants(ed, entry(), list(order)) == 3
        assert ed.doc.sequence() == expected
        undone = 0
        while ed.doc.sequence() != SEQ:
            assert ed.doc.undo(), "nothing left to undo before the original came back"
            undone += 1
        assert undone == 3 and len(ed.doc.history()[1]) == 3   # and three redo steps


def test_nothing_is_applied_while_the_dna_is_locked_or_nothing_can_be():
    """The DNA lock refuses every change (with its message) and leaves the entry alone; an empty or unusable choice is a no-op."""
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


def test_variation_features_cover_the_right_bases():
    """GenBank `variation` features: a SNP on its base, a deletion over the deleted bases with an empty /replace, an
    insertion over the two bases around it (clipped at the ends), with the read counts in the note."""
    kind, loc, quals = variant_feature(variant("snp", 5, "C", "T", count=7, depth=10), 32, "sample.bam")
    assert (kind, loc) == ("variation", "6")
    assert dict(quals) == {"label": "C > T", "replace": "T", "note": "7 of 10 reads (70 %) in sample.bam"}
    assert variant_feature(variant("del", 8, "AAGG"), 32)[1] == "9..12"
    assert dict(variant_feature(variant("del", 8, "AAGG"), 32)[2])["replace"] == ""
    assert variant_feature(variant("ins", 8, alt="CCC"), 32)[1] == "8..9"   # between base 8 and base 9 (1-based)
    assert variant_feature(variant("ins", 0, alt="C"), 32)[1] == "1..1"     # before the first base
    assert variant_feature(variant("ins", 32, alt="C"), 32)[1] == "32..32"  # after the last base
    assert dict(variant_feature(variant("snp", 1, "C", "A", count=0, depth=0), 32)[2])["note"] == "0 of 0 reads (0 %)"


def test_annotating_adds_one_feature_per_variant_in_one_undo_step_and_leaves_the_bases():
    ed, _ = editor()
    chosen = [variant("snp", 5, "C", "T"), variant("del", 8, "AAGG")]
    assert annotate_variants(ed, entry(), chosen) == 2
    features = [f for f in ed.doc.features() if f.kind == "variation"]
    assert sorted(f.location for f in features) == ["6", "9..12"]
    for f in features:   # the bases each feature covers are the reference bases of its variant
        if f.location == "6":
            assert ed.doc.sequence(5, 6) == "C"
        else:
            assert ed.doc.sequence(8, 12) == "AAGG"
    assert ed.doc.sequence() == SEQ
    assert ed.doc.undo()
    assert not [f for f in ed.doc.features() if f.kind == "variation"]
    assert annotate_variants(ed, entry(), []) == 0
