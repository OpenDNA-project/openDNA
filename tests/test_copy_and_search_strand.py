"""Copy / Copy Reverse Complement keep the annotations of the selection (flipped for the reverse complement), and a search
hit is highlighted on the strand it matched only."""
import os
import random
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ["OPENDNA_DATA_DIR"] = tempfile.mkdtemp(prefix="opendna_test_")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "app"))

from PySide6.QtGui import QGuiApplication
from PySide6.QtWidgets import QApplication

from opendna_core import Document

from assembly_panel import FragmentEditor
from editor import reverse_complement
from sequence_view import SequenceView

app = QApplication.instance() or QApplication([])
rng = random.Random(9)
SEQ = "".join(rng.choice("ACGT") for _ in range(300))


def make(circular=False):
    doc = Document.new_sequence("t", SEQ, circular)
    doc.add_features([("gene", "61..100", [("label", "top")]),
                      ("gene", "complement(111..140)", [("label", "bottom")]),
                      ("CDS", "join(51..60,121..130)", [("label", "split")])])
    return FragmentEditor(doc, Path("t.gb"))


def own_bases(ed, label):
    """The bases a feature covers, read along its strand."""
    f = next(f for f in ed.doc.copy_features(0, ed.length)
             if any(q[0] == "label" and q[1] == label for q in f.qualifiers))
    return f


def label_of(f):
    return next(v for k, v in f.qualifiers if k == "label")


def features_by_label(ed):
    return {label_of(f): f for f in ed.doc.copy_features(0, ed.length) if f.kind != "source"}


def test_copy_and_copy_reverse_complement_keep_annotations():
    """Copy Reverse Complement keeps the annotations flipped and pastes them, a plain copy keeps them unflipped, the
    reverse complement works across the origin and twice gives the original, and the selection menu offers both."""
    # --- reverse complement keeps the annotations flipped and pastes them ---
    ed = make()
    ed.select(40, 160)
    mime = ed.copy_reverse_complement_mime()
    assert mime.text() == reverse_complement(SEQ[40:160])
    target = FragmentEditor(Document.new_sequence("u", "A" * 50, False), Path("u.gb"))
    target.set_cursor(10)
    assert target.paste_mime(mime)
    assert target.doc.sequence()[10:130] == reverse_complement(SEQ[40:160])
    pasted = features_by_label(target)
    assert set(pasted) == {"top", "bottom", "split"}
    # the bases of each feature, read along its own strand, are what they were in the original
    for label, strand in (("top", -1), ("bottom", 1), ("split", -1)):
        assert pasted[label].strand == strand, f"{label}: the strand is flipped"
    # "top" covered original 61..100: in the copy of 41..160 reversed it ends at 160 - 61 + 1 = 100 -> 10 + ... shifted by 10
    assert pasted["top"].location == "complement(%d..%d)" % (10 + 160 - 100 + 1, 10 + 160 - 61 + 1)
    # a join keeps both parts, in the new reading order
    assert pasted["split"].location.startswith("complement(join(")

    # --- the plain copy still keeps the annotations unflipped ---
    ed = make()
    ed.select(40, 160)
    mime = ed.copy_mime()
    target = FragmentEditor(Document.new_sequence("u", "A" * 50, False), Path("u.gb"))
    target.set_cursor(0)
    assert target.paste_mime(mime)
    pasted = features_by_label(target)
    assert pasted["top"].strand == 1 and pasted["bottom"].strand == -1

    # --- across the origin, and twice is the original ---
    ed = make(circular=True)
    ed.select(280, 320)   # 280..299 and 0..19
    mime = ed.copy_reverse_complement_mime()
    assert mime.text() == reverse_complement(SEQ[280:] + SEQ[:20])
    # reverse-complementing the reverse complement gives back the plain copy, features included
    middle = FragmentEditor(Document.new_sequence("m", "A" * 10, False), Path("m.gb"))
    middle.set_cursor(0)
    middle.paste_mime(mime)
    middle.select(0, 40)
    again = middle.copy_reverse_complement_mime()
    assert again.text() == SEQ[280:] + SEQ[:20]

    # --- the selection menu offers both copies with their shortcuts ---
    from feature_table import selection_menu
    from PySide6.QtCore import QPoint, QTimer
    ed = make()
    ed.select(40, 80)
    texts = []

    def grab():
        menu = QApplication.activePopupWidget()
        texts.extend(a.text() for a in menu.actions())
        menu.close()
    QTimer.singleShot(50, grab)

    view = SequenceView(ed)   # the menu's parent: it looks up the window it belongs to
    view.show()
    selection_menu(view, ed, QPoint(5, 5))
    view.close()
    assert "Copy	Ctrl+C" in texts and "Copy Reverse Complement	Ctrl+Shift+C" in texts
    QGuiApplication.clipboard().clear()


def _changed_rows(view, base):
    """The y of every pixel row that differs from the picture `base` (taken with no search hit shown)."""
    img = view.viewport().grab().toImage()
    return {y for y in range(img.height()) if any(img.pixel(x, y) != base.pixel(x, y) for x in range(img.width()))}


def test_search_hit_strand():
    """A search hit colours only the strand it matched (none given: both), and the editor keeps the matched strand of
    each hit (0 for no strand or a palindrome found on both strands)."""
    # --- a hit colours only the strand it matched ---
    ed = make()
    view = SequenceView(ed)
    view.resize(900, 400)
    view.show()
    app.processEvents()
    base = view.viewport().grab().toImage()
    ed.set_search_hits([(10, 20, 1)])
    app.processEvents()
    top = _changed_rows(view, base)
    ed.set_search_hits([(10, 20, -1)])
    app.processEvents()
    bottom = _changed_rows(view, base)
    ed.set_search_hits([(10, 20)])
    app.processEvents()
    both = _changed_rows(view, base)
    ed.set_search_hits([])
    app.processEvents()
    assert top and bottom and both
    assert max(top) < min(bottom), "a top-strand hit colours the top letters only, a bottom one the complement only"
    assert both == top | bottom, "no strand: both rows"
    view.close()

    # --- the search bar passes the matched strand to the view ---
    ed = make()
    ed.set_search_hits([(10, 20, -1), (40, 50)])
    assert ed.search_strand[(10, 20)] == -1 and ed.search_strand[(40, 50)] == 0
    ed.set_search_hits([(5, 15, 1), (5, 15, -1)])   # a palindrome: found on both strands
    assert ed.search_strand[(5, 15)] == 0
    ed.set_search_hits([])
    assert ed.search_strand == {} and ed.search_spans == []
