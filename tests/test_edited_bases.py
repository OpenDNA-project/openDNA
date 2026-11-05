"""Edited bases (typed, pasted, replaced, inserted since the last save or load) are painted as blue letters with a thin blue underline in the sequence view (never a background: yellow is a search hit):
the ranges the document keeps, how they follow later edits, undo / redo, saves, and the painting itself."""
import os
import random
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ["GENOMEEDIT_DATA_DIR"] = tempfile.mkdtemp(prefix="genomeedit_test_")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "app"))

from PySide6.QtGui import QColor
from PySide6.QtWidgets import QApplication

from genomeedit_core import Document

import theme
from assembly_panel import FragmentEditor
from sequence_view import SequenceView

app = QApplication.instance() or QApplication([])
rng = random.Random(5)
SEQ = "".join(rng.choice("ACGT") for _ in range(200))


def make(circular=False):
    return FragmentEditor(Document.new_sequence("t", SEQ, circular), Path("t.gb"))


def test_ranges_when_opened_typed_replaced_pasted_deleted_and_shifted_by_later_edits():
    """Nothing is edited when opened; a typed run is one stretch, replace and paste mark their bases; deleting marks nothing and
    later edits shift (or shrink) the ranges."""
    # --- nothing is edited when opened ---
    assert make().doc.edited_ranges() == []

    # --- insert, replace and a typed run ---
    ed = make()
    ed.set_cursor(10)
    ed.type_bases("A")
    ed.type_bases("C")
    ed.type_bases("G")
    assert ed.doc.edited_ranges() == [(10, 13)], "typing a run is one stretch"
    ed.select(50, 54)
    ed.type_bases("TT")           # replace four bases by two
    assert ed.doc.edited_ranges() == [(10, 13), (50, 52)]
    ed.select(100, 100)
    ed.paste("GGGGG")
    assert ed.doc.edited_ranges()[-1] == (100, 105)

    # --- deleting leaves nothing and later edits shift the ranges ---
    ed = make()
    ed.set_cursor(100)
    ed.type_bases("AAAA")
    assert ed.doc.edited_ranges() == [(100, 104)]
    ed.select(0, 10)
    ed.backspace()                # a deletion before: the stretch moves left, nothing new is marked
    assert ed.doc.edited_ranges() == [(90, 94)]
    ed.set_cursor(0)
    ed.type_bases("CC")           # an insertion before: it moves right
    assert ed.doc.edited_ranges() == [(0, 2), (92, 96)]
    ed.select(93, 95)
    ed.backspace()                # a deletion inside shrinks it
    assert ed.doc.edited_ranges() == [(0, 2), (92, 94)]


def test_undo_redo_save_and_the_origin_keep_the_ranges_consistent(tmp_path):
    """Undo / redo keep the ranges consistent; a save clears them but autosave does not; an edit across the circular origin
    survives undo, redo and rotation."""
    # --- undo and redo ---
    ed = make()
    ed.set_cursor(20)
    ed.type_bases("GG")
    ed.select(5, 8)
    ed.backspace()
    assert ed.doc.edited_ranges() == [(17, 19)]
    ed.undo()                     # the deleted bases come back: the stretch is where it was
    assert ed.doc.edited_ranges() == [(20, 22)]
    ed.undo()                     # undoing an insertion removes its highlight
    assert ed.doc.edited_ranges() == []
    ed.redo()
    assert ed.doc.edited_ranges() == [(20, 22)]
    ed.redo()
    assert ed.doc.edited_ranges() == [(17, 19)]

    # --- save clears and autosave does not ---
    ed = make()
    ed.set_cursor(30)
    ed.type_bases("TTT")
    ed.doc.write_genbank(str(tmp_path / "auto.gb"))        # what autosave and the project writer use
    assert ed.doc.edited_ranges() == [(30, 33)], "autosave must not count as a save"
    ed.save(tmp_path / "t.gb")
    assert ed.doc.edited_ranges() == []
    ed.set_cursor(0)
    ed.type_bases("A")
    assert ed.doc.edited_ranges() == [(0, 1)]
    ed.undo()
    ed.undo()                     # back before the save: the removed bases are gone, nothing to paint
    assert ed.doc.edited_ranges() == []
    ed.redo()
    assert ed.doc.edited_ranges() == [], "redo back onto the saved state: nothing differs from the file"
    ed.doc.modified = False       # FASTA and project saves mark the document saved this way
    assert ed.doc.edited_ranges() == []

    # --- circular edit across the origin ---
    ed = make(circular=True)
    ed.select(190, 205)           # 190..199 and 0..4: passes the origin
    ed.type_bases("GGGGGG")
    ranges = ed.doc.edited_ranges()
    assert sum(b - a for a, b in ranges) == 6
    ed.undo()
    assert ed.doc.edited_ranges() == []
    ed.redo()
    assert ed.doc.edited_ranges() == ranges
    ed.set_cursor(0)
    ed.type_bases("A")            # rotation keeps them too
    ed.doc.rotate(50)
    assert sum(b - a for a, b in ed.doc.edited_ranges()) == 7
    ed.doc.undo()
    ed.doc.undo()
    ed.doc.undo()
    assert ed.doc.edited_ranges() == []


def _count(view, color):
    img = view.viewport().grab().toImage()
    want = QColor(color)
    return sum(1 for x in range(img.width()) for y in range(img.height()) if QColor(img.pixel(x, y)) == want)


def _yellow(view):
    """Pixels that read as the yellow of a search hit (translucent gold over white or over the dark background)."""
    img = view.viewport().grab().toImage()
    n = 0
    for x in range(img.width()):
        for y in range(img.height()):
            c = QColor(img.pixel(x, y))
            n += c.red() > 120 and c.green() > 100 and c.blue() < 0.5 * c.red()
    return n


def test_the_view_paints_edited_bases_blue_in_both_themes_and_a_search_hit_stays_a_hit(monkeypatch):
    for dark in (False, True):
        monkeypatch.setattr(theme, "is_dark", lambda dark=dark: dark)
        ed = make()
        view = SequenceView(ed)
        view.resize(900, 400)
        view.show()
        app.processEvents()
        colour = theme.hex("edited_fg")
        assert _count(view, colour) == 0
        ed.set_cursor(10)
        ed.type_bases("ACGTACGT")
        ed.set_cursor(0)              # no selection over the typed bases
        app.processEvents()
        blue = _count(view, colour)
        assert blue > 0, f"edited bases are not painted (dark={dark})"
        # no background: the yellow of a search hit is not on screen, and no wide block of blue either
        assert _yellow(view) == 0
        # a search hit over the edited bases is still plainly a hit (yellow behind them) and their letters stay blue
        ed.set_search_hits([(10, 18)])
        app.processEvents()
        assert _yellow(view) > 100, f"a search hit over edited bases is not visible (dark={dark})"
        assert _count(view, colour) > 0, "the edited letters must stay blue under a hit"
        ed.set_search_hits([])
        ed.undo()
        app.processEvents()
        assert _count(view, colour) == 0
        view.close()
    # blue, and readable on each theme's sequence background
    assert QColor(theme.ROLES["edited_fg"][0]).lightness() < 150 and QColor(theme.ROLES["edited_fg"][1]).lightness() > 128
    for c in theme.ROLES["edited_fg"]:
        assert QColor(c).blue() > QColor(c).red() + 60 and QColor(c).blue() > QColor(c).green() + 40
