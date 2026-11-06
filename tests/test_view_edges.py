"""The sequence view at the edges of a sequence (a primer's tail hanging off the start or the end is drawn in full) and cuts drawn on both strands."""
import os
import random
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ["OPENDNA_DATA_DIR"] = tempfile.mkdtemp(prefix="opendna_test_")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "app"))

from PySide6.QtCore import Qt
from PySide6.QtGui import QColor
from PySide6.QtWidgets import QApplication

from opendna_core import Document

import assembly_panel as P
import sequence_view
from editor import reverse_complement as rc

app = QApplication.instance() or QApplication([])
rng = random.Random(3)
SEQ = "".join(rng.choice("ACGT") for _ in range(120))


def view_with(primers, circular=False):
    doc = Document.new_sequence("t", SEQ, circular)
    editor = P.FragmentEditor(doc, Path("t.gb"))
    editor.set_extra(primers)
    view = P.ReadOnlyView(editor)
    view.resize(1000, 400)
    view.show()
    app.processEvents()
    return view


def primer_pixels(view, x0, x1):
    """How many pixels between x0 and x1 (viewport x) have the colour of a primer tail's bases."""
    img = view.viewport().grab().toImage()
    want = QColor(sequence_view.ALMOST_COLOR)   # the tail is drawn as its own bases, in this orange
    return sum(1 for x in range(int(x0), int(x1)) for y in range(img.height())
               if abs(QColor(img.pixel(x, y)).red() - want.red()) < 6 and abs(QColor(img.pixel(x, y)).green() - want.green()) < 6
               and abs(QColor(img.pixel(x, y)).blue() - want.blue()) < 6)


def test_primer_tails_hanging_off_the_edges_are_drawn():
    """A tail hanging off the start wraps to the row end, one across a plasmid's origin shows on the last row, a reverse primer's tail
    hanging off the end is drawn too, and the rows are as wide at every zoom."""
    # --- a tail hanging off the start wraps to the row end and is drawn ---
    tail = "GGTCTCAAATGC"                                    # 12 nt that match nothing before base 1
    view = view_with([("F", tail + SEQ[:22])])
    (site,) = [s for s in view.editor.primer_sites if s.label == "F"]
    assert (site.start, site.strand, site.tail) == (0, 1, 12)
    assert view.left == view.base_left, "no room is made: the bases never move"
    # The tail goes on at the end of the row, as it would on a row above.
    right = view.left + view.bases_per_line() * view.cw
    assert primer_pixels(view, right - 12 * view.cw, right) > 0
    # --- a tail across the origin of a plasmid shows on the last row ---
    tail = "GGTCTCAAATGC"
    view = view_with([("F", tail + SEQ[:22])], circular=True)
    (site,) = [s for s in view.editor.primer_sites if s.label == "F"]
    assert site.extent[-1] == (len(SEQ) - 12, len(SEQ))
    view.verticalScrollBar().setValue(view.verticalScrollBar().maximum())
    last_row_bases = len(SEQ) % view.bases_per_line() or view.bases_per_line()
    x_end = view.left + last_row_bases * view.cw
    assert primer_pixels(view, x_end - 12 * view.cw, x_end) > 0
    # --- the rows are as wide at every zoom ---
    view = view_with([])
    widths = set()
    for size in (6, 9, 11, 16, 24):
        view.set_size(size)
        widths.add(round(view.left + view.bases_per_line() * view.cw))
    assert len(widths) == 1
    # --- a reverse primer tail hanging off the end is drawn too ---
    # The reverse primer's 5' end is on the right: its tail continues past the last base.
    tail = "CCGGTTAACC"
    view = view_with([("R", tail + rc(SEQ[-22:]))])
    (site,) = [s for s in view.editor.primer_sites if s.label == "R"]
    assert site.strand == -1 and site.tail == 10 and site.end == len(SEQ)
    bpl = view.bases_per_line()
    last_row_bases = len(SEQ) % bpl or bpl
    x_end = view.left + last_row_bases * view.cw
    # It continues right of the last base if the row has room, else at the start of the row (as on the next row).
    room = view.left + bpl * view.cw - x_end
    start, stop = (x_end + 1, x_end + 10 * view.cw) if room >= 10 * view.cw else (view.left, view.left + 10 * view.cw)
    assert primer_pixels(view, start, stop) > 0, "the tail continues after the last base"
    view.viewport().repaint()


# ---- cuts drawn on both strands ---------------------------------------------------------------------------------------


def red_at(view, x, y):
    c = QColor(view.viewport().grab().toImage().pixel(int(x), int(y)))
    return c.red() > 200 and c.green() < 90 and c.blue() < 90


def test_cuts_on_both_strands_are_known_drawn_and_shown_by_the_digest():
    """Enzymes know where each strand is cut, cuts are drawn on the strand they cut (joined when sticky), and the digest shows its cuts."""
    import enzymes
    from tools import DigestDialog
    # --- enzymes know where each strand is cut ---
    # EcoRI G^AATTC: top strand after base 1 of the site, bottom strand after base 5 (a 4 nt 5' overhang).
    cuts = enzymes.find_cuts("TTGAATTCAA", False)
    assert cuts["EcoRI"] == [(3, 7)]
    assert enzymes.find_sites("TTGAATTCAA", False)["EcoRI"] == [3], "find_sites still gives the top-strand cuts"
    # A blunt cutter cuts both strands in the same place.
    assert enzymes.find_cuts("AAGATATCAA", False)["EcoRV"] == [(5, 5)]
    # BsaI GGTCTC(1/5) on either strand: the 4 nt overhang is the same bases between the two cuts.
    fwd = "ATGGTCTCAAAAAAAAAAAAAA"
    (t, b), = enzymes.find_cuts(fwd, False)["BsaI"]
    assert (t, b) == (9, 13) and b - t == 4
    rev = "TTTTTTTTTTTTGAGACCAT"
    (t2, b2), = enzymes.find_cuts(rev, False)["BsaI"]
    assert b2 - t2 == 4 and (t2, b2) == (7, 11), "on the bottom strand the cut is to the left of the site"
    # Across the origin of a circle: positions wrap.
    (tc, bc), = enzymes.find_cuts("ATTCAAAAAAGA", True)["EcoRI"]
    assert (tc, bc) == (11, 3)
    # --- cuts are drawn on the strand they cut and joined when sticky ---
    view = view_with([])
    ed = view.editor
    ed.set_cut_marks("t", [(11, 15, "EcoRI"), (40, 40, "blunt")])
    app.processEvents()
    y = view.offset(0) - view.verticalScrollBar().value() + view.layout_seq_top(0) if hasattr(view, "layout_seq_top") else \
        view.offset(0) - view.verticalScrollBar().value() + view._parts(0)[1]
    lh, left, cw = view.lh, view.left, view.cw
    top_half, bottom_half = y + lh * 0.5, y + lh * 1.5
    # Sticky: the top strand is cut at 11 and the bottom strand at 15 (never the other way round)...
    assert red_at(view, left + 11 * cw, top_half) and not red_at(view, left + 11 * cw, bottom_half)
    assert red_at(view, left + 15 * cw, bottom_half) and not red_at(view, left + 15 * cw, top_half)
    # ...with a dotted line across the middle joining them.
    assert any(red_at(view, left + x * cw, y + lh) for x in (12.5, 13, 13.5, 14))
    # Blunt: both strands in the same place.
    assert red_at(view, left + 40 * cw, top_half) and red_at(view, left + 40 * cw, bottom_half)
    # Removing the marks removes the drawing.
    ed.set_cut_marks("t", [])
    app.processEvents()
    assert not red_at(view, left + 40 * cw, top_half) and not ed.cut_marks
    # --- the digest shows its cuts on the sequence ---
    seq = "ACGT" * 10 + "GAATTC" + "ACGT" * 10 + "GGATCC" + "TTAA" * 5
    doc = Document.new_sequence("t", seq, False)
    ed = P.FragmentEditor(doc, Path("t.gb"))
    dlg = DigestDialog(None, ed)
    assert "EcoRI" in dlg.enzymes.names() and not ed.cut_marks, "nothing is drawn until an enzyme is ticked"
    item = dlg.enzymes.item_of("EcoRI")
    item.setCheckState(0, Qt.Checked)
    assert [(m.top, m.bottom, m.label) for m in ed.cut_marks["digest"]] == [(41, 45, "EcoRI")]
    dlg.enzymes.item_of("BamHI").setCheckState(0, Qt.Checked)
    assert sorted(m.label for m in ed.cut_marks["digest"]) == ["BamHI", "EcoRI"]
    item.setCheckState(0, Qt.Unchecked)
    assert [m.label for m in ed.cut_marks["digest"]] == ["BamHI"]
    dlg.enzymes.item_of("BamHI").setCheckState(0, Qt.Unchecked)
    assert not ed.cut_marks


def test_clicking_the_alignment_view_and_hovering_features_in_the_sequence_view():
    """Clicking a base in the alignment view moves the sequence view to it; the pointer is a hand over a clickable feature."""
    # --- clicking a base in the alignment view moves the sequence view to it ---
    import main as M
    from PySide6.QtCore import QPointF
    from PySide6.QtGui import QMouseEvent
    big = "".join(rng.choice("ACGT") for _ in range(30000))
    win = M.MainWindow()
    win.show()
    win.resize(1500, 900)
    win.show_document(Document.new_sequence("big", big, False), Path("big.gb"), untitled=True)
    win.set_genome_visible(True)
    app.processEvents()
    seq, gv = win.sequence_view, win.genome_view
    canvas, sb = gv.canvas, seq.verticalScrollBar()
    assert sb.value() == 0

    def click(position):
        """A click on the alignment view's canvas over this base (its ruler is framed on the base first)."""
        gv.set_view(position - 1000, position + 1000)
        app.processEvents()
        x = canvas.left + (position - gv.vstart) / (gv.vend - gv.vstart) * canvas.width_px()
        for kind, send in ((QMouseEvent.Type.MouseButtonPress, canvas.mousePressEvent),
                           (QMouseEvent.Type.MouseButtonRelease, canvas.mouseReleaseEvent)):
            send(QMouseEvent(kind, QPointF(x, 90), QPointF(x, 90), Qt.LeftButton, Qt.LeftButton, Qt.NoModifier))
        app.processEvents()

    def shown(position):
        row = position // seq.bases_per_line()
        first, last = seq.row_at(sb.value()), seq.row_at(sb.value() + seq.viewport().height())
        return first <= row <= last

    click(20000)
    assert abs(win.editor.start - 20000) <= 2, "the base was selected"
    assert sb.value() > 0 and shown(win.editor.start), "the sequence view went to it"
    middle = (seq.offset(win.editor.start // seq.bases_per_line()) - sb.value()) / seq.viewport().height()
    assert 0.2 < middle < 0.8, "brought to the middle, not to an edge"
    # A click next to what is on screen does not make the page jump.
    before = sb.value()
    click(win.editor.start + 30)
    assert sb.value() == before
    # And back to the start.
    click(150)
    assert shown(win.editor.start) and sb.value() < before
    win.editor.doc.modified = False
    win.reset_layout()
    win.close()
    # --- the pointer is a hand over a clickable feature in the sequence view ---
    doc = Document.new_sequence("t", SEQ, False)
    doc.add_features([("CDS", "10..40", [("label", "gene")])])
    view = P.ReadOnlyView(P.FragmentEditor(doc, Path("t.gb")))
    view.resize(1000, 400)
    view.show()
    app.processEvents()
    row_y = view.offset(0) - view.verticalScrollBar().value()
    _, seq_top, _, lanes_y = view._parts(0)

    def hover(x, y):
        view.mouseMoveEvent(QMouseEvent(QMouseEvent.Type.MouseMove, QPointF(x, y), QPointF(x, y), Qt.NoButton, Qt.NoButton, Qt.NoModifier))
        return view.viewport().cursor().shape()

    assert hover(view.left + 25 * view.cw, row_y + lanes_y + view.lane_h / 2) == Qt.PointingHandCursor, "over the feature"
    assert hover(view.left + 25 * view.cw, row_y + seq_top + view.lh / 2) == Qt.IBeamCursor, "over the bases: a text cursor"
    assert hover(view.left + 80 * view.cw, row_y + lanes_y + view.lane_h / 2) == Qt.IBeamCursor, "next to the feature"


def test_the_top_strand_alone_is_one_line_of_bases_and_everything_else_stays_where_it_was():
    """Display > Complementary strand off: the rows lose the bottom strand (one line of bases, one line shorter), the primers, the
    features and the hover work as before, the selection, the origin mark and the cuts are as tall as the bases, and the switch is
    in the Display panel and the View menu."""
    from display_panel import DisplayPanel
    doc = Document.new_sequence("t", SEQ, True)
    doc.add_features([("CDS", "10..40", [("label", "gene")])])
    editor = P.FragmentEditor(doc, Path("t.gb"))
    editor.set_extra([("F", SEQ[50:72]), ("R", rc(SEQ[90:112]))])
    view = P.ReadOnlyView(editor)
    view.resize(1000, 400)
    view.show()
    app.processEvents()
    assert editor.show_complement and view.strand_h == 2 * view.lh, "both strands unless asked"
    _, seq_top, rev_top, _ = view._parts(0)
    both_h, both_rev = view.row_h, rev_top - seq_top
    assert both_rev == 2 * view.lh

    def below_the_top_strand():
        """Pixels under the first line of bases (down to where the second line would end) that are not the background."""
        img = view.viewport().grab().toImage()
        base = img.pixel(img.width() - 40, 2)   # (nothing is drawn in the corner)
        _, top, _, _ = view._parts(0)
        y0 = view.offset(0) - view.verticalScrollBar().value() + top
        return sum(1 for x in range(int(view.left + 42 * view.cw), int(view.left + (view.bases_per_line() - 1) * view.cw))   # past the CDS
                   for y in range(int(y0 + view.lh + 1), int(y0 + 2 * view.lh - 1)) if img.pixel(x, y) != base)

    assert below_the_top_strand() > 0, "the complement is drawn"
    editor.set_display(complement=False)
    app.processEvents()
    assert not editor.show_complement and view.strand_h == view.lh
    assert below_the_top_strand() == 0, "no second line of bases"
    _, seq_top, rev_top, lanes_y = view._parts(0)
    assert rev_top - seq_top == view.lh, "the reverse primers sit right under the one line of bases"
    assert abs((both_h - view.row_h) - view.lh) <= 1, "a row is one line shorter"
    # The primers and the feature are still there and hover as before (forward primer above, reverse below, feature in its lane).
    row_y = view.offset(0) - view.verticalScrollBar().value()
    sites = {s.label: s for s in editor.primer_sites}
    assert set(sites) == {"F", "R"} and (sites["F"].strand, sites["R"].strand) == (1, -1)
    bpl = view.bases_per_line()
    for label in ("F", "R"):
        site = sites[label]
        row = site.start // bpl
        y = view.offset(row) - view.verticalScrollBar().value() + view._primer_lane_top(row, site.strand, editor.primer_lanes[site.key]) + view.primer_h / 2
        assert view.hover_target(view.left + (site.start + 3 - row * bpl) * view.cw, y) == site.key
    assert isinstance(view.hover_target(view.left + 25 * view.cw, row_y + lanes_y + view.lane_h / 2), int), "the CDS"
    # A selection is as tall as the bases; the origin mark, the brackets and the cuts follow (nothing here may raise).
    editor.select(20, 30)
    editor.set_cut_marks("t", [(25, 27, "cut")])
    view.viewport().repaint()
    assert view.viewport().grab().toImage().width() > 0
    # The switch: in the Display panel (both ways), in the View menu of the main window, and in what a project remembers.
    panel = DisplayPanel(editor)
    assert not panel.strands.isChecked()
    panel.strands.setChecked(True)
    assert editor.show_complement and view.strand_h == 2 * view.lh and abs(view.row_h - both_h) <= 1
    editor.set_display(complement=False)
    assert not panel.strands.isChecked(), "the panel follows the editor"
    assert below_the_top_strand() == 0, "no second line of bases"
    from exports import project_meta
    assert project_meta(editor, Path("t.odna"))["display"]["complement"] is False
    import main as M
    win = M.MainWindow()
    win.show_document(Document.new_sequence("w", SEQ, False), Path("w.gb"), untitled=True)
    app.processEvents()
    assert win.complement_action.isChecked() and win.editor.show_complement, "two strands unless asked"
    win.complement_action.trigger()   # View > Complementary Strand
    assert not win.editor.show_complement and win.sequence_view.strand_h == win.sequence_view.lh
    win.editor.set_display(complement=True)   # ... and the other way round (the Display panel)
    assert win.complement_action.isChecked() and win.sequence_view.strand_h == 2 * win.sequence_view.lh
    win.editor.doc.modified = False
    win.close()
