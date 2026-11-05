"""The views under the user's hand: keys and wheel steps of the sequence view, its drag selection scrolling on past the edge, the
menus of a right click on empty space, a map click that does not turn into a selection, image exports on a light background, and
an alignment view that does not grow taller than the screen. Headless (Qt offscreen)."""
import os
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("GENOMEEDIT_DATA_DIR", tempfile.mkdtemp(prefix="genomeedit_test_"))
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "app"))

from PySide6.QtCore import QEvent, QPoint, QPointF, Qt
from PySide6.QtGui import QContextMenuEvent, QKeyEvent, QMouseEvent, QWheelEvent
from PySide6.QtWidgets import QApplication, QMenu

from genomeedit_core import Document

import sequence_view
from editor import Editor
from map_view import MapView
from sequence_view import SequenceView

app = QApplication.instance() or QApplication([])
SEQ = "ACGTTGCAAGGCTTAACCGT" * 300   # 6,000 bp: many rows


def view_of(seq=SEQ, circular=False, size=(900, 400)):
    editor = Editor(Document.new_sequence("t", seq, circular), Path("t.gb"))
    view = SequenceView(editor)
    view.resize(*size)
    view.show()
    app.processEvents()
    return view, editor


def key(view, k, mods=Qt.NoModifier, text=""):
    view.keyPressEvent(QKeyEvent(QEvent.KeyPress, k, mods, text))


def wheel(widget, dy, mods=Qt.NoModifier, dx=0, at=QPointF(100, 100)):
    event = QWheelEvent(at, widget.mapToGlobal(at), QPoint(), QPoint(dx, dy), Qt.NoButton, mods, Qt.NoScrollPhase, False)
    widget.wheelEvent(event)


def mouse(widget, kind, at, button=Qt.LeftButton, buttons=Qt.LeftButton, mods=Qt.NoModifier):
    at = QPointF(at)
    event = QMouseEvent(kind, at, widget.mapToGlobal(at), button, buttons, mods)
    {QEvent.MouseButtonPress: widget.mousePressEvent, QEvent.MouseMove: widget.mouseMoveEvent,
     QEvent.MouseButtonRelease: widget.mouseReleaseEvent}[kind](event)


class Recorder(QMenu):
    """A menu that records itself instead of opening (exec would wait for a click)."""
    shown = []

    def exec(self, *args):
        Recorder.shown.append(self)


def menu_texts(menu):
    return {a.text().split("\t")[0]: a.isEnabled() for a in menu.actions() if not a.isSeparator()}


# ---- keys ------------------------------------------------------------------------------------------------------------


def test_ctrl_home_end_and_page_keys():
    view, ed = view_of()
    ed.set_cursor(1234)
    key(view, Qt.Key_End, Qt.ControlModifier)
    assert ed.start == ed.end == ed.length, "Ctrl+End: the end of the sequence (not of the row)"
    assert view.verticalScrollBar().value() == view.verticalScrollBar().maximum()
    key(view, Qt.Key_Home, Qt.ControlModifier | Qt.ShiftModifier)
    assert (ed.start, ed.end) == (0, ed.length), "Ctrl+Shift+Home selects to the start"
    key(view, Qt.Key_Home, Qt.ControlModifier)
    assert ed.start == ed.end == 0 and view.verticalScrollBar().value() == 0
    key(view, Qt.Key_End)
    assert ed.start == view.bases_per_line(), "End alone: the end of the row"
    ed.set_cursor(0)
    view.anchor = 0
    key(view, Qt.Key_PageDown)
    moved = ed.start
    assert moved > 2 * view.bases_per_line() and moved % view.bases_per_line() == 0, "Page Down moves the cursor by a screenful"
    assert view.verticalScrollBar().value() > 0, "and the page scrolls with it"
    key(view, Qt.Key_PageDown, Qt.ShiftModifier)
    assert ed.start == moved and ed.end > moved, "Shift+Page Down extends the selection"
    key(view, Qt.Key_PageUp)
    key(view, Qt.Key_PageUp)
    assert ed.start == ed.end and ed.start < moved


def test_the_cursor_is_faint_in_a_view_without_the_keyboard():
    """With several panes, the cursor is solid only in the view where typing goes."""
    from PySide6.QtGui import QColor
    view, ed = view_of(size=(700, 300))
    ed.set_cursor(5)
    x = view.left + 5 * view.cw
    top, bottom = int(view.offset(0)), int(view.offset(1))
    background = view.palette().color(view.viewport().backgroundRole())

    def ink():
        """How far from the background a strip a few pixels wide around the cursor is, over the whole row. The letters are the
        same with and without the keyboard, so only the cursor makes a difference (one pixel was fragile: with the fractional
        letter widths of another platform's font it fell on the edge of a letter)."""
        image = view.viewport().grab().toImage()
        total = 0
        for px in range(int(x) - 2, int(x) + 3):
            for py in range(top, bottom):
                c = QColor(image.pixel(px, py))
                total += abs(c.red() - background.red()) + abs(c.green() - background.green()) + abs(c.blue() - background.blue())
        return total
    view.activateWindow()
    view.setFocus()
    app.processEvents()
    assert view.hasFocus()
    solid = ink()
    view.clearFocus()
    app.processEvents()
    assert not view.hasFocus()
    faint = ink()
    assert faint < solid, (faint, solid)


# ---- wheel -----------------------------------------------------------------------------------------------------------


def test_wheel_zoom_steps_once_per_notch():
    view, ed = view_of()
    size = view.size
    wheel(view, 60, Qt.ControlModifier)
    assert view.size == size, "half a notch (a touchpad's small step) does not zoom yet"
    wheel(view, 60, Qt.ControlModifier)
    assert view.size == size + 1, "two halves make one step"
    wheel(view, 0, Qt.ControlModifier, dx=120)
    assert view.size == size + 1, "a sideways Ctrl+scroll (no vertical delta) does not zoom out"
    wheel(view, -240, Qt.ControlModifier)
    assert view.size == size - 1, "two notches back: two steps"
    # the map zooms by the same rule (a linear map: it zooms)
    m = MapView(ed)
    m.resize(800, 400)
    m.show()
    app.processEvents()
    wheel(m, 60, at=QPointF(400, 200))
    assert m.view is None
    wheel(m, 60, at=QPointF(400, 200))
    assert m.view is not None, "two half notches: one zoom step"
    span = m.view[1] - m.view[0]
    wheel(m, 0, dx=0, at=QPointF(400, 200))
    assert m.view[1] - m.view[0] == span


# ---- mouse -----------------------------------------------------------------------------------------------------------


def test_a_map_click_that_wobbles_is_not_a_selection():
    view, ed = view_of()
    m = MapView(ed)
    m.resize(800, 400)
    m.show()
    app.processEvents()
    x0, w, y, _ = m._linear()
    at = QPointF(x0 + w / 2, y)
    mouse(m, QEvent.MouseButtonPress, at)
    cursor = ed.start
    mouse(m, QEvent.MouseMove, at + QPointF(1, 0))
    assert ed.start == ed.end == cursor, "a 1 px move after the press does not select"
    mouse(m, QEvent.MouseMove, at + QPointF(QApplication.startDragDistance() + 30, 0))
    assert ed.end > ed.start, "a real drag selects"
    mouse(m, QEvent.MouseButtonRelease, at, buttons=Qt.NoButton)


def test_dragging_past_the_bottom_edge_scrolls_and_extends_the_selection():
    view, ed = view_of()
    vp = view.viewport()
    y = view.offset(0) + view._parts(0)[1] + 2
    mouse(view, QEvent.MouseButtonPress, QPointF(view.left + 2 * view.cw, y))
    mouse(view, QEvent.MouseMove, QPointF(view.left + 5 * view.cw, vp.height() + 40))
    assert view.scroll_timer.isActive(), "out of the bottom with the button held: the timer scrolls"
    end, scroll = ed.end, view.verticalScrollBar().value()
    for _ in range(5):
        view._autoscroll()
    assert view.verticalScrollBar().value() > scroll and ed.end > end, "it scrolled and took the selection along"
    mouse(view, QEvent.MouseButtonRelease, QPointF(view.left, vp.height() + 40), buttons=Qt.NoButton)
    assert not view.scroll_timer.isActive()


# ---- menus -----------------------------------------------------------------------------------------------------------


def test_right_click_on_empty_space_offers_cut_copy_paste_select_all(monkeypatch):
    monkeypatch.setattr(sequence_view, "QMenu", Recorder)
    Recorder.shown.clear()
    view, ed = view_of()
    QApplication.clipboard().setText("ACGT")
    ed.set_cursor(10)
    mouse(view, QEvent.MouseButtonPress, QPointF(view.left + 3 * view.cw, view._parts(0)[1] + 2), Qt.RightButton, Qt.RightButton)
    items = menu_texts(Recorder.shown[-1])
    assert items == {"Cut": False, "Copy": False, "Paste": True, "Select All": True}, items
    # Select All from the menu
    next(a for a in Recorder.shown[-1].actions() if a.text().startswith("Select All")).trigger()
    assert (ed.start, ed.end) == (0, ed.length)
    # the DNA lock greys Paste (and says why); the Menu key opens the same menu
    ed.set_cursor(10)
    ed.dna_locked = True
    view.contextMenuEvent(QContextMenuEvent(QContextMenuEvent.Keyboard, QPoint(5, 5), QPoint(5, 5)))
    items = menu_texts(Recorder.shown[-1])
    assert items["Paste (DNA locked)"] is False and items["Select All"] is True
    # a view that refuses right clicks (the Assembly panel's, which only shows a sequence) refuses the Menu key too
    import assembly_panel
    shown_only = assembly_panel.ReadOnlyView(ed)
    count = len(Recorder.shown)
    shown_only.contextMenuEvent(QContextMenuEvent(QContextMenuEvent.Keyboard, QPoint(5, 5), QPoint(5, 5)))
    assert len(Recorder.shown) == count
    # the map has the same menu on empty space
    m = MapView(ed)
    m.resize(800, 400)
    m.show()
    app.processEvents()
    count = len(Recorder.shown)
    mouse(m, QEvent.MouseButtonPress, QPointF(5, 5), Qt.RightButton, Qt.RightButton)
    assert len(Recorder.shown) == count + 1 and "Select All" in menu_texts(Recorder.shown[-1])


# ---- the alignment view ----------------------------------------------------------------------------------------------


def test_the_alignment_view_scrolls_instead_of_growing_past_the_screen():
    """Tall tracks (many BAM files) scroll inside the view: neither the view nor the window asks for more than the screen, the
    pane is given at most part of it; a long info line does not widen the pane; a click does not take the keyboard."""
    from main import MainWindow
    win = MainWindow()
    win.resize(1400, 850)
    win.show()
    win.load(str(ROOT / "demo" / "pFG051_Chlor.gb"))
    win.set_genome_visible(True)
    for _ in range(3):
        app.processEvents()
    gv = win.genome_view
    for name in ("GC %", "Features"):
        gv.canvas.style(name).height = 3000
    gv.refresh()
    for _ in range(3):
        app.processEvents()
    screen = gv.screen().availableGeometry().height()
    assert gv.canvas.minimumHeight() > 6000
    assert gv.minimumSizeHint().height() < screen and win.minimumSizeHint().height() < screen
    assert gv.tracks_area.verticalScrollBar().maximum() > 0, "the tracks scroll"
    assert gv.height() <= screen * 0.6 + 80, "the pane takes part of the screen, not all of it"
    width = win.width()
    gv.info.setText("read_1234  " * 300)
    app.processEvents()
    assert win.width() == width and gv.info.text().startswith("read_1234") and gv.info.toolTip() == gv.info.text()
    assert gv.canvas.focusPolicy() == Qt.NoFocus and gv.tracks_area.focusPolicy() == Qt.NoFocus
    for name in ("GC %", "Features"):
        gv.canvas.style(name).height = None
    win.editor.doc.modified = False
    win.close()


# ---- panels that are filled anew after an edit -------------------------------------------------------------------------


def test_checkup_tables_keep_their_row_and_scroll_and_fit_their_columns():
    import checkup
    path = ROOT / "demo" / "pFG051_Chlor.gb"
    ed = Editor(Document(str(path)), path)
    panel = checkup.CheckupDialog(None, ed)
    panel.resize(760, 300)
    panel.show()
    app.processEvents()
    table = panel.element_table
    assert table.rowCount() > 6
    table.setCurrentCell(5, 0)
    table.verticalScrollBar().setValue(table.verticalScrollBar().maximum())
    name, scroll = table.item(5, 0).text(), table.verticalScrollBar().value()
    assert scroll > 0, "(the table is scrolled for this test)"
    panel.refresh()
    app.processEvents()
    assert table.currentRow() == 5 and table.item(5, 0).text() == name, "the current row is kept"
    assert table.verticalScrollBar().value() == scroll, "and the scroll position"
    # columns fit what they hold: the details column is not spread over 1000 px, the others are not cut off
    app.processEvents()
    header = table.horizontalHeader()
    for c in range(table.columnCount() - 1):
        need = max(table.sizeHintForColumn(c), header.sectionSizeHint(c))
        assert header.sectionSize(c) >= min(need, 260) - 2, table.horizontalHeaderItem(c).text()
    panel.close()


def test_protein_view_keeps_its_place_through_an_edit():
    """An edit elsewhere shows the same protein again: the text stays where it was scrolled and the letters selected in it stay
    selected. The title is one line (its whole text in the tooltip when cut short); the features list is hidden when empty."""
    import protein
    import random
    rng = random.Random(5)
    coding = "ATG" + "".join(rng.choice(["GCT", "AAA", "TGG", "CAT", "GAA", "CTG"]) for _ in range(600)) + "TAA"
    seq = "TTTTTTTTTT" + coding + "GGGGGGGGGG" * 3
    ed = Editor(Document.new_sequence("long", seq, False), Path("long.gb"))
    ed.add_features([("CDS", f"11..{10 + len(coding)}", [("label", "a CDS with a name long enough to be cut short in a pane")])])
    view = protein.ProteinView(ed)
    view.resize(300, 250)
    view.show()
    cds = next(f for f in ed.features if f.kind == "CDS")
    ed.select_feature(cds.index)
    view.refresh(True)
    app.processEvents()
    assert not view.features.isVisible(), "no feature inside the CDS: no empty list"
    assert view.title.toolTip() and "a CDS with a name" in view.title.text(), "cut short on one line, whole in the tooltip"
    bar = view.text.verticalScrollBar()
    bar.setValue(bar.maximum() // 2)
    cursor = view.text.textCursor()
    cursor.setPosition(200)
    cursor.setPosition(205, cursor.MoveMode.KeepAnchor)
    view.own = True   # (selecting letters selects their codons in the sequence: not what this test is about)
    view.text.setTextCursor(cursor)
    view.own = False
    scroll, selected = bar.value(), view.text.textCursor().selectedText()
    ed.set_feature(cds.index, "CDS", cds.location, list(cds.qualifiers) + [("note", "edited")])   # an edit: a note added
    app.processEvents()
    view.refresh(True)
    assert bar.value() == scroll and view.text.textCursor().selectedText() == selected


# ---- image export ------------------------------------------------------------------------------------------------------


def test_export_is_light_on_a_dark_program_and_leaves_the_controls_out(tmp_path):
    import image_export
    from PySide6.QtGui import QColor, QImage, QPalette
    old = QPalette(app.palette())
    dark = QPalette(old)
    for role, color in ((QPalette.Window, "#2b2d30"), (QPalette.Base, "#1e1f22"), (QPalette.Text, "#e8eaed"),
                        (QPalette.WindowText, "#e8eaed"), (QPalette.PlaceholderText, "#8b9096")):
        dark.setColor(role, QColor(color))
    app.setPalette(dark)
    try:
        view, ed = view_of(size=(700, 300))
        app.processEvents()
        assert view.palette().color(QPalette.Base).lightness() < 60, "(the program is dark for this test)"
        bar = view.verticalScrollBar().parentWidget()
        assert view.zoom_bar.isVisible() and bar.isVisible()
        for light, path in ((True, tmp_path / "light.png"), (False, tmp_path / "dark.png")):
            image_export.write(view, path, "png", 1.0, 96, light=light)
        light_img, dark_img = QImage(str(tmp_path / "light.png")), QImage(str(tmp_path / "dark.png"))
        assert QColor(light_img.pixel(2, 2)).lightness() > 200 and QColor(dark_img.pixel(2, 2)).lightness() < 60
        # the zoom bar and the scroll bar (right edge) are not in the picture: that strip is the background
        zoom_x = view.viewport().mapTo(view, QPoint(view.viewport().width() - 12, 0)).x()
        for x in (zoom_x, light_img.width() - 5):   # the zoom bar's middle, the scroll bar's
            assert {QColor(light_img.pixel(x, y)).name() for y in range(2, 130)} == {"#ffffff"}, \
                f"no zoom bar or scroll bar in the export (x {x})"
        assert view.zoom_bar.isVisible() and bar.isVisible() and view.palette().color(QPalette.Base).lightness() < 60, \
            "all put back afterwards"
    finally:
        app.setPalette(old)
        app.processEvents()


def test_the_export_dialog_does_not_grow_when_it_warns():
    import image_export
    view, ed = view_of(size=(1400, 800))
    d = image_export.ImageExportDialog(None, view)
    d.show()
    app.processEvents()
    before = d.sizeHint().height()
    d.dpi.setValue(2400)
    d.scale.setValue(1000)
    app.processEvents()
    assert "Too large" in d.warning.text() and not d.ok.isEnabled()
    assert d.sizeHint().height() == before, "the warning has its room already"
    assert d.light.isChecked(), "light background: on by default"
    d.close()


# ---- the 3D structure window ---------------------------------------------------------------------------------------------


def test_a_long_local_structure_command_shows_progress_and_can_be_cancelled():
    import time
    from structure_dialog import StructureDialog
    d = StructureDialog(None, "tiny", "MKWVTFISLLLLFSSAYS")
    d.source.setCurrentText("Local command")
    d.command.setText(f'"{sys.executable}" -c "import time; time.sleep(60)" {{fasta}} {{pdb}}')
    d.show()
    d.start()
    t = time.perf_counter()
    while (d.job is None or d.job.proc is None) and time.perf_counter() - t < 10:
        app.processEvents()
        time.sleep(0.02)
    proc = d.job.proc
    assert proc is not None and proc.poll() is None and d.busy.isVisible() and d.cancel_button.isVisible()
    d.cancel()
    proc.wait(timeout=10)
    assert proc.poll() is not None, "the command was stopped"
    assert d.get.isEnabled() and not d.busy.isVisible() and "Cancelled" in d.status.text()
    t = time.perf_counter()
    while time.perf_counter() - t < 0.5:   # the worker's late answer is not shown
        app.processEvents()
    assert "Cancelled" in d.status.text()
    d.settings.clear()
    d.close()
