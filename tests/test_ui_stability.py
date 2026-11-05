"""The interface does not move under the user's hand: no button changes size when it is pressed, toggled or hovered (and so
pushes its neighbours), table columns fit what they hold, counted words read as words. Headless (Qt offscreen)."""
import os
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("GENOMEEDIT_DATA_DIR", tempfile.mkdtemp(prefix="genomeedit_test_"))
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "app"))

import pytest
from PySide6.QtCore import Qt
from PySide6.QtWidgets import QAbstractButton, QApplication, QWidget

from uikit import plural

DEMO = ROOT / "demo" / "pFG051_Chlor.gb"


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


def _window(app, size=(1400, 850)):
    from main import MainWindow
    win = MainWindow()
    win.resize(*size)
    win.show()
    win.load(str(DEMO))
    app.processEvents()
    win.grab()
    return win


def _unstable_buttons(app, root):
    """Every visible button whose size hint changes when it is hovered, pressed or toggled, or whose toggling moves a sibling."""
    moved = []
    for b in [b for b in root.findChildren(QAbstractButton) if b.isVisible()]:
        parent = b.parentWidget()
        siblings = [w for w in parent.findChildren(QWidget, options=Qt.FindDirectChildrenOnly) if w.isVisible()] if parent else []

        def snap():
            app.processEvents()
            return b.sizeHint(), [w.geometry() for w in siblings]
        before = snap()
        b.blockSignals(True)
        b.setAttribute(Qt.WA_UnderMouse, True)
        hovered = snap()
        b.setAttribute(Qt.WA_UnderMouse, False)
        b.setDown(True)
        down = snap()
        b.setDown(False)
        toggled = before
        if b.isCheckable():
            b.setChecked(not b.isChecked())
            toggled = snap()
            b.setChecked(not b.isChecked())
            snap()
        b.blockSignals(False)
        if hovered[0] != before[0] or down[0] != before[0] or toggled != before:
            moved.append(b.text() or b.toolTip()[:40] or type(b).__name__)
    return moved


def test_no_button_changes_size_when_pressed_toggled_or_hovered(app):
    """The toolbar toggles (Lock, Filters, History, Circular) used to shrink or grow by 5 px when pressed: a style sheet styled only
    their checked state. The main window, its side panel and the quick filters are checked; none may move."""
    win = _window(app)
    assert _unstable_buttons(app, win) == []
    win.display_dock.show()
    assert _unstable_buttons(app, win.display_dock) == []
    win.history_dock.show()
    assert _unstable_buttons(app, win.history_dock) == []


def test_feature_table_columns_fit_their_content_in_a_tab_of_its_own(app):
    """A table alone in a wide tab: no column is spread far wider than what it holds (the Strand column used to take 800 px for
    an arrow), and the Name column is not cut off while there is room."""
    import workspace
    win = _window(app)
    ws = win.workspace
    ws.add_tab("Features", workspace.pane("Features"), select=True)
    app.processEvents()
    app.processEvents()
    table = next(v for v in ws.currentWidget().findChildren(QWidget) if type(v).__name__ == "GroupTableView" and v.isVisible())
    header, model = table.horizontalHeader(), table.model()
    fm = table.fontMetrics()
    for c in range(model.columnCount()):
        need = max(fm.horizontalAdvance(str(model.headerData(c, Qt.Horizontal))) + 24,
                   *(fm.horizontalAdvance(str(model.data(model.index(r, c)))) + 16 for r in range(model.rowCount())))
        width = header.sectionSize(c)
        assert width <= max(2.2 * need, need + 120), f"column {model.headerData(c, Qt.Horizontal)!r}: {width} px for {need} px"
    name_need = max(fm.horizontalAdvance(str(model.data(model.index(r, 0)))) for r in range(model.rowCount()))
    assert header.sectionSize(0) >= name_need, "the names are cut off although the tab has room"


def test_counted_words():
    assert [plural(1, "primer"), plural(3, "primer"), plural(1234, "feature"), plural(2, "copy"), plural(0, "mismatch"),
            plural(2, "entry"), plural(5, "bp")] == ["1 primer", "3 primers", "1,234 features", "2 copies", "0 mismatches",
                                                     "2 entries", "5 bp"]


def test_the_new_tab_button_stays_beside_the_tabs_and_in_sight(app):
    """The "+" follows the last tab; with more tabs than room (they scroll) it stays inside the bar, left of the scroll arrows,
    instead of being left after the last tab out of sight. The quick filters take only their buttons' width."""
    import workspace
    win = _window(app)
    ws = win.workspace
    bar = ws.tabBar()
    last = bar.tabRect(bar.count() - 1)
    assert 0 <= bar.add.geometry().left() - last.right() <= 8, "right after the last tab"
    for name in workspace.VIEWS:
        ws.add_tab(name, workspace.pane(name))
    app.processEvents()
    app.processEvents()
    assert bar.rect().contains(bar.add.geometry()), "still in sight when the tabs scroll"
    assert ws.quick_filters.width() <= ws.quick_filters.sizeHint().width() + 4, "the quick filters do not take a quarter of the row"


def test_the_toolbar_fits_without_overflow_down_to_1024_px(app):
    from PySide6.QtWidgets import QToolBar
    win = _window(app, (1024, 768))
    toolbar = win.findChildren(QToolBar)[0]
    hidden = [a.iconText() or a.text() for a in toolbar.actions()
              if not a.isSeparator() and a.isVisible() and toolbar.widgetForAction(a) is not None
              and not toolbar.widgetForAction(a).isVisible()]
    assert hidden == [], f"in the overflow menu: {hidden}"


def test_menus_in_the_usual_order_each_item_its_own_letter(app):
    """File, Edit, View first, Help last; in every menu each underlined letter (Alt+letter) belongs to one item only."""
    import re
    from PySide6.QtWidgets import QMenu
    win = _window(app)
    titles = [a.text().replace("&", "") for a in win.menuBar().actions()]
    assert titles[:3] == ["File", "Edit", "View"] and titles[3:5] == ["Import", "Export"] and titles[-1] == "Help"

    def check(menu, path):
        letters = [m.group(1).lower() for a in menu.actions() if not a.isSeparator()
                   for m in [re.search(r"&([^&])", a.text())] if m]
        assert len(letters) == len(set(letters)), f"{path}: {sorted(letters)}"
        for a in menu.actions():
            if a.menu() is not None:
                check(a.menu(), f"{path} > {a.text()}")
    for a in win.menuBar().actions():
        check(a.menu(), a.text())
    assert len({a.text()[a.text().index('&') + 1].lower() for a in win.menuBar().actions()}) == len(titles)
    assert isinstance(win.menuBar().actions()[0].menu(), QMenu)


def test_the_view_picker_keeps_its_width_whatever_view_it_shows(app):
    import workspace
    win = _window(app)
    pane = win.workspace.panes()[0]
    widths = set()
    for name in workspace.VIEWS:
        pane.picker.set_current(name)
        app.processEvents()
        widths.add(pane.picker.width())
    assert len(widths) == 1, widths


def test_the_window_remembers_its_size_view_switches_and_column_widths(app):
    """Size and place, the View switches and a column width dragged by hand come back in the next window; which side panel was
    open does not (nor the tabs)."""
    from main import MainWindow
    win = _window(app, (790, 560))   # (far from the first-start size on any test screen: 85 % of it)
    size = (win.width(), win.height())   # what the window really took (a platform may keep it wider: its minimum width)
    win.translation_action.setChecked(not win.translation_action.isChecked())
    wanted = win.translation_action.isChecked()
    win.display_dock.show()
    table = win.feature_table.view
    table.horizontalHeader().resizeSection(1, 151)   # as a drag by hand
    win.close()
    again = MainWindow()
    # the same size, give or take the few pixels a platform adds for its minimum width (Windows: 792 came back as 798)
    assert abs(again.width() - size[0]) <= 10 and abs(again.height() - size[1]) <= 10, ((again.width(), again.height()), size)
    assert again.translation_action.isChecked() == wanted
    again.show()
    again.load(str(DEMO))
    app.processEvents()
    app.processEvents()
    assert not again.display_dock.isVisible(), "which side panel was open is not remembered"
    assert again.feature_table.view.horizontalHeader().sectionSize(1) == 151


def test_the_range_box_reads_the_usual_ways_of_writing_a_stretch():
    from main import parse_range
    assert parse_range("101..200", 3780, True) == (100, 200)
    assert parse_range(" 101-200 ", 3780, True) == parse_range("101 to 200", 3780, True) == parse_range("200..101", 3780, False)
    assert parse_range("1,000..2,000", 3780, True) == (999, 2000)
    assert parse_range("3700..50", 3780, True) == (3699, 3830), "through the origin of a circular sequence"
    assert parse_range("3700..50", 3780, False) == (49, 3700), "a linear one: the two ends in order"
    assert parse_range("42", 3780, True) == (42, 42), "one number: the cursor after that base"
    for bad in ("", "abc", "1..2..3", "0..10", "100..4000"):
        with pytest.raises(ValueError):
            parse_range(bad, 3780, True)
