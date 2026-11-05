"""The tool windows behave: they close with Escape, Enter does not press a button nobody chose, one Enter is one NCBI request,
the Digest's "+" sits right after its lanes, tables fit their content, the Primer Database's Delete key and buttons, a mapping
run is not hidden while it runs. Headless (Qt offscreen)."""
import os
import re
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("GENOMEEDIT_DATA_DIR", tempfile.mkdtemp(prefix="genomeedit_test_"))
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "app"))

import pytest
from PySide6.QtCore import Qt, QTimer
from PySide6.QtTest import QTest
from PySide6.QtWidgets import QApplication, QDialog, QMessageBox

DEMO = ROOT / "demo" / "pFG051_Chlor.gb"


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


@pytest.fixture
def win(app):
    from main import MainWindow
    w = MainWindow()
    w.resize(1400, 850)
    w.show()
    w.load(str(DEMO))
    app.processEvents()
    yield w
    w.editor.doc.modified = False
    w.reset_layout()   # closing saves the tabs: leave the default ones for the next window (tests share a data folder)
    w.close()


def _fits(view, model, column, need):
    """No column spread far wider than what it holds."""
    width = view.header().sectionSize(column) if hasattr(view, "header") else view.horizontalHeader().sectionSize(column)
    assert width <= max(2.2 * need, need + 120), f"column {model.headerData(column, Qt.Horizontal)!r}: {width} px for {need} px"


def test_spacer_finder_closes_with_escape(app, win):
    """`self.done = ...` hid QDialog.done(): Escape and the window's close button raised TypeError and the dialog never closed."""
    from cas9 import Cas9Dialog
    dlg = Cas9Dialog(win, win.editor)
    assert callable(dlg.done), "QDialog.done() is hidden by an attribute"   # (before exec(): it would never return)
    seen = {}

    def press():
        QTest.keyClick(dlg, Qt.Key_Escape)
        seen["closed"] = not dlg.isVisible()
        if dlg.isVisible():
            QDialog.done(dlg, QDialog.Accepted)   # (Qt's own done: never leave the test hanging in exec())
    QTimer.singleShot(300, press)
    result = dlg.exec()
    assert seen.get("closed") and result == QDialog.Rejected


def test_spacer_finder_columns_fit_their_content(app, win):
    from cas9 import Cas9Dialog
    dlg = Cas9Dialog(win, win.editor)
    dlg.resize(1800, 600)
    dlg.show()
    for _ in range(3):
        app.processEvents()
    table, model = dlg.table, dlg.table.model()
    assert table.rowCount() > 0 and not table.verticalHeader().isVisible()
    fm = table.fontMetrics()
    for c in range(model.columnCount()):
        need = max([fm.horizontalAdvance(table.horizontalHeaderItem(c).text()) + 24]
                   + [fm.horizontalAdvance(table.item(r, c).text()) + 16 for r in range(table.rowCount())])
        _fits(table, model, c, need)
    dlg.stop_auto(quiet=True)
    dlg.close()


def test_enter_in_a_field_presses_no_stray_button(app, win, monkeypatch):
    """With no button box the first push button was the default: Enter in the ORF length box added the selected ORFs as CDS
    features; in the Digest's name box it added a lane."""
    from tools import DigestDialog, OrfDialog
    dlg = OrfDialog(win, win.editor)
    dlg.min_aa.setValue(30)
    dlg.reload()
    dlg.show()
    new = [r for r, link in enumerate(dlg.links) if link is None]
    assert new, "the demo has an ORF that is not annotated"
    dlg.table.selectRow(new[0])
    assert dlg.add_button.isEnabled()
    before = len(win.editor.features)
    dlg.min_aa.setFocus()
    QTest.keyClick(dlg.min_aa, Qt.Key_Return)
    assert len(win.editor.features) == before, "Enter added CDS features"
    dlg.table.clearSelection()
    assert not dlg.add_button.isEnabled(), "nothing selected: nothing to add"
    dlg.close()

    digest = DigestDialog(win, win.editor)
    digest.show()
    digest.filter.setFocus()
    QTest.keyClick(digest.filter, Qt.Key_Return)
    assert len(digest.lanes) == 1, "Enter added a lane"
    digest.close()


def test_one_enter_is_one_ncbi_search(app, monkeypatch):
    from ncbi_dialog import GenomeDialog, NcbiDialog
    for cls in (NcbiDialog, GenomeDialog):
        d = cls(None)
        sent = []
        monkeypatch.setattr(d, "_run", lambda *a, **k: sent.append(a))
        d.show()
        d.query.setText("pUC19")
        d.query.setFocus()
        QTest.keyClick(d.query, Qt.Key_Return)
        assert len(sent) == 1, f"{cls.__name__}: {len(sent)} requests for one Enter"
        d.close()


def test_ncbi_error_text_is_not_read_as_markup(app):
    from ncbi_dialog import NcbiDialog
    d = NcbiDialog(None)
    d.token = token = object()
    d._done(token, None, "bad <b>reply</b>")
    assert "&lt;b&gt;" in d.status.text()
    d.close()


def test_the_digest_plus_sits_right_after_the_last_lane(app, win):
    from tools import DigestDialog
    dlg = DigestDialog(win, win.editor)
    dlg.resize(1200, 600)
    dlg.show()
    app.processEvents()
    bar = dlg.tabs
    for _ in range(2):
        last = bar.tabRect(bar.count() - 1)
        assert dlg.add.parentWidget() is bar and dlg.add.isVisible()
        assert 0 <= dlg.add.geometry().left() - last.right() <= 12, "the + is not next to the last tab"
        assert dlg.add.toolTip() == "Add another digest (lane)"
        dlg.add.click()
        app.processEvents()
    assert len(dlg.lanes) == 3
    dlg.close()


def _primer_tab(app, win, monkeypatch):
    import db_tabs
    import primer_db
    primer_db.add_primers(None, [("tt_short", "ACGTACGTACGTAC"), ("tt_long_name_primer", "GATCGATCGATCGATCGATCGGATCCAT"),
                                 ("tt_x", "TTGACAGCTAGCTCAGTCCT")], "Imported")
    tab = db_tabs.open_tab(win, db_tabs.PRIMERS)
    win.resize(1800, 850)
    for _ in range(3):
        app.processEvents()
    return tab


def test_primer_database_columns_fit_and_rows_have_no_wasted_indent(app, win, monkeypatch):
    tab = _primer_tab(app, win, monkeypatch)
    tree = tab.view()
    model = tree.model()
    fm = tree.fontMetrics()
    root = tree.invisibleRootItem()
    items = [root.child(i) for i in range(root.childCount())]
    assert items and not tree.rootIsDecorated(), "no primer opens here: no column of room for a ▸"
    for c in range(model.columnCount()):
        need = max([fm.horizontalAdvance(tree.headerItem().text(c)) + 24] + [fm.horizontalAdvance(i.text(c)) + 16 for i in items])
        _fits(tree, model, c, need)
    seq = tree.headerItem().text(1)
    assert all(tree.header().sectionSize(1) >= fm.horizontalAdvance(i.text(1)) for i in items), f"{seq} cut off"


def test_primer_database_delete_key_and_buttons(app, win, monkeypatch):
    import primer_db
    from primer_panel import PrimerDialog
    tab = _primer_tab(app, win, monkeypatch)
    tree = tab.view()
    opened = []
    monkeypatch.setattr(PrimerDialog, "exec", lambda self: opened.append(self) or 0)
    tree.clearSelection()
    tree.setCurrentItem(None)
    app.processEvents()
    assert not tab.edit_button.isEnabled() and not tab.delete_button.isEnabled(), "nothing selected"
    item = next(tree.topLevelItem(i) for i in range(tree.topLevelItemCount()) if tree.topLevelItem(i).text(0) == "tt_x")
    tree.setCurrentItem(item)
    item.setSelected(True)
    assert tab.edit_button.isEnabled() and tab.delete_button.isEnabled()
    # the Delete shortcut needs the list focused in the active window (other tests leave windows open; activateWindow() is
    # asynchronous on Windows): wait until this window is the active one
    win.raise_()
    win.activateWindow()
    QTest.qWaitForWindowActive(win, 2000)   # (best effort: the Windows test machine does not let a window become active)
    tree.setFocus()
    QTest.keyClick(tree, Qt.Key_Return)
    assert not opened, "Enter in the list opened a primer dialog"
    asked = []
    monkeypatch.setattr(QMessageBox, "question", staticmethod(lambda *a, **k: asked.append(a) or QMessageBox.Yes))
    QTest.keyClick(tree, Qt.Key_Delete)
    assert asked and primer_db.db().find("tt_x") is None, "the Delete key deletes the primer selected"


def test_delete_questions_default_to_no(app, monkeypatch):
    """A destructive question has No as its default: Enter or a stray click on the default button keeps the primers."""
    import db_tabs
    seen = []

    def question(*args):
        seen.append(args)
        return args[-1]   # the default button
    monkeypatch.setattr(QMessageBox, "question", staticmethod(question))
    db_tabs.PrimersDatabase().delete(None, ["anything"])
    assert seen and seen[0][-1] == QMessageBox.No


def _python_steps(monkeypatch):
    """Run the dialog's steps with Python instead of bash (through WSL on Windows, which a CI machine does not have): what is
    tested is the dialog, not the shell. "sleep N" sleeps, "exit N" ends with that code."""
    import runner

    def command_line(text):
        sleep, code = re.search(r"sleep (\d+)", text), re.search(r"exit (\d+)", text)
        body = f"import time; time.sleep({sleep.group(1)})" if sleep else f"import sys; sys.exit({code.group(1) if code else 0})"
        return [sys.executable, "-c", body]
    monkeypatch.setattr(runner, "command_line", command_line)


def test_a_running_mapping_is_not_hidden_by_escape(app, monkeypatch):
    """Escape went to reject(), which skips closeEvent: the window hid while bwa / minimap2 kept running."""
    import runner
    _python_steps(monkeypatch)
    run = runner.RunDialog(None, "sleep", ["sleep 30"])
    run.show()
    QTest.qWait(300)
    assert run.running() and run.busy.isVisible()
    monkeypatch.setattr(QMessageBox, "question", staticmethod(lambda *a, **k: QMessageBox.No))
    QTest.keyClick(run, Qt.Key_Escape)
    assert run.isVisible() and run.running(), "No: it goes on"
    monkeypatch.setattr(QMessageBox, "question", staticmethod(lambda *a, **k: QMessageBox.Yes))
    QTest.keyClick(run, Qt.Key_Escape)
    assert not run.isVisible() and not run.running(), "Yes: stopped, then closed"


def test_a_failed_step_says_its_exit_code(app, monkeypatch):
    import runner
    _python_steps(monkeypatch)
    run = runner.RunDialog(None, "fail", ["exit 3"])
    for _ in range(100):
        if not run.running():
            break
        QTest.qWait(50)
    QTest.qWait(100)
    assert run.status.text() == "Failed (exit code 3)."
    run.close()
