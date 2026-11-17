"""File tabs: the files of one window.

A window is a FileHost; each file in it is a *page* (a main.MainWindow embedded in the host, with its own editor, view tabs,
side panels, search box and timers: nothing is shared between pages, so one file's input can never reach another). FileSet
is the bookkeeping, without any widget: which pages, in which order, which one is shown, which are unsaved. Each page draws
the strip of file tabs itself (FileStrip) from the host's FileSet, so no widget is shared either.

The program-wide choice "open files in tabs, or in windows of their own" is `open_in_tabs()`.
"""
import os
from pathlib import Path

from PySide6.QtCore import QObject, QSettings, QSize, Qt, QTimer, QUrl, Signal
from PySide6.QtGui import QDesktopServices, QGuiApplication
from PySide6.QtWidgets import QMenu, QSizePolicy, QStackedLayout, QTabBar, QWidget

from primer_db import data_dir


class FileSet:
    """Pages in tab order and the shown one. A page is any object with `file_path()` (a Path, None while untitled) and
    `has_unsaved()`."""

    def __init__(self):
        self.pages, self.current = [], None

    def add(self, page):
        """Append `page` and show it."""
        self.pages.append(page)
        self.current = page
        return page

    def find(self, path):
        """The page showing the file at `path`, or None: the same file whatever the spelling (case on Windows, relative path)."""
        key = os.path.normcase(os.path.realpath(path))
        return next((p for p in self.pages if p.file_path() is not None
                     and os.path.normcase(os.path.realpath(p.file_path())) == key), None)

    def close(self, page):
        """Remove `page`. If it was shown, the page that took its place (the one to its right, else to its left) is."""
        i = self.pages.index(page)
        self.pages.remove(page)
        if self.current is page:
            self.current = self.pages[min(i, len(self.pages) - 1)] if self.pages else None
        return self.current

    def move(self, src, dst):
        """Drag a tab from position `src` to `dst`."""
        self.pages.insert(dst, self.pages.pop(src))

    def step(self, by):
        """The page `by` places after the shown one, round the end."""
        return self.pages[(self.pages.index(self.current) + by) % len(self.pages)]

    def unsaved(self):
        return [p for p in self.pages if p.has_unsaved()]


HOSTS = []   # every open window: keeps them referenced so they stay open


def all_pages():
    """Every file page of every window."""
    return [p for h in HOSTS for p in h.files.pages]


def find_open(path):
    """The page of any window that shows the file at `path`, or None."""
    return next((p for p in (h.files.find(path) for h in HOSTS) if p is not None), None)


def page_of(widget):
    """The file page (the main.MainWindow) a widget is inside: its window() is the host, which has none of a page's attributes."""
    while widget is not None and not getattr(widget, "is_file_page", False):
        widget = widget.parentWidget()
    return widget


def window_settings():
    return QSettings(str(data_dir() / "window.ini"), QSettings.IniFormat)


class _Mode(QObject):
    changed = Signal()


mode = _Mode()   # told when open_in_tabs() changes, so every window shows the choice


def open_in_tabs():
    """View > Open Files in Tabs: a file opens as a tab of the window (else in a window of its own). The same for every
    window, kept between sessions; on unless switched off."""
    return str(window_settings().value("open_in_tabs", "true")).lower() == "true"


def set_open_in_tabs(on):
    window_settings().setValue("open_in_tabs", bool(on))
    mode.changed.emit()


class FileStrip(QTabBar):
    """The file tabs as one page shows them. It takes the room the toolbar has left (never more than the tabs need, and it
    scrolls when they do not fit), so the search box stays at the right."""

    def __init__(self, page):
        super().__init__()
        self.page = page
        self.setTabsClosable(True)
        self.setMovable(True)
        self.setExpanding(False)
        self.setDocumentMode(True)
        self.setElideMode(Qt.ElideMiddle)
        self.setUsesScrollButtons(True)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.currentChanged.connect(lambda i: self.page.host.select(self.page.host.files.pages[i]) if i >= 0 else None)
        self.tabCloseRequested.connect(lambda i: self.page.host.close_page(self.page.host.files.pages[i]))
        self.tabMoved.connect(lambda a, b: self.page.host.moved(self, a, b))

    def sizeHint(self):
        # No wish of its own: the toolbar gives it what the search box and the range box leave (it would otherwise push them into
        # the overflow » of a narrow window before the tabs themselves do)
        return QSize(1, super().sizeHint().height())

    def minimumSizeHint(self):
        return QSize(1, super().minimumSizeHint().height())

    def sync(self):
        """The tabs the host's FileSet says, the shown one current."""
        host = self.page.host
        self.blockSignals(True)
        while self.count():
            self.removeTab(0)
        for p in host.files.pages:
            i = self.addTab(p.tab_title())
            self.setTabToolTip(i, str(p.file_path() or p.tab_title()))
        self.setCurrentIndex(host.files.pages.index(host.files.current))
        self.blockSignals(False)

    def contextMenuEvent(self, event):
        i = self.tabAt(event.pos())
        if i < 0:
            return
        host = self.page.host
        page = host.files.pages[i]
        path = page.file_path()
        menu = QMenu(self)
        menu.addAction("Close", lambda: host.close_page(page))
        menu.addAction("Close Others", lambda: host.close_others(page)).setEnabled(len(host.files.pages) > 1)
        menu.addSeparator()
        menu.addAction("Copy Path", lambda: QGuiApplication.clipboard().setText(str(path))).setEnabled(path is not None)
        menu.addAction("Reveal in Folder", lambda: QDesktopServices.openUrl(QUrl.fromLocalFile(str(Path(path).parent)))
                       ).setEnabled(path is not None)
        menu.addSeparator()
        menu.addAction("Move to New Window", page.move_to_new_window).setEnabled(len(host.files.pages) > 1)
        menu.exec(event.globalPos())


class FileHost(QWidget):
    """The window itself: size and place (remembered between sessions), title, closing. It shows its current page.
    `make_page` makes an empty page (main.MainWindow)."""

    def __init__(self, make_page):
        super().__init__()
        self.make_page = make_page
        self.files = FileSet()
        self.stack = QStackedLayout(self)
        self.stack.setContentsMargins(0, 0, 0, 0)
        HOSTS.append(self)
        self._restore()
        mode.changed.connect(self.sync)

    # ---- the pages ----------------------------------------------------------------------------------------------------
    def add_page(self, page):
        self.files.add(page)
        page.host = self
        self.stack.addWidget(page)
        page.windowTitleChanged.connect(lambda title, page=page: self._retitled(page))
        self.select(page)

    def new_page(self):
        """A page with nothing in it, added after the others and shown."""
        page = self.make_page()
        self.add_page(page)
        return page

    def select(self, page):
        self.files.current = page
        self.stack.setCurrentWidget(page)
        self.setWindowTitle(page.windowTitle())
        self.sync()
        page.focus_view()

    def _retitled(self, page):
        if page is self.files.current:
            self.setWindowTitle(page.windowTitle())
        self.sync()   # (the tab's name and its "unsaved" mark follow the title)

    def sync(self):
        """Every page's strip shows the same tabs; the strips show when files are in tabs (or more than one is open)."""
        pages = self.files.pages
        shown = (open_in_tabs() or len(pages) > 1) and any(p.editor is not None for p in pages)
        for p in pages:
            p.file_strip.sync()
            p.file_strip.setVisible(shown)
            p.file_mode_changed()

    def moved(self, strip, src, dst):
        """A tab was dragged in `strip`: the order of the files follows, in the other pages' strips too."""
        self.files.move(src, dst)
        for p in self.files.pages:
            if p.file_strip is not strip:
                p.file_strip.sync()

    def cycle(self, by):
        if len(self.files.pages) > 1:
            self.select(self.files.step(by))

    def jump(self, n):
        """Ctrl+1..9: the n-th tab (the 9th is the last one, as in a browser)."""
        pages = self.files.pages
        i = len(pages) - 1 if n == 9 else n - 1
        if 0 <= i < len(pages):
            self.select(pages[i])

    def close_page(self, page, ask=True):
        """Close a file tab. False when the user cancels (unsaved changes). The last file's tab leaves an empty page (the empty
        window is shown); closing the empty page itself closes the window."""
        if page.editor is None and len(self.files.pages) == 1:
            self.close()
            return True
        if ask and page.has_unsaved():
            self.select(page)   # (so the question is asked over the file it is about)
        if ask and not page.confirm_discard():
            return False
        page.shut_down()
        was_shown = page is self.files.current
        self.stack.removeWidget(page)
        page.hide()
        page.deleteLater()
        after = self.files.close(page)
        if after is None:
            self.new_page()
        elif was_shown:
            self.select(after)
        else:
            self.sync()
        return True

    def close_others(self, page):
        for p in list(self.files.pages):
            if p is not page and not self.close_page(p):
                return

    # ---- the window ---------------------------------------------------------------------------------------------------
    def _restore(self):
        geometry = window_settings().value("geometry")
        if geometry is None or not self.restoreGeometry(geometry):
            # first start: most of the screen, centred (a fixed 1400 x 850 was taller than a 1366 x 768 laptop screen)
            screen = (self.screen() or QGuiApplication.primaryScreen()).availableGeometry()
            w, h = min(1400, int(screen.width() * 0.85)), min(850, int(screen.height() * 0.85))
            self.resize(w, h)
            self.move(screen.center().x() - w // 2, screen.center().y() - h // 2)

    def closeEvent(self, event):
        pages = list(self.files.pages)
        for p in pages:   # each unsaved file asks, shown first; one Cancel keeps the window and every file
            if p.has_unsaved():
                self.select(p)
            if not p.confirm_discard():
                event.ignore()
                return
        for p in pages:
            p.shut_down()
        s = window_settings()
        s.setValue("geometry", self.saveGeometry())   # size, place and maximised (restoreGeometry keeps it on a screen)
        s.remove("panel_width")   # (the side panels were once as wide as dragged: they all have one fixed width now)
        event.accept()
        # Not now: this object is still being closed. (The list's reference is what keeps a window alive.)
        QTimer.singleShot(0, lambda host=self: HOSTS.remove(host) if host in HOSTS else None)
