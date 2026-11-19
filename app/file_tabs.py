"""File tabs: the files of one window.

A window is a FileHost; each file in it is a *page* (a main.MainWindow embedded in the host, with its own editor, view tabs,
side panels, search box and timers: nothing is shared between pages, so one file's input can never reach another). FileSet
is the bookkeeping, without any widget: which pages, in which order, which one is shown, which are unsaved. Each page draws
the strip of file tabs itself (FileStrip) from the host's FileSet, so no widget is shared either.

The program-wide choice "open files in tabs, or in windows of their own" is `open_in_tabs()`.
"""
import os
from pathlib import Path

from PySide6.QtCore import QEvent, QObject, QSettings, QSize, Qt, QTimer, QUrl, Signal
from PySide6.QtGui import QDesktopServices, QGuiApplication
from PySide6.QtWidgets import QMenu, QSizePolicy, QStackedLayout, QTabBar, QVBoxLayout, QWidget

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


_last = None   # the window that was in front last


def last_host():
    """The window in use last (the newest one if none was ever in front), None with no window."""
    return _last if _last in HOSTS else HOSTS[-1] if HOSTS else None


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
_in_tabs = None   # open_in_tabs(), read from the ini file once


def open_in_tabs():
    """View > Open Files in Tabs: a file opens as a tab of the window (else in a window of its own). The same for every
    window, kept between sessions; on unless switched off."""
    global _in_tabs
    if _in_tabs is None:
        _in_tabs = str(window_settings().value("open_in_tabs", "true")).lower() == "true"
    return _in_tabs


def set_open_in_tabs(on):
    global _in_tabs
    _in_tabs = bool(on)
    window_settings().setValue("open_in_tabs", _in_tabs)
    mode.changed.emit()


MIN_CHARS = 12   # the least of a file's name a tab shows, before the tabs go on a second row (then scroll)


class _Row(QTabBar):
    """One row of the file tabs (FileStrip has two; the second shows only when the tabs do not fit in the first)."""

    def __init__(self, strip):
        super().__init__()
        self.strip = strip
        self.setTabsClosable(True)
        self.setMovable(True)
        self.setExpanding(False)
        self.setDocumentMode(True)
        self.setElideMode(Qt.ElideRight)
        self.setUsesScrollButtons(True)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.currentChanged.connect(lambda i: strip.picked(self, i))
        self.tabCloseRequested.connect(lambda i: strip.closed(self, i))
        self.tabMoved.connect(lambda a, b: strip.moved(self, a, b))

    def min_width(self, text, i=0):
        """Width of a tab showing the first MIN_CHARS letters of `text` (and "…"), like the room around the text in tab `i`."""
        fm = self.fontMetrics()
        chrome = self.tabSizeHint(i).width() - fm.horizontalAdvance(self.tabText(i))
        return chrome + fm.horizontalAdvance(text if len(text) <= MIN_CHARS else text[:MIN_CHARS] + "…")

    def minimumTabSizeHint(self, i):   # Qt shrinks the tabs to fit down to this (then it scrolls)
        return QSize(self.min_width(self.tabText(i), i), super().minimumTabSizeHint(i).height())

    def sizeHint(self):
        # No wish of its own: the toolbar gives it what the search box and the range box leave (it would otherwise push them into
        # the overflow » of a narrow window before the tabs themselves do)
        return QSize(1, super().sizeHint().height())

    def minimumSizeHint(self):
        return QSize(1, super().minimumSizeHint().height())

    def contextMenuEvent(self, event):
        self.strip.context(self, event)


class FileStrip(QWidget):
    """The file tabs as one page shows them. It takes the room the toolbar has left (never more than the tabs need), so the search
    box stays at the right. Each tab keeps at least MIN_CHARS letters of its name; when the tabs do not fit like that they go on
    two rows, and only past that do they scroll."""

    def __init__(self, page):
        super().__init__()
        self.page = page
        self.rows = [_Row(self), _Row(self)]
        self.split = 0   # tabs in the first row
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        for row in self.rows:
            layout.addWidget(row)
        self.rows[1].hide()
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)

    def sizeHint(self):
        return QSize(1, self.rows[0].sizeHint().height() * (2 if self.rows[1].isVisible() else 1))

    def minimumSizeHint(self):
        return self.sizeHint()

    def _titles(self):
        return [p.tab_title() for p in self.page.host.files.pages]

    def _split(self, titles):
        """How many tabs go in the first row: all of them when they fit at MIN_CHARS letters, else half (or what fits, if less)."""
        row, room, used, fits = self.rows[0], self.width() - 70, 0, 0   # (70: the scroll arrows)
        if row.count() == 0:
            return len(titles)
        for t in titles:
            used += row.min_width(t)
            if used > room:
                break
            fits += 1
        return len(titles) if fits >= len(titles) else max(1, min(fits, (len(titles) + 1) // 2))

    def _shown(self):
        return [r.tabText(i) for r in self.rows for i in range(r.count())]

    def sync(self):
        """The tabs the host's FileSet says, the shown one current. The tabs are only rebuilt when the files changed: a click on a
        tab comes here from a row's currentChanged, which must not tear the tabs down under the mouse."""
        host = self.page.host
        pages, titles = host.files.pages, self._titles()
        split = self._split(titles)
        for r in self.rows:
            r.blockSignals(True)
        if self._shown() != titles or split != self.split:
            for r in self.rows:
                while r.count():
                    r.removeTab(0)
            for i, p in enumerate(pages):
                r = self.rows[0 if i < split else 1]
                r.setTabToolTip(r.addTab(titles[i]), str(p.file_path() or titles[i]))
            self.split = split
        two = split < len(pages)
        if self.rows[1].isVisible() != two:
            self.rows[1].setVisible(two)
            self.updateGeometry()
        now = pages.index(host.files.current)
        self.rows[0].setCurrentIndex(now if now < split else -1)
        self.rows[1].setCurrentIndex(now - split if now >= split else -1)
        for r in self.rows:
            r.blockSignals(False)

    def resizeEvent(self, event):
        super().resizeEvent(event)
        host = getattr(self.page, "host", None)   # (not yet set while the page is being built)
        if host is not None and host.files.pages and self._split(self._titles()) != self.split:
            self.sync()

    def retitle(self):
        """A file's name or "unsaved" mark changed: only the text of the tabs, not the tabs (an edit does this at typing speed)."""
        titles = self._titles()
        for i, title in enumerate(titles):
            row, j = self._at(i)
            if row.tabText(j) != title:
                row.setTabText(j, title)
        if self._split(titles) != self.split:
            self.sync()

    def _at(self, i):
        """(row, index in it) of the file tab `i`."""
        return (self.rows[0], i) if i < self.split else (self.rows[1], i - self.split)

    def _global(self, row, i):
        return i if row is self.rows[0] else i + self.split

    def picked(self, row, i):
        if i >= 0:
            self.page.host.select(self.page.host.files.pages[self._global(row, i)])

    def closed(self, row, i):
        self.page.host.close_page(self.page.host.files.pages[self._global(row, i)])

    def moved(self, row, a, b):
        self.page.host.moved(self, self._global(row, a), self._global(row, b))

    def context(self, row, event):
        i = row.tabAt(event.pos())
        if i < 0:
            return
        host = self.page.host
        page = host.files.pages[self._global(row, i)]
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
        for p in self.files.pages:   # (the tab's name and its "unsaved" mark follow the title)
            p.file_strip.retitle()

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
        if ask and not self._confirm(page):
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

    def _confirm(self, page):
        """False if the user cancels closing `page`: an unsaved file is shown first, so the question is asked over it."""
        if page.has_unsaved():
            self.select(page)
        return page.confirm_discard()

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

    def changeEvent(self, event):
        if event.type() == QEvent.ActivationChange and self.isActiveWindow():
            global _last
            _last = self
        super().changeEvent(event)

    def closeEvent(self, event):
        pages = list(self.files.pages)
        if not all(self._confirm(p) for p in pages):   # each unsaved file asks; one Cancel keeps the window and every file
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
