"""Small pieces every panel of the interface shares, so the same thing looks and behaves the same everywhere: counted words
("3 primers", never "3 primer(s)"), tables whose columns fit what they hold, and buttons that Enter does not press by surprise."""
import threading
import time

import shiboken6

from PySide6.QtCore import QEvent, QEventLoop, QModelIndex, QObject, Qt, QTimer, Signal
from PySide6.QtWidgets import (
    QApplication, QHeaderView, QMessageBox, QProgressDialog, QPushButton, QStyle, QStyleOptionViewItem, QTableView, QTreeView,
)

# Words whose plural is not the word + "s".
_IRREGULAR = {"copy": "copies", "match": "matches", "mismatch": "mismatches", "ORF": "ORFs", "is": "are", "was": "were",
              "this": "these", "it": "they", "has": "have", "entry": "entries", "library": "libraries", "sequence": "sequences",
              "assembly": "assemblies", "index": "indexes", "bp": "bp", "nt": "nt", "aa": "aa", "kb": "kb"}


class _TallBoxes(QObject):
    """A message box whose text would be taller than the screen keeps its first lines and puts the whole text in the
    scrollable details area (opened), so its buttons stay reachable."""
    HEAD = 4

    def eventFilter(self, obj, event):
        if event.type() == QEvent.Show and isinstance(obj, QMessageBox) and not obj.detailedText():
            text = obj.text()
            screen = obj.screen() or QApplication.primaryScreen()
            lines = text.count("\n") + 1 + len(text) // 90   # (a rough count of the lines it takes)
            if screen is not None and lines * obj.fontMetrics().lineSpacing() > screen.availableGeometry().height() * 0.6:
                obj.setDetailedText(text)
                obj.setText("\n".join(text.split("\n")[:self.HEAD]) + "\n…")
                for b in obj.findChildren(QPushButton):
                    if "etails" in b.text():
                        b.click()
                        break
        return False


def limit_tall_boxes(app):
    """Install once on the QApplication: every message box then fits the screen."""
    app._tall_boxes = _TallBoxes(app)
    app.installEventFilter(app._tall_boxes)


def plural(n, word, many=None):
    """'1 primer', '3 primers', '1,234 features' (thousands separated); `many` for a plural the table above does not know."""
    if n == 1:
        form = word
    else:
        form = many or _IRREGULAR.get(word) or (word[:-1] + "ies" if word.endswith("y") and word[-2:-1] not in "aeiou"
                                               else word + "es" if word.endswith(("s", "x", "ch", "sh")) else word + "s")
    return f"{n:,} {form}"


def no_auto_default(*widgets):
    """Enter in a field of a dialog presses its default button, and a dialog with no button box makes the first push button the
    default: Enter in a filter box opened "New…", added a lane, added every ticked feature. Called on the dialog (or any parent),
    no push button below it is pressed by Enter unless one is made the default on purpose afterwards."""
    for w in widgets:
        buttons = [w] if isinstance(w, QPushButton) else w.findChildren(QPushButton)
        for b in buttons:
            b.setAutoDefault(False)
            b.setDefault(False)


def _column_settings():
    from primer_db import window_state
    return window_state()


class _Fitter(QObject):
    """Keeps a table's columns fitted (see fit_columns): again when the rows change or the table is resized, never undoing a width
    the user dragged by hand. A dragged width is remembered for the next session (per table, known by its column names);
    double-clicking a column's edge gives the column back to the fitting."""

    def __init__(self, view, stretch, cap, min_stretch):
        super().__init__(view)
        self.view, self.stretch, self.cap, self.min_stretch = view, stretch, cap, min_stretch
        self.by_hand = set()
        self.busy = False
        self.header.setStretchLastSection(False)   # (a tree's header stretches its last column by default: not a drag by the user)
        self.timer = QTimer(self, singleShot=True, interval=0, timeout=self.fit)
        view.viewport().installEventFilter(self)
        self.header.sectionResized.connect(self._resized)
        self.header.sectionHandleDoubleClicked.connect(self._unpin)
        self.loaded = False
        self.watch(view.model())

    @property
    def header(self):
        return self.view.horizontalHeader() if isinstance(self.view, QTableView) else self.view.header()

    def watch(self, model):
        if model is None:
            return
        for signal in (model.modelReset, model.rowsInserted, model.rowsRemoved, model.layoutChanged, model.dataChanged):
            signal.connect(self._changed)

    UNSIZED_ROLES = {Qt.ForegroundRole, Qt.BackgroundRole, Qt.ToolTipRole, Qt.StatusTipRole, Qt.CheckStateRole}

    def _changed(self, *args):
        """The rows changed: fit again, unless the change cannot need it: only colours changed (dataChanged with those roles: a
        feature hidden, a display switch), or the model says so (`fit_quiet`: a gene opened or closed in the Features view)."""
        roles = args[2] if len(args) == 3 and isinstance(args[2], (list, tuple)) else None   # (dataChanged; rowsInserted has ints)
        if roles and all(r in self.UNSIZED_ROLES for r in roles):
            return
        try:
            quiet = getattr(self.view.model(), "fit_quiet", False)
        except RuntimeError:   # the table is being deleted: its own model (a QTableWidget's) clears its rows on the way out
            return
        if not quiet and not self.timer.isActive():   # (restarting it for every cell of a big refill cost seconds: one fit follows anyway)
            self.timer.start()

    def _names(self):
        model = self.view.model()
        clean = lambda c: str(model.headerData(c, Qt.Horizontal) or "").rstrip(" ▾▴*•·").strip() or str(c)
        return [clean(c) for c in range(model.columnCount())] if model is not None else []

    def _key(self):
        return "columns/" + "|".join(self._names()).replace("/", "-")

    def _resized(self, column, _old, new):
        if not self.busy:   # dragged by hand: that width is the user's from now on, and next time too
            self.by_hand.add(column)
            names = self._names()
            if column < len(names):
                _column_settings().setValue(f"{self._key()}/{names[column]}", new)

    def _unpin(self, column):
        """Double-click on a column's edge: back to fitting what it holds (and forgotten)."""
        self.by_hand.discard(column)
        names = self._names()
        if column < len(names):
            _column_settings().remove(f"{self._key()}/{names[column]}")
        self.timer.start()

    def _load(self):
        """The widths dragged by hand in an earlier session, put back once the table has its columns."""
        names = self._names()
        if not names:
            return
        self.loaded = True
        s = _column_settings()
        for c, name in enumerate(names):
            width = s.value(f"{self._key()}/{name}")
            if width is not None:
                try:
                    self.header.resizeSection(c, max(30, int(width)))
                    self.by_hand.add(c)
                except (TypeError, ValueError):
                    pass

    def eventFilter(self, obj, event):
        if event.type() == QEvent.Resize:
            self.timer.start()
        return False

    ARROW = 14   # room for a sort arrow beside a header's text
    tight = frozenset()   # columns as wide as their cells only (a checkbox, a short number): the header's text is cut

    def content_width(self, column):
        """Its widest cell, or its header's text with the style's margins and room for a sort arrow (and a filter mark).
        (Not Qt's own measure of the header section: with sorting on, it gives every column the room of the sort arrow at the
        style's size, which the Windows style makes large, so every column came out tens of pixels wider than it needs.)"""
        v, h = self.view, self.header
        if column in self.tight:
            return self.cells_width(column) + 8
        margin = v.style().pixelMetric(QStyle.PM_HeaderMargin, None, h)
        text = str(v.model().headerData(column, Qt.Horizontal) or "")
        title = h.fontMetrics().horizontalAdvance(text) + 2 * margin + self.ARROW
        return max(title, self.cells_width(column) + 8)

    MEASURED_ROWS = 2000   # up to this many rows, every row is measured; above, Qt's measure of the rows in sight

    def cells_width(self, column):
        """The widest cell of a column, every row counted. Qt's own sizeHintForColumn only measures the rows in sight, so in a
        short window (a laptop, the Windows test machine's screen) a long name further down was cut off. A view that measures
        its column itself (a sizeHintForColumn of its own, e.g. one that skips spanned group rows) keeps its way."""
        v = self.view
        model = v.model()
        own = any("sizeHintForColumn" in vars(k) for k in type(v).__mro__ if not k.__module__.startswith(("PySide6", "Shiboken")))
        if own or model.rowCount() > self.MEASURED_ROWS:
            return v.sizeHintForColumn(column)
        tree = isinstance(v, QTreeView)
        option = QStyleOptionViewItem()
        v.initViewItemOption(option)
        delegate = v.itemDelegateForColumn(column) or v.itemDelegate()
        widest = 0
        for r in range(model.rowCount()):
            if tree:
                if v.isRowHidden(r, QModelIndex()) or v.isFirstColumnSpanned(r, QModelIndex()):
                    continue
            elif v.isRowHidden(r) or v.columnSpan(r, column) > 1:
                continue
            widest = max(widest, delegate.sizeHint(option, model.index(r, column)).width())
        if tree and column == 0 and v.rootIsDecorated():
            widest += v.indentation()
        return widest

    def fit(self):
        v, h = self.view, self.header
        model = v.model()
        if model is None or not v.isVisible():
            return
        self.busy = True
        try:
            if not self.loaded:
                self._load()
            h.setStretchLastSection(False)
            columns = [c for c in range(model.columnCount()) if not h.isSectionHidden(c)]
            stretch = self.stretch if self.stretch in columns else None
            used = 0
            for c in columns:
                h.setSectionResizeMode(c, QHeaderView.Interactive)
                if c == stretch:
                    continue
                if c not in self.by_hand:
                    h.resizeSection(c, min(self.content_width(c), self.cap))
                used += h.sectionSize(c)
            if stretch is not None and stretch not in self.by_hand:
                room = v.viewport().width() - used
                # the descriptive column takes the room left, up to what it holds (any rest stays empty at the right: no column
                # is spread wider than its content), and never less than min_stretch (then the table scrolls sideways)
                h.resizeSection(stretch, max(self.min_stretch, min(room, self.content_width(stretch))))
        finally:
            self.busy = False


def fit_columns(view, stretch=None, cap=320, min_stretch=120, row_numbers=False):
    """Columns as wide as what they hold (header included, at most `cap` px, longer text is elided by the view), the one
    descriptive column `stretch` (a name, a note, a source) taking the room left up to its own content, and no column spread
    wider than its content: what is left stays empty at the right, instead of a "Strand" column 800 px wide in a tab of its own.
    A width the user drags is kept. Rows are as tall as the font needs; Qt's row numbers are hidden unless `row_numbers`."""
    if isinstance(view, QTableView):
        rows = view.verticalHeader()
        rows.setDefaultSectionSize(view.fontMetrics().height() + 8)
        rows.setVisible(row_numbers)
    elif hasattr(view, "setUniformRowHeights"):   # a tree: rows of one height too (and quicker to lay out)
        view.setUniformRowHeights(True)
    view.setTextElideMode(Qt.ElideRight)
    fitter = getattr(view, "_fitter", None)
    if fitter is None:
        fitter = view._fitter = _Fitter(view, stretch, cap, min_stretch)
        # a model set later (setModel) is watched too
        set_model = view.setModel

        def setModel(model, _set=set_model, _view=view):
            _set(model)
            _view._fitter.watch(model)
            _view._fitter.timer.start()
        view.setModel = setModel
    else:
        fitter.stretch, fitter.cap, fitter.min_stretch = stretch, cap, min_stretch
    fitter.timer.start()
    return fitter


def right_aligned(item):
    """A table item holding a number: right-aligned (digits line up), vertically centred."""
    item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
    return item


class _BusyDialog(QProgressDialog):
    """The "busy" dialog of run_busy: neither Escape nor its close button hides it while the work goes on (hiding it gave the
    window back to the user in the middle of a save or a load). Once the work is done it closes like any window, so even a dialog
    that came back somehow can never stay on screen for ever."""
    working = True

    def reject(self):
        if not self.working:
            super().reject()

    def closeEvent(self, event):
        if self.working:
            event.ignore()
        else:
            super().closeEvent(event)


def dispose_progress(dialog):
    """Put a progress dialog away for good, now (run_busy's busy dialog, the background genome search).

    QProgressDialog shows itself from a timer of its own (minimumDuration), and hide() / close() do not stop that timer;
    deleteLater() is not enough either, as Qt does not delete an object while a dialog opened after it runs its own event loop
    (the next question, the Import as New Chromosome dialog...). The timer then fired inside that dialog and showed the busy
    dialog again, which ignores Escape and its close button: it stayed on screen for ever, over a window that no longer answered.
    reset() stops the timer, and the dialog is deleted at once."""
    dialog.working = False   # (_BusyDialog: it may close now)
    dialog.reset()
    dialog.hide()
    try:
        shiboken6.delete(dialog)
    except RuntimeError:   # already deleted with its window
        pass


def show_after(dialog, ms):
    """A progress dialog with no steps (a busy bar) shows itself after `ms`. setMinimumDuration alone does not do it: Qt starts
    a 4 s show timer when the dialog is made and restarts it only once the value is the minimum, which a busy bar's (-1) is not;
    setValue(0) restarts it with `ms`."""
    dialog.setMinimumDuration(ms)
    dialog.setValue(dialog.minimum())


def run_cancellable(parent, title, work, steps=0):
    """Run `work(state)` in a thread with a progress dialog that has a Cancel button (shown after 400 ms), keeping the window
    responsive. `work` reads `state["cancel"]` to stop early and may set `state["text"]` and `state["done"]` (out of `steps`; 0
    steps: a busy bar). Returns what `work` returned, or re-raises its error."""
    state = {"text": title + "…", "cancel": False, "done": 0}
    result = {}

    def target():
        try:
            result["value"] = work(state)
        except BaseException as e:   # handed back to the GUI thread below
            result["error"] = e

    thread = threading.Thread(target=target, daemon=True)
    thread.start()
    dialog = QProgressDialog(state["text"], "Cancel", 0, steps, parent)
    dialog.setWindowTitle(title)
    dialog.setWindowModality(Qt.WindowModal)
    show_after(dialog, 400)
    dialog.canceled.connect(lambda: state.update(cancel=True))
    while thread.is_alive():
        QApplication.processEvents(QEventLoop.AllEvents, 30)
        dialog.setLabelText(state["text"])
        if steps:
            dialog.setValue(state["done"])
        thread.join(0.02)
    dispose_progress(dialog)   # (close() left its show timer running: it could come back over the next dialog)
    if "error" in result:
        raise result["error"]
    return result.get("value")


BUSY_DELAY_MS = 1500   # run_busy's dialog shows only for work longer than this: a shorter one (opening a chromosome) would just flash


def run_busy(parent, text, work):
    """Run `work()` in a thread with a "busy" dialog (shown only if it takes more than BUSY_DELAY_MS), keeping the window responsive. Returns its
    result or re-raises its error. The Rust calls used here release Python's lock, so the GUI really keeps running.

    Until the dialog shows, clicks and keys are held back (repaints still happen): otherwise a click in that first moment reached
    the window and could start a second save, open or close in the middle of this one."""
    result = {}

    def target():
        try:
            result["value"] = work()
        except BaseException as e:  # handed back to the GUI thread below
            result["error"] = e

    thread = threading.Thread(target=target, daemon=True)
    thread.start()
    dialog = _BusyDialog(text, None, 0, 0, parent)
    dialog.setWindowTitle(text.rstrip("…. ") or "Working")
    dialog.setWindowModality(Qt.WindowModal)
    show_after(dialog, BUSY_DELAY_MS)
    if QApplication.overrideCursor() is None:
        QApplication.setOverrideCursor(Qt.WaitCursor)
        restore = True
    else:
        restore = False
    try:
        while thread.is_alive():
            flags = QEventLoop.AllEvents if dialog.isVisible() else QEventLoop.ExcludeUserInputEvents
            QApplication.processEvents(flags, 30)
            thread.join(0.02)
    finally:
        if restore:
            QApplication.restoreOverrideCursor()
        dispose_progress(dialog)
    if "error" in result:
        raise result["error"]
    return result["value"]


class _Delivery(QObject):
    """Hands a background job's answer to the GUI thread (a queued signal). Lives forever, so a late answer never hits a deleted object."""
    done = Signal(object, object, object)   # done callback, value, error

    def __init__(self):
        super().__init__()
        self.done.connect(self._call)   # a slot of this object (made on the GUI thread): queued, so it runs there

    def _call(self, done, value, error):
        if error is not None:
            raise error
        done(value)


_delivery = _Delivery()
_jobs = {}   # job -> (label, time.monotonic() at its start): what the status bar says is running


SYNC_LIMIT = 50_000   # bp: on a longer sequence the analyses that follow every edit run in a background thread (`Latest`)


def in_background(label, work, done):
    """Run `work()` in a worker thread without blocking anything; then `done(value)` on the GUI thread (an exception `work` raised is
    raised there instead). `label` is what the status bar says while it runs (`background_jobs`). Returns the thread (`.join()` waits)."""
    job = object()

    def target():
        value = error = None
        try:
            value = work()
        except Exception as e:   # raised on the GUI thread
            error = e
        finally:
            _jobs.pop(job, None)
        try:
            _delivery.done.emit(done, value, error)
        except RuntimeError:   # the program is closing
            pass

    _jobs[job] = (label, time.monotonic())
    thread = threading.Thread(target=target, daemon=True)
    thread.start()
    return thread


def background_jobs(after=1.0):
    """The labels of the `in_background` jobs running for more than `after` seconds, oldest first, each once."""
    now = time.monotonic()
    return list(dict.fromkeys(label for label, start in list(_jobs.values()) if now - start > after))   # (a dict keeps its start order)


class Latest:
    """An analysis that follows the document (`run` after every edit): now on a short sequence, else `in_background`, one job at a time.
    A request made while one runs waits, only the newest of them; `take(value)` gets an answer only if nothing was asked since (nor
    `cancel`), and only while `owner` (a QObject) exists. So a stale answer is never shown, and fast typing starts no pile of threads.
    An error is raised only from the current request (an older one may have read a document being edited)."""

    def __init__(self, owner, label):
        self.owner, self.label = owner, label
        self.asked, self.running, self.waiting = 0, False, None

    def run(self, work, take, now=False):
        """`now`: on the GUI thread, at once (a short sequence: see SYNC_LIMIT). `work` must not read what an edit changes (give it copies), or
        `take` must check that the document is the one it read (ToolDialog.compute)."""
        self.cancel()
        if now:
            take(work())
        elif self.running:
            self.waiting = (work, take)
        else:
            self._start(work, take)

    def cancel(self):
        self.asked += 1
        self.waiting = None

    def _start(self, work, take):
        self.running, ticket = True, self.asked

        def guarded():
            try:
                return work(), None
            except Exception as e:
                return None, e

        def done(answer):
            self.running = False
            if self.waiting:
                self._start(*self.waiting)
                self.waiting = None
            elif ticket == self.asked and shiboken6.isValid(self.owner):
                value, error = answer
                if error is not None:
                    raise error
                take(value)
        in_background(self.label, guarded, done)


def ask(parent, title, text, choices, default=None, detail=""):
    """A question whose buttons say what they do ("Delete Gene and Parts", "Delete Gene Only", Cancel) instead of Yes / No /
    Cancel with their meaning spelled out in the text. `choices`: [(key, button text)]; returns the key clicked, None for Cancel
    (and Escape). `default`: the key of the button Enter presses (None: Cancel)."""
    box = QMessageBox(QMessageBox.Question, title, text, parent=parent)
    if detail:
        box.setInformativeText(detail)
    buttons = {box.addButton(label, QMessageBox.AcceptRole): key for key, label in choices}
    cancel = box.addButton(QMessageBox.Cancel)
    box.setEscapeButton(cancel)
    box.setDefaultButton(next((b for b, k in buttons.items() if k == default), cancel))
    box.exec()
    return buttons.get(box.clickedButton())
