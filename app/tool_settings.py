"""Per-tool settings: each tool declares its parameters once (register), reads them with
value(tool, key), and Tools → Settings shows one page per view or tool, listed on the left in collapsible groups that follow
the new-tab picker (Views, Annotations, Tools) and a General group, each page with a form built from the declarations (a page
without options says so). Stored in settings.ini in the app data folder (tests use a temp folder).

A change emits notifier.changed(tool), so open tools can redo their work. Every setting is declared once (its default, limits
and tip live in the Param and nowhere else); a tool whose settings limit each other gives `register` a `check` that brings
them back in range, applied both on OK and on every read, so a hand-edited settings.ini cannot break the tool.
"""
import os
from contextlib import contextmanager
from types import SimpleNamespace

from PySide6.QtCore import QObject, QSettings, Qt, Signal
from PySide6.QtGui import QPalette
from PySide6.QtWidgets import (
    QCheckBox, QComboBox, QFileDialog, QFrame, QGroupBox, QLabel, QLineEdit, QDialog, QDialogButtonBox, QDoubleSpinBox,
    QFormLayout, QHBoxLayout, QMessageBox, QPushButton, QScrollArea, QSpinBox, QStackedWidget, QTreeWidget, QTreeWidgetItem,
    QVBoxLayout, QWidget,
)

from user_error import UserError

TOOLS = {}  # tool name -> [Param], in declaration order
CHECKS = {}  # tool name -> check(values): brings settings that limit each other back in range (see register)
# The name a tool's page shows when it differs from its key (the key stays: it is where users' saved values are stored)
LABELS = {"Detect Features": "Detect Common Features"}

GENERAL = "General"
GENERAL_PAGES = ["Genetic code", "Display", "History", "Updates", "Performance"]   # settings that belong to no view; the Display page is the hover delay
# A registered tool whose settings are shown on the page of a view (by default a tool has the page of its own name).
PAGE_OF = {"Alignments": "Alignment", "Primer checks": "Primers"}
# A page without options of its own that depends on another page's: it says so and leads there.
SEE_ALSO = {"Sequence": "Display", "Map": "Display", "Protein": "Genetic code", "Detect ORFs": "Genetic code",
            "Spacer Finder": "Primers"}
NO_OPTIONS = "No options for this tool."
CUSTOM = {}   # page name -> factory(parent) of a widget with save() -> bool and defaults(), shown above the page's forms


def register_page(page, factory):
    CUSTOM[page] = factory


def page_layout():
    """[(group, [page names])] as the Settings list shows them: the new-tab picker's groups, then General (and any tool whose
    page is nowhere else)."""
    import workspace
    groups = [(g, [n for n in names if n]) for g, names in workspace.grouped_views()]
    known = {n for _, names in groups for n in names} | set(GENERAL_PAGES)
    extra = [PAGE_OF.get(t, t) for t in TOOLS if PAGE_OF.get(t, t) not in known]
    return groups + [(GENERAL, GENERAL_PAGES + list(dict.fromkeys(extra)))]


def Param(key, label, default, minimum=None, maximum=None, step=None, suffix="", tip="", choices=None, info=None, browse=None):
    """One setting. Its type follows `default`: bool -> checkbox, int / float -> spin box,
    str with `choices` -> drop-down list, other str -> a line of text (with `browse`, a file filter: a Browse… button beside it).
    `info`: optional function of the chosen value giving a small read-only line
    shown under the drop-down (e.g. the stop and start codons of the selected genetic code)."""
    return SimpleNamespace(key=key, label=label, default=default, minimum=minimum, maximum=maximum,
                           step=step, suffix=suffix, tip=tip, choices=choices, info=info, browse=browse)


class _Notifier(QObject):
    changed = Signal(str)  # tool name


notifier = _Notifier()


def store():
    from primer_db import data_dir   # (here, not at the top: primer_db declares its settings with this module)
    return QSettings(str(data_dir() / "settings.ini"), QSettings.IniFormat)


# The layout of settings.ini. Bump it when a saved key is renamed or changes meaning, and convert the older layout in stamp_format.
FORMAT, FORMAT_KEY = 1, "openDNA/format"


def stamp_format():
    """Write the layout number at start (never lowering one a newer openDNA wrote)."""
    s = store()
    if not s.contains(FORMAT_KEY):
        s.setValue(FORMAT_KEY, FORMAT)


def register(tool, params, check=None):
    """Declare a tool's settings. `check(values)`: for settings that limit each other; it changes `values` ({key: value}) in
    place to bring them back in range and returns a sentence saying what it changed ("" when nothing), or raises UserError
    for a value it cannot fix (Settings then stays open)."""
    TOOLS[tool] = list(params)
    if check is not None:
        CHECKS[tool] = check
    _cache.pop(tool, None)


_cache = {}  # tool -> {key: value}: reading settings.ini takes ~ms, and drawing asks per base


_overrides = {}   # (tool, key) -> value of a run that was given one (never saved): see `override`


@contextmanager
def override(values):
    """For the length of the block, `value` answers `values` ({(tool, key): value}) instead of what is saved; nothing is written to the settings file."""
    before = {k: _overrides[k] for k in values if k in _overrides}
    _overrides.update(values)
    try:
        yield
    finally:
        for k in values:
            _overrides.pop(k, None)
        _overrides.update(before)


def value(tool, key):
    if (tool, key) in _overrides:
        return _overrides[tool, key]
    return values(tool)[key]


def values(tool):
    """{key: value} of every setting of `tool` as saved (each within its limits, and passed through the tool's check)."""
    if tool not in _cache:
        saved = store()
        got = {p.key: _read(saved, tool, p) for p in TOOLS[tool]}
        if tool in CHECKS:
            try:
                CHECKS[tool](got)
            except UserError:   # a value the check cannot fix (a missing file...): it is the reader's to handle
                pass
        _cache[tool] = got
    given = {k: v for (t, k), v in _overrides.items() if t == tool} if _overrides else None
    return {**_cache[tool], **given} if given else _cache[tool]


def _read(saved, tool, p):
    raw = saved.value(f"{tool}/{p.key}", p.default)
    # QSettings gives strings back from an .ini file: convert to the default's type.
    if isinstance(p.default, bool):
        return raw if isinstance(raw, bool) else str(raw).lower() in ("true", "1")
    if p.choices:
        return raw if raw in p.choices else p.default
    try:
        v = type(p.default)(raw)
    except (TypeError, ValueError):
        return p.default
    if isinstance(v, (int, float)):   # a hand-edited file stays within the limits the Settings box allows
        v = max(v, p.minimum) if p.minimum is not None else v
        v = min(v, p.maximum) if p.maximum is not None else v
    return v


PERFORMANCE = "Performance"
CPUS = os.cpu_count() or 1
DEFAULT_CORES = max(1, min(8, CPUS - 1))   # one core left free so the window and the rest of the PC stay responsive
register(PERFORMANCE, [
    Param("cores", "CPU cores openDNA may use", DEFAULT_CORES, 1, CPUS, suffix=f" of {CPUS}",
          tip="For the analyses that run on several cores at once (Local Align, Folder Search). The installer asks it too."),
])


def cores():
    """How many threads a multithreaded analysis may use (Settings > General > Performance)."""
    return value(PERFORMANCE, "cores")


def set_value(tool, key, v):
    """Save one setting (the caller emits notifier.changed(tool) once it has set all it changes)."""
    store().setValue(f"{tool}/{key}", v)
    _cache.pop(tool, None)


class Pages(QWidget):
    """The pages of the settings: the groups and their pages as a tree on the left (a click on a group opens or folds it), the
    chosen page on the right (the usual Preferences layout)."""

    def __init__(self):
        super().__init__()
        self.tree = QTreeWidget(headerHidden=True, indentation=14)
        self.stack = QStackedWidget()
        self.items = {}   # page name -> its tree item
        self.tree.currentItemChanged.connect(self._shown)
        self.tree.itemClicked.connect(lambda item, _c: item.parent() is None and item.setExpanded(not item.isExpanded()))
        row = QHBoxLayout(self)
        row.setContentsMargins(0, 0, 0, 0)
        row.addWidget(self.tree)
        row.addWidget(self.stack, 1)
        self.wide = 0

    def add_group(self, name):
        item = QTreeWidgetItem([name])
        font = item.font(0)
        font.setBold(True)
        item.setFont(0, font)
        item.setFlags(Qt.ItemIsEnabled)   # a heading: it folds, it is not a page
        self.tree.addTopLevelItem(item)
        return item

    def add_page(self, group, name, widget, label=None):
        item = QTreeWidgetItem([label or name])
        item.setData(0, Qt.UserRole, self.stack.addWidget(widget))
        group.addChild(item)
        self.items[name] = item
        # the list as wide as its longest name (plus the indent and room for a scroll bar), no wider
        self.wide = max(self.wide, self.tree.fontMetrics().horizontalAdvance(label or name) + 60)
        self.tree.setFixedWidth(self.wide)

    def _shown(self, item, _previous=None):
        if item is not None and item.parent() is not None:
            self.stack.setCurrentIndex(item.data(0, Qt.UserRole))

    def show_page(self, name):
        """Open `name`: its group unfolds, the other groups fold."""
        item = self.items.get(name)
        if item is None:
            return False
        for i in range(self.tree.topLevelItemCount()):
            self.tree.topLevelItem(i).setExpanded(self.tree.topLevelItem(i) is item.parent())
        self.tree.setCurrentItem(item)
        return True

    def widget(self, name):
        return self.stack.widget(self.items[name].data(0, Qt.UserRole))

    def currentWidget(self):
        return self.stack.currentWidget()

    def currentPage(self):
        item = self.tree.currentItem()
        return next((n for n, i in self.items.items() if i is item), None)


class SettingsDialog(QDialog):
    """One page per view or tool; `tool` (a page name, or the key of a registered tool) is the page shown first. OK saves and
    notifies."""

    def __init__(self, parent, tool=None):
        super().__init__(parent)
        self.setWindowTitle("Settings")
        self.resize(780, 580)
        self.tabs = Pages()
        self.widgets = {}  # (tool, key) -> (Param, widget)
        self.custom = {}   # page name -> its custom widget
        self.pages = {}    # page name -> the tools whose forms it shows
        for group, names in page_layout():
            heading = self.tabs.add_group(group)
            for page in names:
                self.pages[page] = [t for t in TOOLS if PAGE_OF.get(t, t) == page]
                self.tabs.add_page(heading, page, self._page(page), LABELS.get(page, page))
        self.tabs.show_page(PAGE_OF.get(tool, tool)) or self.tabs.show_page(next(iter(self.tabs.items)))
        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout = QVBoxLayout(self)
        layout.addWidget(self.tabs)
        layout.addWidget(buttons)

    def page_for(self, tool):
        """The page widget that shows the settings of the registered tool `tool`."""
        return self.tabs.widget(PAGE_OF.get(tool, tool))

    def _page(self, page):
        body = QWidget()
        column = QVBoxLayout(body)
        tools = self.pages[page]
        column.addWidget(QLabel(f"<b>{LABELS.get(page, page)}</b>"))
        if page in SEE_ALSO:
            see = SEE_ALSO[page]
            row = QHBoxLayout()
            row.addWidget(QLabel(f"This one uses the settings of {LABELS.get(see, see)}."))
            row.addWidget(QPushButton(f"Open {LABELS.get(see, see)}", clicked=lambda _=False, n=see: self.tabs.show_page(n)))
            row.addStretch()
            column.addLayout(row)
        if page in CUSTOM:
            self.custom[page] = CUSTOM[page](body)
            column.addWidget(self.custom[page])
        for name in tools:
            form = QFormLayout()
            form.setFieldGrowthPolicy(QFormLayout.FieldsStayAtSizeHint)   # a number box as wide as a number, not the whole page
            for p in TOOLS[name]:
                w = self.param_widget(p, value(name, p.key))
                w.setToolTip(p.tip)
                self.widgets[(name, p.key)] = (p, w)
                form.addRow(p.label, w)
                if p.info and isinstance(w, QComboBox):   # a small read-only line under the choice, kept in step with it
                    line = QLabel(p.info(w.currentText()), wordWrap=True)
                    line.setForegroundRole(QPalette.PlaceholderText)   # faint, in either theme
                    w.currentTextChanged.connect(lambda text, line=line, info=p.info: line.setText(info(text)))
                    form.addRow("", line)
            if len(tools) + (page in CUSTOM) > 1:   # several forms on one page: each under its own name
                box = QGroupBox(LABELS.get(name, name))
                box.setLayout(form)
                column.addWidget(box)
            else:
                column.addLayout(form)
        if not tools and page not in CUSTOM:
            note = QLabel(NO_OPTIONS)
            note.setForegroundRole(QPalette.PlaceholderText)
            column.addWidget(note)
        else:
            reset = QPushButton("Restore Defaults", toolTip=f"Put every setting of {LABELS.get(page, page)} back to its default value")
            reset.clicked.connect(lambda _=False, p=page: self._defaults(p))
            column.addWidget(reset, alignment=Qt.AlignLeft)
        column.addStretch()
        scroll = QScrollArea(widgetResizable=True, frameShape=QFrame.NoFrame)
        scroll.setWidget(body)
        return scroll

    @staticmethod
    def param_widget(p, v):
        """The Settings page's box for the Param `p` showing `v` (also Help > Run Setup Again's cores box): its limits, steps and suffix."""
        if isinstance(p.default, bool):
            w = QCheckBox()
            w.setChecked(v)
        elif p.choices:
            w = QComboBox()
            w.addItems(p.choices)
            w.setSizeAdjustPolicy(QComboBox.AdjustToMinimumContentsLengthWithIcon)   # not as wide as the longest choice (a long list)
            w.setMinimumContentsLength(32)
            w.setCurrentText(v)
        elif isinstance(p.default, str):
            w = _PathBox(p.browse) if p.browse else QLineEdit()
            w.setText(v)
        elif isinstance(p.default, int):
            w = QSpinBox(minimum=p.minimum if p.minimum is not None else 0,
                         maximum=p.maximum if p.maximum is not None else 1_000_000, suffix=p.suffix)
            w.setSingleStep(p.step or 1)
            w.setValue(v)
        else:
            w = QDoubleSpinBox(minimum=p.minimum if p.minimum is not None else 0.0,
                               maximum=p.maximum if p.maximum is not None else 1e9, suffix=p.suffix, decimals=2)
            w.setSingleStep(p.step or 0.1)
            w.setValue(v)
        return w

    @staticmethod
    def _get(w):
        if isinstance(w, QCheckBox):
            return w.isChecked()
        if isinstance(w, (QLineEdit, _PathBox)):
            return w.text().strip()
        return w.currentText() if isinstance(w, QComboBox) else w.value()

    @staticmethod
    def _put(w, v):
        if isinstance(w, QCheckBox):
            w.setChecked(v)
        elif isinstance(w, QComboBox):
            w.setCurrentText(v)
        elif isinstance(w, (QLineEdit, _PathBox)):
            w.setText(v)
        else:
            w.setValue(v)

    def _defaults(self, page):
        if page in self.custom:
            self.custom[page].defaults()
        for (name, _), (p, w) in self.widgets.items():
            if name in self.pages[page]:
                self._put(w, p.default)

    def accept(self):
        for w in self.custom.values():   # a custom page that refuses (a wrong file name...) keeps the dialog open
            if not w.save():
                return
        chosen = {}   # tool -> {key: value} as the boxes show them
        for (name, key), (_, w) in self.widgets.items():
            chosen.setdefault(name, {})[key] = self._get(w)
        notes = []
        for name, got in chosen.items():
            if name in CHECKS:
                try:
                    note = CHECKS[name](got)
                except UserError as e:
                    QMessageBox.warning(self, LABELS.get(name, name), str(e))
                    self.tabs.show_page(PAGE_OF.get(name, name))
                    return
                if note:
                    notes.append(note)
                    for key, v in got.items():   # the boxes show what is saved
                        self._put(self.widgets[name, key][1], v)
        changed = set()
        for name, got in chosen.items():
            before = dict(values(name))   # (set_value empties the tool's cache: read it once)
            for key, v in got.items():
                if v != before[key]:
                    set_value(name, key, v)
                    changed.add(name)
        if notes:
            QMessageBox.information(self, "Settings", "\n".join(notes))
        super().accept()
        for name in changed:
            notifier.changed.emit(name)


class _PathBox(QWidget):
    """A file name and a Browse… button (a Param with `browse`)."""

    def __init__(self, file_filter):
        super().__init__()
        self.line = QLineEdit()
        self.line.setMinimumWidth(self.line.fontMetrics().horizontalAdvance("M") * 24)   # room for a file name
        row = QHBoxLayout(self)
        row.setContentsMargins(0, 0, 0, 0)
        row.addWidget(self.line)
        row.addWidget(QPushButton("Browse…", clicked=lambda: self._browse(file_filter)))

    def _browse(self, file_filter):
        path, _ = QFileDialog.getOpenFileName(self, "Choose a file", self.line.text(), file_filter)
        if path:
            self.line.setText(path)

    def text(self):
        return self.line.text()

    def setText(self, v):
        self.line.setText(v)
