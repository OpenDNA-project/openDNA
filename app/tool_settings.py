"""Per-tool settings: each tool declares its parameters once (register), reads them with
value(tool, key), and Tools → Settings shows one page per view or tool, listed on the left in collapsible groups that follow
the new-tab picker (Views, Annotations, Tools) and a General group, each page with a form built from the declarations (a page
without options says so). Stored in settings.ini in the app data folder (tests use a temp folder).

A change emits notifier.changed(tool), so open tools can redo their work.
"""
from types import SimpleNamespace

from PySide6.QtCore import QObject, QSettings, Qt, Signal
from PySide6.QtGui import QPalette
from PySide6.QtWidgets import (
    QCheckBox, QComboBox, QFrame, QGroupBox, QLabel, QDialog, QDialogButtonBox, QDoubleSpinBox, QFormLayout, QHBoxLayout,
    QPushButton, QScrollArea, QSpinBox, QStackedWidget, QTreeWidget, QTreeWidgetItem, QVBoxLayout, QWidget,
)

from primer_db import data_dir

TOOLS = {}  # tool name -> [Param], in declaration order
# The name a tool's page shows when it differs from its key (the key stays: it is where users' saved values are stored)
LABELS = {"Detect Features": "Detect Common Features"}

GENERAL = "General"
GENERAL_PAGES = ["Genetic code", "Display", "History", "Updates"]   # settings that belong to no view; the Display page is the hover delay
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


def Param(key, label, default, minimum=None, maximum=None, step=None, suffix="", tip="", choices=None, info=None):
    """One setting. Its type follows `default`: bool -> checkbox, int / float -> spin box,
    str with `choices` -> drop-down list. `info`: optional function of the chosen value giving a small read-only line
    shown under the drop-down (e.g. the stop and start codons of the selected genetic code)."""
    return SimpleNamespace(key=key, label=label, default=default, minimum=minimum, maximum=maximum,
                           step=step, suffix=suffix, tip=tip, choices=choices, info=info)


class _Notifier(QObject):
    changed = Signal(str)  # tool name


notifier = _Notifier()


def store():
    return QSettings(str(data_dir() / "settings.ini"), QSettings.IniFormat)


def register(tool, params):
    TOOLS[tool] = list(params)


_cache = {}  # (tool, key) -> value: reading settings.ini takes ~ms, and drawing asks per base


def value(tool, key):
    if (tool, key) not in _cache:
        _cache[tool, key] = _read(tool, key)
    return _cache[tool, key]


def _read(tool, key):
    p = next(p for p in TOOLS[tool] if p.key == key)
    raw = store().value(f"{tool}/{key}", p.default)
    # QSettings gives strings back from an .ini file: convert to the default's type.
    if isinstance(p.default, bool):
        return raw if isinstance(raw, bool) else str(raw).lower() in ("true", "1")
    if p.choices:
        return raw if raw in p.choices else p.default
    try:
        return type(p.default)(raw)
    except (TypeError, ValueError):
        return p.default


def set_value(tool, key, v):
    store().setValue(f"{tool}/{key}", v)
    _cache.pop((tool, key), None)


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
                w = self._widget(p, value(name, p.key))
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
    def _widget(p, v):
        if isinstance(p.default, bool):
            w = QCheckBox()
            w.setChecked(v)
        elif p.choices:
            w = QComboBox()
            w.addItems(p.choices)
            w.setSizeAdjustPolicy(QComboBox.AdjustToMinimumContentsLengthWithIcon)   # not as wide as the longest choice (a long list)
            w.setMinimumContentsLength(32)
            w.setCurrentText(v)
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
        return w.currentText() if isinstance(w, QComboBox) else w.value()

    @staticmethod
    def _put(w, v):
        if isinstance(w, QCheckBox):
            w.setChecked(v)
        elif isinstance(w, QComboBox):
            w.setCurrentText(v)
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
        changed = set()
        for (name, key), (p, w) in self.widgets.items():
            v = self._get(w)
            if v != value(name, key):
                set_value(name, key, v)
                changed.add(name)
        super().accept()
        for name in changed:
            notifier.changed.emit(name)
