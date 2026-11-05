"""Per-tool settings: each tool declares its parameters once (register), reads them with
value(tool, key), and Tools → Settings shows one page per tool (listed on the left) with a form built from the
declarations. Stored in settings.ini in the app data folder (tests use a temp folder).

A change emits notifier.changed(tool), so open tools can redo their work.
"""
from types import SimpleNamespace

from PySide6.QtCore import QObject, QSettings, Qt, Signal
from PySide6.QtGui import QPalette
from PySide6.QtWidgets import (
    QCheckBox, QComboBox, QLabel, QDialog, QDialogButtonBox, QDoubleSpinBox, QFormLayout, QHBoxLayout, QListWidget, QPushButton,
    QSpinBox, QStackedWidget, QVBoxLayout, QWidget,
)

from primer_db import data_dir

TOOLS = {}  # tool name -> [Param], in declaration order
# The name a tool's page shows when it differs from its key (the key stays: it is where users' saved values are stored)
LABELS = {"Detect Features": "Detect Common Features"}


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
    """The pages of the settings: a list of their names on the left, the chosen page on the right (the usual Preferences layout;
    nine tabs did not fit across the dialog and scrolled, the last ones cut off). Answers like a QTabWidget (addTab, widget,
    currentIndex, setCurrentIndex, count, tabText), so the callers need not know."""

    def __init__(self):
        super().__init__()
        self.names = QListWidget(spacing=2)
        self.stack = QStackedWidget()
        self.names.currentRowChanged.connect(self.stack.setCurrentIndex)
        row = QHBoxLayout(self)
        row.setContentsMargins(0, 0, 0, 0)
        row.addWidget(self.names)
        row.addWidget(self.stack, 1)

    def addTab(self, page, name):
        self.names.addItem(name)
        self.stack.addWidget(page)
        # the list as wide as its longest name (plus room for a scroll bar), no wider
        self.names.setFixedWidth(max(self.names.sizeHintForColumn(0), 100) + 28)
        if self.names.currentRow() < 0:
            self.names.setCurrentRow(0)
        return self.stack.count() - 1

    def widget(self, i):
        return self.stack.widget(i)

    def count(self):
        return self.stack.count()

    def tabText(self, i):
        return self.names.item(i).text()

    def currentIndex(self):
        return self.stack.currentIndex()

    def currentWidget(self):
        return self.stack.currentWidget()

    def setCurrentIndex(self, i):
        self.names.setCurrentRow(i)


class SettingsDialog(QDialog):
    """One page per tool (or `tool`'s page first). OK saves and notifies."""

    def __init__(self, parent, tool=None):
        super().__init__(parent)
        self.setWindowTitle("Settings")
        self.tabs = Pages()
        self.widgets = {}  # (tool, key) -> widget
        for name, params in TOOLS.items():
            page = QWidget()
            form = QFormLayout(page)
            form.setFieldGrowthPolicy(QFormLayout.FieldsStayAtSizeHint)   # a number box as wide as a number, not the whole page
            for p in params:
                w = self._widget(p, value(name, p.key))
                w.setToolTip(p.tip)
                self.widgets[(name, p.key)] = (p, w)
                form.addRow(p.label, w)
                if p.info and isinstance(w, QComboBox):   # a small read-only line under the choice, kept in step with it
                    line = QLabel(p.info(w.currentText()), wordWrap=True)
                    line.setForegroundRole(QPalette.PlaceholderText)   # faint, in either theme
                    w.currentTextChanged.connect(lambda text, line=line, info=p.info: line.setText(info(text)))
                    form.addRow("", line)
            reset = QPushButton("Restore Defaults", toolTip=f"Put every setting of {name} back to its default value")
            reset.clicked.connect(lambda _=False, n=name: self._defaults(n))
            form.addRow("", reset)
            form.setAlignment(reset, Qt.AlignLeft)
            i = self.tabs.addTab(page, LABELS.get(name, name))
            if name == tool:
                self.tabs.setCurrentIndex(i)
        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout = QVBoxLayout(self)
        layout.addWidget(self.tabs)
        layout.addWidget(buttons)

    @staticmethod
    def _widget(p, v):
        if isinstance(p.default, bool):
            w = QCheckBox()
            w.setChecked(v)
        elif p.choices:
            w = QComboBox()
            w.addItems(p.choices)
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

    def _defaults(self, tool):
        for (name, _), (p, w) in self.widgets.items():
            if name == tool:
                self._put(w, p.default)

    def accept(self):
        changed = set()
        for (name, key), (p, w) in self.widgets.items():
            v = self._get(w)
            if v != value(name, key):
                set_value(name, key, v)
                changed.add(name)
        super().accept()
        for name in changed:
            notifier.changed.emit(name)
