"""Filter by values: right-click a column header (or click the ▼ it shows once filtered) to tick which of that column's values
are shown, like a spreadsheet's filter. Works on any QTableView / QTableWidget / QTreeWidget, by hiding rows, so the table's own
code (which reads its rows by index) does not change. A table is refilled often (an edit, a new search): the filter is applied
again after every refill, and a value it has never seen (a new row) is shown, not hidden.

    install(table, ["Type", "Strand"])   # the columns, by header text, that can be filtered (categories, not numbers)

The tables it is installed on, and why the others are left out, are in docs/FILTER_AUDIT.md.
"""
from PySide6.QtCore import QObject, Qt, QTimer, Signal
from PySide6.QtWidgets import (
    QAbstractItemView, QHBoxLayout, QLineEdit, QListWidget, QListWidgetItem, QMenu, QPushButton, QTreeWidget, QTreeWidgetItem,
    QVBoxLayout, QWidget, QWidgetAction,
)

MARK = " ▼"   # after a filtered column's header text
EMPTY = "(empty)"


class ValueFilter(QObject):
    """The filter of one table: {column: set of the values hidden}. Stored as what is hidden, so values that appear later show."""
    changed = Signal()

    def __init__(self, view, columns, keys=None):
        super().__init__(view)
        self.view, self.names, self.hidden = view, list(columns), {}
        self.keys = dict(keys or {})   # column name -> function(text) -> the value it is filtered by (a part of the text)
        self.tree = isinstance(view, QTreeWidget)
        self.header = view.header() if self.tree else view.horizontalHeader()
        self.header.setContextMenuPolicy(Qt.CustomContextMenu)
        self.header.customContextMenuRequested.connect(lambda pos: self.popup(self.header.logicalIndexAt(pos), pos))
        self._later = QTimer(self, singleShot=True, interval=0, timeout=self.apply)
        model = view.model()
        for signal in (model.modelReset, model.layoutChanged, model.rowsInserted, model.dataChanged):
            signal.connect(self._refilled)
        self._applying = False

    # ---- what the table holds ----------------------------------------------------------------------------------------------
    def columns(self):
        """{column number: header text without the mark} of the columns that can be filtered."""
        out = {}
        model = self.view.model()
        for c in range(model.columnCount()):
            text = str(model.headerData(c, Qt.Horizontal) or "").removesuffix(MARK)
            if text in self.names:
                out[c] = text
        return out

    def rows(self):
        """Top-level rows: (row number, {column: value}), the value being the text, or the part of it the column is filtered by."""
        model, cols = self.view.model(), self.columns()
        for r in range(model.rowCount()):
            yield r, {c: self._value(name, _text(model.index(r, c))) for c, name in cols.items()}

    def _value(self, name, text):
        key = self.keys.get(name)
        return (key(text) or EMPTY) if key is not None and text != EMPTY else text

    def values(self, column):
        """[(value, rows having it)] of a column, the most common first, then by name."""
        counts = {}
        for _, texts in self.rows():
            v = texts[column]
            counts[v] = counts.get(v, 0) + 1
        return sorted(counts.items(), key=lambda kv: (kv[0] == EMPTY, kv[0].lower()))

    # ---- applying ----------------------------------------------------------------------------------------------------------
    def _refilled(self, *_):
        if self.hidden and not self._applying:
            self._later.start()

    def apply(self):
        """Hide the rows with a hidden value in any filtered column, show the others; mark the filtered headers."""
        self._applying = True
        try:
            cols = self.columns()
            hidden = {c: self.hidden.get(name, set()) for c, name in cols.items()}
            root = self.view.invisibleRootItem() if self.tree else None
            for r, texts in self.rows():
                hide = any(texts[c] in values for c, values in hidden.items())
                if self.tree:
                    root.child(r).setHidden(hide)
                else:
                    self.view.setRowHidden(r, hide)
            self._mark(cols)
        finally:
            self._applying = False
        self.changed.emit()

    def _mark(self, cols):
        for c, name in cols.items():
            text = name + (MARK if self.hidden.get(name) else "")
            if self.tree:
                self.view.headerItem().setText(c, text)
            elif hasattr(self.view, "horizontalHeaderItem") and self.view.horizontalHeaderItem(c) is not None:
                self.view.horizontalHeaderItem(c).setText(text)
            elif hasattr(self.view.model(), "marked"):   # a model of its own (the feature table) adds the mark itself
                self.view.model().marked = {n for n, v in self.hidden.items() if v}
                self.view.model().headerDataChanged.emit(Qt.Horizontal, 0, self.view.model().columnCount() - 1)

    def set_hidden(self, name, values):
        """Hide the rows whose column `name` holds one of `values` (an empty set clears that column's filter)."""
        if values:
            self.hidden[name] = set(values)
        else:
            self.hidden.pop(name, None)
        self.apply()

    def clear(self):
        self.hidden.clear()
        self.apply()

    def shown_rows(self):
        """The row numbers not hidden (for the tests)."""
        if self.tree:
            root = self.view.invisibleRootItem()
            return [r for r in range(root.childCount()) if not root.child(r).isHidden()]
        return [r for r in range(self.view.model().rowCount()) if not self.view.isRowHidden(r)]

    # ---- the popup ---------------------------------------------------------------------------------------------------------
    def popup(self, column, pos):
        name = self.columns().get(column)
        if name is None:
            return
        menu = QMenu(self.view)
        picker = ValuePicker(self, name, column)
        action = QWidgetAction(menu)
        action.setDefaultWidget(picker)
        menu.addAction(action)
        picker.done.connect(menu.close)
        picker.search.setFocus()
        menu.exec(self.header.mapToGlobal(pos))


class ValuePicker(QWidget):
    """The values of one column, each with a tick (shown) and how many rows have it; a box narrows the list."""
    done = Signal()

    def __init__(self, vf, name, column):
        super().__init__()
        self.vf, self.name = vf, name
        self.search = QLineEdit(placeholderText=f"Find a value of {name}…", clearButtonEnabled=True)
        self.list = QListWidget()
        self.list.setSelectionMode(QAbstractItemView.NoSelection)
        hidden = vf.hidden.get(name, set())
        for value, n in vf.values(column):
            item = QListWidgetItem(f"{value}  ({n})")
            item.setData(Qt.UserRole, value)
            item.setFlags(item.flags() | Qt.ItemIsUserCheckable)
            item.setCheckState(Qt.Unchecked if value in hidden else Qt.Checked)
            self.list.addItem(item)
        self.list.setMinimumWidth(260)
        self.list.setMaximumHeight(320)
        self.search.textChanged.connect(self._narrow)
        row = QHBoxLayout()
        for text, slot in (("All", lambda: self._tick(True)), ("None", lambda: self._tick(False)),
                           ("Clear filter", self._clear), ("OK", self._ok)):
            row.addWidget(QPushButton(text, clicked=slot))
        layout = QVBoxLayout(self)
        layout.setContentsMargins(6, 6, 6, 6)
        layout.addWidget(self.search)
        layout.addWidget(self.list)
        layout.addLayout(row)

    def items(self):
        return [self.list.item(i) for i in range(self.list.count())]

    def _narrow(self, text):
        needle = text.strip().lower()
        for item in self.items():
            item.setHidden(bool(needle) and needle not in item.data(Qt.UserRole).lower())

    def _tick(self, on):
        """All / None act on the values the search box leaves listed."""
        for item in self.items():
            if not item.isHidden():
                item.setCheckState(Qt.Checked if on else Qt.Unchecked)

    def _clear(self):
        self.vf.set_hidden(self.name, set())
        self.done.emit()

    def _ok(self):
        self.vf.set_hidden(self.name, {i.data(Qt.UserRole) for i in self.items() if i.checkState() != Qt.Checked})
        self.done.emit()


def _text(index):
    text = index.data(Qt.DisplayRole)
    text = "" if text is None else str(text).strip()
    return text or EMPTY


def install(view, columns, keys=None):
    """Filter by values on these columns (header texts) of `view`; the ValueFilter is also `view.value_filter`. `keys`: {column:
    function(text) -> value} for a column filtered by a part of its text."""
    vf = view.value_filter = ValueFilter(view, columns, keys)
    vf.header.setToolTip(f"Right-click a column header ({', '.join(columns)}) to show only some of its values")
    return vf


class SortItem(QTreeWidgetItem):
    """A tree row that sorts by number where it shows one (`keys`: {column: value}, e.g. a length or a Tm), by text elsewhere."""
    keys = {}

    def __lt__(self, other):
        column = self.treeWidget().sortColumn() if self.treeWidget() is not None else 0
        mine, theirs = self.keys.get(column), getattr(other, "keys", {}).get(column)
        if mine is not None and theirs is not None and type(mine) is type(theirs):
            return mine < theirs
        return self.text(column).lower() < other.text(column).lower()
