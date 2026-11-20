"""The one table of the application: a table of rows that stand for things the program holds (hits, guides, records, fragments), with
what every such table offers, written once:

  * click a column header to sort (numbers as numbers: "1,234-5,678", "94 %" and "12 aa" sort by their first number; a column with a
    better order gives its own key); the order is kept when the table is refilled; before any click the rows stay as they were filled;
  * right-click a header, or a cell, to filter by values (value_filter.py): only for the columns named in `filtered`, the categories;
  * the row a user sees is not the thing's number: `data_row(row)` is the number the row was filled as, and `row_of(number)` is the
    reverse. Code that reads its list of things by the row number (`self.hits[row]`) must read it through `data_row` once the table
    can be sorted. `current_data_row()` and `selected_data_rows()` are the two usual ways in.

To fill it:

    with table.refilling(len(things)):
        for i, thing in enumerate(things):
            table.put(i, 0, thing.name)                          # sorts by its text (or its first number)
            table.put(i, 1, f"{thing.n:,}", key=thing.n)         # a key of its own; right=True for a number's alignment

Never write a cell of a sorted table outside `refilling` unless the sort column cannot change (a cell written in the sort column moves
its row under the caller). Not for fixed lists that are not data (shortcuts, the self-test, timings, the codon grid): they keep a plain
QTableWidget. Every other table that lists things the program holds is a DataTable.
"""
import re
from contextlib import contextmanager

from PySide6.QtCore import Qt
from PySide6.QtWidgets import QTableWidget, QTableWidgetItem

import value_filter

INDEX = Qt.UserRole + 1   # in the first cell of a row: the number the row was filled as
LEADING_NUMBER = re.compile(r"\s*([-+]?\d[\d,]*\.?\d*)")


def sort_key(text):
    """What a cell sorts by when it is given no key: a number it begins with (numbers first, in order), else its text, without case."""
    m = LEADING_NUMBER.match(text)
    if m:
        try:
            return (0, float(m.group(1).replace(",", "")), text.lower())
        except ValueError:
            pass
    return (1, 0.0, text.lower())


class KeyedItem(QTableWidgetItem):
    """A cell that sorts by a key (a number, a tuple) kept in its UserRole data instead of by its text."""

    def __lt__(self, other):
        a, b = self.data(Qt.UserRole), other.data(Qt.UserRole)
        if a is None or b is None:
            return super().__lt__(other)
        try:
            return a < b
        except TypeError:   # keys of different kinds in one column (a number and a dash)
            return str(a) < str(b)


class DataTable(QTableWidget):
    """`columns`: the header texts. `filtered`: those of them that hold categories, to filter by value (right-click). `keys`: for
    value_filter, {column: function(text) -> the value it is filtered by}. `cell_menu`: a right-click on a cell of a filtered column offers
    Show only / Hide its value; False for a table with a context menu of its own, which adds `self.value_filter.cell_actions(menu, pos)`."""

    def __init__(self, columns, filtered=(), keys=None, cell_menu=True, parent=None):
        super().__init__(0, len(columns), parent)
        self.setHorizontalHeaderLabels(list(columns))
        header = self.horizontalHeader()
        header.setSortIndicator(-1, Qt.AscendingOrder)   # no column sorts until one is clicked (Qt's own default is the first, descending)
        self.setSortingEnabled(True)
        if filtered:
            value_filter.install(self, list(filtered), keys=keys, cells=cell_menu)

    @contextmanager
    def refilling(self, rows):
        """Empty the table and give it `rows` rows to write with put(); the sort the user chose is kept (applied when the block ends)."""
        header = self.horizontalHeader()
        section, order = header.sortIndicatorSection(), header.sortIndicatorOrder()
        self.setSortingEnabled(False)   # (a row written in a sorted table would be moved before its other cells are written)
        self.setRowCount(0)
        self.setRowCount(rows)
        try:
            yield self
        finally:
            header.setSortIndicator(section, order)
            self.setSortingEnabled(True)

    def put(self, row, column, text, key=None, right=False, tip=None, number=None):
        """Write a cell while refilling and return its item. `key`: what the column sorts by (default: see sort_key); `right`: a number.
        `number`: the thing the row stands for, when it is not the row's number at filling (rewriting one row of a sorted table)."""
        item = KeyedItem(text)
        item.setData(Qt.UserRole, sort_key(text) if key is None else key)
        if right:
            item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
        if tip:
            item.setToolTip(tip)
        if column == 0:
            item.setData(INDEX, row if number is None else number)
        self.setItem(row, column, item)
        return item

    # ---- the row a user sees, and the thing it stands for ----------------------------------------------------------------------
    def data_row(self, row):
        """The number row `row` of the table was filled as (the index of its thing), or -1 for no such row."""
        item = self.item(row, 0) if 0 <= row < self.rowCount() else None
        number = None if item is None else item.data(INDEX)
        return -1 if number is None else number

    def current_data_row(self):
        return self.data_row(self.currentRow())

    def selected_data_rows(self):
        """The things of the selected rows, by number."""
        return sorted({self.data_row(i.row()) for i in self.selectionModel().selectedRows()} - {-1})

    def row_of(self, number):
        """The row now showing thing `number` (sorting moves it), or None."""
        return next((r for r in range(self.rowCount()) if self.data_row(r) == number), None)

    def select_data_row(self, number):
        row = self.row_of(number)
        if row is not None:
            self.selectRow(row)
