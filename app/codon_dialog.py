"""The codon editor (Tools > Settings > General > Genetic code > "Codon usage tables…", or the toolbar switch's "Manage codon tables…"):
two tabs, the codon usage tables and the genetic codes, and the table choice on the Genetic code page.

Left, your tables (none ships with openDNA: New, or Import a CSV). Right, the 64 codons of the selected table with their frequency
per thousand codons, editable, and the same frequency relative to the best synonymous codon."""
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QCheckBox, QComboBox, QDialog, QDialogButtonBox, QTabWidget, QFileDialog, QHBoxLayout, QHeaderView, QInputDialog, QLabel, QLineEdit, QListWidget,
    QListWidgetItem, QMessageBox, QPushButton, QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget,
)

import codon_usage
import genetic_code
import tool_settings


class UsagePage(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.tables = codon_usage.all_tables()
        self.current = None   # the table shown
        self.edited = False   # its numbers were changed and not saved yet
        self.names = QListWidget()
        self.names.setFixedWidth(230)
        self.names.currentRowChanged.connect(lambda _row: self._show(self.names.currentItem().data(Qt.UserRole) if self.names.currentItem() else None))
        buttons = QVBoxLayout()
        for text, tip, slot in (("New", "A new table: all codons equal", self.new), ("Duplicate", "A copy you can edit", self.duplicate),
                                ("Import…", "Read a CSV (codon, amino acid, per thousand)", self.import_csv),
                                ("Export…", "Write this table as CSV", self.export_csv), ("Delete", "Delete this table of yours", self.delete)):
            buttons.addWidget(QPushButton(text, toolTip=tip, clicked=slot))
        self.buttons = buttons
        left = QVBoxLayout()
        left.addWidget(self.names, 1)
        left.addLayout(buttons)
        self.title = QLabel()
        self.source = QLineEdit(placeholderText="Source or note")
        self.source.textEdited.connect(lambda _t: setattr(self, "edited", True))
        self.grid = QTableWidget(64, 4)
        self.grid.setHorizontalHeaderLabels(["Codon", "Amino acid", "Per thousand", "Relative to the best synonym"])
        self.grid.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.grid.verticalHeader().hide()
        self.grid.itemChanged.connect(self._edited)
        self.total = QLabel()
        right = QVBoxLayout()
        for w in (self.title, self.source, self.grid, self.total):
            right.addWidget(w)
        use = QPushButton("Use this table", toolTip="Score coding sequences against it (Plasmid Checkup)", clicked=self.use)
        row = QHBoxLayout()
        row.addWidget(use)
        row.addStretch()
        right.addLayout(row)
        top = QHBoxLayout(self)
        top.addLayout(left)
        top.addLayout(right, 1)
        self.use_button = use
        self._fill(codon_usage.active_name())

    def _fill(self, select=None):
        self.names.blockSignals(True)
        self.names.clear()
        for name, t in self.tables.items():
            item = QListWidgetItem(name)
            item.setData(Qt.UserRole, name)
            self.names.addItem(item)
        self.names.blockSignals(False)
        for i in range(self.names.count()):
            if self.names.item(i).data(Qt.UserRole) == select:
                self.names.setCurrentRow(i)
                return
        self.names.setCurrentRow(0)
        if not self.names.count():   # no row, so no signal: say it is empty
            self._show(None)

    def _show(self, name):
        self._save()
        self.current = self.tables.get(name)
        t = self.current
        if t is None:
            self.title.setText("No codon usage table yet. <b>New</b> makes one; <b>Import…</b> reads a CSV (codon, amino acid, per thousand).")
            return
        editable = True
        self.title.setText(f"<b>{t.name}</b>")
        self.source.setText(t.source if not t.note else (t.source + " " + t.note).strip())
        self.source.setReadOnly(not editable)
        w = codon_usage.adaptiveness(t)
        code = genetic_code.codon_map(1)
        self.grid.blockSignals(True)
        for i, c in enumerate(codon_usage.CODONS):
            cells = [c, code[c], f"{t.per_thousand[c]:.2f}", f"{w[c]:.2f}" if c in w else ""]
            for j, text in enumerate(cells):
                item = QTableWidgetItem(text)
                item.setFlags(item.flags() | Qt.ItemIsEditable if (editable and j == 2) else item.flags() & ~Qt.ItemIsEditable)
                self.grid.setItem(i, j, item)
        self.grid.blockSignals(False)
        self.edited = False
        self._sum()
        for b in range(self.buttons.count()):
            if self.buttons.itemAt(b).widget().text() == "Delete":
                self.buttons.itemAt(b).widget().setEnabled(editable)

    def _sum(self):
        total = sum(float(self.grid.item(i, 2).text() or 0) for i in range(64))
        self.total.setText(f"Total {total:.1f} per thousand" + ("" if abs(total - 1000) < 5 else "  (a table usually adds up to 1000)"))

    def _edited(self, item):
        if item.column() == 2:
            try:
                float(item.text())
            except ValueError:
                item.setText(f"{self.current.per_thousand[codon_usage.CODONS[item.row()]]:.2f}")
                return
            self.edited = True
            self._sum()

    def _save(self):
        """Store the numbers edited on a table of yours (when leaving it)."""
        t = self.current
        if t is None or not self.edited:
            return
        for i, c in enumerate(codon_usage.CODONS):
            t.per_thousand[c] = float(self.grid.item(i, 2).text())
        t.source, t.note = self.source.text(), ""
        try:
            codon_usage.save_user(t)
        except ValueError as e:   # the file could not be written: said, not lost in silence
            QMessageBox.warning(self, "Codon table", str(e))

    # -- the buttons
    def _ask_name(self, title, start=""):
        name, ok = QInputDialog.getText(self, title, "Name of the table:", text=start)
        return name.strip() if ok and name.strip() else None

    def _add(self, table):
        try:
            codon_usage.save_user(table)
        except ValueError as e:
            QMessageBox.warning(self, "Codon table", str(e))
            return
        self.tables = codon_usage.all_tables()
        self._fill(table.name)

    def new(self):
        name = self._ask_name("New codon table")
        if name:
            self._add(codon_usage._table(name, {c: 1000 / 64 for c in codon_usage.CODONS}))

    def duplicate(self):
        if self.current is not None:
            name = self._ask_name("Duplicate codon table", f"{self.current.name} (copy)")
            if name:
                self._save()
                self._add(codon_usage.duplicate(self.current, name))

    def import_csv(self):
        path, _ = QFileDialog.getOpenFileName(self, "Import a codon table", "", "Codon tables (*.csv *.tsv *.txt);;All files (*)")
        if path:
            try:
                table = codon_usage.read_csv(path)
            except (OSError, ValueError) as e:
                QMessageBox.warning(self, "Codon table", str(e))
                return
            self._add(table)

    def export_csv(self):
        if self.current is not None:
            path, _ = QFileDialog.getSaveFileName(self, "Export the codon table", f"{self.current.name}.csv", "CSV (*.csv)")
            if path:
                codon_usage.write_csv(self.current, path)

    def delete(self):
        t = self.current
        if t is not None and QMessageBox.question(self, "Delete", f"Delete the table '{t.name}'?") == QMessageBox.Yes:
            self.current, self.edited = None, False
            codon_usage.delete_user(t.name)
            self.tables = codon_usage.all_tables()
            self._fill(codon_usage.active_name())

    def use(self):
        self._save()
        if self.current is not None:
            codon_usage.set_active(self.current.name)
            tool_settings.notifier.changed.emit(genetic_code.TOOL)

    def leave(self):
        self._save()


LETTERS = "ACDEFGHIKLMNPQRSTVWY*" + genetic_code.REMOVED   # the amino acids, the stop, and "X": the codon is not used any more


class CodesPage(QWidget):
    """The genetic codes: the NCBI tables and the recoded ones (read-only: Duplicate to change one) and yours. For each codon, the
    amino acid it codes for (X: the organism does not use it any more) and whether it can start a protein."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.current = None
        self.names = QListWidget()
        self.names.setFixedWidth(260)
        self.names.currentRowChanged.connect(lambda _row: self._show(self.names.currentItem().data(Qt.UserRole) if self.names.currentItem() else None))
        self.delete_button = QPushButton("Delete", toolTip="Delete this code of yours", clicked=self.delete)
        left = QVBoxLayout()
        left.addWidget(self.names, 1)
        for b in (QPushButton("Duplicate", toolTip="A copy you can edit", clicked=self.duplicate), self.delete_button):
            left.addWidget(b)
        self.title = QLabel()
        self.name = QLineEdit(placeholderText="Name")
        self.grid = QTableWidget(64, 3)
        self.grid.setHorizontalHeaderLabels(["Codon", "Amino acid (X: not used)", "Can start"])
        self.grid.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.grid.verticalHeader().hide()
        self.note = QLabel(wordWrap=True)
        self.save_button = QPushButton("Save", toolTip="Keep the changes to this code", clicked=self.save)
        right = QVBoxLayout()
        for w in (self.title, self.name, self.grid, self.note, self.save_button):
            right.addWidget(w)
        top = QHBoxLayout(self)
        top.addLayout(left)
        top.addLayout(right, 1)
        self.boxes = []
        self._fill(genetic_code.table_id())

    def _fill(self, select):
        self.names.blockSignals(True)
        self.names.clear()
        for id, t in sorted(genetic_code.TABLES.items()):
            item = QListWidgetItem(f"{id} – {t.name}" + ("  🔒" if id < genetic_code.FIRST_USER else ""))
            item.setData(Qt.UserRole, id)
            self.names.addItem(item)
        self.names.blockSignals(False)
        for i in range(self.names.count()):
            if self.names.item(i).data(Qt.UserRole) == select:
                self.names.setCurrentRow(i)
                return
        self.names.setCurrentRow(0)

    def _show(self, id):
        t = genetic_code.TABLES.get(id)
        self.current = t
        if t is None:
            return
        editable = id >= genetic_code.FIRST_USER
        self.title.setText(f"<b>{t.id} – {t.name}</b>" + ("" if editable else "  (built in: Duplicate to change it)"))
        self.name.setText(t.name)
        self.name.setReadOnly(not editable)
        self.boxes = []
        for i, c in enumerate(codon_usage.CODONS):
            self.grid.setItem(i, 0, QTableWidgetItem(c))
            self.grid.item(i, 0).setFlags(Qt.ItemIsEnabled)
            aa = QComboBox()
            aa.addItems(list(LETTERS))
            aa.setCurrentText(t.codons[c])
            start = QCheckBox(checked=t.starts[i] != "-")
            for w in (aa, start):
                w.setEnabled(editable)
            aa.currentTextChanged.connect(self._note)
            self.grid.setCellWidget(i, 1, aa)
            self.grid.setCellWidget(i, 2, start)
            self.boxes.append((aa, start))
        self.save_button.setEnabled(editable)
        self.delete_button.setEnabled(editable)
        self._note()

    def _note(self, *_):
        gone = [c for c, (aa, _s) in zip(codon_usage.CODONS, self.boxes) if aa.currentText() == genetic_code.REMOVED]
        self.note.setText(f"{64 - len(gone)} codons used" + (f"; not used: {', '.join(gone)}" if gone else "")
                          + (f"\n{self.current.source}" if self.current is not None and self.current.id >= genetic_code.FIRST_CUSTOM else ""))

    def _strings(self):
        aas = "".join(aa.currentText() for aa, _s in self.boxes)
        return aas, "".join("M" if start.isChecked() else "-" for _a, start in self.boxes)

    def save(self):
        try:
            id = genetic_code.save_code(self.name.text(), *self._strings(), id=self.current.id)
        except ValueError as e:
            QMessageBox.warning(self, "Genetic code", str(e))
            return
        self._fill(id)

    def duplicate(self):
        t = self.current
        if t is None:
            return
        name, ok = QInputDialog.getText(self, "Duplicate genetic code", "Name of the code:", text=f"{t.name} (copy)")
        if ok and name.strip():
            try:
                id = genetic_code.save_code(name, t.aas, t.starts, note=f"Copy of {t.name}.")
            except ValueError as e:
                QMessageBox.warning(self, "Genetic code", str(e))
                return
            self._fill(id)

    def delete(self):
        t = self.current
        if t is not None and t.id >= genetic_code.FIRST_USER and QMessageBox.question(self, "Delete", f"Delete the code '{t.name}'?") == QMessageBox.Yes:
            genetic_code.delete_code(t.id)
            self._fill(genetic_code.table_id())


class CodonDialog(QDialog):
    """The codon usage tables and the genetic codes, one tab each."""

    def __init__(self, parent=None, tab=0):
        super().__init__(parent)
        self.setWindowTitle("Codon tables")
        self.resize(860, 640)
        self.usage, self.codes = UsagePage(), CodesPage()
        self.tabs = QTabWidget()
        self.tabs.addTab(self.usage, "Codon usage")
        self.tabs.addTab(self.codes, "Genetic codes")
        self.tabs.setCurrentIndex(tab)
        close = QDialogButtonBox(QDialogButtonBox.Close)
        close.rejected.connect(self.reject)
        column = QVBoxLayout(self)
        column.addWidget(self.tabs)
        column.addWidget(close)

    def done(self, result):
        self.usage.leave()
        super().done(result)


class UsageBox(QWidget):
    """The "Codon usage" line of the Genetic code page of the Settings: the table in use, and the editor."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.choice = QComboBox()
        self.choice.setToolTip("The table Plasmid Checkup scores coding sequences against (codon adaptation index, rare codons)")
        self.choice.setSizeAdjustPolicy(QComboBox.AdjustToMinimumContentsLengthWithIcon)   # not as wide as its longest name
        self.choice.setMinimumContentsLength(28)
        column = QVBoxLayout(self)
        column.setContentsMargins(0, 0, 0, 0)
        row = QHBoxLayout()
        row.addWidget(QLabel("Codon usage table"))
        row.addWidget(self.choice, 1)
        column.addLayout(row)
        buttons = QHBoxLayout()
        buttons.addWidget(QPushButton("Codon usage tables…", toolTip="Make, edit, import and export codon usage tables", clicked=self.edit))
        buttons.addWidget(QPushButton("Genetic codes…", toolTip="The recoded and your own genetic codes", clicked=lambda: self.edit(1)))
        buttons.addStretch()
        column.addLayout(buttons)
        self.fill()

    def fill(self):
        self.choice.clear()
        self.choice.addItems(list(codon_usage.all_tables()))
        current = codon_usage.active()
        if current is not None:
            self.choice.setCurrentText(current.name)

    def edit(self, tab=0):
        CodonDialog(self, tab if isinstance(tab, int) else 0).exec()
        self.fill()

    def save(self):
        codon_usage.set_active(self.choice.currentText())
        return True

    def defaults(self):
        self.choice.setCurrentText(codon_usage.active_name())

