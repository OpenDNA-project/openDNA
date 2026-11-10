"""The codon editor (Tools > Settings > General > Genetic code > "Codon tables…", or the toolbar switch's "Manage codon tables…"):
one table for the genetic code and the codon usage together.

The 64 codons are laid out as the familiar chart — sixteen rows, one per first and second base, four columns for the third — so the
whole code is on screen at once instead of a list you scroll. Each cell carries the codon, the amino acid it codes for (X: the organism
does not use it any more), whether it can start a protein, and how often the organism uses it (per thousand codons).

The letters and the start flags belong to the genetic code, the frequencies to the codon usage table, and the two are still chosen
separately above the chart: one code is read with any organism's usage, and a table imported from a CSV or carried in a project stays
usable with every code. Every table can be edited, a built-in one into a copy of itself: changing a letter or a start codon on an NCBI
table or a recoded one renames it "… (edited)" and Save writes a code of your own, so the tables themselves never change and a file
translates today to what it translated to yesterday. Only the usual NCBI tables are listed, with the recoded codes and yours; "All
tables" brings back the other twenty (a CDS with /transl_table=N is read with table N whatever the list shows).

Import and Export carry both at once: one CSV holds the frequencies and the code (its name in `# code:`, its letters in the `aa` column,
its start codons in `start`). Importing one adds the table and, unless a code already here says the very same thing, the code too.
"""
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QCheckBox, QComboBox, QDialog, QDialogButtonBox, QFileDialog, QHBoxLayout, QHeaderView, QInputDialog, QLabel, QLineEdit,
    QMessageBox, QPushButton, QTableWidget, QVBoxLayout, QWidget,
)

import codon_usage
import genetic_code
import tool_settings

LETTERS = genetic_code.LETTERS   # the amino acids, the stop, and "X": the codon is not used any more
# The NCBI tables the chart offers: the standard one, the two mitochondrial ones most asked for, and the bacterial table the recoded
# codes are built from. The other twenty are a tick away ("Show all NCBI tables"), never gone: a file with /transl_table=5 is still
# read with table 5 whatever this list shows.
SHORT_CODES = (1, 2, 3, 11)
PER_ROW = 4            # the chart: four columns, one per third base
ROWS = 64 // PER_ROW


def row_label(row):
    """The first two bases of the codons of `row` ("TT", "TC"...)."""
    return codon_usage.CODONS[row * PER_ROW][:2]


class CodonPage(QWidget):
    """The chart: the genetic code and the codon usage of one organism, side by side on the same 64 codons."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.tables = codon_usage.all_tables()
        self.current = None        # the usage table shown
        self.current_code = None   # the genetic code shown
        self.edited = False        # its frequencies were changed and not saved yet
        self.dirty = False         # its letters or start codons were changed and not saved yet
        self.boxes = []            # per codon, in codon_usage.CODONS order: (amino acid, can start, per thousand)

        self.code = QComboBox(toolTip="The genetic code: what each codon means")
        self.code.setSizeAdjustPolicy(QComboBox.AdjustToMinimumContentsLengthWithIcon)   # not as wide as the longest NCBI table name
        self.code.setMinimumContentsLength(24)
        self.code.currentIndexChanged.connect(lambda _i: self._show_code())
        self.show_all = QCheckBox("All tables", toolTip="List every NCBI translation table, not just the ones most used")
        self.show_all.toggled.connect(lambda _on: self._fill_codes(self.code.currentData()))
        self.name = QLineEdit(placeholderText="Name of the code")
        self.save_button = QPushButton("Save", toolTip="Keep the changes to this code", clicked=self.save)
        self.duplicate_code = QPushButton("Duplicate", toolTip="A copy you can edit", clicked=self.duplicate)
        self.delete_code = QPushButton("Delete", toolTip="Delete this code of yours", clicked=self.delete)

        self.usage = QComboBox(toolTip="The codon usage table: how often the organism uses each codon")
        self.usage.setSizeAdjustPolicy(QComboBox.AdjustToMinimumContentsLengthWithIcon)
        self.usage.setMinimumContentsLength(16)
        self.usage.currentIndexChanged.connect(lambda _i: self._show_usage())
        self.source = QLineEdit(placeholderText="Source or note")
        self.source.textEdited.connect(lambda _t: setattr(self, "edited", True))
        self.usage_buttons = {}
        for text, tip, slot in (("New", "A new table: all codons equal", self.new), ("Duplicate", "A copy you can edit", self.duplicate_table),
                                ("Import…", "Read a CSV: its frequencies and the genetic code in it", self.import_csv),
                                ("Export…", "Write this table and the code shown as one CSV", self.export_csv),
                                ("Delete", "Delete this table of yours", self.delete_table),
                                ("Use", "Score coding sequences against it (Plasmid Checkup)", self.use)):
            self.usage_buttons[text] = QPushButton(text, toolTip=tip, clicked=slot)

        self.grid = QTableWidget(ROWS, PER_ROW)
        self.grid.setHorizontalHeaderLabels(list(codon_usage.BASES))
        self.grid.setVerticalHeaderLabels([row_label(r) for r in range(ROWS)])
        self.grid.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.grid.setSelectionMode(QTableWidget.NoSelection)
        for i, c in enumerate(codon_usage.CODONS):
            self.grid.setCellWidget(i // PER_ROW, i % PER_ROW, self._cell(i, c))
        self.grid.resizeRowsToContents()
        self.title = QLabel()
        self.note = QLabel(wordWrap=True)
        self.total = QLabel()

        column = QVBoxLayout(self)
        column.addWidget(self.title)
        code_row = QHBoxLayout()
        code_row.addWidget(QLabel("Genetic code"))
        code_row.addWidget(self.code, 1)
        code_row.addWidget(self.show_all)
        code_row.addWidget(self.name, 1)
        for b in (self.save_button, self.duplicate_code, self.delete_code):
            code_row.addWidget(b)
        column.addLayout(code_row)
        usage_row = QHBoxLayout()
        usage_row.addWidget(QLabel("Codon usage"))
        usage_row.addWidget(self.usage, 1)
        usage_row.addWidget(self.source, 1)
        for b in self.usage_buttons.values():
            usage_row.addWidget(b)
        column.addLayout(usage_row)
        column.addWidget(self.grid, 1)
        column.addWidget(self.note)
        column.addWidget(self.total)

        self._fill_codes()
        self._fill_usage()

    def _cell(self, i, codon):
        """One cell of the chart: the codon, its amino acid, whether it can start, and its frequency."""
        aa = QComboBox(toolTip=f"The amino acid {codon} codes for ({genetic_code.REMOVED}: not used any more)")
        aa.addItems(list(LETTERS))
        aa.currentTextChanged.connect(self._code_edited)
        start = QCheckBox(toolTip=f"{codon} can start a protein")
        start.toggled.connect(self._code_edited)
        freq = QLineEdit(alignment=Qt.AlignRight, toolTip=f"How often the organism uses {codon}, per thousand codons")
        freq.setFixedWidth(54)
        freq.textChanged.connect(lambda text, k=i: self._usage_edited(k, text))
        self.boxes.append((aa, start, freq))
        box = QWidget()
        row = QHBoxLayout(box)
        row.setContentsMargins(4, 1, 4, 1)
        row.setSpacing(3)
        row.addWidget(QLabel(f"<b>{codon}</b>"))
        row.addWidget(aa)
        row.addWidget(start)
        row.addWidget(freq)
        return box

    # ---- the genetic code -------------------------------------------------------------------------------------------------
    def _listed(self, select):
        """The code ids the dropdown offers: the usual NCBI tables, the recoded ones, yours — and whatever is selected, so a code
        chosen elsewhere is never silently swapped for another one."""
        return sorted(t.id for t in genetic_code.TABLES.values()
                      if self.show_all.isChecked() or t.id in SHORT_CODES or t.id >= genetic_code.FIRST_CUSTOM or t.id == select)

    def _fill_codes(self, select=None):
        want = genetic_code.table_id() if select is None else select
        self.code.blockSignals(True)
        self.code.clear()
        previous = None
        for id in self._listed(want):
            if previous is not None and previous < genetic_code.FIRST_CUSTOM <= id:
                self.code.insertSeparator(self.code.count())   # the NCBI tables, then the recoded ones and yours
            self.code.addItem(f"{id} – {genetic_code.TABLES[id].name}", id)
            previous = id
        self.code.blockSignals(False)
        self.code.setCurrentIndex(max(0, self.code.findData(want)))
        self._show_code()

    def _show_code(self):
        # shortcut: changing the dropdown drops an unsaved letter or start change without asking, as it did before any table could be
        # edited; ask to keep it (or save it like the frequencies are) if anyone loses work this way.
        t = genetic_code.TABLES.get(self.code.currentData())
        self.current_code = t
        self.dirty = False
        if t is None:
            return
        self.name.setText(t.name)
        for i, c in enumerate(codon_usage.CODONS):
            aa, start, _freq = self.boxes[i]
            for w in (aa, start):
                w.blockSignals(True)
            aa.setCurrentText(t.codons[c])
            start.setChecked(t.starts[i] != "-")
            for w in (aa, start):
                w.blockSignals(False)
                w.setEnabled(True)      # every table can be edited; a built-in one into a copy of itself (see _code_edited)
        self.delete_code.setEnabled(t.id >= genetic_code.FIRST_USER)
        self._update_save()
        self._titles()
        self._note()

    def _code_edited(self, *_):
        """A letter or a start codon was changed. On a built-in table the edit starts a copy — NCBI's tables and the recoded ones are
        never altered, so what a file translated to yesterday it still translates to today."""
        t = self.current_code
        if t is not None and not self.dirty:
            self.dirty = True
            if t.id < genetic_code.FIRST_USER:
                self.name.setText(f"{t.name} (edited)")
            self._update_save()
            self._titles()
        self._note()

    def _update_save(self):
        own = self.current_code is not None and self.current_code.id >= genetic_code.FIRST_USER
        self.save_button.setEnabled(own or self.dirty)
        self.save_button.setToolTip("Keep the changes to this code" if own else "Save these changes as a genetic code of your own")

    def _strings(self):
        aas = "".join(aa.currentText() for aa, _s, _f in self.boxes)
        return aas, "".join("M" if start.isChecked() else "-" for _a, start, _f in self.boxes)

    def save(self):
        """Keep the changes: on a code of yours in place, on a built-in table as a new code of your own."""
        t = self.current_code
        if t is None:
            return
        own = t.id >= genetic_code.FIRST_USER
        try:
            id = genetic_code.save_code(self.name.text(), *self._strings(), id=t.id if own else None,
                                        note="" if own else f"Edited copy of {t.id} – {t.name}.")
        except ValueError as e:
            QMessageBox.warning(self, "Genetic code", str(e))
            return
        self._fill_codes(id)

    def duplicate(self):
        t = self.current_code
        if t is None:
            return
        name, ok = QInputDialog.getText(self, "Duplicate genetic code", "Name of the code:", text=f"{t.name} (copy)")
        if ok and name.strip():
            try:
                id = genetic_code.save_code(name, t.aas, t.starts, note=f"Copy of {t.name}.")
            except ValueError as e:
                QMessageBox.warning(self, "Genetic code", str(e))
                return
            self._fill_codes(id)

    def delete(self):
        t = self.current_code
        if t is not None and t.id >= genetic_code.FIRST_USER and QMessageBox.question(self, "Delete", f"Delete the code '{t.name}'?") == QMessageBox.Yes:
            genetic_code.delete_code(t.id)
            self._fill_codes(genetic_code.table_id())

    def _note(self, *_):
        gone = [c for c, (aa, _s, _f) in zip(codon_usage.CODONS, self.boxes) if aa.currentText() == genetic_code.REMOVED]
        t = self.current_code
        self.note.setText(f"{64 - len(gone)} codons used" + (f"; not used: {', '.join(gone)}" if gone else "")
                          + (f"\n{t.source}" if t is not None and t.id >= genetic_code.FIRST_CUSTOM else ""))

    # ---- the codon usage table --------------------------------------------------------------------------------------------
    def _fill_usage(self, select=None):
        self.usage.blockSignals(True)
        self.usage.clear()
        for name in self.tables:
            self.usage.addItem(name, name)
        self.usage.blockSignals(False)
        found = self.usage.findData(codon_usage.active_name() if select is None else select)
        self.usage.setCurrentIndex(max(0, found))
        self._show_usage()

    def _show_usage(self):
        self._save()
        self.current = self.tables.get(self.usage.currentData())
        t = self.current
        for _aa, _start, freq in self.boxes:
            freq.blockSignals(True)
            freq.setText("")
            freq.setEnabled(t is not None)
            freq.blockSignals(False)
        if t is not None:
            for i, c in enumerate(codon_usage.CODONS):
                freq = self.boxes[i][2]
                freq.blockSignals(True)
                freq.setText(f"{t.per_thousand[c]:.2f}")
                freq.blockSignals(False)
            self.source.setText(t.source if not t.note else (t.source + " " + t.note).strip())
        else:
            self.source.setText("")
        self.source.setEnabled(t is not None)
        for text in ("Duplicate", "Export…", "Delete", "Use"):
            self.usage_buttons[text].setEnabled(t is not None)
        self.edited = False
        self._titles()
        self._sum()

    def _titles(self):
        t, code = self.current, self.current_code
        if code is None:
            code_text = "no genetic code"
        else:
            code_text = f"<b>{code.id} – {code.name}</b>"
            if code.id >= genetic_code.FIRST_USER:
                pass                      # a code of yours: Save writes it back
            elif self.dirty:
                code_text += " · built in, so this is a copy: <b>Save</b> makes it yours"
            else:
                code_text += " · built in: editing it starts a copy of your own"
        if t is None:
            self.title.setText(code_text + " · No codon usage table yet. <b>New</b> makes one; <b>Import…</b> reads a CSV "
                                           "(codon, amino acid, per thousand).")
        else:
            self.title.setText(f"{code_text} · <b>{t.name}</b>")

    def _value(self, i):
        try:
            return float(self.boxes[i][2].text())
        except ValueError:
            return 0.0

    def _sum(self):
        if self.current is None:
            self.total.setText("")
            return
        total = sum(self._value(i) for i in range(64))
        self.total.setText(f"Total {total:.1f} per thousand" + ("" if abs(total - 1000) < 5 else "  (a table usually adds up to 1000)"))

    def _usage_edited(self, i, text):
        """A frequency typed in: a number is kept, anything else brings the old one back."""
        if self.current is None:
            return
        try:
            float(text)
        except ValueError:
            freq = self.boxes[i][2]
            freq.blockSignals(True)
            freq.setText(f"{self.current.per_thousand[codon_usage.CODONS[i]]:.2f}")
            freq.blockSignals(False)
            return
        self.edited = True
        self._sum()

    def _save(self):
        """Store the frequencies edited on a table of yours (when leaving it)."""
        t = self.current
        if t is None or not self.edited:
            return
        for i, c in enumerate(codon_usage.CODONS):
            t.per_thousand[c] = self._value(i)
        t.source, t.note = self.source.text(), ""
        try:
            codon_usage.save_user(t)
        except ValueError as e:   # the file could not be written: said, not lost in silence
            QMessageBox.warning(self, "Codon table", str(e))

    def _ask_name(self, title, start=""):
        name, ok = QInputDialog.getText(self, title, "Name of the table:", text=start)
        return name.strip() if ok and name.strip() else None

    def _add(self, table):
        try:
            codon_usage.save_user(table)
        except ValueError as e:
            QMessageBox.warning(self, "Codon table", str(e))
            return False
        self.tables = codon_usage.all_tables()
        self._fill_usage(table.name)
        return True

    def new(self):
        name = self._ask_name("New codon table")
        if name:
            self._add(codon_usage._table(name, {c: 1000 / 64 for c in codon_usage.CODONS}))

    def duplicate_table(self):
        if self.current is not None:
            name = self._ask_name("Duplicate codon table", f"{self.current.name} (copy)")
            if name:
                self._save()
                self._add(codon_usage.duplicate(self.current, name))

    def import_csv(self):
        """Read a CSV: its frequencies become a table of yours, and the genetic code it carries (if any) a code of yours."""
        path, _ = QFileDialog.getOpenFileName(self, "Import a codon table", "", "Codon tables (*.csv *.tsv *.txt);;All files (*)")
        if path:
            try:
                table = codon_usage.read_csv(path)
            except (OSError, ValueError) as e:
                QMessageBox.warning(self, "Codon table", str(e))
                return
            if self._add(table):   # a table that could not be saved: its code is not imported on its own
                self._adopt_code(table)

    def _adopt_code(self, table):
        """The genetic code a file carried: the code of yours that already says the same thing, else a new one."""
        if not table.aas:
            return
        same = next((t for t in genetic_code.TABLES.values() if t.aas == table.aas and t.starts == table.starts), None)
        if same is not None:        # the standard table, or a code already here: nothing to add
            self._fill_codes(same.id)
            return
        wanted = table.code_name or f"{table.name} code"
        note = f"Imported with the codon usage table {table.name}."
        for name in (wanted, f"{wanted} (imported)"):   # the second try when a code of that name says something else
            try:
                self._fill_codes(genetic_code.save_code(name, table.aas, table.starts, note=note))
                return
            except ValueError as e:
                last = e
        QMessageBox.warning(self, "Genetic code", str(last))

    def export_csv(self):
        """Write the table and the code shown to one CSV, so the pair travels together."""
        if self.current is not None:
            path, _ = QFileDialog.getSaveFileName(self, "Export the codon table", f"{self.current.name}.csv", "CSV (*.csv)")
            if path:
                self._save()
                codon_usage.write_csv(self.current, path, self.current_code)

    def delete_table(self):
        t = self.current
        if t is not None and QMessageBox.question(self, "Delete", f"Delete the table '{t.name}'?") == QMessageBox.Yes:
            self.current, self.edited = None, False
            codon_usage.delete_user(t.name)
            self.tables = codon_usage.all_tables()
            self._fill_usage(codon_usage.active_name())

    def use(self):
        self._save()
        if self.current is not None:
            codon_usage.set_active(self.current.name)
            tool_settings.notifier.changed.emit(genetic_code.TOOL)

    def leave(self):
        self._save()


class CodonDialog(QDialog):
    """The one chart: the genetic code and the codon usage of one organism together."""

    def __init__(self, parent=None, tab=0):
        super().__init__(parent)
        self.setWindowTitle("Codon tables")
        self.resize(1060, 700)   # the chart's own minimum is ~1030 wide: opening at it would leave nothing to spare
        self.page = CodonPage()
        close = QDialogButtonBox(QDialogButtonBox.Close)
        close.rejected.connect(self.reject)
        column = QVBoxLayout(self)
        column.addWidget(self.page)
        column.addWidget(close)

    def done(self, result):
        self.page.leave()
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
        buttons.addWidget(QPushButton("Codon tables…", toolTip="The genetic codes and the codon usage tables, on one chart",
                                      clicked=self.edit))
        buttons.addStretch()
        column.addLayout(buttons)
        self.fill()

    def fill(self):
        self.choice.clear()
        self.choice.addItems(list(codon_usage.all_tables()))
        current = codon_usage.active()
        if current is not None:
            self.choice.setCurrentText(current.name)

    def edit(self, *_):
        CodonDialog(self).exec()
        self.fill()

    def save(self):
        codon_usage.set_active(self.choice.currentText())
        return True

    def defaults(self):
        self.choice.setCurrentText(codon_usage.active_name())
