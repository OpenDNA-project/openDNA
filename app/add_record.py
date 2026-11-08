"""Import > DNA > Sequence File as New Chromosome: the sequences of another file are added to the open file as new records (a
GenBank file then holds several records, like a genome and its plasmids; a FASTA file gets more sequences).

The new records are written into the file at once (the records already in it are copied unchanged, see multigb.append_genbank /
append_fasta), because every tool that works across the records of a file (the map's list, the whole-file search, the primer
sites on the other records, Detect Features) reads them from the file. The open record keeps its unsaved edits until Save.
"""
from pathlib import Path
from types import SimpleNamespace

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (QAbstractItemView, QDialog, QDialogButtonBox, QHeaderView, QLabel, QTableWidget, QTableWidgetItem,
                               QVBoxLayout)

from opendna_core import Document, Fasta

import exports
import multigb
from reference import is_fasta_file
import theme
from uikit import plural, right_aligned


def source_records(path):
    """The sequences of the file `path` that can be added: [SimpleNamespace(name, length, load)] in file order, `load()` giving
    the Document of that one sequence (a big file is not read whole: only the records chosen are loaded, one at a time)."""
    path = Path(path)
    if is_fasta_file(path):
        fasta = Fasta(str(path))
        out = [SimpleNamespace(name=n, length=int(l), load=lambda n=n: Document.from_fasta(str(path), n)) for n, l in fasta.sequences()]
    elif multigb.is_multi_genbank(path):
        out = [SimpleNamespace(name=n, length=int(l), load=lambda n=n: multigb.read_record(path, n))
               for n, l, _ in multigb.records(path)]
    else:
        doc = exports.open_project(path)[0] if path.suffix.lower() == exports.PROJECT_EXT else Document(str(path))
        out = [SimpleNamespace(name=doc.name or path.stem, length=doc.length, load=lambda: doc)]
    if not out:
        raise OSError("no sequence in this file")
    return out


def current_names(editor):
    """The record names the open file already has (the LOCUS names of a GenBank file, the sequence names of a FASTA file)."""
    if editor.fasta_records:
        return [n for n, _ in editor.fasta_records]
    try:
        return [n for n, *_ in multigb.records(editor.path)] or [editor.doc.name]
    except OSError:
        return [editor.doc.name]


class AddRecordsDialog(QDialog):
    """Which sequences of the chosen file to add, and under which names (each must be new in the open file)."""

    def __init__(self, parent, target, fasta, records, taken):
        super().__init__(parent)
        self.setWindowTitle("Import as New Chromosome")
        self.records, self.taken = records, set(taken)
        kind = "sequences" if fasta else "records"
        intro = QLabel(wordWrap=True)
        intro.setText(
            f"The ticked sequences are added to <b>{target.name}</b> as new {kind}, after the "
            f"{plural(len(self.taken), kind[:-1])} it already has. The file is written now: the new {kind} do not need Save, and Undo does not "
            "remove them (the sequence open now keeps its unsaved edits until you save it)."
            + ("<br><br>A FASTA file holds only bases: features and primers of the added sequences are not kept. Save the file as "
               "GenBank first to keep them." if fasta else ""))
        self.table = QTableWidget(len(records), 2)
        self.table.setHorizontalHeaderLabels(["Name in the file", "Length"])
        self.table.verticalHeader().hide()
        self.table.setSelectionMode(QAbstractItemView.NoSelection)
        self.table.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
        self.table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeToContents)
        names = set(self.taken)
        for row, r in enumerate(records):
            name = multigb.record_name(r.name, names)   # a name already in the file gets "_2"
            names.add(name)
            item = QTableWidgetItem(name)
            item.setFlags(item.flags() | Qt.ItemIsUserCheckable | Qt.ItemIsEditable)
            item.setCheckState(Qt.Checked)
            item.setToolTip(f"{r.name} in the chosen file; double-click to rename")
            self.table.setItem(row, 0, item)
            length = QTableWidgetItem(f"{r.length:,} bp")
            length.setFlags(Qt.ItemIsEnabled)
            right_aligned(length)
            self.table.setItem(row, 1, length)
        self.problem = QLabel(styleSheet=f"color: {theme.hex('status_bad')}", wordWrap=True)
        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.button(QDialogButtonBox.Ok).setText(f"&Add to {target.name}")
        buttons.accepted.connect(self._accept)
        buttons.rejected.connect(self.reject)
        layout = QVBoxLayout(self)
        layout.addWidget(intro)
        layout.addWidget(self.table)
        layout.addWidget(self.problem)
        layout.addWidget(buttons)
        self.resize(520, min(560, 260 + 30 * len(records)))

    def chosen(self):
        """[(record, name)] of the ticked rows (the names as typed, spaces made "_")."""
        out = []
        for row, r in enumerate(self.records):
            item = self.table.item(row, 0)
            if item.checkState() == Qt.Checked:
                out.append((r, multigb.record_name(item.text())))
        return out

    def check(self):
        """'' if the ticked names can be added, else what is wrong."""
        chosen = self.chosen()
        if not chosen:
            return "Tick at least one sequence."
        seen = set()
        for _, name in chosen:
            if name in self.taken:
                return f"The file already has a sequence named {name}: rename it."
            if name in seen:
                return f"Two sequences are named {name}: rename one."
            seen.add(name)
        return ""

    def _accept(self):
        problem = self.check()
        self.problem.setText(problem)
        if not problem:
            self.accept()
