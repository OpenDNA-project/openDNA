"""The 3D structure window: choose where the structure comes from (see structure.py), get it in a worker thread (the Protein panel
then shows it in 3D: structure_view.py), open it in the browser as a self-contained 3Dmol.js page, or save it as a PDB file. It says before asking which service the sequence goes to.
With the PDB and Linked entries sources it first lists the entries (title, identity, method, resolution, organism) and you pick the one
to show; Linked entries is chosen to begin with when the feature links PDB entries."""
import threading
from html import escape
from pathlib import Path
from types import SimpleNamespace

from PySide6.QtCore import QObject, QSettings, Qt, QUrl, Signal
from PySide6.QtGui import QDesktopServices
from PySide6.QtWidgets import (
    QAbstractItemView, QComboBox, QDialog, QFileDialog, QFormLayout, QHBoxLayout, QHeaderView, QLabel, QLineEdit, QMessageBox,
    QProgressBar, QPushButton, QTableWidget, QTableWidgetItem, QVBoxLayout,
)

import structure
from primer_db import data_dir
import theme
import wip


HIT_COLUMNS = ["Entry", "Identity", "Covers", "Title", "Organism", "Method", "Resolution", "Chains"]


class _Relay(QObject):
    """Brings the answer from the worker thread to the GUI thread (lives forever: an answer arriving after the window closed is harmless)."""
    done = Signal(object, object, object)   # token, Structure / [PdbHit] / None, error text


_relay = _Relay()


class StructureDialog(QDialog):
    found = Signal(str, object)   # the amino acids, the Structure got for them (the Protein panel shows it in 3D)

    def __init__(self, parent, name, protein, accession="", links=()):
        super().__init__(parent)
        self.name, self.protein, self.result, self.token = name, protein, None, None
        self.links = list(links)   # the PDB codes the feature links (structure.pdb_links)
        self.job = None   # the request under way: its local command can be stopped (see cancel)
        self.settings = QSettings(str(data_dir() / "structure.ini"), QSettings.IniFormat)
        self.setWindowTitle(f"3D Structure of {name}")
        self.resize(560, 330)
        self.source = QComboBox()
        self.source.addItems(structure.SOURCES)
        self.source.setCurrentText("Linked entries" if self.links else self.settings.value("source", "Auto"))
        self.accession = QLineEdit(accession, placeholderText="UniProt accession, e.g. P69905 (for AlphaFold DB)")
        self.command = QLineEdit(self.settings.value("command", ""),
                                 placeholderText="e.g. colabfold_batch {fasta} outdir …  (must write {pdb})")
        self.notice = QLabel(wordWrap=True)
        self.status = QLabel(wordWrap=True)
        self.get = QPushButton("Get structure", clicked=self.start, default=True)
        self.open = QPushButton("Open 3D view in the browser", clicked=self.open_view, enabled=False)
        self.save = QPushButton("Save PDB…", clicked=self.save_pdb, enabled=False)
        # While a structure is being got: a busy bar (a prediction can take a minute, a local command up to an hour) and a
        # Cancel button that stops a local command and forgets the answer of a request.
        self.busy = QProgressBar(minimum=0, maximum=0, textVisible=False, visible=False)
        # The PDB entries found (PDB source): pick one, then "Show this entry" (or double-click it)
        self.hits = []
        self.table = QTableWidget(0, len(HIT_COLUMNS), visible=False, editTriggers=QAbstractItemView.NoEditTriggers,
                                  selectionBehavior=QAbstractItemView.SelectRows, selectionMode=QAbstractItemView.SingleSelection)
        self.table.setHorizontalHeaderLabels(HIT_COLUMNS)
        self.table.verticalHeader().hide()
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.table.horizontalHeader().setSectionResizeMode(HIT_COLUMNS.index("Title"), QHeaderView.Stretch)
        self.table.itemSelectionChanged.connect(lambda: self.pick.setEnabled(bool(self.table.selectedItems()) and self.job is None))
        self.table.itemDoubleClicked.connect(lambda _: self.show_hit())
        self.pick = QPushButton("Show this entry", clicked=self.show_hit, enabled=False, visible=False)
        self.cancel_button = QPushButton("Cancel", clicked=self.cancel, visible=False)
        form = QFormLayout()
        title = QHBoxLayout()   # the protein, and the "WIP" tag while the protein interface is unfinished (wip.FEATURES)
        title.addWidget(QLabel(f"{name} · {len(protein):,} amino acids"))
        title.addStretch()
        title.addWidget(wip.Badge("protein"))
        form.addRow("Protein", title)
        form.addRow("Source", self.source)
        form.addRow("Accession", self.accession)
        form.addRow("Local command", self.command)
        row = QHBoxLayout()
        for w in (self.get, self.open, self.save):
            row.addWidget(w)
        working = QHBoxLayout()
        working.addWidget(self.busy, 1)
        working.addWidget(self.cancel_button)
        layout = QVBoxLayout(self)
        layout.addLayout(form)
        layout.addWidget(self.notice)
        layout.addLayout(row)
        layout.addLayout(working)
        layout.addWidget(self.status)
        layout.addWidget(self.table, 1)
        layout.addWidget(self.pick, 0, Qt.AlignRight)
        self.source.currentTextChanged.connect(self._source_changed)
        _relay.done.connect(self._done)
        self._source_changed(self.source.currentText())

    def _source_changed(self, source):
        linked = (f" The feature links {len(self.links)} entr{'y' if len(self.links) == 1 else 'ies'} ({', '.join(self.links[:4])}"
                  f"{'…' if len(self.links) > 4 else ''})." if source == "Linked entries" and self.links else "")
        self.notice.setText(f"<i>{structure.NOTICE[source]}{escape(linked)}</i>")
        self.accession.setEnabled(source in ("Auto", "AlphaFold DB", "Linked entries"))
        self.command.setEnabled(source == "Local command")
        self.get.setText({"PDB": "Search PDB", "Linked entries": "List entries"}.get(source, "Get structure"))

    # ---- getting it -------------------------------------------------------------------------------------------------------------
    def start(self):
        source = self.source.currentText()
        self.settings.setValue("source", source)
        self.settings.setValue("command", self.command.text())
        protein, accession, command = self.protein, self.accession.text(), self.command.text()
        if source == "PDB":   # a list to pick from: see show_hit
            self._run(lambda job: structure.search_pdb(protein), "Searching the PDB…")
        elif source == "Linked entries":
            links = self.links
            if not links and not accession.strip():
                self.status.setText("The feature links no PDB entry (db_xref PDB:…) and there is no UniProt accession: type one, or "
                                    "use PDB to search by sequence.")
                return
            self._run(lambda job: structure.linked_pdb(protein, links, accession), "Looking up the linked entries…")
        else:
            self._run(lambda job: structure.get(protein, source, accession, command, job=job), "Working… (a prediction can take a minute)")

    def show_hit(self):
        """Get the PDB entry picked in the list (and keep it in the cache: Auto then finds this one for this protein)."""
        rows = self.table.selectionModel().selectedRows()
        if not rows or self.job is not None:
            return
        hit, protein = self.hits[rows[0].row()], self.protein

        def work(job):
            found = structure.pdb_entry(hit)
            structure.store(protein, "", found)
            return found
        self._run(work, f"Getting PDB {hit.entry}…")

    def _run(self, task, message):
        """Run `task(job)` in a worker thread; its answer (a Structure, or the [PdbHit] of a search) comes to _done."""
        self.token = token = object()
        self.job = job = SimpleNamespace(proc=None, cancelled=False)
        self.get.setEnabled(False)
        self.pick.setEnabled(False)
        self._working(True)
        self.status.setText(message)

        def work():
            try:
                found, error = task(job), None
            except structure.StructureError as e:
                found, error = None, str(e)
            except Exception as e:   # a bug must not leave the button disabled for good
                found, error = None, f"Unexpected error: {e}"
            try:
                _relay.done.emit(token, found, error)
            except RuntimeError:
                pass   # the program is closing: the relay is already gone
        threading.Thread(target=work, daemon=True).start()

    def _working(self, on):
        self.busy.setVisible(on)
        self.cancel_button.setVisible(on)

    def cancel(self):
        """Stop waiting: a local command is stopped, the answer of a request to a service is not used when it comes."""
        job, self.job, self.token = self.job, None, None
        if job is not None:
            job.cancelled = True
            if job.proc is not None and job.proc.poll() is None:
                job.proc.terminate()
        self.get.setEnabled(True)
        self.pick.setEnabled(bool(self.table.selectedItems()))
        self._working(False)
        self.status.setText("Cancelled.")

    def closeEvent(self, event):
        if self.job is not None:   # a local command does not run on for an hour behind a closed window
            self.cancel()
        super().closeEvent(event)

    def _done(self, token, found, error):
        if token is not self.token:
            return
        self.job = None
        self.get.setEnabled(True)
        self.pick.setEnabled(bool(self.table.selectedItems()))
        self._working(False)
        if isinstance(found, list):
            self._list_hits(found)
            return
        if found is None:
            self.status.setText(f"<span style='color:{theme.hex('status_bad')}'>{escape(error).replace(chr(10), '<br>')}</span>")
            return
        self.result = found
        self.found.emit(self.protein, found)
        self.open.setEnabled(True)
        self.save.setEnabled(True)
        confidence = ("Colours show the confidence of each part (pLDDT): the blue parts are reliable, the orange ones are not." if found.plddt
                      else "")
        self.status.setText(f"<b>{found.detail}</b>{' (from the cache)' if found.cached else ''}<br>{confidence}")

    def _list_hits(self, hits):
        self.hits = hits
        self.table.setRowCount(0)
        self.table.setVisible(bool(hits))
        self.pick.setVisible(bool(hits))
        if not hits:
            self.status.setText("No PDB entry is linked to this protein." if self.source.currentText() == "Linked entries" else
                                "No entry of the PDB has a chain like this protein (50 % identity or more).")
            return
        for h in hits:
            row = self.table.rowCount()
            self.table.insertRow(row)
            cells = [h.entry, f"{h.identity:.0%}", f"{h.covered:.0%}", h.title, h.organism, h.method.lower(),
                     f"{h.resolution:.2f} Å" if h.resolution else "", h.chains]
            for column, text in enumerate(cells):
                item = QTableWidgetItem(text)
                item.setToolTip(text)
                self.table.setItem(row, column, item)
        self.table.selectRow(0)
        self.resize(max(self.width(), 900), max(self.height(), 560))
        self.status.setText(f"{len(hits)} entries, the closest first. <b>Identity</b> is over the aligned part, <b>Covers</b> how much of "
                            "the protein the chain spans: pick the entry to show (double-click, or <i>Show this entry</i>).")

    # ---- using it ---------------------------------------------------------------------------------------------------------------
    def view_file(self):
        """Write the 3D page of the structure in hand (next to the cache) and return its path."""
        path = structure.cache_key(self.protein, "", "view").with_suffix(".html")
        path.write_text(structure.viewer_html(self.result, f"{self.name} · {len(self.protein):,} aa"), encoding="utf-8")
        return path

    def open_view(self):
        if self.result is not None:
            try:
                page = self.view_file()
            except OSError as e:   # a full disk, a cache folder that cannot be written: said, not a traceback
                QMessageBox.warning(self, "3D structure", f"The 3D page could not be written:\n{e}")
                return
            QDesktopServices.openUrl(QUrl.fromLocalFile(str(page)))

    def save_pdb(self):
        if self.result is None:
            return
        path, _ = QFileDialog.getSaveFileName(self, "Save structure", f"{self.name}.pdb", "PDB files (*.pdb)")
        if path:
            try:
                Path(path).write_text(self.result.pdb)
            except OSError as e:
                QMessageBox.warning(self, "Save structure", f"{path}\n\n{e}")
                return
            self.status.setText(f"Saved {escape(path)}")
