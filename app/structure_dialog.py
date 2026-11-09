"""The 3D structure window: choose where the structure comes from (see structure.py), get it in a worker thread, open it in the browser
as a self-contained 3Dmol.js page, or save it as a PDB file. It says before asking which service the sequence goes to."""
import threading
from html import escape
from pathlib import Path
from types import SimpleNamespace

from PySide6.QtCore import QObject, QSettings, QUrl, Signal
from PySide6.QtGui import QDesktopServices
from PySide6.QtWidgets import (
    QComboBox, QDialog, QFileDialog, QFormLayout, QHBoxLayout, QLabel, QLineEdit, QMessageBox, QProgressBar, QPushButton,
    QVBoxLayout,
)

import structure
from primer_db import data_dir
import theme
import wip


class _Relay(QObject):
    """Brings the answer from the worker thread to the GUI thread (lives forever: an answer arriving after the window closed is harmless)."""
    done = Signal(object, object, object)   # token, Structure or None, error text


_relay = _Relay()


class StructureDialog(QDialog):
    def __init__(self, parent, name, protein, accession=""):
        super().__init__(parent)
        self.name, self.protein, self.result, self.token = name, protein, None, None
        self.job = None   # the request under way: its local command can be stopped (see cancel)
        self.settings = QSettings(str(data_dir() / "structure.ini"), QSettings.IniFormat)
        self.setWindowTitle(f"3D Structure of {name}")
        self.resize(560, 330)
        self.source = QComboBox()
        self.source.addItems(structure.SOURCES)
        self.source.setCurrentText(self.settings.value("source", "Auto"))
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
        layout.addWidget(self.status, 1)
        self.source.currentTextChanged.connect(self._source_changed)
        _relay.done.connect(self._done)
        self._source_changed(self.source.currentText())

    def _source_changed(self, source):
        self.notice.setText(f"<i>{structure.NOTICE[source]}</i>")
        self.accession.setEnabled(source in ("Auto", "AlphaFold DB"))
        self.command.setEnabled(source == "Local command")

    # ---- getting it -------------------------------------------------------------------------------------------------------------
    def start(self):
        source = self.source.currentText()
        self.settings.setValue("source", source)
        self.settings.setValue("command", self.command.text())
        self.token = token = object()
        self.job = job = SimpleNamespace(proc=None, cancelled=False)
        self.get.setEnabled(False)
        self._working(True)
        self.status.setText("Working… (a prediction can take a minute)")
        protein, accession, command = self.protein, self.accession.text(), self.command.text()

        def work():
            try:
                found, error = structure.get(protein, source, accession, command, job=job), None
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
        self._working(False)
        if found is None:
            self.status.setText(f"<span style='color:{theme.hex('status_bad')}'>{escape(error).replace(chr(10), '<br>')}</span>")
            return
        self.result = found
        self.open.setEnabled(True)
        self.save.setEnabled(True)
        confidence = ("Colours show the confidence of each part (pLDDT): the blue parts are reliable, the orange ones are not." if found.plddt
                      else "")
        self.status.setText(f"<b>{found.detail}</b>{' (from the cache)' if found.cached else ''}<br>{confidence}")

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
