"""File > Import from NCBI: search the Nucleotide or Protein database, pick a record, open it (a nucleotide record in a new window with
its annotation; a protein record in the 3D structure window) or save it as a file. The network work runs in a worker thread.

File > Import > Background Genome from NCBI (`GenomeDialog`): search NCBI Datasets for an organism's assemblies, download one
(sequences + GFF3 annotation), open it (the FASTA with its annotation as a track) or make it the background genome of the primer
and spacer off-target checks."""
import threading
from html import escape
from pathlib import Path

from PySide6.QtCore import QObject, QSettings, Qt, Signal
from PySide6.QtWidgets import (
    QAbstractItemView, QCheckBox, QComboBox, QDialog, QFileDialog, QFormLayout, QHBoxLayout, QLabel, QLineEdit, QMessageBox, QProgressBar,
    QPushButton, QTableWidget, QTableWidgetItem, QVBoxLayout,
)

import ncbi
import online
import validate
from primer_db import data_dir, notifier
from structure_dialog import StructureDialog
from uikit import fit_columns, no_auto_default, plural, right_aligned
import value_filter
import theme

COLUMNS = ["Accession", "Title", "Length", "Organism", "Updated"]


class _Relay(QObject):
    """Brings the answer of the worker thread to the GUI thread (lives forever: an answer arriving after the window closed is harmless)."""
    done = Signal(object, object, object)   # token, result, error text
    progress = Signal(object, object)       # token, bytes downloaded so far


_relay = _Relay()


class Cancelled(Exception):
    """The user stopped a download (Cancel, or the window was closed)."""


def error_html(text):
    """An error for the status line (rich text): in red, its own text escaped (a "<" in an NCBI message is not markup)."""
    return f"<span style='color:{theme.hex('status_bad')}'>{escape(text)}</span>"


class NcbiDialog(QDialog):
    def __init__(self, window):
        super().__init__(window)
        self.window, self.records, self.token, self.opened, self.pending = window, [], None, None, None
        self.settings = QSettings(str(data_dir() / "ncbi.ini"), QSettings.IniFormat)
        self.setWindowTitle("Import from NCBI")
        self.resize(900, 520)
        self.database = QComboBox()
        self.database.addItems(ncbi.DATABASES)
        self.query = QLineEdit(placeholderText="An accession (NC_000913.3, U00096.3, P69905) or words (pUC19 cloning vector)",
                               clearButtonEnabled=True)
        self.query.returnPressed.connect(self.search)   # the one way Enter searches (a default button too made it search twice)
        self.find = QPushButton("Search", clicked=self.search)
        online.mark(self.find, "your search words (or accession) to NCBI (eutils.ncbi.nlm.nih.gov)")
        self.table = QTableWidget(0, len(COLUMNS))
        self.table.setHorizontalHeaderLabels(COLUMNS)
        value_filter.install(self.table, ["Organism"])   # right-click a header: filter by values
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SingleSelection)
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        fit_columns(self.table, stretch=1)   # the title takes the room; a long one is elided, whole in its tooltip
        self.table.itemSelectionChanged.connect(self._selected)
        self.table.cellDoubleClicked.connect(lambda *_: self.open())
        self.key = QLineEdit(self.settings.value("api_key", ""), placeholderText="optional: 10 requests a second instead of 3")
        self.key.setEchoMode(QLineEdit.Password)
        self.key.setMaxLength(100)
        self.query.setMaxLength(500)
        self.status = QLabel("Your search words (or accession) are sent to NCBI (eutils.ncbi.nlm.nih.gov); nothing else leaves this computer.",
                             wordWrap=True)
        self.open_button = QPushButton("Open", clicked=self.open, enabled=False)
        self.save_button = QPushButton("Save as file…", clicked=self.save, enabled=False)
        top = QHBoxLayout()
        top.addWidget(self.database)
        top.addWidget(self.query, 1)
        top.addWidget(self.find)
        keys = QFormLayout()
        keys.addRow("NCBI API key", self.key)
        row = QHBoxLayout()
        row.addWidget(self.status, 1)
        row.addWidget(self.open_button)
        row.addWidget(self.save_button)
        layout = QVBoxLayout(self)
        layout.addLayout(top)
        layout.addWidget(self.table, 1)
        layout.addLayout(keys)
        layout.addLayout(row)
        no_auto_default(self)
        _relay.done.connect(self._done)
        # A wrong character in the key is named as soon as it is typed or pasted (and nothing is sent while it is there).
        self.key.textChanged.connect(lambda text: self.status.setText(validate.api_key_problem(text.strip()) or self.privacy))
        self.privacy = self.status.text()

    # ---- searching --------------------------------------------------------------------------------------------------------------
    def _run(self, work, what):
        """Run `work()` in a thread; `_done` gets its result (or the error text) back."""
        self.key.setText(self.key.text().strip())      # a pasted key often carries a space or line break
        self.settings.setValue("api_key", self.key.text())
        self.settings.remove("email")   # an address typed in an older version: no longer asked for, not kept
        self.token = token = object()
        self.find.setEnabled(False)
        self.open_button.setEnabled(False)
        self.save_button.setEnabled(False)
        self.status.setText(what)

        def target():
            try:
                result, error = work(), None
            except ncbi.NcbiError as e:
                result, error = None, str(e)
            except Exception as e:   # a bug must not leave the buttons disabled for good
                result, error = None, f"The NCBI request failed: {e or type(e).__name__}"
            try:
                _relay.done.emit(token, result, error)
            except RuntimeError:
                pass   # the program is closing: the relay is already gone
        threading.Thread(target=target, daemon=True).start()

    def _key_ok(self):
        """The API key, trimmed (a pasted key often carries a space or line break) and checked: letters and digits only. A bad one
        is named in the status line and nothing is sent."""
        self.key.setText(self.key.text().strip())
        problem = validate.api_key_problem(self.key.text())
        if problem:
            self.status.setText(problem)
            self.key.setFocus()
        return not problem

    def search(self):
        if not self._key_ok():
            return
        db, term, key = ncbi.DATABASES[self.database.currentText()], self.query.text(), self.key.text()
        self._run(lambda: ("search", ncbi.search(db, term, key)), "Searching NCBI…")

    def _done(self, token, result, error):
        if token is not self.token:
            return
        self.find.setEnabled(True)
        if error:
            self.status.setText(error_html(error))
            self._selected()
            return
        kind, value = result
        if kind == "search":
            self._show(value)
        else:
            self._fetched(value)

    def _show(self, records):
        self.records = records
        self.table.setRowCount(len(records))
        for i, r in enumerate(records):
            unit = "aa" if r.db == "protein" else "bp"
            for c, v in enumerate([r.accession, r.title, f"{r.length:,} {unit}", r.organism, r.updated]):
                item = QTableWidgetItem(v)
                if c == 2:
                    right_aligned(item)
                if c == 1:
                    item.setToolTip(v)   # (elided when long)
                self.table.setItem(i, c, item)
        self.status.setText(f"{plural(len(records), 'record')}." if records else "Nothing found: try other words or the exact accession.")
        if records:
            self.table.selectRow(0)
        self._selected()

    def current(self):
        row = self.table.currentRow()
        return self.records[row] if 0 <= row < len(self.records) else None

    def _selected(self):
        r = self.current()
        self.open_button.setEnabled(r is not None)
        self.save_button.setEnabled(r is not None)
        self.open_button.setText("3D structure…" if r is not None and r.db == "protein" else "Open")

    # ---- using a record ---------------------------------------------------------------------------------------------------------
    def _fetch(self, then):
        """Fetch the selected record (from the cache if it is there), then call `then(record, path)`."""
        r = self.current()
        if r is None:
            return
        if r.db == "nuccore" and r.length > ncbi.WARN_BASES and not ncbi.path_of(r).exists():
            ask = QMessageBox.question(self, "Large record", f"{r.accession} has {r.length:,} bases: fetching it can take minutes. Continue?")
            if ask != QMessageBox.Yes:
                return
        if not self._key_ok():
            return
        self.pending = then
        key = self.key.text()
        self._run(lambda: ("fetched", (r, *ncbi.get(r, key))), f"Fetching {r.accession}…")

    def _fetched(self, value):
        record, path, cached = value
        self.status.setText(f"{record.accession}: {'from the cache' if cached else 'downloaded'} ({path.stat().st_size:,} bytes).")
        self._selected()
        self.pending(record, path)

    def open(self):
        self._fetch(self._open_record)

    def _open_record(self, record, path):
        from genomeedit_core import Document
        if record.db == "protein":
            name, protein = ncbi.protein_of(path)
            self.opened = StructureDialog(self.window, record.title.split(";")[0][:60] or name, protein, record.uniprot)
            self.opened.show()
            return
        try:
            doc = Document(str(path))
        except (OSError, ValueError) as e:
            QMessageBox.warning(self, "Cannot open", f"{record.accession}: {e}")
            return
        target = Path(f"{record.accession}.gb")
        self.opened = self.window.open_untitled(doc, target)
        self.accept()

    def save(self):
        self._fetch(self._save_record)

    def _save_record(self, record, path):
        name, _ = QFileDialog.getSaveFileName(self, "Save record", f"{record.accession}{path.suffix}",
                                              "GenBank (*.gb)" if path.suffix == ".gb" else "FASTA (*.fasta)")
        if name:
            try:
                Path(name).write_text(path.read_text(encoding="utf-8"), encoding="utf-8")
            except OSError as e:   # a read-only folder, a full disk: said, not a crash
                QMessageBox.critical(self, "Save record", f"{name}\n\n{e}")
                return
            self.status.setText(f"Saved {escape(name)}")


GENOME_COLUMNS = ["Accession", "Organism", "Assembly", "Level", "Size", "Sequences", "Released"]


class GenomeDialog(QDialog):
    """Background genomes from NCBI Datasets: search, download (kept in the data folder), open or use for off-target checks."""

    def __init__(self, window):
        super().__init__(window)
        self.window, self.assemblies, self.token, self.opened = window, [], None, None
        self.settings = QSettings(str(data_dir() / "ncbi.ini"), QSettings.IniFormat)
        self.setWindowTitle("Background Genome from NCBI")
        self.resize(980, 520)
        self.query = QLineEdit(placeholderText="An organism (Escherichia coli, Saccharomyces cerevisiae, 562) or an assembly accession "
                                               "(GCF_000005845.2)", clearButtonEnabled=True)
        self.query.returnPressed.connect(self.search)   # the one way Enter searches (a default button too made it search twice)
        self.reference_only = QCheckBox("NCBI reference genomes only", checked=True,
                                        toolTip="Only the assembly NCBI picked as the reference of each species")
        self.find = QPushButton("Search", clicked=self.search)
        self.table = QTableWidget(0, len(GENOME_COLUMNS))
        self.table.setHorizontalHeaderLabels(GENOME_COLUMNS)
        value_filter.install(self.table, ["Organism", "Assembly", "Level"])   # right-click a header: filter by values
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SingleSelection)
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        fit_columns(self.table, stretch=1)
        self.table.itemSelectionChanged.connect(self._selected)
        self.table.cellDoubleClicked.connect(lambda *_: self.open())
        self.status = QLabel("Your search words are sent to NCBI (api.ncbi.nlm.nih.gov); genomes are kept in "
                             f"{data_dir() / 'genomes'}.", wordWrap=True)
        self.open_button = QPushButton("Download and open", clicked=self.open, enabled=False)
        self.reference_button = QPushButton("Use for off-target checks", clicked=self.use_as_reference, enabled=False,
                                            toolTip="The genome the Primers and Spacer Finder tools check for other binding sites "
                                                    "(Primers > Options); one-sequence genomes only")
        # While NCBI answers: a moving bar (a download's size is not known beforehand), and for a download a Cancel.
        self.progress = QProgressBar(minimum=0, maximum=0, textVisible=False, maximumWidth=160)
        self.cancel_button = QPushButton("Cancel", clicked=self.cancel, toolTip="Stop the download (nothing half-written is kept)")
        for w in (self.progress, self.cancel_button):
            keep = w.sizePolicy()
            keep.setRetainSizeWhenHidden(True)   # the buttons beside them do not move when they show
            w.setSizePolicy(keep)
            w.hide()
        self.stop = threading.Event()
        top = QHBoxLayout()
        top.addWidget(self.query, 1)
        top.addWidget(self.reference_only)
        top.addWidget(self.find)
        row = QHBoxLayout()
        row.addWidget(self.status, 1)
        row.addWidget(self.progress)
        row.addWidget(self.cancel_button)
        row.addWidget(self.reference_button)
        row.addWidget(self.open_button)
        layout = QVBoxLayout(self)
        layout.addLayout(top)
        layout.addWidget(self.table, 1)
        layout.addLayout(row)
        no_auto_default(self)
        _relay.done.connect(self._done)
        _relay.progress.connect(self._progress)

    def _key(self):
        return self.settings.value("api_key", "")   # set in Import from NCBI

    def _run(self, work, what, cancellable=False):
        self.token = token = object()
        self.stop.set()   # a download still running for an older request stops
        stop = self.stop = threading.Event()
        for b in (self.find, self.open_button, self.reference_button):
            b.setEnabled(False)
        self.status.setText(what)
        self.progress.show()
        self.cancel_button.setVisible(cancellable)

        def target():
            try:
                result, error = work(token, stop), None
            except Cancelled:
                result, error = ("cancelled", None), None
            except ncbi.NcbiError as e:
                result, error = None, str(e)
            except OSError as e:
                result, error = None, f"Could not write the genome: {e}"
            except Exception as e:   # a bug must not leave the buttons disabled for good
                result, error = None, f"The NCBI request failed: {e or type(e).__name__}"
            try:
                _relay.done.emit(token, result, error)
            except RuntimeError:
                pass
        threading.Thread(target=target, daemon=True).start()

    def cancel(self):
        """Stop the download: the download loop sees the flag at its next block and removes what it wrote."""
        self.stop.set()
        self.cancel_button.setEnabled(False)
        self.status.setText("Cancelling the download…")

    def done(self, result):   # closed (Close, Escape, the window's ✕): a download still running stops too
        self.stop.set()
        super().done(result)

    def search(self):
        term, only, key = self.query.text(), self.reference_only.isChecked(), self._key()
        self._run(lambda _token, _stop: ("search", ncbi.search_genomes(term, only, key)), "Searching NCBI…")

    def _progress(self, token, got):
        if token is self.token:
            self.status.setText(f"Downloading… {got / 1e6:,.1f} MB")

    def _done(self, token, result, error):
        if token is not self.token:
            return
        self.find.setEnabled(True)
        self.progress.hide()
        self.cancel_button.hide()
        self.cancel_button.setEnabled(True)
        if error:
            self.status.setText(error_html(error))
            self._selected()
            return
        kind, value = result
        if kind == "search":
            self._show(value)
        elif kind == "cancelled":
            self.status.setText("Download cancelled.")
            self._selected()
        else:
            value[0](*value[1:])   # then(assembly, fasta, gff, cached)

    def _show(self, assemblies):
        self.assemblies = assemblies
        self.table.setRowCount(len(assemblies))
        for i, a in enumerate(assemblies):
            values = [a.accession, a.organism, a.name, a.level, f"{a.length / 1e6:,.2f} Mb", f"{a.sequences:,}", a.released]
            for c, v in enumerate(values):
                item = QTableWidgetItem(v)
                if c in (4, 5):
                    item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
                if a.category:
                    item.setToolTip(a.category)
                self.table.setItem(i, c, item)
        hint = " Untick “NCBI reference genomes only” to see every assembly." if self.reference_only.isChecked() else ""
        self.status.setText(f"{plural(len(assemblies), 'assembly')}." if assemblies
                            else "Nothing found: check the spelling, or try the taxonomy number." + hint)
        if assemblies:
            self.table.selectRow(0)
        self._selected()

    def current(self):
        row = self.table.currentRow()
        return self.assemblies[row] if 0 <= row < len(self.assemblies) else None

    def _selected(self):
        a = self.current()
        self.open_button.setEnabled(a is not None)
        self.reference_button.setEnabled(a is not None)

    def _download(self, then, index=False):
        a = self.current()
        if a is None:
            return
        folder = ncbi.genome_folder(a)
        if a.length > ncbi.GENOME_WARN_BASES and not (folder / f"{folder.name}.fna").exists():
            if QMessageBox.question(self, "Large genome", f"{a.accession} has {a.length / 1e6:,.0f} Mb: the download is large "
                                    "and can take a long time. Continue?") != QMessageBox.Yes:
                return
        key = self._key()

        def work(token, stop):   # `stop`: this request's flag (Cancel, closing the window or a newer request sets it)
            def progress(got):
                if stop.is_set():
                    raise Cancelled()   # download_genome removes the half-written files on the way out
                _relay.progress.emit(token, got)
            fasta, gff, cached = ncbi.download_genome(a, key, progress)
            if index:   # index every sequence here, on the worker thread (a big genome takes a while)
                from genomeedit_core import Fasta
                Fasta(str(fasta))
            return "downloaded", (then, a, fasta, gff, cached)
        self._run(work, f"Downloading {a.accession}…", cancellable=True)

    def open(self):
        self._download(self._open_genome)

    def _open_genome(self, assembly, fasta, gff, cached):
        self.status.setText(f"{assembly.accession}: {'already downloaded' if cached else 'downloaded'} ({fasta.parent}).")
        self._selected()
        win = self.window.open_path(fasta)
        if win is not None and gff is not None:
            win.load_track(str(gff))
        self.opened = win
        if win is not None:
            self.accept()

    def use_as_reference(self):
        self._download(self._set_reference, index=True)

    def _set_reference(self, assembly, fasta, gff, cached):
        """The off-target checks read the downloaded FASTA itself: every sequence of it (all chromosomes) is searched."""
        from primer_db import settings
        from reference import Reference
        try:
            ref = Reference(str(fasta))   # the .fai index was built by the download step
        except (OSError, ValueError) as e:
            QMessageBox.warning(self, "Background Genome", f"{assembly.accession}: {e}")
            return
        s = settings()
        s.setValue("primers/genome", str(fasta))
        s.setValue("primers/genome_check", "true")
        import reference_picker
        reference_picker.remember(fasta)
        notifier.changed.emit()   # every window: the toolbar's reference box, and the primers re-checked
        if getattr(self.window, "editor", None) is not None:
            self.window.primer_options_changed()
        self._selected()
        self.status.setText(f"{assembly.organism} ({assembly.accession}) is now the background genome of the off-target checks: "
                            f"{ref.label()}.")
