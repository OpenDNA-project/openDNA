"""Import > Tracks > Track Hub: the dialog that reads a UCSC track hub (a URL or a local folder), lists its tracks in the hub's own
groups, all ticked (the hidden and grouped ones too: those the hub hides at first load hidden), with where the files will be
copied, and then copies them: several at once, one row each with its bar and its own cancel, Cancel stopping all of them at once
even when a server does not answer. main.import_hub loads what was copied; the hub itself is read by track_hub."""
from pathlib import Path

from PySide6.QtCore import Qt, QTimer
from PySide6.QtWidgets import (
    QComboBox, QDialog, QDialogButtonBox, QFileDialog, QFormLayout, QHBoxLayout, QLabel, QLineEdit, QMessageBox, QProgressBar,
    QPushButton, QRadioButton, QTreeWidget, QTreeWidgetItem, QVBoxLayout, QWidget,
)

import track_hub
from uikit import ask, plural, run_busy
from user_error import UserError


class HubDialog(QDialog):
    def __init__(self, parent, editor):
        super().__init__(parent)
        self.setWindowTitle("Import a track hub")
        self.editor, self.hub, self.tracks = editor, None, []
        self.transfers, self.genome_transfer = [], None   # the copies of the files, once begun
        self.copied, self.genome_got = [], None   # what main.import_hub loads: [(track, path or None, error)], and copy_genome's answer
        self.source = QLineEdit(placeholderText="https://…/hub.txt, or a folder holding hub.txt", clearButtonEnabled=True)
        self.source.returnPressed.connect(self.read)
        folder_btn = QPushButton("Folder…", clicked=self._pick_folder)
        read_btn = QPushButton("Read", clicked=self.read)
        row = QHBoxLayout()
        row.addWidget(self.source, 1)
        row.addWidget(folder_btn)
        row.addWidget(read_btn)
        self.genome = QComboBox()
        self.genome.currentTextChanged.connect(self._list)
        self.tree = QTreeWidget(headerHidden=False)
        self.tree.setHeaderLabels(["Track", "Type", "When it opens"])
        self.tree.itemChanged.connect(lambda *_: self._count())
        self.summary = QLabel(wordWrap=True)
        # Where the tracks go: on the sequence that is open, or on the hub's own genome, fetched and opened with them.
        self.on_file = QRadioButton("Show the tracks on the open file", checked=True)
        self.with_genome = QRadioButton("Import the hub's genome with the tracks (opens in a new tab, saved as a new .oDNA project)")
        self.on_file.toggled.connect(lambda *_: self._folder_state())
        self.folder = QLineEdit()
        self.folder.textChanged.connect(lambda *_: self._count())
        browse = QPushButton("Browse…", clicked=self._pick_target)
        target = QHBoxLayout()
        target.addWidget(self.folder, 1)
        target.addWidget(browse)
        form = QFormLayout()
        form.addRow("Hub:", row)
        form.addRow("Genome:", self.genome)
        form.addRow("Tracks:", self.on_file)
        form.addRow("", self.with_genome)
        form.addRow("Copy the files to:", target)
        self.buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        self.buttons.button(QDialogButtonBox.Ok).setText("Load Tracks")
        self.on_file.toggled.connect(lambda on: self.buttons.button(QDialogButtonBox.Ok).setText("Load Tracks" if on else "Import Genome and Tracks"))
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)
        self.retry = self.buttons.addButton("Retry Failed", QDialogButtonBox.ActionRole)
        self.retry.clicked.connect(self._retry)
        self.retry.hide()
        self.rows = QTreeWidget(headerHidden=False, rootIsDecorated=False)
        self.rows.setHeaderLabels(["File", "State", "", ""])
        self.rows.hide()
        self.timer = QTimer(self, interval=100, timeout=self._tick)
        self.top = QWidget()
        self.top.setLayout(form)
        form.setContentsMargins(0, 0, 0, 0)
        layout = QVBoxLayout(self)
        layout.addWidget(self.top)
        layout.addWidget(self.tree, 1)
        layout.addWidget(self.rows, 1)
        layout.addWidget(self.summary)
        layout.addWidget(self.buttons)
        self.resize(720, 560)
        self._count()

    # ---- copying ------------------------------------------------------------------------------------------------------------
    def accept(self):
        if self.transfers:   # the copies are over (some failed or were stopped): load what came
            return self._finish()
        folder, ticked = self.target(), self.chosen()
        genome = self.genome.currentText() if self.with_genome.isChecked() else None
        if genome and not self._confirm_genome(genome, len(ticked), folder):
            return
        if folder is None:   # a hub on this computer: its files are read where they are
            self.copied = [(t, Path(t.url), "") for t in ticked]
            return super().accept()
        try:
            tracks = track_hub.track_transfers(ticked, folder)
        except OSError as e:
            self.summary.setText(f"cannot copy to {folder}: {e}")
            return
        if genome:
            self.genome_transfer = track_hub.Transfer(f"Genome {genome}", lambda t: track_hub.copy_genome(
                self.hub, genome, folder, lambda d, n, name: t.update(d, n, f"converting {name}" if name else "")))
        self.transfers = ([self.genome_transfer] if genome else []) + tracks
        self.top.setEnabled(False)
        self.tree.hide()
        self.rows.show()
        self.rows.clear()
        for i, t in enumerate(self.transfers):
            item = QTreeWidgetItem([t.label, "waiting"])
            self.rows.addTopLevelItem(item)
            self.rows.setItemWidget(item, 2, QProgressBar(textVisible=False, maximumHeight=12))
            self.rows.setItemWidget(item, 3, QPushButton("✕", maximumWidth=28, toolTip="Cancel this file", clicked=lambda _=False, i=i: self.transfers[i].cancel()))
        self.rows.setColumnWidth(0, 300)
        self.rows.setColumnWidth(1, 220)
        self.rows.setColumnWidth(2, 140)
        self._running()
        track_hub.start(self.transfers)

    def _running(self):
        self.buttons.button(QDialogButtonBox.Ok).setEnabled(False)
        self.buttons.button(QDialogButtonBox.Cancel).setText("Cancel All")
        self.retry.hide()
        self.timer.start()

    def reject(self):
        if self.timer.isActive():   # copying: stop every file at once (a silent server is left behind, unseen)
            for t in self.transfers:
                t.cancel()
            return
        super().reject()

    def _confirm_genome(self, genome, count, folder):
        """Say how big the genome download is, and ask."""
        source = track_hub.genome_source(self.hub, genome)
        try:
            size = run_busy(self, "Looking up the genome…", lambda: track_hub.genome_size(source))
        except OSError as e:
            QMessageBox.warning(self, "Import a track hub", f"The sequence of {genome} cannot be reached at\n{source}\n\n{e}")
            return False
        how_big = f"about {size / 1e6:,.0f} MB" if size else "of a size the server does not say"
        return ask(self, "Import a track hub", f"Import the genome {genome} with {plural(count, 'track')}?", [("go", "Import")],
                   detail=f"It is {how_big}, from {source}. It is written as a FASTA file (larger) in {folder}, "
                          "with the tracks, and opened in a new tab.") == "go"

    def _retry(self):
        for i, t in enumerate(self.transfers):
            if t.state in ("failed", "cancelled"):
                self.transfers[i] = t.again()
                if t is self.genome_transfer:
                    self.genome_transfer = self.transfers[i]
        self._running()
        track_hub.start([t for t in self.transfers if t.state == "waiting"])

    def _tick(self):
        """Show every row; when all are over, load them if all came, else say what failed (Retry Failed, or load the rest)."""
        for i, t in enumerate(self.transfers):
            item = self.rows.topLevelItem(i)
            bar = self.rows.itemWidget(item, 2)
            busy = t.state == "copying" and not t.total   # (a size the server does not say: a moving bar)
            bar.setRange(0, 0 if busy else t.total if t.state == "copying" else 1)
            bar.setValue(t.done if t.state == "copying" else 1 if t.state == "done" else 0)
            mb = f"{t.done / 1e6:,.1f} of {t.total / 1e6:,.1f} MB" if t.total else f"{t.done / 1e6:,.1f} MB" if t.done else "starting…"
            item.setText(1, {"waiting": "waiting", "copying": t.text or mb, "done": "done",
                             "failed": f"failed: {t.error}", "cancelled": "cancelled"}[t.state])
            item.setToolTip(1, t.error)
            self.rows.itemWidget(item, 3).setEnabled(not t.finished())
        if not all(t.finished() for t in self.transfers):
            return
        self.timer.stop()
        bad = [t for t in self.transfers if t.state != "done"]
        if not bad:
            return self._finish()
        got = len(self.transfers) - len(bad)
        self.summary.setText(f"{got} of {len(self.transfers)} copied; {len(bad)} failed or stopped (see State). Retry Failed copies them "
                             "again; Load Tracks keeps those that came.")
        self.retry.show()
        self.buttons.button(QDialogButtonBox.Cancel).setText("Close")
        genome_ok = self.genome_transfer is None or self.genome_transfer.state == "done"
        self.buttons.button(QDialogButtonBox.Ok).setEnabled(genome_ok and got > 0)

    def _finish(self):
        self.copied = track_hub.copied([t for t in self.transfers if t.track is not None])
        if self.genome_transfer is not None:
            self.genome_got = self.genome_transfer.result
        super().accept()

    def _pick_folder(self):
        folder = QFileDialog.getExistingDirectory(self, "Folder of a track hub (holding hub.txt)")
        if folder:
            self.source.setText(folder)
            self.read()

    def _pick_target(self):
        folder = QFileDialog.getExistingDirectory(self, "Copy the hub's files to", self.folder.text())
        if folder:
            self.folder.setText(folder)

    def read(self):
        text = self.source.text().strip()
        if not text:
            return
        try:
            self.hub = run_busy(self, "Reading the hub…", lambda: track_hub.open_hub(text))
        except UserError as e:
            self.hub = None
            self.summary.setText(str(e))
            self.tree.clear()
            self._count()
            return
        ed = self.editor
        self.folder.setText("")
        self._folder_state()
        self.genome.blockSignals(True)
        self.genome.clear()
        self.genome.addItems(list(self.hub.genomes))
        pick = track_hub.pick_genome(self.hub, [ed.fasta_record, ed.doc.name])
        if pick:
            self.genome.setCurrentText(pick)
        self.genome.blockSignals(False)
        self._list(self.genome.currentText())

    def _folder_state(self):
        """The folder is asked for when files are copied: a hub on the internet (copied beside the open file, in a folder named
        after the hub), or the genome too (its sequence and the tracks are kept there, with the project). A hub on this
        computer whose tracks go on the open file is read where it is."""
        if self.hub is None:
            return
        ed = self.editor
        need = track_hub.is_url(self.hub.source) or self.with_genome.isChecked()
        here = Path(ed.path).parent if not getattr(ed, "untitled", False) and ed.path else Path.home()
        if need and not self.folder.text().strip():
            self.folder.setText(str(here / self.hub.name))
        elif not need:
            self.folder.setText("")
        self.folder.setEnabled(need)
        self.folder.setPlaceholderText("" if need else "(a hub on this computer: its files are read where they are)")
        self._count()

    def _list(self, genome):
        self.tree.clear()
        self.tracks = []
        if self.hub is None or not genome:
            self._count()
            return
        try:
            self.tracks = run_busy(self, "Reading its tracks…", lambda: track_hub.hub_tracks(self.hub, genome))
        except UserError as e:
            self.summary.setText(str(e))
            return
        self.tree.blockSignals(True)
        groups = {}
        for i, t in enumerate(self.tracks):
            parent = groups.get(t.group)
            if parent is None:
                parent = groups[t.group] = QTreeWidgetItem([t.group])
                parent.setFlags(Qt.ItemIsEnabled | Qt.ItemIsUserCheckable | Qt.ItemIsAutoTristate)
                self.tree.addTopLevelItem(parent)
            item = QTreeWidgetItem([t.label, t.kind, "hidden" if t.hidden else "shown"])
            item.setData(0, Qt.UserRole, i)
            item.setToolTip(0, t.url)
            if t.why:   # a type openDNA does not read: listed, never ticked
                item.setFlags(Qt.NoItemFlags)
                item.setText(2, "not read")
                item.setToolTip(2, t.why)
            else:
                item.setFlags(Qt.ItemIsEnabled | Qt.ItemIsUserCheckable)
                item.setCheckState(0, Qt.Checked)
            parent.addChild(item)
        self.tree.expandAll()
        self.tree.resizeColumnToContents(0)
        self.tree.blockSignals(False)
        self._count()

    def chosen(self):
        """The tracks ticked (those openDNA reads)."""
        out = []
        for i in range(self.tree.topLevelItemCount()):
            top = self.tree.topLevelItem(i)
            for j in range(top.childCount()):
                item = top.child(j)
                if item.flags() & Qt.ItemIsUserCheckable and item.checkState(0) == Qt.Checked:
                    out.append(self.tracks[item.data(0, Qt.UserRole)])
        return out

    def target(self):
        """The folder to copy the files to, None for a hub on this computer (read where it is)."""
        return self.folder.text().strip() or None if self.folder.isEnabled() else None

    def _count(self):
        chosen = self.chosen()
        ok = self.buttons.button(QDialogButtonBox.Ok)
        ok.setEnabled(bool(chosen) and (not self.folder.isEnabled() or bool(self.folder.text().strip())))
        if self.hub is None:
            if not self.summary.text():
                self.summary.setText("Type the address of a hub's hub.txt (UCSC track hub), or pick a folder holding one, then Read.")
            return
        unread = [t for t in self.tracks if t.why]
        self.summary.setText(
            f"{plural(len(self.tracks), 'track')} with data, {len(chosen)} ticked; {sum(t.hidden for t in chosen)} of them start hidden "
            "(as the hub has them): tick them in the Tracks panel to show them."
            + (f" {plural(len(unread), 'track')} of a type openDNA does not read: listed, not loaded." if unread else ""))
