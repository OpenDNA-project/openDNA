"""Import > Features and Tracks > Track Hub: the dialog that reads a UCSC track hub (a URL or a local folder) and lists its tracks
in the hub's own groups, all ticked (the hidden and grouped ones too: those the hub hides at first load hidden), with where the
files will be copied. main.import_hub copies and loads them; the hub itself is read by track_hub."""
from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QComboBox, QDialog, QDialogButtonBox, QFileDialog, QFormLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton, QTreeWidget,
    QTreeWidgetItem, QVBoxLayout,
)

import track_hub
from uikit import plural, run_busy
from user_error import UserError


class HubDialog(QDialog):
    def __init__(self, parent, editor):
        super().__init__(parent)
        self.setWindowTitle("Import a track hub")
        self.editor, self.hub, self.tracks = editor, None, []
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
        self.folder = QLineEdit()
        self.folder.textChanged.connect(lambda *_: self._count())
        browse = QPushButton("Browse…", clicked=self._pick_target)
        target = QHBoxLayout()
        target.addWidget(self.folder, 1)
        target.addWidget(browse)
        form = QFormLayout()
        form.addRow("Hub:", row)
        form.addRow("Genome:", self.genome)
        form.addRow("Copy the files to:", target)
        self.buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        self.buttons.button(QDialogButtonBox.Ok).setText("Load Tracks")
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)
        layout = QVBoxLayout(self)
        layout.addLayout(form)
        layout.addWidget(self.tree, 1)
        layout.addWidget(self.summary)
        layout.addWidget(self.buttons)
        self.resize(720, 560)
        self._count()

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
        # A hub on the internet is copied beside the open file (a folder named after the hub); a hub on this computer is read
        # where it is.
        remote = track_hub.is_url(self.hub.source)
        here = Path(ed.path).parent if not getattr(ed, "untitled", False) and ed.path else Path.home()
        self.folder.setText(str(here / self.hub.name) if remote else "")
        self.folder.setEnabled(remote)
        self.folder.setPlaceholderText("" if remote else "(a hub on this computer: its files are read where they are)")
        self.genome.blockSignals(True)
        self.genome.clear()
        self.genome.addItems(list(self.hub.genomes))
        pick = track_hub.pick_genome(self.hub, [ed.fasta_record, ed.doc.name])
        if pick:
            self.genome.setCurrentText(pick)
        self.genome.blockSignals(False)
        self._list(self.genome.currentText())

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
