"""The background genome box of the toolbar: the genome the primers (and the Spacer Finder's guides) are checked against, chosen
among the last ones used (RECENT of them), or "Add genome file…" (the box's last entry) for another. It is the same setting as Primers → Options → Genome file, so
choosing one here, there, or downloading one from NCBI (Import → DNA → Background Genome) is the same thing; every window's box
follows (primer_db.notifier)."""
from pathlib import Path

from PySide6.QtCore import QEvent, Qt, Signal
from PySide6.QtWidgets import QComboBox, QFileDialog, QHBoxLayout, QWidget

import icons
import primer_db

RECENT = 3
KEY = "primers/recent_genomes"
OFF = "No background genome"
ADD = "Add genome file…"
FILTER = "Sequence files (*.gb *.gbk *.gbff *.genbank *.dna *.fa *.fasta *.fna *.fas);;All files (*)"


def recent():
    """The last background genomes used, newest first (only files that still exist)."""
    raw = primer_db.settings().value(KEY, [])
    paths = [raw] if isinstance(raw, str) and raw else list(raw or [])
    return [p for p in paths if p and Path(p).is_file()][:RECENT]


def remember(path):
    """Put `path` first in the list of the last genomes used."""
    if not path:
        return
    path = str(path)
    primer_db.settings().setValue(KEY, [path] + [p for p in recent() if Path(p) != Path(path)][:RECENT - 1])


def current():
    """The genome in use, or None when the genome check is off or there is no genome at all."""
    from primer_panel import genome_path
    if primer_db.settings().value("primers/genome_check", "true") != "true":
        return None
    return genome_path() or None


def use(path):
    """Check primers against `path` from now on (None: against no genome)."""
    s = primer_db.settings()
    if path:
        s.setValue("primers/genome", str(path))
        s.setValue("primers/genome_check", "true")
        remember(path)
    else:
        s.setValue("primers/genome_check", "false")
    primer_db.notifier.changed.emit()   # every window re-checks its primers (and its box shows the choice)


class ReferencePicker(QWidget):
    """A compact box in the toolbar: a genome icon and the genome's file name (its whole path and what it is for in the tooltip);
    its list holds the last genomes used, "No background genome", then "Add genome file…"."""
    changed = Signal()   # the genome was changed from this box

    def __init__(self, parent=None):
        super().__init__(parent)
        self.box = QComboBox(minimumContentsLength=9, sizeAdjustPolicy=QComboBox.AdjustToMinimumContentsLengthWithIcon)
        self.box.setMaximumWidth(self.box.fontMetrics().horizontalAdvance("M" * 9) + 40)
        self.box.activated.connect(self._chosen)
        row = QHBoxLayout(self)
        row.setContentsMargins(0, 0, 0, 0)
        row.addWidget(self.box)
        primer_db.notifier.changed.connect(self.refresh)   # (a QObject slot: the link goes with the box)
        self.refresh()

    def refresh(self):
        """The last genomes used (the one in use first if it is not among them, e.g. the default one), "no genome", "add"."""
        now = current()
        paths = recent()
        if now and not any(Path(p) == Path(now) for p in paths):
            paths = [now] + paths[:RECENT - 1]
        icon = icons.genome_icon(self.palette())
        self.box.blockSignals(True)
        self.box.clear()
        for p in paths:
            self.box.addItem(icon, Path(p).name, p)
            self.box.setItemData(self.box.count() - 1, p, Qt.ToolTipRole)   # the whole path
        self.box.addItem(OFF, None)
        self.box.insertSeparator(self.box.count())
        self.box.addItem(ADD, ADD)
        match = next((i for i, p in enumerate(paths) if now and Path(p) == Path(now)), None)
        self.box.setCurrentIndex(match if match is not None else len(paths))
        self.box.blockSignals(False)
        shown = self.box.currentData()
        self.box.setToolTip(f"Background genome: {shown}\nPrimers and guides are checked against it for other binding sites."
                            if shown else "No background genome: primers and guides are checked against this sequence only.")

    def _chosen(self, i):
        path = self.box.itemData(i)
        if path == ADD:
            if not self.browse():
                self.refresh()   # cancelled: the box shows the genome in use again
            return
        if path != current():
            use(path)
            self.changed.emit()

    def browse(self):
        path, _ = QFileDialog.getOpenFileName(self, "Background genome", "", FILTER)
        if path:
            use(path)
            self.changed.emit()
        return path

    def changeEvent(self, event):
        if event.type() == QEvent.PaletteChange:   # the icon follows the theme
            self.refresh()
        super().changeEvent(event)
