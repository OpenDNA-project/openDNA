"""Folder Search (the "Folder" scope of Local Align): "which of my plasmids carries this sequence?" Aligns the query
against every sequence record of every GenBank / .dna / FASTA file in a folder (the Rust `align`, both strands,
mismatches and gaps allowed) and lists the hits; a click opens that file at that position.

No index: `align` takes about 0.02 s per plasmid, so a freezer of 500 files takes seconds. It runs in a background
thread with progress and a Cancel button, and a file that cannot be read is skipped and counted, never fatal."""
import threading
from pathlib import Path
from types import SimpleNamespace

from PySide6.QtCore import QSettings, QTimer
from PySide6.QtWidgets import (
    QAbstractItemView, QApplication, QFileDialog, QHBoxLayout, QLabel, QLineEdit, QProgressBar, QPushButton,
    QWidget,
)

from opendna_core import Document, Fasta, align

import multigb
import tool_settings
from local_align import TOOL, strand_text, time_limit_note, warn_time_limit
from primer_db import data_dir
from reference import is_fasta_file
from sourcefiles import FASTA_EXT
from tools import empty_hint
from data_table import DataTable
from uikit import fit_columns, plural

SEQUENCE_EXT = multigb.GENBANK_EXT | FASTA_EXT | {".dna"}
PER_FILE = 5          # hits kept per file: one repetitive genome must not flood the list
FOLDER_KEY = "folder_search/folder"


def sequence_files(folder):
    """Every sequence file in `folder` and its sub-folders (freezer boxes are often one folder each), sorted by path."""
    return sorted(p for p in Path(folder).rglob("*") if p.is_file() and p.suffix.lower() in SEQUENCE_EXT)


def records_of(path):
    """(record name, bases, circular) for each record of a file; raises OSError / ValueError when it cannot be read."""
    path = Path(path)
    if is_fasta_file(path):
        fasta = Fasta(str(path))
        for name, length in fasta.sequences():
            yield name, fasta.fetch(name, 0, int(length)), False   # (FASTA has no topology: searched as linear)
    elif multigb.is_multi_genbank(path):
        for name, _length, _circular in multigb.records(path):
            doc = multigb.read_record(path, name)
            yield name, doc.sequence(), doc.circular
    else:
        doc = Document(str(path))
        yield doc.name or path.stem, doc.sequence(), doc.circular


def search_folder(folder, query, progress=None, cancelled=lambda: False):
    """(hits, skipped files, files searched, notes). A hit has path, record, start, end (0-based, half open), strand, identity and
    length (bases of the query that matched). Sorted by identity, then length, best first. The notes name each file whose hits were
    cut (more than PER_FILE) or where `align` reached a search limit. `progress(done, total)` is called after each file;
    `cancelled()` is asked before each one, and the hits so far are returned when it says yes."""
    files = sequence_files(folder)
    hits, skipped, notes = [], 0, []
    for i, path in enumerate(files):
        if cancelled():
            break
        try:
            found = []
            for name, text, circular in records_of(path):
                pieces, limits = align(text, query, circular, PER_FILE + 1, tool_settings.value(TOOL, "time_limit"), tool_settings.cores())   # (one more than kept: tells whether the file has more)
                found += [SimpleNamespace(path=path, record=name, start=h[0], end=h[1], strand=h[2], identity=h[4],
                                          length=h[6] - h[5]) for h in pieces]
                notes += [f"{path.name}: record {name}: {n}" for n in limits]
        except (OSError, ValueError):
            skipped += 1
            continue
        finally:
            if progress:
                progress(i + 1, len(files))
        found.sort(key=lambda h: (-h.identity, -h.length))
        if len(found) > PER_FILE:
            notes.append(f"{path.name}: more than {PER_FILE} hits, only the {PER_FILE} best are listed (open the file and use Local Align for all).")
        hits += found[:PER_FILE]
    hits.sort(key=lambda h: (-h.identity, -h.length, str(h.path), h.start))
    return hits, skipped, len(files), notes


def notes_summary(notes):
    """The status-line words for the notes of `search_folder`: how many files they concern and the first one."""
    if not notes:
        return ""
    files = {n.split(": ", 1)[0] for n in notes}   # (each note starts with its file name)
    return f"; {plural(len(files), 'file')} reached a search limit (Copy lists them all), e.g. {notes[0]}"


def hits_text(hits, notes=()):
    """The hits as tab-separated lines with a header, to paste in a lab notebook or a spreadsheet; then the notes, one per line."""
    lines = ["File\tRecord\tPosition\tStrand\tIdentity\tLength (bp)"]
    lines += [row_text(h, "\t") for h in hits]
    lines += [f"Note: {n}" for n in notes]
    return "\n".join(lines)


def row_text(h, sep):
    return sep.join([Path(h.path).name, h.record, f"{h.start + 1}-{h.end}", strand_text(h.strand), f"{h.identity:.1%}", str(h.length)])


class FolderSearchPanel(QWidget):
    """The folder row of Local Align: folder box and Browse. The owner places `self.table` (the hits), `self.cancel_button` and
    `self.bar` (shown while a search runs), and runs `start()` from its own Align button; `on_busy(True / False)` says when a
    search starts and ends. `get_query()` returns the bases to search for or raises ValueError (its text is shown by `say`);
    `win` opens a hit's file."""
    COLUMNS = ["File", "Record", "Position", "Strand", "Identity", "Length"]

    def __init__(self, win, get_query, say, on_start=lambda: None, on_busy=lambda busy: None):
        super().__init__()
        self.win, self.get_query, self.say, self.on_start, self.on_busy, self.hits, self.notes = win, get_query, say, on_start, on_busy, [], []
        self.thread, self.state = None, None
        self.folder = QLineEdit(str(QSettings(str(data_dir() / "window.ini"), QSettings.IniFormat).value(FOLDER_KEY, "")),
                                placeholderText="Folder of sequence files")
        self.bar = QProgressBar(maximumWidth=160, visible=False)
        self.table = DataTable(self.COLUMNS, filtered=["Strand"])
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        fit_columns(self.table, stretch=0)
        self.table.setWordWrap(False)
        empty_hint(self.table, "Pick a folder and press Align")
        self.table.cellDoubleClicked.connect(lambda row, _c: self.open_hit(self.table.data_row(row)))
        self.table.itemActivated.connect(lambda item: self.open_hit(self.table.data_row(item.row())))
        self.cancel_button = QPushButton("Cancel", clicked=self.cancel, visible=False)
        row = QHBoxLayout(self)
        row.setContentsMargins(0, 0, 0, 0)
        row.addWidget(QLabel("Folder"))
        row.addWidget(self.folder, 1)
        row.addWidget(QPushButton("Browse…", clicked=self.browse))
        self.timer = QTimer(self, interval=50, timeout=self._poll)

    def browse(self):
        path = QFileDialog.getExistingDirectory(self, "Folder of sequence files", self.folder.text())
        if path:
            self.folder.setText(path)

    def start(self):
        folder = self.folder.text().strip()
        try:
            query = self.get_query()
        except ValueError as e:
            self.say(str(e))
            return
        if not Path(folder).is_dir():
            self.say("Pick an existing folder")
            return
        QSettings(str(data_dir() / "window.ini"), QSettings.IniFormat).setValue(FOLDER_KEY, folder)
        self.state = SimpleNamespace(done=0, total=0, cancel=False, result=None, error=None)
        s = self.state

        def work():
            try:
                s.result = search_folder(folder, query, lambda d, t: (setattr(s, "done", d), setattr(s, "total", t)),
                                         lambda: s.cancel)
            except BaseException as e:   # shown in the status line by _poll
                s.error = e

        self.thread = threading.Thread(target=work, daemon=True)
        self.thread.start()
        self.on_busy(True)
        self.cancel_button.setVisible(True)
        self.bar.setVisible(True)
        self.say("Searching…")
        self.on_start()
        self.timer.start()

    def cancel(self):
        if self.state:
            self.state.cancel = True
            self.say("Cancelling…")

    def _poll(self):
        s = self.state
        self.bar.setMaximum(max(s.total, 1))
        self.bar.setValue(s.done)
        if self.thread.is_alive():
            return
        self.timer.stop()
        self.on_busy(False)
        self.cancel_button.setVisible(False)
        self.bar.setVisible(False)
        if s.error:
            self.say(f"Search failed: {s.error}")
            return
        self.hits, skipped, total, self.notes = s.result
        self.fill()
        note = f"; {plural(skipped, 'file')} could not be read and skipped" if skipped else ""
        stopped = f" (cancelled after {s.done:,} of {total:,} files)" if s.cancel else f" in {plural(total, 'file')}"
        self.say((plural(len(self.hits), "hit") if self.hits else "No hit") + stopped + note + notes_summary(self.notes))
        if time_limit_note(self.notes):
            warn_time_limit(self, time_limit_note(self.notes), len({n.split(": ", 1)[0] for n in self.notes if time_limit_note([n])}))

    def stop(self):
        """The dialog is closing: stop the search too."""
        if self.thread and self.thread.is_alive():
            self.state.cancel = True
            self.thread.join()

    def fill(self):
        with self.table.refilling(len(self.hits)):
            for i, h in enumerate(self.hits):
                values = [Path(h.path).name, h.record, f"{h.start + 1:,}–{h.end:,}", strand_text(h.strand), f"{h.identity:.1%}", f"{h.length:,}"]
                for c, v in enumerate(values):
                    self.table.put(i, c, v, right=c in (4, 5), tip=str(h.path) if c == 0 else None)

    def copy_results(self):
        QApplication.clipboard().setText(hits_text(self.hits, self.notes))
        self.say(f"Copied {plural(len(self.hits), 'hit')}")

    def open_hit(self, row):
        """Open the hit's file (its record, for a multi-record file) and select the matching bases."""
        if not 0 <= row < len(self.hits):
            return None
        h = self.hits[row]
        win = self.win.open_path(str(h.path))   # (an empty window is used, else a new one)
        if win is not None and win.editor is not None:
            if win.editor.fasta_record and win.editor.fasta_record != h.record:
                win.load(str(h.path), h.record)
            win.editor.select(h.start, h.end)
            win.editor.reveal.emit(h.start)
        return win
