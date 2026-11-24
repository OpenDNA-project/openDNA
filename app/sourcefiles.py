"""Files a tool can take DNA from: the recently opened list and a reader for any sequence file."""
from pathlib import Path

from PySide6.QtCore import QSettings
from PySide6.QtWidgets import QInputDialog, QMessageBox

from opendna_core import Document, Fasta

from primer_db import data_dir

KEY, PIN_KEY, MAX_RECENT = "files/recent", "files/pinned", 12
FASTA_EXT = {".fa", ".fasta", ".fna", ".fas", ".ffn", ".frn", ".mfa"}
OPEN_FILTER = "Sequences (*.gb *.gbk *.gbff *.genbank *.ape *.embl *.dna *.fa *.fasta *.fna *.fas);;All files (*)"


def recent_store():
    """recent.ini: the files and background genomes used lately."""
    return QSettings(str(data_dir() / "recent.ini"), QSettings.IniFormat)


def _existing(key):
    value = recent_store().value(key, [])
    paths = [value] if isinstance(value, str) else list(value or [])
    return [p for p in paths if Path(p).exists()]


def recent_files():
    """Paths opened lately, newest first, that still exist."""
    return _existing(KEY)


def add_recent(path):
    path = str(Path(path).resolve())
    paths = [path] + [p for p in recent_files() if p != path]
    recent_store().setValue(KEY, paths[:MAX_RECENT])


def pinned_files():
    """Paths pinned on the start page (kept on top, never pushed out by newer files), that still exist."""
    return _existing(PIN_KEY)


def set_pinned(path, pinned):
    paths = [p for p in pinned_files() if p != path]
    recent_store().setValue(PIN_KEY, paths + [path] if pinned else paths)


def remove_recent(path):
    """Off the start page: out of the recent list and unpinned (the file itself is not touched)."""
    set_pinned(path, False)
    recent_store().setValue(KEY, [p for p in recent_files() if p != path])


def read_document(parent, path):
    """The Document of a GenBank / .dna / FASTA file (a multi-sequence FASTA asks which one),
    or None after telling the user why not."""
    path = Path(path)
    try:
        if path.suffix.lower() in FASTA_EXT:
            names = [n for n, _ in Fasta(str(path)).sequences()]
            if not names:
                raise OSError("no sequence in this FASTA file")
            name = names[0]
            if len(names) > 1:
                name, ok = QInputDialog.getItem(parent, "Which sequence?", f"{path.name} holds {len(names)} sequences:",
                                                names, 0, False)
                if not ok:
                    return None
            return Document.from_fasta(str(path), name)
        return Document(str(path))
    except (OSError, ValueError) as e:
        QMessageBox.warning(parent, "Cannot Open File", f"{path}\n\n{e}")
        return None
