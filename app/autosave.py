"""Recovery copies: what the program writes while a file has unsaved changes, so a crash does not lose them. They live in the
program's data folder (autosave/), not beside the user's files; each has a small `.src` file next to it saying which file it is for."""
import hashlib
import os
import shutil
import time
from pathlib import Path

from primer_db import data_dir

SWEEP_DAYS = 30   # a copy whose file is gone is deleted after this long


def folder():
    path = data_dir() / "autosave"
    path.mkdir(exist_ok=True)
    return path


def _name(path, record):
    return Path(path).name + (f".{record}" if record else "")


def autosave_path(path, record=None):
    """The recovery copy of a file: plasmid.dna -> autosave/plasmid.dna-3fa9c1d2.autosave.gb. The 8 hex digits come from the file's
    full path, so two files of one name in different folders never share a copy. FASTA: one per sequence (genome.fa.chr2-...)."""
    key = hashlib.sha1(os.path.normcase(os.path.realpath(path)).encode("utf-8")).hexdigest()[:8]
    return folder() / f"{_name(path, record)}-{key}.autosave.gb"


def old_path(path, record=None):
    """Where earlier versions wrote it: beside the file (plasmid.dna.autosave.gb)."""
    return Path(path).with_name(_name(path, record) + ".autosave.gb")


def remember_source(target, path):
    """Note which file the copy `target` is for (sweep needs it)."""
    target.with_suffix(".src").write_text(str(path), encoding="utf-8")


def drop(target):
    target.unlink(missing_ok=True)
    target.with_suffix(".src").unlink(missing_ok=True)


def drop_all(path, record=None):
    """Every recovery copy of this file (record), wherever a version of the program put it."""
    drop(autosave_path(path, record))
    old_path(path, record).unlink(missing_ok=True)


def adopt_old(path, record=None):
    """A copy an earlier version left beside the file moves into the data folder, to be offered like any other (once: declining or
    saving deletes it). The newer one is kept if both exist."""
    old = old_path(path, record)
    if not old.exists():
        return
    new = autosave_path(path, record)
    try:
        if new.exists() and new.stat().st_mtime >= old.stat().st_mtime:
            old.unlink()
        else:
            shutil.move(str(old), str(new))
            remember_source(new, path)
    except OSError:
        pass   # (a read-only folder: the copy stays where it is and is simply not offered)


def sweep():
    """At start: delete the copies of files that no longer exist and are older than SWEEP_DAYS. A copy of a file that is still
    there is never deleted: it waits until the file is opened."""
    limit = time.time() - SWEEP_DAYS * 86400
    for copy in folder().glob("*.autosave.gb"):
        try:
            source = copy.with_suffix(".src")
            gone = not source.exists() or not Path(source.read_text(encoding="utf-8")).exists()
            if gone and copy.stat().st_mtime < limit:
                drop(copy)
        except OSError:
            pass
