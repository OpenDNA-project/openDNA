"""The primer database: every primer you create or import, independent of any sequence
file. Stored as a tab-separated file in the app's data folder
(%APPDATA%/openDNA/primers.tsv on Windows; OPENDNA_DATA_DIR overrides
it, e.g. for tests). Export/import as TSV, CSV or FASTA.

Names are unique: the same name with a different sequence is refused, never overwritten.
"""
import csv
import datetime
import os
import shutil
from pathlib import Path

from PySide6.QtCore import QObject, QSettings, QStandardPaths, Signal

FIELDS = ["name", "sequence", "note", "added"]
MAX_PRIMER, MAX_NAME, MAX_NOTE = 200, 60, 200   # a longer "primer" is a sequence (and would stall the checks); a name this long is a mistake
PRIMER_BASES = set("ACGTRYSWKMBDHVNacgtryswkmbdhvn")   # DNA and the IUPAC ambiguity codes, either case (a lowercase 5' tail is fine)


def primer_problem(name, seq, note=""):
    """What is wrong with a primer name / sequence / note, in a sentence, or "" when it can be saved. Every way into the database
    (the primer dialog, Design, Assembly, an imported list, the primers of an opened file) goes through it."""
    if not name:
        return "Give the primer a name."
    if any(c.isspace() and c != " " for c in name) or not name.isprintable():
        return "The name has a tab, line break or other control character: use plain letters, digits and spaces."
    if len(name) > MAX_NAME:
        return f"The name is too long: at most {MAX_NAME} characters ({len(name)} given)."
    if not seq:
        return "There is no sequence."
    bad = sorted({c for c in seq if c not in PRIMER_BASES})
    if bad:
        shown = " ".join(repr(c)[1:-1] if not c.isprintable() else c for c in bad[:8])
        return f"The sequence has characters that are not DNA ({shown}). Use A C G T (or IUPAC codes such as N, R, Y)."
    if len(seq) < 8:
        return f"A primer needs at least 8 bases ({len(seq)} given)."
    if len(seq) > MAX_PRIMER:
        return f"A primer is at most {MAX_PRIMER} bases ({len(seq):,} given)."
    if note and (len(note) > MAX_NOTE or any(c.isspace() and c != " " for c in note) or not note.isprintable()):
        return f"The note must be one line of at most {MAX_NOTE} characters, with no tab or control character."
    return ""


class _Notifier(QObject):
    changed = Signal()  # the database (or a primer option) changed: re-match primers


notifier = _Notifier()


def binding_mismatches():
    """Mismatches allowed in a primer's last 15 bases (never the 3 at the 3' end) for a
    site to be shown: 0 = exact sites only; 1+ also shows where it almost binds (default 2)."""
    return int(settings().value("primers/binding_mismatches", 2))


def settings():
    """The primer and PCR settings: an .ini file in the data folder like every other setting, not the Windows registry."""
    return QSettings(str(data_dir() / "primers.ini"), QSettings.IniFormat)


def data_dir():
    folder = os.environ.get("OPENDNA_DATA_DIR") or QStandardPaths.writableLocation(QStandardPaths.AppDataLocation)
    path = Path(folder)
    path.mkdir(parents=True, exist_ok=True)
    return path


OLD_DATA_NAME = "OpenGenEdit"   # the program's name before openDNA


def copy_old_data(new=None):
    """Copies the data folder of the old name (%APPDATA%/OpenGenEdit: settings, primers, recent files) to the new one, once: only
    while the new folder does not exist. The old folder stays as it was. Returns True when it copied."""
    if new is None:
        if os.environ.get("OPENDNA_DATA_DIR"):
            return False
        new = Path(QStandardPaths.writableLocation(QStandardPaths.AppDataLocation))
    old = new.with_name(OLD_DATA_NAME)
    if new.exists() or not old.is_dir():
        return False
    try:
        shutil.copytree(old, new)
    except OSError:   # a half copy: the program starts with what arrived, and the old folder still holds everything
        return False
    return True


class PrimerDB:
    def __init__(self, path):
        self.path = Path(path)
        self.rows = []  # dicts with FIELDS
        if self.path.exists():
            with open(self.path, encoding="utf-8", newline="") as fh:
                self.rows = [{k: r.get(k) or "" for k in FIELDS} for r in csv.DictReader(fh, delimiter="\t")]

    def save(self):
        tmp = self.path.with_suffix(".tmp")  # write then rename: never a half-written DB
        with open(tmp, "w", encoding="utf-8", newline="") as fh:
            writer = csv.DictWriter(fh, FIELDS, delimiter="\t")
            writer.writeheader()
            writer.writerows(self.rows)
        os.replace(tmp, self.path)
        notifier.changed.emit()

    def find(self, name):
        return next((r for r in self.rows if r["name"] == name), None)

    def add(self, name, sequence, note=""):
        """Returns "added", "exists" (same name and sequence), "conflict" (name taken) or "invalid" (see primer_problem: nothing is
        stored)."""
        result = self._add(name, sequence, note)
        if result == "added":
            self.save()
        return result

    def _add(self, name, sequence, note=""):
        name, sequence = name.strip(), "".join(sequence.split())
        if primer_problem(name, sequence, note):
            return "invalid"
        old = self.find(name)
        if old:
            return "exists" if old["sequence"].upper() == sequence.upper() else "conflict"
        self.rows.append({"name": name, "sequence": sequence, "note": note,
                          "added": datetime.date.today().isoformat()})
        return "added"

    def add_many(self, primers):
        """[(name, sequence)] -> {"added": n, "exists": n, "conflict": [names], "invalid": [(name, why)]}; one save."""
        report = {"added": 0, "exists": 0, "conflict": [], "invalid": []}
        for name, sequence in primers:
            result = self._add(name, sequence)
            if result == "conflict":
                report["conflict"].append(name)
            elif result == "invalid":
                report["invalid"].append((name, primer_problem(name.strip(), "".join(sequence.split()))))
            else:
                report[result] += 1
        if report["added"]:
            self.save()
        return report

    def update(self, old_name, name, sequence, note):
        """Edit a primer. Returns False if `name` belongs to another primer; raises ValueError if the new values are not a primer."""
        name, sequence = name.strip(), "".join(sequence.split())
        problem = primer_problem(name, sequence, note)
        if problem:
            raise ValueError(problem)
        if name != old_name and self.find(name):
            return False
        row = self.find(old_name)
        row.update(name=name, sequence=sequence, note=note)
        self.save()
        return True

    def remove(self, names):
        self.rows = [r for r in self.rows if r["name"] not in set(names)]
        self.save()

    def restore(self, indexed_rows):
        """Put back rows [(index, row)] (as they were, with their date) where they were. Returns the names that could not
        come back because the name is taken now. One save."""
        skipped = []
        for index, row in sorted(indexed_rows, key=lambda x: x[0]):
            if self.find(row["name"]):
                skipped.append(row["name"])
            else:
                self.rows.insert(min(index, len(self.rows)), dict(row))
        if len(skipped) < len(indexed_rows):
            self.save()
        return skipped

    def export(self, path):
        """FASTA for .fa/.fasta, CSV for .csv, otherwise TSV (same columns as the DB)."""
        path = Path(path)
        if path.suffix.lower() in (".fa", ".fasta"):
            path.write_text("".join(f">{r['name']}\n{r['sequence']}\n" for r in self.rows), encoding="utf-8")
            return
        with open(path, "w", encoding="utf-8", newline="") as fh:
            writer = csv.DictWriter(fh, FIELDS, delimiter="," if path.suffix.lower() == ".csv" else "\t")
            writer.writeheader()
            writer.writerows(self.rows)


_db = None


def db():
    """The shared database (opened on first use)."""
    global _db
    if _db is None:
        _db = PrimerDB(data_dir() / "primers.tsv")
    return _db


# ---- changes that go in a document's History, and so can be undone with it ------------------------------------------------------
# The database is shared by every open document, but a change is made from one: it is a step in that document's history (see
# Editor.record_external). Undoing it does the opposite to the database, unless the primer was changed again since: then it is
# left as it is and the status bar says so. Without an editor (`None`) the change is made and not recorded.


def _names(names, limit=3):
    return ", ".join(names[:limit]) + (f" and {len(names) - limit} more" if len(names) > limit else "")


def _noun(n):
    return "primer" if n == 1 else "primers"


def _take_out(rows):
    """Remove `rows` from the database if they are still as they were. Returns a message for the ones that are not."""
    changed = [r["name"] for r in rows if db().find(r["name"]) != r]
    still = [r["name"] for r in rows if r["name"] not in changed]
    if still:
        db().remove(still)
    return f"Not removed, changed since: {_names(changed)}." if changed else None


def _put_back(indexed):
    skipped = db().restore(indexed)
    return f"Not restored, the name is used now: {_names(skipped)}." if skipped else None


def _in_file(editor, rows, present):
    """Mark `rows` as primers of the document `editor` shows (or not): what Save writes into the file."""
    if editor is not None and rows:
        editor.set_in_file([(r["name"], r["sequence"]) for r in rows], present)


def add_primer(editor, name, sequence, note=""):
    """db().add, as a history step: "added", "exists" or "conflict". A primer made from this document's sequence also
    becomes one of the file's own primers (drawn as a regular primer, written on Save)."""
    result = db().add(name, sequence, note)
    if result == "added" and editor is not None:
        row = dict(db().find(name.strip()))
        index = db().rows.index(row)
        _in_file(editor, [row], True)

        def undo():
            _in_file(editor, [row], False)
            return _take_out([row])

        def redo():
            _in_file(editor, [row], True)
            return _put_back([(index, row)])
        editor.record_external(f"Added primer {row['name']}", undo, redo)
    return result


def add_primers(editor, primers, what="Added", in_file=False):
    """db().add_many, as one history step. Returns its report. `in_file`: the primers belong to this document too (designed
    for its sequence), not only to the database."""
    before = {r["name"] for r in db().rows}
    report = db().add_many(primers)
    if report["added"] and editor is not None:
        indexed = [(i, dict(r)) for i, r in enumerate(db().rows) if r["name"] not in before]
        rows = [r for _, r in indexed]
        if in_file:
            _in_file(editor, rows, True)

        def undo():
            if in_file:
                _in_file(editor, rows, False)
            return _take_out(rows)

        def redo():
            if in_file:
                _in_file(editor, rows, True)
            return _put_back(indexed)
        editor.record_external(f"{what} {len(rows)} {_noun(len(rows))}: {_names([r['name'] for r in rows])}", undo, redo)
    return report


def update_primer(editor, old_name, name, sequence, note):
    """db().update, as a history step. False if `name` belongs to another primer."""
    before = dict(db().find(old_name))
    if not db().update(old_name, name, sequence, note):
        return False
    after = dict(db().find(name.strip()))
    if editor is not None and after != before:
        def switch(source, target):
            if db().find(source["name"]) != source:
                return f"Primer {source['name']} was changed since: left as it is."
            if not db().update(source["name"], target["name"], target["sequence"], target["note"]):
                return f"The name {target['name']} is used now: primer {source['name']} left as it is."
            editor.rename_in_file(source, target)
            return None
        editor.rename_in_file(before, after)
        what = f"Renamed primer {before['name']} to {after['name']}" if before["name"] != after["name"] else f"Edited primer {after['name']}"
        editor.record_external(what, lambda: switch(after, before), lambda: switch(before, after))
    return True


def remove_primers(editor, names):
    """db().remove, as one history step."""
    indexed = [(i, dict(r)) for i, r in enumerate(db().rows) if r["name"] in set(names)]
    db().remove(names)
    if indexed and editor is not None:
        rows = [r for _, r in indexed]
        stored = [r for r in rows if editor.in_file(r["name"], r["sequence"])]   # deleted from the file's own primers too
        _in_file(editor, stored, False)

        def undo():
            _in_file(editor, stored, True)
            return _put_back(indexed)

        def redo():
            _in_file(editor, stored, False)
            return _take_out(rows)
        editor.record_external(f"Deleted {_noun(len(rows))} {_names([r['name'] for r in rows])}", undo, redo)
