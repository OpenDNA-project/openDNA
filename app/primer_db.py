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

from opendna_core import BINDING_SEED, EXACT_3P
import fileio
import tool_settings
from tool_settings import Param
from user_error import UserError

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


GENOME_FILTER = "Sequence files (*.gb *.gbk *.gbff *.genbank *.ape *.embl *.dna *.fa *.fasta *.fna *.fas);;All files (*)"

# The Primers page of Tools > Settings: declared here once (defaults, limits), read everywhere through tool_settings.
PCR = "PCR conditions"   # for every Tm (Primer3's defaults); no sodium or no primer at all makes the Tm meaningless, hence the minimums
tool_settings.register(PCR, [
    Param("na_mm", "Na⁺/K⁺ (mM)", 50.0, 1.0, 2000),
    Param("mg_mm", "Mg²⁺ (mM)", 1.5, 0.0, 100),
    Param("dntp_mm", "dNTPs (mM)", 0.6, 0.0, 100),
    Param("primer_nm", "Primer (nM)", 50.0, 1.0, 100000),
])


def _binding_check(v):
    """The binding settings limit each other: the exact 3′ bases are part of the compared ones, and the mismatches fit in the rest.
    A background genome that is switched on must be a file."""
    notes = []
    if v["exact_3p"] > v["seed"]:
        v["exact_3p"] = v["seed"]
        notes.append(f"exact 3′ bases lowered to {v['seed']} (the bases compared)")
    if v["mismatches"] > v["seed"] - v["exact_3p"]:
        v["mismatches"] = v["seed"] - v["exact_3p"]
        notes.append(f"mismatches lowered to {v['mismatches']} (the bases compared minus the exact 3′ ones)")
    if v["genome_check"] and v["genome"] and not Path(v["genome"]).is_file():
        raise UserError(f"There is no file '{v['genome']}'. Browse for the genome (FASTA or GenBank), or clear the box for none.")
    return "Primer binding: " + "; ".join(notes) + "." if notes else ""


BINDING = "Primer binding"   # where a primer is drawn, its "Binds here" check, the 3′ mis-priming check and the background genome
tool_settings.register(BINDING, [
    Param("seed", "3′ bases compared", BINDING_SEED, 8, 30, suffix=" nt",
          tip="A site is where the primer's last bases match the sequence (then the match extends toward the 5′ end)."),
    Param("exact_3p", "Of which exact at the 3′ end", EXACT_3P, 1, 30, suffix=" nt",
          tip="The polymerase starts at the 3′ end: these bases never tolerate a mismatch."),
    Param("mismatches", "Mismatches allowed", 2, 0, 29,
          tip="Mismatches allowed in the compared bases (never in the exact 3′ ones). 0 = exact sites only."),
    Param("min_tm", "Minimum Tm of a site", 40.0, 0.0, 100.0, suffix=" °C",
          tip="A site counts only if the primer's bases that anneal there melt at this temperature or above (PCR conditions,\n"
              "mismatches are not taken off). 0 = no minimum."),
    Param("best_only", "Show only the best binding site by default", False,
          tip="Fewest mismatches in the last 5 bases, then fewest mismatches, then the longest annealed part; ties are all shown.\n"
              "Also in the Display filter."),
    Param("homology_seed", "3′ check: bases", 12, 6, 20, suffix=" nt", tip="Bases at the primer's 3′ end compared for the mis-priming check"),
    Param("homology_mismatches", "3′ check: mismatches", 0, 0, 3),
    Param("genome_check", "Also check a genome", True),
    Param("genome", "Genome file", "", browse=GENOME_FILTER,
          tip="A FASTA or GenBank file; one with several sequences (chromosomes) is searched in full. Empty: no background genome."),
], check=_binding_check)
tool_settings.PAGE_OF.update({PCR: "Primers", BINDING: "Primers"})


def conditions():
    """The PCR conditions of every Tm (Settings > Primers), as keyword arguments of tm()."""
    return dict(tool_settings.values(PCR))


def tm_rule():
    """The minimum Tm of a site and the PCR conditions, as keyword arguments of find_primers / Document.binding_sites: the core
    keeps only the sites that melt that high, their mismatches counted. Every search for a primer's sites passes them."""
    return {"min_tm": tool_settings.value(BINDING, "min_tm"), **conditions()}


def binding():
    """(bases compared at the 3′ end, of which exact, mismatches allowed): the one rule for where a primer binds."""
    v = tool_settings.values(BINDING)
    return v["seed"], v["exact_3p"], v["mismatches"]


def data_dir():
    folder = os.environ.get("OPENDNA_DATA_DIR") or QStandardPaths.writableLocation(QStandardPaths.AppDataLocation)
    path = Path(folder)
    path.mkdir(parents=True, exist_ok=True)
    return path


def window_state(group=None):
    """window.ini in the data folder: what the interface remembers between sessions that is not a setting (window size, column
    widths, the search bar's and dialogs' last choices, the updater's and the sponsor banner's dates), one group per part."""
    state = QSettings(str(data_dir() / "window.ini"), QSettings.IniFormat)
    if group:
        state.beginGroup(group)
    return state


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
        fileio.replace(tmp, self.path)
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


ORDER_FORMATS = ["IDT plate (CSV)", "Name and sequence (CSV)"]
SCALES = ["25nm", "100nm", "250nm", "1um", "2um", "5um", "10um"]   # IDT's synthesis scales
PURIFICATIONS = ["STD", "PAGE", "HPLC", "IEHPLC", "RNASE", "PAGE_HPLC"]
PLATE_SIZE = 96


def well(i):
    """Well of the i-th oligo (from 0) on a 96-well plate, filled down the columns: A1, B1 … H1, A2 … (the order IDT expects)."""
    return f"{'ABCDEFGH'[i % 8]}{i // 8 % 12 + 1}"


def order_rows(rows):
    """[(name, sequence)] ready for a supplier's importer. A tab, line break or other control character in a name becomes a space
    (it would split the field); commas and quotes are left alone, the CSV writer quotes them. Degenerate bases are kept."""
    return [(" ".join("".join(c if c.isprintable() else " " for c in r["name"]).split()), "".join(r["sequence"].split()).upper())
            for r in rows]


def export_order(rows, path, fmt=ORDER_FORMATS[0], scale="25nm", purification="STD"):
    """Writes `rows` (primer-database rows) as an order file. Returns the files written: one for the plain CSV; for the IDT format
    one per 96 oligos (a single file keeps `path`, several become name_plate1.csv, name_plate2.csv …). ValueError when empty."""
    primers = order_rows(rows)
    if not primers:
        raise ValueError("There are no primers to export.")
    path = Path(path)
    if fmt == ORDER_FORMATS[1]:
        chunks = [primers]
    else:
        chunks = [primers[i:i + PLATE_SIZE] for i in range(0, len(primers), PLATE_SIZE)]
    written = []
    for n, chunk in enumerate(chunks, 1):
        target = path if len(chunks) == 1 else path.with_name(f"{path.stem}_plate{n}{path.suffix}")
        with open(target, "w", encoding="utf-8", newline="") as fh:
            writer = csv.writer(fh)
            if fmt == ORDER_FORMATS[1]:
                writer.writerows(chunk)
            else:
                writer.writerow(["Well Position", "Name", "Sequence", "Scale", "Purification"])
                writer.writerows((well(i), name, seq, scale, purification) for i, (name, seq) in enumerate(chunk))
        written.append(target)
    return written


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


def pair(row):
    """A database row's (name, sequence)."""
    return row["name"], row["sequence"]


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
            editor.rename_in_file(pair(source), pair(target))
            return None
        editor.rename_in_file(pair(before), pair(after))
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
