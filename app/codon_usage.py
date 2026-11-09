"""Codon usage tables: how often an organism uses each of the 64 codons (per thousand codons). No GUI code here.

openDNA ships none: the published tables (Kazusa, HIVE-CUTs ...) have no licence that lets a program redistribute them. You bring
your own: import a CSV (`codon,aa,per_thousand`, or just codon and frequency), or make a table by hand in the codon editor. Each one is
a JSON file in `codon_usage/` of the data folder (the one place where all the app's data lives, see primer_db.data_dir).

A table is a SimpleNamespace(name, per_thousand {codon: float}, source, note). A codon the organism does not use at all has frequency 0;
which codons are *removed* from a genetic code is the code's business (genetic_code.py), not the table's.

`cai` and `rare_codons` score a coding sequence against a table: the usual codon adaptation index (Sharp & Li 1987) from each codon's
frequency relative to the most used synonymous codon."""
import csv
import json
import math
import re
from pathlib import Path
from types import SimpleNamespace

import fileio
import genetic_code
import tool_settings
from primer_db import data_dir

BASES = "TCAG"
CODONS = [a + b + c for a in BASES for b in BASES for c in BASES]
RARE = 0.1   # a codon used less than this fraction of the most used synonymous codon is "rare"


def _table(name, per_thousand, source="", note=""):
    return SimpleNamespace(name=name, per_thousand={c: float(per_thousand.get(c, 0.0)) for c in CODONS}, source=source, note=note)


# ---- reading and writing a CSV ---------------------------------------------------------------------------------------------
def read_csv(path, name=None):
    """A table from a CSV file: `# key: value` comment lines (name, source...), then `codon,aa,per_thousand[,count]`. Also accepts
    a plain list of codon,frequency pairs. Raises ValueError with a sentence when the file is not a codon table."""
    meta, rows = {}, {}
    text = Path(path).read_text(encoding="utf-8-sig")
    for line in text.splitlines():
        line = line.strip()
        if line.startswith("#"):
            if ":" in line:
                key, _, val = line[1:].partition(":")
                meta[key.strip().lower()] = val.strip()
            continue
        if not line:
            continue
        cells = [c.strip() for c in re.split(r"[,;\t]", line)]
        codon = cells[0].upper().replace("U", "T")
        if codon not in CODONS:
            if rows or any(re.fullmatch(r"\d+(\.\d+)?", c) for c in cells[:1]):
                raise ValueError(f"'{cells[0]}' is not a codon (three letters A, C, G, T).")
            continue   # the header line
        try:
            rows[codon] = float(cells[2] if len(cells) > 2 else cells[1])
        except (IndexError, ValueError):
            raise ValueError(f"No frequency for codon {codon}.") from None
    if len(rows) != 64:
        raise ValueError(f"A codon table has 64 codons; this file has {len(rows)}.")
    source = " · ".join(x for x in (meta.get("source"), meta.get("citation")) if x)   # a file may say where it comes from
    return _table(name or meta.get("name") or meta.get("organism") or Path(path).stem, rows, source, meta.get("note", ""))


def write_csv(table, path):
    code = genetic_code.codon_map(1)
    with open(path, "w", encoding="utf-8", newline="") as f:
        f.write(f"# name: {table.name}\n")
        if table.source:
            f.write(f"# source: {table.source}\n")
        if table.note:
            f.write(f"# note: {table.note}\n")
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["codon", "aa", "per_thousand"])
        for c in CODONS:
            w.writerow([c, code[c], f"{table.per_thousand[c]:.2f}"])


# ---- the tables ------------------------------------------------------------------------------------------------------------
def folder():
    f = data_dir() / "codon_usage"
    f.mkdir(parents=True, exist_ok=True)
    return f


def _file(name):
    return folder() / (re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("_") + ".json")


def user_tables():
    out = {}
    for path in sorted(folder().glob("*.json")):
        try:
            d = json.loads(path.read_text(encoding="utf-8"))
            out[d["name"]] = _table(d["name"], d["per_thousand"], d.get("source", ""), d.get("note", ""))
        except (OSError, ValueError, KeyError, TypeError):
            continue   # a file that is not a table is left alone
    return out


def all_tables():
    """name -> table: the tables you made or imported."""
    return user_tables()


def save_user(table):
    """Save `table` under its name (a table of that name is replaced). Raises ValueError for an empty name or a file that cannot be written."""
    if not table.name.strip():
        raise ValueError("A codon table needs a name.")
    tmp = _file(table.name).with_suffix(".tmp")
    try:
        tmp.write_text(json.dumps({"name": table.name, "source": table.source, "note": table.note,
                                   "per_thousand": table.per_thousand}, indent=1), encoding="utf-8")
        fileio.replace(tmp, _file(table.name))
    except OSError as e:
        raise ValueError(f"The table '{table.name}' could not be saved ({e}). Nothing was changed on disk.") from None
    tool_settings.notifier.changed.emit(genetic_code.TOOL)


def delete_user(name):
    _file(name).unlink(missing_ok=True)
    if active_name() == name:
        set_active("")
    tool_settings.notifier.changed.emit(genetic_code.TOOL)


def duplicate(table, name):
    return _table(name, table.per_thousand, table.source, f"Copy of {table.name}." + (f" {table.note}" if table.note else ""))


# ---- the table in use ------------------------------------------------------------------------------------------------------
_KEY = f"{genetic_code.TOOL}/usage_table"


def active_name():
    return str(tool_settings.store().value(_KEY, ""))


def set_active(name):
    tool_settings.store().setValue(_KEY, name)


def project_entry():
    """The table in use, for an Open DNA Project (None when there is none)."""
    t = active()
    return None if t is None else {"name": t.name, "per_thousand": t.per_thousand, "source": t.source, "note": t.note}


def adopt(entry):
    """Use the codon usage table a project was saved with; one that is missing here is added, and one of the same name with other
    numbers is kept beside it as "name (from project)". Returns a sentence for the status bar ('' when nothing changed)."""
    if not entry:
        return ""
    try:
        name = entry["name"]
        mine = _table(name, entry["per_thousand"], entry.get("source", ""), entry.get("note", ""))
        have = all_tables().get(name)
        if have is not None and any(abs(have.per_thousand[c] - mine.per_thousand[c]) > 1e-6 for c in CODONS):
            mine.name = name = f"{name} (from project)"
            have = all_tables().get(name)
        if have is None:
            save_user(mine)
    except (KeyError, TypeError, ValueError) as e:
        return f"The project's codon table could not be used ({e})." if isinstance(e, ValueError) else ""
    if name == active_name():
        return ""
    set_active(name)
    tool_settings.notifier.changed.emit(genetic_code.TOOL)
    return f"Codon usage table set to {name}, as in the project."


def active():
    """The table in use: the one chosen in Settings, or the first one you have when it is gone; None when you have none."""
    tables = all_tables()
    return tables.get(active_name()) or next(iter(tables.values()), None)


# ---- scoring a coding sequence ---------------------------------------------------------------------------------------------
def adaptiveness(table, code=None):
    """codon -> w, its frequency over the highest frequency of the codons of the same amino acid (0 for a codon never used, and the
    stops, Met and Trp, which have no synonym, are left out of the index)."""
    code = code or genetic_code.codon_map()
    best = {}
    for c, aa in code.items():
        best[aa] = max(best.get(aa, 0.0), table.per_thousand[c])
    return {c: (table.per_thousand[c] / best[aa] if best[aa] else 0.0) for c, aa in code.items() if aa != "*"}


def _scored(seq, table, code):
    code = code or genetic_code.codon_map()
    w = adaptiveness(table, code)
    synonymous = {aa: sum(1 for a in code.values() if a == aa) for aa in set(code.values())}
    seq = seq.upper().replace("U", "T")
    for i in range(0, len(seq) - len(seq) % 3, 3):
        c = seq[i:i + 3]
        if c in code and code[c] != "*" and synonymous[code[c]] > 1:
            yield i // 3, c, w[c]


def cai(seq, table, code=None):
    """The codon adaptation index of the coding sequence `seq` (0-1, geometric mean of w over its codons that have synonyms), or
    None when it has none. A codon the table never uses counts as 0.01 (Sharp & Li's floor) so one such codon does not give 0."""
    ws = [max(w, 0.01) for _, _, w in _scored(seq, table, code)]
    return math.exp(sum(math.log(w) for w in ws) / len(ws)) if ws else None


def rare_codons(seq, table, code=None, limit=RARE):
    """[(codon number from 1, codon)] of the codons of `seq` the table uses less than `limit` times as much as the best synonym."""
    return [(i + 1, c) for i, c, w in _scored(seq, table, code) if w < limit]
