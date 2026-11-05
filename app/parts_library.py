"""The Features Database: the known features (parts) Detect Common Features and the Plasmid Checkup look for.

It is one library made of *sources*, one file each, so a source can be left out whole:
- `app/parts/*.tsv`: bundled with the program, one file per published source (named by its `#title`);
- `parts.tsv` in the app data folder: your own features, shown as the source "Local" (right-click a feature → Add to
  Features Database, Import → Features and Tracks → Import Features to Features Database). Each keeps the file it came from.

Databases → Features Database opens it as a tab (db_tabs.py): the Features view on it, source a column of its own.

A source file is tab-separated, with `#key: value` lines before the header row:

    #title: Terminators (Chen 2013)
    #reference: Chen YJ, et al. Nat Methods 10, 659-664 (2013).
    #doi: 10.1038/nmeth.2515
    #note: Termination strength (TS) measured in E. coli; see Help > Citations.
    name<TAB>type<TAB>sequence<TAB>organism<TAB>strength<TAB>units<TAB>source

`strength` is the number the source measured and `units` says what it means (TS, RPU, relative),
so values are never compared across sources. `organism` is the host it was measured in; a part with
no measurement leaves both empty. `source` is the record or file the feature was taken from. Only `name`,
`type` and `sequence` are required.

Which sources are used is remembered in settings.ini (`Parts Library/off`): the checkup and
Detect Common Features both read `load_parts()`, which skips the sources that are switched off.
"""
import csv
from pathlib import Path
from types import SimpleNamespace

from PySide6.QtCore import QObject, Signal

import tool_settings
from primer_db import data_dir

BUNDLED_DIR = Path(__file__).with_name("parts")
YOURS = "yours"          # the collection id of the user's own parts
LOCAL = "Local"          # ... and the source name they are shown under
SETTING = ("Parts Library", "off")   # settings.ini: the ids that are switched off, comma separated
FIELDS = ("name", "type", "sequence", "organism", "strength", "units", "source")


class _Notifier(QObject):
    changed = Signal()   # a source was switched on / off, or features were added: what is found and listed changes


notifier = _Notifier()


def source_name(cid, meta=None):
    """The name a source is listed under: "Local" for your own features, else the title of its file."""
    if cid == YOURS:
        return LOCAL
    return (meta or {}).get("title") or cid


SOURCE_SEP = " · "   # between a source's name and the record or file a feature was taken from, in one field


def source_label(part):
    """Where a feature comes from, in one field: "Cloning vectors (NCBI records) · pBR322 (J01749)", "Local · pMine.gb"."""
    return part.source_name + (SOURCE_SEP + part.source if part.source else "")


def source_of(label):
    """The source's name of a `source_label` (what the Source column is filtered by)."""
    return label.split(SOURCE_SEP)[0]


def user_library():
    return data_dir() / "parts.tsv"


def _meta_and_rows(path):
    """(`#key: value` header, csv.DictReader rows) of a collection file."""
    meta, lines = {}, []
    with open(path, newline="", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                key, _, v = line[1:].partition(":")
                meta[key.strip().lower()] = v.strip()
            else:
                lines.append(line)
    return meta, list(csv.DictReader(lines, delimiter="\t"))


def _part(row, collection):
    """One library part, or None when the row has no usable sequence."""
    seq = (row.get("sequence") or "").upper().strip()
    if not seq or set(seq) - set("ACGTN"):
        return None
    strength = row.get("strength") or ""
    try:
        strength = float(strength) if strength else None
    except ValueError:
        strength = None
    return SimpleNamespace(name=row.get("name") or "part", kind=row.get("type") or "misc_feature", sequence=seq,
                           organism=(row.get("organism") or "").strip(), strength=strength,
                           units=(row.get("units") or "").strip(), source=(row.get("source") or "").strip(),
                           collection=collection)


def files():
    """(collection id, path) of every collection there is, the user's own last (its names win)."""
    out = [(p.stem, p) for p in sorted(BUNDLED_DIR.glob("*.tsv"))] if BUNDLED_DIR.is_dir() else []
    if user_library().exists():
        out.append((YOURS, user_library()))
    return out


def collections():
    """What the Parts Library dialog lists: SimpleNamespace(id, title, reference, doi, note, organism,
    counts {type: n}, measured, on) per collection, in file order."""
    out = []
    for cid, path in files():
        try:
            meta, rows = _meta_and_rows(path)
        except OSError:
            continue
        parts = [p for p in (_part(r, cid) for r in rows) if p is not None]
        counts = {}
        for p in parts:
            counts[p.kind] = counts.get(p.kind, 0) + 1
        out.append(SimpleNamespace(
            id=cid, path=path, title=source_name(cid, meta),
            reference=meta.get("reference", ""), doi=meta.get("doi", ""), note=meta.get("note", ""),
            organism=meta.get("organism", ""), counts=counts, measured=sum(p.strength is not None for p in parts),
            on=is_on(cid)))
    return out


def all_parts():
    """Every usable part of every source, switched on or not, none merged away (what the Features Database tab lists):
    SimpleNamespace (see `_part`) with `source_name` too, in file order."""
    out = []
    for cid, path in files():
        try:
            meta, rows = _meta_and_rows(path)
        except (OSError, UnicodeDecodeError, csv.Error):
            continue
        name = source_name(cid, meta)
        for row in rows:
            p = _part(row, cid)
            if p is not None:
                p.source_name = name
                out.append(p)
    return out


def off_ids():
    raw = tool_settings.store().value(f"{SETTING[0]}/{SETTING[1]}", "")
    return {s for s in str(raw or "").split(",") if s}


def is_on(cid):
    return cid not in off_ids()


def set_on(cid, on):
    ids = off_ids() - {cid} if on else off_ids() | {cid}
    if ids != off_ids():
        tool_settings.store().setValue(f"{SETTING[0]}/{SETTING[1]}", ",".join(sorted(ids)))
        notifier.changed.emit()


def load_parts(kinds=None):
    """The parts of the collections that are switched on, as SimpleNamespace (see `_part`), the user's
    own last so its names win where the same sequence is in two collections. `kinds`: only these types."""
    parts = []
    for cid, path in files():
        if not is_on(cid):
            continue
        try:
            _, rows = _meta_and_rows(path)
        except OSError:
            continue
        for row in rows:
            p = _part(row, cid)
            if p is not None and (kinds is None or p.kind in kinds):
                parts.append(p)
    by_seq = {p.sequence: p for p in parts}   # a dict keeps the last one per sequence
    return list(by_seq.values())


def describe(part):
    """The measurement after a part's name, e.g. "TS 98 (E. coli)", or "" when it has none."""
    if part.strength is None:
        return ""
    number = f"{part.strength:g}"
    unit = f"{part.units} {number}" if part.units else number
    return f"{unit} ({part.organism})" if part.organism else unit


def part_problem(name, kind, seq, source=""):
    """Why a part cannot go in the library (a sentence), or "". The file is tab-separated, so a tab or line break in a text would
    break it; the search only knows A, C, G, T and N."""
    for what, text in (("name", name), ("type", kind), ("source", source)):
        if any(not c.isprintable() for c in text or ""):
            return f"its {what} has a tab, line break or other control character"
    if not name:
        return "it has no name"
    if not seq:
        return "it has no sequence"
    bad = sorted(set(seq.upper()) - set("ACGTN"))
    if bad:
        return f"its sequence has characters other than A, C, G, T, N ({' '.join(bad[:6])})"
    return ""


def add_parts(new, rejected=None):
    """Append (name, kind, sequence, source) rows to the user's collection, skipping sequences any
    collection already has. Returns how many were added. A part that cannot be stored (part_problem) is left out, and
    `rejected` (a list) receives (name, why) for it."""
    path = user_library()
    known = {p.sequence for p in load_parts()}
    rows = []
    for name, kind, seq, source in new:
        seq = seq.upper()
        problem = part_problem(name, kind, seq, source)
        if problem:
            if rejected is not None:
                rejected.append((name, problem))
            continue
        if seq not in known:
            known.add(seq)
            rows.append((name, kind, seq, "", "", "", source))
    if rows:
        fresh = not path.exists()
        with open(path, "a", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            if fresh:
                fh.write(f"#title: {LOCAL}\n")
                w.writerow(FIELDS)
            w.writerows(rows)
        notifier.changed.emit()
    return len(rows)


def remove_parts(keys):
    """Take your own (Local) features out of the library: `keys` = [(name, sequence)]. The file keeps its header lines and every
    other row as they were. Returns how many rows went."""
    path = user_library()
    keys = {(n, s.upper()) for n, s in keys}
    if not keys or not path.exists():
        return 0
    with open(path, newline="", encoding="utf-8") as fh:
        lines = fh.readlines()
    header = next((i for i, line in enumerate(lines) if not line.startswith("#")), len(lines))
    columns = lines[header].rstrip("\r\n").split("\t") if header < len(lines) else list(FIELDS)
    kept, gone = lines[:header + 1], 0
    for line in lines[header + 1:]:
        row = dict(zip(columns, line.rstrip("\r\n").split("\t")))
        if ((row.get("name") or "part"), (row.get("sequence") or "").upper().strip()) in keys:
            gone += 1
        else:
            kept.append(line)
    if gone:
        tmp = path.with_suffix(".tmp")   # write then rename: never a half-written library
        tmp.write_text("".join(kept), encoding="utf-8", newline="")
        tmp.replace(path)
        notifier.changed.emit()
    return gone
