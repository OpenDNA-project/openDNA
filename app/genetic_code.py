"""The genetic code in use: the Tools > Settings "Genetic code" tab, and everything that needs a table to translate.

The tables themselves (NCBI translation tables) are data in the Rust core (core/src/code.rs); this module reads them
once, offers them as dictionaries for the drawing code, and keeps these persistent settings:

- the table used for DNA that is not a CDS with its own /transl_table (3-frame translation, selection in the Protein
  pane, ORFs, protein-motif search) and for CDSs without /transl_table: default 1 (standard);
- "Allow alternative start codons" for the ORF finder (default off: ATG only, as before);
- "Keep /translation up to date" (Feature menu): see `keep_translation`;
- "Show the genetic code in the toolbar": the quick switch (code_switch.py).

Besides the NCBI tables there are codes of id 100 and up, registered in the Rust core when this module loads: the recoded ones
below and the ones you write in the codon editor (codon_dialog.py), one JSON file each in `genetic_codes/` of the data folder. In
such a code a codon that is no longer used (a recoded organism) has the letter "X" and is translated as X.

A CDS that carries /transl_table=N always uses table N; /codon_start is honoured; a simple /transl_except
(selenocysteine ...) is applied. All of that is worked out by Rust (`Document.cds_plan`), so the drawing and the
/translation refresh cannot disagree.
"""
import json
from types import SimpleNamespace

from opendna_core import genetic_codes, register_genetic_code

import fileio
import tool_settings
from primer_db import data_dir

TOOL = "Genetic code"

# id -> table: name, the 64 codon letters, a codon -> amino acid dict, and the start codons (all of the table's).
TABLES = {}
_BASES = "TCAG"
_CODONS = [a + b + c for a in _BASES for b in _BASES for c in _BASES]   # NCBI order, same as the Rust strings
for _id, _name, _aas, _starts in genetic_codes():
    TABLES[_id] = SimpleNamespace(
        id=_id, name=_name, aas=_aas, starts=_starts, builtin=True, source="NCBI translation table", note="",
        codons=dict(zip(_CODONS, _aas)),
        start_codons=[c for c, s in zip(_CODONS, _starts) if s == "M"],
        stop_codons=[c for c, a in zip(_CODONS, _aas) if a == "*"])

LABELS = [f"{t.id} – {t.name}" for t in TABLES.values()]
DEFAULT_LABEL = LABELS[0]   # "1 – Standard"

REMOVED = "X"        # the letter of a codon the organism no longer uses
LETTERS = "ACDEFGHIKLMNPQRSTVWY*" + REMOVED   # what a codon may code for: an amino acid, the stop, or nothing any more
FIRST_CUSTOM = 100   # ids of the codes that are not NCBI's: 100-109 are built in, yours start at 110
FIRST_USER = 110
LAST_ID = 254
# Built in: id -> (name, NCBI table it starts from, its removed codons, where it is described). Every other codon keeps its meaning,
# the removed ones are read as X.
RECODED = {
    101: ("Syn61 (recoded E. coli)", 11, ["TCG", "TCA", "TAG"],
          "Fredens et al., Nature 569, 514-518 (2019), doi:10.1038/s41586-019-1192-5"),
    102: ("Syn57 (recoded E. coli)", 11, ["TCG", "TCA", "TCC", "TCT", "GCA", "GCG", "TAG"],
          "Robertson et al., Science 390, eady4368 (2025), doi:10.1126/science.ady4368 (read from the preprint, bioRxiv 2025.05.02.651837)"),
}


def _label(t):
    return f"{t.id} – {t.name}"


def _add(id, name, aas, starts, builtin, source="", note=""):
    """Register a code of id 100+ here and in the Rust core, and list it (replacing the code of that id if there is one)."""
    register_genetic_code(id, name, aas, starts)
    old = TABLES.get(id)
    TABLES[id] = SimpleNamespace(
        id=id, name=name, aas=aas, starts=starts, builtin=builtin, source=source, note=note,
        codons=dict(zip(_CODONS, aas)),
        start_codons=[c for c, s in zip(_CODONS, starts) if s == "M"],
        stop_codons=[c for c, a in zip(_CODONS, aas) if a == "*"])
    if old is not None and _label(old) in LABELS:
        LABELS[LABELS.index(_label(old))] = _label(TABLES[id])
    else:
        LABELS.append(_label(TABLES[id]))   # (the list itself, not a copy: the Settings page and the quick switch hold it)
    return TABLES[id]


def removed_codons(table):
    """The codons the code `table` (an id) does not use any more, in TCAG order."""
    return [c for c, a in TABLES.get(table, TABLES[1]).codons.items() if a == REMOVED]


def is_custom(table):
    return table >= FIRST_CUSTOM


def recoded(base, removed):
    """(aas, starts) of NCBI table `base` with the `removed` codons not used any more."""
    t = TABLES[base]
    aas = "".join(REMOVED if c in removed else a for c, a in zip(_CODONS, t.aas))
    return aas, "".join("-" if c in removed else s for c, s in zip(_CODONS, t.starts))


def folder():
    f = data_dir() / "genetic_codes"
    f.mkdir(parents=True, exist_ok=True)
    return f


def _file(id):
    return folder() / f"{id}.json"


def save_code(name, aas, starts, id=None, note=""):
    """Add a code of yours (or change it, given its `id`) and use its id; ValueError with a sentence when it cannot be saved."""
    name = name.strip()
    if not name:
        raise ValueError("A genetic code needs a name.")
    if len(aas) != 64 or len(starts) != 64 or not aas.isascii() or not starts.isascii():
        raise ValueError("A genetic code is 64 amino-acid letters and 64 start marks, in TCAG order.")
    for other in TABLES.values():
        if other.name.lower() == name.lower() and other.id != id:
            raise ValueError(f"There is already a genetic code called '{name}'.")
    if id is None:
        id = next((i for i in range(FIRST_USER, LAST_ID + 1) if i not in TABLES), None)
        if id is None:
            raise ValueError("No room for another genetic code of your own (145 at most).")
    elif id < FIRST_USER or id > LAST_ID or (id in TABLES and TABLES[id].builtin):
        raise ValueError("A built-in genetic code cannot be changed: duplicate it first.")
    tmp = _file(id).with_suffix(".tmp")
    try:
        tmp.write_text(json.dumps({"id": id, "name": name, "aas": aas, "starts": starts, "note": note}, indent=1), encoding="utf-8")
        fileio.replace(tmp, _file(id))
    except OSError as e:
        raise ValueError(f"The genetic code '{name}' could not be saved ({e}). Nothing was changed on disk.") from None
    _add(id, name, aas, starts, False, note=note)
    tool_settings.notifier.changed.emit(TOOL)
    return id


def delete_code(id):
    """Delete a code of yours; the setting goes back to the standard code if it was the one in use."""
    if id < FIRST_USER or id not in TABLES:
        return
    if table_id() == id:
        tool_settings.set_value(TOOL, "table", DEFAULT_LABEL)
    _file(id).unlink(missing_ok=True)
    gone = TABLES.pop(id)
    LABELS.remove(_label(gone))
    tool_settings.notifier.changed.emit(TOOL)


def project_entry():
    """The genetic code in use, for an Open DNA Project: its id and name, and its letters when it is one of yours (the others are
    the same on every machine)."""
    t = TABLES[table_id()]
    return {"id": t.id, "name": t.name, **({"aas": t.aas, "starts": t.starts} if t.id >= FIRST_USER else {})}


def adopt(entry):
    """Use the genetic code a project was saved with; a code of yours that is missing here is added (matched by name: the ids differ
    between machines). Returns a sentence for the status bar ('' when nothing changed or the entry is unusable)."""
    try:
        id, name = entry["id"], entry["name"]
        if id < FIRST_USER:
            if id not in TABLES:
                return ""
        else:
            same = next((t for t in TABLES.values() if t.id >= FIRST_USER and t.name.lower() == name.lower()), None)
            if same is not None and (same.aas, same.starts) != (entry["aas"], entry["starts"]):
                name += " (from project)"   # a code of that name that is not this one: both are kept
                same = next((t for t in TABLES.values() if t.id >= FIRST_USER and t.name.lower() == name.lower()), None)
            id = same.id if same is not None else save_code(name, entry["aas"], entry["starts"], note="From a project.")
    except (KeyError, TypeError, ValueError) as e:
        return f"The project's genetic code could not be used ({e})." if isinstance(e, ValueError) else ""
    if id == table_id():
        return ""
    tool_settings.set_value(TOOL, "table", _label(TABLES[id]))
    tool_settings.notifier.changed.emit(TOOL)
    return f"Genetic code set to {_label(TABLES[id])}, as in the project."


def _load_custom():
    for id, (name, base, removed, source) in RECODED.items():
        _add(id, name, *recoded(base, removed), True, source=source)
    for path in sorted(folder().glob("*.json")):
        try:
            d = json.loads(path.read_text(encoding="utf-8"))
            if FIRST_USER <= d["id"] <= LAST_ID and len(d["aas"]) == 64 and len(d["starts"]) == 64 and d["id"] not in TABLES:
                _add(d["id"], d["name"], d["aas"], d["starts"], False, note=d.get("note", ""))
        except (OSError, ValueError, KeyError, TypeError):
            continue   # a file that is not a code is left alone


def describe(label):
    """The read-only line under the table choice: its stop and start codons (RNA letters are not used: this is DNA)."""
    t = TABLES.get(_id_of(label), TABLES[1])
    return f"Stop codons: {', '.join(t.stop_codons) or 'none'}  ·  Start codons: {', '.join(t.start_codons)}"


def _id_of(label):
    try:
        return int(str(label).split()[0])
    except (ValueError, IndexError):
        return 1


tool_settings.register(TOOL, [
    tool_settings.Param("table", "Genetic code", DEFAULT_LABEL, choices=LABELS, info=describe,
                        tip="The NCBI translation table used to translate DNA: the 3-frame translation, a selection in the Protein pane, "
                            "ORFs, protein-motif search. A CDS with its own /transl_table keeps that table."),
    tool_settings.Param("alt_starts", "Allow alternative start codons", False,
                        tip="Detect ORFs: also start at the table's other start codons (GTG, TTG, CTG, ATT...), not only ATG. "
                            "Off by default."),
    tool_settings.Param("show_switch", "Show the genetic code in the toolbar", True,
                        tip="A quick switch left of the Select box that shows the table in use and changes it."),
])


_load_custom()


def table_id():
    """NCBI id of the table chosen in the settings (1 = standard)."""
    return _id_of(tool_settings.value(TOOL, "table"))


def alt_starts():
    return bool(tool_settings.value(TOOL, "alt_starts"))


def codon_map(table=None):
    """codon -> amino acid ('*' = stop) of `table` (default: the setting). Unknown ids give the standard table."""
    return TABLES.get(table or table_id(), TABLES[1]).codons


def name_of(table):
    t = TABLES.get(table, TABLES[1])
    return f"Table {t.id}"


def table_note(table=None):
    """' · Table 11' for the headers when the table is not the standard one, else ''."""
    table = table or table_id()
    return f" · {name_of(table)}" if table != 1 else ""


# ---- one CDS ---------------------------------------------------------------------------------------------------------------
def cds_plan(doc, f):
    """How CDS feature `f` is translated: SimpleNamespace(table, frame, first_met, exceptions {codon number: letter},
    why). `why` is '' normally; otherwise the reason the annotation could not be followed in full (an unusual location, a
    /transl_except we cannot reproduce): the translation is then drawn plainly, from /codon_start and the table only."""
    table, frame, first_met, exceptions, why = doc.cds_plan(f.index, table_id())
    if why:
        quals = dict(f.qualifiers)
        start = quals.get("codon_start")
        frame = int(start) - 1 if start in ("1", "2", "3") else 0
        own = quals.get("transl_table")
        table = int(own) if own and own.isdigit() and int(own) in TABLES else table_id()
        first_met = False
    return SimpleNamespace(table=table, frame=frame, first_met=first_met, exceptions=dict(exceptions), why=why)


def cds_note(doc, f):
    """One line for the tooltip of a CDS: its table when not the standard one, and 'translation not reproducible' when it is."""
    plan = cds_plan(doc, f)
    parts = []
    if plan.table != 1:
        parts.append(name_of(plan.table))
    if plan.why:
        parts.append(f"translation not reproducible: {plan.why}")
    return " · ".join(parts)


# ---- Features > Keep /translation up to date -------------------------------------------------------------------------------
def keep_translation():
    """Whether saving a GenBank file refreshes the /translation of its CDSs (default on)."""
    raw = tool_settings.store().value("features/keep_translation", True)
    return raw if isinstance(raw, bool) else str(raw).lower() in ("true", "1")


def set_keep_translation(on):
    tool_settings.store().setValue("features/keep_translation", bool(on))


def summarize(report):
    """Status-bar text for a refresh report [(index, status, reason)]: 'Updated /translation of 3 CDSs, skipped 1 pseudo'."""
    updated = sum(s == "updated" for _, s, _ in report)
    skipped = {}
    for _, status, _ in report:
        if status not in ("updated", "unchanged"):
            kind = {"other": "not reproducible"}.get(status, status)
            skipped[kind] = skipped.get(kind, 0) + 1
    text = f"Updated /translation of {updated} CDS{'s' if updated != 1 else ''}"
    if skipped:
        text += ", skipped " + ", ".join(f"{n} {kind}" for kind, n in sorted(skipped.items()))
    return text
