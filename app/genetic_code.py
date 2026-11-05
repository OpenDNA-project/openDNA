"""The genetic code in use: the Tools > Settings "Genetic code" tab, and everything that needs a table to translate.

The tables themselves (NCBI translation tables) are data in the Rust core (core/src/code.rs); this module reads them
once, offers them as dictionaries for the drawing code, and keeps three persistent settings:

- the table used for DNA that is not a CDS with its own /transl_table (3-frame translation, selection in the Protein
  pane, ORFs, protein-motif search) and for CDSs without /transl_table: default 1 (standard);
- "Allow alternative start codons" for the ORF finder (default off: ATG only, as before);
- "Keep /translation up to date" (Feature menu): see `keep_translation`.

A CDS that carries /transl_table=N always uses table N; /codon_start is honoured; a simple /transl_except
(selenocysteine ...) is applied. All of that is worked out by Rust (`Document.cds_plan`), so the drawing and the
/translation refresh cannot disagree.
"""
from types import SimpleNamespace

from genomeedit_core import genetic_codes

import tool_settings

TOOL = "Genetic code"

# id -> table: name, the 64 codon letters, a codon -> amino acid dict, and the start codons (all of the table's).
TABLES = {}
_BASES = "TCAG"
_CODONS = [a + b + c for a in _BASES for b in _BASES for c in _BASES]   # NCBI order, same as the Rust strings
for _id, _name, _aas, _starts in genetic_codes():
    TABLES[_id] = SimpleNamespace(
        id=_id, name=_name, aas=_aas, starts=_starts,
        codons=dict(zip(_CODONS, _aas)),
        start_codons=[c for c, s in zip(_CODONS, _starts) if s == "M"],
        stop_codons=[c for c, a in zip(_CODONS, _aas) if a == "*"])

LABELS = [f"{t.id} – {t.name}" for t in TABLES.values()]
DEFAULT_LABEL = LABELS[0]   # "1 – Standard"


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
])


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
