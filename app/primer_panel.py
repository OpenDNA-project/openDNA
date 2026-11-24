"""Primers: the Primers tab (sites on the open sequence, grouped by primer), the primer
dialog (type a sequence, see live where it binds), the primer database window, options.

A primer is just a name + sequence in the primer database (primer_db.py). Where it binds
is never stored: Editor.match_primers recomputes it on the current sequence after every
edit, including sites with a few mismatches ("almost binds", orange, the mismatched bases written on the line).
"""
from html import escape
import threading
from pathlib import Path

from PySide6.QtCore import QEventLoop, QObject, QPointF, QSize, Qt, QTimer, Signal
from PySide6.QtGui import QColor, QFont, QFontMetrics, QPainter, QPalette, QPolygonF
from PySide6.QtWidgets import (
    QAbstractItemView, QCheckBox, QComboBox, QDialog, QDialogButtonBox, QFileDialog, QFrame,
    QApplication, QFormLayout, QHBoxLayout, QLabel, QLineEdit, QMenu, QMessageBox, QProgressDialog, QPushButton,
    QStyle, QStyledItemDelegate, QStyleOptionViewItem, QTextBrowser, QTreeWidget, QTreeWidgetItem, QVBoxLayout, QWidget,
)

from opendna_core import dimer, hairpin, site_tm as duplex_tm, tm, tm_range

import tool_settings
from editor import MAX_SITES, dna_field, reverse_complement
from file_tabs import page_of
from value_filter import SortItem
from tool_settings import Param, SettingsDialog
from reference import Reference, as_reference
from sequence_view import TEMP_COLOR
from uikit import dispose_progress, fit_columns, plural, show_after
from primer_db import (ORDER_FORMATS, PLATE_SIZE, PURIFICATIONS, SCALES, export_order, conditions, tm_rule, MAX_NAME, MAX_NOTE, MAX_PRIMER, add_primer, add_primers, BINDING, binding, db,  # noqa: F401 (re-exported)
                       primer_problem, update_primer)
import theme

MAX_LISTED = 8     # positions listed in a report before "... and N more"


def melting_temp(seq):
    """Tm in °C under the current PCR conditions, or None (non-ACGT or too short)."""
    return tm(seq, **conditions())


def is_degenerate(seq):
    """True when the primer has IUPAC codes (R, Y, N...): it is a mixture of oligos, not one."""
    return any(c not in "ACGTUacgtu" for c in seq)


def tm_text(seq):
    """The Tm to show: "59.3", or for a degenerate primer the span of its mixture, "58.1–61.4"; "" when there is none."""
    if is_degenerate(seq):
        r = tm_range(seq, **conditions())
        return "" if r is None else f"{r[0]:.1f}–{r[1]:.1f}"
    t = melting_temp(seq)
    return "" if t is None else f"{t:.1f}"


def site_tm_text(bound, template):
    """(text, value) of the Tm of a site: the primer's annealed bases `bound` on the template's bases there `template` (the
    primer's sense; None: unknown), its mismatches counted. A degenerate primer, or a template that is not plain bases: the
    primer's own Tm (a range for a mixture), value None."""
    t = duplex_tm(bound, template, **conditions()) if template is not None and len(template) == len(bound) else None
    return (tm_text(bound), None) if t is None else (f"{t:.1f}", t)


def min_tm_text():
    """", Tm ≥ 40 °C" when a minimum Tm of a site applies (what the checks count as a match), else ""."""
    low = tool_settings.value(BINDING, "min_tm")
    return f", Tm ≥ {low:g} °C" if low > 0 else ""


def tm_sort_key(seq):
    """Number to sort a Tm column by (the lowest Tm of a degenerate primer)."""
    r = tm_range(seq, **conditions()) if is_degenerate(seq) else None
    t = r[0] if r else melting_temp(seq)
    return -1e9 if t is None else t


def structure(seq):
    """Secondary structure ΔG (kcal/mol, 37 °C) and warnings for one primer.
    Limits are the usual rules of thumb: hairpin 3' < -2 / any < -3, self-dimer 3' < -5 /
    any < -9."""
    c = conditions()
    salt = {k: c[k] for k in ("na_mm", "mg_mm", "dntp_mm")}
    hp, hp3 = hairpin(seq, **salt)
    sd, sd3 = dimer(seq, seq, **salt)
    warnings = []
    if hp3 < -2 or hp < -3:
        warnings.append(f"hairpin {min(hp, hp3):.1f}")
    if sd3 < -5:
        warnings.append(f"3′ self-dimer {sd3:.1f}")
    elif sd < -9:
        warnings.append(f"self-dimer {sd:.1f}")
    return {"hairpin": hp, "self": sd, "self3": sd3, "warnings": warnings}


def structure_text(seq):
    if is_degenerate(seq):   # the pairing model only knows A C G T: a made-up number would mislead
        return "hairpin and dimer ΔG not calculated for a degenerate primer"
    s = structure(seq)
    text = f"hairpin {s['hairpin']:.1f}, self-dimer {s['self']:.1f} (3′ end {s['self3']:.1f}) kcal/mol"
    return text + ("  ⚠ " + "; ".join(s["warnings"]) if s["warnings"] else "  ✓")


def homology_settings():
    """(3'-end bases that must match, mismatches allowed in them) for the mis-priming check.
    Default 12 nt exact: ~0.5 chance hits per primer in a 5 Mb genome (1 mismatch: ~15)."""
    v = tool_settings.values(BINDING)
    return v["homology_seed"], v["homology_mismatches"]


def genome_path():
    """Off-target genome file, or "" when the check is off / none is set (the default: no background genome)."""
    v = tool_settings.values(BINDING)
    return v["genome"] if v["genome_check"] else ""


_genomes = {}  # (path, modified time) -> Reference, so the genome file is opened (and indexed) once per session


def genome():
    """The background genome (a `reference.Reference`: one or many sequences), or None when off / unreadable."""
    path = genome_path()
    if not path:
        return None
    try:
        key = (path, Path(path).stat().st_mtime_ns)
    except OSError:
        return None
    if key not in _genomes:
        try:
            _genomes[key] = Reference(path)
        except (OSError, ValueError):
            return None
    return _genomes[key]


def genome_label_text(g):
    """'<b>BW25113</b> (4.6 Mb)' or '<b>24 sequences, 3.1 Gb</b>': what the Primers and Spacer Finder tools show."""
    return f"<b>{escape(g.label())}</b>"


def annealed(editor, f):
    """Bases a primer_bind feature covers, written as primer sequence."""
    bases = "".join(editor.doc.sequence(s, e) for s, e in f.spans)
    return reverse_complement(bases) if f.strand < 0 else bases


def location_text(start, end, strand, length):
    """'2,001–2,022 →' (1-based; end may pass the origin of a circular sequence)."""
    return f"{start + 1:,}–{(end - 1) % length + 1:,} {'→' if strand > 0 else '←'}"


def site_text(start, end, strand, mismatches, length):
    """'2,001–2,022 → exact'"""
    match = "exact" if mismatches == 0 else plural(mismatches, "mismatch")
    return f"{location_text(start, end, strand, length)} {match}"


def three_prime_end(site):
    """Position of the primer's 3' end at a site (same for a short or a full match)."""
    start, end, strand, _ = site
    return (end if strand > 0 else start, strand)


def _listing(items):
    more = f" … and {len(items) - MAX_LISTED} more" if len(items) > MAX_LISTED else ""
    return "; ".join(items[:MAX_LISTED]) + more


CHECKS = "Primer checks"  # the thresholds: Tools → Settings (also from Primers → Options)
tool_settings.register(CHECKS, [
    Param("tm_low", "Tm ✓ from", 52.0, 30, 90, step=0.5, suffix=" °C"),
    Param("tm_high", "Tm ✓ up to", 65.0, 30, 90, step=0.5, suffix=" °C"),
    Param("tm_margin", "Tm ? within this much outside (else ✗)", 7.0, 0, 30, step=0.5, suffix=" °C"),
    Param("gc_low", "GC ✓ from", 40, 0, 100, suffix=" %"),
    Param("gc_high", "GC ✓ up to", 60, 0, 100, suffix=" %"),
    Param("gc_margin", "GC ? within this much outside (else ✗)", 10, 0, 50, suffix=" %"),
    Param("clamp_min", "3′ clamp ✓: at least … G/C in the last 5 bases", 1, 0, 5),
    Param("clamp_max", "3′ clamp ✓: at most … G/C in the last 5 bases", 3, 0, 5),
    Param("hairpin_bad", "Hairpin ✗ below", -6.0, -20, 0, step=0.5, suffix=" kcal/mol",
          tip="? from the usual limits: 3′ end below -2 or anywhere below -3 kcal/mol."),
    Param("dimer3_bad", "Self-dimer ✗: 3′ end below", -8.0, -20, 0, step=0.5, suffix=" kcal/mol",
          tip="? from the usual limits: 3′ end below -5 or anywhere below -9 kcal/mol."),
    Param("dimer_bad", "Self-dimer ✗: anywhere below", -12.0, -30, 0, step=0.5, suffix=" kcal/mol"),
])


def limit(key):
    return tool_settings.value(CHECKS, key)


def band(value, low, high, margin):
    """"ok" inside [low, high], "warn" within `margin` outside it, else "bad"."""
    return "ok" if low <= value <= high else "warn" if low - margin <= value <= high + margin else "bad"


# Check results: "ok" ✓ green, "warn" ? yellow, "bad" ✗ red, "info" no verdict.
# (✓ / ✗ rather than ✔ / ✖: Windows draws those as coloured emoji, ignoring the colour.)
MARKS = theme.Live({"ok": ("✓", "status_ok"), "warn": ("?", "status_warn"), "bad": ("✗", "status_bad"), "info": ("·", "status_info")})


def genome_query(editor, sequence):
    """(reference, annealing part, mismatches, seed, record to leave out) for the genome check of `sequence`, or None."""
    g = as_reference(genome())   # a Document (a sequence built in code) works as a one-sequence reference
    if g is None:
        return None
    seed, mm = homology_settings()
    return g, annealing_part(sequence), mm, seed, (editor.doc.name, editor.length)


def annealing_part(sequence):
    """The annealing part of a primer: a lowercase 5' tail does not need to match anywhere."""
    return sequence[len(sequence) - len(sequence.lstrip("acgtnryswkmbdhv")):] or sequence


def genome_search(editor, sequence, progress=None, cancel=None):
    """Run (or fetch from the cache) the search of `sequence` in every sequence of the reference. None = no genome."""
    q = genome_query(editor, sequence)
    if q is None:
        return None
    g, core, mm, seed, skip = q
    return g.search(sequence, core, mm, seed, skip, progress, cancel, warm=tm_rule())


def genome_cached(editor, sequence):
    q = genome_query(editor, sequence)
    if q is None:
        return None
    g, core, mm, seed, skip = q
    return g._cache.get(g.cache_key(sequence, core, mm, seed, skip, tm_rule()))


def site_in_record(h):
    """'chr3:1,200–1,221 →' for a genome hit (name, start, end, strand, mismatches, length, circular)."""
    name, start, end, strand, _, length, circular = h
    text = location_text(start, end, strand, length) if circular else f"{start + 1:,}–{end:,} {'→' if strand > 0 else '←'}"
    return f"{name}:{text}"


def genome_check(editor, sequence, mode="wait"):
    """The ("Non-specific (genome)", level, text) check: red when the whole primer matches somewhere in the reference, yellow
    when only its 3' end does, green when nowhere; the open sequence is never compared with itself."""
    name = "Non-specific (genome)"
    q = genome_query(editor, sequence)
    if q is None:
        return (name, "info", "no genome set (Primers → Options)")
    g, core, mm, seed, skip = q
    if all((n, l) == skip for n, l, _ in g.records):
        return (name, "info", "this sequence is the genome")
    r = genome_cached(editor, sequence) if mode == "cached" else genome_search(editor, sequence)
    if r is None:
        return (name, "info", f"searching {plural(g.count, 'sequence')}…")
    rule = f"last {seed} nt, {'exact' if mm == 0 else f'≤{mm} mismatch'}" + min_tm_text()
    note = ""
    if not r.complete:
        note = (f"  Only {r.searched} of {r.total} sequences were searched ("
                f"{'time limit' if r.stopped == 'time' else 'stopped'}): not conclusive.")
    where = g.name if g.single else "the genome"
    if r.whole:
        return (name, "bad", f"{where}: the whole primer binds " + _listing([site_in_record(h) for h in r.whole]) + note)
    if r.hits:
        return (name, "warn", f"{where}: 3′ end ({rule}) matches {plural(r.n_hits, 'place')}: "
                + _listing([site_in_record(h) for h in r.hits]) + note)
    if not r.complete:
        return (name, "warn", f"no 3′ match ({rule}) in the part searched." + note)
    return (name, "ok", f"{where}: no 3′ match ({rule})" + ("" if g.single else f" in {r.total} sequences"))


def warm_genome(parent, editor, sequences):
    """Search the background genome for each of `sequences` before their checks are drawn (a big genome takes a while): a
    progress dialog names the sequence being read and has a Cancel button. The finished searches are cached, so
    `primer_checks(..., "cached")` then answers at once; after Cancel the genome check just says "searching…"."""
    todo = [q for q in sequences if len(q) >= 8 and genome_query(editor, q) is not None and genome_cached(editor, q) is None]
    if not todo:
        return
    state = {"text": "Searching the background genome…", "cancel": False}

    def work():
        for i, q in enumerate(todo, 1):
            if state["cancel"]:
                return
            prefix = f"Primer {i} of {len(todo)}: " if len(todo) > 1 else ""
            genome_search(editor, q, lambda done, total, name, n, count, p=prefix: state.update(
                text=f"{p}searching {name} ({n} of {count}), {100 * done // max(total, 1)} %"), lambda: state["cancel"])
    thread = threading.Thread(target=work, daemon=True)
    thread.start()
    dialog = QProgressDialog(state["text"], "Cancel", 0, 0, parent)
    dialog.setWindowTitle("Background Genome")
    dialog.setWindowModality(Qt.WindowModal)
    show_after(dialog, 400)
    dialog.canceled.connect(lambda: state.update(cancel=True))
    while thread.is_alive():
        QApplication.processEvents(QEventLoop.AllEvents, 30)
        dialog.setLabelText(state["text"])
        thread.join(0.02)
    dispose_progress(dialog)   # (close() left its show timer running: it could come back over the next dialog)


class _GenomeSignals(QObject):
    progress = Signal(int, str)   # token, text
    done = Signal(int)            # token
    failed = Signal(int, str)     # token, what went wrong


class ReportView(QTextBrowser):
    """A check report (rich text) in a box of a fixed number of lines: a progress tick or another row clicked changes its text,
    never its height, so the table or the fields above it stay where they are (a longer report scrolls). Reads like a label:
    no frame, the window's background."""

    def __init__(self, lines):
        super().__init__(frameShape=QFrame.NoFrame, openLinks=False)
        self.viewport().setAutoFillBackground(False)
        self.setFixedHeight(lines * self.fontMetrics().lineSpacing() + 2 * int(self.document().documentMargin()) + 4)
        self._html = ""

    def setText(self, html):
        self._html = html
        self.setHtml(html)

    def text(self):
        return self._html


def primer_checks(editor, sequence, genome_mode="wait"):
    """[(check, level, text)] for one primer: Tm, GC, 3′ clamp, hairpin, self-dimer, where
    it binds here, and non-specific priming in this sequence and in the genome (red: the
    whole primer matches elsewhere; yellow: only its 3′ end does). The genome is searched record by record (every chromosome);
    `genome_mode` "wait": search now (a big genome takes a while), "cached": use the finished search or show "searching…"."""
    if len(sequence) < 8:
        return [("Primer", "info", "Type at least 8 bases.")]
    checks = []
    t = tm_sort_key(sequence)
    if t > -1e8:   # a degenerate primer is judged on its lowest Tm
        shown = tm_text(sequence)
        checks.append(("Tm", band(t, limit("tm_low"), limit("tm_high"), limit("tm_margin")),
                       f"{shown} °C ({len(sequence)} nt)"))
    up = sequence.upper()
    gc = 100 * sum(c in "GCS" for c in up) / len(up)
    checks.append(("GC", band(gc, limit("gc_low"), limit("gc_high"), limit("gc_margin")), f"{gc:.0f} %"))
    clamp, lo, hi = sum(c in "GC" for c in up[-5:]), limit("clamp_min"), limit("clamp_max")
    checks.append(("3′ clamp", "ok" if lo <= clamp <= hi else "warn",
                   f"{clamp} G/C in the last 5 bases" + ("" if lo <= clamp <= hi else f" ({lo}–{hi} is best)")))
    if is_degenerate(sequence):   # the pairing model only knows A C G T
        checks.append(("Hairpin / dimer", "info", "not calculated for a degenerate primer"))
    else:
        s = structure(sequence)
        hp = min(s["hairpin"], 0.0)
        checks.append(("Hairpin", "bad" if s["hairpin"] < limit("hairpin_bad") else
                       "warn" if any(w.startswith("hairpin") for w in s["warnings"]) else "ok", f"ΔG {hp:.1f} kcal/mol"))
        checks.append(("Self-dimer", "bad" if s["self3"] < limit("dimer3_bad") or s["self"] < limit("dimer_bad") else
                       "warn" if any("dimer" in w for w in s["warnings"]) else "ok",
                       f"ΔG {s['self']:.1f} (3′ end {s['self3']:.1f}) kcal/mol"))

    n = editor.length
    seed3, exact, mismatches = binding()
    warm = tm_rule()   # (every site below the minimum Tm left out, as where they are drawn)
    sites = editor.doc.binding_sites(sequence, mismatches, seed3, exact, **warm)
    checks.append(("Binds here", "ok" if sites else "warn",
                   _listing([site_text(*x, n) for x in sites]) if sites else "nowhere in this sequence"))
    seed, mm = homology_settings()
    rule = f"last {seed} nt, {'exact' if mm == 0 else f'≤{mm} mismatch'}" + min_tm_text()
    core = annealing_part(sequence)
    full = editor.doc.binding_sites(core, 0, len(core), **warm)  # the whole annealing part, exact
    bound = {three_prime_end(x) for x in sites}
    others = [h for h in editor.doc.binding_sites(sequence, mm, seed, **warm) if three_prime_end(h) not in bound]
    if len(full) > 1:
        checks.append(("Non-specific (sequence)", "bad", f"the whole primer binds {len(full)} places"))
    elif others:
        checks.append(("Non-specific (sequence)", "warn",
                       f"3′ end ({rule}) also matches: " + _listing([site_text(*h, n) for h in others])))
    else:
        checks.append(("Non-specific (sequence)", "ok", f"no other 3′ match ({rule})"))
    checks.append(genome_check(editor, sequence, genome_mode))
    return checks


def checks_html(checks):
    """The checks as a small table: mark, check, details."""
    rows = []
    for name, level, text in checks:
        mark, color = MARKS[level]
        rows.append(f"<tr><td><b><font color='{color}'>{mark}</font></b>&nbsp;</td>"
                    f"<td><b>{escape(name)}</b>&nbsp;&nbsp;</td><td>{escape(text)}</td></tr>")
    return "<table cellspacing='1' cellpadding='1'>" + "".join(rows) + "</table>"


def worst(checks):
    """The most serious level among checks (for one mark per primer)."""
    order = ["bad", "warn", "ok", "info"]
    return min((level for _, level, _ in checks), key=order.index, default="info")


def binding_report(editor, sequence):
    """Plain-text version of the checks, one per line ("Binds here: …")."""
    return "\n".join(f"{name}: {text}" for name, _, text in primer_checks(editor, sequence))


HEADER_WORDS = {"name", "primer", "primer name", "id", "oligo", "sequence", "seq", "sequence (5'-3')", "sequence 5'-3'",
                "sequence (5'→3')", "primer sequence"}


def parse_primer_list(text, rejected=None):
    """[(name, sequence)] from a primer list, every primer checked (primer_problem). Two formats:

    - FASTA: ">name optional words" then the sequence on one or more lines (spaces and line breaks are removed);
    - a table, one primer per line: name, then sequence, separated by a tab, a comma or a semicolon (the first of these found in the
      line); more columns are ignored, quotes around a value are removed. A first line naming the columns (name, sequence...) is
      skipped, and so are blank lines and lines starting with #.

    Every other line is refused, never guessed: `rejected` (a list) receives (line number, the line, why) for each."""
    out = []
    rejected = [] if rejected is None else rejected
    lines = text.splitlines()
    if text.lstrip().startswith(">"):
        name, start, body = None, 0, []

        def close():
            if name is None:
                return
            seq = "".join("".join(body).split())
            problem = primer_problem(name, seq)
            if problem:
                rejected.append((start, ">" + name, problem))
            else:
                out.append((name, seq))
        for number, line in enumerate(lines, 1):
            if line.startswith(">"):
                close()
                words = line[1:].split()
                name, start, body = (words[0] if words else ""), number, []
            elif name is not None:
                body.append(line)
            elif line.strip():
                rejected.append((number, line, "Text before the first '>' header of a FASTA file."))
        close()
        return out
    first = True
    for number, line in enumerate(lines, 1):
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        sep = next((s for s in ("\t", ",", ";") if s in line), None)
        if sep is None:
            first = False
            rejected.append((number, line, "No tab, comma or semicolon between a name and a sequence."))
            continue
        name, seq = (x.strip().strip('"').strip() for x in line.split(sep)[:2])
        if first and (name.lower() in HEADER_WORDS or seq.lower() in HEADER_WORDS):
            first = False
            continue   # the column names
        first = False
        seq = "".join(seq.split())
        problem = primer_problem(name, seq)
        if problem:
            rejected.append((number, line, problem))
        else:
            out.append((name, seq))
    return out


def import_primers(path, editor=None, rejected=None):
    """Add primers from a FASTA/CSV/TSV file to the database (one step in `editor`'s history). Returns a report in words, which names
    each refused line and why; `rejected` (a list) also receives them as (line number, line, why)."""
    rejected = [] if rejected is None else rejected
    with open(path, encoding="utf-8-sig") as fh:
        text = fh.read()
    if "\0" in text:
        raise ValueError(f"{Path(path).name} is not a text file (it holds binary data): a primer list is a TSV, CSV or FASTA file.")
    primers = parse_primer_list(text, rejected)
    report = add_primers(editor, primers, "Imported")
    out = f"Imported {len(primers)} primers: {report['added']} new, {report['exists']} already in the database."
    if report["conflict"]:   # the first few names: a long list made the message wider than the screen
        names = report["conflict"]
        out += (f" Name already used for another sequence (skipped): {', '.join(names[:5])}"
                + (f" and {len(names) - 5} more" if len(names) > 5 else ""))
    if rejected:
        out += f"\n\n{plural(len(rejected), 'line')} not imported:\n" + "\n".join(
            f"• line {n}: {line.strip()[:40]!r}: {why}" for n, line, why in rejected[:15])
        if len(rejected) > 15:
            out += f"\n… and {len(rejected) - 15} more."
    return out


class PrimerDialog(QDialog):
    """New primer (prefilled from the selection) or edit one from the database. While you
    type, shows where it binds here and where its 3' end could mis-prime."""

    def __init__(self, parent, editor, name=None):
        super().__init__(parent)
        self.editor, self.old = editor, db().find(name) if name else None
        # A primer drawn here that is not in the database: one of the file's own is edited in the file; a temporary one (an
        # assembly's, designed and not saved) keeps its name and sequence until Save puts it in the database (a rename here
        # added a second primer beside it).
        drawn = None if self.old or not name else next((p for p in editor.primer_list() if p[0] == name), None)
        self.in_file = drawn if drawn and editor.in_file(*drawn) else None
        unsaved = drawn and not self.in_file
        if self.old:
            self.setWindowTitle(f"Primer {name}")
        elif drawn:
            self.setWindowTitle(f"Primer {name} ({'in this file only' if self.in_file else 'not saved yet'})")
        else:
            self.setWindowTitle("New Primer")
        self.setMinimumWidth(640)
        self.name = QLineEdit(self.old["name"] if self.old else drawn[0] if drawn else self._free_name())
        self.name.setMaxLength(max(MAX_NAME, len(self.name.text())))   # not shorter than a name already in the database
        start = self.old["sequence"] if self.old else drawn[1] if drawn else editor.selected_sequence()
        if len(start) > MAX_PRIMER:   # a long selection is not a primer: leave the box empty rather than cut it
            start = ""
        self.sequence = dna_field(QLineEdit(start))
        self.sequence.setMaxLength(max(MAX_PRIMER, len(start)))
        self.sequence.setPlaceholderText(f"5'→3'; a lowercase 5' tail is fine (up to {MAX_PRIMER} letters)")
        self.note = QLineEdit(self.old["note"] if self.old else "")
        self.note.setMaxLength(max(MAX_NOTE, len(self.note.text())))
        self.reverse = QCheckBox("Reverse complement")
        self.reverse.toggled.connect(lambda: self.sequence.setText(reverse_complement(self.sequence.text())))
        if self.in_file:
            self.note.setEnabled(False)
            self.note.setToolTip("A primer of this file only has no note (the Primer Database keeps notes)")
        if unsaved:
            for box in (self.name, self.sequence):
                box.setReadOnly(True)
                box.setToolTip("Save puts this primer in the Primer Database; rename or change it afterwards")
            self.reverse.setEnabled(False)
        self.info = QLabel()
        # Every check with a mark: ✓ fine, ? worth a look, ✗ a problem.
        self.report = ReportView(10)
        self.report.setToolTip("Hairpin / dimer: ΔG at 37 °C of the most stable perfectly paired structure\n"
                               "(no internal loops or mismatches: a quick check, not a full fold).\n"
                               "Non-specific: ✗ the whole primer also binds there, ? only its 3′ end does.")
        self.timer = QTimer(self, singleShot=True, interval=250, timeout=self._update_report)
        # The genome check runs on a worker thread (every chromosome is read in turn); a newer sequence supersedes it.
        self._token, self._progress, self._error = 0, "", ""
        self._signals = _GenomeSignals(self)
        self._signals.progress.connect(self._genome_progress)
        self._signals.done.connect(self._genome_done)
        self._signals.failed.connect(self._genome_failed)
        self.sequence.textChanged.connect(self._changed)
        self._changed()
        buttons = QDialogButtonBox(QDialogButtonBox.Save | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        form = QFormLayout(self)
        form.addRow("Name", self.name)
        form.addRow("Sequence", self.sequence)
        form.addRow("", self.reverse)
        form.addRow("Note", self.note)
        form.addRow("Length", self.info)
        form.addRow("Checks", self.report)
        form.addRow(buttons)
        self.setTabOrder(self.sequence, self.reverse)   # the order of the rows, not of creation
        self.setTabOrder(self.reverse, self.note)
        self._update_report()

    def _free_name(self):
        taken = {r["name"] for r in db().rows}
        i = 1
        while f"primer{i}" in taken:
            i += 1
        return f"primer{i}"

    def seq(self):
        return "".join(self.sequence.text().split())

    def _changed(self):
        seq = self.seq()
        t = tm_text(seq)
        self.info.setText(f"{len(seq)} nt, Tm {t} °C (whole primer)" if t else f"{len(seq)} nt")
        self.timer.start()  # re-check once typing pauses

    def _update_report(self):
        seq = self.seq()
        self._token += 1   # stops a search still running for an older sequence
        token = self._token
        self._progress, self._error = "", ""
        self._render()
        if len(seq) < 8 or genome_query(self.editor, seq) is None or genome_cached(self.editor, seq) is not None:
            return
        editor, signals = self.editor, self._signals

        def work():
            try:
                genome_search(editor, seq, lambda done, total, name, n, count: signals.progress.emit(
                    token, f"searching {name} ({n} of {count}), {100 * done // max(total, 1)} %"), lambda: token != self._token)
                signal, args = signals.done, (token,)
            except Exception as e:   # said in the check (it stayed on "searching…" for ever), never lost in a thread
                signal, args = signals.failed, (token, str(e) or type(e).__name__)
            try:
                signal.emit(*args)
            except RuntimeError:
                pass   # the dialog is gone
        threading.Thread(target=work, daemon=True).start()

    def _render(self):
        checks = primer_checks(self.editor, self.seq(), "cached")
        say = f"the genome search failed: {self._error}" if self._error else self._progress
        if say:
            checks = [(n, lv, say) if n == "Non-specific (genome)" and lv == "info" else (n, lv, t) for n, lv, t in checks]
        self.report.setText(checks_html(checks))

    def _genome_progress(self, token, text):
        if token == self._token:
            self._progress = text
            self._render()

    def _genome_done(self, token):
        if token == self._token:
            self._progress = ""
            self._render()

    def _genome_failed(self, token, error):
        if token == self._token:
            self._progress, self._error = "", error
            self._render()

    def done(self, result):
        self._token += 1   # cancels a running genome search
        super().done(result)

    def accept(self):
        name, seq, note = self.name.text().strip(), self.seq(), self.note.text().strip()
        problem = primer_problem(name, seq, note)
        if problem:
            QMessageBox.warning(self, "Invalid primer", problem)
            (self.name if not name or any(c.isspace() and c != " " for c in name) else self.sequence).setFocus()
            return
        if self.in_file:
            if (name, seq) != self.in_file:
                if any(p[0] == name and p != self.in_file for p in self.editor.primer_list()):
                    QMessageBox.warning(self, "Name in use", f"There is already a primer '{name}' here.")
                    return
                self.editor.change_file_primer(self.in_file, (name, seq))
            super().accept()
            return
        ok = update_primer(self.editor, self.old["name"], name, seq, note) if self.old else add_primer(self.editor, name, seq, note) in ("added", "exists")
        if not ok:
            QMessageBox.warning(self, "Name in use", f"The primer database already has a primer '{name}'.")
            return
        super().accept()


def open_options(parent):
    """Tools > Settings on the Primers page; True when OK was pressed (the primer lists then redo their matching)."""
    return SettingsDialog(parent, "Primers").exec() == QDialog.Accepted


def import_file(parent, editor):
    """Import a primer list into the database (the New / Import of the Primer Database tab). Returns the report, or None."""
    path, _ = QFileDialog.getOpenFileName(parent, "Import primers", "", "Primer lists (*.tsv *.csv *.txt *.fa *.fasta);;All files (*)")
    if not path:
        return None
    try:
        rejected = []
        report = import_primers(path, editor, rejected)
        (QMessageBox.warning if rejected else QMessageBox.information)(parent, "Import primers", report)
        return report
    except (OSError, UnicodeDecodeError, ValueError) as e:
        QMessageBox.critical(parent, "Cannot read file", f"{path}\n\n{e}")
        return None


def export_file(parent):
    path, _ = QFileDialog.getSaveFileName(parent, "Export primer database", "primers.tsv",
                                          "Tab-separated (*.tsv);;CSV (*.csv);;FASTA (*.fasta)")
    if path:
        try:
            db().export(path)
        except OSError as e:
            QMessageBox.critical(parent, "Export failed", f"{path}\n\n{e}")


class OrderDialog(QDialog):
    """Export for Ordering: the file format, and for a plate the synthesis scale and purification."""

    def __init__(self, parent, count):
        super().__init__(parent)
        self.setWindowTitle("Export primers for ordering")
        self.format = QComboBox()
        self.format.addItems(ORDER_FORMATS)
        self.scale = QComboBox()
        self.scale.addItems(SCALES)
        self.purification = QComboBox()
        self.purification.addItems(PURIFICATIONS)
        self.format.currentIndexChanged.connect(self._format_changed)
        form = QFormLayout(self)
        form.addRow("Format", self.format)
        form.addRow("Scale", self.scale)
        form.addRow("Purification", self.purification)
        plates = -(-count // PLATE_SIZE)
        form.addRow(QLabel(f"{count} {'primer' if count == 1 else 'primers'}" + (f", on {plates} plates: one file per plate" if plates > 1 else ""),
                           wordWrap=True))
        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.button(QDialogButtonBox.Ok).setText("Export…")
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        form.addRow(buttons)

    def _format_changed(self):
        plate = self.format.currentIndex() == 0   # the plain two-column file has no scale or purification
        self.scale.setEnabled(plate)
        self.purification.setEnabled(plate)


def export_for_ordering(parent, names=None):
    """Export for Ordering…: the primers called `names` (the ones selected in the Primer Database), all of them when none."""
    rows = [r for r in db().rows if not names or r["name"] in set(names)]
    if not rows:
        QMessageBox.information(parent, "Export for ordering", "The primer database is empty: there is nothing to export.")
        return []
    dialog = OrderDialog(parent, len(rows))
    if dialog.exec() != QDialog.Accepted:
        return []
    path, _ = QFileDialog.getSaveFileName(parent, "Export primers for ordering", "primer_order.csv", "CSV (*.csv)")
    if not path:
        return []
    try:
        written = export_order(rows, path, dialog.format.currentText(), dialog.scale.currentText(), dialog.purification.currentText())
    except (OSError, ValueError) as e:
        QMessageBox.critical(parent, "Export failed", f"{path}\n\n{e}")
        return []
    if len(written) > 1:
        QMessageBox.information(parent, "Export for ordering", "Written, one file per plate:\n" + "\n".join(str(w) for w in written))
    return written


MISMATCHES = Qt.UserRole + 2   # on a site's sequence cell: the indices of its bases that do not pair (drawn red)
# The rows of the background genome under a primer (its sites, and a note: searching…, no site…), the only ones redone alone when
# a primer is opened or its search ends; they come after the primer is opened, so they never set a column's width.
GENOME_KIND, GENOME_NOTE = "genome", "genome note"
GROUP_ROLE = Qt.UserRole + 3   # GENOME_KIND or GENOME_NOTE on those rows
SITE_ROLE = Qt.UserRole + 4    # on a primer with one site (no row of its own under it): what a click on that site's row would select
MISMATCH_COLOR, MISMATCH_COLOR_SELECTED = "#d32f2f", "#ffcdd2"
_PAIRS = {"A": "A", "C": "C", "G": "G", "T": "T", "U": "T", "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
          "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}


def mismatch_indices(primer, template):
    """Indices of `primer` (5'→3') whose base does not match `template` (the same strand, as long): an ambiguity code matches
    any of its bases."""
    return [k for k, (p, t) in enumerate(zip(primer.upper(), template.upper())) if t not in _PAIRS.get(p, p)]


def primer_item(item):
    """The primer's row of any row under it (a group of sites, a site)."""
    while item.parent() is not None:
        item = item.parent()
    return item


def _below(item):
    """Every row under `item`, at any depth."""
    out = []
    for i in range(item.childCount()):
        out.append(item.child(i))
        out += _below(item.child(i))
    return out


class _PanelSignals(QObject):
    done = Signal(str, str)   # a genome search finished: the primer's sequence (upper case), and what went wrong or ""


class PrimerTree(QTreeWidget):
    """The primers' tree. Its indent is small (a site sits just under its primer), so its ▸ / ▾ are drawn here at a size of their
    own, and a click on one, or just beside it, opens or closes the primer instead of starting a selection. Each column is as wide as the primers and their sites on this file, opened or not, so opening or
    closing a primer never changes a width; the background genome's rows (added when a primer is opened) do not count: a long
    one is cut, its tooltip says it whole."""

    INDENT = 10     # px per level (Qt's 20 pushed the sites far right)
    ARROW = 8       # the ▸ / ▾, px
    SLACK = 6       # a click this far right of the arrow still opens or closes the primer

    def _opens(self, item):
        return item is not None and (item.childCount() > 0 or item.childIndicatorPolicy() == QTreeWidgetItem.ShowIndicator)

    def drawBranches(self, painter, rect, index):
        item = self.itemFromIndex(index)
        if not self._opens(item):
            return
        c = QPointF(rect.right() - self.indentation() / 2 + 1, rect.center().y())
        h = self.ARROW / 2
        points = ([QPointF(c.x() - h, c.y() - h / 2), QPointF(c.x() + h, c.y() - h / 2), QPointF(c.x(), c.y() + h / 2 + 1)]
                  if item.isExpanded() else
                  [QPointF(c.x() - h / 2, c.y() - h), QPointF(c.x() + h / 2 + 1, c.y()), QPointF(c.x() - h / 2, c.y() + h)])
        painter.save()
        painter.setRenderHint(QPainter.Antialiasing)
        painter.setPen(Qt.NoPen)
        painter.setBrush(self.palette().color(QPalette.Text))
        painter.drawPolygon(QPolygonF(points))
        painter.restore()

    def mousePressEvent(self, event):
        pos = event.position().toPoint()
        index = self.indexAt(pos)
        item = self.itemFromIndex(index) if index.isValid() and index.column() == 0 else None
        if event.button() == Qt.LeftButton and self._opens(item):
            left = self.visualRect(index).left()
            if left - self.indentation() <= pos.x() <= left + self.SLACK:
                item.setExpanded(not item.isExpanded())
                event.accept()
                return
        super().mousePressEvent(event)

    def mouseDoubleClickEvent(self, event):   # (a quick second click on the arrow: opened and closed again, not an edit)
        index = self.indexAt(event.position().toPoint())
        item = self.itemFromIndex(index) if index.isValid() and index.column() == 0 else None
        if self._opens(item) and event.position().x() <= self.visualRect(index).left() + self.SLACK:
            return self.mousePressEvent(event)
        super().mouseDoubleClickEvent(event)

    def sizeHintForColumn(self, column):
        option = QStyleOptionViewItem()
        self.initViewItemOption(option)
        widest = 0
        for i in range(self.topLevelItemCount()):
            top = self.topLevelItem(i)
            for depth, item in [(0, top)] + [(1, top.child(k)) for k in range(top.childCount())]:
                if item.data(0, GROUP_ROLE):
                    continue
                index = self.indexFromItem(item, column)
                width = self.itemDelegateForIndex(index).sizeHint(option, index).width()
                if column == 0:
                    width += (depth + self.rootIsDecorated()) * self.indentation()
                widest = max(widest, width)
        return widest


class MismatchDelegate(QStyledItemDelegate):
    """The sequence of a binding site, letter by letter: the bases that do not pair with the template in red, bold. Every letter
    takes the room it has in the normal font (a bold one centred in it), so the bases after a mismatch stay in line with the
    primer's sequence on the row above and with the exact sites."""

    def paint(self, painter, option, index):
        marks = set(index.data(MISMATCHES) or [])
        if not marks:
            return super().paint(painter, option, index)
        opt = QStyleOptionViewItem(option)
        self.initStyleOption(opt, index)
        text, opt.text = opt.text, ""
        style = opt.widget.style() if opt.widget is not None else QApplication.style()
        style.drawControl(QStyle.CE_ItemViewItem, opt, painter, opt.widget)   # background, selection, focus
        rect = style.subElementRect(QStyle.SE_ItemViewItemText, opt, opt.widget).adjusted(3, 0, -3, 0)
        selected = bool(opt.state & QStyle.State_Selected)
        plain = opt.palette.color(QPalette.HighlightedText if selected else QPalette.Text)
        red = QColor(MISMATCH_COLOR_SELECTED if selected else MISMATCH_COLOR)
        bold = QFont(opt.font)
        bold.setBold(True)
        painter.save()
        painter.setClipRect(rect)
        metrics, bold_metrics = QFontMetrics(opt.font), QFontMetrics(bold)
        y = rect.top() + (rect.height() + metrics.ascent() - metrics.descent()) // 2
        x = rect.left()
        for k, letter in enumerate(text):
            mismatch = k in marks
            width = metrics.horizontalAdvance(letter)
            painter.setFont(bold if mismatch else opt.font)
            painter.setPen(red if mismatch else plain)
            shift = (width - bold_metrics.horizontalAdvance(letter)) / 2 if mismatch else 0
            painter.drawText(QPointF(x + shift, y), letter)
            x += width
        painter.restore()


class PrimerPanel(QWidget):
    """Primers that bind the open sequence, one group per primer (collapsed by default),
    one child row per binding site.

    The Primer Database tab (db_tabs.py) uses this same view on a document made of every primer of the database: its editor's
    `database` lists every primer (with its note and the day it was added), and gives the buttons and right-click menu of the
    database instead of those of a file."""
    COLUMNS = ["Name", "Sequence (5'→3')", "Length", "Tm °C", "Match", "Structure"]
    DATABASE_COLUMNS = ["In file", "Added", "Note"]   # the note last: it takes the room left, the date stays in view
    FILTERED = ["Match", "Structure", "Added"]   # filter by values (value_filter.py)

    def __init__(self, editor):
        super().__init__()
        self.editor = editor
        self.database = getattr(editor, "database", None)
        # A file's list ticks what is drawn, a whole primer or one site of it, in a Visible column shown first (added last, so
        # every other column keeps its index)
        self.columns = self.COLUMNS + (self.DATABASE_COLUMNS if self.database else ["Visible"])
        self.visible = None if self.database else len(self.columns) - 1
        # The database's list: a box ticked for each primer that is part of the file open in the window (tick: import it)
        self.in_file = self.columns.index("In file") if self.database else None
        self.tree = PrimerTree()
        self.tree.setHeaderLabels(self.columns)
        if self.visible is not None:   # an eye, not a word: the column is only as wide as its checkboxes
            self.tree.headerItem().setText(self.visible, "")
            self.tree.header().moveSection(self.visible, 0)
            self._paint_eye()
        if self.in_file is not None:   # first, like a file's eye
            self.tree.header().moveSection(self.in_file, 0)
        self.tree.setUniformRowHeights(True)
        # (before fit_columns: a tree's header stretches its last column by default, and that width would count as one the user
        # dragged, never fitted again)
        self.tree.header().setStretchLastSection(False)
        # Columns as wide as what they hold; the sequence is cut at 24 bases (drag the column wider to read more) so this panel
        # is no wider than the Features table in the same pane: switching between the two keeps the pane's width. The note
        # (the database), else the structure, takes the room left.
        fitter = fit_columns(self.tree, stretch=self.columns.index("Note" if self.database else "Structure"),
                             cap=self.tree.fontMetrics().horizontalAdvance("G" * 24) + 24)
        # Visible, Length, Tm: as wide as a checkbox or a number (their header is cut; its tooltip names it)
        fitter.tight = frozenset(c for c in (self.visible, 2, 3) if c is not None)
        for c in fitter.tight:
            self.tree.headerItem().setToolTip(c, self.columns[c])
        self.tree.setIndentation(PrimerTree.INDENT)   # a site sits just under its primer
        self.tree.setToolTip("Tm: of the bases that anneal at that site (5' tail excluded)")
        self.tree.setSelectionMode(QAbstractItemView.ExtendedSelection)   # several primers: Import to This Sequence for all
        self.tree.currentItemChanged.connect(self._select)
        self.tree.itemDoubleClicked.connect(lambda item, _: self.edit(item))
        self.tree.itemChanged.connect(self._ticked)
        self.tree.setContextMenuPolicy(Qt.CustomContextMenu)
        self.tree.customContextMenuRequested.connect(self._menu)
        self.expanded = set()  # primer names the user opened; kept across refreshes
        self.tree.itemExpanded.connect(self._expanded)
        self.tree.itemCollapsed.connect(lambda item: item.parent() is None and self.expanded.discard(item.data(0, Qt.UserRole)))
        self.tree.setItemDelegateForColumn(1, MismatchDelegate(self.tree))
        self._searching = set()   # primers whose genome search is running
        self._waiting = []        # primers (their sequence) whose genome search waits its turn (one runs at a time)
        self._failed = {}         # primer sequence (upper case) -> why its genome search failed
        self._genome_signals = _PanelSignals(self)
        self._genome_signals.done.connect(self._genome_found)
        self._reload_later = QTimer(self, singleShot=True, interval=0, timeout=self.reload)

        # Click a header to sort by it (again: the other way). Until then the primers are in the order they bind along the sequence.
        self.sort = None   # (column, order) chosen by a click, kept across refreshes
        header = self.tree.header()
        header.setSectionsClickable(True)
        header.setSortIndicatorShown(False)
        header.sectionClicked.connect(self._sort_by)
        import value_filter
        value_filter.install(self.tree, self.FILTERED)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        self.count = QLabel()
        self.count.setEnabled(False)   # greyed: a detail
        self.count.setMinimumWidth(1)   # cut rather than widen the pane (see minimumSizeHint)
        self.edit_button = None
        if not self.database:   # the Primer Database tab has these buttons (and its own) in its header
            # The same header as the database tabs': the count on the left, the actions as one group on the right.
            head = QHBoxLayout()
            head.setContentsMargins(4, 2, 4, 2)
            head.addWidget(self.count)
            head.addStretch()
            for text, tip, slot in (("New…", "A new primer (from the selection, if any)", lambda: PrimerDialog(self, editor).exec()),
                                    ("Edit…", "Edit the primer clicked", lambda: self.edit(self.tree.currentItem())),
                                    ("Database…", "Databases → Primer Database, as a tab", self.open_database),
                                    ("Options…", "PCR conditions and the primer checks", self.edit_options)):
                button = QPushButton(text, clicked=slot, toolTip=tip)
                head.addWidget(button)
                if text == "Edit…":
                    self.edit_button = button
            layout.addLayout(head)
            self.tree.itemSelectionChanged.connect(self._can_edit)
        layout.addWidget(self.tree)
        editor.primers_changed.connect(self.reload)
        editor.elsewhere_changed.connect(self.reload)
        self.sources = self._sources()
        editor.display_changed.connect(self._display_changed)
        self.reload()

    def minimumSizeHint(self):
        """As narrow as the Features view (its table's minimum): switching a pane between the two never widens it; the buttons
        above the list are squeezed instead."""
        return QSize(self.tree.minimumSizeHint().width(), super().minimumSizeHint().height())

    def _sources(self):
        """The primer sources switched on in the Display filter: (From this sequence, From the database)."""
        return self.editor.show_file_primers, self.editor.show_db_primers

    def _display_changed(self):
        """The Display filter changed: the list follows its primer switches (only those: the other switches change nothing here)."""
        if not self.database and self._sources() != self.sources:
            self.sources = self._sources()
            self.reload()

    def _can_edit(self):
        if self.edit_button is not None:
            self.edit_button.setEnabled(self.tree.currentItem() is not None and bool(self.tree.selectedItems()))

    def _paint_eye(self):
        import icons
        self.tree.headerItem().setIcon(self.visible, icons.eye_icon(self.palette()))

    def theme_changed(self):
        """Light / dark switch (theme.refresh_views): the greyed (hidden) and temporary primers, and the eye, take the new colours."""
        if self.visible is not None:
            self._paint_eye()
        self.reload()

    def reload(self):
        ed = self.editor
        groups, elsewhere = {}, {}
        lengths = dict(ed.fasta_records)
        for site in ed.primer_sites:
            groups.setdefault(site.label, []).append(site)
        for site in ed.elsewhere:   # a multi-record file: the sites on the other records are listed too
            elsewhere.setdefault(site.label, []).append(site)
            groups.setdefault(site.label, [])
        rows, filtered = {}, 0
        here = self.database.file_editor() if self.database else None   # the file open in the window (the database's Match)
        genome_on = genome() is not None
        cut = ed.cut_primers()   # primers with more sites than are listed (the first MAX_SITES of each sequence are)
        if self.database:   # every primer of the database, also one that binds nowhere (an ambiguity code...)
            # Match: where each binds the file open in this window (not the list's own document); those that bind it first.
            rows = {r["name"]: r for r in db().rows}
            on = {}
            for s in here.primer_sites if here is not None else ():
                if s.label in rows and s.sequence.upper() == rows[s.label]["sequence"].upper():
                    on.setdefault(s.label, []).append(s)
            groups = {name: on.get(name, []) for name in sorted(rows, key=lambda n: n not in on)}
            cut = here.cut_primers() if here is not None else {}
            self.tree.headerItem().setToolTip(self.columns.index("Match"), f"Where the primer binds {here.path.name}, the file open in this window"
                                              if here is not None else "Open a file in this window to see where the primers bind it")
        elif not all(self._sources()):
            # Only the primers of the sources switched on in the Display filter (From this sequence / From the database), as
            # the views draw them; a hidden primer of a source shown stays listed, greyed.
            in_db = ed._db_primer_keys()
            shown = {name: sites for name, sites in groups.items()
                     if ed.primer_source_shown(name, (sites or elsewhere[name])[0].sequence, in_db)}
            filtered, groups = len(groups) - len(shown), shown
        self.tree.blockSignals(True)  # don't select sites while refilling
        current = self.tree.currentItem()
        keep = primer_item(current).data(0, Qt.UserRole) if current is not None else None
        self.tree.clear()
        for name, sites in groups.items():
            away = elsewhere.get(name, [])
            seq = (sites or away)[0].sequence if sites or away else rows[name]["sequence"]
            exact = sum(s.mismatches == 0 for s in sites + away)
            t = tm_text(seq)
            total = len(sites) + len(away)
            summary = f"{total} site{'s' * (total > 1)}" + (f", {exact} exact" if exact != total else "") + \
                (f" ({len(away)} on other sequences)" if away else "")
            if name in cut:
                summary = f"{cut[name]:,} sites, the first {MAX_SITES} on each sequence listed"
            warnings = [] if is_degenerate(seq) else structure(seq)["warnings"]
            mark, color = MARKS["warn" if warnings else "ok"]
            group = SortItem([name, seq, str(len(seq)), t or "–", summary if total else "" if self.database else "no site",
                              f"{mark} " + "; ".join(warnings) if warnings else mark]
                             + (["", rows[name]["added"], rows[name]["note"]] if self.database else []))
            group.keys = {0: name.lower(), 2: len(seq), 3: tm_sort_key(seq), 4: total}
            group.setForeground(5, QColor(color))
            group.setData(0, Qt.UserRole, name)
            if sites and all(s.provisional for s in sites) and not self.database:   # drawn in the temporary colour: not in this file
                group.setForeground(0, QColor(TEMP_COLOR))
                group.setToolTip(0, "Not part of this file: from the primer database (or designed, not saved). "
                                    "Right-click > Import to This Sequence keeps it in the file.")
            if here is not None:
                self._in_file_box(group, here, name, seq, bool(sites))
            if self.visible is not None:   # ticked: drawn in this file (set before the row is in the tree: no signal)
                hidden = name in ed.hidden_primers
                group.setFlags(group.flags() | Qt.ItemIsUserCheckable)
                group.setCheckState(self.visible, Qt.Unchecked if hidden else Qt.Checked)
                group.setToolTip(self.visible, "Untick to stop drawing this primer in this file; it stays in the primer database")
                if hidden:   # still listed, greyed
                    for col in range(len(self.columns)):
                        group.setForeground(col, self.palette().color(QPalette.PlaceholderText))
            self._add_sites(group, seq, sites, away, lengths)
            if not group.childCount() and genome_on:
                self._genome_arrow(group, seq)
            self._align(group)
            self.tree.addTopLevelItem(group)
            group.setExpanded(name in self.expanded)  # collapsed unless the user opened it
            self._span_notes(group)
            if name == keep:   # the primer picked stays picked (its row only: the sites may have moved)
                self.tree.setCurrentItem(group)
        self._apply_sort()
        self._decorate()
        self.tree.blockSignals(False)
        n = len(groups)
        text = plural(n, "primer") if self.database else f"{plural(n, 'primer')} on this sequence"
        if filtered:
            text += f" ({filtered} more not shown: Display filter)"
        self.count.setText(text)
        self.count.setToolTip(text)   # (cut in a narrow pane)
        self._can_edit()

    # ---- the sites of a primer: Match says where (Sequence: this file; Background: the background genome) -----------------------
    def _add_sites(self, group, seq, sites, away, lengths):
        """Under a primer: its sites on this sequence, on the other sequences of the file (their name first), then in the
        background genome; in each site's sequence the bases that do not pair with the template are red (MismatchDelegate)."""
        ed = self.editor
        if sites and not self.database:   # (the database's own list document is not a sequence of yours)
            for site in sites:
                bound = seq[-sum(e - s for s, e in site.spans):]
                template = "".join(ed.doc.sequence(a, b) for a, b in site.spans)
                row = self._site_row(location_text(site.start, site.end, site.strand, ed.length), bound,
                                     template if site.strand > 0 else reverse_complement(template), ("", site.start), site.key)
                row.setFlags(row.flags() | Qt.ItemIsUserCheckable)   # ticked: drawn (set before the row is in the tree: no signal)
                row.setCheckState(self.visible, Qt.Checked if ed.site_shown(site) else Qt.Unchecked)
                row.setToolTip(self.visible, "Untick to stop drawing this site (the primer's other sites stay)"
                               + (" · best site" if site.best else ""))
                group.addChild(row)
        for site in away:
            bases = getattr(site, "bases", None)
            group.addChild(self._site_row(
                f"{site.record}: " + location_text(site.start, site.end, site.strand, lengths.get(site.record, site.end)),
                site.sequence[-(site.end - site.start):], bases if site.strand > 0 or bases is None else reverse_complement(bases),
                (site.record, site.start), (site.record, site.start, site.end)))   # a click opens that record
        ticked = lambda row: not row.flags() & Qt.ItemIsUserCheckable or row.checkState(self.visible) == Qt.Checked
        if group.childCount() == 1 and ticked(group.child(0)):
            self._fold(group, len(seq))
        if group.data(0, Qt.UserRole) in self.expanded:   # the genome is searched only for a primer opened
            self._genome_rows(group, seq)

    def _in_file_box(self, group, here, name, seq, binds):
        """The database's In file box of a primer: ticked when it is part of `here`'s file. Only a primer with a site on it has
        one (Save writes a primer as its sites)."""
        inside = here.in_file(name, seq)
        group.keys = {**group.keys, self.in_file: not inside}   # (sorted: those in the file first)
        if not binds:   # (no check state at all: Qt draws a box for any row that has one)
            return
        group.setFlags(group.flags() | Qt.ItemIsUserCheckable)
        group.setCheckState(self.in_file, Qt.Checked if inside else Qt.Unchecked)
        group.setToolTip(self.in_file, f"Part of {here.path.name}: untick to take it out (it stays in the database)" if inside
                         else f"Tick to import it to {here.path.name}")

    @staticmethod
    def _fold(group, length):
        """A primer with one site is one row: the site's match on it (its mismatches red in the primer's sequence), its place in
        the name's tooltip, a click selects it. (A site unticked by hand keeps its row: its box ticks it again.)"""
        row = group.takeChild(0)
        group.setText(4, row.text(4))
        if row.data(1, MISMATCHES):   # (its mismatches lower the Tm of the site)
            group.setText(3, row.text(3))
            group.keys = {**group.keys, 3: row.keys[3]}
        group.setData(1, MISMATCHES, [k + length - len(row.text(1)) for k in row.data(1, MISMATCHES) or []])
        group.setData(0, SITE_ROLE, row.data(0, Qt.UserRole))
        group.setToolTip(0, "\n".join(x for x in (row.text(0), group.toolTip(0)) if x))

    def _genome_rows(self, group, seq):
        q = genome_query(self.editor, seq)
        if q is None:
            return
        result = genome_cached(self.editor, seq)
        if result is None and seq.upper() in self._failed:
            return self._note(group, f"{q[0].name}: the search failed ({self._failed[seq.upper()]})")
        if result is None:
            self._search_genome(seq, first=True)
            return self._note(group, f"{q[0].name}: searching…")
        listed = 0
        for whole, hits in ((True, result.whole), (False, result.hits)):
            for h in hits:
                name, start, end, strand, _, length, circular = h
                bases = result.bases.get((name, start, end, strand)) if hasattr(result, "bases") else None
                where = site_in_record(h)
                row = self._site_row(where, seq[-(end - start):],
                                     None if bases is None else bases if strand > 0 else reverse_complement(bases),
                                     (name, start), None, "whole primer" if whole else "3′ end", "Background")
                row.setData(0, GROUP_ROLE, GENOME_KIND)
                row.setToolTip(0, where)   # (its column may be narrower: these rows do not set its width)
                group.addChild(row)
                listed += 1
        more = result.n_whole + result.n_hits - listed
        if not listed:
            self._note(group, f"{q[0].name}: no site" + ("" if result.complete else " in the part searched"))
        elif more > 0 or not result.complete:
            self._note(group, f"{q[0].name}: " + "; ".join(
                x for x in (f"{more} more not listed" if more > 0 else "", "" if result.complete else "not searched to the end") if x))

    def _site_row(self, where, bound, template, sort_key, data, kind="", source="Sequence"):
        """One binding site: where, the primer's annealing bases (mismatches red when the template is known), length, Tm, and
        Match: the map it is on (Sequence, Background), how it binds there."""
        marks = mismatch_indices(bound, template) if template is not None and len(template) == len(bound) else None
        tb, t = site_tm_text(bound, template)
        match = "" if marks is None else "exact" if not marks else plural(len(marks), 'mismatch')
        row = SortItem([where, bound, str(len(bound)), tb or "–", ", ".join(x for x in (source, kind, match) if x)])
        row.keys = {0: sort_key, 2: len(bound), 3: tm_sort_key(bound) if t is None else t, 4: (source, len(marks or []))}
        row.setData(0, Qt.UserRole, data)
        row.setData(1, MISMATCHES, marks or [])
        if marks:
            row.setToolTip(1, "Red: the primer's base does not pair with the template there")
        return row

    def _note(self, parent, label):
        """A line about the background genome's search under a primer (searching…, no site…), in italics across the row."""
        g = QTreeWidgetItem([label])
        f = g.font(0)
        f.setItalic(True)
        g.setFont(0, f)
        g.setData(0, GROUP_ROLE, GENOME_NOTE)
        parent.addChild(g)

    def _genome_arrow(self, group, seq):
        """A primer with no row under it (one site, folded) opens only to show its sites in the background genome: its ▸ is there
        only once the genome is known to have some (or the search failed: opened, it says why). The search runs in the
        background; the Primer Database's list (every primer) searches only the primer opened, so its ▸ is always there."""
        if self.database:
            return group.setChildIndicatorPolicy(QTreeWidgetItem.ShowIndicator)
        result = genome_cached(self.editor, seq)
        if result is None and seq.upper() not in self._failed:
            self._search_genome(seq)
        found = seq.upper() in self._failed or (result is not None and result.n_whole + result.n_hits > 0)
        group.setChildIndicatorPolicy(QTreeWidgetItem.ShowIndicator if found else QTreeWidgetItem.DontShowIndicator)

    def _decorate(self):
        """No row opens: no column of empty room for the ▸ on every row."""
        items = [self.tree.topLevelItem(i) for i in range(self.tree.topLevelItemCount())]
        self.tree.setRootIsDecorated(any(i.childCount() or i.childIndicatorPolicy() == QTreeWidgetItem.ShowIndicator for i in items))

    def _search_genome(self, seq, first=False):
        """Search the background genome for `seq` on a worker thread, then show the sites. One search runs at a time: the others
        wait their turn, a primer just opened (`first`) ahead of the rest."""
        key = seq.upper()
        if key in self._searching:
            return
        if self._searching:
            self._waiting = [w for w in self._waiting if w.upper() != key]   # (as typed: a lowercase 5' tail is not searched)
            self._waiting.insert(0, seq) if first else self._waiting.append(seq)
            return
        self._searching.add(key)
        editor, signals = self.editor, self._genome_signals

        def work():
            try:
                genome_search(editor, seq)
                error = ""
            except Exception as e:   # said in the group (it stayed on "searching…"), never lost in a thread
                error = str(e) or type(e).__name__
            try:
                signals.done.emit(seq.upper(), error)
            except RuntimeError:
                pass   # the panel is gone
        threading.Thread(target=work, daemon=True).start()

    def _genome_found(self, key, error):
        self._searching.discard(key)
        if error:
            self._failed[key] = error   # not searched again at each refresh: the group says why
        for i in range(self.tree.topLevelItemCount()):   # the primers that have this sequence: only their genome rows change
            item = self.tree.topLevelItem(i)
            if item.text(1).upper() != key:
                continue
            if item.data(0, Qt.UserRole) in self.expanded:
                self._refresh_genome(item)
            elif not item.childCount():
                self._genome_arrow(item, item.text(1))
        self._decorate()
        while self._waiting:   # the next one waiting (unless it was found meanwhile)
            seq = self._waiting.pop(0)
            if genome_cached(self.editor, seq) is None and genome_query(self.editor, seq) is not None:
                self._search_genome(seq)
                break

    def _expanded(self, item):
        if item.parent() is None:
            name = item.data(0, Qt.UserRole)
            if name not in self.expanded:
                self.expanded.add(name)
                # Its sites in the reference genome: searched now. Only that group is added to the row (the list is not rebuilt: it
                # would scroll back to the top), after this event, not inside it.
                QTimer.singleShot(0, lambda: self._refresh_genome_of(name))

    def _refresh_genome_of(self, name):
        for i in range(self.tree.topLevelItemCount()):
            if self.tree.topLevelItem(i).data(0, Qt.UserRole) == name:
                self._refresh_genome(self.tree.topLevelItem(i))
                return

    def _refresh_genome(self, item):
        """Redo the genome rows under the primer `item` alone, the rest of the list (and its scroll position) untouched."""
        self.tree.blockSignals(True)   # no site gets selected by the change
        for i in reversed(range(item.childCount())):
            if item.child(i).data(0, GROUP_ROLE):
                item.removeChild(item.child(i))
        self._genome_rows(item, item.text(1))
        self._align(item)
        self._span_notes(item)
        self.tree.setRootIsDecorated(True)
        self.tree.blockSignals(False)

    @staticmethod
    def _align(group):
        """Numbers (length, Tm) on the right, in every row of a primer, the genome's added later too."""
        for item in [group] + _below(group):
            for col in (2, 3):
                item.setTextAlignment(col, Qt.AlignRight | Qt.AlignVCenter)

    def _span_notes(self, group):
        for i in range(group.childCount()):   # a note's text across the row (only once it is in the tree)
            if group.child(i).data(0, GROUP_ROLE) == GENOME_NOTE:
                group.child(i).setFirstColumnSpanned(True)

    def _sort_by(self, column):
        """A header clicked: sort by it, ascending; clicked again: descending."""
        order = Qt.DescendingOrder if self.sort == (column, Qt.AscendingOrder) else Qt.AscendingOrder
        self.sort = (column, order)
        self._apply_sort()

    def _apply_sort(self):
        if self.sort is not None:
            header = self.tree.header()
            header.setSortIndicatorShown(True)
            header.setSortIndicator(*self.sort)
            self.tree.sortItems(*self.sort)

    def _select(self, item, _previous=None):
        if item is None:
            return
        key = item.data(0, Qt.UserRole)
        if item.parent() is None:  # a primer: select its site, else its first site on this sequence
            key = item.data(0, SITE_ROLE) or next(
                (row.data(0, Qt.UserRole) for row in _below(item) if isinstance(row.data(0, Qt.UserRole), str)), None)
        elif item.data(0, GROUP_ROLE):
            return
        if isinstance(key, tuple):   # a site on another record: open it (after this click: the panel is replaced)
            ed = self.editor
            QTimer.singleShot(0, lambda: ed.record_hit_requested.emit(*key))
        elif key:
            self.editor.select_primer(key)

    def chosen(self):
        """The names of the primers selected (a site counts for its primer)."""
        return sorted({primer_item(i).data(0, Qt.UserRole) for i in self.tree.selectedItems()})

    def _menu(self, pos):
        item = self.tree.itemAt(pos)
        if item is not None and self.database:
            menu = QMenu(self)
            self.database.menu(menu, self, self.chosen() or [primer_item(item).data(0, Qt.UserRole)])
            self.tree.value_filter.cell_actions(menu, pos)   # show only / hide the value of the cell
            menu.exec(self.tree.viewport().mapToGlobal(pos))
        elif item is not None:   # the same menu as on a primer in the sequence view
            from feature_table import primer_actions, primers_actions  # here: feature_table imports the tools
            menu = QMenu(self)
            name = primer_item(item).data(0, Qt.UserRole)
            chosen = {primer_item(i).data(0, Qt.UserRole) for i in self.tree.selectedItems()}
            if len(chosen) > 1 and name in chosen:   # several primers selected: every action applies to all of them
                primers_actions(menu, self.editor, chosen)
            else:
                primer_actions(menu, self.editor, name, self.editor.site_by_key.get(item.data(0, Qt.UserRole)))
            self.tree.value_filter.cell_actions(menu, pos)   # show only / hide the value of the cell
            menu.exec(self.tree.viewport().mapToGlobal(pos))

    def _ticked(self, item, column):
        """A Visible checkbox: a primer's shows or hides it in this file, a site's that site only."""
        if self.in_file is not None and column == self.in_file and item.parent() is None:
            here, row = self.database.file_editor(), db().find(item.data(0, Qt.UserRole))
            on = item.checkState(column) == Qt.Checked
            if here is not None and row is not None and on != here.in_file(row["name"], row["sequence"]):
                QTimer.singleShot(0, lambda: here.import_primers([(row["name"], row["sequence"])], on))
            return
        if column != self.visible or self.visible is None:
            return
        ed, on = self.editor, item.checkState(column) == Qt.Checked
        # (after this signal: the change rebuilds the tree, and the clicked row must not be deleted while Qt still handles it)
        if item.parent() is None:
            name = item.data(0, Qt.UserRole)
            if on == (name in ed.hidden_primers):
                QTimer.singleShot(0, lambda: ed.set_primer_hidden(name, not on))
            return
        site = ed.site_by_key.get(item.data(0, Qt.UserRole))
        if site is not None and on != ed.site_shown(site):
            QTimer.singleShot(0, lambda: ed.set_sites_shown([site], on))

    def edit(self, item):
        if item is not None:
            top = primer_item(item)
            if self.database:   # the database's tab: checked against the primer on view
                self.database.edit(top.data(0, Qt.UserRole))
            else:
                PrimerDialog(self, self.editor, top.data(0, Qt.UserRole)).exec()

    def open_database(self):
        """Databases → Primer Database, as a tab of this window."""
        import db_tabs
        db_tabs.open_tab(page_of(self), "primers")

    def edit_options(self):
        if open_options(self):
            self._failed.clear()   # another genome, other rules: searched again
            self.reload()
            self.editor.selection_changed.emit()  # status bar Tm
