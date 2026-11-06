"""Primers: the Primers tab (sites on the open sequence, grouped by primer), the primer
dialog (type a sequence, see live where it binds), the primer database window, options.

A primer is just a name + sequence in the primer database (primer_db.py). Where it binds
is never stored: Editor.match_primers recomputes it on the current sequence after every
edit, including sites with a few mismatches ("almost binds", orange, the mismatched bases written on the line).
"""
from html import escape
import sys
import threading
from pathlib import Path

from PySide6.QtCore import QEventLoop, QObject, QPointF, Qt, QTimer, Signal
from PySide6.QtGui import QColor, QFont, QFontMetrics, QPalette
from PySide6.QtWidgets import (
    QAbstractItemView, QCheckBox, QDialog, QDialogButtonBox, QDoubleSpinBox, QFileDialog, QFrame,
    QApplication, QFormLayout, QHBoxLayout, QLabel, QLineEdit, QMenu, QMessageBox, QProgressDialog, QPushButton,
    QSpinBox, QStyle, QStyledItemDelegate, QStyleOptionViewItem, QTextBrowser, QTreeWidget, QTreeWidgetItem, QVBoxLayout, QWidget,
)

from opendna_core import dimer, hairpin, tm

import tool_settings
from editor import dna_field, reverse_complement
from value_filter import SortItem
from tool_settings import Param, SettingsDialog
from reference import Reference, as_reference
from sequence_view import TEMP_COLOR
from uikit import dispose_progress, fit_columns, plural, show_after
from primer_db import (MAX_NAME, MAX_NOTE, MAX_PRIMER, add_primer, add_primers, binding_mismatches, db, notifier,  # noqa: F401 (re-exported)
                       primer_problem, settings, update_primer)
import theme

# Where a demo genome may sit: the repository folder, or (the one-file exe unpacks itself into a temp folder) next to the exe.
ROOT = Path(sys.executable).resolve().parent if getattr(sys, "frozen", False) else Path(__file__).resolve().parent.parent
# key, label, default (Primer3's defaults), max
CONDITIONS = [("na_mm", "Na⁺/K⁺ (mM)", 50.0, 2000), ("mg_mm", "Mg²⁺ (mM)", 1.5, 100),
              ("dntp_mm", "dNTPs (mM)", 0.6, 100), ("primer_nm", "Primer (nM)", 50.0, 100000)]
MAX_LISTED = 8     # positions listed in a report before "... and N more"


# Smallest value that still gives a meaningful Tm: no sodium or no primer at all makes the salt / concentration terms
# undefined (the Tm came out as -273 °C). Mg and dNTPs may be 0.
CONDITION_MINIMUM = {"na_mm": 1.0, "primer_nm": 1.0}


def conditions():
    s = settings()
    out = {}
    for key, _, default, _ in CONDITIONS:
        try:
            out[key] = max(float(s.value(f"pcr/{key}", default)), CONDITION_MINIMUM.get(key, 0.0))
        except (TypeError, ValueError):   # a hand-edited settings file
            out[key] = default
    return out


def melting_temp(seq):
    """Tm in °C under the current PCR conditions, or None (non-ACGT or too short)."""
    return tm(seq, **conditions())


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
    s = structure(seq)
    text = f"hairpin {s['hairpin']:.1f}, self-dimer {s['self']:.1f} (3′ end {s['self3']:.1f}) kcal/mol"
    return text + ("  ⚠ " + "; ".join(s["warnings"]) if s["warnings"] else "  ✓")


def homology_settings():
    """(3'-end bases that must match, mismatches allowed in them) for the mis-priming check.
    Default 12 nt exact: ~0.5 chance hits per primer in a 5 Mb genome (1 mismatch: ~15)."""
    s = settings()
    return int(s.value("primers/homology_seed", 12)), int(s.value("primers/homology_mismatches", 0))


def default_genome():
    """The example E. coli BW25113 genome (demo folder, or next to the app), if present. Import > DNA > Background Genome from
    NCBI sets another one."""
    for path in (ROOT / "demo" / "BW25113_2026.gb", ROOT / "BW25113_2026.gb"):
        if path.exists():
            return str(path)
    return ""


def genome_path():
    """Off-target genome file, or "" when the check is off / no file."""
    s = settings()
    if s.value("primers/genome_check", "true") != "true":
        return ""
    return s.value("primers/genome", "") or default_genome()


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
    # The annealing part: a lowercase 5' tail does not need to match anywhere.
    core = sequence[len(sequence) - len(sequence.lstrip("acgtnryswkmbdhv")):] or sequence
    return g, core, mm, seed, (editor.doc.name, editor.length)


def genome_search(editor, sequence, progress=None, cancel=None):
    """Run (or fetch from the cache) the search of `sequence` in every sequence of the reference. None = no genome."""
    q = genome_query(editor, sequence)
    if q is None:
        return None
    g, core, mm, seed, skip = q
    return g.search(sequence, core, mm, seed, skip, progress, cancel)


def genome_cached(editor, sequence):
    q = genome_query(editor, sequence)
    if q is None:
        return None
    g, core, mm, seed, skip = q
    return g._cache.get((sequence.upper(), core.upper(), mm, seed, skip))


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
    rule = f"last {seed} nt, {'exact' if mm == 0 else f'≤{mm} mismatch'}"
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
    t = melting_temp(sequence)
    if t is not None:
        checks.append(("Tm", band(t, limit("tm_low"), limit("tm_high"), limit("tm_margin")),
                       f"{t:.1f} °C ({len(sequence)} nt)"))
    up = sequence.upper()
    gc = 100 * sum(c in "GCS" for c in up) / len(up)
    checks.append(("GC", band(gc, limit("gc_low"), limit("gc_high"), limit("gc_margin")), f"{gc:.0f} %"))
    clamp, lo, hi = sum(c in "GC" for c in up[-5:]), limit("clamp_min"), limit("clamp_max")
    checks.append(("3′ clamp", "ok" if lo <= clamp <= hi else "warn",
                   f"{clamp} G/C in the last 5 bases" + ("" if lo <= clamp <= hi else f" ({lo}–{hi} is best)")))
    s = structure(sequence)
    hp = min(s["hairpin"], 0.0)
    checks.append(("Hairpin", "bad" if s["hairpin"] < limit("hairpin_bad") else
                   "warn" if any(w.startswith("hairpin") for w in s["warnings"]) else "ok", f"ΔG {hp:.1f} kcal/mol"))
    checks.append(("Self-dimer", "bad" if s["self3"] < limit("dimer3_bad") or s["self"] < limit("dimer_bad") else
                   "warn" if any("dimer" in w for w in s["warnings"]) else "ok",
                   f"ΔG {s['self']:.1f} (3′ end {s['self3']:.1f}) kcal/mol"))

    n = editor.length
    sites = editor.doc.binding_sites(sequence, binding_mismatches())
    checks.append(("Binds here", "ok" if sites else "warn",
                   _listing([site_text(*x, n) for x in sites]) if sites else "nowhere in this sequence"))
    seed, mm = homology_settings()
    rule = f"last {seed} nt, {'exact' if mm == 0 else f'≤{mm} mismatch'}"
    # The annealing part: a lowercase 5' tail does not need to match anywhere.
    core = sequence[len(sequence) - len(sequence.lstrip("acgtnryswkmbdhv")):] or sequence
    full = editor.doc.binding_sites(core, 0, len(core))  # the whole annealing part, exact
    bound = {three_prime_end(x) for x in sites}
    others = [h for h in editor.doc.binding_sites(sequence, mm, seed) if three_prime_end(h) not in bound]
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
        self.setWindowTitle(f"Primer {name}" if self.old else "New Primer")
        self.setMinimumWidth(640)
        self.name = QLineEdit(self.old["name"] if self.old else self._free_name())
        self.name.setMaxLength(max(MAX_NAME, len(self.name.text())))   # not shorter than a name already in the database
        start = self.old["sequence"] if self.old else editor.selected_sequence()
        if len(start) > MAX_PRIMER:   # a long selection is not a primer: leave the box empty rather than cut it
            start = ""
        self.sequence = dna_field(QLineEdit(start))
        self.sequence.setMaxLength(max(MAX_PRIMER, len(start)))
        self.sequence.setPlaceholderText(f"5'→3'; a lowercase 5' tail is fine (up to {MAX_PRIMER} letters)")
        self.note = QLineEdit(self.old["note"] if self.old else "")
        self.note.setMaxLength(max(MAX_NOTE, len(self.note.text())))
        self.reverse = QCheckBox("Reverse complement")
        self.reverse.toggled.connect(lambda: self.sequence.setText(reverse_complement(self.sequence.text())))
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
        t = melting_temp(seq)
        self.info.setText(f"{len(seq)} nt, Tm {t:.1f} °C (whole primer)" if t is not None else f"{len(seq)} nt")
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
        ok = update_primer(self.editor, self.old["name"], name, seq, note) if self.old else add_primer(self.editor, name, seq, note) in ("added", "exists")
        if not ok:
            QMessageBox.warning(self, "Name in use", f"The primer database already has a primer '{name}'.")
            return
        super().accept()


class OptionsDialog(QDialog):
    """PCR conditions for Tm, which near-matches to show, and the 3' mis-priming check."""

    def __init__(self, parent):
        super().__init__(parent)
        self.setWindowTitle("Primer Options")
        self.setMinimumWidth(540)
        form = QFormLayout(self)
        current = conditions()
        self.boxes = {}
        for key, label, _, maximum in CONDITIONS:
            box = QDoubleSpinBox(decimals=2, minimum=CONDITION_MINIMUM.get(key, 0.0), maximum=maximum, value=current[key])
            form.addRow(label, box)
            self.boxes[key] = box
        self.binding_mm = QSpinBox(minimum=0, maximum=3, value=binding_mismatches())
        self.binding_mm.setToolTip("Show sites where the primer's last 15 bases match with up to this many\n"
                                   "mismatches (never in the 3 bases at the 3' end). 0 = exact sites only.")
        form.addRow("Show sites with up to … mismatches", self.binding_mm)
        seed, mm = homology_settings()
        self.seed = QSpinBox(minimum=6, maximum=20, value=seed)
        self.seed.setToolTip("Bases at the primer's 3' end compared for the mis-priming check")
        self.homology_mm = QSpinBox(minimum=0, maximum=3, value=mm)
        form.addRow("3′ check: bases", self.seed)
        form.addRow("3′ check: mismatches", self.homology_mm)
        s = settings()
        self.genome_check = QCheckBox("Also check a genome", checked=s.value("primers/genome_check", "true") == "true")
        self.genome = QLineEdit(s.value("primers/genome", ""))
        self.genome.setPlaceholderText(f"default: {default_genome() or 'none found'}")
        row = QHBoxLayout()
        row.addWidget(self.genome)
        row.addWidget(QPushButton("Browse…", clicked=self._browse))
        form.addRow("", self.genome_check)
        form.addRow("Genome file", row)
        self.genome_note = QLabel(wordWrap=True)
        ref = genome()
        self.genome_note.setText(f"Background genome now: {genome_label_text(ref)}. A FASTA or GenBank file with several sequences "
                                 "(chromosomes) is searched in full, one sequence at a time." if ref is not None
                                 else "No readable background genome now. A FASTA or GenBank file with several sequences "
                                 "(chromosomes) is searched in full.")
        form.addRow("", self.genome_note)
        form.addRow("✓ / ? / ✗ limits", QPushButton("Check thresholds (Tm, GC, hairpin…)…",
                                                   clicked=lambda: SettingsDialog(self, CHECKS).exec()))
        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel | QDialogButtonBox.RestoreDefaults)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        buttons.button(QDialogButtonBox.RestoreDefaults).clicked.connect(self._defaults)
        form.addRow(buttons)

    def _browse(self):
        path, _ = QFileDialog.getOpenFileName(self, "Genome for off-target check", "",
                                              "Sequence files (*.gb *.gbk *.gbff *.genbank *.dna *.fa *.fasta *.fna *.fas);;All files (*)")
        if path:
            self.genome.setText(path)

    def _defaults(self):
        for key, _, default, _ in CONDITIONS:
            self.boxes[key].setValue(default)
        self.binding_mm.setValue(1)
        self.seed.setValue(12)
        self.homology_mm.setValue(0)
        self.genome_check.setChecked(True)
        self.genome.clear()

    def accept(self):
        path = self.genome.text().strip()
        if self.genome_check.isChecked() and path and not Path(path).is_file():
            QMessageBox.warning(self, "Genome file", f"There is no file '{path}'. Browse for the genome (FASTA or GenBank), "
                                "or clear the box to use the default.")
            self.genome.setFocus()
            return
        s = settings()
        for key, box in self.boxes.items():
            s.setValue(f"pcr/{key}", box.value())
        s.setValue("primers/binding_mismatches", self.binding_mm.value())
        s.setValue("primers/homology_seed", self.seed.value())
        s.setValue("primers/homology_mismatches", self.homology_mm.value())
        s.setValue("primers/genome_check", "true" if self.genome_check.isChecked() else "false")
        s.setValue("primers/genome", self.genome.text().strip())
        if self.genome_check.isChecked() and self.genome.text().strip():
            import reference_picker
            reference_picker.remember(self.genome.text().strip())   # the toolbar's box lists it
        notifier.changed.emit()  # re-match primers with the new rules
        super().accept()


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


MISMATCHES = Qt.UserRole + 2   # on a site's sequence cell: the indices of its bases that do not pair (drawn red)
GROUP_ROLE = Qt.UserRole + 3   # on the rows that group a primer's sites ("On this sequence", "On the reference genome")
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
    """The primers' tree. Its Name column is as wide as the names and sites: the rows that group a primer's sites ("On the
    background genome, …") run across the whole row, so their long label does not count."""

    def sizeHintForColumn(self, column):
        if column != 0:
            return super().sizeHintForColumn(column)
        fm, widest = self.fontMetrics(), 0

        def walk(item, depth):
            nonlocal widest
            if not item.isFirstColumnSpanned():
                widest = max(widest, fm.horizontalAdvance(item.text(0)) + depth * self.indentation())
            if item.isExpanded():
                for k in range(item.childCount()):
                    walk(item.child(k), depth + 1)
        for i in range(self.topLevelItemCount()):
            walk(self.topLevelItem(i), 1 if self.rootIsDecorated() else 0)
        return widest + 12


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
    DATABASE_COLUMNS = ["Added", "Note"]   # the note last: it takes the room left, the date stays in view
    FILTERED = ["Match", "Structure", "Added"]   # filter by values (value_filter.py)

    def __init__(self, editor):
        super().__init__()
        self.editor = editor
        self.database = getattr(editor, "database", None)
        self.columns = self.COLUMNS + (self.DATABASE_COLUMNS if self.database else [])
        self.tree = PrimerTree()
        self.tree.setHeaderLabels(self.columns)
        self.tree.setUniformRowHeights(True)
        # (before fit_columns: a tree's header stretches its last column by default, and that width would count as one the user
        # dragged, never fitted again)
        self.tree.header().setStretchLastSection(False)
        # Columns as wide as what they hold (a long primer's whole sequence, up to 480 px); the note (the database), else the
        # structure, takes the room left.
        fit_columns(self.tree, stretch=self.columns.index("Note" if self.database else "Structure"),
                    cap=self.tree.fontMetrics().horizontalAdvance("G" * 64) + 24)   # a long primer (60 nt and more) read whole
        self.tree.setToolTip("Tm: of the bases that anneal at that site (5' tail excluded)")
        self.tree.setSelectionMode(QAbstractItemView.ExtendedSelection)   # several primers: Import to This Sequence for all
        self.tree.currentItemChanged.connect(self._select)
        self.tree.itemDoubleClicked.connect(lambda item, _: self.edit(item))
        self.tree.setContextMenuPolicy(Qt.CustomContextMenu)
        self.tree.customContextMenuRequested.connect(self._menu)
        self.expanded = set()  # primer names the user opened; kept across refreshes
        self.tree.itemExpanded.connect(self._expanded)
        self.tree.itemCollapsed.connect(lambda item: item.parent() is None and self.expanded.discard(item.data(0, Qt.UserRole)))
        self.tree.setItemDelegateForColumn(1, MismatchDelegate(self.tree))
        self._searching = set()   # primers whose genome search is running
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
        self.reload()

    def _can_edit(self):
        if self.edit_button is not None:
            self.edit_button.setEnabled(self.tree.currentItem() is not None and bool(self.tree.selectedItems()))

    def theme_changed(self):
        """Light / dark switch (theme.refresh_views): the greyed (hidden) and temporary primers take the new colours."""
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
        rows = {}
        if self.database:   # every primer of the database, also one that binds nowhere (an ambiguity code...)
            rows = {r["name"]: r for r in db().rows}
            groups = {name: groups.get(name, []) for name in rows}
        self.tree.blockSignals(True)  # don't select sites while refilling
        current = self.tree.currentItem()
        keep = primer_item(current).data(0, Qt.UserRole) if current is not None else None
        self.tree.clear()
        for name, sites in groups.items():
            away = elsewhere.get(name, [])
            seq = (sites or away)[0].sequence if sites or away else rows[name]["sequence"]
            exact = sum(s.mismatches == 0 for s in sites + away)
            t = melting_temp(seq)
            total = len(sites) + len(away)
            summary = f"{total} site{'s' * (total > 1)}" + (f", {exact} exact" if exact != total else "") + \
                (f" ({len(away)} on other sequences)" if away else "")
            warnings = structure(seq)["warnings"]
            mark, color = MARKS["warn" if warnings else "ok"]
            group = SortItem([name, seq, str(len(seq)), f"{t:.1f}" if t is not None else "–", summary if total else "no site",
                              f"{mark} " + "; ".join(warnings) if warnings else mark]
                             + ([rows[name]["added"], rows[name]["note"]] if self.database else []))
            group.keys = {0: name.lower(), 2: len(seq), 3: t if t is not None else -1e9, 4: total}
            group.setForeground(5, QColor(color))
            group.setData(0, Qt.UserRole, name)
            if sites and all(s.provisional for s in sites) and not self.database:   # drawn in the temporary colour: not in this file
                group.setForeground(0, QColor(TEMP_COLOR))
                group.setToolTip(0, "Not part of this file: from the primer database (or designed, not saved). "
                                    "Right-click > Import to This Sequence keeps it in the file.")
            if name in ed.hidden_primers:   # still listed, greyed: right-click > Show in This File
                group.setText(0, f"{name} (hidden)")
                for col in range(len(self.columns)):
                    group.setForeground(col, self.palette().color(QPalette.PlaceholderText))
                group.setToolTip(0, "Hidden in this file (right-click > Show in This File); still in the primer database")
            self._add_sites(group, seq, sites, away, lengths)
            for item in [group] + _below(group):
                for col in (2, 3):
                    item.setTextAlignment(col, Qt.AlignRight | Qt.AlignVCenter)
            self.tree.addTopLevelItem(group)
            group.setExpanded(name in self.expanded)  # collapsed unless the user opened it
            for i in range(group.childCount()):   # the groups of sites (this sequence, the genome) open with it
                group.child(i).setExpanded(True)
                group.child(i).setFirstColumnSpanned(True)   # their label across the row (only once it is in the tree)
            if name == keep:   # the primer picked stays picked (its row only: the sites may have moved)
                self.tree.setCurrentItem(group)
        self._apply_sort()
        # No row opens: no column of empty room for the ▸ on every row.
        self.tree.setRootIsDecorated(any(self.tree.topLevelItem(i).childCount() for i in range(self.tree.topLevelItemCount())))
        self.tree.blockSignals(False)
        n = len(groups)
        self.count.setText(plural(n, "primer") if self.database else f"{plural(n, 'primer')} on this sequence")
        self._can_edit()

    # ---- the sites of a primer, grouped by where they are -------------------------------------------------------------------
    def _add_sites(self, group, seq, sites, away, lengths):
        """Under a primer: the sites on this sequence, on the other sequences of the file, and in the background genome, each a
        group of its own; in each site's sequence the bases that do not pair with the template are red (MismatchDelegate)."""
        ed = self.editor
        if sites and not self.database:   # (the database's own list document is not a sequence of yours)
            rows = []
            for site in sites:
                bound = seq[-sum(e - s for s, e in site.spans):]
                template = "".join(ed.doc.sequence(a, b) for a, b in site.spans)
                rows.append(self._site_row(location_text(site.start, site.end, site.strand, ed.length), bound,
                                           template if site.strand > 0 else reverse_complement(template),
                                           ("", site.start), site.key))
            self._group(group, f"On this sequence ({len(rows)})", rows)
        if away:
            rows = [self._site_row(f"{site.record}: " + location_text(site.start, site.end, site.strand, lengths.get(site.record, site.end)),
                                   site.sequence[-(site.end - site.start):],
                                   getattr(site, "bases", None) if site.strand > 0 or getattr(site, "bases", None) is None
                                   else reverse_complement(site.bases),
                                   (site.record, site.start), (site.record, site.start, site.end))   # a click opens that record
                    for site in away]
            self._group(group, f"On the other sequences of this file ({len(rows)})", rows)
        if group.data(0, Qt.UserRole) in self.expanded:   # the genome is searched only for a primer opened
            self._genome_group(group, seq)

    def _genome_group(self, group, seq):
        q = genome_query(self.editor, seq)
        if q is None:
            return
        result = genome_cached(self.editor, seq)
        if result is None and seq.upper() in self._failed:
            self._group(group, f"On the background genome: the search failed ({self._failed[seq.upper()]})", [])
            return
        if result is None:
            self._search_genome(seq)
            self._group(group, "On the background genome: searching…", [])
            return
        rows = []
        for whole, hits in ((True, result.whole), (False, result.hits)):
            for h in hits:
                name, start, end, strand, _, length, circular = h
                bases = result.bases.get((name, start, end, strand)) if hasattr(result, "bases") else None
                rows.append(self._site_row(site_in_record(h), seq[-(end - start):],
                                           None if bases is None else bases if strand > 0 else reverse_complement(bases),
                                           (name, start), None, "whole primer" if whole else "3′ end"))
        more = result.n_whole + result.n_hits - len(rows)
        label = f"On the background genome, {q[0].name} ({result.n_whole + result.n_hits}" + (f", {more} not listed" if more > 0 else "") + ")"
        if not result.complete:
            label += ": not searched to the end"
        if not rows:
            label = f"On the background genome, {q[0].name}: no site" + ("" if result.complete else " in the part searched")
        self._group(group, label, rows)

    def _site_row(self, where, bound, template, sort_key, data, kind=""):
        """One binding site: where, the primer's annealing bases (mismatches red when the template is known), length, Tm, match."""
        tb = melting_temp(bound)
        marks = mismatch_indices(bound, template) if template is not None and len(template) == len(bound) else None
        match = "" if marks is None else "exact" if not marks else plural(len(marks), 'mismatch')
        row = SortItem([where, bound, str(len(bound)), f"{tb:.1f}" if tb is not None else "–",
                        ", ".join(x for x in (kind, match) if x)])
        row.keys = {0: sort_key, 2: len(bound), 3: tb if tb is not None else -1e9, 4: len(marks or [])}
        row.setData(0, Qt.UserRole, data)
        row.setData(1, MISMATCHES, marks or [])
        if marks:
            row.setToolTip(1, "Red: the primer's base does not pair with the template there")
        return row

    def _group(self, parent, label, rows):
        g = QTreeWidgetItem([label])
        f = g.font(0)
        f.setItalic(True)
        g.setFont(0, f)
        g.setData(0, GROUP_ROLE, True)
        parent.addChild(g)
        g.addChildren(rows)
        return g

    def _search_genome(self, seq):
        """Search the background genome for `seq` on a worker thread, then show the sites (one search at a time per sequence)."""
        if seq.upper() in self._searching:
            return
        self._searching.add(seq.upper())
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
        self._reload_later.start()

    def _expanded(self, item):
        if item.parent() is None:
            name = item.data(0, Qt.UserRole)
            if name not in self.expanded:
                self.expanded.add(name)
                # Its sites in the reference genome: searched now. The list is redone after this event, not inside it (Qt may be
                # opening the primer to show a site being selected: deleting that row now would crash).
                self._reload_later.start()

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
        if item.parent() is None:  # a primer: select its first site on this sequence
            key = next((row.data(0, Qt.UserRole) for row in _below(item) if isinstance(row.data(0, Qt.UserRole), str)), None)
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
            menu.exec(self.tree.viewport().mapToGlobal(pos))
        elif item is not None:   # the same menu as on a primer in the sequence view
            from feature_table import primer_actions, primers_actions  # here: feature_table imports the tools
            menu = QMenu(self)
            name = primer_item(item).data(0, Qt.UserRole)
            chosen = {primer_item(i).data(0, Qt.UserRole) for i in self.tree.selectedItems()}
            if len(chosen) > 1 and name in chosen:   # several primers selected: every action applies to all of them
                primers_actions(menu, self.editor, chosen)
            else:
                primer_actions(menu, self.editor, name)
            menu.exec(self.tree.viewport().mapToGlobal(pos))

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
        db_tabs.open_tab(self.window(), "primers")

    def edit_options(self):
        if OptionsDialog(self).exec() == QDialog.Accepted:
            self._failed.clear()   # another genome, other rules: searched again
            self.reload()
            self.editor.selection_changed.emit()  # status bar Tm
