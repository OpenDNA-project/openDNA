"""Primer pair design around a selected target (our own scorer; no Primer3, GPL).

Forward primers end (3') at most `window` bases before the target, reverse primers
start (3') at most `window` bases after it, so the product always contains the whole
target. Each candidate gets a penalty (lower is better):
- distance from the target Tm (1 per °C),
- GC outside 40-60 % (1 per 5 %),
- no G/C at the 3' end (+1), more than 3 G/C in the last 5 bases (+1),
- runs of 5+ identical bases (+2),
- hairpin / 3' self-dimer ΔG beyond the usual limits (1 per kcal/mol),
- 3' end (12 nt) matching elsewhere in the sequence (+5 per site).
A pair adds the Tm difference and their 3' cross-dimer ΔG beyond -5 kcal/mol.
"""
import re

from PySide6.QtCore import Qt
from PySide6.QtGui import QColor
from PySide6.QtWidgets import (
    QAbstractItemView, QDialogButtonBox, QDoubleSpinBox, QFormLayout, QHBoxLayout,
    QInputDialog, QLabel, QMessageBox, QPushButton, QSpinBox, QTableWidget, QTableWidgetItem,
    QVBoxLayout,
)

from opendna_core import dimer, hairpin, tm

from editor import reverse_complement
from primer_db import add_primers
from primer_panel import (MARKS, MAX_NAME, ReportView, checks_html, conditions, homology_settings, primer_checks, primer_problem,
                          warm_genome, worst)
from tools import ToolDialog
from uikit import fit_columns, no_auto_default, plural, right_aligned, run_busy

KEEP = 25  # best candidates per side that get the (slower) uniqueness check and pairing


def salt(cond):
    return {k: cond[k] for k in ("na_mm", "mg_mm", "dntp_mm")}


def candidate(seq, tm_target, cond):
    """Penalty and properties of one primer sequence (uppercase), or None if hopeless.
    `cond`: PCR conditions, read once per design run (settings reads are slow)."""
    t = tm(seq, **cond)
    if t is None or abs(t - tm_target) > 5:
        return None
    gc = 100 * sum(c in "GC" for c in seq) / len(seq)
    penalty = abs(t - tm_target) + max(0.0, abs(gc - 50) - 10) / 5
    notes = []
    if seq[-1] not in "GC":
        penalty += 1
        notes.append("no 3′ G/C")
    if sum(c in "GC" for c in seq[-5:]) > 3:
        penalty += 1
        notes.append("GC-rich 3′")
    if re.search(r"(A{5,}|C{5,}|G{5,}|T{5,})", seq):
        penalty += 2
        notes.append("run")
    hp, hp3 = hairpin(seq, **salt(cond))
    sd, sd3 = dimer(seq, seq, **salt(cond))
    penalty += max(0.0, -2 - hp3) + max(0.0, -3 - hp) + max(0.0, -5 - sd3) + max(0.0, -9 - sd)
    if hp3 < -2 or hp < -3:
        notes.append(f"hairpin {min(hp, hp3):.1f}")
    if sd3 < -5:
        notes.append(f"3′ dimer {sd3:.1f}")
    return {"seq": seq, "tm": t, "gc": gc, "penalty": penalty, "notes": notes}


def design(editor, start, end, tm_target=60.0, min_len=18, max_len=25, window=150, count=5):
    """Best primer pairs whose product contains [start, end). Returns a list of dicts:
    fwd/rev candidates, their positions, product size, pair penalty."""
    n = editor.length
    seq = editor.doc.sequence().upper()
    cond = conditions()
    fwds, revs = [], []
    # Forward: 3' end at p (bases [p - L, p)), p in [start - window, start].
    for p in range(max(min_len, start - window), start + 1):
        for length in range(min_len, max_len + 1):
            if p - length < 0:
                break
            c = candidate(seq[p - length:p], tm_target, cond)
            if c:
                c.update(start=p - length, end=p)
                fwds.append(c)
    # Reverse: anneals [q, q + L) with its 3' end at q, q in [end, end + window].
    for q in range(end, min(n - min_len, end + window) + 1):
        for length in range(min_len, max_len + 1):
            if q + length > n:
                break
            c = candidate(reverse_complement(seq[q:q + length]), tm_target, cond)
            if c:
                c.update(start=q, end=q + length)
                revs.append(c)
    seed, _ = homology_settings()
    for side in (fwds, revs):
        side.sort(key=lambda c: c["penalty"])
        del side[KEEP:]
        for c in side:  # other places the 3' end could prime, in this sequence
            extra = len(editor.doc.binding_sites(c["seq"], 0, seed)) - 1
            if extra > 0:
                c["penalty"] += 5 * extra
                c["notes"].append(f"3′ matches {plural(extra, 'other site')}")
    pairs = []
    for f in fwds:
        for r in revs:
            _, cross3 = dimer(f["seq"], r["seq"], **salt(cond))
            penalty = f["penalty"] + r["penalty"] + abs(f["tm"] - r["tm"]) + max(0.0, -5 - cross3)
            notes = [f"F: {x}" for x in f["notes"]] + [f"R: {x}" for x in r["notes"]]
            if cross3 < -5:
                notes.append(f"cross-dimer {cross3:.1f}")
            pairs.append({"fwd": f, "rev": r, "size": r["end"] - f["start"], "penalty": penalty, "notes": notes})
    pairs.sort(key=lambda p: p["penalty"])
    return pairs[:count]


class DesignDialog(ToolDialog):
    COLUMNS = ["Forward 5'→3'", "Reverse 5'→3'", "Tm F / R °C", "Product bp", "Penalty", "Notes"]

    def __init__(self, parent, editor):
        super().__init__(parent)
        self.editor, self.pairs = editor, []
        self.start, self.end = editor.start, editor.end
        self.setWindowTitle("Design Primers")
        self.resize(1000, 420)
        self.tm = QDoubleSpinBox(minimum=45, maximum=75, value=60, decimals=1, suffix=" °C")
        self.min_len = QSpinBox(minimum=12, maximum=35, value=18)
        self.max_len = QSpinBox(minimum=12, maximum=40, value=25)
        self.min_len.valueChanged.connect(lambda v: self.max_len.setValue(max(v, self.max_len.value())))   # shortest <= longest
        self.max_len.valueChanged.connect(lambda v: self.min_len.setValue(min(v, self.min_len.value())))
        self.window = QSpinBox(minimum=20, maximum=2000, value=150, suffix=" bp")
        self.window.setToolTip("How far outside the target the primers may sit, on each side")
        form = QFormLayout()
        self.target = QLabel()
        self._set_target()
        form.addRow("Target", self.target)
        form.addRow("Tm", self.tm)
        lengths = QHBoxLayout()
        lengths.addWidget(self.min_len)
        lengths.addWidget(QLabel("to"))
        lengths.addWidget(self.max_len)
        form.addRow("Length", lengths)
        form.addRow("Search window", self.window)
        self.table = QTableWidget(0, len(self.COLUMNS))
        self.table.setHorizontalHeaderLabels(self.COLUMNS)
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SingleSelection)
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        fit_columns(self.table, stretch=5)
        self.table.currentCellChanged.connect(lambda row, *_: self._show(row))
        self.destroyed.connect(lambda *_, ed=editor: self._forget(ed))   # the pane is closed: its primers go with it
        self.details = ReportView(10)   # the same height whatever the pair clicked: the table above does not move
        run = QPushButton("Design", clicked=self.run)
        self.save_button = save = QPushButton("Save pair to database…", clicked=self.save, enabled=False)   # once there is a pair
        buttons = QDialogButtonBox(QDialogButtonBox.Close)
        buttons.rejected.connect(self.reject)
        row = QHBoxLayout()
        row.addWidget(run)
        row.addWidget(save)
        row.addStretch()
        row.addWidget(buttons)
        layout = QVBoxLayout(self)
        layout.addLayout(form)
        layout.addWidget(self.table)
        layout.addWidget(self.details)
        layout.addLayout(row)
        no_auto_default(self)
        run.setDefault(True)   # Enter in a field designs (never closes the window or saves a pair)

    def _set_target(self):
        self.target.setText(f"{self.start + 1:,}–{self.end:,} ({self.end - self.start:,} bp, the selection)")

    def run(self):
        if self.embedded:  # in a pane: design for whatever is selected now
            self.start, self.end = self.editor.start, self.editor.end
            self._set_target()
        if self.end <= self.start:
            QMessageBox.warning(self, "No target", "Select the region to amplify first.")
            return
        if self.end > self.editor.length:
            # shortcut: design() reads one linear stretch; rotate the template to support this.
            QMessageBox.warning(self, "Across the origin", "Designing primers across the origin is not supported "
                                "yet: select a region that does not pass base 1.")
            return
        self._forget(self.editor)
        args = (self.editor, self.start, self.end, self.tm.value(), self.min_len.value(),
                max(self.min_len.value(), self.max_len.value()), self.window.value())
        # (thousands of candidates: in a thread behind a busy dialog, so the window is drawn and stays alive meanwhile)
        self.pairs = run_busy(self, "Designing primers…", lambda: design(*args))
        self.save_button.setEnabled(bool(self.pairs))
        self.table.setRowCount(len(self.pairs))
        # The genome search (every chromosome) first, with a progress dialog; the checks below then read its results.
        warm_genome(self, self.editor, [p[k]["seq"] for p in self.pairs for k in ("fwd", "rev")])
        for i, p in enumerate(self.pairs):
            values = [p["fwd"]["seq"], p["rev"]["seq"], f"{p['fwd']['tm']:.1f} / {p['rev']['tm']:.1f}",
                      f"{p['size']:,}", f"{p['penalty']:.1f}", "; ".join(p["notes"])]
            # The pair's mark: the worst of all checks of both primers (details below the table).
            level = worst(primer_checks(self.editor, p["fwd"]["seq"], "cached") + primer_checks(self.editor, p["rev"]["seq"], "cached"))
            mark, color = MARKS["warn" if level == "ok" and p["notes"] else level]
            values[5] = f"{mark} {values[5]}".strip()
            for c, v in enumerate(values):
                item = QTableWidgetItem(v)
                if c in (2, 3, 4):
                    right_aligned(item)
                if c == 5:
                    item.setForeground(QColor(color))
                self.table.setItem(i, c, item)
        if not self.pairs:
            self.details.setText("No pair found: widen the window, the length range or the Tm.")
        else:
            self.table.selectRow(0)

    @staticmethod
    def _forget(editor):
        """Stop drawing the designed primers (they are not in the database)."""
        try:
            editor.set_extra([])
        except RuntimeError:
            pass   # the editor is already being deleted with its window

    def temporary(self, row):
        """The pair on `row` as temporary primers: [(name, sequence)], named pair 1 F / pair 1 R..."""
        p = self.pairs[row]
        return [(f"pair{row + 1}_F", p["fwd"]["seq"]), (f"pair{row + 1}_R", p["rev"]["seq"])]

    def _show(self, row):
        """The checks of both primers side by side (✓ / ? / ✗), sequence and genome hits apart; the pair is drawn as temporary
        primers where it binds."""
        if 0 <= row < len(self.pairs):
            p = self.pairs[row]
            self.editor.set_extra([], self.temporary(row))
            self.editor.reveal.emit((p["fwd"]["start"] + p["rev"]["end"]) // 2)
            fwd, rev = (primer_checks(self.editor, p[k]["seq"], "cached") for k in ("fwd", "rev"))
            self.details.setText("<table width='100%'><tr>"
                                 f"<td valign='top'><b>Forward</b>{checks_html(fwd)}</td>"
                                 f"<td valign='top'><b>Reverse</b>{checks_html(rev)}</td></tr></table>")

    def save(self, base=None):
        row = self.table.currentRow()
        if not 0 <= row < len(self.pairs):
            return None
        if base is None:
            base, ok = QInputDialog.getText(self, "Save pair", "Name (saved as <name>_F and <name>_R):",
                                            text=f"{self.editor.doc.name}_{self.start + 1}")
            if not ok or not base.strip():
                return None
        p = self.pairs[row]
        names = [f"{base.strip()}_F", f"{base.strip()}_R"]
        problem = primer_problem(names[0], p["fwd"]["seq"]) or primer_problem(names[1], p["rev"]["seq"])
        if problem or len(names[0]) > MAX_NAME:   # a tab in the name, or one too long for the database's name box
            QMessageBox.warning(self, "Save pair", problem or f"The name is too long (at most {MAX_NAME - 2} characters).")
            return None
        if add_primers(self.editor, [(names[0], p["fwd"]["seq"]), (names[1], p["rev"]["seq"])], "Saved", in_file=True)["conflict"]:
            QMessageBox.warning(self, "Name in use", "A primer with that name and another sequence exists.")
            return None
        self._forget(self.editor)   # now regular primers, under their new names
        self.details.setText(f"Saved {names[0]} and {names[1]} to the primer database.")
        if self.editor is not None and not self.embedded:
            self.editor.error.emit(f"Saved {names[0]} and {names[1]} to the primer database.")   # status bar: the window closes
            self.accept()
        return names


def design_for_selection(window, editor):
    """Right-click a selection > Design Primers for This Sequence (or Primers > Design Primers for Selection): the designer in
    its own window over the main one (not modal: the views show where each pair binds). It closes once a pair is saved to the
    primer database; closing it otherwise drops the pairs it drew."""
    if editor.end <= editor.start:
        QMessageBox.information(window, "Design Primers", "Select the region to amplify first.")
        return None
    dialog = DesignDialog(window, editor)
    dialog.setAttribute(Qt.WA_DeleteOnClose)
    dialog.finished.connect(lambda _result, ed=editor: DesignDialog._forget(ed))
    dialog.show()
    dialog.run()
    return dialog
