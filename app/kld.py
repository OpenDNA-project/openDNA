"""KLD: the primers that make the unsaved edits of a plasmid (point mutations, insertions, deletions) by inverse PCR + KLD.

Edit the sequence in place as usual; the KLD tab compares it with the sequence as last saved, takes the first change (and the next
ones while they all fit within `group` bp), and gives a back-to-back primer pair: the forward primer starts at the change and runs
right (new bases first, then bases that anneal), the reverse primer starts just left of it and runs left. Bases to insert or
substitute ride on the 5' ends (split between both primers when long). `design` is pure logic; the panel is at the end.
"""
from difflib import SequenceMatcher
from html import escape
from types import SimpleNamespace

from PySide6.QtGui import QGuiApplication
from PySide6.QtWidgets import QDoubleSpinBox, QHBoxLayout, QLabel, QPushButton, QSpinBox, QTextBrowser, QVBoxLayout, QWidget

from assembly import TAIL_COLOR, anneal_len, oligo_checks
from editor import reverse_complement as rc
from primer_panel import melting_temp

MAX_INSERT = 500    # bp: new bases beyond this are not made with primers
ONE_TAIL = 40       # bp: up to here the new bases ride on the forward primer alone; beyond, they are split between the two
MIN_LENGTH = 100    # bp: a shorter plasmid cannot take two annealing parts and a change
MAX_MIDDLE = 20_000  # bp: a bigger stretch between the first and last difference is taken as one change, not aligned
LONG_PRIMER = 70    # nt: a longer primer is flagged


def changes(old, new):
    """[(start, end, text)]: how `old` became `new`, as replacements of old[start:end] by `text` (in order). Equal bases at both
    ends are skipped first, so a change across the origin of a circle shows as two."""
    p = 0
    while p < min(len(old), len(new)) and old[p] == new[p]:
        p += 1
    s = 0
    while s < min(len(old), len(new)) - p and old[len(old) - 1 - s] == new[len(new) - 1 - s]:
        s += 1
    a, b = old[p:len(old) - s], new[p:len(new) - s]
    if not a and not b:
        return []
    if len(a) > MAX_MIDDLE or len(b) > MAX_MIDDLE:   # shortcut: one change, not an alignment; split it if that ever matters
        return [(p, p + len(a), b)]
    ops = SequenceMatcher(None, a, b, autojunk=False).get_opcodes()
    return [(p + i1, p + i2, b[j1:j2]) for tag, i1, i2, j1, j2 in ops if tag != "equal"]


def design(old, new, circular=True, tm=60.0, group=30, warn_insert=120):
    """The KLD primers that turn `old` (the sequence as saved) into `new` (as edited). Returns SimpleNamespace(rows, notes,
    product, change, skipped): rows [(name, sequence, annealing length)], notes [(level, text)] ("warn" is the orange one), the
    sequence the primers make, the change designed (start, end, text) on `old`, and the changes left out. Raises ValueError
    when it cannot be done."""
    old, new = old.upper(), new.upper()
    if not circular:
        raise ValueError("KLD needs a circular sequence: the product is closed into a plasmid.")
    if len(old) < MIN_LENGTH:
        raise ValueError(f"The sequence is too short for KLD (under {MIN_LENGTH} bp).")
    found = changes(old, new)
    if not found:
        raise ValueError("No unsaved change: edit the sequence first (mutate, insert or delete bases).")
    taken = [c for c in found if c[1] - found[0][0] <= group] or found[:1]   # the first change, then those near enough
    skipped = found[len(taken):]
    a, b = found[0][0], taken[-1][1]
    text = new[a:len(new) - (len(old) - b)]   # what replaces old[a:b]: bases after b are the same count in both
    if len(text) > MAX_INSERT:
        raise ValueError(f"{len(text):,} new bp are too many for primers (the limit is {MAX_INSERT} bp).")
    notes = []
    if len(text) > warn_insert:
        notes.append(("warn", f"{len(text):,} new bp (over {warn_insert}): the primers are long and may be hard or expensive to order."))
    if not text:
        notes.append(("info", f"Deletes {b - a:,} bp."))
    n = len(old)
    wide = old * 3
    right = wide[b:b + 60]
    left = rc(wide[a + n - 60:a + n])   # the bases just left of the change, read away from it
    nf, nr = anneal_len(right, tm), anneal_len(left, tm)
    k = 0 if len(text) <= ONE_TAIL else len(text) // 2   # bases of the text carried by the reverse primer
    rows = [("KLD-F", text[k:].lower() + right[:nf], nf), ("KLD-R", rc(text[:k]).lower() + left[:nr], nr)]
    # The product, read from the primers alone, must be the edited plasmid (up to where the circle is opened).
    gap = n - (b - a) - nf - nr
    if gap < 0:
        raise ValueError("The sequence is too short for these primers.")
    product = rows[0][1].upper() + wide[b + nf:b + nf + gap] + rc(rows[1][1].upper())
    expected = old[:a] + text + old[b:]
    if len(product) != len(expected) or product not in expected + expected:
        raise ValueError("The primers do not rebuild the edited sequence (unexpected: please report it).")
    for name, seq, an in rows:
        for label, level, what in oligo_checks((name, seq, an, "")):
            if level in ("warn", "bad"):
                notes.append(("warn", f"{name}: {label} {what}"))
        if len(seq) > LONG_PRIMER:
            notes.append(("warn", f"{name} is {len(seq)} nt long."))
    if skipped:
        notes.append(("info", f"{len(skipped)} other change{'s' if len(skipped) > 1 else ''} not included (more than {group} bp "
                              "from the first): save, then run KLD again."))
    return SimpleNamespace(rows=rows, notes=notes, product=expected, change=(a, b, text), skipped=skipped)


# ---- the panel ---------------------------------------------------------------------------------------------


class KldPanel(QWidget):
    """Tools > KLD (also an Assembly choice): the primers for the unsaved edits of the open sequence."""

    def __init__(self, win):
        super().__init__()
        self.win = win
        self.tm = QDoubleSpinBox(minimum=45, maximum=72, value=60, suffix=" °C")
        self.group = QSpinBox(minimum=1, maximum=200, value=30, suffix=" bp")
        self.warn = QSpinBox(minimum=1, maximum=MAX_INSERT, value=120, suffix=" bp")
        self.group.setToolTip("Changes closer than this are made by the same primer pair.")
        self.warn.setToolTip(f"More new bases than this get an orange warning (primers may be hard or expensive to order); "
                             f"over {MAX_INSERT} bp the change is refused.")
        top = QHBoxLayout()
        for label, w in (("Annealing Tm", self.tm), ("Group changes within", self.group), ("Warn for insertions above", self.warn)):
            top.addWidget(QLabel(label))
            top.addWidget(w)
            w.valueChanged.connect(self.refresh)
        top.addStretch()
        top.addWidget(QPushButton("Refresh", clicked=self.refresh, toolTip="Read the unsaved edits of the open sequence again"))
        top.addWidget(QPushButton("Copy primers", clicked=self.copy_primers))
        self.page = QTextBrowser()
        lay = QVBoxLayout(self)
        lay.addLayout(top)
        lay.addWidget(self.page)
        self.result = None

    def showEvent(self, e):
        super().showEvent(e)
        self.refresh()

    def copy_primers(self):
        if self.result:
            QGuiApplication.clipboard().setText("".join(f"{n}\t{s}\n" for n, s, _ in self.result.rows))

    def refresh(self, *_):
        ed = getattr(self.win, "editor", None)
        self.result = None
        if ed is None:
            self.page.setHtml("<p>Open a sequence, edit it (mutate, insert or delete bases), then come back here.</p>")
            return
        old = ed.doc.saved_sequence()
        try:
            if old is None:
                raise ValueError("There is no saved sequence to compare with (a new sequence, or one rotated or "
                                 "reverse-complemented since the last save): save it, then edit.")
            self.result = r = design(old, ed.doc.sequence(), ed.circular, self.tm.value(), self.group.value(), self.warn.value())
        except ValueError as e:
            self.page.setHtml(f"<p style='color:#c62828'>{escape(str(e))}</p>")
            return
        a, b, text = r.change
        what = "Mutation" if len(text) == b - a else "Insertion" if a == b else "Deletion" if not text else "Replacement"
        html = [f"<p><b>{what}</b> at {a + 1:,}–{b:,} ({b - a:,} bp → {len(text):,} bp)</p>", "<table cellspacing=6>"]
        for name, seq, n in r.rows:
            tail = len(seq) - n
            t = melting_temp(seq[tail:])
            html.append(f"<tr><td><b>{name}</b></td><td style='font-family:monospace'><span style='color:{TAIL_COLOR}'>"
                        f"{seq[:tail]}</span>{seq[tail:]}</td><td>{len(seq)} nt</td><td>{f'{t:.1f} °C' if t else ''}</td></tr>")
        html.append("</table>")
        for level, note in r.notes:
            html.append(f"<p style='color:{'#ef6c00' if level == 'warn' else 'gray'}'>{'⚠ ' if level == 'warn' else ''}{escape(note)}</p>")
        self.page.setHtml("".join(html))
