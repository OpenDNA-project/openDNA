"""KLD: the primers that make the unsaved edits of a plasmid (point mutations, insertions, deletions) by inverse PCR + KLD.

Edit the sequence in place as usual; the KLD method of the Assembly panel compares it with the sequence as last saved, takes the first change (and the next
ones while they all fit within `group` bp), and gives a back-to-back primer pair: the forward primer starts at the change and runs
right (new bases first, then bases that anneal), the reverse primer starts just left of it and runs left. Bases to insert or
substitute ride on the 5' ends (split between both primers when long). `design` is pure logic; the panel is at the end.
"""
from difflib import SequenceMatcher
from pathlib import Path
from html import escape
from types import SimpleNamespace

from PySide6.QtCore import QTimer, Qt
from PySide6.QtGui import QColor, QGuiApplication
from PySide6.QtWidgets import (
    QComboBox, QDoubleSpinBox, QHBoxLayout, QLabel, QListWidget, QListWidgetItem, QPushButton, QSpinBox, QSplitter, QTextBrowser,
    QVBoxLayout, QWidget,
)

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


def auto(found, group=30, n=None):
    """The changes selected to begin with: the first, then the next ones while all of them fit within `group` bp. With `n`
    (the length of a circular sequence), the last changes join too when they are within `group` bp of the first across the origin."""
    if not found:
        return []
    chosen = [found[0]] + [c for c in found[1:] if c[1] - found[0][0] <= group]   # (the first one whatever its size)
    if n:
        end = chosen[-1][1]
        for c in reversed(found[len(chosen):]):
            if n - c[0] + end > group:
                break
            chosen.insert(0, c)
    return chosen


def design(old, new, circular=True, tm=60.0, group=30, warn_insert=120):
    """`design_selected` for the changes `auto` picks; `.skipped` lists the others."""
    found = changes(old.upper(), new.upper())
    chosen = auto(found, group, len(old) if circular else None)
    r = design_selected(old, chosen, circular, tm, group, warn_insert)
    r.skipped = [c for c in found if c not in chosen]
    if r.skipped:
        r.notes.append(("info", f"{len(r.skipped)} other change{'s' if len(r.skipped) > 1 else ''} not included (more than {group} bp "
                                "from the first): save, then run KLD again."))
    return r


def design_selected(old, chosen, circular=True, tm=60.0, group=30, warn_insert=120):
    """The KLD primers that make the `chosen` changes ((start, end, text) on `old`, as `changes` gives them) in `old`, the sequence
    as saved. Returns SimpleNamespace(rows, notes, product, change, skipped): rows [(name, sequence, annealing length)], notes
    [(level, text)] ("warn" is the orange one), the sequence the primers make, the stretch replaced (start, end, text) and []. Bases
    between chosen changes stay as they are in `old`. Raises ValueError when it cannot be done."""
    old = old.upper()
    if not circular:
        raise ValueError("KLD needs a circular sequence: the product is closed into a plasmid.")
    if len(old) < MIN_LENGTH:
        raise ValueError(f"The sequence is too short for KLD (under {MIN_LENGTH} bp).")
    if not chosen:
        raise ValueError("No change selected: edit the sequence (mutate, insert or delete bases) and tick a change.")
    chosen = sorted(chosen)
    n = len(old)
    # The chosen changes are close if they are close round the circle: the origin sits in the biggest stretch between them.
    gaps = [(chosen[i + 1][0] - chosen[i][1], chosen[i][1]) for i in range(len(chosen) - 1)]
    wrap = chosen[0][0] + n - chosen[-1][1]
    if len(chosen) > 1 and n - max(max(g for g, _ in gaps), wrap) > group:
        raise ValueError(f"The selected changes span {n - max(max(g for g, _ in gaps), wrap):,} bp (more than {group}): untick some, "
                         "or raise \"Group changes within\".")
    linear = old   # the sequence after the changes, bases as they stand in the file: the origin stays where it is
    for s_, e_, t_ in reversed(chosen):
        linear = linear[:s_] + t_.upper() + linear[e_:]
    spin = 0
    if gaps and max(g for g, _ in gaps) > wrap:
        g, at = max(gaps)
        spin = at + g // 2
    old = old[spin:] + old[:spin]   # (designed in a frame where the changes do not meet the origin)
    chosen = sorted(((s_ - spin) % n, (s_ - spin) % n + e_ - s_, t_) for s_, e_, t_ in chosen)
    a, b = chosen[0][0], chosen[-1][1]
    text, at = "", a
    for s_, e_, t_ in chosen:   # the new text of old[a:b]
        text += old[at:s_] + t_.upper()
        at = e_
    if len(text) > MAX_INSERT:
        raise ValueError(f"{len(text):,} new bp are too many for primers (the limit is {MAX_INSERT} bp).")
    notes = []
    if len(text) > warn_insert:
        notes.append(("warn", f"{len(text):,} new bp (over {warn_insert}): the primers are long and may be hard or expensive to order."))
    if not text:
        notes.append(("info", f"Deletes {b - a:,} bp."))
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
    expected = linear
    if len(product) != len(expected) or product not in expected + expected:
        raise ValueError("The primers do not rebuild the edited sequence (unexpected: please report it).")
    for name, seq, an in rows:
        for label, level, what in oligo_checks((name, seq, an, "")):
            if level in ("warn", "bad"):
                notes.append(("warn", f"{name}: {label} {what}"))
        if len(seq) > LONG_PRIMER:
            notes.append(("warn", f"{name} is {len(seq)} nt long."))
    return SimpleNamespace(rows=rows, notes=notes, product=expected, change=((a + spin) % n, b + spin if b + spin <= n else b + spin - n, text), skipped=[])


# ---- the page of the Assembly panel -------------------------------------------------------------------------


def describe(old, c, warn):
    """One change in words, for the list."""
    s, e, t = c
    show = lambda x: x if len(x) <= 12 else x[:12] + "…"
    if s == e:
        what = f"Insertion at {s + 1:,}: +{len(t):,} bp {show(t)}"
    elif not t:
        what = f"Deletion {s + 1:,}–{e:,}: −{e - s:,} bp"
    elif len(t) == e - s:
        what = f"Mutation at {s + 1:,}: {show(old[s:e])} → {show(t)}"
    else:
        what = f"Replacement {s + 1:,}–{e:,}: {e - s:,} → {len(t):,} bp"
    return what + (" ⚠" if len(t) > warn else "")


class KldPage(QWidget):
    """The KLD method of the Assembly panel. Left: the open sequence, editable (the edits are the changes). Middle: the changes
    KLD can make, ticked or not. Right: the sequence after the ticked changes, read only, with the primers drawn on it.
    Everything follows the edits and the ticks."""

    def __init__(self, panel):
        super().__init__()
        from assembly_panel import FragmentEditor, ReadOnlyView
        self.panel, self.win = panel, panel.win
        self._editor_class, self._view_class = FragmentEditor, ReadOnlyView
        self.editor = self.left = self.right = self.right_editor = None
        self.found, self.old, self.manual, self.result, self.product_doc = [], "", {}, None, None
        self.tm = QDoubleSpinBox(minimum=45, maximum=72, value=60, suffix=" °C")
        self.group = QSpinBox(minimum=1, maximum=MAX_INSERT, value=30, suffix=" bp")
        self.warn = QSpinBox(minimum=1, maximum=MAX_INSERT, value=120, suffix=" bp")
        self.group.setToolTip("Changes closer than this are ticked together and made by the same primer pair.")
        self.warn.setToolTip(f"More new bases than this get an orange warning (primers may be hard or expensive to order); "
                             f"over {MAX_INSERT} bp the change is refused.")
        self.method = QComboBox()   # the same list as the panel's: the way back to the other methods
        self.method.setModel(panel.method.model())
        self.method.setCurrentIndex(panel.method.currentIndex())
        self.method.activated.connect(panel.method.setCurrentIndex)
        panel.method.currentIndexChanged.connect(self.method.setCurrentIndex)
        top = QHBoxLayout()
        for label, w in (("Method", self.method), ("Annealing Tm", self.tm), ("Group changes within", self.group),
                         ("Warn for insertions above", self.warn)):
            top.addWidget(QLabel(label))
            top.addWidget(w)
        top.addStretch()
        top.addWidget(QPushButton("Copy primers", clicked=self.copy_primers))
        self.save = QPushButton("Save as new file…", clicked=self.save_file,
                                toolTip="Write the sequence after the ticked changes, with its annotations and the two primers, to a new GenBank file")
        top.addWidget(self.save)
        self.changes = QListWidget()
        self.changes.itemChanged.connect(self._ticked)
        self.changes.itemClicked.connect(self._show_change)
        self.info = QTextBrowser()
        middle = QSplitter(Qt.Vertical)
        middle.addWidget(self.changes)
        middle.addWidget(self.info)
        middle.setSizes([260, 340])
        self.left_title = QLabel("<b>Current sequence</b> (editable)")
        self.right_title = QLabel("<b>After KLD</b> (read only), with the primers")
        cols, boxes = [], []
        for title in (self.left_title, QLabel("<b>Changes KLD can make</b>"), self.right_title):
            w, box = QWidget(), QVBoxLayout()
            box.setContentsMargins(0, 0, 0, 0)
            box.addWidget(title)
            w.setLayout(box)
            cols.append(w)
            boxes.append(box)
        self.left_box, self.right_box = boxes[0], boxes[2]
        boxes[1].addWidget(middle)
        split = QSplitter(Qt.Horizontal)
        for w in cols:
            split.addWidget(w)
        split.setSizes([450, 260, 450])
        lay = QVBoxLayout(self)
        lay.addLayout(top)
        lay.addWidget(split, 1)
        self.timer = QTimer(self, singleShot=True, interval=150, timeout=self.refresh)   # edits come fast: one refresh after a pause
        for w in (self.tm, self.group, self.warn):
            w.valueChanged.connect(self.refresh)

    def showEvent(self, e):
        super().showEvent(e)
        self.refresh()

    def copy_primers(self):
        if self.result:
            QGuiApplication.clipboard().setText("".join(f"{n}\t{s}\n" for n, s, _ in self.result.rows))

    # -- the left view: the open sequence of the window, so editing it here edits the file
    def _follow(self, ed):
        if ed is self.editor:
            return
        if self.editor is not None:
            try:
                self.editor.changed.disconnect(self.timer.start)
            except (RuntimeError, TypeError):
                pass
        if self.left is not None:
            self.left.setParent(None)
            self.left.deleteLater()
            self.left = None
        self.editor, self.manual = ed, {}
        if ed is not None:
            from sequence_view import SequenceView
            self.left = SequenceView(ed)
            self.left_box.addWidget(self.left, 1)
            ed.changed.connect(self.timer.start)

    def _show_right(self, widget):
        if self.right is not None:
            self.right.setParent(None)
            self.right.deleteLater()
        self.right = widget
        self.right_box.addWidget(widget, 1)

    def _say(self, text, error=True):
        """Nothing to show on the right: say why (in the changes' place too)."""
        self.right_editor = self.product_doc = None
        self.save.setEnabled(False)
        note = QLabel(text, wordWrap=True, alignment=Qt.AlignCenter)
        self._show_right(note)
        self.info.setHtml(f"<p style='color:{'#c62828' if error else 'gray'}'>{escape(text)}</p>")

    # -- the changes and the design
    def refresh(self, *_):
        ed = getattr(self.win, "editor", None)
        self._follow(ed)
        self.result = None
        self.changes.blockSignals(True)
        self.changes.clear()
        self.changes.blockSignals(False)
        if ed is None:
            self.left_title.setText("<b>Current sequence</b> (editable)")
            return self._say("Open a sequence, then edit it here (mutate, insert or delete bases).", error=False)
        self.left_title.setText(f"<b>{escape(ed.doc.name)}</b> (editable)")
        self.old = ed.doc.saved_sequence()
        if self.old is None:
            return self._say("There is no saved sequence to compare with (a new sequence, or one rotated or reverse-complemented "
                             "since the last save): save it, then edit.")
        self.found = changes(self.old.upper(), ed.doc.sequence().upper())
        auto_keys = set(auto(self.found, self.group.value(), len(self.old) if self.editor.circular else None))
        self.changes.blockSignals(True)
        for c in self.found:
            item = QListWidgetItem(describe(self.old, c, self.warn.value()))
            item.setData(Qt.UserRole, c)
            item.setFlags(item.flags() | Qt.ItemIsUserCheckable)
            item.setCheckState(Qt.Checked if self.manual.get(c, c in auto_keys) else Qt.Unchecked)
            if len(c[2]) > self.warn.value():
                item.setForeground(QColor("#ef6c00"))
            self.changes.addItem(item)
        self.changes.blockSignals(False)
        self._design()

    def _ticked(self, item):
        self.manual[item.data(Qt.UserRole)] = item.checkState() == Qt.Checked
        self._design()

    def _show_change(self, item):
        """Click a change: select it on the left (its place in the current sequence)."""
        s, e, t = item.data(Qt.UserRole)
        shift = sum(len(x[2]) - (x[1] - x[0]) for x in self.found if x[0] < s)
        if self.editor is not None and self.editor.length:
            self.editor.select(s + shift, min(max(s + shift + len(t), s + shift + 1), self.editor.length))
            self.left.ensure_visible(s + shift)

    def _product(self, r):
        """The sequence after the ticked changes as a Document: the open sequence with its annotations (they follow the edits), the
        unticked changes put back, and the two primers added where they bind."""
        from assembly import snapshot, transfer_primers
        doc = snapshot(self.editor.doc, self.editor.doc.name)
        shift = 0
        spots = []   # (start, end, old text) of every change in the current sequence
        for s, e, t in self.found:
            spots.append((s + shift, s + shift + len(t), self.old[s:e], (s, e, t)))
            shift += len(t) - (e - s)
        ticked = set(self._checked())
        for cs, ce, was, c in reversed(spots):   # last first: the places before it do not move
            if not self.manual.get(c, c in ticked):
                doc.replace(cs, ce, was)
        if doc.sequence().upper() != r.product:   # (never expected: the plain sequence still shows the primers)
            from opendna_core import Document
            doc = Document.new_sequence(self.editor.doc.name, r.product, True)
        transfer_primers(doc, [SimpleNamespace(fwd=(n, s) if n.endswith("F") else None, rev=(n, s) if n.endswith("R") else None)
                               for n, s, _ in r.rows])
        return doc

    def _checked(self):
        """The ticked changes of the list, in order."""
        return [self.changes.item(i).data(Qt.UserRole) for i in range(self.changes.count())
                if self.changes.item(i).checkState() == Qt.Checked]

    def save_file(self):
        """Save as new file…: the product (annotations and primers included) as a GenBank file, next to the open one."""
        if self.product_doc is None:
            return
        from PySide6.QtWidgets import QFileDialog, QMessageBox
        ed = self.editor
        folder = ed.path.parent if not ed.untitled else Path.cwd()
        path, _ = QFileDialog.getSaveFileName(self, "Save the sequence after KLD", str(folder / f"{ed.doc.name}_KLD.gb"),
                                              "GenBank (*.gb *.gbk)")
        if path:
            try:
                self.product_doc.write_genbank(path, None)
            except OSError as e:
                QMessageBox.critical(self, "Save Failed", f"{path}\n\n{e}")
            else:
                self.win.statusBar().showMessage(f"Saved {path}", 5000)

    def _design(self):
        try:
            self.result = r = design_selected(self.old, self._checked(), self.editor.circular, self.tm.value(), self.group.value(),
                                              self.warn.value())
        except ValueError as e:
            self.result = None
            return self._say(str(e))
        a, b, text = r.change
        keep = self.right.verticalScrollBar().value() if isinstance(self.right, self._view_class) else 0
        self.product_doc = doc = self._product(r)
        self.save.setEnabled(True)
        self.right_editor = ed = self._editor_class(doc, "After KLD")
        view = self._view_class(ed)
        self._show_right(view)
        if text and a < b:   # (a region across the origin is not selected as a whole)
            ed.select(a, a + len(text))
        else:
            ed.set_cursor(a)
        view.ensure_visible(a)
        if keep:
            view.verticalScrollBar().setValue(keep)
        html = ["<table cellspacing=6>"]
        for name, seq, n in r.rows:
            tail = len(seq) - n
            t = melting_temp(seq[tail:])
            html.append(f"<tr><td><b>{name}</b></td><td style='font-family:monospace'><span style='color:{TAIL_COLOR}'>"
                        f"{seq[:tail]}</span>{seq[tail:]}</td><td>{len(seq)} nt</td><td>{f'{t:.1f} °C' if t else ''}</td></tr>")
        html.append("</table>")
        for level, note in r.notes:
            html.append(f"<p style='color:{'#ef6c00' if level == 'warn' else 'gray'}'>{'⚠ ' if level == 'warn' else ''}{escape(note)}</p>")
        self.info.setHtml("".join(html))
