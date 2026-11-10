"""Codon tools for moving a gene into a new host: recode a CDS for a codon usage table, and read it under another genetic code.

Tools > Recode CDS for a Codon Usage Table… rewrites the codons of a selected CDS (or of a selected stretch of whole codons) so that
the protein stays exactly the same while the DNA suits the new host. Tools > Compare CDS Under Another Genetic Code… shows the
protein under two genetic codes side by side, stops included: it tells you whether a gene survives a move into an amber-suppressor
strain, a mitochondrion or Mycoplasma (where TGA is Trp, not a stop).

The recoding is plain Python (`recode`, `compare`): no GUI is needed to use it. The dialogs write the result back with
`Editor.overwrite_bases`, one normal edit: undo, redo and the features all work as for any other edit.
"""
import random
import re
from html import escape
from types import SimpleNamespace

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QCheckBox, QComboBox, QDialogButtonBox, QFormLayout, QHBoxLayout, QLabel, QListWidget, QListWidgetItem, QPushButton,
    QSpinBox, QTextBrowser, QVBoxLayout,
)

import codon_usage
import genetic_code
import theme
from editor import mono_font, reverse_complement
from enzymes import COMMON, site_spans
from tools import ToolDialog
from uikit import no_auto_default, plural

RBS = "AGGAGG"          # the Shine-Dalgarno core: inside a CDS it makes the ribosome start (or stall) in the wrong place
MAX_RUN = 5             # no more than this many identical bases in a row
GC_WINDOW = 50          # the GC limits hold for every stretch of this many bases
CONTEXT = 12            # bases of the neighbours read too, so a site made across the CDS's edge is seen (longest site is 8)
MAX_PASSES = 30         # rounds of repair before the remaining problems are reported as unmet
MATCHED, HIGHEST = "Frequency matched", "Highest frequency"
BASES = set("ACGT")


# ---- the CDS to work on ---------------------------------------------------------------------------------------------------
def cds_source(editor):
    """The coding sequence the selection stands for: SimpleNamespace(start, end, strand, seq, left, right, table, frame, label)
    with `seq` read 5'→3' on the coding strand, `left`/`right` the bases on either side of it (same orientation), or a sentence
    saying why there is none (a string). A selected CDS feature is read as annotated (its table, /codon_start); bare DNA selected
    is read from its first base on the top strand with the table chosen in Settings."""
    doc = editor.doc
    f = editor.feature_by_index(editor.current_feature) if editor.current_feature is not None else None
    if f is not None and f.kind == "CDS":
        if len(f.spans) != 1 or editor.wraps:
            return "This CDS is made of several parts (a join, or it crosses the origin): select one piece of DNA to work on instead."
        (start, end), strand, label = f.spans[0], f.strand, f.label
        plan = genetic_code.cds_plan(doc, f)
        table, frame = plan.table, plan.frame
        if plan.exceptions or plan.why:
            return "This CDS has a /transl_except or an annotation openDNA cannot follow in full (see its tooltip)."
    elif editor.end > editor.start and not editor.wraps:
        start, end, strand, label = editor.start, editor.end, 1, "Selection"
        table, frame = genetic_code.table_id(), 0
    else:
        return "Select a CDS feature, or a stretch of DNA that is a whole number of codons."
    seq = doc.sequence(start, end)
    left, right = doc.sequence(max(0, start - CONTEXT), start), doc.sequence(end, min(editor.length, end + CONTEXT))
    if strand < 0:
        seq, left, right = reverse_complement(seq), reverse_complement(right), reverse_complement(left)
    seq = seq[frame:]
    if len(seq) % 3 or not seq:
        return f"{label} is {len(seq):,} bases: not a whole number of codons."
    if not set(seq.upper()) <= BASES:
        return f"{label} holds bases other than A, C, G and T: fix them first."
    return SimpleNamespace(start=start + (frame if strand > 0 else 0), end=end - (0 if strand > 0 else frame), strand=strand,
                           seq=seq.upper(), left=left.upper(), right=right.upper(), table=table, frame=frame, label=label)


# ---- recoding ---------------------------------------------------------------------------------------------------------------
def translate(seq, code):
    """The protein of `seq` under `code` (a codon → letter dict): * for a stop, X for what the code does not know."""
    return "".join(code.get(seq[i:i + 3], "X") for i in range(0, len(seq) - len(seq) % 3, 3))


def gc(seq):
    return 100 * sum(c in "GC" for c in seq) / max(len(seq), 1)


def problems(text, lo, hi, avoid, max_run, rbs, gc_range, window):
    """What is wrong in the part [lo, hi) of `text` (the rest of `text` is only the neighbours, read so that a site across the edge
    is seen): [(what, start, end)] with positions in `text`. `avoid` are Enzyme tuples, `gc_range` (low, high) per cent or None."""
    out = []
    for e in avoid:
        out += [(f"{e.name} site", a, b) for a, b, _ in site_spans(text, False, e) if a < hi and b > lo]
    for m in re.finditer(r"A+|C+|G+|T+", text):
        if m.end() - m.start() > max_run and m.start() < hi and m.end() > lo:
            out.append((f"run of {m.end() - m.start()} {m.group()[0]}", m.start(), m.end()))
    if rbs:
        out += [("ribosome binding site", m.start(), m.start() + len(RBS)) for m in re.finditer(f"(?={RBS})", text)
                if m.start() < hi and m.start() + len(RBS) > lo]
    if gc_range:
        out += _gc_problems(text[lo:hi], gc_range, window, lo)
    return out


def _gc_problems(seq, gc_range, window, offset):
    """One problem per stretch of consecutive windows outside the limits."""
    low, high = gc_range
    n = len(seq)
    if n <= window:
        return [] if low <= gc(seq) <= high else [("GC content", offset, offset + n)]
    prefix = [0]
    for c in seq:
        prefix.append(prefix[-1] + (c in "GC"))
    out, run = [], None
    for i in range(n - window + 1):
        bad = not low <= 100 * (prefix[i + window] - prefix[i]) / window <= high
        if bad and run is None:
            run = i
        if run is not None and (not bad or i == n - window):
            out.append(("GC content", offset + run, offset + i + (window if bad else window - 1)))
            run = None
    return out


def _score(text, lo, hi, *rules):
    """How much is wrong: the number of problems, with a GC stretch counting by its length so a repair that shortens it scores better."""
    return sum((b - a) // 10 + 1 if what == "GC content" else 1 for what, a, b in problems(text, lo, hi, *rules))


def recode(seq, code, table, method=MATCHED, avoid=(), gc_range=None, window=GC_WINDOW, max_run=MAX_RUN, rbs=True,
           left="", right="", seed=0):
    """Rewrite the codons of `seq` (a whole number of codons, A/C/G/T) for the codon usage `table`, keeping the protein that `code`
    (a codon → letter dict) reads in it. Returns SimpleNamespace(seq, codons, changed, gc_before, gc_after, unmet): `unmet` is a
    list of sentences for the rules that could not be met (empty when all were). `left` and `right` are the neighbouring bases,
    used only to see problems across the CDS's edges.

    method: MATCHED picks each codon at random in proportion to its frequency in `table` (no long run of the best codon, no
    depleted tRNA); HIGHEST always takes the most used synonym. A codon the table never uses is not chosen unless it is the
    only choice. The same `seed` gives the same result. The rules are then enforced by changing one codon at a time, always
    to a synonym, keeping a change only when it leaves fewer problems."""
    if len(seq) % 3 or not set(seq) <= BASES:
        raise ValueError("A coding sequence is a whole number of codons made of A, C, G and T.")
    rng = random.Random(seed)
    freq = table.per_thousand
    synonyms = {}
    for c, aa in code.items():
        synonyms.setdefault(aa, []).append(c)
    old = [seq[i:i + 3] for i in range(0, len(seq), 3)]
    # a codon the code reads as X (removed in a recoded organism) has no synonym to offer: it stays as it is
    options = [sorted((c for c in synonyms.get(code.get(o, "X"), []) if freq[c] > 0), key=lambda c: -freq[c]) or [o] for o in old]
    new = [opts[0] if method == HIGHEST or len(opts) == 1 else rng.choices(opts, [freq[c] for c in opts])[0] for opts in options]
    rules = (avoid, max_run, rbs, gc_range, window)
    lo = len(left)
    full = lambda cods: left + "".join(cods) + right
    for _ in range(MAX_PASSES):
        bad = problems(full(new), lo, lo + len(seq), *rules)
        if not bad:
            break
        best, improved = _score(full(new), lo, lo + len(seq), *rules), False
        for _what, a, b in bad:
            first, last = max(0, (a - lo) // 3), min(len(new) - 1, (b - 1 - lo) // 3)
            for i in rng.sample(range(first, last + 1), last - first + 1):
                for alt in options[i]:
                    if alt == new[i]:
                        continue
                    trial = new[:i] + [alt] + new[i + 1:]
                    score = _score(full(trial), lo, lo + len(seq), *rules)
                    if score < best:
                        new, best, improved = trial, score, True
                        break
        if not improved:
            break
    result = "".join(new)
    if translate(result, code) != translate(seq, code):   # the one thing that must never be wrong
        raise RuntimeError("Recoding changed the protein: nothing was written.")
    unmet = [f"{what} at bases {a - lo + 1:,}–{b - lo:,}" for what, a, b in problems(full(new), lo, lo + len(seq), *rules)]
    return SimpleNamespace(seq=result, codons=len(old), changed=sum(a != b for a, b in zip(old, new)),
                           gc_before=gc(seq), gc_after=gc(result), unmet=unmet)


# ---- the same CDS under two genetic codes -----------------------------------------------------------------------------------
def compare(seq, code_a, code_b):
    """The protein of `seq` under two codes: SimpleNamespace(a, b, differ [codon numbers from 1], gained [..], lost [..]).
    `gained`: codons that are an internal stop only under `code_b` (the gene would be cut short there), `lost`: internal stops of
    `code_a` that `code_b` reads as an amino acid (read-through). The last codon is not internal: a final stop is expected."""
    a, b = translate(seq, code_a), translate(seq, code_b)
    last = len(a) - 1
    differ = [i + 1 for i, (x, y) in enumerate(zip(a, b)) if x != y]
    return SimpleNamespace(a=a, b=b, differ=differ,
                           gained=[i for i in differ if b[i - 1] == "*" and i - 1 != last],
                           lost=[i for i in differ if a[i - 1] == "*" and i - 1 != last])


def comparison_html(res, name_a, name_b, per_line=50, block=10):
    """The two proteins as rows of `per_line` letters in groups of `block`, numbered like the Protein view; a letter that
    differs is bold and coloured, a stop is drawn *. Returns HTML for a monospace text browser."""
    mark = theme.hex("cut")
    width = max(len(name_a), len(name_b), 6)
    rows = []
    for start in range(0, len(res.a), per_line):
        lines = []
        for label, text in ((name_a, res.a), (name_b, res.b), ("", None)):
            cells = []
            for i in range(start, min(start + per_line, len(res.a))):
                if text is None:
                    cells.append(f'<span style="color:{mark}">▲</span>' if (i + 1) in res.differ else " ")
                else:
                    ch = escape(text[i])
                    cells.append(f'<b style="color:{mark}">{ch}</b>' if (i + 1) in res.differ else ch)
                if (i + 1 - start) % block == 0 and i + 1 < min(start + per_line, len(res.a)):
                    cells.append(" ")
            lines.append(f"{escape(label):<{width}} " + "".join(cells) + (f"  {min(start + per_line, len(res.a))}" if label == name_a else ""))
        rows.append("\n".join(lines))
    return f'<pre style="margin:0">{(chr(10) * 2).join(rows)}</pre>'


# ---- the dialogs ------------------------------------------------------------------------------------------------------------
class _CdsDialog(ToolDialog):
    """What both tools share: the CDS read from the selection, named in a line at the top, or the reason there is none."""

    def __init__(self, parent, editor, title):
        super().__init__(parent)
        self.editor = editor
        self.setWindowTitle(title)
        self.src = cds_source(editor)
        self.layout_ = QVBoxLayout(self)
        self.head = QLabel(wordWrap=True)
        self.layout_.addWidget(self.head)
        if isinstance(self.src, str):
            self.head.setText(self.src)
            return
        s = self.src
        strand = "−" if s.strand < 0 else "+"
        self.head.setText(f"<b>{escape(s.label)}</b> · {len(s.seq):,} bp · {plural(len(s.seq) // 3, 'codon')} · strand {strand} · "
                          f"{genetic_code.name_of(s.table)}")

    @property
    def usable(self):
        return not isinstance(self.src, str)


class RecodeDialog(_CdsDialog):
    def __init__(self, parent, editor):
        super().__init__(parent, editor, "Recode CDS for a Codon Usage Table")
        self.result, self.attempt = None, 0
        if not self.usable:
            self._buttons(has_apply=False)
            return
        self.setMinimumWidth(560)
        self.tables = codon_usage.all_tables()
        if not self.tables:
            self.head.setText(self.head.text() + "<br>There is no codon usage table yet: import or make one in Tools > Settings "
                              "> Codon usage, then come back.")
            self._buttons(has_apply=False)
            return
        self.usage = QComboBox()
        self.usage.addItems(list(self.tables))
        active = codon_usage.active()
        if active:
            self.usage.setCurrentText(active.name)
        self.method = QComboBox()
        self.method.addItems([MATCHED, HIGHEST])
        self.method.setToolTip("Frequency matched picks each codon in proportion to the table: no long run of the best codon and "
                               "no tRNA depletion. Highest frequency always takes the most used codon.")
        self.enzymes = QListWidget(maximumHeight=110)
        self.enzymes.setToolTip("Tick the restriction enzymes whose sites must not be in the recoded gene.")
        for e in COMMON:
            item = QListWidgetItem(f"{e.name}  {e.site}")
            item.setFlags(item.flags() | Qt.ItemIsUserCheckable)
            item.setCheckState(Qt.Unchecked)
            item.setData(Qt.UserRole, e)
            self.enzymes.addItem(item)
        self.gc_on = QCheckBox("Keep GC between", checked=True)
        self.gc_low = QSpinBox(minimum=0, maximum=100, value=35, suffix=" %")
        self.gc_high = QSpinBox(minimum=0, maximum=100, value=65, suffix=" %")
        self.window = QSpinBox(minimum=10, maximum=1000, value=GC_WINDOW, suffix=" bases")
        gc_row = QHBoxLayout()
        for w in (self.gc_on, self.gc_low, QLabel("and"), self.gc_high, QLabel("in every"), self.window):
            gc_row.addWidget(w)
        gc_row.addStretch()
        self.run = QSpinBox(minimum=3, maximum=20, value=MAX_RUN, suffix=" identical bases")
        self.rbs = QCheckBox(f"No internal ribosome binding site ({RBS})", checked=True)
        form = QFormLayout()
        form.addRow("Codon usage table", self.usage)
        form.addRow("Method", self.method)
        form.addRow("Avoid sites", self.enzymes)
        form.addRow("", gc_row)
        form.addRow("No run of more than", self.run)
        form.addRow("", self.rbs)
        self.layout_.addLayout(form)
        self.preview = QPushButton("Preview recoding", clicked=self.run_recode)
        self.layout_.addWidget(self.preview)
        self.report = QLabel(wordWrap=True, textFormat=Qt.RichText)
        self.report.setMinimumHeight(70)
        self.report.setAlignment(Qt.AlignTop)
        self.layout_.addWidget(self.report)
        for w in (self.usage, self.method):
            w.currentIndexChanged.connect(self._stale)
        for w in (self.gc_low, self.gc_high, self.window, self.run):
            w.valueChanged.connect(self._stale)
        for w in (self.gc_on, self.rbs):
            w.toggled.connect(self._stale)
        self.enzymes.itemChanged.connect(self._stale)
        self._buttons(has_apply=True)
        no_auto_default(self)

    def _buttons(self, has_apply):
        box = QDialogButtonBox(QDialogButtonBox.Close)
        if has_apply:
            self.apply_button = box.addButton("Apply to sequence", QDialogButtonBox.AcceptRole)
            self.apply_button.setEnabled(False)
        box.accepted.connect(self.apply)
        box.rejected.connect(self.reject)
        self.layout_.addWidget(box)

    def _stale(self, *_):
        """Any option changed: the preview no longer says what Apply would do."""
        self.result = None
        self.apply_button.setEnabled(False)

    def rules(self):
        avoid = [self.enzymes.item(i).data(Qt.UserRole) for i in range(self.enzymes.count())
                 if self.enzymes.item(i).checkState() == Qt.Checked]
        return dict(method=self.method.currentText(), avoid=avoid, max_run=self.run.value(), rbs=self.rbs.isChecked(),
                    gc_range=(self.gc_low.value(), self.gc_high.value()) if self.gc_on.isChecked() else None,
                    window=self.window.value())

    def run_recode(self):
        s = self.src
        if self.gc_on.isChecked() and self.gc_low.value() > self.gc_high.value():
            self.report.setText("The lowest GC is above the highest.")
            return
        self.attempt += 1   # each press of the button draws a new sample (frequency matched)
        table = self.tables[self.usage.currentText()]
        try:
            self.result = recode(s.seq, genetic_code.codon_map(s.table), table, left=s.left, right=s.right, seed=self.attempt,
                                 **self.rules())
        except (ValueError, RuntimeError) as e:
            self.report.setText(str(e))
            return
        r = self.result
        text = (f"<b>{plural(r.changed, 'codon')}</b> of {r.codons:,} changed · GC {r.gc_before:.1f} % → {r.gc_after:.1f} % · "
                f"the protein is identical.")
        if r.unmet:
            shown = "; ".join(r.unmet[:6]) + (f"; and {len(r.unmet) - 6} more" if len(r.unmet) > 6 else "")
            text += f"<br><b>Could not be met:</b> {escape(shown)}."
        self.report.setText(text)
        self.apply_button.setEnabled(True)

    def apply(self):
        """Write the recoded bases into the sequence: one edit, so one Undo restores them."""
        r, s = self.result, self.src
        if r is None:
            return
        bases = r.seq if s.strand > 0 else reverse_complement(r.seq)
        if self.editor.overwrite_bases(s.start, bases):
            self.accept()


class CompareDialog(_CdsDialog):
    def __init__(self, parent, editor):
        super().__init__(parent, editor, "Compare CDS Under Another Genetic Code")
        box = QDialogButtonBox(QDialogButtonBox.Close)
        box.rejected.connect(self.reject)
        if not self.usable:
            self.layout_.addWidget(box)
            return
        self.resize(760, 520)
        self.code_a, self.code_b = QComboBox(), QComboBox()
        for combo in (self.code_a, self.code_b):
            for t in genetic_code.TABLES.values():
                combo.addItem(f"{t.id} – {t.name}", t.id)
        self.code_a.setCurrentIndex(self.code_a.findData(self.src.table))
        self.code_b.setCurrentIndex(self.code_b.findData(4 if self.src.table != 4 else 1))   # 4: Mycoplasma, TGA = Trp
        top = QHBoxLayout()
        for text, w in (("Read as", self.code_a), ("and as", self.code_b)):
            top.addWidget(QLabel(text))
            top.addWidget(w, 1)
        self.summary = QLabel(wordWrap=True, textFormat=Qt.RichText)
        self.view = QTextBrowser()
        self.view.setFont(mono_font(10))
        self.layout_.addLayout(top)
        self.layout_.addWidget(self.summary)
        self.layout_.addWidget(self.view, 1)
        self.layout_.addWidget(box)
        no_auto_default(self)
        for combo in (self.code_a, self.code_b):   # only now: the slot needs the widgets above
            combo.currentIndexChanged.connect(self.show_comparison)
        self.show_comparison()

    def show_comparison(self, *_):
        a, b = self.code_a.currentData(), self.code_b.currentData()
        res = compare(self.src.seq, genetic_code.codon_map(a), genetic_code.codon_map(b))
        if not res.differ:
            text = "The protein is the same under both codes."
        else:
            text = f"<b>{plural(len(res.differ), 'residue')}</b> differ."
            if res.gained:
                text += (f" <b>Under code {b} the gene gains an internal stop</b> at {plural(len(res.gained), 'codon')} "
                         f"({', '.join(map(str, res.gained[:8]))}): it would be cut short.")
            if res.lost:
                text += f" Code {b} reads {plural(len(res.lost), 'internal stop')} of code {a} as an amino acid (read-through)."
        self.summary.setText(text)
        self.view.setHtml(comparison_html(res, f"Code {a}", f"Code {b}"))
