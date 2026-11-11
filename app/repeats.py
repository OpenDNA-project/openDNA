"""Repeats: stretches of a construct that occur more than once, and the dot plot of a sequence against another (or itself).

It is the paired view of Local Align: with two sequences in the boxes, `pair_repeats` lists the stretches they share and
`pair_dots` / `DotPlot` draw them (A across, B down); the same sequence in both boxes finds the repeats of one sequence.

Why it matters: two copies of the same 20+ bases let the cell's own recombination delete what lies between them, an inverted
repeat folds into a hairpin that stalls polymerases, and synthesis companies refuse sequences full of either.

Engine (no GUI): `find_repeats` lists direct repeats (the same stretch twice, same orientation), inverted repeats (a stretch
and its reverse complement, arms apart) and tandem arrays (a unit repeated head to tail, three times or more), each with all
its copies. `pair_dots` gives the dots of the plot. Both work the way core/src/align.rs does: every k-mer of the sequence is
put in a dictionary, a k-mer seen twice is a seed, and a seed is stretched along its diagonal while mismatches stay within the
allowance. Pure Python, no numpy (it is not a dependency of the program): 10 kb takes well under a second, 100 kb a few.

Work is capped (MAX_LENGTH, MAX_OCCURRENCES, MAX_SEEDS, MAX_REPEATS); `Result.capped` says in words what was left out, so a
capped answer is never mistaken for "no repeats".
"""
from collections import defaultdict, namedtuple

from PySide6.QtCore import QPointF, QRectF, Qt, Signal
from PySide6.QtGui import QColor, QFontMetrics, QImage, QPainter, QPen
from PySide6.QtWidgets import (
    QAbstractItemView, QHBoxLayout, QLabel, QPushButton, QSizePolicy, QSpinBox, QSplitter, QTreeWidget, QTreeWidgetItem,
    QVBoxLayout, QWidget,
)

import theme
from editor import reverse_complement
from genome_view import nice_step
from uikit import fit_columns, plural, run_busy

MAX_LENGTH = 1_000_000     # bp looked at: a longer sequence is cut here (and the result says so)
MAX_OCCURRENCES = 200      # a k-mer seen more often than this is skipped (poly-A, microsatellites: they would be most of the work)
MAX_SEEDS = 2_000_000      # seed pairs examined in one search
MAX_REPEATS = 500          # repeats listed
MAX_DOTS = 300_000         # dots drawn
MAX_UNIT = 200             # longest unit of a tandem array
MIN_SEED = 8               # shortest seed: below it random matches swamp the real ones
WRAP = 1000                # on a circle, a repeat that runs over the origin is found if it is shorter than this

# One repeat: kind ("direct", "inverted", "tandem"), length (of one copy), copies [(start, end, strand)] (0-based, end exclusive,
# on a circle end may pass the length; strand -1 = the copy reads as the reverse complement of the first), mismatches (the
# most between any two copies seen, 0 = exact), unit (tandem only: the unit's length).
Repeat = namedtuple("Repeat", "kind length copies mismatches unit")
Result = namedtuple("Result", "repeats capped")


def _clean(seq):
    """Upper case, anything but A, C, G, T as N (positions must not shift); a k-mer with an N is never a seed."""
    return "".join(c if c in "ACGT" else "N" for c in seq.upper())


def _seeds(a, b, k, same, budget):
    """{diagonal (position in b - position in a): [position in a, ...]} of every k-mer of `a` that b holds, `same`: a is b (only
    pairs where b is after a). Returns (diagonals, capped note or '')."""
    index = defaultdict(list)
    for q in range(len(b) - k + 1):
        if "N" not in b[q:q + k]:
            index[b[q:q + k]].append(q)
    diagonals, count, note = defaultdict(list), 0, ""
    for i in range(len(a) - k + 1):
        places = index.get(a[i:i + k], ())
        if len(places) > MAX_OCCURRENCES:
            note = "some very repetitive stretches were skipped"
            continue
        for q in places:
            if same and q <= i:
                continue
            count += 1
            if count > budget:
                return diagonals, "the search was stopped at its limit of seed matches"
            diagonals[q - i].append(i)
    return diagonals, note


def _segments(a, b, d, starts, k, mismatches):
    """Stretches of the diagonal d (a[i] against b[i + d]) around the seeds `starts`: [(i, length, mismatches)], each beginning
    and ending on a match, with at most `mismatches` in all."""
    out, cursor = [], 0
    la, lb = len(a), len(b)
    for s in starts:
        if s < cursor:
            continue
        x, mm, end = s + k, 0, s + k
        while x < la and x + d < lb:
            if a[x] == b[x + d]:
                end = x + 1
            else:
                mm += 1
                if mm > mismatches:
                    break
            x += 1
        start, x = s, s - 1
        mm = sum(a[y] != b[y + d] for y in range(s + k, end))   # the left side gets what the right side left over
        while x >= cursor and x + d >= 0:
            if a[x] == b[x + d]:
                start = x
            else:
                mm += 1
                if mm > mismatches:
                    break
            x -= 1
        cursor = end
        out.append((start, end - start, sum(a[y] != b[y + d] for y in range(start, end))))
    return out


def _overlap(a, b, n):
    """Bases two copies (start, length) share, on a circle of n (n huge for a line)."""
    best = 0
    for (s1, l1), (s2, l2) in (((a), (b)), ((b), (a))):
        d = (s2 - s1) % n
        if d < l1:
            best = max(best, min(l1 - d, l2))
    return best


def _group(pairs, n):
    """Copies that are the same repeat into one group: [(length, [(start, strand)], mismatches)]. A pair whose copies overlap
    (by half or more) copies already in a group joins it; the copy it adds that the group lacks is added to the group."""
    groups = []
    for (s1, st1, s2, st2, length, mm) in pairs:
        found = None
        for g in groups:
            for c in g["copies"]:
                if _overlap((c[0], g["length"]), (s1, length), n) * 2 >= min(length, g["length"]) or \
                        _overlap((c[0], g["length"]), (s2, length), n) * 2 >= min(length, g["length"]):
                    found = g
                    break
            if found:
                break
        if not found:
            found = {"length": length, "copies": [], "mm": 0}
            groups.append(found)
        found["mm"] = max(found["mm"], mm)
        found["length"] = max(found["length"], length)
        for start, strand in ((s1, st1), (s2, st2)):
            if not any(_overlap((c[0], found["length"]), (start, length), n) * 2 >= min(length, found["length"]) for c in found["copies"]):
                found["copies"].append((start, strand))
    return groups


def find_repeats(seq, circular=False, min_direct=20, min_inverted=12, min_tandem=20, mismatches=0, max_loop=200):
    """A Result: the repeats of `seq`, direct ones of at least `min_direct` bases, inverted of `min_inverted`, tandem arrays
    of `min_tandem` bases in all (and at least three copies), each allowing `mismatches` differences between its copies.
    The two arms of an inverted repeat are at most `max_loop` bases apart (None: any distance): a stem of 12 pairs up by chance
    somewhere in a long enough sequence, a hairpin needs its arms close. Longest first. Shortcut: with mismatches, only matches that hold an exact stretch of `seed` bases (at least 8) are found."""
    seq = _clean(seq)
    notes = []
    if len(seq) > MAX_LENGTH:
        seq = seq[:MAX_LENGTH]
        notes.append(f"only the first {MAX_LENGTH:,} bp were searched")
    n = len(seq)
    wide = n if circular else 10 ** 12
    t = seq + seq[:min(n - 1, WRAP)] if circular else seq
    seed = max(MIN_SEED, min(min_direct, min_inverted, min_tandem) // (mismatches + 1)) if n else MIN_SEED
    seed = min(seed, min(min_direct, min_inverted))   # (never a seed longer than the shortest repeat asked for)
    if n < 2 * seed:
        return Result([], "; ".join(notes))
    # ---- same orientation: direct repeats and tandem arrays
    diagonals, note = _seeds(t, t, seed, True, MAX_SEEDS)
    notes += [note] if note else []
    direct, tandem = [], []
    for d, starts in diagonals.items():
        if circular and d == n:
            continue   # the sequence against its own copy in the tail
        for i, length, mm in _segments(t, t, d, starts, seed, mismatches):
            if circular and i >= n:
                continue   # entirely in the copied tail: the same as one that starts at the front
            if d >= length and length >= min_direct:
                direct.append((i % n, 1, (i + d) % n, 1, length, mm))
            elif d < length and d <= MAX_UNIT and (length + d) // d >= 3 and length + d >= min_tandem:
                tandem.append((d, i, length + d, mm))
    # tandem arrays, smallest unit first: ACACAC is a unit of 2, not also of 4 and 6
    arrays = []
    for d, i, total, mm in sorted(tandem):
        if any(a0 <= i and i + total <= a1 for _, a0, a1, _ in arrays):
            continue
        arrays.append((d, i, i + total, mm))
    found = []
    for d, i, end, mm in arrays:
        copies = (end - i) // d
        found.append(Repeat("tandem", d, [(i + c * d, i + (c + 1) * d, 1) for c in range(copies)], mm, d))
    # a direct or inverted pair inside one array is the array's own doing
    def inside_array(s1, s2, length):
        return any(a0 <= s1 and s1 + length <= a1 and a0 <= s2 and s2 + length <= a1 for _, a0, a1, _ in arrays)
    direct = [p for p in direct if not inside_array(p[0], p[2], p[4])]
    # ---- opposite orientation: inverted repeats (the sequence against its reverse complement)
    r = reverse_complement(t)
    m = len(t)
    inverted = []
    diagonals, note = _seeds(t, r, seed, False, MAX_SEEDS)
    notes += [note] if note and note not in notes else []
    for d, starts in diagonals.items():
        for i, length, mm in _segments(t, r, d, starts, seed, mismatches):
            j = m - (i + d) - length   # where the matching stretch lies in t, read as its reverse complement
            if length < min_inverted or (circular and (i >= n or j >= n)):
                continue
            if i > j or j < i + length:
                continue   # the same pair seen from the other arm, or arms that overlap (a palindrome, not a hairpin)
            loop = j - (i + length)
            if circular:
                loop = min(loop, n - (j + length) + i)   # (the way round the origin may be the short one)
            if max_loop is not None and loop > max_loop:
                continue
            if not inside_array(i, j, length):
                inverted.append((i % n, 1, j % n, -1, length, mm))
    for kind, pairs in (("direct", direct), ("inverted", inverted)):
        for g in _group(pairs, wide):
            copies = sorted(((s, s + g["length"], st) for s, st in g["copies"]), key=lambda c: c[0])
            found.append(Repeat(kind, g["length"], copies, g["mm"], 0))
    found.sort(key=lambda r: (-r.length, r.copies[0][0]))
    if len(found) > MAX_REPEATS:
        notes.append(f"only the longest {MAX_REPEATS} of {len(found):,} repeats are listed")
        found = found[:MAX_REPEATS]
    return Result(found, "; ".join(notes))


def pair_dots(a, b, k=12):
    """(direct, inverted, note): the dots of `a` (x) against `b` (y): a direct dot is where the window at x equals the window at y, an inverted dot where it equals the reverse complement of the one at y (the window middles, in bases). With the same sequence twice only the
    pairs off the main diagonal are listed (both ways round)."""
    a, b = _clean(a)[:MAX_LENGTH], _clean(b)[:MAX_LENGTH]
    same, half, direct, inverted, notes = a == b, k / 2, [], [], []
    if len(a) < k or len(b) < k:
        return direct, inverted, ""
    diagonals, note = _seeds(a, b, k, same, MAX_SEEDS)
    notes += [note] if note else []
    for d, starts in diagonals.items():
        for i in starts:
            direct.append((i + half, i + d + half))
            if same:
                direct.append((i + d + half, i + half))
    diagonals, _ = _seeds(a, reverse_complement(b), k, False, MAX_SEEDS)
    for d, starts in diagonals.items():
        for i in starts:
            inverted.append((i + half, len(b) - (i + d) - k + half))
    if len(direct) + len(inverted) > MAX_DOTS:
        direct, inverted = direct[:MAX_DOTS // 2], inverted[:MAX_DOTS // 2]
        notes.append(f"only {MAX_DOTS:,} dots are drawn")
    return direct, inverted, "; ".join(notes)


def pair_repeats(a, b, min_len=12, mismatches=0):
    """A Result whose repeats are the stretches `a` and `b` share, (kind, start in a, start in b, length, differences), longest
    first: "direct" (same orientation) and "inverted" (one is the reverse complement of the other; start in b is the left end of
    the stretch on b's top strand). With the same sequence twice each pair is listed once and the main diagonal is left out."""
    a, b = _clean(a)[:MAX_LENGTH], _clean(b)[:MAX_LENGTH]
    seed = min(max(MIN_SEED, min_len // (mismatches + 1)), min_len)
    same, found, notes = a == b, [], []
    for kind, target in (("direct", b), ("inverted", reverse_complement(b))):
        diagonals, note = _seeds(a, target, seed, same and kind == "direct", MAX_SEEDS)
        notes += [note] if note and note not in notes else []
        for d, starts in diagonals.items():
            for i, length, mm in _segments(a, target, d, starts, seed, mismatches):
                j = i + d if kind == "direct" else len(b) - (i + d) - length
                if length >= min_len and not (same and kind == "inverted" and i > j):
                    found.append((kind, i, j, length, mm))
    found.sort(key=lambda r: (-r[3], r[1]))
    if len(found) > MAX_REPEATS:
        notes.append(f"only the longest {MAX_REPEATS} of {len(found):,} stretches are listed")
        found = found[:MAX_REPEATS]
    return Result(found, "; ".join(notes))


# ---------------------------------------------------------------------------------------------------- the dot plot

class DotPlot(QWidget):
    """One sequence against another (or itself): x is a position in the first (`length` bases), y in the second (`rows` bases), a
    blue dot where the window at x equals the window at y, a red one where it equals its reverse complement (lines running the
    other way). Click a point to say its place; `picked` gives (x, y) in bases. Colours come from the theme's roles, and the picture is kept until the size, data or theme changes."""
    picked = Signal(int, int)
    MARGIN = (52, 14, 14, 36)   # left, top, right, bottom: room for the axis numbers

    def __init__(self, parent=None):
        super().__init__(parent)
        self.length, self.rows, self.direct, self.inverted, self.window = 0, 0, [], [], 12
        self.cache, self.hover = None, None
        self.setMinimumSize(260, 260)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.setMouseTracking(True)

    def set_data(self, length, rows, direct, inverted, window):
        self.length, self.rows, self.direct, self.inverted, self.window = length, rows, direct, inverted, window
        self.theme_changed()

    def theme_changed(self):
        self.cache = None
        self.update()

    def plot_rect(self):
        l, t, r, b = self.MARGIN
        side = max(10, min(self.width() - l - r, self.height() - t - b))
        return QRectF(l, t, side, side)

    def to_bases(self, pos):
        rect = self.plot_rect()
        x = (pos.x() - rect.left()) / rect.width() * self.length
        y = (pos.y() - rect.top()) / rect.height() * self.rows
        return (int(x), int(y)) if 0 <= x < self.length and 0 <= y < self.rows else None

    def _picture(self, rect):
        """The dots on their own transparent image, drawn once per size / data / theme."""
        key = (self.size(), theme.is_dark())
        if self.cache and self.cache[0] == key:
            return self.cache[1]
        image = QImage(self.size(), QImage.Format_ARGB32_Premultiplied)
        image.fill(Qt.transparent)
        p = QPainter(image)
        sx, sy = rect.width() / max(self.length, 1), rect.height() / max(self.rows, 1)
        for points, role in ((self.direct, "dot_direct"), (self.inverted, "dot_inverted")):
            p.setPen(QPen(theme.color(role), 2, Qt.SolidLine, Qt.SquareCap))
            for x, y in points:
                p.drawPoint(QPointF(rect.left() + x * sx, rect.top() + y * sy))
        p.end()
        self.cache = (key, image)
        return image

    def paintEvent(self, event):
        p = QPainter(self)
        p.fillRect(self.rect(), self.palette().window())
        rect = self.plot_rect()
        p.fillRect(rect, self.palette().base())
        text, faint = self.palette().windowText().color(), QColor(theme.hex("status_info"))
        if self.length:
            p.drawImage(0, 0, self._picture(rect))
            fm = QFontMetrics(self.font())
            label = lambda v, step: f"{v / 1000:g} kb" if step >= 1000 else f"{v:,}"
            step = nice_step(self.length / 6)
            for v in range(0, self.length + 1, step):   # across: the first sequence
                x = rect.left() + v / self.length * rect.width()
                p.setPen(QPen(faint, 1))
                p.drawLine(QPointF(x, rect.bottom()), QPointF(x, rect.bottom() + 4))
                p.setPen(text)
                p.drawText(QPointF(x - fm.horizontalAdvance(label(v, step)) / 2, rect.bottom() + 6 + fm.ascent()), label(v, step))
            step = nice_step(self.rows / 6)
            for v in range(0, self.rows + 1, step):   # down: the second
                y = rect.top() + v / self.rows * rect.height()
                p.setPen(QPen(faint, 1))
                p.drawLine(QPointF(rect.left() - 4, y), QPointF(rect.left(), y))
                p.setPen(text)
                p.drawText(QPointF(rect.left() - 8 - fm.horizontalAdvance(label(v, step)), y + fm.ascent() / 2 - 1), label(v, step))
        p.setPen(QPen(theme.color("group_border"), 1))
        p.setBrush(Qt.NoBrush)
        p.drawRect(rect)
        if self.hover:
            p.setPen(QPen(theme.color("highlight"), 1, Qt.DotLine))
            f = self.hover
            x, y = rect.left() + f[0] / self.length * rect.width(), rect.top() + f[1] / self.rows * rect.height()
            p.drawLine(QPointF(x, rect.top()), QPointF(x, rect.bottom()))
            p.drawLine(QPointF(rect.left(), y), QPointF(rect.right(), y))

    def mouseMoveEvent(self, event):
        self.hover = self.to_bases(event.position()) if self.length else None
        self.setToolTip("" if not self.hover else f"A {self.hover[0] + 1:,} (across) against B {self.hover[1] + 1:,} (down)")
        self.update()

    def leaveEvent(self, event):
        self.hover = None
        self.update()

    def mousePressEvent(self, event):
        at = self.to_bases(event.position()) if self.length and event.button() == Qt.LeftButton else None
        if at:
            self.picked.emit(*at)


# ---------------------------------------------------------------------------------------------------- the paired view

class PairRepeats(QWidget):
    """The "Dot plot" tab of Local Align's paired mode: the stretches the two sequences share (list) and their dot plot."""

    def __init__(self):
        super().__init__()
        self.a = self.b = ""
        self.result = Result([], "")
        spin = lambda lo, hi, v, tip: QSpinBox(minimum=lo, maximum=hi, value=v, toolTip=tip)
        self.shortest = spin(8, 500, 12, "Shortest shared stretch listed (20 bases and more can recombine in the cell)")
        self.mismatch = spin(0, 3, 0, "Differences allowed inside a shared stretch")
        self.window = spin(6, 32, 12, "Window of the dot plot: a dot is a window that matches in the other sequence")
        top = QHBoxLayout()
        for text, w in (("Shared ≥", self.shortest), ("Mismatches", self.mismatch), ("Plot window", self.window)):
            top.addWidget(QLabel(text))
            top.addWidget(w)
        top.addWidget(QPushButton("Find", clicked=self.reload))
        top.addStretch()
        self.tree = QTreeWidget()
        self.tree.setHeaderLabels(["Kind", "Length", "In A", "In B", "Differences"])
        self.tree.setRootIsDecorated(False)
        self.tree.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.tree.setEditTriggers(QAbstractItemView.NoEditTriggers)
        fit_columns(self.tree, stretch=2)
        self.plot = DotPlot()
        self.plot.picked.connect(lambda x, y: self.status.setText(f"A {x + 1:,} (across) against B {y + 1:,} (down)."))
        split = QSplitter(Qt.Horizontal)
        split.addWidget(self.tree)
        split.addWidget(self.plot)
        split.setSizes([420, 420])
        self.status = QLabel(wordWrap=True, textFormat=Qt.PlainText)
        self.status.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Preferred)
        layout = QVBoxLayout(self)
        layout.addLayout(top)
        layout.addWidget(split, 1)
        layout.addWidget(self.status)

    def show_pair(self, a, b):
        self.a, self.b = a, b
        self.reload()

    def reload(self):
        a, b, window = self.a, self.b, self.window.value()
        args = (self.shortest.value(), self.mismatch.value())
        # a long pair takes seconds: in a thread behind the busy dialog, so the window stays alive
        self.result, (direct, inverted, note) = run_busy(
            self, "Comparing…", lambda: (pair_repeats(a, b, *args), pair_dots(a, b, window)))
        self.tree.clear()
        for kind, i, j, length, mm in self.result.repeats:
            self.tree.addTopLevelItem(QTreeWidgetItem([KIND_TEXT[kind], f"{length:,} bp", f"{i + 1:,}–{i + length:,}",
                                                       f"{j + 1:,}–{j + length:,}" + ("←" if kind == "inverted" else ""),
                                                       str(mm) if mm else "exact"]))
        self.plot.set_data(min(len(a), MAX_LENGTH), min(len(b), MAX_LENGTH), direct, inverted, window)
        notes = [x for x in (self.result.capped, note) if x]
        self.status.setText(plural(len(self.result.repeats), "shared stretch", "shared stretches") + f" between {len(a):,} bp (A) and {len(b):,} bp (B)"
                            + (". Not complete: " + "; ".join(notes) if notes else "") + ".")


KIND_TEXT = {"direct": "Direct", "inverted": "Inverted"}
