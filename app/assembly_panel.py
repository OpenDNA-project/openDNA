"""The Assembly panel: a tab of its own (Tools > Assembly; not a pane, so it is not in the pane drop-down).
Top: the method and its options, then the map: small fragment boxes with a compact pill for every junction between
them (what joins, Tm, GC, and the verdict of its checks), an arrow that loops back to Fragment 1 when the product
is a circle, and, beside it, the product as a ring (circle) or a bar (line). Bottom: the program's own sequence view,
read-only, on the file behind the selected fragment: select the stretch to use, or click the two primers that
amplify it. The primers designed for the fragment are drawn there too, with their tails. A click on a junction opens the
junction inspector (JunctionDialog): the same sequence view on the product around the joint, the top strand alone, with
what each fragment and each primer contributes marked on it and the checks of that junction under it.

Flow: choose a method and the number of fragments; click a gray box and say where its DNA comes from;
select a stretch (or two primers) in the sequence view and press Use; the junctions are checked at
once; "Design primers" designs primers that carry what the method needs; when every junction
checks out, "Perform assembly" makes the exact product. The references are in Help > Citations.
"""
import math
from html import escape
from pathlib import Path
from types import SimpleNamespace

from PySide6.QtCore import QPointF, QRectF, Qt, Signal
from PySide6.QtGui import (
    QColor, QFont, QFontDatabase, QFontMetrics, QKeySequence, QPainter, QPainterPath, QPalette, QPen, QPolygonF,
)
from PySide6.QtWidgets import (
    QAbstractItemView, QCheckBox, QComboBox, QDialog, QDoubleSpinBox, QFileDialog, QFrame, QGridLayout,
    QHBoxLayout, QLabel, QMenu, QMessageBox, QPushButton, QRadioButton, QScrollArea, QSpinBox, QSplitter,
    QTableWidget, QStackedWidget, QTableWidgetItem, QTabWidget, QTextBrowser, QVBoxLayout, QWidget,
)

import runner
from assembly import (
    FRAG_COLORS, GG_ENZYMES, JOINED_BY, JOINT_COLOR, TAIL_COLOR, Frag, analyse, auto_primers, joint_tm, joints_html, junction_cards,
    junction_html, junction_summary, junction_window, oligo_checks, piece_doc, simulate, snapshot,
)
from editor import Editor, feature_color, primer_sequence, reverse_complement as rc
from primer_db import db
from primer_panel import MARKS, checks_html, melting_temp, worst
from sequence_view import SequenceView
from sourcefiles import OPEN_FILTER, add_recent, read_document, recent_files
from uikit import fit_columns, plural, right_aligned
import theme
import wip

MAX_FRAGMENTS = 12
WHOLE_FILE_MAX = 200_000    # bp: a file up to this size is used whole as soon as it is chosen
LEVEL_COLORS = theme.Live({"ok": "status_ok", "warn": "status_warn", "bad": "status_bad", "info": "status_info"})


def mono(points):
    font = QFontDatabase.systemFont(QFontDatabase.FixedFont)
    font.setPointSizeF(points)
    return font


# ---- the map ----------------------------------------------------------------------------------------------


class AssemblyMap(QWidget):
    """The assembly in order: small fragment boxes and, between them, a junction. An empty fragment is a gray
    box named Fragment N; a filled one is an arrow (the way it is read: left = reverse-complemented) with a mini
    map of its features. A junction is a compact pill: what joins the two neighbours and its Tm and GC, in the
    colour of its verdict (green ✓ passed; orange ? and red ✗ name the checks that did not pass, the reasons are in
    its tooltip and its detail). A circular assembly has one junction more, the last, and an arrow that leaves it,
    runs back over the boxes and enters Fragment 1: the loop is closed where the eye follows the fragments, not only
    in the product. Click a gray box (or the ▾ of a filled one) to choose where its DNA comes from, click a filled
    one to select it, drag to reorder, double-click to reverse-complement, click a junction for its full detail.
    The product itself is the ProductTile beside the map."""

    picked = Signal(int)                 # a fragment was selected
    source_asked = Signal(int, object)   # (index, global position): choose where this fragment's DNA comes from
    menu_asked = Signal(int, object)     # right-click: everything that can be done to this fragment
    reordered = Signal(int, int)         # drag: move fragment i to position j
    flip_asked = Signal(int)             # double-click
    junction_clicked = Signal(int)

    BOX_W, BOX_H, TOP = 124, 78, 8
    LANE, LANES = 8, 2                   # height of a feature lane in a box, lanes per box
    GAP, PILL_MIN = 46, 84
    LOOP_H, LOOP_X, LOOP_R = 30, 34, 7   # the return arrow: the band above the boxes, the room left of Fragment 1, its corner radius
    EMPTY_HINT = "click to choose a source"
    BOTTOM = 8                           # room under the boxes
    # How to use it, as the map's tooltip where nothing else is under the pointer (it was a permanent line under the map).
    HELP = ("Click a gray box to choose its source · drag a fragment to reorder · right-click a fragment for everything else · "
            "click a junction for details")

    def __init__(self):
        super().__init__()
        self.slots, self.circular, self.cards, self.hints, self.fchecks = [], False, [], [], {}
        self.current, self.dragging, self.tracks = -1, False, []
        self.boxes, self.card_rects = [], []   # card_rects: where each junction (its pill) is
        self.loop = None                       # the return arrow of a circular assembly: (leave x, y, right x, top y, left x)
        self.mono_font, self.small = mono(8.5), QFont(self.font())
        self.small.setPointSizeF(max(7.0, self.small.pointSizeF() - 1.5))
        self.setMouseTracking(True)
        self._compute()

    def refresh(self, slots, circular, cards=(), select=None, hints=(), fchecks=None):
        """`cards`: junction_cards(design) after an analysis, else none. `hints`: a line for each empty box.
        `fchecks`: {fragment: [(name, level, text)]} drawn as a mark on the box."""
        self.slots, self.circular, self.cards = list(slots), circular, list(cards)
        self.hints, self.fchecks = list(hints), dict(fchecks or {})
        if select is not None:
            self.current = select
        self.current = min(self.current, len(self.slots) - 1)
        self.tracks = [self._track(f) if f else [] for f in self.slots]
        self._compute()
        self.update()

    # -- the features drawn inside a box
    def _track(self, f):
        """[(parts, lane, colour, strand, label)], parts as fractions of the fragment's length; features that would
        overlap go in the next lane (2 lanes, 300 features)."""
        length, lanes, out = max(len(f.seq), 1), [], []
        for ft in sorted(f.doc.features(), key=lambda x: x.spans[0][0])[:300]:
            parts = [(s / length, e / length) for s, e in ft.spans]
            first, last = min(a for a, _ in parts), max(b for _, b in parts)
            lane = next((k for k, end in enumerate(lanes) if end <= first), len(lanes))
            if lane >= self.LANES:
                continue
            if lane == len(lanes):
                lanes.append(last)
            else:
                lanes[lane] = last
            out.append((parts, lane, feature_color(ft), ft.strand, ft.label))
        return out

    def _inner(self, box, flipped):
        head = min(12.0, box.width() / 4)
        return box.left() + 5 + (head if flipped else 0), box.width() - 10 - head, box.top() + 20

    def _feature_rect(self, inner, part, lane):
        left, width, top = inner
        a, b = part
        x0 = left + a * width
        return QRectF(x0, top + lane * self.LANE, max(2.0, left + b * width - x0), self.LANE - 1.5)

    def _draw_features(self, p, box, i):
        inner = self._inner(box, self.slots[i].flipped)
        left, width, top = inner
        p.setPen(Qt.NoPen)
        p.setBrush(QColor(255, 255, 255, 90))
        p.drawRoundedRect(QRectF(left - 2, top - 2, width + 4, self.LANE * self.LANES + 2), 3, 3)
        p.setPen(QPen(QColor(0, 0, 0, 110), 0.7))
        for parts, lane, color, strand, _ in self.tracks[i]:
            p.setBrush(color)
            for part in parts:
                r = self._feature_rect(inner, part, lane)
                tip, mid = min(4.0, r.width() / 2), r.center().y()
                pts = ([(r.right(), r.top()), (r.left() + tip, r.top()), (r.left(), mid), (r.left() + tip, r.bottom()), (r.right(), r.bottom())]
                       if strand < 0 else
                       [(r.left(), r.top()), (r.right() - tip, r.top()), (r.right(), mid), (r.right() - tip, r.bottom()), (r.left(), r.bottom())])
                p.drawPolygon(QPolygonF([QPointF(x, y) for x, y in pts]))

    def _feature_at(self, i, pos):
        inner = self._inner(self.boxes[i], self.slots[i].flipped)
        for parts, lane, _, _, label in self.tracks[i]:
            if any(self._feature_rect(inner, part, lane).contains(pos) for part in parts):
                return label

    # -- geometry: boxes with a junction between each pair, left to right
    @staticmethod
    def _failed(card):
        """The short names of the checks of a junction that did not pass (an orange ? or a red ✗), worst first."""
        return [name for name, level, _ in sorted(card.chips, key=lambda c: c[1] != "bad") if level in ("warn", "bad")]

    def _pill_lines(self, card):
        """(bold line, second line, level of the second): what joins, then Tm and GC when it is fine, else what is wrong."""
        main = ("↻ " if card.closing else "") + MARKS.get(card.level, ("", ""))[0] + " " + card.brief[0]
        if card.level not in ("warn", "bad"):
            return main, card.brief[1], "ok"
        failed = self._failed(card)
        return main, (" · ".join(failed[:3]) + (" …" if len(failed) > 3 else "")) if card.brief[1] else "click for why", card.level

    def _pill_size(self, card):
        bold = QFont(self.small)
        bold.setBold(True)
        main, second, _ = self._pill_lines(card)
        sm = QFontMetrics(self.small)
        return max(self.PILL_MIN, max(QFontMetrics(bold).horizontalAdvance(main), sm.horizontalAdvance(second)) + 22), 2 * sm.height() + 10

    def _compute(self):
        n = len(self.slots)
        sizes = [self._pill_size(c) for c in self.cards]
        # The band above the boxes is kept in a line too: ticking Circular moves nothing below the map.
        mid = self.TOP + self.LOOP_H + (self.BOX_H + 16) / 2
        self.boxes, self.card_rects, self.loop = [], [], None
        x = float(self.LOOP_X if self.circular else 10)
        for i in range(n):
            self.boxes.append(QRectF(x, mid - self.BOX_H / 2, self.BOX_W, self.BOX_H))
            x += self.BOX_W
            if i < len(sizes):
                w, h = sizes[i]
                self.card_rects.append(QRectF(x + 8, mid - h / 2, w, h))
                x += 16 + w
            elif i < n - 1:
                x += self.GAP
        if self.circular and self.boxes:
            # It leaves the last junction (the last box while there are no junctions yet), runs back above the boxes, enters Fragment 1.
            leave = self.card_rects[-1].right() if len(self.card_rects) == n else self.boxes[-1].right()
            x = leave + 2 * self.LOOP_R + 6
            self.loop = (leave, mid, x - self.LOOP_R, self.TOP + 20, self.LOOP_R + 1)
        self.setMinimumSize(int(x + 10), self.typical_height())
        self.setFixedHeight(self.typical_height())

    def typical_height(self):
        """The map's height, always the same (the boxes, the band of the return arrow, a little room under): the panel keeps
        it whatever the junctions, so the sequence view below does not jump as they come and go."""
        return int(self.TOP + self.LOOP_H + self.BOX_H + 16 + self.BOTTOM)

    def _at(self, pos):
        if not self.boxes:
            return -1
        return min(range(len(self.boxes)), key=lambda i: abs(self.boxes[i].center().x() - pos.x()))

    @staticmethod
    def _arrow(box, flipped):
        head = min(12.0, box.width() / 4)
        t, b, m = box.top(), box.bottom(), box.center().y()
        if flipped:
            pts = [(box.right(), t), (box.left() + head, t), (box.left(), m), (box.left() + head, b), (box.right(), b)]
        else:
            pts = [(box.left(), t), (box.right() - head, t), (box.right(), m), (box.right() - head, b), (box.left(), b)]
        return QPolygonF([QPointF(x, y) for x, y in pts])

    @staticmethod
    def _source_rect(box):
        return QRectF(box.right() - 24, box.top() + 1, 18, 15)

    # -- painting
    def paintEvent(self, event):
        p = QPainter(self)
        p.setRenderHint(QPainter.Antialiasing)
        dim = self.palette().color(QPalette.PlaceholderText)
        for i, (box, f) in enumerate(zip(self.boxes, self.slots)):
            self._draw_box(p, i, box, f)
            if i < len(self.card_rects):
                self._draw_pill(p, i)
            elif i < len(self.boxes) - 1:   # nothing to say yet: an arrow
                y = box.center().y()
                p.setPen(QPen(dim, 1.4))
                p.drawLine(QPointF(box.right() + 6, y), QPointF(box.right() + self.GAP - 6, y))
                p.drawLine(QPointF(box.right() + self.GAP - 12, y - 4), QPointF(box.right() + self.GAP - 6, y))
                p.drawLine(QPointF(box.right() + self.GAP - 12, y + 4), QPointF(box.right() + self.GAP - 6, y))
        if self.loop:
            self._draw_loop(p, dim)

    def _draw_loop(self, p, dim):
        """The arrow that closes a circle, in the colour of the closing junction's verdict (gray until there is one)."""
        leave, y, right, top, left = self.loop
        n, r, first = len(self.slots), self.LOOP_R, self.boxes[0].left()
        closing = self.cards[-1] if len(self.cards) == n else None
        color = QColor(LEVEL_COLORS.get(closing.level, LEVEL_COLORS["info"])) if closing else dim
        path = QPainterPath(QPointF(leave, y))
        path.lineTo(right - r, y)
        path.quadTo(right, y, right, y - r)
        path.lineTo(right, top + r)
        path.quadTo(right, top, right - r, top)
        path.lineTo(left + r, top)
        path.quadTo(left, top, left, top + r)
        path.lineTo(left, y - r)
        path.quadTo(left, y, left + r, y)
        path.lineTo(first - 8, y)
        p.setBrush(Qt.NoBrush)
        p.setPen(QPen(color, 1.8))
        p.drawPath(path)
        p.setBrush(color)
        p.setPen(Qt.NoPen)
        p.drawPolygon(QPolygonF([QPointF(first - 1, y), QPointF(first - 10, y - 5), QPointF(first - 10, y + 5)]))
        p.setFont(self.small)
        p.setPen(color)
        p.drawText(QRectF(left + 12, self.TOP, max(right - left - 24, 0), 14), Qt.AlignLeft | Qt.AlignVCenter,
                   f"junction {n} → 1 closes the circle" if n > 1 else "the end closes onto the start")

    def _draw_pill(self, p, j):
        card, rect = self.cards[j], self.card_rects[j]
        color = QColor(LEVEL_COLORS.get(card.level, LEVEL_COLORS["info"]))
        base = self.palette().color(QPalette.Base)   # the theme's background with a tenth of the verdict's colour in it
        p.setBrush(QColor(int(base.red() * 0.9 + color.red() * 0.1), int(base.green() * 0.9 + color.green() * 0.1),
                          int(base.blue() * 0.9 + color.blue() * 0.1)))
        p.setPen(QPen(color, 2.4 if card.level == "bad" else 1.5))
        p.drawRoundedRect(rect, 12, 12)
        main, second, level = self._pill_lines(card)
        bold = QFont(self.small)
        bold.setBold(True)
        line = QFontMetrics(self.small).height()
        inner = rect.adjusted(6, 4, -6, -4)
        p.setFont(bold)
        p.setPen(color)
        p.drawText(QRectF(inner.left(), inner.top(), inner.width(), line), Qt.AlignCenter, main)
        p.setFont(self.small)
        p.setPen(self.palette().color(QPalette.Text) if level == "ok" else QColor(LEVEL_COLORS[level]))
        p.drawText(QRectF(inner.left(), inner.top() + line + 2, inner.width(), line), Qt.AlignCenter, second)

    # The gray boxes (an empty fragment) take the theme's grays: fixed light grays glared in the dark.
    def _gray(self, current=False):
        return self.palette().color(QPalette.Button if current else QPalette.Midlight)

    def _faint(self):
        c = QColor(self.palette().color(QPalette.Text))
        c.setAlpha(170)
        return c

    def _draw_box(self, p, i, box, f):
        text = self.palette().color(QPalette.Text)
        if f is None:
            p.setBrush(self._gray(i == self.current))
            p.setPen(QPen(text if i == self.current else self.palette().color(QPalette.PlaceholderText), 2.5 if i == self.current else 1.2,
                          Qt.DashLine))
            p.drawRoundedRect(box, 6, 6)
            p.setPen(text)
            p.setFont(self.font())
            p.drawText(box.adjusted(4, 8, -4, -34), Qt.AlignCenter, f"Fragment {i + 1}")
            p.setFont(self.small)
            p.setPen(self._faint())
            hint = self.hints[i] if i < len(self.hints) and self.hints[i] else self.EMPTY_HINT
            p.drawText(box.adjusted(4, 40, -4, -4), Qt.AlignCenter | Qt.TextWordWrap, hint)
            return
        color = QColor(FRAG_COLORS[i % len(FRAG_COLORS)])
        p.setBrush(color.lighter(135) if i == self.current else color)
        p.setPen(QPen(text if i == self.current else color.darker(150), 2.5 if i == self.current else 1.2,
                      Qt.DashLine if f.vector else Qt.SolidLine))
        p.drawPolygon(self._arrow(box, f.flipped))
        p.setPen(QColor("#1b1b1b"))
        p.setFont(self.small)
        flags = self.fchecks.get(i, [])
        level = "bad" if any(c[1] == "bad" for c in flags) else "warn" if any(c[1] == "warn" for c in flags) else None
        p.drawText(box.adjusted(6, 2, -38 if level else -26, -(box.height() - 18)), Qt.AlignLeft | Qt.AlignVCenter,
                   p.fontMetrics().elidedText(f"{i + 1} · {f.name}", Qt.ElideMiddle, int(box.width()) - (44 if level else 32)))
        p.drawText(self._source_rect(box), Qt.AlignCenter, "▾")
        if level:
            p.setPen(QColor(LEVEL_COLORS[level]))
            p.drawText(QRectF(box.right() - 40, box.top() + 1, 16, 15), Qt.AlignCenter, MARKS[level][0])
            p.setPen(QColor("#1b1b1b"))
        self._draw_features(p, box, i)
        tails = f.tails[0] + f.tails[1]
        p.drawText(box.adjusted(6, box.height() - 18, -6, -2), Qt.AlignCenter,
                   f"{len(f.seq):,} bp" + (f" +{tails}" if tails else "") + ("  ◂" if f.flipped else "") + ("  [vector]" if f.vector else ""))

    # -- mouse
    def mousePressEvent(self, event):
        for j, r in enumerate(self.card_rects):
            if r.contains(event.position()):
                self.junction_clicked.emit(j)
                return
        i = self._at(event.position())
        if i < 0:
            return
        f = self.slots[i]
        self.current, self.dragging = i, False
        self.picked.emit(i)
        if event.button() == Qt.RightButton:
            self.menu_asked.emit(i, event.globalPosition().toPoint())
            return
        if f is None or self._source_rect(self.boxes[i]).contains(event.position()):
            self.source_asked.emit(i, event.globalPosition().toPoint())
            return
        self.dragging = event.button() == Qt.LeftButton
        self.update()

    def mouseMoveEvent(self, event):
        over = any(b.contains(event.position()) for b in self.boxes) or any(r.contains(event.position()) for r in self.card_rects)
        self.setCursor(Qt.PointingHandCursor if over else Qt.ArrowCursor)   # boxes and junctions can be clicked
        for r, card in zip(self.card_rects, self.cards):
            if r.contains(event.position()):
                # The pill shows two short lines, so its tooltip is the rest in words: what joins, then every check.
                self.setToolTip(f"{card.title}\n\n{card.text}".strip())
                return
        i = self._at(event.position())
        if self.dragging and i >= 0 and i != self.current and self.slots[i] is not None:
            self.reordered.emit(self.current, i)
        elif not self.dragging and i >= 0 and self.boxes[i].contains(event.position()):
            f = self.slots[i]
            if f is None:
                self.setToolTip(f"Fragment {i + 1}: click to choose where its DNA comes from")
                return
            feature = self._feature_at(i, event.position())
            flags = "\n".join(f"{MARKS[lv][0]} {t}" for _, lv, t in self.fchecks.get(i, []))
            self.setToolTip((f"{feature}\n\n" if feature else "") + f"{f.name}\n{len(f.seq):,} bp, {plural(len(f.doc.features()), 'feature')}"
                            + ("\nreverse-complemented" if f.flipped else "") + ("\nGolden Gate destination vector" if f.vector else "")
                            + (f"\n\n{flags}" if flags else ""))
        else:
            self.setToolTip(self.HELP)

    def mouseReleaseEvent(self, event):
        self.dragging = False

    def mouseDoubleClickEvent(self, event):
        i = self._at(event.position())
        if i >= 0 and self.slots[i] is not None and self.boxes[i].contains(event.position()):
            self.flip_asked.emit(i)


class ProductTile(QFrame):
    """The product, at the right end of the map, which scrolls (when there are too many fragments to fit) while this does not: a ring when the assembly is a circle, a bar when it is
    a line, built from the fragments' own colours and lengths, with a dot on each junction in the colour of its verdict (the one
    at the top of a ring closes the circle). A dashed gray outline until every junction fits; then a click shows the exact
    sequence below."""

    clicked = Signal()
    WIDTH, RING, BAND = 132, 60, 9
    MIN_SHARE = 0.05   # a very short fragment still gets this much of the ring (or bar), so that it can be seen

    def __init__(self):
        super().__init__()
        self.setFrameShape(QFrame.StyledPanel)   # the same frame as the map beside it
        self.parts, self.levels, self.circular, self.text, self.selected = [], [], False, None, False
        self.inset = 0   # room kept free at the bottom (the map keeps some for its scroll bar), so the ring lines up with the boxes
        self.small = QFont(self.font())
        self.small.setPointSizeF(max(7.0, self.small.pointSizeF() - 1.5))
        self.setFixedWidth(self.WIDTH)
        self.refresh(False)

    def refresh(self, circular, parts=(), levels=(), text=None):
        """`parts`: [(bases, colour)] of the fragments in product order; `levels`: the verdict of each junction (a circle's last one
        closes it); `text`: the size line ("4,323 bp circular") when the product can be previewed, else None."""
        self.circular, self.parts, self.levels, self.text = circular, list(parts), list(levels), text
        self.setCursor(Qt.PointingHandCursor if self.ready() else Qt.ArrowCursor)
        self.setToolTip("The exact sequence this assembly makes: click to see it below" if self.ready()
                        else "The product can be shown once every junction fits")
        self.update()

    def ready(self):
        return self.text is not None and bool(self.parts)

    def _ends(self):
        """Where each fragment ends, in turns of the whole (0..1), each given at least MIN_SHARE."""
        sizes = [max(b, 1) for b, _ in self.parts]
        shares = [max(s / sum(sizes), self.MIN_SHARE) for s in sizes]
        out, end = [], 0.0
        for s in shares:
            end += s / sum(shares)
            out.append(end)
        return out

    def mousePressEvent(self, event):
        if self.ready() and event.button() == Qt.LeftButton:
            self.clicked.emit()

    def paintEvent(self, event):
        super().paintEvent(event)   # the frame
        p = QPainter(self)
        p.setRenderHint(QPainter.Antialiasing)
        pal = self.palette()
        dim, text, base = pal.color(QPalette.PlaceholderText), pal.color(QPalette.Text), pal.color(QPalette.Base)
        ready, size, band, cx = self.ready(), self.RING, self.BAND, self.width() / 2
        line = QFontMetrics(self.small).height()
        group = size + 8 + 2 * line
        top = max(6.0, (self.height() - self.inset - group) / 2)
        ring = QRectF(cx - size / 2 + band / 2, top + band / 2, size - band, size - band)
        mid = QPointF(cx, top + size / 2)
        if self.selected:
            p.setPen(QPen(text, 2))
            p.setBrush(Qt.NoBrush)
            p.drawRoundedRect(QRectF(3, top - 5, self.width() - 6, group + 10), 8, 8)
        ends = self._ends() if ready else []
        if self.circular:
            if not ready:
                p.setPen(QPen(dim, 1.6, Qt.DashLine))
                p.setBrush(Qt.NoBrush)
                p.drawEllipse(ring)
            for k, end in enumerate(ends):
                start = ends[k - 1] if k else 0.0
                p.setPen(QPen(QColor(self.parts[k][1]), band, Qt.SolidLine, Qt.FlatCap))
                # Qt measures angles counter-clockwise from 3 o'clock, in sixteenths of a degree; the ring starts at 12 and runs clockwise.
                p.drawArc(ring, int((90 - (start * 360 + 1)) * 16), int(-((end - start) * 360 - 2) * 16))
            for k, end in enumerate(ends):
                if k < len(self.levels):
                    turn = math.radians(end * 360)
                    self._dot(p, QPointF(mid.x() + ring.width() / 2 * math.sin(turn), mid.y() - ring.height() / 2 * math.cos(turn)),
                              self.levels[k], base)
        else:
            width = self.width() - 36
            bar = QRectF(18, mid.y() - band / 2, width, band)
            if not ready:
                p.setPen(QPen(dim, 1.6, Qt.DashLine))
                p.setBrush(Qt.NoBrush)
                p.drawRoundedRect(bar, 4, 4)
            for k, end in enumerate(ends):
                start = ends[k - 1] if k else 0.0
                p.setPen(Qt.NoPen)
                p.setBrush(QColor(self.parts[k][1]))
                p.drawRect(QRectF(bar.left() + start * width + (1 if k else 0), bar.top(), (end - start) * width - (1 if k else 0), band))
            for k, end in enumerate(ends[:-1]):
                if k < len(self.levels):
                    self._dot(p, QPointF(bar.left() + end * width, bar.center().y()), self.levels[k], base)
        y = top + size + 8
        if ready:
            size_line, _, shape = self.text.rpartition(" ")
            p.setPen(text)
            p.setFont(self.font())
            p.drawText(QRectF(0, y, self.width(), line + 2), Qt.AlignHCenter | Qt.AlignTop, size_line)
            p.setFont(self.small)
            p.drawText(QRectF(0, y + line + 2, self.width(), line), Qt.AlignHCenter | Qt.AlignTop, shape)
        else:
            p.setPen(dim)
            p.setFont(self.small)
            p.drawText(QRectF(6, y - 2, self.width() - 12, 2 * line + 2), Qt.AlignHCenter | Qt.AlignTop | Qt.TextWordWrap,
                       "Product, when every junction fits")

    @staticmethod
    def _dot(p, at, level, base):
        p.setPen(QPen(base, 1.6))
        p.setBrush(QColor(LEVEL_COLORS.get(level, LEVEL_COLORS["info"])))
        p.drawEllipse(at, 4.6, 4.6)


# ---- the sequence view ----------------------------------------------------------------------------------


class FragmentEditor(Editor):
    """An Editor on a copy of a source file, behind the panel's sequence view. Besides the primer database it looks
    for the file's own primers and the ones designed for the fragment (`set_extra`), so they are drawn like any primer, tail included."""


class ReadOnlyView(SequenceView):
    """The sequence view of the program with its editing taken away: select, scroll, zoom and copy only."""

    def keyPressEvent(self, event):
        key = event.key()
        typing = event.text() and not event.modifiers() & (Qt.ControlModifier | Qt.AltModifier)
        if event.matches(QKeySequence.Cut) or event.matches(QKeySequence.Paste) or key in (Qt.Key_Backspace, Qt.Key_Delete) or typing:
            return
        super().keyPressEvent(event)

    def mousePressEvent(self, event):
        if event.button() != Qt.RightButton:   # its menus edit features and the sequence
            super().mousePressEvent(event)


# ---- details, on request ------------------------------------------------------------------------------------


class ReportDialog(QDialog):
    """Everything about the current design in text: every junction in detail, every check, every oligo."""

    def __init__(self, parent, d):
        super().__init__(parent)
        self.setWindowTitle("Assembly Report")
        self.resize(900, 700)
        junctions, report = QTextBrowser(), QTextBrowser()
        junctions.setHtml(junction_html(d))
        how = "".join(f"<li>{escape(j)}</li>" for j in d.junctions)
        report.setHtml(f"<b>Checks</b>{checks_html(d.checks)}<br><b>Junctions</b>{joints_html(d)}<br><b>How each junction is made</b><ul>{how}</ul>")
        table = QTableWidget(len(d.rows), 5)
        table.setHorizontalHeaderLabels(["Name", "Sequence 5′→3′ (tail in lowercase)", "Anneal Tm", "Length", "Checks"])
        table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        fit_columns(table, stretch=1, cap=480)   # an oligo's whole sequence, up to 480 px
        for r, row in enumerate(d.rows):
            checks = [c for c in oligo_checks(row) if row[2] or c[0] in ("Hairpin", "Self-dimer")]
            t = melting_temp(row[1][len(row[1]) - row[2]:]) if row[2] else None
            for c, v in enumerate([row[0], row[1], f"{t:.1f} °C" if t else "", f"{len(row[1])} nt", MARKS[worst(checks)][0]]):
                item = QTableWidgetItem(v)
                if c in (2, 3):
                    right_aligned(item)
                item.setToolTip(checks_html(checks))
                table.setItem(r, c, item)
        tabs = QTabWidget()
        for widget, title in ((junctions, "Junctions"), (report, "Checks"), (table, "Oligos")):
            tabs.addTab(widget, title)
        QVBoxLayout(self).addWidget(tabs)


class JunctionDialog(QDialog):
    """The junction inspector: one junction in the program's own sequence view, on the top strand alone. It shows the product
    around the joint (the joint selected), so what is where is read the way it is everywhere else: each fragment a bar in its
    colour, the joint a slate bar, what a primer's tail added an orange bar, the primers drawn as primers, and every base, primer
    and feature answers to the pointer as in any other view. A sentence says how much of the joint a primer adds and how much the
    two fragments already shared; the checks of the junction are listed under the sequence. A junction whose ends do not fit has
    no product to show: the text of its two ends and the reason stand in the view's place."""

    def __init__(self, parent, slots, d, j):
        super().__init__(parent)
        n, jt = len(d.names), d.joints[j]
        a, b = j % n, (j + 1) % n
        self.setWindowTitle(f"Junction {j + 1}: {jt.left} → {jt.right}")
        self.resize(1000, 600)
        level = d.jlevels[j] if j < len(d.jlevels) else "ok"
        colour = LEVEL_COLORS.get(level, LEVEL_COLORS["info"])
        if jt.kind == "none":
            what = "✗ the two ends do not fit"
        else:
            unit, word = ("nt", "toehold") if jt.kind == "toehold" else ("bp", JOINED_BY[jt.kind])
            tm, gc = joint_tm(jt.seq), 100 * sum(c in "GC" for c in jt.seq.upper()) / max(len(jt.seq), 1)
            what = f"{MARKS.get(level, ('', ''))[0]} {jt.kind}, {len(jt.seq)} {unit} {word}" + (f" · Tm {tm:.1f} °C" if tm is not None else "") + f" · GC {gc:.0f} %"
        head = QLabel(f"<b>Junction {j + 1}</b> &nbsp;Fragment {a + 1} → Fragment {b + 1} &nbsp;<span style='color:{colour};font-weight:600'>{escape(what)}</span>")
        orient = lambda i: f"Fragment {i + 1} = {escape(d.names[i])} {'◀ reverse-complemented' if d.flipped[i] else '▶ as given'}"
        summary = QLabel(f"{escape(junction_summary(d, j))}<br><span style='color:gray'>{orient(a)} &nbsp;·&nbsp; {orient(b)}</span>", wordWrap=True)
        box = QVBoxLayout(self)
        self.view, labels = None, [head, summary]
        if jt.kind == "none":   # no product to show: the two ends as they are, under the reason above
            text = QTextBrowser()
            text.setHtml(f"<tt>Fragment {a + 1} ends&nbsp;&nbsp;…{escape(d.seqs[a][-30:])}<br>Fragment {b + 1} starts {escape(d.seqs[b][:30])}…</tt>")
            body = text
        else:
            doc = simulate(slots, d, preview=True)
            window, at, first = junction_window(doc, d, j)
            swatch = lambda colour: f"<span style='background-color:{colour}'>&nbsp;&nbsp;&nbsp;&nbsp;</span>"
            last = first + window.length
            where = f"{first + 1:,}–{last:,}" if last <= doc.length else f"{first + 1:,}–{doc.length:,} and 1–{last - doc.length:,}"
            legend = QLabel(f"{swatch(FRAG_COLORS[a % len(FRAG_COLORS)])} Fragment {a + 1} &nbsp;{swatch(FRAG_COLORS[b % len(FRAG_COLORS)])} Fragment {b + 1}"
                            f" &nbsp;&nbsp;{swatch(TAIL_COLOR)} added by a primer's tail &nbsp;{swatch(JOINT_COLOR)} the joint"
                            f" &nbsp;<span style='color:{theme.hex('primer')}'><b>━▸</b></span> a primer (hover it)"
                            f" &nbsp;·&nbsp; top strand, 5′ → 3′, bases {where} of the product (the numbers count from the first shown)")
            labels.append(legend)
            editor = Editor(window, f"Junction {j + 1}")
            editor.set_display(complement=False)   # the second strand would only repeat the first
            self.view = body = ReadOnlyView(editor)
            body.resize(960, 280)
            editor.select(at, at + len(jt.seq))
        for label in labels:
            label.setTextFormat(Qt.RichText)
            label.setWordWrap(True)
            label.setTextInteractionFlags(Qt.TextSelectableByMouse)
            box.addWidget(label)
        box.addWidget(body, 3)
        box.addWidget(QLabel("<b>Checks</b>"))
        checks = QTextBrowser()
        checks.setHtml(checks_html(d.jchecks.get(j, [])))
        box.addWidget(checks, 2)


# ---- the panel ---------------------------------------------------------------------------------------------


METHOD_WIP = [None, "golden_gate", "sidewinder"]   # per method of the list: its feature in wip.FEATURES, if unfinished


class AssemblyPanel(QWidget):
    """See the module docstring. `win`: the main window (for the open files and for opening the product)."""

    def __init__(self, win):
        super().__init__()
        self.win = win
        self.slots = [None] * 3          # Frag or None, in assembly order
        self.sources = [None] * 3        # per slot: SimpleNamespace(doc, label, path, editor, view, file_primers) or None
        self.design, self.error, self.dialogs = None, "", []
        self.product_view, self.product_shown = None, False   # the preview of the product, in the sequence stack
        self.method = QComboBox()
        self.method.addItems(["Gibson / HiFi", "Golden Gate", "Sidewinder"])
        # The unfinished methods carry the "WIP" tag: in the list, and beside the box while one is chosen (wip.FEATURES).
        self.method.setItemDelegate(wip.Delegate(lambda index: METHOD_WIP[index.row()], self.method))
        self.method_wip = wip.Badge(None)
        if any(wip.is_wip(k) for k in METHOD_WIP if k):   # its room kept when hidden: the options on its right stay put
            keep = self.method_wip.sizePolicy()
            keep.setRetainSizeWhenHidden(True)
            self.method_wip.setSizePolicy(keep)
        method_box = QWidget()
        method_row = QHBoxLayout(method_box)
        method_row.setContentsMargins(0, 0, 0, 0)
        method_row.setSpacing(wip.GAP)
        method_row.addWidget(self.method)
        method_row.addWidget(self.method_wip)
        self.method.currentIndexChanged.connect(lambda m: self.method_wip.set_key(METHOD_WIP[m]))
        self.count = QSpinBox(minimum=1, maximum=MAX_FRAGMENTS, value=3)
        self.enzyme = QComboBox()
        self.enzyme.addItems(list(GG_ENZYMES))
        self.circular = QCheckBox("Circular product (plasmid)", checked=True)
        self.min_ov = QSpinBox(minimum=15, maximum=60, value=20, suffix=" bp")
        self.max_ov = QSpinBox(minimum=15, maximum=80, value=40, suffix=" bp")
        self.ov_tm = QDoubleSpinBox(minimum=40, maximum=75, value=52, suffix=" °C")
        self.anneal_tm = QDoubleSpinBox(minimum=45, maximum=72, value=60, suffix=" °C")
        self.report_button = QPushButton("Full report…", clicked=self.show_report, toolTip="Every junction, check and oligo in detail")
        self.options = []   # (label, widget, method it belongs to) shown only for that method
        grid = QGridLayout()
        grid.setHorizontalSpacing(8)
        for col, (label, w, only) in enumerate((("Method", method_box, None), ("Fragments", self.count, None), ("Enzyme", self.enzyme, 1),
                                                ("Overlap", self.min_ov, 0), ("to", self.max_ov, 0), ("Overlap Tm ≥", self.ov_tm, 0),
                                                ("Annealing Tm", self.anneal_tm, None))):
            lab = QLabel(label)
            grid.addWidget(lab, 0, col)
            grid.addWidget(w, 1, col)
            for x in (lab, w):   # a hidden option keeps its room: "Annealing Tm" and "Circular" stay put when the method changes
                keep = x.sizePolicy()
                keep.setRetainSizeWhenHidden(True)
                x.setSizePolicy(keep)
            self.options.append((lab, w, only))
        grid.addWidget(self.circular, 1, 7)
        grid.setColumnStretch(8, 1)
        grid.addWidget(self.report_button, 1, 9)
        self.strip = AssemblyMap()
        self.strip.picked.connect(self._pick)
        self.strip.source_asked.connect(self._ask_source)
        self.strip.menu_asked.connect(self._context_menu)
        self.strip.reordered.connect(self._move_to)
        self.strip.flip_asked.connect(self._flip_at)
        self.strip.junction_clicked.connect(self._show_junction)
        self.map_area = QScrollArea()
        self.map_area.setWidget(self.strip)
        self.map_area.setWidgetResizable(True)   # the map fills the whole width; a scroll bar only when the fragments do not fit
        self.map_area.setToolTip(AssemblyMap.HELP)
        # One height whatever the junctions, so the splitter below never jumps; the map keeps room for its scroll bar.
        room = self.map_area.horizontalScrollBar().sizeHint().height() + 6
        self.map_area.setFixedHeight(self.strip.typical_height() + room)
        # The product stays in view while the fragments scroll: the end of the line is always on screen.
        self.product_tile = ProductTile()
        self.product_tile.setFixedHeight(self.strip.typical_height() + room)
        self.product_tile.inset = room
        self.product_tile.clicked.connect(self._show_product)
        # The actions are a column on the right of the map: the verdict first, the buttons under it, and the one that
        # matters, "Perform assembly", last and highlighted.
        self.buttons = {}
        for key, text, tip, fn in (
                ("homology", "Design primers", "Design new primers that carry what the method needs at each junction "
                 "(an overlap tail for Gibson, a site + spacer + overhang for Golden Gate)", self.add_homology),
                ("save", "Save fragments…", "Write each fragment as it is after PCR to its own GenBank file", self.save_fragments),
                ("primers", "Save primers to database", "Keep the primers of this design in the primer database", self.save_primers),
                ("perform", "Perform assembly", "Make the exact sequence the assembly would give (needs every junction to check out)",
                 self.perform)):
            self.buttons[key] = QPushButton(text, clicked=fn, toolTip=tip)
        self.buttons["perform"].setObjectName("primary")
        self.status = QLabel(wordWrap=True)
        self.notes = QLabel(wordWrap=True, textFormat=Qt.RichText)
        actions = QWidget()
        actions.setFixedWidth(250)
        column = QVBoxLayout(actions)
        column.setContentsMargins(4, 0, 0, 0)
        column.addWidget(self.status)
        column.addWidget(self.notes)
        column.addStretch()
        for key in ("homology", "save", "primers"):
            column.addWidget(self.buttons[key])
        column.addSpacing(10)
        column.addWidget(self.buttons["perform"])
        upper = QWidget()
        box = QVBoxLayout(upper)
        box.setContentsMargins(6, 4, 6, 0)
        box.addLayout(grid)
        row = QHBoxLayout()
        row.setSpacing(4)
        row.addWidget(self.map_area, 1, Qt.AlignTop)
        row.addWidget(self.product_tile, 0, Qt.AlignTop)
        row.addWidget(actions)
        box.addLayout(row, 1)
        # the sequence view: one per source, the selected fragment's on top; its controls are a column on its right
        self.seqstack = QStackedWidget()
        self.nothing = QLabel("Choose a source for the selected fragment.", alignment=Qt.AlignCenter)
        self.seqstack.addWidget(self.nothing)
        self.which = QLabel()
        self.primer_info = QTextBrowser(maximumHeight=90)
        self.by_region = QRadioButton("Selected stretch", checked=True)
        self.by_primers = QRadioButton("PCR with two primers")
        self.forward, self.reverse = QComboBox(minimumContentsLength=14), QComboBox(minimumContentsLength=14)
        self.use = QPushButton("Use as fragment", clicked=self.use_selection)
        self.use.setObjectName("primary")
        self.chosen = QLabel(wordWrap=True, styleSheet="color: gray")
        controls = QWidget()
        controls.setFixedWidth(250)
        side = QVBoxLayout(controls)
        side.setContentsMargins(4, 0, 0, 0)
        for w in (self.by_region, self.by_primers, QLabel("Forward primer"), self.forward, QLabel("Reverse primer"), self.reverse):
            side.addWidget(w)
        side.addSpacing(8)
        side.addWidget(self.use)
        side.addWidget(self.chosen)
        side.addStretch()
        lower = QWidget()
        lbox = QVBoxLayout(lower)
        lbox.setContentsMargins(6, 0, 6, 4)
        lbox.addWidget(self.which)
        middle = QHBoxLayout()
        left = QVBoxLayout()
        left.addWidget(self.seqstack, 1)
        left.addWidget(self.primer_info)
        middle.addLayout(left, 1)
        middle.addWidget(controls)
        lbox.addLayout(middle, 1)
        # The two buttons that matter are drawn in the highlight colour of the theme.
        self.setStyleSheet("QPushButton#primary { background: palette(highlight); color: palette(highlighted-text); "
                           "font-weight: 700; padding: 8px 14px; border: 1px solid palette(highlight); border-radius: 4px; }"
                           "QPushButton#primary:disabled { background: palette(button); color: palette(mid); "
                           "border: 1px solid palette(mid); }")
        split = QSplitter(Qt.Vertical)
        split.setChildrenCollapsible(False)
        split.addWidget(upper)
        split.addWidget(lower)
        split.setSizes([360, 600])
        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.addWidget(split)
        self.method.currentIndexChanged.connect(self._method_changed)
        self.count.valueChanged.connect(self._count_changed)
        self.enzyme.currentIndexChanged.connect(self.analyse)
        self.circular.toggled.connect(self.analyse)
        # The shortest overlap never exceeds the longest: moving one pushes the other along.
        self.min_ov.valueChanged.connect(lambda v: self.max_ov.setValue(max(v, self.max_ov.value())))
        self.max_ov.valueChanged.connect(lambda v: self.min_ov.setValue(min(v, self.min_ov.value())))
        for w in (self.min_ov, self.max_ov, self.ov_tm, self.anneal_tm):
            w.valueChanged.connect(self.analyse)
        self._method_changed()

    @property
    def view(self):
        """The sequence view on show (None when no source is chosen)."""
        w = self.seqstack.currentWidget()
        return w if isinstance(w, ReadOnlyView) else None

    # -- state helpers
    @property
    def current(self):
        return self.strip.current

    def is_circular(self):
        return self.circular.isChecked() and self.method.currentIndex() != 2

    def _method_changed(self):
        m = self.method.currentIndex()
        for lab, w, only in self.options:
            lab.setVisible(only is None or only == m)
            w.setVisible(only is None or only == m)
        self.circular.setEnabled(m != 2)
        self.buttons["homology"].setToolTip(
            "Design new primers that carry what the method needs at each junction (an overlap tail for Gibson, a site + "
            "spacer + overhang for Golden Gate)" if m != 2 else "Sidewinder oligos are designed for you")
        self.analyse()

    def _count_changed(self, n):
        chosen = [i + 1 for i, src in enumerate(self.sources[n:], n) if src is not None]
        if chosen and QMessageBox.question(   # No by default: the sources chosen would be lost
                self, "Fewer fragments", f"{'Fragment' if len(chosen) == 1 else 'Fragments'} {', '.join(map(str, chosen))} "
                f"{'has' if len(chosen) == 1 else 'have'} a source chosen. Drop {'it' if len(chosen) == 1 else 'them'}?",
                QMessageBox.Yes | QMessageBox.No, QMessageBox.No) != QMessageBox.Yes:
            self.count.blockSignals(True)
            self.count.setValue(len(self.slots))
            self.count.blockSignals(False)
            return
        for src in self.sources[n:]:
            self._drop_source(src)
        self.slots = (self.slots + [None] * n)[:n]
        self.sources = (self.sources + [None] * n)[:n]
        self.strip.current = min(max(self.strip.current, 0), n - 1)
        self.analyse()

    def _free_name(self, name, skip=None):
        names = {f.name for i, f in enumerate(self.slots) if f and i != skip}
        return name if name not in names else next(f"{name} ({k})" for k in range(2, 999) if f"{name} ({k})" not in names)

    # -- choosing a source
    def _source_entries(self, menu, i):
        """Where fragment i's DNA can come from: open files, recent files, or a file to browse to."""
        opened = menu.addMenu("Currently open files")
        docs = self.win.open_documents()
        for label, doc in docs:
            opened.addAction(label, lambda label=label, doc=doc: self._set_source(i, snapshot(doc), label, None))
        opened.setEnabled(bool(docs))
        recent = menu.addMenu("Recently opened files")
        paths = recent_files()
        for path in paths:
            recent.addAction(Path(path).name, lambda path=path: self._open_path(i, path)).setToolTip(path)
        recent.setEnabled(bool(paths))
        menu.addAction("Browse…", lambda: self._browse(i))

    def _ask_source(self, i, pos):
        """A click on a gray box (or the arrow of a filled one): choose where its DNA comes from."""
        menu = QMenu(self)
        self._source_entries(menu, i)
        menu.exec(pos)

    def _context_menu(self, i, pos):
        """Right-click: show the menu of fragment i."""
        self._build_context_menu(i).exec(pos)

    def _build_context_menu(self, i):
        """Everything that can be done to fragment i."""
        menu = QMenu(self)
        self._source_entries(menu.addMenu("Source"), i)
        f = self.slots[i]
        if f is not None:
            menu.addSeparator()
            menu.addAction("Move left", lambda: self._move_to(i, i - 1)).setEnabled(i > 0)
            menu.addAction("Move right", lambda: self._move_to(i, i + 1)).setEnabled(i < len(self.slots) - 1)
            menu.addAction("Reverse complement", lambda: self._flip_at(i))
            if self.method.currentIndex() == 1:
                vector = menu.addAction("Destination vector (cut at its two sites)")
                vector.setCheckable(True)
                vector.setChecked(f.vector)
                vector.triggered.connect(lambda _checked=False: self._toggle_vector_at(i))
            menu.addSeparator()
            menu.addAction("Clear this fragment", lambda: self._clear_at(i))
        return menu

    def _browse(self, i):
        path, _ = QFileDialog.getOpenFileName(self, f"Source of fragment {i + 1}", "", OPEN_FILTER)
        if path:
            self._open_path(i, path)

    def _open_path(self, i, path):
        doc = read_document(self, path)
        if doc is not None:
            add_recent(path)
            self._set_source(i, snapshot(doc), Path(path).name, path)

    def _drop_source(self, src):
        """Forget a source: its view goes (its editor, linked weakly to the primer database, goes with it)."""
        if src is not None:
            self.seqstack.removeWidget(src.view)
            src.view.deleteLater()

    def _set_source(self, i, doc, label, path):
        """The DNA of fragment i comes from `doc`: show it below. A small file is used whole at once."""
        self._drop_source(self.sources[i])
        editor = FragmentEditor(doc, path or label)
        view = ReadOnlyView(editor)
        file_primers = []   # the primers stored in the file are looked for too, like any primer of the database
        for f in doc.features():
            if f.kind == "primer_bind":
                bases = "".join(doc.sequence(s, e) for s, e in f.spans)
                file_primers.append((f.label, primer_sequence(f) or (rc(bases) if f.strand < 0 else bases)))
        editor.set_extra(file_primers)
        self.sources[i] = SimpleNamespace(doc=doc, label=label, path=path, editor=editor, view=view, file_primers=file_primers)
        self.seqstack.addWidget(view)
        # Bound methods, found again by the sender: a lambda holding `src` would keep the editor and the panel alive for ever.
        editor.selection_changed.connect(self._on_selection_changed)
        editor.primers_changed.connect(self._on_primers_changed)
        self.strip.current = i
        self._show_source(i)
        if doc.length <= WHOLE_FILE_MAX:
            self._use(i, 0, doc.length)
            self.chosen.setText(f"Fragment {i + 1} uses the whole file ({doc.length:,} bp). Select a stretch, or click two primers, "
                                "below and press Use to change that.")
        else:
            self.slots[i] = None
            self.analyse()
            self.chosen.setText(f"{label} is {doc.length:,} bp: select the stretch to use, then press Use.")

    # -- the sequence view
    def _pick(self, i):
        self.strip.current = i
        self._show_source(i)

    def _source_of(self, editor):
        return next((s for s in self.sources if s is not None and s.editor is editor), None)

    def _on_selection_changed(self):
        src = self._source_of(self.sender())
        if src is not None:
            self._selection_changed(src)

    def _on_primers_changed(self):
        src = self._source_of(self.sender())
        if src is not None:
            self._fill_primer_lists(src)

    def _exact_sites(self, src, strand):
        return [s for s in src.editor.primer_sites if s.strand == strand and not s.mismatches]

    def _fill_primer_lists(self, src):
        """The forward and reverse lists: every primer that binds this sequence exactly."""
        if self.sources[self.current if self.current >= 0 else 0] is not src:
            return
        for box, strand in ((self.forward, 1), (self.reverse, -1)):
            keep = box.currentData().key if box.currentData() is not None else None
            box.clear()
            for s in self._exact_sites(src, strand):
                box.addItem(f"{s.label}  {s.start + 1:,}–{s.end:,}", s)
                if s.key == keep:
                    box.setCurrentIndex(box.count() - 1)
        self.by_primers.setEnabled(bool(self.forward.count() and self.reverse.count()))

    def _drop_product(self):
        if self.product_view is not None:
            self.seqstack.removeWidget(self.product_view)
            self.product_view.deleteLater()
        self.product_view, self.product_shown, self.product_tile.selected = None, False, False
        self.product_tile.update()

    def _product_ready(self):
        """The preview needs a design whose every junction fits (the product is then well defined, checks aside)."""
        d = self.design
        return d is not None and all(j.kind != "none" for j in d.joints)

    def _show_product(self):
        """The exact product below, like a fragment: each fragment a coloured bar, each homology / overhang / toehold amber."""
        if not self._product_ready():
            return
        self._drop_product()
        doc = simulate(self.slots, self.design, preview=True, tails=False)   # (the tails are for the junction inspector)
        self.product_view = ReadOnlyView(Editor(doc, "Assembly product"))
        self.seqstack.addWidget(self.product_view)
        self.seqstack.setCurrentWidget(self.product_view)
        self.product_shown, self.product_tile.selected = True, True
        self.product_tile.update()
        self.which.setText(f"<b>Product</b> · {doc.length:,} bp {'circular' if doc.circular else 'linear'} · {len(self.slots)} fragments")
        for w in (self.use, self.by_region, self.by_primers, self.forward, self.reverse):
            w.setEnabled(False)   # they are about choosing a fragment; _show_source turns them back on
        self.primer_info.setVisible(False)
        self.chosen.setText("Each fragment is a coloured bar named after it; the homologies (Gibson), overhangs (Golden Gate) or toeholds "
                            "(Sidewinder) that join them are slate. The primers that made the fragments are drawn as primers: hover or "
                            "click one to see where it is. Click a junction on the map to see what its primers added.")

    def _show_source(self, i):
        if self.product_shown:
            self._drop_product()
            for w in (self.by_region, self.forward, self.reverse):
                w.setEnabled(True)
        src = self.sources[i] if 0 <= i < len(self.sources) else None
        f = self.slots[i] if 0 <= i < len(self.slots) else None
        self.which.setText(f"<b>Fragment {i + 1}</b>" + (f" · {src.label} · {src.doc.length:,} bp" if src else " · choose a source"))
        self.use.setEnabled(src is not None)
        if src is None:
            self.seqstack.setCurrentWidget(self.nothing)
            self.forward.clear()
            self.reverse.clear()
            self.by_primers.setEnabled(False)
            self.chosen.setText("")
        else:
            self.seqstack.setCurrentWidget(src.view)
            self._refresh_primers()
            self._fill_primer_lists(src)
            if f is not None:
                ed = src.editor
                ed.select(f.start, f.end)
                ed.reveal.emit(f.start)
                for box, p in ((self.forward, f.fwd), (self.reverse, f.rev)):
                    for k in range(box.count()):
                        if p and box.itemData(k).label == p[0]:
                            box.setCurrentIndex(k)
        self.buttons_update()

    def _refresh_primers(self):
        """The primers of the selected fragment are drawn by the view (arrow, tail beside it; teal and "temporary" until they are in
        the database, details on hover). Only a Sidewinder fragment's oligos, which bind nowhere, are listed in a box."""
        i = self.current
        f = self.slots[i] if 0 <= i < len(self.slots) else None
        src = self.sources[i] if 0 <= i < len(self.sources) else None
        if src is not None:
            known = {r["name"]: r["sequence"].upper() for r in db().rows}
            designed = [(p[0], p[1]) for p in ((f.fwd, f.rev) if f else ()) if p]
            src.editor.set_extra(src.file_primers, [p for p in designed if known.get(p[0]) != p[1].upper()])
        oligos = self.method.currentIndex() == 2
        self.primer_info.setVisible(oligos)
        self.primer_info.setHtml(self._primers_html(f) if oligos else "")

    def _primers_html(self, f):
        if f is None or not (self.design and self.method.currentIndex() == 2):
            return ""
        rows = [r for r in self.design.rows if r[0].startswith(f.name + "-")]
        note = "Sidewinder oligos of this fragment (the product is then amplified with its own primer pair)."
        out = [f"<span style='color:gray'>{escape(note)}</span><br>"] if note else []
        for name, seq, anneal, _ in rows:
            # A Sidewinder oligo is not a primer: only its hairpin and self-dimer mean something.
            checks = [c for c in oligo_checks((name, seq, anneal, "")) if anneal or c[0] in ("Hairpin", "Self-dimer")]
            tail, core = (seq[:len(seq) - anneal], seq[len(seq) - anneal:]) if anneal else ("", seq)
            t = melting_temp(core) if anneal else None
            level = worst(checks)
            bad = "; ".join(f"{n} {tx}" for n, lv, tx in checks if lv in ("warn", "bad"))
            out.append(f"<tt><b>{escape(name)}</b> 5′-<span style='color:#bf360c'>{escape(tail.lower())}</span>{escape(core.upper())}-3′</tt> "
                       f"&nbsp;{len(seq)} nt" + (f" · Tm {t:.1f} °C" if t else "") +
                       f" <b><font color='{MARKS[level][1]}'>{MARKS[level][0]}</font></b>"
                       + (f" <span style='color:gray'>{escape(bad)}</span>" if bad else "") + "<br>")
        return "".join(out)

    def _selection_changed(self, src):
        """Something was selected in a source's view: a stretch, or primers (click one, or shift-click a second one on the
        other strand to select what they amplify, as in the editor)."""
        if self.sources[self.current if self.current >= 0 else 0] is not src:
            return
        ed, n = src.editor, src.doc.length
        site = ed.site_by_key.get(ed.current_primer)
        pair = None
        if ed.end > ed.start:
            fw = next((s for s in self._exact_sites(src, 1) if s.start == ed.start % n), None)
            rv = next((s for s in self._exact_sites(src, -1) if s.end == ed.end or s.end == ed.end - n * (ed.end > n)), None)
            pair = (fw, rv) if fw and rv else None
        if pair:
            for box, s in ((self.forward, pair[0]), (self.reverse, pair[1])):
                for k in range(box.count()):
                    if box.itemData(k).key == s.key:
                        box.setCurrentIndex(k)
            self.by_primers.setChecked(True)
            self.chosen.setText(f"Primers {pair[0].label} + {pair[1].label} selected: press Use to make the PCR product.")
        elif site is not None and not site.mismatches:
            box = self.forward if site.strand > 0 else self.reverse
            for k in range(box.count()):
                if box.itemData(k).key == site.key:
                    box.setCurrentIndex(k)
            self.chosen.setText(f"Primer {site.label} selected: shift-click one on the other strand for the PCR product.")
        else:
            self.by_region.setChecked(True)
            self.chosen.setText(f"Selected {ed.start + 1:,}–{ed.end:,} ({ed.end - ed.start:,} bp)" if ed.end > ed.start else
                                "Nothing selected: drag over the bases, or click a feature.")

    def use_selection(self):
        i = self.current
        src = self.sources[i] if 0 <= i < len(self.sources) else None
        if src is None:
            return
        L = src.doc.length
        if self.by_primers.isChecked():
            fw, rv = self.forward.currentData(), self.reverse.currentData()
            if fw is None or rv is None:
                self.chosen.setText("Pick a forward and a reverse primer.")
                return
            end = rv.end if rv.end > fw.start else rv.end + L if src.doc.circular else None
            if end is None:
                self.chosen.setText("The reverse primer lies before the forward primer on this linear sequence.")
                return
            self._use(i, fw.start, end, (fw.label, fw.sequence, fw.end - fw.start), (rv.label, rv.sequence, rv.end - rv.start))
        else:
            a, b = src.editor.start, src.editor.end
            if b <= a:
                self.chosen.setText("Select a stretch first (drag over the bases, or click a feature).")
                return
            self._use(i, a, b)

    def _use(self, i, start, end, fwd=None, rev=None):
        src, old = self.sources[i], self.slots[i]
        stem = Path(src.label).stem
        name = self._free_name(stem if (start, end) == (0, src.doc.length) else f"{stem} {start + 1}-{end}", i)
        f = Frag(name, src.doc, start, end, fwd, rev, vector=bool(old and old.vector), source=src.label)
        self.slots[i] = f
        self.analyse()
        self.chosen.setText(f"Fragment {i + 1} = {src.label} {start + 1:,}–{end:,} ({f.doc.length:,} bp"
                            + (", PCR with " + " + ".join(p[0] for p in (fwd, rev)) if fwd else "") + ")")

    # -- editing the list
    def move(self, step):
        self._move_to(self.current, self.current + step)

    def _move_to(self, i, j):
        if 0 <= i < len(self.slots) and 0 <= j < len(self.slots) and i != j:
            for lst in (self.slots, self.sources):
                lst.insert(j, lst.pop(i))
            self.strip.current = j
            self.analyse()
            self._show_source(j)

    def flip(self):
        self._flip_at(self.current)

    def _flip_at(self, i):
        if 0 <= i < len(self.slots) and self.slots[i]:
            self.slots[i].flip()
            self.analyse()

    def toggle_vector(self):
        self._toggle_vector_at(self.current)

    def _toggle_vector_at(self, i):
        f = self.slots[i] if 0 <= i < len(self.slots) else None
        if f:
            f.vector = not f.vector
            if f.vector:
                f.fwd = f.rev = None   # a vector is cut, not amplified
            f.rebuild()
            self.analyse()

    def clear(self):
        self._clear_at(self.current)

    def _clear_at(self, i):
        if 0 <= i < len(self.slots):
            self._drop_source(self.sources[i])
            self.slots[i] = self.sources[i] = None
            self._show_source(i)
            self.analyse()

    # -- analysis
    def params(self):
        return dict(method=self.method.currentIndex(), circular=self.is_circular(), enzyme=self.enzyme.currentText(),
                    min_ov=self.min_ov.value(), min_tm=self.ov_tm.value(), tm=self.anneal_tm.value())

    def analyse(self):
        """Check the junctions of the fragments as they are now and show everything that follows from it."""
        self.design, self.error = None, ""
        empty = [i + 1 for i, f in enumerate(self.slots) if f is None]
        hints = ["" for _ in self.slots]
        for i, s in enumerate(self.sources):
            if s and self.slots[i] is None:
                hints[i] = f"{s.label}: select a stretch below"
        cards, fchecks = [], {}
        if not empty:
            try:
                self.design = analyse(self.slots, **self.params())
            except ValueError as e:
                self.error = str(e)
            if self.design is not None:
                cards, fchecks = junction_cards(self.design), self.design.fchecks
        ready = self._product_ready()
        self.strip.refresh(self.slots, self.is_circular(), cards, self.current if self.current >= 0 else 0, hints, fchecks)
        self._refresh_tile(ready)
        self._status(empty)
        self._refresh_primers()
        if self.product_shown:   # something changed under the preview: redraw it, or go back to the fragment if it can no longer be made
            self._show_product() if ready else self._show_source(self.current)
        self.buttons_update()

    def _refresh_tile(self, ready):
        """The product tile: each fragment's share of it, in its own colour, and the verdict of each junction."""
        d = self.design
        parts = [(len(d.seqs[i]) - a - b, FRAG_COLORS[i % len(FRAG_COLORS)]) for i, a, b in d.pieces] if ready else []
        self.product_tile.refresh(self.is_circular(), parts, d.jlevels if ready else [],
                                  f"{len(d.product()):,} bp {'circular' if d.circular else 'linear'}" if ready else None)

    def _status(self, empty):
        d = self.design
        notes = ""
        if self.error:
            text, color = self.error, "bad"
        elif empty:
            text, color = (f"{plural(len(empty), 'fragment')} still empty ({', '.join(map(str, empty))}): choose a source for each gray box.",
                           None)
        else:
            missing = sum(j.kind == "none" for j in d.joints)
            product = d.product()
            if d.ok:
                text, color = (f"✓ Every junction checks out: {len(product):,} bp {'circular' if d.circular else 'linear'} product, "
                               f"{len(d.rows)} oligos. Perform assembly is ready."), "ok"
            else:
                bad = [c for c, level, _ in d.checks if level == "bad"]
                text, color = (f"✗ {plural(missing, 'junction')} {'does' if missing == 1 else 'do'} not fit" if missing else "✗ A check failed") + \
                    (": press \"Design primers\" to give the ends what the method needs, or change the fragments."
                     if missing and self.method.currentIndex() != 2 else f" ({', '.join(bad[:3])})."), "bad"
            # What is not about one junction or one fragment (the circle, leftover sites...); the rest is on the cards and boxes.
            notes = "<br>".join(f"<span style='color:{LEVEL_COLORS[lv]}'>{MARKS[lv][0]} <b>{escape(n)}</b>: {escape(t)}</span>"
                                for n, lv, t in d.global_checks if lv in ("warn", "bad"))
        self.status.setText(text)
        self.status.setStyleSheet(f"color: {LEVEL_COLORS[color]}; font-weight: 600;" if color else "")
        self.notes.setText(notes)

    def buttons_update(self):
        filled = [f for f in self.slots if f]
        m = self.method.currentIndex()
        self.buttons["homology"].setEnabled(m != 2 and len(filled) == len(self.slots))
        self.buttons["perform"].setEnabled(bool(self.design and self.design.ok))
        self.buttons["save"].setEnabled(bool(filled))
        self.buttons["primers"].setEnabled(bool(self.design and self.design.rows))
        self.report_button.setEnabled(bool(self.design))

    # -- actions
    def add_homology(self):
        if None in self.slots or self.method.currentIndex() == 2:
            return
        p = self.params()
        n = auto_primers(self.slots, p["method"], p["circular"], p["enzyme"], self.min_ov.value(), self.max_ov.value(), p["min_tm"], p["tm"])
        self.analyse()
        self.chosen.setText(f"New primers designed for {plural(n, 'fragment')} (any primers chosen before were replaced). They are drawn in teal on "
                            "each fragment's sequence as temporary primers: hover one for its details; Save primers to database keeps them.")
        self._show_source(self.current)

    def perform(self):
        d = self.design
        if d and d.ok:
            self.win.open_product(simulate(self.slots, d))

    def show_report(self):
        if self.design:
            dialog = ReportDialog(self, self.design)
            self.dialogs.append(dialog)
            dialog.show()

    def _show_junction(self, j):
        if self.design and j < len(self.design.joints):
            dialog = JunctionDialog(self, self.slots, self.design, j)
            self.dialogs.append(dialog)
            dialog.show()

    def save_fragments(self):
        folder = QFileDialog.getExistingDirectory(self, "Save the fragments in")
        if not folder:
            return
        d, saved = self.design, 0
        todo = [(i, f, Path(folder) / f"{runner.safe_name(f.name)}{'_backbone' if f.vector else ''}.gb")
                for i, f in enumerate(self.slots) if f is not None and not (f.vector and not d)]
        there = [path.name for _, _, path in todo if path.exists()]
        if there and QMessageBox.question(   # asked once, for all of them; No by default (the files would be lost)
                self, "Save fragments", f"{plural(len(there), 'file')} already {'exists' if len(there) == 1 else 'exist'} in {folder}:"
                f"\n\n{', '.join(there[:8])}{' …' if len(there) > 8 else ''}\n\nReplace {'it' if len(there) == 1 else 'them'}?",
                QMessageBox.Yes | QMessageBox.No, QMessageBox.No) != QMessageBox.Yes:
            return 0
        for i, f, path in todo:
            doc = piece_doc(self.slots, d, i) if f.vector else f.doc
            try:
                doc.save_genbank(path)
            except OSError as e:   # a read-only folder, a full disk: said, and the files written so far are counted
                QMessageBox.critical(self, "Save fragments", f"{path}\n\n{e}")
                break
            saved += 1
        self.status.setText(f"{plural(saved, 'fragment')} written to {folder} (PCR products with their tails; a vector as its backbone).")
        self.status.setStyleSheet("")
        return saved

    def save_primers(self):
        d = self.design
        if not d:
            return
        counts = {"added": 0, "exists": 0, "conflict": 0, "invalid": 0}
        for name, seq, _, kind in d.rows:
            counts[db().add(name, seq, kind)] += 1
        self.status.setText(f"Database: {counts['added']} added, {counts['exists']} already there, "
                            f"{counts['conflict']} skipped (name taken by another sequence)"
                            + (f", {counts['invalid']} refused (not a valid primer: a name or bases the database cannot take)"
                               if counts["invalid"] else ""))
        self.status.setStyleSheet("")
        self._refresh_primers()   # the ones now in the database are regular primers
