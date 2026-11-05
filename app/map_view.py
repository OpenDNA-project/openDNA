"""Plasmid map: circular for circular records, linear otherwise (a linear one is centred in the pane).

A FASTA file of several sequences (a genome's chromosomes) adds a list above the map, one line per sequence with a bar
as long as the sequence is; a click on a line opens that sequence, the wheel scrolls a long list (hundreds of contigs).

Features are drawn once into a cached pixmap (re-drawn only when the document changes
or the widget is resized); selection and cursor are painted on top each time.

A crowded map keeps its largest features: the smallest ones are left out until the rest are few
enough to read (one per SPACING px of backbone) and fit in lanes still thick enough to see (a note
says how many). Enlarging the pane brings them back; the sequence view always has them all.

A linear map zooms: the wheel zooms at the pointer (a slider with + / - at the right edge does too), Shift+wheel, a sideways
wheel or a middle-button drag pans. Its zoom is the map's own: the sequence view's letter size is not touched. Only the features
in view are laid out, so a zoomed-in region shows the small ones the whole map had to leave out.

A position scale runs inside the features (under them on a linear map): ticks every 1, 2 or 5 x 10^k bp, the step chosen
from the map's size so its labels never crowd, from a 100 bp oligo (every 10 bp) to a 10 Mb genome (every 1 or 2 Mb).
"""
import math

from PySide6.QtCore import QPoint, QPointF, QRectF, Qt, QTimer
from PySide6.QtGui import (
    QColor, QContextMenuEvent, QFont, QFontMetricsF, QPainter, QPainterPath, QPalette, QPen, QPixmap, QPolygonF,
)
from PySide6.QtWidgets import QApplication, QSlider, QToolButton, QVBoxLayout, QWidget

import image_export
import units
import theme
from editor import arrow, pack_lanes, text_color_on
from hover import HoverTips
from perf import timed
from sequence_view import ZOOM_BAR_W, WheelSteps, _Magnifier, edit_menu, primer_color
from uikit import plural

MAX_MARKS = 3000   # search hits drawn on the map
MIN_BAND = 5.0     # px: the thinnest feature lane worth drawing; fewer features rather than thinner lanes
CIRCLE = 0.8       # the circle's diameter, as a share of the room it has (width or height, the smaller)
MIN_VIEW = 50       # bp: the deepest zoom of a linear map
SPACING = 20.0     # px of backbone per feature drawn, on average: a map holds at most its length / SPACING features


def tick_step(bp_per_px, spacing):
    """The smallest 1, 2 or 5 x 10^k bp (at least 1) that puts ticks `spacing` px or more apart."""
    raw = max(1.0, bp_per_px * spacing)
    k = 10 ** math.floor(math.log10(raw))
    return next(int(m * k) for m in (1, 2, 5, 10) if m * k >= raw)


def minor_step(step):
    """Unlabelled ticks between the labelled ones: 1000 -> 200, 2000 -> 500, 5000 -> 1000 (1 bp at least)."""
    return max(1, step // (4 if str(step)[0] == "2" else 5))


def tick_label(pos, step):
    """A tick's position, in a unit that suits the step: 750, 2,500, 15 kb, 2.5 Mb (and 1, the first base, always 1)."""
    if pos == 1:
        return "1"
    if step >= 100_000:
        return f"{pos / 1e6:g} Mb"
    if step >= 5_000:
        return f"{pos / 1e3:g} kb"
    return f"{pos:,}"


def ticks(length, length_px, fm, last=None):
    """(labelled step, minor step or None) for a backbone of `length` bp drawn `length_px` px long: labels a little
    more than their own width apart, minor ticks only when 4 px or more apart. `last`: the biggest position labelled
    when that is not `length` (a zoomed map shows bases a..b: length b - a, last b)."""
    bp_per_px = max(1, length) / max(1.0, length_px)
    step = tick_step(bp_per_px, 60)
    widest = fm.horizontalAdvance(tick_label((last or length) // step * step, step))
    step = tick_step(bp_per_px, max(40.0, widest + 18))
    minor = minor_step(step)
    return step, (minor if minor < step and minor / bp_per_px >= 4 else None)


def fm_font(widget):
    """The map's small font: labels and the position scale."""
    font = QFont(widget.font())
    font.setPointSizeF(max(7.0, widget.font().pointSizeF() - 1))
    return font


def note_text(count):
    """The note of a crowded map: how many small features it leaves out."""
    return f"{plural(count, 'small feature')} not shown (enlarge the map to see more)"


def QLineF_length(a, b):
    return math.hypot(a.x() - b.x(), a.y() - b.y())


class MapView(QWidget):
    def __init__(self, editor):
        super().__init__()
        self.editor = editor
        self.cache = None
        self.layout_cache = None   # (key, layout): see _layout
        self.anchor = None
        self.view = None     # zoomed linear map: (first base, end) in view; None = the whole sequence
        self.pan = None      # middle-button drag: (x at the press, first base at the press)
        self.scroll = None   # first chromosome line shown (None: put the open sequence in view at the first paint)
        self.site_hits = []   # where the last paint drew each enzyme (View > Restriction Sites): (segment or rect, key)
        self.press = None    # where the left button went down on the backbone; a selection starts once it moves far enough
        self._zoom_wheel, self._list_wheel = WheelSteps(), WheelSteps()
        self.setMinimumSize(250, 200)
        # Big sequences: re-render the (many-feature) map once typing pauses, not per key.
        self.timer = QTimer(self, singleShot=True, interval=400, timeout=self._invalidate)
        editor.changed.connect(self._edited)   # a bound method, not a lambda: Qt drops it when this pane is closed
        editor.display_changed.connect(self._invalidate)
        editor.primers_changed.connect(self._invalidate)
        editor.selection_changed.connect(self.update)
        editor.search_changed.connect(self.update)
        editor.reveal.connect(self._reveal)
        editor.navigate.connect(self._reveal)
        self._build_zoom_bar()
        self.hover = HoverTips(self, editor, lambda pt: self._site_at(QPointF(pt)) or self._hit(QPointF(pt))[1])

    # ---- zoom (linear maps) ----------------------------------------------------------
    def _build_zoom_bar(self):
        """Magnifier, "+", slider, "-" down the right edge, as in the sequence view; the slider goes up to zoom in."""
        bar = self.zoom_bar = QWidget(self)
        bar.setFixedWidth(ZOOM_BAR_W)
        column = QVBoxLayout(bar)
        column.setContentsMargins(0, 4, 0, 4)
        column.setSpacing(2)
        column.addWidget(_Magnifier(), 0, Qt.AlignHCenter)
        self.zoom_slider = QSlider(Qt.Vertical)
        self.zoom_slider.setRange(0, 100)
        self.zoom_slider.setFixedHeight(90)
        self.zoom_slider.setToolTip("Zoom (wheel)")
        self.zoom_slider.valueChanged.connect(self._slider_zoom)
        buttons = []
        for text, factor in (("+", 1 / 1.5), ("\u2212", 1.5)):
            button = QToolButton()
            button.setText(text)
            button.setAutoRaise(True)
            button.setToolTip("Zoom in" if factor < 1 else "Zoom out")
            button.clicked.connect(lambda _=False, f=factor: self._zoom(f))
            buttons.append(button)
        column.addWidget(buttons[0], 0, Qt.AlignHCenter)
        column.addWidget(self.zoom_slider, 0, Qt.AlignHCenter)
        column.addWidget(buttons[1], 0, Qt.AlignHCenter)
        bar.adjustSize()
        self._place_zoom_bar()

    def _place_zoom_bar(self):
        self.zoom_bar.setVisible(not self.editor.doc.circular and self.editor.length > MIN_VIEW and not image_export.ACTIVE)
        self.zoom_bar.move(self.width() - ZOOM_BAR_W, self._top() + 4)

    def export_hidden(self):
        """What an exported image leaves out (image_export): the zoom bar."""
        return [self.zoom_bar]

    def _win(self):
        """(first base, end) of the part of a linear map in view, always inside the sequence."""
        n = max(1, self.editor.length)
        if self.view is None or self.editor.doc.circular:
            return 0.0, float(n)
        span = min(self.view[1] - self.view[0], n)
        a = max(0.0, min(self.view[0], n - span))
        return a, a + span

    def _set_view(self, a, b):
        n = max(1, self.editor.length)
        span = max(min(MIN_VIEW, n), min(b - a, n))
        a = max(0.0, min(a, n - span))
        self.view = None if span >= n else (a, a + span)
        self._sync_slider()
        self._invalidate()

    def _sync_slider(self):
        a, b = self._win()
        n = max(1, self.editor.length)
        value = 0 if n <= MIN_VIEW else round(100 * math.log(n / (b - a)) / math.log(n / MIN_VIEW))
        self.zoom_slider.blockSignals(True)
        self.zoom_slider.setValue(max(0, min(100, value)))
        self.zoom_slider.blockSignals(False)

    def _zoom(self, factor, at=None):
        """Zoom by `factor` (< 1 in) keeping the base `at` (default the middle of the view) where it is."""
        a, b = self._win()
        at = (a + b) / 2 if at is None else at
        keep = (at - a) / (b - a)
        span = (b - a) * factor
        self._set_view(at - keep * span, at - keep * span + span)

    def _slider_zoom(self, value):
        a, b = self._win()
        n = max(1, self.editor.length)
        span = n * (min(MIN_VIEW, n) / n) ** (value / 100)
        mid = (a + b) / 2
        self._set_view(mid - span / 2, mid + span / 2)

    def _reveal(self, pos):
        """Another view went to base `pos`: bring it into a zoomed map."""
        a, b = self._win()
        if self.view is not None and not a <= pos < b:
            self._set_view(pos - (b - a) / 2, pos + (b - a) / 2)

    def wheelEvent(self, event):
        if self._rows() and event.position().y() < self._top():
            self._first_row()
            self.scroll -= self._list_wheel(event.angleDelta().y())
            self.update()
        elif not self.editor.doc.circular and self.editor.length > 1:
            a, b = self._win()
            x0, w, _, _ = self._linear()
            delta = event.angleDelta()
            if event.modifiers() & Qt.ShiftModifier or abs(delta.x()) > abs(delta.y()):   # pan: a tenth of the view per notch
                step = (delta.x() or delta.y()) / 120 * (b - a) * 0.1
                self._set_view(a - step, b - step)
            elif steps := self._zoom_wheel(delta.y()):   # one step per notch, however finely a touchpad sends it
                self._zoom(1.25 ** -steps, self._pos_at(event.position().x()))
            event.accept()
        else:
            event.ignore()

    def _edited(self):
        if self.editor.length > 200_000:
            self.timer.start()
        else:
            self._invalidate()

    def _invalidate(self):
        self.cache = None
        self._place_zoom_bar()
        self.update()

    def resizeEvent(self, event):
        self._invalidate()
        self._place_zoom_bar()

    def theme_changed(self):
        """Light / dark switch (theme.refresh_views): the cached map was drawn in the old colours."""
        self._invalidate()

    # ---- which features fit -------------------------------------------------------
    def _layout(self):
        """(features drawn, {feature index: lane}, lane count, how many left out) at this size. The largest features are
        kept: the smallest go until there are at most one per SPACING px of backbone and they pack into lanes of at least
        MIN_BAND px. A map that is not crowded draws them all."""
        ed = self.editor
        circular = ed.doc.circular
        if circular:
            r = self._radius()
            length_px, room = 2 * math.pi * r, r * 0.45      # lanes fill the outer 45 % of the radius
        else:
            length_px, room = self.width() - 40.0, self._room() - 80.0
        view = self._win() if not circular else None
        key = (self.width(), self._room(), circular, ed.length, view)
        # ed.map_features is a new list whenever the features or the filter change: the cache holds on to it.
        if self.layout_cache and self.layout_cache[0] == key and self.layout_cache[1] is ed.map_features:
            return self.layout_cache[2]
        features = ed.map_features
        count = ed.map_lane_count
        if view is not None and view != (0.0, float(max(1, ed.length))):   # zoomed in: the features in view only
            features = [f for f in features if any(e > view[0] and s < view[1] for s, e in f.spans)]
            count = 1 + max((ed.map_lanes[f.index] for f in features), default=-1)
        most = max(1, int(room / (MIN_BAND + 3)))
        budget = max(1, int(length_px / SPACING))

        if len(features) <= budget and count <= most:
            layout = (features, ed.map_lanes, count, 0)   # everything fits: the editor's packing
        else:
            big = sorted(features, key=lambda f: sum(e - s for s, e in f.spans), reverse=True)
            lo, hi = 0, min(len(big), budget)   # the most of the largest features that fit `most` lanes (binary search)
            while lo < hi:
                k = (lo + hi + 1) // 2
                if pack_lanes(big[:k])[1] <= most:
                    lo = k
                else:
                    hi = k - 1
            lanes, count = pack_lanes(big[:lo])
            layout = ([f for f in features if f.index in lanes], lanes, count, len(features) - lo)
        self.layout_cache = (key, ed.map_features, layout)
        return layout

    # ---- the chromosome list (a FASTA file of several sequences) --------------------
    def _records(self):
        return self.editor.fasta_records if len(self.editor.fasta_records) > 1 else []

    def row_h(self):
        """px per line of the chromosome list: its font's height and a little room (a fixed 17 px clipped larger fonts)."""
        return round(QFontMetricsF(fm_font(self)).height()) + 4

    def _rows(self):
        """How many lines of the list are shown: all of them, up to 40 % of the pane (at least two)."""
        n = len(self._records())
        return min(n, max(2, int(self.height() * 0.4 // self.row_h()))) if n and not image_export.ACTIVE else 0

    def _top(self):
        """Height of the list, 0 when the file holds one sequence (and in an exported image)."""
        return self._rows() * self.row_h() + 6 if self._rows() else 0

    def _room(self):
        """Height left for the map itself under the list."""
        return max(1, self.height() - self._top())

    def _first_row(self):
        rows, n = self._rows(), len(self._records())
        if self.scroll is None:
            names = [name for name, _ in self._records()]
            cur = names.index(self.editor.fasta_record) if self.editor.fasta_record in names else 0
            self.scroll = cur - rows // 2
        self.scroll = max(0, min(self.scroll, n - rows))
        return self.scroll

    def _paint_records(self, p):
        recs = self._records()
        rows = self._rows()
        if not rows:
            return
        first = self._first_row()
        font = fm_font(self)
        bold = QFont(font)
        bold.setBold(True)
        fm = QFontMetricsF(bold)
        name_w = min(150.0, max(fm.horizontalAdvance(n) for n, _ in recs) + 8)
        size_w = fm.horizontalAdvance(f"{max(l for _, l in recs):,} bp") + 8
        x0, x1 = 12 + name_w, self.width() - 22 - size_w
        longest = max(1, max(l for _, l in recs))
        dim = self.palette().color(QPalette.PlaceholderText)
        text = self.palette().color(QPalette.Text)
        hl = QColor(theme.selection_color())
        row_h = self.row_h()
        p.setClipRect(QRectF(0, 0, self.width(), self._top() - 2))
        p.fillRect(QRectF(0, 0, self.width(), self._top() - 2), self.palette().color(QPalette.Base))
        for i in range(rows):
            name, length = recs[first + i]
            y = 2 + i * row_h
            current = name == self.editor.fasta_record
            if current:
                shade = QColor(hl)
                shade.setAlpha(50)
                p.fillRect(QRectF(0, y, self.width(), row_h), shade)
            p.setFont(bold if current else font)
            p.setPen(text)
            p.drawText(QRectF(8, y, name_w, row_h), Qt.AlignVCenter | Qt.AlignLeft,
                       fm.elidedText(name, Qt.ElideRight, name_w - 4))
            p.setPen(Qt.NoPen)
            p.setBrush(hl if current else dim)
            p.drawRect(QRectF(x0, y + 4, max(2.0, (x1 - x0) * length / longest), row_h - 8))
            p.setFont(font)
            p.setPen(dim)
            p.drawText(QRectF(x1 + 6, y, size_w, row_h), Qt.AlignVCenter | Qt.AlignRight, f"{length:,} bp")
        if len(recs) > rows:   # a thin scroll bar at the right edge
            h = rows * row_h
            p.setPen(Qt.NoPen)
            p.setBrush(dim)
            p.drawRect(QRectF(self.width() - 5, 2 + h * first / len(recs), 3, max(8.0, h * rows / len(recs))))
        p.setClipping(False)
        p.setPen(QPen(dim, 1))
        p.drawLine(0, self._top() - 3, self.width(), self._top() - 3)

    def _record_at(self, pt):
        """Name of the sequence whose line of the list is under a widget point, else None."""
        rows, row_h = self._rows(), self.row_h()
        if not rows or pt.y() < 2 or pt.y() >= 2 + rows * row_h:
            return None
        return self._records()[self._first_row() + int((pt.y() - 2) // row_h)][0]

    # ---- geometry (circular) -------------------------------------------------------
    def _radius(self):
        return max(40.0, CIRCLE * min(self.width(), self._room()) / 2)

    def _circle(self):
        """Centre, backbone radius, and per-lane band width."""
        c = QPointF(self.width() / 2, self._top() + self._room() / 2)
        r = self._radius()
        lanes = max(1, self._layout()[2])
        band = max(3.0, min(10.0, r * 0.45 / lanes - 3))
        return c, r, band

    def _lane_radii(self, lane, r, band):
        outer = r - 6 - lane * (band + 3)
        return outer, outer - band

    def _angle(self, pos):
        """Qt angle in degrees (0 = 3 o'clock, counter-clockwise) with base 0 at 12 o'clock."""
        return 90 - 360 * pos / max(1, self.editor.length)

    # ---- geometry (linear) ---------------------------------------------------------
    def _linear(self):
        margin = 20
        _, _, lanes, left_out = self._layout()
        lanes = max(1, lanes)
        band = max(3.0, min(14.0, (self._room() - 80) / lanes - 3))
        block = 8 + lanes * (band + 3) + 2 + 28          # features, then the scale with its labels
        note = self._note_rect(left_out)   # the note of left-out features keeps its strip at the bottom
        room = self._room() - (note.height() + 4 if note is not None else 0)
        y = self._top() + max(40.0, (room - block) / 2)   # the whole map centred in what is left
        return margin, self.width() - 2 * margin - ZOOM_BAR_W, y, band  # x0, width, backbone y, band

    def _x(self, pos):
        """Widget x of base `pos` (kept near the map's edges: a zoomed-in map has huge coordinates off screen)."""
        x0, w, _, _ = self._linear()
        a, b = self._win()
        return max(x0 - 50.0, min(x0 + w + 50.0, x0 + w * (pos - a) / (b - a)))

    def _pos_at(self, x):
        """Base under widget x of a linear map (a float, inside the sequence)."""
        x0, w, _, _ = self._linear()
        a, b = self._win()
        return a + max(0.0, min(1.0, (x - x0) / w)) * (b - a)

    # ---- painting ------------------------------------------------------------------
    def paintEvent(self, event):
        if image_export.ACTIVE:   # an image export: draw at the export's resolution, not the screen-sized cache
            p = QPainter(self)
            p.fillRect(self.rect(), self.palette().color(QPalette.Base))
            p.setRenderHint(QPainter.Antialiasing)
            (self._paint_circular if self.editor.doc.circular else self._paint_linear)(p)
            self._paint_selection(p)
            p.end()
            return
        if self.cache is None or self.cache.size() != self.size() * self.devicePixelRatioF():
            with timed("map: draw (cached until an edit or resize)"):
                self.cache = self._render_static()
        p = QPainter(self)
        p.drawPixmap(0, 0, self.cache)
        p.setRenderHint(QPainter.Antialiasing)
        self._paint_selection(p)
        self._paint_records(p)
        p.end()

    def _render_static(self):
        dpr = self.devicePixelRatioF()
        pix = QPixmap(self.size() * dpr)
        pix.setDevicePixelRatio(dpr)
        pix.fill(self.palette().color(QPalette.Base))
        p = QPainter(pix)
        p.setRenderHint(QPainter.Antialiasing)
        if self.editor.doc.circular:
            self._paint_circular(p)
        else:
            self._paint_linear(p)
        p.end()
        return pix

    def _paint_circular(self, p):
        ed = self.editor
        c, r, band = self._circle()
        text = self.palette().color(QPalette.Text)
        p.setPen(QPen(text, 2))
        p.drawEllipse(c, r, r)

        label_font = QFont(self.font())
        label_font.setPointSizeF(max(7.0, self.font().pointSizeF() - 1))
        fm = QFontMetricsF(label_font)
        p.setFont(label_font)
        features, lanes, _, left_out = self._layout()
        note = self._note_rect(left_out)
        taken = [note] if note is not None else []  # label rectangles already drawn (to skip overlapping labels)
        # Bigger features first, so they win the label space.
        for f in sorted(features, key=lambda f: -sum(e - s for s, e in f.spans)):
            lane = lanes[f.index]
            outer, inner = self._lane_radii(lane, r, band)
            if inner < 10:
                continue  # too many lanes for this size; skip the innermost
            color = ed.colors[f.index]
            g = units.geometry(f, ed.length, ed.doc.circular)
            if g is not None:
                self._unit_arc(p, c, outer, inner, g, color.darker(130))
            else:
                p.setPen(color.darker(140))
                p.setBrush(color)
            for i, (s, e) in enumerate(f.spans if g is None else ()):
                # Arrow head at the 3' end: last part going clockwise (+), first part (−).
                way = arrow(f)
                tip = way if (way > 0 and i == len(f.spans) - 1) or (way < 0 and i == 0) else 0
                p.drawPath(self._arc_band(c, outer, inner, s, e, tip))
            # Label outside the backbone, at the middle of the first span.
            s, e = f.spans[0]
            mid = math.radians(self._angle((s + e) / 2))
            direction = QPointF(math.cos(mid), -math.sin(mid))
            anchor = c + direction * (r + 12)
            w, h = fm.horizontalAdvance(f.label), fm.height()
            rect = QRectF(anchor.x() if direction.x() >= 0 else anchor.x() - w, anchor.y() - h / 2, w, h)
            if self.rect().contains(rect.toRect()) and not any(rect.intersects(t) for t in taken):
                taken.append(rect.adjusted(-2, -1, 2, 1))
                p.setPen(QPen(color.darker(150), 1))
                p.drawLine(c + direction * (r + 2), c + direction * (r + 10))
                p.setPen(text)
                p.drawText(rect, Qt.AlignVCenter | Qt.AlignLeft, f.label)

        # Primer sites: a thin arc just outside the backbone (orange = with mismatches).
        p.setPen(Qt.NoPen)
        for site in ed.drawn_primers:
            p.setBrush(QColor(primer_color(site)))
            for s, e in site.spans:
                p.drawPath(self._arc_band(c, r + 7, r + 3, s, e))

        # Unique cutters: tick across the backbone, name outside where there is room.
        dim = self.palette().color(QPalette.PlaceholderText)
        hits = []   # becomes self.site_hits below (not while exporting an image: that is another size)
        for pos, name in ed.shown_sites():
            a = math.radians(self._angle(pos))
            d = QPointF(math.cos(a), -math.sin(a))
            p.setPen(QPen(dim, 1))
            p.drawLine(c + d * (r - 4), c + d * (r + 4))
            hits.append(((c + d * (r - 4), c + d * (r + 4)), ("site", name, pos)))
            anchor = c + d * (r + 30)
            w, h = fm.horizontalAdvance(name), fm.height()
            rect = QRectF(anchor.x() if d.x() >= 0 else anchor.x() - w, anchor.y() - h / 2, w, h)
            if self.rect().contains(rect.toRect()) and not any(rect.intersects(t) for t in taken):
                taken.append(rect.adjusted(-2, -1, 2, 1))
                p.drawLine(c + d * (r + 4), c + d * (r + 28))
                p.drawText(rect, Qt.AlignVCenter | Qt.AlignLeft, name)
                hits.append((rect, ("site", name, pos)))

        if not image_export.ACTIVE:
            self.site_hits = hits

        title = QFont(self.font())
        title.setBold(True)
        box = QRectF(c.x() - r, c.y() - 30, 2 * r, 60)
        name = f"{ed.doc.name}\n{ed.length:,} bp"
        title_rect = QFontMetricsF(title).boundingRect(box, Qt.AlignCenter, name)
        self._scale_circular(p, c, r - 6 - self._layout()[2] * (band + 3) - 1, fm, title_rect)
        p.setFont(title)
        p.setPen(text)
        p.drawText(box, Qt.AlignCenter, name)
        self._paint_left_out(p, left_out)

    def _scale_circular(self, p, c, radius, fm, keep_clear):
        """A thin ring just inside the innermost lane, ticks pointing in, the labelled ones named further in (the names
        that would touch another or `keep_clear`, the title, are left out). Position 0 is at 12 o'clock."""
        n = self.editor.length
        if radius < 20 or n < 2:
            return
        dim = self.palette().color(QPalette.PlaceholderText)
        p.setPen(QPen(dim, 1))
        p.setBrush(Qt.NoBrush)
        p.drawEllipse(c, radius, radius)
        step, minor = ticks(n, 2 * math.pi * radius, fm)
        for pos in range(minor, n, minor) if minor else ():
            if pos % step:
                a = math.radians(self._angle(pos))
                d = QPointF(math.cos(a), -math.sin(a))
                p.drawLine(c + d * radius, c + d * (radius - 3))
        p.setFont(fm_font(self))
        taken = [keep_clear]
        h = fm.height()
        for pos in range(0, n, step):
            a = math.radians(self._angle(pos))
            d = QPointF(math.cos(a), -math.sin(a))
            p.drawLine(c + d * radius, c + d * (radius - 6))
            if pos == 0:
                continue   # the origin: where the scale starts (12 o'clock), its "0" not needed
            label = tick_label(pos, step)
            w = fm.horizontalAdvance(label)
            # Centre of the label: inward from the tick by half the label's extent along that direction.
            centre = c + d * (radius - 9 - abs(d.x()) * w / 2 - abs(d.y()) * h / 2)
            rect = QRectF(centre.x() - w / 2, centre.y() - h / 2, w, h)
            corners = (rect.topLeft(), rect.topRight(), rect.bottomLeft(), rect.bottomRight())
            inside = max(QLineF_length(q, c) for q in corners) < radius - 4   # clear of the ring and its ticks
            if inside and not any(rect.intersects(t) for t in taken):
                taken.append(rect.adjusted(-3, -1, 3, 1))
                p.drawText(rect, Qt.AlignCenter, label)

    def _note_rect(self, count):
        """Where the note of _paint_left_out goes (bottom left corner), or None when there is none. The labels of a
        circular map keep out of it, and a linear map is centred above it: it is never written over a label."""
        if not count or image_export.ACTIVE:
            return None
        fm = QFontMetricsF(fm_font(self))
        return QRectF(4, self.height() - 4 - fm.height() - 2, fm.horizontalAdvance(note_text(count)) + 6, fm.height() + 4)

    def _paint_left_out(self, p, count):
        """A note in the bottom corner: how many small features the map leaves out (not in an exported image)."""
        rect = self._note_rect(count)
        if rect is None:
            return
        p.setFont(fm_font(self))
        p.setPen(self.palette().color(QPalette.PlaceholderText))
        p.drawText(rect.adjusted(2, 0, 0, -2), Qt.AlignLeft | Qt.AlignBottom, note_text(count))

    def _unit_arc(self, p, c, outer, inner, g, color):
        """A transcription unit on the circle: thin arc over the promoter, thick from +1, a radial T bar at the terminator
        (a small arrow head when it runs on). See units.py."""
        mid, band = (outer + inner) / 2, outer - inner
        n = max(1, self.editor.length)
        rect = QRectF(c.x() - mid, c.y() - mid, 2 * mid, 2 * mid)
        p.setBrush(Qt.NoBrush)
        for a, b, heavy in g.pieces:
            p.setPen(QPen(color, max(2.0, band * 0.55) if heavy else 1.0, Qt.SolidLine, Qt.FlatCap))
            p.drawArc(rect, round(self._angle(a) * 16), round(-360 * (b - a) / n * 16))

        def radial(pos, r0, r1):
            t = math.radians(self._angle(pos))
            d = QPointF(math.cos(t), -math.sin(t))
            return c + d * r0, c + d * r1

        p.setPen(QPen(color, 1.0))
        p.drawLine(*radial(g.tss, mid, outer + 1))   # +1
        if g.stopped:
            p.setPen(QPen(color, 2.0, Qt.SolidLine, Qt.FlatCap))
            p.drawLine(*radial(g.end, inner - 2, outer + 2))
        else:
            t = math.radians(self._angle(g.end))
            tangent = QPointF(math.sin(t), math.cos(t)) * (1 if g.strand > 0 else -1)   # clockwise for +
            tip, (o, i) = c + QPointF(math.cos(t), -math.sin(t)) * mid + tangent * 7, radial(g.end, outer + 2, inner - 2)
            p.setPen(Qt.NoPen)
            p.setBrush(color)
            p.drawPolygon(QPolygonF([tip, o, i]))

    def _arc_band(self, c, outer, inner, s, e, tip=0):
        """Closed path for the ring segment between two radii covering bases [s, e).
        `tip` 1 / -1: an arrow head at the clockwise / counter-clockwise end (the strand)."""
        a = self._angle(s)
        sweep = -max(0.3, 360 * (e - s) / max(1, self.editor.length))  # negative = clockwise
        out_rect = QRectF(c.x() - outer, c.y() - outer, 2 * outer, 2 * outer)
        in_rect = QRectF(c.x() - inner, c.y() - inner, 2 * inner, 2 * inner)
        mid = (outer + inner) / 2
        # Head length: about the band's thickness, at most 60 % of the segment (in degrees).
        head = min(-sweep * 0.6, math.degrees((outer - inner) * 0.9 / max(1.0, mid))) if tip else 0

        def point(angle, radius):
            t = math.radians(angle)
            return QPointF(c.x() + radius * math.cos(t), c.y() - radius * math.sin(t))

        path = QPainterPath()
        if tip < 0:  # head at the start (counter-clockwise end)
            path.moveTo(point(a, mid))
            path.lineTo(point(a - head, outer))
            path.arcTo(out_rect, a - head, sweep + head)
            path.arcTo(in_rect, a + sweep, -(sweep + head))
        else:
            body = sweep + head  # (head = 0 when there is no tip)
            path.arcMoveTo(out_rect, a)
            path.arcTo(out_rect, a, body)
            if tip:
                path.lineTo(point(a + sweep, mid))
            path.arcTo(in_rect, a + body, -body)
        path.closeSubpath()
        return path

    def _paint_linear(self, p):
        ed = self.editor
        x0, w, y, band = self._linear()
        text = self.palette().color(QPalette.Text)
        p.setPen(QPen(text, 2))
        p.drawLine(QPointF(x0, y), QPointF(x0 + w, y))
        a, b = self._win()
        shown = "" if self.view is None else f"   showing {int(a) + 1:,} \u2013 {int(b):,}"
        p.drawText(QRectF(x0, y - 36, w, 20), Qt.AlignLeft, f"{ed.doc.name}  ({ed.length:,} bp){shown}")
        p.setClipRect(QRectF(x0, self._top(), w, self._room()))   # a zoomed-in map: what is off to the side is not drawn
        fm = QFontMetricsF(self.font())
        p.setPen(self.palette().color(QPalette.PlaceholderText))
        free_x = -1.0
        hits = []   # becomes self.site_hits below (not while exporting an image: that is another size)
        for pos, name in ed.shown_sites():
            if not a - 1 <= pos <= b + 1:
                continue
            x = self._x(pos)
            p.drawLine(QPointF(x, y - 5), QPointF(x, y + 5))
            hits.append(((QPointF(x, y - 5), QPointF(x, y + 5)), ("site", name, pos)))
            if x >= free_x:  # skip names that would overlap
                p.drawText(QPointF(x + 2, y - 7), name)
                width = fm.horizontalAdvance(name)
                hits.append((QRectF(x + 2, y - 7 - fm.ascent(), width, fm.height()), ("site", name, pos)))
                free_x = x + width + 6
        if not image_export.ACTIVE:
            self.site_hits = hits
        p.setPen(Qt.NoPen)
        for site in ed.drawn_primers:  # thin marks just above the backbone
            p.setBrush(QColor(primer_color(site)))
            for s, e in site.spans:
                p.drawRect(QRectF(self._x(s), y - 5, max(2.0, self._x(e) - self._x(s)), 3))
        features, lanes, _, left_out = self._layout()
        p.setFont(self.font())
        for f in features:
            top = y + 8 + lanes[f.index] * (band + 3)
            color = ed.colors[f.index]
            g = units.geometry(f, ed.length, ed.doc.circular)
            if g is not None:
                self._unit_line(p, top, band, g, color.darker(130), f.label, fm)
                continue
            p.setPen(color.darker(140))
            p.setBrush(color)
            for i, (s, e) in enumerate(f.spans):
                rect = QRectF(self._x(s), top, max(1.0, self._x(e) - self._x(s)), band)
                head = min(band * 0.7, rect.width() * 0.6)  # arrow head at the 3' end (strand)
                x0, x1, y0, y1, ym = rect.left(), rect.right(), rect.top(), rect.bottom(), rect.center().y()
                if arrow(f) > 0 and i == len(f.spans) - 1:
                    pts = [(x0, y0), (x1 - head, y0), (x1, ym), (x1 - head, y1), (x0, y1)]
                elif arrow(f) < 0 and i == 0:
                    pts = [(x0 + head, y0), (x1, y0), (x1, y1), (x0 + head, y1), (x0, ym)]
                else:
                    pts = [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]
                p.drawPolygon(QPolygonF([QPointF(*pt) for pt in pts]))
                if band >= fm.height() - 2 and rect.width() > fm.horizontalAdvance(f.label) + 4:
                    p.setPen(text_color_on(color))
                    p.drawText(rect, Qt.AlignCenter, f.label)
                    p.setPen(color.darker(140))
        p.setClipping(False)
        self._paint_left_out(p, left_out)   # (after the clip: it starts left of the backbone)
        self._scale_linear(p, x0, w, y + 8 + self._layout()[2] * (band + 3) + 2)

    def _scale_linear(self, p, x0, w, top):
        """A thin rule under the lowest lane, ticks pointing down, the labelled ones named under them (from 1 at the left)."""
        n = self.editor.length
        if n < 2:
            return
        a, b = self._win()
        font = fm_font(self)
        fm = QFontMetricsF(font)
        p.setFont(font)
        dim = self.palette().color(QPalette.PlaceholderText)
        p.setPen(QPen(dim, 1))
        p.drawLine(QPointF(x0, top), QPointF(x0 + w, top))
        step, minor = ticks(max(1, int(b - a)), w, fm, int(b))
        first = max(minor or 1, -int(-a // (minor or 1)) * (minor or 1))
        for pos in range(first, int(b), minor) if minor else ():
            if pos % step:
                p.drawLine(QPointF(self._x(pos), top), QPointF(self._x(pos), top + 3))
        free_x = -1e9
        for pos in ([1] if a <= 1 else []) + list(range(max(step, -int(-a // step) * step), int(b) + 1, step)):
            x = self._x(pos - 1 if pos == 1 else pos)
            p.drawLine(QPointF(x, top), QPointF(x, top + 6))
            label = tick_label(pos, step)
            lw = fm.horizontalAdvance(label)
            left = min(max(x - lw / 2, 0.0), self.width() - lw)
            if left >= free_x:
                p.drawText(QPointF(left, top + 7 + fm.ascent()), label)
                free_x = left + lw + 8

    def _unit_line(self, p, top, band, g, color, label, fm):
        """A transcription unit on the linear map (see _unit_arc)."""
        ym = top + band * 0.6
        p.setBrush(Qt.NoBrush)
        for a, b, heavy in g.pieces:
            p.setPen(QPen(color, max(2.0, band * 0.45) if heavy else 1.0, Qt.SolidLine, Qt.FlatCap))
            p.drawLine(QPointF(self._x(a), ym), QPointF(self._x(b), ym))
        p.setPen(QPen(color, 1.0))
        p.drawLine(QPointF(self._x(g.tss), ym), QPointF(self._x(g.tss), top))
        ex = self._x(g.end)
        if g.stopped:
            p.setPen(QPen(color, 2.0, Qt.SolidLine, Qt.FlatCap))
            p.drawLine(QPointF(ex, top), QPointF(ex, top + band))
        else:
            way = 1 if g.strand > 0 else -1
            p.setPen(Qt.NoPen)
            p.setBrush(color)
            p.drawPolygon(QPolygonF([QPointF(ex + way * 4, ym), QPointF(ex - way * 2, top), QPointF(ex - way * 2, top + band)]))
        if band >= fm.height() - 2:
            p.setPen(color)
            p.drawText(QPointF(self._x(g.tss) + 3 if g.strand > 0 else self._x(g.tss) - 3 - fm.horizontalAdvance(label),
                               top + fm.ascent() - 2), label)

    def _paint_hits(self, p):
        """What the search bar found: yellow marks over the backbone (the first MAX_MARKS; the counter has the total)."""
        ed = self.editor
        p.setPen(Qt.NoPen)
        p.setBrush(theme.hit_color())
        if ed.doc.circular:
            c, r, _ = self._circle()
            for s, e in ed.search_spans[:MAX_MARKS]:
                p.drawPath(self._arc_band(c, r + 12, r - 12, s, max(e, s + 1)))
        else:
            x0, w, y, _ = self._linear()
            for s, e in ed.search_spans[:MAX_MARKS]:
                p.drawRect(QRectF(self._x(s), y - 12, max(2.0, self._x(e) - self._x(s)), 24))

    def _paint_selection(self, p):
        ed = self.editor
        self._paint_hits(p)
        hl = QColor(theme.selection_color())
        if ed.doc.circular:
            c, r, _ = self._circle()
            if ed.end > ed.start:
                hl.setAlpha(90)
                p.setPen(Qt.NoPen)
                p.setBrush(hl)
                p.drawPath(self._arc_band(c, r + 8, r - 8, ed.start, ed.end))
            else:
                a = math.radians(self._angle(ed.start))
                d = QPointF(math.cos(a), -math.sin(a))
                p.setPen(QPen(hl, 2))
                p.drawLine(c + d * (r - 12), c + d * (r + 12))
        else:
            x0, w, y, _ = self._linear()
            hl.setAlpha(90)
            if ed.end > ed.start:
                p.fillRect(QRectF(self._x(ed.start), y - 10, max(2.0, self._x(ed.end) - self._x(ed.start)), 20), hl)
            else:
                hl.setAlpha(255)
                p.setPen(QPen(hl, 2))
                p.drawLine(QPointF(self._x(ed.start), y - 12), QPointF(self._x(ed.start), y + 12))

    # ---- mouse ---------------------------------------------------------------------
    def _hit(self, pt):
        """(base position, feature index or None) under a widget point."""
        ed = self.editor
        if ed.doc.circular:
            c, r, band = self._circle()
            dx, dy = pt.x() - c.x(), c.y() - pt.y()
            theta = math.degrees(math.atan2(dy, dx))
            pos = int(((90 - theta) / 360 % 1.0) * ed.length)
            dist = math.hypot(dx, dy)
            lane = int((r - 6 - dist) // (band + 3))
            inside_band = lane >= 0 and dist >= self._lane_radii(lane, r, band)[1]
        else:
            x0, w, y, band = self._linear()
            pos = int(self._pos_at(pt.x()))
            lane = int((pt.y() - y - 8) // (band + 3))
            inside_band = lane >= 0 and (pt.y() - y - 8) - lane * (band + 3) <= band
        pos = min(pos, ed.length)
        if inside_band:
            base = min(pos, ed.length - 1)
            for f in ed.visible_features(base, base + 1):
                if self._layout()[1].get(f.index) == lane:  # hidden gene parts and left-out features are not there
                    return pos, f.index
        return pos, None

    def _site_at(self, pt):
        """Key ("site", enzyme name, cut position) of the enzyme tick or name drawn under a widget point, else None."""
        for shape, key in self.site_hits:
            if isinstance(shape, QRectF):
                if shape.adjusted(-1, -1, 1, 1).contains(pt):
                    return key
            else:   # a tick: a few pixels either side of its segment count
                a, b = shape
                ab = b - a
                length2 = ab.x() ** 2 + ab.y() ** 2 or 1.0
                t = max(0.0, min(1.0, ((pt.x() - a.x()) * ab.x() + (pt.y() - a.y()) * ab.y()) / length2))
                if math.hypot(pt.x() - a.x() - t * ab.x(), pt.y() - a.y() - t * ab.y()) <= 3:
                    return key
        return None

    def _primer_at(self, pt):
        """Key of the primer site whose thin mark is under a widget point, else None."""
        ed = self.editor
        if ed.doc.circular:
            c, r, _ = self._circle()
            dx, dy = pt.x() - c.x(), c.y() - pt.y()
            if not r + 1 <= math.hypot(dx, dy) <= r + 9:
                return None
            pos, slack = ((90 - math.degrees(math.atan2(dy, dx))) / 360 % 1.0) * ed.length, ed.length * 3 / (2 * math.pi * r)
        else:
            x0, w, y, _ = self._linear()
            if not y - 8 <= pt.y() <= y - 1:
                return None
            a, b = self._win()
            pos, slack = self._pos_at(pt.x()), (b - a) * 3 / w
        for site in ed.drawn_primers:   # the marks are a few pixels: a click close to one counts
            if any(s - slack <= pos < e + slack for s, e in site.spans):
                return site.key
        return None

    def _menu_at(self, pt, global_pos):
        """The right-click menu at widget point `pt` (None: the Menu key): the primer's or the feature's under it, else that
        of the selected bases, else Cut / Copy / Paste / Select All."""
        site, (_, feature) = (self._primer_at(pt), self._hit(pt)) if pt is not None else (None, (None, None))
        if site is not None:
            from feature_table import primer_menu  # here: feature_table imports the tools
            primer_menu(self, self.editor, site, global_pos)
        elif feature is not None:
            from feature_table import feature_menu  # here: feature_table imports the tools
            feature_menu(self, self.editor, feature, global_pos)
        elif self.editor.end > self.editor.start:  # selected bases: primers, feature, gene, copy
            from feature_table import selection_menu
            selection_menu(self, self.editor, global_pos)
        else:
            edit_menu(self, self.editor, global_pos)

    def contextMenuEvent(self, event):
        """The Menu key: the menu of the selection (or Cut / Copy / Paste / Select All). A right click opened its menu when the
        button went down (mousePressEvent)."""
        if event.reason() != QContextMenuEvent.Mouse:
            self._menu_at(None, self.mapToGlobal(QPoint(self.width() // 2, self.height() // 2)))
        event.accept()

    def mousePressEvent(self, event):
        if event.button() == Qt.RightButton and event.position().y() < self._top():
            return
        if event.button() == Qt.RightButton:
            self._menu_at(event.position(), event.globalPosition().toPoint())
            return
        if event.button() == Qt.MiddleButton and not self.editor.doc.circular:
            self.pan = (event.position().x(), self._win()[0])
            self.setCursor(Qt.ClosedHandCursor)
            return
        if event.button() != Qt.LeftButton:
            return
        if event.position().y() < self._top():   # the chromosome list (switching replaces this pane: do it after the click)
            name = self._record_at(event.position())
            if name is not None and name != self.editor.fasta_record:
                QTimer.singleShot(0, lambda: self.editor.record_requested.emit(name))
            return
        pos, feature = self._hit(event.position())
        site = self._primer_at(event.position())
        enzyme = self._site_at(event.position())
        if enzyme is not None:   # an enzyme wins over a primer drawn under it, as in the sequence view
            self.editor.select_site(enzyme[1], enzyme[2])
            self.anchor = None
        elif site is not None:
            self.editor.select_primer(site, extend=bool(event.modifiers() & Qt.ShiftModifier))
            self.anchor = None
        elif feature is not None:
            self.editor.select_feature(feature)
            self.anchor = None
        else:
            self.anchor, self.press = pos, event.position()
            self.editor.set_cursor(pos)
            self.editor.reveal.emit(pos)

    def mouseDoubleClickEvent(self, event):
        _, feature = self._hit(event.position())
        if feature is not None:
            self.editor.edit_feature_requested.emit(feature)

    def mouseMoveEvent(self, event):
        if self.pan is not None:
            a, b = self._win()
            first = self.pan[1] - (event.position().x() - self.pan[0]) / self._linear()[1] * (b - a)
            self._set_view(first, first + b - a)
        elif self.anchor is not None and event.buttons() & Qt.LeftButton:
            # A click that wobbles by a pixel is still a click: the selection starts once the pointer has really moved.
            if self.press is not None and (event.position() - self.press).manhattanLength() < QApplication.startDragDistance():
                return
            self.press = None
            pos, _ = self._hit(event.position())
            self.editor.select(self.anchor, pos)
        elif not event.buttons():   # over a feature or a line of the list: a hand, as it can be clicked
            if event.position().y() < self._top():
                self.setCursor(Qt.PointingHandCursor)
                return
            self.setCursor(Qt.PointingHandCursor if self._hit(event.position())[1] is not None or self._primer_at(event.position()) is not None
                           or self._site_at(event.position()) is not None else Qt.ArrowCursor)

    def mouseReleaseEvent(self, event):
        if self.pan is not None:
            self.pan = None
            self.unsetCursor()   # the hover cursor comes back with the next move
        if self.anchor is not None and self.editor.end > self.editor.start:
            self.editor.reveal.emit(self.editor.start)
        self.anchor = self.press = None
