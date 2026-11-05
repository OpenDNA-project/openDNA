"""IGV-style genome view: a horizontal, zoomable strip from the whole sequence down to
single bases, with tracks: overview, ruler, GC %, bases (zoomed in), features, primer
sites, loaded annotation tracks (GFF3/GTF/BED/VCF) and BAM alignments (coverage + reads).

Mouse: wheel = zoom at the pointer, Shift+wheel or a sideways wheel = pan, drag = pan, Shift+drag = select, click = cursor
or feature, double-click = zoom in. The locus box takes "12,000-15,000", "12000" or a name. Many tracks scroll up and down
under the overview and the ruler, which stay at the top (the wheel over the track names scrolls them).

Everything is drawn from what is visible: GC comes from Rust in screen-sized bins (a
cached block sum when zoomed out), crowded tracks switch to a density plot. BAM coverage
of wide windows needs one pass over the file's reads (done in a background thread);
reads are drawn in narrower windows (Alignments settings), like IGV's visibility window.
Read colours, variant colouring and flags: see alignments.py.

Each track can be expanded, condensed (thin lanes, no labels) or hidden (click its name to
hide/show, right-click it for the menu), resized by dragging the thin line under it (taller
= more lanes of reads), and bar tracks (coverage, dense annotations) drawn on a log scale.
"""
import math
import re
from html import escape
from types import SimpleNamespace

from PySide6.QtCore import QObject, QPointF, QRectF, Qt, QThreadPool, QTimer, Signal
from PySide6.QtGui import QBrush, QColor, QFont, QFontMetricsF, QPainter, QPalette, QPen, QPixmap, QPolygonF
from PySide6.QtWidgets import (
    QFrame, QHBoxLayout, QLabel, QLineEdit, QMenu, QPushButton, QScrollArea, QSizePolicy, QVBoxLayout, QWidget,
)

import image_export
import theme
import tool_settings
from hover import HoverTips, Text
from perf import Laps, timed
from alignments import (COLORS, INVERTED_COLOR, LABELS, MATE_UNMAPPED_COLOR, TOOL as ALIGNMENTS, classify, flag_level, insert_limits,
                        relay as scan_relay, setting, start_scan)
from editor import arrow, text_color_on
from variants import annotate_variants, apply_variants, can_apply, change_text
from sequence_view import TEMP_COLOR, primer_color

LEFT = 74         # px for track names, at least (longer names widen the column: see GenomeCanvas.left)
MAX_LEFT = 160    # px: the widest the names column gets (longer names are cut short, the whole name in their tooltip)
MIN_SPAN = 40     # bp: deepest zoom
DENSE = 2500      # more items than this in view -> density plot
BASE_COLORS = {"A": "#2e7d32", "C": "#1565c0", "G": "#ef6c00", "T": "#c62828"}
TRACK_COLORS = ["#00897b", "#6d4c41", "#3949ab", "#c0ca33", "#8e24aa"]
READ_LANES = 20     # rows of reads at the default height (drag the track taller for more)
READ_COLORS = COLORS         # by category, see alignments.classify
INSERTION_COLOR = "#7b1fa2"
FLAG_COLOR = "#d32f2f"
GAP, HIDDEN_H = 6, 12        # px between tracks; height of a hidden track (its name only)
HITS_HEAD = 14               # px above the Local Align rows: the track's name (click: hide / show)
# kind -> (lanes, lane px expanded, lane px condensed, lanes condensed); 0 lanes = a bar plot
# of that height. Condensed reads: twice the rows in less room (IGV's "squished").
KINDS = {"gc": (0, 30, 16, 0), "features": (4, 13, 6, 4), "primers": (2, 9, 5, 2), "track": (3, 11, 5, 3),
         "coverage": (0, 30, 16, 0), "reads": (READ_LANES, 7, 3, 2 * READ_LANES), "hits": (1, 14, 8, 1),
         "variants": (2, 12, 5, 2)}
TEXT_KINDS = ("features", "track", "variants", "hits")   # lanes that write names in them: taller with a larger font
FONT_PX = 15.0      # the font height (px) the sizes above are made for; a larger font scales the lanes that hold text
MAX_SHARE = 0.6     # the alignment view asks for at most this share of the screen's height (its tracks scroll beyond)
MODES = ["expanded", "condensed", "hidden"]
SCALES = {"linear": lambda v: v, "log2": lambda v: math.log2(1 + v), "log10": lambda v: math.log10(1 + v)}
SCALED = ("coverage", "track")  # kinds with a scale choice


class _Relay(QObject):
    """Tells canvases (GUI thread) that a background coverage pass finished. Lives forever."""
    depth_ready = Signal()


_relay = _Relay()


def nice_step(raw):
    """1, 2 or 5 x 10^k just above `raw`."""
    k = 10 ** math.floor(math.log10(max(raw, 1)))
    return next(m * k for m in (1, 2, 5, 10) if m * k >= raw)


def format_bp(v, step=1):
    if step >= 1_000_000:
        return f"{v / 1e6:g} Mb"
    if step >= 1000:
        return f"{v / 1e3:g} kb"
    return f"{v:,}"


def pack(boxes, max_lanes):
    """Greedy lanes in pixel space: [(x0, x1, payload)] -> [(lane, x0, x1, payload)]."""
    ends, out = [], []
    for x0, x1, payload in sorted(boxes, key=lambda b: b[0]):
        lane = next((i for i, e in enumerate(ends) if e + 2 <= x0), None)
        if lane is None:
            if len(ends) < max_lanes:
                ends.append(x1)
                lane = len(ends) - 1
            else:
                lane = max_lanes - 1  # crowded: overlap in the last lane
        ends[lane] = max(ends[lane], x1)
        out.append((lane, x0, x1, payload))
    return out


class GenomeCanvas(QWidget):
    def __init__(self, panel):
        super().__init__()
        self.panel, self.editor = panel, panel.editor
        self.hits = []      # (QRectF, action) of clickable boxes from the last paint
        self.drag = None    # ("pan" | "select", x at press, view at press, moved?)
        self.resizing = None  # (track name, y at press, height at press)
        self.styles = {}      # track name -> SimpleNamespace(mode, height, scale)
        self.kinds = {}       # track name -> kind (see KINDS), from the last layout
        self.setMouseTracking(True)  # resize cursor over the lines between tracks
        # No keyboard focus: nothing here takes keys, and a click that took the focus left typed bases going nowhere (they
        # go on to the sequence view that had it).
        self.setFocusPolicy(Qt.NoFocus)
        self.left = LEFT        # px of the names column (see layout_rows)
        self.header_h = 44      # px of the overview and the ruler, kept at the top when the tracks scroll (see layout_rows)
        self.ruler_base = 30    # y of the ruler's labels (their baseline)
        self.layer, self.layer_key, self.dirty = None, None, True  # cached image of the tracks (see paintEvent)
        self.pan_timer = QTimer(self, singleShot=True, interval=120, timeout=self.invalidate)
        self.hover = HoverTips(self, self.editor, self.hover_target)
        _relay.depth_ready.connect(self.invalidate)
        scan_relay.variants_ready.connect(self.invalidate)
        tool_settings.notifier.changed.connect(self._settings_changed)
        self.read_kinds = {}

    def theme_changed(self):
        """Light / dark switch (theme.refresh_views): the cached tracks and the reads' colours are made again."""
        for a in self.editor.alignments:
            a.cache = None
        self.invalidate()

    def _settings_changed(self, tool):
        if tool == ALIGNMENTS:
            for a in self.editor.alignments:
                a.cache = None      # colours, insert limits
                a.variants = None   # thresholds: scan again
            self.invalidate()

    # ---- coordinates ---------------------------------------------------------------
    def width_px(self):
        return max(10.0, self.width() - self.left - 8)

    def x(self, pos):
        p = self.panel
        return self.left + (pos - p.vstart) / (p.vend - p.vstart) * self.width_px()

    def pos(self, x):
        p = self.panel
        return p.vstart + (x - self.left) / self.width_px() * (p.vend - p.vstart)

    def view_top(self):
        """How far the tracks are scrolled down (the canvas sits in a scroll area): the overview and the ruler are drawn at
        this y, so they stay in view."""
        return max(0, -self.y())

    def scale(self):
        """1 at the font the sizes are made for; more with a larger font (or screen scaling), so names fit their lanes."""
        return max(1.0, QFontMetricsF(self.font()).height() / FONT_PX)

    def name_h(self):
        """px of a track's name line (14 at the usual font)."""
        return max(14, int(QFontMetricsF(self.font()).height()))

    # ---- tracks: display mode, height, scale ---------------------------------------
    def style(self, name):
        """How a track is shown (kept per view, by track name, so a reloaded file keeps it)."""
        return self.styles.setdefault(name, SimpleNamespace(mode="expanded", height=None, scale="linear"))

    def lane_px(self, kind, mode):
        """px per lane of a kind of track in a mode (bar tracks: their whole height)."""
        _, expanded, condensed, _ = KINDS[kind]
        px = condensed if mode == "condensed" else expanded
        return round(px * self.scale()) if kind in TEXT_KINDS else px

    def lane_h(self, name):
        """px per lane (bar tracks: their whole height)."""
        return self.lane_px(self.kinds[name], self.style(name).mode)

    def lanes(self, name, h):
        return max(1, int(h // self.lane_h(name)))

    def track_height(self, name, kind):
        st = self.style(name)
        if st.mode == "hidden":
            return max(HIDDEN_H, self.name_h() - 2)
        lanes, _, _, condensed_lanes = KINDS[kind]
        px = self.lane_px(kind, st.mode)
        if kind == "hits":   # Local Align: one row per alignment (hit), under a line with the track's name
            rows = max(1, len(self.editor.query_hits))
            return st.height or self.hits_head() + rows * px
        default = px * max(1, condensed_lanes if st.mode == "condensed" else lanes)
        return st.height or default

    def hits_head(self):
        """px above the Local Align rows (the track's name)."""
        return round(HITS_HEAD * self.scale())

    def layout_rows(self):
        """y of each track: list of (name, top, height), and the bottom of the last one. Also sets the ruler's place and the
        names column's width from the font and the names (`ruler_base`, `header_h`, `left`)."""
        ed = self.editor
        fm = QFontMetricsF(self.font())
        self.ruler_base = 18 + max(12, math.ceil(fm.ascent()))   # labels on the ruler; ticks and its line under them
        self.header_h = self.ruler_base + 14
        rows, y = [("", 2, 12), ("", 18, self.header_h - 22)], self.header_h
        self.kinds = {}

        def add(name, kind):
            nonlocal y
            self.kinds[name] = kind
            h = self.track_height(name, kind)
            rows.append((name, y, h))
            y += h + GAP

        add("GC %", "gc")
        if self.panel.vend - self.panel.vstart <= self.width_px() / 7:
            rows.append(("Bases", y, self.name_h()))
            y += self.name_h() + GAP
        add("Features", "features")
        if ed.drawn_primers:
            add("Primers", "primers")
        if ed.query_hits:
            add("Query hits", "hits")
        for t in ed.tracks:
            add(t.name, "track")
        for a in ed.alignments:
            add(a.name, "coverage")
            add(a.name + " variants", "variants")
            add(a.name + " reads", "reads")
        # The names column: as wide as the longest name needs (at least LEFT, at most MAX_LEFT; a longer one is cut short)
        widest = max((fm.horizontalAdvance("▾ " + self.track_label(name)) for name, _, _ in rows if name), default=0)
        if ed.query_hits:
            widest = max(widest, fm.horizontalAdvance(self.hit_label(len(ed.query_hits))))
        self.left = int(max(LEFT, min(MAX_LEFT, widest + 12)))
        return rows, y

    def track_label(self, name):
        """What the names column says for a track (a BAM's reads: "reads", its file name is just above)."""
        return "reads" if self.kinds.get(name) == "reads" else name

    def shown(self, tops, name):
        """(top, height) of a track that is drawn, None if hidden."""
        return None if self.style(name).mode == "hidden" else tops[name]

    def in_header(self, y):
        """Is y on the overview or the ruler (they stay at the top of what shows, over the tracks scrolled under them)?"""
        return y < self.view_top() + self.header_h - 2

    def row_at(self, y):
        """Named track whose rows (or the gap under them) contain y: (name, top, height)."""
        rows = self.layout_rows()[0]
        if self.in_header(y):
            return None
        for name, top, h in rows:
            if name and top <= y < top + h + GAP:
                return name, top, h
        return None

    def edge_at(self, y):
        """Track whose bottom line is under y (drag it to resize), or None."""
        row = self.row_at(y)
        if row and row[0] in self.kinds and self.style(row[0]).mode != "hidden" and \
                abs(y - (row[1] + row[2] + GAP / 2)) <= 3:
            return row
        return None

    def set_style(self, name, **changes):
        st = self.style(name)
        if "mode" in changes and changes["mode"] != st.mode:
            st.height = None  # each mode has its own natural height
        for k, v in changes.items():
            setattr(st, k, v)
        self.panel.refresh()

    def track_menu(self, name, global_pos):
        kind = self.kinds[name]
        st = self.style(name)
        menu = QMenu(self)
        menu.addSection(name)
        for mode in MODES:
            a = menu.addAction(mode.capitalize(), lambda m=mode: self.set_style(name, mode=m))
            a.setCheckable(True)
            a.setChecked(st.mode == mode)
        if kind in SCALED:
            scale = menu.addMenu("Scale")
            for key, label in (("linear", "Linear"), ("log2", "Log2"), ("log10", "Log10")):
                a = scale.addAction(label, lambda k=key: self.set_style(name, scale=k))
                a.setCheckable(True)
                a.setChecked(st.scale == key)
        menu.addAction("Reset Height", lambda: self.set_style(name, height=None)).setEnabled(st.height is not None)
        if kind == "variants":
            a = next((a for a in self.editor.alignments if a.name + " variants" == name), None)
            consensus = [v for v in (a.variants or []) if v.level == "consensus"] if a else []
            menu.addSeparator()
            menu.addAction(f"Apply All Consensus Changes to Sequence ({len(consensus)})",
                           lambda: apply_variants(self.editor, a, consensus)).setEnabled(bool(consensus))
            menu.addAction(f"Annotate All Consensus Changes ({len(consensus)})",
                           lambda: annotate_variants(self.editor, a, consensus)).setEnabled(bool(consensus))
        menu.exec(global_pos)

    def variant_menu(self, a, v, global_pos):
        """Right-click on a variant: write it into the sequence, or annotate it."""
        menu = QMenu(self)
        menu.addSection(f"{v.pos + 1:,}: {change_text(v)}")
        apply = menu.addAction("Change the Sequence to Match", lambda: apply_variants(self.editor, a, [v]))
        apply.setEnabled(can_apply(v) and not self.editor.dna_locked)
        if self.editor.dna_locked:
            apply.setText(apply.text() + " (DNA locked)")
        menu.addAction("Annotate This Variant", lambda: annotate_variants(self.editor, a, [v]))
        menu.exec(global_pos)

    # ---- painting ------------------------------------------------------------------
    # ---- painting ------------------------------------------------------------------
    # The tracks are drawn into a cached image (`layer`), redrawn only when what they show
    # changes (view window, size, edits, tracks, settings: `invalidate`). Cursor, selection
    # and the sequence view's frame are drawn on top at every paint, so moving the cursor
    # or scrolling the sequence view costs ~1 ms even with many BAM files open.
    # While dragging to pan, the cached image just slides along, and the tracks are redrawn
    # once the drag pauses (redrawing many BAM files at every mouse move would stutter).
    def invalidate(self):
        self.dirty = True
        self.update()

    def paintEvent(self, event):
        if image_export.ACTIVE:   # an image export: draw at the export's resolution, not into the cached layer
            p = QPainter(self)
            self._paint_layers(p)
            self._paint_pinned(p)
            self._paint_overlay(p)
            p.end()
            return
        dpr = self.devicePixelRatioF()
        key = (self.panel.vstart, self.panel.vend, self.width(), self.height(), dpr)
        old = self.layer_key
        if not self.dirty and self.drag and self.drag[0] == "pan" and old and old[2:] == key[2:] \
                and abs((old[1] - old[0]) - (key[1] - key[0])) < 1e-6:
            dx = (old[0] - key[0]) / (key[1] - key[0]) * self.width_px()
            p = QPainter(self)
            p.fillRect(self.rect(), self.palette().color(QPalette.Base))
            p.drawPixmap(QPointF(dx, 0), self.layer)
            names = QRectF(0, 0, self.left - 4, self.height())  # the names column does not move
            p.drawPixmap(names, self.layer, QRectF(0, 0, (self.left - 4) * dpr, self.height() * dpr))
            self._paint_pinned(p)
            self._paint_overlay(p)
            p.end()
            self.pan_timer.start()
            return
        if self.dirty or key != old:
            self.layer = QPixmap(self.size() * dpr)
            self.layer.setDevicePixelRatio(dpr)
            p = QPainter(self.layer)
            self._paint_layers(p)
            p.end()
            self.layer_key, self.dirty = key, False
        p = QPainter(self)
        p.drawPixmap(0, 0, self.layer)
        self._paint_pinned(p)
        with timed("genome view: cursor + selection (on the cached tracks)"):
            self._paint_overlay(p)
        p.end()

    def _paint_pinned(self, p):
        """Tracks scrolled down: the overview and the ruler again at the top of what shows, over the tracks under them."""
        top = self.view_top()
        if not top:
            return
        p.save()
        p.translate(0, top)
        p.fillRect(QRectF(0, 0, self.width(), self.header_h - 2), self.palette().color(QPalette.Base))
        self._paint_header(p)
        line = QColor(self.palette().color(QPalette.PlaceholderText))
        line.setAlpha(60)
        p.setPen(line)
        p.drawLine(QPointF(0, self.header_h - 2), QPointF(self.width(), self.header_h - 2))
        p.restore()

    def _paint_overlay(self, p):
        ed, panel = self.editor, self.panel
        hl = QColor(theme.selection_color())
        vs, ve = panel.vstart, panel.vend
        _, bottom = self.layout_rows()
        top = self.view_top() + 18   # the ruler's top, wherever the tracks are scrolled to
        if panel.seq_window:  # what the sequence view shows, framed on the ruler
            a, b = panel.seq_window
            p.setPen(QPen(hl, 1, Qt.DashLine))
            p.drawRect(QRectF(self.x(a), top + 1, max(2.0, self.x(b) - self.x(a)), self.ruler_base - 10))
        p.setClipRect(QRectF(self.left - 4, 0, self.width(), self.height()))  # not over the track names
        for s, e in ed.search_spans[:3000]:   # what the search bar found
            if e > vs and s < ve:
                p.fillRect(QRectF(self.x(max(s, vs)), top, max(2.0, self.x(min(e, ve)) - self.x(max(s, vs))), bottom - top), theme.hit_color())
        if ed.end > ed.start:
            sel = QColor(hl)
            sel.setAlpha(60)
            for s, e in ed.selected_spans():  # two parts when it passes the origin
                if e > vs and s < ve:
                    p.fillRect(QRectF(self.x(max(s, vs)), top, self.x(min(e, ve)) - self.x(max(s, vs)) + 1,
                                      bottom - top), sel)
        elif vs <= ed.start <= ve:
            p.setPen(QPen(hl, 1))
            p.drawLine(QPointF(self.x(ed.start), top), QPointF(self.x(ed.start), bottom))

    def _paint_layers(self, p):
        ed, panel = self.editor, self.panel
        pal = self.palette()
        dim = pal.color(QPalette.PlaceholderText)
        p.fillRect(self.rect(), pal.color(QPalette.Base))
        fm = QFontMetricsF(self.font())
        self.hits = []
        rows, _ = self.layout_rows()
        tops = {name: (top, h) for name, top, h in rows if name}
        W = self.width_px()
        vs, ve = panel.vstart, panel.vend
        left = self.left   # the names column's width
        name_h = self.name_h()
        lap = Laps("genome view: tracks")
        self._paint_header(p)
        lap("overview + ruler")

        # GC % (20-80 % range, 50 % midline).
        if row := self.shown(tops, "GC %"):
            top, h = row
            bins = max(1, int(min(W // 2, (ve - vs) // 50)))  # >= 50 bp per bin: no single-base noise
            values = ed.doc.gc_bins(int(vs), int(math.ceil(ve)), bins)
            p.setPen(QPen(dim, 1, Qt.DotLine))
            p.drawLine(QPointF(left, top + h / 2), QPointF(left + W, top + h / 2))
            p.setPen(Qt.NoPen)
            p.setBrush(QColor("#7986cb"))
            for i, v in enumerate(values):
                if v == v:  # not NaN
                    y = top + h * (1 - min(1.0, max(0.0, (v - 0.2) / 0.6)))
                    p.drawRect(QRectF(left + i * W / bins, y, W / bins + 0.5, top + h - y))
        lap("GC %")

        # Bases at deep zoom.
        if "Bases" in tops:
            top, h = tops["Bases"]
            seq = ed.doc.sequence(int(vs), int(math.ceil(ve)))
            cw = W / (ve - vs)
            for i, c in enumerate(seq):
                x = self.x(int(vs) + i)
                p.fillRect(QRectF(x, top, cw - 0.5, h), QColor(BASE_COLORS.get(c.upper(), "#9e9e9e")))
                if cw >= fm.horizontalAdvance("W") + 2:
                    p.setPen(QColor("white"))
                    p.drawText(QRectF(x, top, cw, h), Qt.AlignCenter, c.upper())
            lap("bases")

        # Features.
        if row := self.shown(tops, "Features"):
            top, h = row
            lh = self.lane_h("Features")
            labels = lh >= 11 * self.scale()  # condensed: thin lanes, no labels
            feats = [f for f in ed.shown if any(s < ve and e > vs for s, e in f.spans)]
            if len(feats) > DENSE:
                self._density(p, [(s, e) for f in feats for s, e in f.spans], top, h, QColor("#90a4ae"))
            else:
                boxes = [(self.x(min(s for s, _ in f.spans)), self.x(max(e for _, e in f.spans)), f) for f in feats]
                for lane, x0, x1, f in pack(boxes, self.lanes("Features", h)):
                    y = top + lane * lh
                    color = ed.colors[f.index]
                    for s, e in f.spans:
                        if s < ve and e > vs:
                            r = QRectF(max(left, self.x(s)), y,
                                       max(1.5, min(left + W, self.x(e)) - max(left, self.x(s))), lh - 2)
                            p.fillRect(r, color)
                            self.hits.append((r, ("feature", f.index)))
                    r = QRectF(max(left, x0), y, min(left + W, x1) - max(left, x0), lh - 2)
                    if labels and r.width() > fm.horizontalAdvance(f.label) + 6:
                        p.setPen(text_color_on(color))
                        way = arrow(f)
                        p.drawText(r, Qt.AlignCenter, ("‹ " if way < 0 else "") + f.label + (" ›" if way > 0 else ""))
        lap("features")

        # Primer sites.
        if "Primers" in tops and (row := self.shown(tops, "Primers")):
            top, h = row
            lh = self.lane_h("Primers")
            boxes = [(self.x(s.extent[0][0]), self.x(s.extent[-1][1]), s) for s in ed.drawn_primers
                     if any(a < ve and b > vs for a, b in s.spans)]
            for lane, x0, x1, site in pack(boxes, self.lanes("Primers", h)):
                y = top + lane * lh + lh / 2
                p.setPen(QPen(QColor(primer_color(site)), 2))
                for a, b in site.spans:
                    p.drawLine(QPointF(self.x(a), y), QPointF(self.x(b), y))
                    self.hits.append((QRectF(self.x(a), y - 4, max(4.0, self.x(b) - self.x(a)), 8), ("primer", site.key)))
        lap("primers")

        # Loaded tracks.
        for i, t in enumerate(ed.tracks):
            if not (row := self.shown(tops, t.name)):
                continue
            top, h = row
            lh = self.lane_h(t.name)
            color = QColor(TRACK_COLORS[i % len(TRACK_COLORS)])
            items, more = t.track.query(t.seqid, int(vs), int(math.ceil(ve)), DENSE)
            if more:
                counts = t.track.density(t.seqid, int(vs), int(math.ceil(ve)), max(1, int(W // 2)))
                self._bars(p, counts, top, h, color, self.style(t.name).scale)
                continue
            boxes = [(self.x(s), self.x(e), (s, e, strand, name)) for s, e, strand, name, _ in items]
            for lane, x0, x1, (s, e, strand, name) in pack(boxes, self.lanes(t.name, h)):
                r = QRectF(max(left, x0), top + lane * lh, max(1.5, min(left + W, x1) - max(left, x0)), lh - 2)
                p.fillRect(r, color)
                self.hits.append((r, ("range", s, e)))
                if lh >= 11 * self.scale() and r.width() > fm.horizontalAdvance(name) + 6:
                    p.setPen(text_color_on(color))
                    p.drawText(r, Qt.AlignCenter, name)
        lap("annotation tracks")

        if row := self.shown(tops, "Query hits") if ed.query_hits else None:
            self._draw_hits(p, row, pal.color(QPalette.Base), dim, fm)
            lap("query hits")

        # BAM alignments: coverage, then reads.
        for a in ed.alignments:
            if row := self.shown(tops, a.name):
                self._coverage(p, a, row, dim)
            lap("BAM coverage")
            if row := self.shown(tops, a.name + " variants"):
                self._variants(p, a, row, dim, fm)
            lap("BAM variants")
            if not (row := self.shown(tops, a.name + " reads")):
                continue
            top, h = row
            if a.error:
                p.setPen(dim)
                p.drawText(QRectF(left + 4, top, W, name_h), Qt.AlignLeft, f"Cannot read: {a.error}")
            elif ve - vs > setting("reads_span"):
                p.setPen(dim)
                p.drawText(QRectF(left + 4, top, W, name_h), Qt.AlignLeft,
                           f"Zoom in under {format_bp(setting('reads_span'), 1000)} to see reads")
            elif self._draw_reads(p, a, row, pal.color(QPalette.Base), dim) and (cov := self.shown(tops, a.name)):
                p.setPen(dim)  # at the right end of the coverage row (the reads rows are full)
                p.drawText(QRectF(left, cov[0], W - 4, name_h), Qt.AlignRight,
                           "more reads than shown: drag the line under the reads down")
            lap("BAM reads")

        # Track names (▾ shown / ▸ hidden: click to toggle) and a thin line under each track.
        p.fillRect(QRectF(0, 0, left - 4, self.height()), pal.color(QPalette.Base))
        line = QColor(dim)
        line.setAlpha(60)
        for name, top, h in rows:
            if not name:
                continue
            p.setPen(line)
            y = top + h + GAP / 2
            p.drawLine(QPointF(0, y), QPointF(self.width(), y))
            p.setPen(dim)
            toggle = "" if name not in self.kinds else ("▸ " if self.style(name).mode == "hidden" else "▾ ")
            p.drawText(QRectF(2, top, left - 8, min(h, name_h)), Qt.AlignLeft | Qt.AlignVCenter,
                       fm.elidedText(toggle + self.track_label(name), Qt.ElideRight, left - 8))
            if self.kinds.get(name) == "hits" and self.style(name).mode != "hidden":
                self._hit_labels(p, top, h, fm)
        lap("track names")
        lap.done()

    def _paint_header(self, p):
        """The overview (the whole sequence, the window in view) and the ruler, at the top."""
        ed, panel = self.editor, self.panel
        dim = self.palette().color(QPalette.PlaceholderText)
        left, W, n = self.left, self.width_px(), max(1, ed.length)
        vs, ve = panel.vstart, panel.vend
        p.fillRect(QRectF(left, 6, W, 4), dim)
        hl = QColor(theme.selection_color())
        p.setPen(QPen(hl, 2))
        p.setBrush(Qt.NoBrush)
        p.drawRect(QRectF(left + vs / n * W, 3, max(3.0, (ve - vs) / n * W), 10))
        # Ruler: labels on their baseline, the ticks and the line under them.
        base = self.ruler_base
        step = nice_step((ve - vs) / (W / 110))
        p.setPen(dim)
        p.drawLine(QPointF(left, base + 8), QPointF(left + W, base + 8))
        tick = math.ceil(vs / step) * step
        while tick <= ve:
            x = self.x(tick)
            p.drawLine(QPointF(x, base + 2), QPointF(x, base + 8))
            p.drawText(QPointF(x + 2, base), format_bp(int(tick), step))
            tick += step

    def _hit_labels(self, p, top, height, fm):
        """One label per alignment row, in the names column (it stays put while the view pans): number, strand, identity. The
        inverted ones are in the purple of their bars; the picked one is bold."""
        lh = self.lane_h("Query hits")
        dim = self.palette().color(QPalette.PlaceholderText)
        font = QFont(self.font())
        for i, hit in enumerate(self.editor.query_hits):
            y = top + self.hits_head() + i * lh
            if y + lh > top + height + 1:
                break
            font.setBold(i == self.editor.current_query_hit)
            p.setFont(font)
            p.setPen(QColor(INVERTED_COLOR) if hit.strand < 0 else dim)
            p.drawText(QRectF(2, y, self.left - 8, lh), Qt.AlignLeft | Qt.AlignVCenter,
                       QFontMetricsF(font).elidedText(self.hit_label(i + 1), Qt.ElideRight, self.left - 8))
        p.setFont(self.font())

    def _density(self, p, spans, top, h, color):
        bins = max(1, int(self.width_px() // 2))
        counts = [0] * bins
        vs, ve = self.panel.vstart, self.panel.vend
        width = (ve - vs) / bins
        for s, e in spans:
            a, b = max(0, int((s - vs) / width)), min(bins, int(math.ceil((e - vs) / width)))
            for i in range(a, max(a + 1, b)):
                if i < bins:
                    counts[i] += 1
        self._bars(p, counts, top, h, color)

    def _coverage(self, p, a, row, dim):
        """Depth bars. Narrow windows: counted per base from the reads, and where enough reads
        disagree with the reference the bar is split by base (IGV-style); wide windows: per-
        block totals that need one background pass over the file first. Consensus changes
        found by the variant scan are flagged with a red mark above the bar at any zoom."""
        top, h = row
        vs, ve = self.panel.vstart, self.panel.vend
        W, span = self.width_px(), ve - vs
        scale = self.style(a.name).scale
        left = self.left
        p.setPen(dim)
        if span <= setting("reads_span"):
            try:
                peak = self._pileup_bars(p, a, top, h, scale)
            except OSError as e:
                a.error = str(e)
                return
        else:
            if not a.bam.has_depth(a.seqid):
                self._start_depth(a)
                p.drawText(QRectF(left + 4, top, W, self.name_h()), Qt.AlignLeft, a.error or "Computing coverage…")
                return
            # Bins of >= 512 bp, so Rust answers from its block totals, not from reads.
            bins = max(1, min(int(W // 2), int(span // 512)))
            try:
                values = a.bam.coverage(a.seqid, int(vs), int(math.ceil(ve)), bins)
            except OSError as e:
                a.error = str(e)
                return
            self._bars(p, values, top, h, QColor("#9e9e9e"), scale)
            peak = max(values, default=0)
        p.setPen(dim)
        baseline = top + max(10, math.ceil(QFontMetricsF(self.font()).ascent()) - 2)
        p.drawText(QPointF(left + 2, baseline), f"{peak:,.0f}×" + ("" if scale == "linear" else f" ({scale})"))
        self._flags(p, a, top)

    def _pileup_bars(self, p, a, top, h, scale="linear"):
        """Per-base coverage of the view; returns the peak depth."""
        f = SCALES[scale]
        ed, vs, ve = self.editor, int(self.panel.vstart), int(math.ceil(self.panel.vend))
        with timed("genome view: BAM pileup (Rust)"):
            counts = a.bam.pileup(a.seqid, vs, ve, setting("min_baseq"), setting("min_mapq"))
        ref = ed.doc.sequence(vs, min(ve, ed.length)).upper()
        depths = [sum(c[:5]) for c in counts]  # bases (deleted ones not drawn, like IGV)
        W = self.width_px()
        bins = max(1, min(int(W // 2), len(depths)))
        per = len(depths) / bins if depths else 1
        means = [sum(depths[int(i * per):max(int(i * per) + 1, int((i + 1) * per))]) /
                 max(1, int((i + 1) * per) - int(i * per)) for i in range(bins)] if depths else []
        peak = f(max(means, default=0)) or 1
        self._bars(p, means, top, h, QColor("#9e9e9e"), scale)
        # Positions where reads disagree: the bar split by base, height to scale.
        cw = W / max(1, ve - vs)
        for i, c in enumerate(counts):
            want = ref[i] if i < len(ref) else "N"
            slot = "ACGT".find(want)
            depth = depths[i]
            differ = depth - (c[slot] if slot >= 0 else 0)
            if slot < 0 or not flag_level(differ, depth):
                continue
            x, y = self.x(vs + i), top + h
            total = h * f(depth) / peak  # the bar's height on the chosen scale, split by base share
            for b, n in zip("ACGT", c[:4]):
                if n:
                    bh = total * n / depth
                    y -= bh
                    p.fillRect(QRectF(x, y, max(2.0, cw), bh), QColor(BASE_COLORS[b]))
        return max(depths, default=0)

    def _flags(self, p, a, top):
        """Red marks for consensus changes from the variant scan (started on demand)."""
        if a.variants is None:
            start_scan(a, self.editor)
            return
        vs, ve = self.panel.vstart, self.panel.vend
        p.setPen(Qt.NoPen)
        p.setBrush(QColor(FLAG_COLOR))
        for v in a.variants:
            if v.level == "consensus" and vs <= v.pos <= ve:
                x = self.x(v.pos + 0.5)
                p.drawPolygon(QPolygonF([QPointF(x - 4, top - 4), QPointF(x + 4, top - 4), QPointF(x, top + 2)]))
                self.hits.append((QRectF(x - 5, top - 5, 10, 8), ("read", v.pos, v.pos + max(1, len(v.ref)), v.text)))

    def _variants(self, p, a, row, dim, fm):
        """Every position where the reads differ (variant scan): solid = consensus change, faded = minority variant.
        SNPs in the colour of the sample's base, deletions dark grey, insertions purple. Right-click one to apply or
        annotate it."""
        top, h = row
        if a.variants is None:
            start_scan(a, self.editor)
            p.setPen(dim)
            p.drawText(QRectF(self.left + 4, top, self.width_px(), self.name_h()), Qt.AlignLeft, "Scanning the reads…")
            return
        vs, ve = self.panel.vstart, self.panel.vend
        seen = [v for v in a.variants if v.pos + max(1, len(v.ref)) > vs and v.pos < ve]
        if len(seen) > DENSE:
            self._density(p, [(v.pos, v.pos + max(1, len(v.ref))) for v in seen], top, h, QColor(FLAG_COLOR))
            return
        lh = self.lane_h(a.name + " variants")
        boxes = [(self.x(v.pos), self.x(v.pos + max(1, len(v.ref))), v) for v in seen]
        for lane, x0, x1, v in pack(boxes, self.lanes(a.name + " variants", h)):
            color = QColor(BASE_COLORS.get(v.alt[:1], "#757575") if v.kind == "snp" else
                           "#424242" if v.kind == "del" else INSERTION_COLOR)
            if v.level != "consensus":
                color.setAlpha(110)
            left = self.left
            r = QRectF(max(left, x0 - 1), top + lane * lh, max(3.0, min(left + self.width_px(), x1) - max(left, x0) + 1), lh - 2)
            p.fillRect(r, color)
            self.hits.append((r, ("variant", a, v)))
            label = change_text(v)
            if lh >= 11 * self.scale() and r.width() > fm.horizontalAdvance(label) + 6:
                p.setPen(text_color_on(color))
                p.drawText(r, Qt.AlignCenter, label)

    def _start_depth(self, a):
        if a.depth_pending:
            return
        a.depth_pending = True
        bam, seqid, length = a.bam, a.seqid, a.length

        def work():  # runs in a pool thread; Rust releases Python's lock meanwhile
            try:
                bam.coverage(seqid, 0, length, 1)  # one wide bin = compute the block totals
            except OSError as e:
                a.error = str(e)
            _relay.depth_ready.emit()

        QThreadPool.globalInstance().start(work)

    def _reads(self, a, max_lanes=READ_LANES):
        """Reads around the view, packed into rows. Cached for a window 1.5x the view on each
        side, so panning reuses it; rows are packed in bp so they do not jump while panning."""
        ed, vs, ve = self.editor, self.panel.vstart, self.panel.vend
        c = a.cache
        if c and c[0] <= vs and ve <= c[1] and c[2] == ed.version and c[6] == max_lanes:
            return c[3], c[4]
        span = ve - vs
        s, e = max(0, int(vs - 1.5 * span)), min(ed.length, int(ve + 1.5 * span) + 1)
        with timed("genome view: BAM reads fetch (Rust)"):
            reads, more = a.bam.reads(a.seqid, s, e, ed.doc.sequence(s, e))
        limits = insert_limits(reads)  # from all reads in the window, drawn or not
        ends, laid = [], []
        for r in reads:  # sorted by start
            lane = next((i for i, end in enumerate(ends) if end < r[0]), None)
            if lane is None:
                if len(ends) == max_lanes:
                    more = True  # too deep to draw
                    continue
                ends.append(0)
                lane = len(ends) - 1
            ends[lane] = r[1]
            laid.append((lane, r))
        a.cache = (s, e, ed.version, laid, more, limits, max_lanes)
        return laid, more

    def _draw_reads(self, p, a, row, background, dim):
        top, height = row
        name = a.name + " reads"
        lane_h = self.lane_h(name)
        try:
            laid, more = self._reads(a, self.lanes(name, height))
        except OSError as e:
            a.error = str(e)
            return
        limits, mode = a.cache[5], setting("color_by")
        vs, ve = self.panel.vstart, self.panel.vend
        cw = self.width_px() / (ve - vs)
        h = max(1, lane_h - 2) if lane_h > 3 else lane_h - 1  # condensed: 2-px reads, 1-px gap
        colors = {k: QColor(v) for k, v in READ_COLORS.items()}
        faint = {k: QColor(c.red(), c.green(), c.blue(), 70) for k, c in colors.items()}
        base_colors = {b: QColor(c) for b, c in BASE_COLORS.items()}
        unmapped_mate = QPen(QColor(MATE_UNMAPPED_COLOR), 1)
        insertion = QColor(INSERTION_COLOR)
        x0 = self.left - vs * cw  # x of position p = x0 + p * cw (self.x inlined: called per mismatch)
        mark_w = max(1.0, cw)
        self.read_kinds = {}  # category -> count drawn (for the info line and tests)
        for lane, read in laid:
            s, e, strand, mapq, name, mism, ins, dels, flags, _, _, insert = read
            if e <= vs or s >= ve:
                continue
            kind = classify(read, limits, mode)
            self.read_kinds[kind] = self.read_kinds.get(kind, 0) + 1
            y = top + lane * lane_h
            r = QRectF(x0 + s * cw, y, max(1.0, (e - s) * cw), h)
            p.fillRect(r, (faint if mapq == 0 else colors)[kind])  # MAPQ 0 = could map elsewhere
            if flags & 0x1 and flags & 0x8:  # paired, mate unmapped
                p.setPen(unmapped_mate)
                p.setBrush(Qt.NoBrush)
                p.drawRect(r)
            p.setPen(dim)
            for ds, de in dels:
                if de > vs and ds < ve:
                    a, b = x0 + ds * cw, x0 + de * cw
                    p.fillRect(QRectF(a, y, b - a, h), background)
                    p.drawLine(QPointF(a, y + h / 2), QPointF(b, y + h / 2))
            # Only marks in view, and one per pixel column when zoomed out (long noisy reads
            # have thousands; more in one pixel would not show anyway).
            last = None
            for pos, b in mism:
                if vs <= pos < ve:
                    x = x0 + pos * cw
                    if int(x) != last:
                        last = int(x)
                        p.fillRect(QRectF(x, y, mark_w, h), base_colors.get(b, dim))
            last = None
            for pos, _ in ins:
                if vs <= pos < ve:
                    x = x0 + pos * cw
                    if int(x) != last:
                        last = int(x)
                        p.fillRect(QRectF(x - 1, y - 1, 2, h + 2), insertion)
            info = f"{name}  {s + 1:,}–{e:,} ({'+' if strand > 0 else '−'})  MAPQ {mapq}"
            if mism or ins or dels:
                info += f"  {len(mism)} mismatch, {len(ins)} ins, {len(dels)} del"
            if insert:
                info += f"  insert {abs(insert):,} bp"
            if LABELS.get(kind):
                info += f"  ({LABELS[kind]})"
            if flags & 0x1 and flags & 0x8:
                info += "  (mate unmapped)"
            self.hits.append((r, ("read", s, e, info)))
        return more

    def hit_label(self, rank):
        """The text of alignment row `rank` (1-based) in the names column: its number, strand and identity."""
        hit = self.editor.query_hits[rank - 1]
        return f"#{rank} {'↺' if hit.strand < 0 else '→'} {hit.identity:.1%}"

    def hit_row_at(self, y):
        """Index of the Local Align row under y (a row of the "Query hits" track), or -1."""
        row = self.row_at(y)
        if not row or self.kinds.get(row[0]) != "hits" or self.style(row[0]).mode == "hidden":
            return -1
        head = self.hits_head()
        i = int((y - row[1] - head) // self.lane_h(row[0])) if y >= row[1] + head else -1
        return i if 0 <= i < len(self.editor.query_hits) else -1

    def _draw_hits(self, p, row, background, dim, fm):
        """Local Align hits, drawn like reads, each alignment on its own row (so overlapping or repeated hits never sit on top
        of each other): a bar over the template span, coloured mismatches, insertion marks, deletion gaps. A hit across the
        origin is drawn twice. The picked row is tinted; its label is in the names column (see _hit_labels)."""
        top, height = row
        bottom = top + height
        top += self.hits_head()
        lh = self.lane_h("Query hits")
        left = self.left
        vs, ve, n = self.panel.vstart, self.panel.vend, self.editor.length
        cw = self.width_px() / (ve - vs)
        h = max(1, lh - 2)
        base_colors = {b: QColor(c) for b, c in BASE_COLORS.items()}
        for rank, hit in enumerate(self.editor.query_hits, 1):
            y0 = top + (rank - 1) * lh
            if y0 + lh > bottom + 1:   # the track was made shorter than its rows
                break
            if rank - 1 == self.editor.current_query_hit:
                tint = QColor(theme.selection_color())
                tint.setAlpha(45)
                p.fillRect(QRectF(left, y0 - 1, self.width_px() + 8, lh), tint)
            inverted = hit.strand < 0
            color = QColor(INVERTED_COLOR if inverted else COLORS["forward"])
            info = (f"Hit {rank}  {hit.start + 1:,}–{(hit.end - 1) % n + 1:,} "
                    f"({'↺ inverted, reverse complement' if inverted else 'forward'})  "
                    f"{hit.identity:.1%} identity, score {hit.score}  {len(hit.mismatches)} mismatch, "
                    f"{len(hit.insertions)} ins, {len(hit.deletions)} del")
            for shift in (0, -n) if hit.end > n else (0,):
                s, e = hit.start + shift, hit.end + shift
                if e <= vs or s >= ve:
                    continue
                x0 = left - vs * cw + shift * cw  # x of template position q = x0 + q * cw
                r = QRectF(left - vs * cw + s * cw, y0, max(1.0, (e - s) * cw), h)
                p.fillRect(r, color)
                if inverted:  # a light hatch marks the inverted pieces even without colour vision
                    p.fillRect(r, QBrush(QColor(255, 255, 255, 80), Qt.BDiagPattern))
                for ds, de in hit.deletions:
                    a, b = x0 + ds * cw, x0 + de * cw
                    p.fillRect(QRectF(a, y0, b - a, h), background)
                    p.setPen(dim)
                    p.drawLine(QPointF(a, y0 + h / 2), QPointF(b, y0 + h / 2))
                for pos, b in hit.mismatches:
                    if vs <= pos + shift < ve:
                        p.fillRect(QRectF(x0 + pos * cw, y0, max(1.0, cw), h), base_colors.get(b, dim))
                for pos, _ in hit.insertions:
                    if vs <= pos + shift < ve:
                        p.fillRect(QRectF(x0 + pos * cw - 1, y0 - 1, 2, h + 2), QColor(INSERTION_COLOR))
                label = f"{'↺ ' if inverted else ''}#{rank}  {hit.identity:.1%}"
                if r.width() > fm.horizontalAdvance(label) + 8:
                    p.setPen(QColor("black"))
                    p.drawText(r, Qt.AlignCenter, label)
                self.hits.append((r, ("qhit", rank - 1, info)))

    def _bars(self, p, counts, top, h, color, scale="linear"):
        f = SCALES[scale]
        peak = f(max(counts, default=0)) or 1
        W = self.width_px()
        p.setPen(Qt.NoPen)
        p.setBrush(color)
        for i, c in enumerate(counts):
            if c:
                bh = max(1.0, h * f(c) / peak)
                p.drawRect(QRectF(self.left + i * W / len(counts), top + h - bh, W / len(counts) + 0.5, bh))

    # ---- mouse -----------------------------------------------------------------------
    def name_at(self, x, y):
        """Track whose name (the line that hides / shows it) is under (x, y), else None."""
        row = self.row_at(y)
        if x < self.left - 4 and row and row[0] in self.kinds and y < row[1] + self.name_h():
            return row[0]
        return None

    def hover_target(self, pt):
        """What a tooltip is for under a point: a feature (its index), or a track's name (its whole name and what a click
        does: the column cuts long names short)."""
        x, y = pt.x(), pt.y()
        if self.in_header(y):
            return None
        name = self.name_at(x, y)
        if name is not None:
            hidden = self.style(name).mode == "hidden"
            return Text(f"<b>{escape(name)}</b><br>Click: {'show' if hidden else 'hide'} this track. "
                        "Right-click: expanded, condensed or hidden" + (", scale" if self.kinds[name] in SCALED else "") + ".")
        return next((a[1] for r, a in reversed(self.hits) if a[0] == "feature" and r.contains(QPointF(pt))), None)

    def mousePressEvent(self, event):
        x, y = event.position().x(), event.position().y()
        row = self.row_at(y)
        header = self.in_header(y)
        if event.button() == Qt.RightButton:
            pt = event.position()
            hit = next((a for r, a in reversed(self.hits) if a[0] in ("feature", "primer", "variant")
                        and r.adjusted(-1, -2, 1, 2).contains(pt)), None) if not header else None
            if hit and hit[0] == "variant":
                self.variant_menu(hit[1], hit[2], event.globalPosition().toPoint())
            elif hit:  # a feature or a primer: its menu; elsewhere in a track: the track's menu
                from feature_table import feature_menu, primer_menu  # here: feature_table imports the tools
                (primer_menu if hit[0] == "primer" else feature_menu)(self, self.editor, hit[1], event.globalPosition().toPoint())
            elif row and row[0] in self.kinds:
                self.track_menu(row[0], event.globalPosition().toPoint())
            return
        if event.button() != Qt.LeftButton:
            return
        if self.on_overview(x, y):  # overview bar: jump there
            self.panel.center_on(int((x - self.left) / self.width_px() * self.editor.length))
            return
        if edge := self.edge_at(y):
            self.resizing = (edge[0], y, edge[2])
            return
        if x < self.left - 4 and (i := self.hit_row_at(y)) >= 0:   # an alignment's label: pick its row
            self.editor.select_query_hit(i)
            self.editor.navigate.emit(self.editor.start)
            return
        if name := self.name_at(x, y):  # the name: hide / show
            self.set_style(name, mode="expanded" if self.style(name).mode == "hidden" else "hidden")
            return
        mode = "select" if event.modifiers() & Qt.ShiftModifier else "pan"
        self.drag = (mode, x, self.panel.vstart, self.panel.vend, False)

    def on_overview(self, x, y):
        """Is (x, y) on the overview bar at the top (a click there jumps)?"""
        return x >= self.left and y - self.view_top() < 16

    def mouseMoveEvent(self, event):
        y = event.position().y()
        if self.resizing:
            name, y0, h0 = self.resizing
            bar = KINDS[self.kinds[name]][0] == 0
            self.style(name).height = max(10 if bar else self.lane_h(name), int(h0 + y - y0))  # >= one lane
            self.panel.refresh()
            return
        if not self.drag:
            pt = event.position()
            x = pt.x()
            # a hand over what a click does something with: a feature, a primer, a read..., a track's name, the overview
            clickable = self.on_overview(x, y) or (not self.in_header(y) and (
                any(r.adjusted(-1, -2, 1, 2).contains(pt) for r, _ in self.hits) or self.name_at(x, y) is not None
                or (x < self.left - 4 and self.hit_row_at(y) >= 0)))
            self.setCursor(Qt.SplitVCursor if self.edge_at(y) else Qt.PointingHandCursor if clickable else Qt.ArrowCursor)
            return
        mode, x0, vs, ve, _ = self.drag
        x = event.position().x()
        self.drag = (mode, x0, vs, ve, abs(x - x0) > 3)
        if not self.drag[4]:
            return
        if mode == "pan":
            self.setCursor(Qt.ClosedHandCursor)
            shift = (x0 - x) / self.width_px() * (ve - vs)
            self.panel.set_view(vs + shift, ve + shift)
        else:
            self.editor.select(int(self.pos(x0)), int(self.pos(x)))

    def mouseReleaseEvent(self, event):
        if self.resizing:
            self.resizing = None
            return
        if not self.drag:
            return
        mode, x0, *_, moved = self.drag
        self.drag = None
        self.setCursor(Qt.ArrowCursor)
        if moved:
            if mode == "select":
                self.editor.reveal.emit(self.editor.start)
                self.editor.navigate.emit(self.editor.start)
            else:
                self.invalidate()  # end of a pan: draw the tracks for where it stopped
            return
        pt = event.position()
        for rect, action in reversed(self.hits) if not self.in_header(pt.y()) else ():   # (tracks under the ruler: hidden)
            if rect.adjusted(-1, -2, 1, 2).contains(pt):
                if action[0] == "feature":
                    self.editor.select_feature(action[1])
                elif action[0] == "primer":
                    self.editor.select_primer(action[1], extend=bool(event.modifiers() & Qt.ShiftModifier))
                elif action[0] == "qhit":   # an alignment row: pick it (the dialog's table follows)
                    self.editor.select_query_hit(action[1])
                    self.panel.info.setText(action[2])
                elif action[0] == "read":
                    self.editor.select(action[1], action[2])
                    self.panel.info.setText(action[3])
                elif action[0] == "variant":
                    v = action[2]
                    self.editor.select(v.pos, min(self.editor.length, v.pos + max(1, len(v.ref))))
                    self.panel.info.setText(v.text + "  (right-click: apply or annotate)")
                else:
                    self.editor.select(action[1], action[2])
                    self.editor.reveal.emit(action[1])
                self.editor.navigate.emit(self.editor.start)   # whatever was picked, the other views go to it
                return
        pos = int(round(self.pos(pt.x())))
        self.editor.set_cursor(pos)
        self.editor.reveal.emit(pos)
        self.editor.navigate.emit(pos)

    def mouseDoubleClickEvent(self, event):
        self.panel.zoom(0.25, self.pos(event.position().x()))

    def wheelEvent(self, event):
        delta = event.angleDelta()
        if event.position().x() < self.left - 4 and not event.modifiers() & Qt.ShiftModifier:
            event.ignore()   # over the track names: the scroll area scrolls the tracks up and down
            return
        vs, ve = self.panel.vstart, self.panel.vend
        if event.modifiers() & Qt.ShiftModifier or abs(delta.x()) > abs(delta.y()):   # pan: a tenth of the view per notch
            step = (delta.x() or delta.y()) / 120 * (ve - vs) * 0.1
            self.panel.set_view(vs - step, ve - step)
        elif delta.y():   # zoom, in proportion to the delta (a touchpad's small steps zoom a little)
            self.panel.zoom(0.8 ** (delta.y() / 120), self.pos(event.position().x()))
        event.accept()


def legend_html():
    """What each colour and mark means (the "?" of the header), with the current thresholds."""
    def swatch(color, text, cell=None):
        cell = cell or f'<td bgcolor="{color}" width="22">&nbsp;</td>'
        return f"<tr>{cell}<td>&nbsp;{text}</td></tr>"

    def mark(char, color, text):
        return swatch("", text, f'<td width="22" align="center"><b><font color="{color}">{char}</font></b></td>')

    def section(title, *items):  # one table each: a title across both columns would widen the swatches
        return f"<p><b>{title}</b></p><table cellspacing='2' cellpadding='0'>{''.join(items)}</table>"

    faint = QColor(COLORS["larger"]).lighter(160).name()
    reads = [swatch(COLORS[k], LABELS[k] or "normal pair, or unpaired read")
             for k in ("normal", "larger", "smaller", "LL", "RR", "RL", "elsewhere")]
    return section(
        "Reads, coloured by pair orientation &amp; insert size", *reads,
        mark("□", MATE_UNMAPPED_COLOR, "purple outline: mate unmapped"),
        swatch(faint, "faded colour: mapping quality 0 (could map elsewhere)"),
    ) + section(
        "Reads, coloured by strand (Settings → Alignments)",
        swatch(COLORS["forward"], "forward"), swatch(COLORS["reverse"], "reverse"),
    ) + section(
        "Inside reads",
        *[swatch(c, f"mismatch: {b}") for b, c in BASE_COLORS.items()],
        swatch(INSERTION_COLOR, "insertion (thin bar)"),
        mark("—", "#757575", "deletion (gap with a line)"),
    ) + section(
        "Coverage",
        swatch("#9e9e9e", "depth (peak written at the left; right-click for a log scale)"),
        swatch(BASE_COLORS["G"], f"split by base: over {setting('variant_pct')} % of reads differ"),
        mark("▼", FLAG_COLOR, f"consensus change: over {setting('consensus_pct')} % of reads, "
                              f"depth ≥ {setting('min_depth')}"),
    ) + section(
        "Variants track (under each BAM's coverage)",
        *[swatch(c, f"substitution to {b}") for b, c in BASE_COLORS.items()],
        swatch("#424242", "deletion"), swatch(INSERTION_COLOR, "insertion"),
        swatch("#bdbdbd", "faded: a minority of the reads (solid: consensus change)"),
        mark("☰", "#757575", "right-click a variant: change the sequence to match, or annotate it"),
    ) + section(
        "Other tracks",
        swatch("#7986cb", "GC % (20–80 %, dotted line = 50 %)"),
        swatch("#5e35b1", "primer site (in this file)"), swatch("#ef6c00", "primer site with mismatches"),
        swatch(TEMP_COLOR, "primer not in this file: designed, or only in the primer database"),
    ) + section(
        "Tracks",
        mark("▾", "#757575", "click a track name to hide / show it"),
        mark("☰", "#757575", "right-click a track: expanded / condensed / hidden, scale"),
        mark("↕", "#757575", "drag the thin line under a track to resize it"),
    )


class InfoLabel(QLabel):
    """The header's info line. Its text changes at every pan and zoom, and a click on a read puts a long line in it: it takes
    the room the header leaves it and never widens the pane (the text is cut short with "…", the whole of it in its tooltip)."""

    def __init__(self):
        super().__init__()
        self.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Preferred)
        self.setMinimumWidth(0)

    def setText(self, text):
        super().setText(text)
        self._tip()

    def resizeEvent(self, event):
        super().resizeEvent(event)
        self._tip()

    def _tip(self):
        cut = self.fontMetrics().horizontalAdvance(self.text()) > self.contentsRect().width()
        self.setToolTip(self.text() if cut else "")

    def paintEvent(self, event):
        p = QPainter(self)
        r = self.contentsRect()
        p.setPen(self.palette().color(QPalette.WindowText))
        p.drawText(r, Qt.AlignLeft | Qt.AlignVCenter, self.fontMetrics().elidedText(self.text(), Qt.ElideRight, r.width()))
        p.end()


class LegendWindow(QScrollArea):
    """The legend in its own small window (stays on top of the app), until closed or its
    button is clicked again."""
    closed = Signal()

    def __init__(self, parent):
        super().__init__(parent)
        self.setWindowFlags(Qt.Tool)  # its own small window, kept above the app
        self.setWindowTitle("Alignment view legend")
        self.label = QLabel(legend_html(), margin=10)
        self.label.setTextFormat(Qt.RichText)
        self.setWidget(self.label)
        self.setWidgetResizable(True)
        self.resize(self.label.sizeHint().width() + 30, min(720, self.label.sizeHint().height() + 20))

    def refresh(self):
        self.label.setText(legend_html())  # the thresholds are in it

    def closeEvent(self, event):
        self.closed.emit()
        super().closeEvent(event)


class GenomeView(QWidget):
    """Header (locus box, zoom buttons, legend) + canvas."""

    def __init__(self, editor):
        super().__init__()
        self.editor = editor
        self.vstart, self.vend = 0.0, float(max(1, editor.length))
        self.seq_window = None
        self.locus = QLineEdit(placeholderText="Go to: 12,000-15,000 · 12000 · gene name", clearButtonEnabled=True)
        self.locus.returnPressed.connect(self.go)
        self.info = InfoLabel()
        header = QHBoxLayout()
        header.setContentsMargins(4, 2, 4, 0)
        header.addWidget(self.locus, 2)
        for text, slot, tip in (("−", lambda: self.zoom(2.0), "Zoom out"), ("+", lambda: self.zoom(0.5), "Zoom in"),
                                ("Whole", self.show_all, "Whole sequence"),
                                ("Selection", self.show_selection, "Zoom to the selection")):
            b = QPushButton(text, clicked=slot, toolTip=tip)
            b.setMaximumWidth(max(80, b.sizeHint().width()))   # (a larger font: the text still fits)
            header.addWidget(b)
        header.addWidget(self.info, 3)
        self.legend = QPushButton("Legend", checkable=True, toolTip="What the colours and marks mean (own window)")
        self.legend.setMaximumWidth(max(80, self.legend.sizeHint().width()))
        self.legend.toggled.connect(self.show_legend)
        self.legend_window = None
        tool_settings.notifier.changed.connect(self._settings_changed)
        header.addWidget(self.legend)
        self.header = header
        self.canvas = GenomeCanvas(self)
        # The tracks scroll up and down when they need more room than the pane has (the overview and the ruler stay at the
        # top: see GenomeCanvas._paint_pinned): a few BAM files no longer push the window past the bottom of the screen.
        self.tracks_area = QScrollArea(widgetResizable=True, frameShape=QFrame.NoFrame)
        self.tracks_area.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.tracks_area.setFocusPolicy(Qt.NoFocus)   # (a click on the canvas would give it the keyboard: see GenomeCanvas)
        self.tracks_area.setWidget(self.canvas)
        self.tracks_area.verticalScrollBar().valueChanged.connect(self.canvas.update)   # the pinned ruler is drawn again
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        layout.addLayout(header)
        layout.addWidget(self.tracks_area, 1)
        # After edits: redraw once typing pauses (GC/feature data may be large to refetch).
        self.timer = QTimer(self, singleShot=True, interval=250, timeout=self._after_edit)
        editor.changed.connect(self.timer.start)
        editor.display_changed.connect(self.refresh)
        editor.primers_changed.connect(self.timer.start)  # re-matched after every edit: wait for a pause
        editor.selection_changed.connect(self.canvas.update)  # cursor / selection only: tracks stay cached
        editor.search_changed.connect(self.canvas.update)
        editor.reveal.connect(self.follow)
        # The pane is given the room its tracks need once the layout has settled (see _fit_splitter), and again only when
        # that changes: a pane made smaller by hand is not grown back at every pan.
        self.fit_timer = QTimer(self, singleShot=True, interval=0, timeout=self._fit)
        self.fitted = None
        self.refresh()

    def show_legend(self, on):
        if on and self.legend_window is None:
            self.legend_window = LegendWindow(self)
            self.legend_window.closed.connect(lambda: self.legend.setChecked(False))
        if self.legend_window is not None:
            self.legend_window.setVisible(on)

    def _settings_changed(self, tool):
        if tool == ALIGNMENTS and self.legend_window is not None:
            self.legend_window.refresh()

    def refresh(self, redraw=True):
        """Layout and header; `redraw`: the tracks' content changed (not just the view
        window, which the canvas notices by itself)."""
        try:
            _, height = self.canvas.layout_rows()
            if self.canvas.minimumHeight() != height + 4:
                self.canvas.setMinimumHeight(height + 4)
            if self.fitted != height + 4:
                self.fit_timer.start()
            span = self.vend - self.vstart
            self.info.setText(f"{int(self.vstart) + 1:,}–{int(self.vend):,}  ({format_bp(int(span), span / 10)})")
            if redraw:
                self.canvas.invalidate()
            else:
                self.canvas.update()
        except RuntimeError:
            pass   # the pane is gone (closed with its window) while the editor still signals

    def export_hidden(self):
        """What an exported image leaves out (image_export): the header row (locus box, buttons, info line) and the scroll
        bar of the tracks."""
        bar = self.tracks_area.verticalScrollBar()
        return [self.header.itemAt(i).widget() for i in range(self.header.count()) if self.header.itemAt(i).widget()] + \
            [bar.parentWidget() if bar.parentWidget() is not self.tracks_area else bar]

    def _fit_splitter(self, canvas_height):
        """A splitter keeps its pane sizes: give this view (or its workspace pane) the room
        its tracks need, taken from the neighbour below (or above), but never more than MAX_SHARE of the
        screen's height: beyond that the tracks scroll."""
        box = self.parentWidget() if self.parentWidget() and not hasattr(self.parentWidget(), "setSizes") else self
        splitter = box.parentWidget()
        if not (hasattr(splitter, "setSizes") and self.isVisible() and splitter.orientation() == Qt.Vertical
                and splitter.count() > 1):
            return
        # Header rows (fixed content) + the canvas. Not sizeHint(): Qt updates it later.
        needed = self.layout().itemAt(0).sizeHint().height() + canvas_height
        if box is not self:
            needed += box.layout().itemAt(0).sizeHint().height()
        screen = self.screen()
        if screen is not None:
            needed = min(needed, int(screen.availableGeometry().height() * MAX_SHARE))
        sizes = splitter.sizes()
        i = splitter.indexOf(box)
        if sizes[i] >= needed:
            return
        j = i + 1 if i + 1 < len(sizes) else i - 1
        give = min(needed - sizes[i], max(0, sizes[j] - 120))  # leave the neighbour usable
        sizes[i] += give
        sizes[j] -= give
        splitter.setSizes(sizes)

    def _fit(self):
        if self.isVisible():
            self.fitted = self.canvas.minimumHeight()
            self._fit_splitter(self.fitted)

    def showEvent(self, event):
        super().showEvent(event)
        self.fit_timer.start()  # tracks may have been added while hidden

    def _after_edit(self):
        if self.vend > self.editor.length:  # sequence got shorter
            self.set_view(self.vstart, self.editor.length)
        self.refresh()

    # ---- view window -------------------------------------------------------------------
    def set_view(self, start, end):
        n = max(1, self.editor.length)
        span = min(max(end - start, MIN_SPAN), n)
        start = min(max(0.0, start), n - span)
        self.vstart, self.vend = start, start + span
        self.refresh(redraw=False)

    def zoom(self, factor, around=None):
        around = (self.vstart + self.vend) / 2 if around is None else around
        self.set_view(around - (around - self.vstart) * factor, around + (self.vend - around) * factor)

    def center_on(self, pos):
        half = (self.vend - self.vstart) / 2
        self.set_view(pos - half, pos + half)

    def show_all(self):
        self.set_view(0, self.editor.length)

    def show_selection(self):
        ed = self.editor
        if ed.end > ed.start:
            margin = (ed.end - ed.start) * 0.1
            self.set_view(ed.start - margin, ed.end + margin)

    def follow(self, pos):
        """Keep `pos` (cursor, found item) on screen."""
        if not self.vstart <= pos <= self.vend:
            self.center_on(pos)

    def set_sequence_window(self, start, end):
        """Region shown by the sequence view: framed on the ruler, followed when the sequence
        view scrolls out of it (not when it is only resized: that would undo a pan here)."""
        scrolled = self.seq_window is None or self.seq_window[0] != start
        self.seq_window = (start, end)
        if scrolled and (end < self.vstart or start > self.vend):
            self.center_on((start + end) / 2)
        else:
            self.canvas.update()

    def go(self):
        """Locus box: "a-b", "a" or a name (features, primers, tracks)."""
        text = self.locus.text().strip()
        # Digits only: a lone "," or "-" holds no number, and a number past the sequence is said so, not passed on to Qt
        numbers = [int(x.replace(",", "")) for x in re.findall(r"\d[\d,]*", text)] if re.fullmatch(r"[\d,\s\-–:]+", text) else []
        ed = self.editor
        if numbers and (min(numbers[:2]) < 1 or max(numbers[:2]) > ed.length):
            ed.error.emit(f"'{text}': positions run from 1 to {ed.length:,}")
            return
        if not numbers and text and re.fullmatch(r"[\d,\s\-–:]+", text):
            ed.error.emit(f"'{text}' is not a position: type 12000, 12,000-15,000 or a name")
            return
        if len(numbers) >= 2:
            a, b = sorted(numbers[:2])
            self.set_view(a - 1, b)
            ed.select(a - 1, b)
            ed.reveal.emit(a - 1)
        elif len(numbers) == 1:
            ed.set_cursor(numbers[0] - 1)
            ed.reveal.emit(numbers[0] - 1)
            self.center_on(numbers[0] - 1)
        elif text:
            hit = ed.lookup(text)
            if hit is None:
                ed.error.emit(f"'{text}' not found")
                return
            start, end, feature = hit
            if feature is not None:
                ed.select_feature(feature)
            else:
                ed.select(start, end)
                ed.reveal.emit(start)
            margin = max(MIN_SPAN, (end - start) * 0.5)
            self.set_view(start - margin, end + margin)
