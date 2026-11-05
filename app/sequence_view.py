"""Sequence view: wrapped rows of top strand, complement, and feature lanes.

Only the rows on screen are fetched from Rust and painted, so a 5 Mb genome scrolls
as smoothly as a plasmid.

Each row is as tall as its content: an enzyme line if a cut site falls in it, the two
strands, an optional 3-frame translation, and only the feature lanes used in that row.
Most rows have the same height (`row_h`); only rows with extra content are stored, with
running sums, so a 250 Mb chromosome (2.5 million rows) costs nothing to lay out.
`offset(row)` gives a row's top y; `row_at(y)` binary-searches back.

Circular sequences are shown twice in a row (a red mark at the origin), so the end runs
on into the start and a selection can be dragged through the origin. Screen rows of the
second copy map back to the same "real" rows (`real_row`); positions on screen ("virtual",
up to 2 x length) are brought back to the sequence by the editor (Editor.select).
"""
import math
from bisect import bisect_left
from itertools import accumulate

from PySide6.QtCore import QEvent, QPoint, QPointF, QRect, QRectF, Qt, QTimer, Signal
from PySide6.QtGui import (
    QColor, QContextMenuEvent, QFont, QFontMetricsF, QGuiApplication, QKeyEvent, QKeySequence, QMouseEvent, QPainter,
    QPainterPath, QPalette, QPen, QPolygonF,
)
from PySide6.QtWidgets import (
    QAbstractScrollArea, QApplication, QMenu, QSlider, QToolButton, QToolTip, QVBoxLayout, QWidget,
)

from cut_sites import CAS9_KEY, feature_layout, runs
import genetic_code
import tool_settings
from editor import arrow, complement, label_font, mono_font, reverse_complement, text_color_on
from hover import HoverTips
from perf import timed
import theme
import units

DEFAULT_SIZE = 11          # letter size, points
MIN_SIZE, MAX_SIZE = 6, 32
ZOOM_BAR_W = 24            # strip kept free at the right edge for the zoom bar
ORIGIN_COLOR = "#d50000"   # the origin mark between the copies of a circular sequence
PRIMER_COLOR = "#5e35b1"   # exact binding site, light theme (the views use theme "primer", which has a dark-theme twin)
ALMOST_COLOR = "#ef6c00"   # site with mismatches (dashed)
TEMP_COLOR = "#00897b"   # a primer that is not part of this file: designed here (not saved yet) or only in the primer database


def primer_color(site):
    """The colour of a primer site, in every view: temporary style (designed here, or only in the database and not imported to
    this file), else orange when it binds with mismatches, else the regular primer colour."""
    return TEMP_COLOR if getattr(site, "provisional", site.temporary) else ALMOST_COLOR if site.mismatches else theme.hex("primer")
STOP_COLOR = "#e53935"     # stop codons in translations
HIT_COLOR = QColor(255, 214, 0, 130)   # what the search bar found, light theme (drawn with theme.hit_color(), which follows the theme)

CODONS = genetic_code.TABLES[1].codons   # the standard code (table 1); other tables: genetic_code.codon_map(id)


def codons(doc, spans, strand, frame=0, table=1, exceptions=None, first_met=False):
    """[(middle base, amino acid, first base, last base + 1)] of the codons of the bases in
    `spans` (in that order; minus strand read backwards), starting `frame` bases in, read with NCBI genetic code `table`.
    `exceptions` {codon number: letter} replaces single codons (/transl_except); `first_met`: a first codon that is a
    start codon of the table is shown as M, as GenBank translates a CDS. Positions are on the sequence."""
    positions = [p for s, e in spans for p in range(s, e)]
    seq = "".join(doc.sequence(s, e) for s, e in spans).upper()
    if strand < 0:
        positions.reverse()
        seq = reverse_complement(seq)
    code = genetic_code.codon_map(table)
    out = [(positions[i + 1], code.get(seq[i:i + 3], "X"), min(positions[i:i + 3]), max(positions[i:i + 3]) + 1)
           for i in range(frame, len(seq) - 2, 3)]
    if first_met and out and seq[frame:frame + 3] in genetic_code.TABLES[table].start_codons:
        out[0] = (out[0][0], "M", out[0][2], out[0][3])
    for k, aa in (exceptions or {}).items():
        if k < len(out):
            out[k] = (out[k][0], aa, out[k][2], out[k][3])
    return out


def cds_codons(doc, f):
    """The codons of CDS feature `f`, read the way its annotation says (/transl_table, /codon_start, /transl_except), or with
    the table chosen in Tools > Settings when it has no /transl_table."""
    plan = genetic_code.cds_plan(doc, f)
    return codons(doc, f.spans, f.strand, plan.frame, plan.table, plan.exceptions, plan.first_met)


class WheelSteps:
    """Whole notches of a wheel. A mouse wheel sends 120 per notch, a touchpad many small deltas, and a sideways scroll 0 on
    the vertical axis: one zoom step per event made a touchpad zoom far too fast, and a sideways Ctrl+scroll zoomed out. Here
    a step is taken for every 120 gathered, in the direction they were gathered."""

    def __init__(self):
        self.rest = 0

    def __call__(self, delta):
        """Steps (+ up / away, - down) for a wheel event's angleDelta().y()."""
        if not delta:
            return 0
        if self.rest and (delta > 0) != (self.rest > 0):   # turned the other way: what was gathered no longer counts
            self.rest = 0
        self.rest += delta
        steps = int(self.rest / 120)
        self.rest -= steps * 120
        return steps


def edit_menu(parent, editor, global_pos, editable=True, act=None):
    """The menu of a right-click with nothing under the pointer (and of the Menu key): Cut, Copy, Paste, Select All, greyed
    when they cannot be done (nothing selected, the DNA locked, a view that cannot change the sequence). `act(name)` does
    the action ("cut", "copy", "paste", "select_all"); by default straight on the editor."""
    menu = QMenu(parent)
    add_edit_actions(menu, editor, editable, act)
    menu.exec(global_pos)
    return menu


def add_edit_actions(menu, editor, editable=True, act=None):
    """Cut, Copy, Paste, Select All added to `menu` (greyed as edit_menu says): also the end of the right-click menu on
    selected bases (feature_table.selection_menu), so both menus offer the same four."""
    if act is None:
        def act(name):
            if name in ("cut", "copy") and editor.end > editor.start:
                QGuiApplication.clipboard().setMimeData(editor.copy_mime())
                if name == "cut":
                    editor.backspace()
            elif name == "paste":
                editor.paste_mime(QGuiApplication.clipboard().mimeData())
            elif name == "select_all":
                editor.select(0, editor.length)
    selected, changes = editor.end > editor.start, editable and not editor.dna_locked
    clip = QGuiApplication.clipboard().mimeData()
    for name, text, key, enabled in (
            ("cut", "Cut", QKeySequence.Cut, selected and changes),
            ("copy", "Copy", QKeySequence.Copy, selected),
            ("paste", "Paste", QKeySequence.Paste, changes and clip is not None and clip.hasText()),
            (None, None, None, None),
            ("select_all", "Select All", QKeySequence.SelectAll, editor.length > 0)):
        if name is None:
            menu.addSeparator()
            continue
        if name in ("cut", "paste") and editor.dna_locked and editable:
            text += " (DNA locked)"   # says why it is greyed, as the variant menu of the alignment view does
        shortcut = QKeySequence(key).toString(QKeySequence.NativeText)
        action = menu.addAction(f"{text}\t{shortcut}", lambda n=name: act(n))
        action.setEnabled(bool(enabled))


class _Magnifier(QWidget):
    """Black-and-white magnifying glass (drawn in the text colour, so it follows the theme)."""
    def __init__(self):
        super().__init__()
        self.setFixedSize(18, 18)

    def paintEvent(self, event):
        p = QPainter(self)
        p.setRenderHint(QPainter.Antialiasing)
        p.setPen(QPen(self.palette().color(QPalette.Text), 1.8, Qt.SolidLine, Qt.RoundCap))
        p.drawEllipse(QRectF(2, 2, 9, 9))
        p.drawLine(QPointF(10, 10), QPointF(15.5, 15.5))


class SequenceView(QAbstractScrollArea):
    window_changed = Signal(int, int)  # bases currently on screen (for the genome view)
    editable = True   # False in a view that only shows a sequence: its menu greys Cut and Paste

    def __init__(self, editor):
        super().__init__()
        self._window = None
        self.editor = editor
        self._codon_cache = {}  # (feature index, editor version, table setting) -> [(middle base, amino acid, ...)]
        self._set_metrics(DEFAULT_SIZE)
        self.anchor = None  # fixed end of the selection (set by click/drag/shift+arrows)
        self.row_h, self.rows = 1, 1
        self.rpc, self.copies = 1, 1  # rows per copy of the sequence; 3+ copies when looping
        self.loop = False  # circular: endless scrolling, switched by clicking the red origin bar
        self.special, self.cum = [], [0]  # rows taller than row_h, running sum of their extra height
        self.lane_counts = {}    # row -> number of feature lanes used in that row
        self.row_fwd, self.row_rev = {}, {}  # same for forward / reverse primer lanes
        self.enzyme_rows = set()
        self._wheel = WheelSteps()
        # Dragging a selection past the top or bottom edge scrolls on while the button is held (see _autoscroll).
        self._drag_at = None   # where the pointer is while selecting with the left button, else None
        self.scroll_timer = QTimer(self, interval=40, timeout=self._autoscroll)
        self._origin_tip = False   # the origin bar's tooltip is up

        self.setFocusPolicy(Qt.StrongFocus)
        self.viewport().setCursor(Qt.IBeamCursor)
        editor.changed.connect(self._relayout)
        editor.display_changed.connect(self._relayout)
        editor.primers_changed.connect(self._relayout)
        tool_settings.notifier.changed.connect(self._code_changed)   # a bound method: dropped with this view
        editor.selection_changed.connect(self.viewport().update)
        editor.search_changed.connect(self.viewport().update)
        editor.reveal.connect(self.ensure_visible)
        editor.navigate.connect(self.center_on_position)
        self.hovered_site = None   # (start, end, top cut, bottom cut) of the restriction site being hovered, drawn on the bases
        self.hover = HoverTips(self.viewport(), editor, lambda pt: self.hover_target(pt.x(), pt.y()),
                               self._hover_changed, self._tip_anchor)
        self._build_zoom_bar()
        self._relayout()

    def _build_zoom_bar(self):
        """Zoom control down the right edge, in the strip no bases use: magnifier, "+", slider, "-"
        (same as Ctrl+wheel). The slider goes up for bigger letters."""
        bar = self.zoom_bar = QWidget(self.viewport())
        bar.setFixedWidth(ZOOM_BAR_W)
        column = QVBoxLayout(bar)
        column.setContentsMargins(0, 4, 0, 4)
        column.setSpacing(2)
        column.addWidget(_Magnifier(), 0, Qt.AlignHCenter)
        self.zoom_slider = QSlider(Qt.Vertical)
        self.zoom_slider.setRange(MIN_SIZE, MAX_SIZE)
        self.zoom_slider.setValue(self.size)
        self.zoom_slider.setFixedHeight(90)
        self.zoom_slider.setToolTip("Zoom (Ctrl+wheel)")
        self.zoom_slider.valueChanged.connect(self.set_size)
        buttons = []
        for text, step in (("+", 1), ("−", -1)):
            button = QToolButton()
            button.setText(text)
            button.setAutoRaise(True)
            button.setToolTip("Zoom in" if step > 0 else "Zoom out")
            button.clicked.connect(lambda _=False, s=step: self.set_size(self.size + s))
            buttons.append(button)
        column.addWidget(buttons[0], 0, Qt.AlignHCenter)
        column.addWidget(self.zoom_slider, 0, Qt.AlignHCenter)
        column.addWidget(buttons[1], 0, Qt.AlignHCenter)
        bar.adjustSize()
        self._place_zoom_bar()
        bar.show()

    def _place_zoom_bar(self):
        self.zoom_bar.move(self.viewport().width() - ZOOM_BAR_W, 0)

    @property
    def strand_h(self):
        """Height of the bases of a row: both strands, or the top one alone (Display > Complementary strand off), which is what
        a view shows when the second strand would only repeat the first (the Assembly's junction inspector)."""
        return self.lh * (2 if self.editor.show_complement else 1)

    def export_hidden(self):
        """What an exported image leaves out (image_export): the zoom bar and the scroll bar."""
        bar = self.verticalScrollBar()
        return [self.zoom_bar, bar.parentWidget() if bar.parentWidget() is not self else bar]

    # ---- letter size (Ctrl+wheel, Ctrl+= / Ctrl+- / Ctrl+0) -----------------------------
    def _set_metrics(self, size):
        self.size = size
        self.font = mono_font(size)
        self.labels = label_font(QApplication.font())
        self.labels.setPointSizeF(self.labels.pointSizeF() * size / DEFAULT_SIZE)
        self.label_metrics = QFontMetricsF(self.labels)
        fm = QFontMetricsF(self.font)
        self.cw0, self.lh, self.ascent = fm.horizontalAdvance("A"), fm.height(), fm.ascent()
        self.cw = self.cw0  # letter pitch: cw0 stretched so the rows fill the width (_update_left)
        self.digit_w = fm.horizontalAdvance("0")
        # Row numbers in a copy of the font without the letter spacing of the bases (see _update_left): spaced out, a long
        # number ran into the first base.
        self.number_font = mono_font(size)
        self.number_metrics = QFontMetricsF(self.number_font)
        self.base_left = 0  # set in _update_left (room for the position numbers)
        self.left = self.base_left  # x of the first base: base_left, plus room for a primer tail hanging off the start
        self.primer_h = self.label_metrics.height() + 8  # name + line + raised tail
        self.lane_h = self.lh
        # Text that shares a lane with other text (amino acids over a CDS's name, "PAM" over a guide's name) takes half
        # the lane.
        self.half_font = QFont(self.labels)
        self.half_font.setPixelSize(max(6, int(self.lane_h / 2)))
        self.verticalScrollBar().setSingleStep(int(self.lh))

    def set_size(self, size):
        """Letter size in points (6-32); the same bases stay at the top."""
        size = max(MIN_SIZE, min(MAX_SIZE, size))
        if self.zoom_slider.value() != size:  # (keeps the slider in step with Ctrl+wheel and keys)
            self.zoom_slider.blockSignals(True)
            self.zoom_slider.setValue(size)
            self.zoom_slider.blockSignals(False)
        if size != self.size:
            self._set_metrics(size)
            self._relayout()

    def wheelEvent(self, event):
        if event.modifiers() & Qt.ControlModifier:
            steps = self._wheel(event.angleDelta().y())
            if steps:
                self.set_size(self.size + steps)
            event.accept()
        else:
            super().wheelEvent(event)

    # ---- geometry ------------------------------------------------------------------
    def bases_per_line(self):
        n = int((self.viewport().width() - self.base_left - 4 - ZOOM_BAR_W) / self.cw0)
        return max(10, n // 10 * 10)

    def row_count(self):
        """Rows of one copy of the sequence."""
        return max(1, math.ceil(self.editor.length / self.bases_per_line()))

    def real_row(self, row):
        """Screen row -> row of the sequence (rows of the second copy repeat the first)."""
        return row % self.rpc

    def near(self, pos, v):
        """The copy of sequence position `pos` closest to screen position `v`."""
        n = self.editor.length
        return pos + n * round((v - pos) / n) if self.copies > 1 else pos

    def _update_left(self):
        """Room for the row numbers on the left, and the letter pitch: a row has a whole number of tens of bases, so
        the pitch is stretched a little (never squeezed) to make the rows exactly as wide as the free space, whatever
        the zoom. Tails that do not fit in a row go on in the next / previous one (see _draw_primer)."""
        self.base_left = self.number_metrics.horizontalAdvance(f"{self.editor.length + 1:,}") + 10  # the widest row number
        self.left = self.base_left
        free = self.viewport().width() - self.base_left - 4 - ZOOM_BAR_W
        self.cw = max(self.cw0, free / self.bases_per_line())
        self.font.setLetterSpacing(QFont.AbsoluteSpacing, self.cw - self.cw0)

    def _relayout(self):
        self._update_left()
        ed, bpl = self.editor, self.bases_per_line()
        # Keep the same bases at the top (rows change height and width with edits/resizes).
        sb = self.verticalScrollBar()
        top_row = self.row_at(sb.value())
        top = (self.real_row(top_row) * self._layout_bpl, sb.value() - self.offset(top_row)) \
            if getattr(self, "_layout_bpl", None) else (0, 0)
        self._layout_bpl = bpl
        rows = self.rpc = self.row_count()
        self.copies = 2 if ed.circular and ed.length and self.loop else 1  # endless scrolling: opt-in
        # Lanes used per row, by counting: per lane, merge the row ranges its features
        # cover, then a difference array (+1 at the first row, -1 after the last) summed
        # once. O(features + rows); which lanes exactly is only worked out for rows on
        # screen (row_lanes).
        by_lane = {}
        for f in ed.shown:
            rows_of = by_lane.setdefault(ed.lanes[f.index], [])
            for s, e in f.spans:
                rows_of.append((s // bpl, (max(e, s + 1) - 1) // bpl))
        diff = [0] * (rows + 1)
        for ranges in by_lane.values():
            ranges.sort()
            a, b = ranges[0]
            for c, d in ranges[1:] + [(rows + 1, rows + 1)]:
                if c > b + 1:  # gap: close the merged range
                    diff[min(a, rows)] += 1
                    diff[min(b + 1, rows)] -= 1
                    a, b = c, d
                else:
                    b = max(b, d)
        self.lane_counts, running = {}, 0
        for row in range(rows):
            running += diff[row]
            if running:
                self.lane_counts[row] = running
        self.enzyme_rows = {min(pos, ed.length - 1) // bpl for pos, _ in ed.shown_sites()}
        # a pinned cut (the Digest's selected fragment) writes its enzyme's name on the line above the bases too
        self.enzyme_rows |= {min(m.top, ed.length - 1) // bpl for marks in ed.cut_marks.values() for m in marks if m.pinned}
        # Primer lanes per row, forward (above the bases) and reverse (below), tails included.
        fwd, rev = {}, {}
        for site in ed.drawn_primers:
            used = fwd if site.strand > 0 else rev
            for s, e in site.extent:
                for row in range(s // bpl, (max(e, s + 1) - 1) // bpl + 1):
                    used.setdefault(row, set()).add(ed.primer_lanes[site.key])
        self.row_fwd = {row: sorted(used) for row, used in fwd.items()}
        self.row_rev = {row: sorted(used) for row, used in rev.items()}

        # Rounded up so the drawn parts (float sizes) always fit inside their row.
        self.row_h = math.ceil(0.8 * self.lh + self.strand_h + (3 * self.lh if ed.show_translation else 0))
        extra = {}
        for row, count in self.lane_counts.items():
            extra[row] = count * (self.lane_h + 2)
        for used_by_row in (self.row_fwd, self.row_rev):
            for row, used in used_by_row.items():
                extra[row] = extra.get(row, 0) + len(used) * self.primer_h
        for row in self.enzyme_rows - set(self.row_fwd):   # a row with a forward primer writes the names over its lane
            extra[row] = extra.get(row, 0) + self.lh
        if self.copies > 1:
            # Endless scrolling: a copy above and below the one in use, plus enough for a
            # window taller than the whole sequence (see scrollContentsBy).
            one = rows * self.row_h + sum(math.ceil(extra[r]) for r in extra if r < rows)
            self.copies = 2 + max(1, math.ceil(self.viewport().height() / max(1, one)))
        self.rows = rows * self.copies
        self.special = sorted(r + k * rows for r in extra if r < rows for k in range(self.copies))
        self.cum = [0] + list(accumulate(math.ceil(extra[r % rows]) for r in self.special))

        sb.setPageStep(self.viewport().height())
        sb.setRange(0, max(0, self.offset(self.rows) - self.viewport().height()))
        row = min(top[0] // bpl, rows - 1) + (rows if self.copies > 1 else 0)  # middle copy when circular
        anchor, self.anchor = self.anchor, None  # (not a scroll by the user: no anchor shift)
        sb.setValue(self.offset(row) + min(top[1], self.offset(row + 1) - self.offset(row) - 1))
        self._wrap_scroll()
        self.anchor = anchor
        self.viewport().update()

    def _wrap_scroll(self):
        """Circular: keep the scroll position within the second copy. The copies are
        identical, so moving by exactly one copy's height changes nothing on screen: the
        sequence scrolls on forever, both ways."""
        if self.copies < 2:
            return False
        sb, one = self.verticalScrollBar(), self.offset(self.rpc)
        v = sb.value()
        if v < one or v >= 2 * one:
            new = one + (v - one) % one
            if self.anchor is not None:  # a selection being made: its fixed end moves along
                self.anchor += round((new - v) / one) * self.editor.length
            sb.setValue(new)
            return True
        return False

    def _parts(self, row):
        """y offsets inside a row, top to bottom: (forward primers, bases, reverse
        primers, feature lanes). The enzyme line sits above the bases (its own line, or the forward primer lane nearest them)
        and the 3-frame translation between the reverse primers and the features."""
        fwd_top = self.lh if row in self.enzyme_rows and row not in self.row_fwd else 0   # names write over a primer lane
        seq_top = fwd_top + len(self.row_fwd.get(row, ())) * self.primer_h
        rev_top = seq_top + self.strand_h
        lanes = rev_top + len(self.row_rev.get(row, ())) * self.primer_h
        lanes += (3 * self.lh if self.editor.show_translation else 0) + 2
        return fwd_top, seq_top, rev_top, lanes

    def _name_y(self, row, pos, width):
        """y (inside the row) of the top of an enzyme name: always the line just above the bases, over a forward primer if
        there is one (`pos` and `width` are kept so the callers read the same)."""
        return max(0, self._parts(row)[1] - self.lh)

    def _primer_lane_top(self, row, strand, lane):
        """y (inside the row) of a primer lane; lane 0 is the one nearest the bases."""
        fwd_top, _, rev_top, _ = self._parts(row)
        if strand > 0:
            used = self.row_fwd.get(row, [])
            return fwd_top + (len(used) - 1 - used.index(lane)) * self.primer_h
        return rev_top + self.row_rev.get(row, []).index(lane) * self.primer_h

    def primer_at(self, x, y):
        """Key of the primer site drawn under a viewport point, or None."""
        pos, row, y_in_row = self.pos_at(x, y)
        row = self.real_row(row)
        base = row * self.bases_per_line() + int((x - self.left) // self.cw)
        ed = self.editor
        for site in ed.drawn_primers:
            lane = ed.primer_lanes[site.key]
            used = (self.row_fwd if site.strand > 0 else self.row_rev).get(row, [])
            if lane not in used or not any(s <= base < e for s, e in site.extent):
                continue
            top = self._primer_lane_top(row, site.strand, lane)
            if top <= y_in_row < top + self.primer_h:
                return site.key
        return None

    def _tip_anchor(self, pt, key):
        """(above, below) viewport points for the tooltip of a cut `key`: the top left corner of its enzyme name (the tip's
        bottom left corner goes there) and the bottom of the bases, for when there is no room above."""
        _, row, y_in_row = self.pos_at(pt.x(), pt.y())
        row = self.real_row(row)
        top = pt.y() - y_in_row
        _, seq_top, _, _ = self._parts(row)
        if key[0] == "site":
            pos, name = key[2], key[1]
        else:
            m = self.editor.cut_marks[key[1]][key[2]]
            pos, name = m.top, m.label
        name_top = self._name_y(row, pos, len(name)) if row in self.enzyme_rows else seq_top
        x = int(self.left + (pos - row * self.bases_per_line()) * self.cw)
        return QPoint(x, top + name_top), QPoint(x, top + seq_top + self.strand_h)

    def _hover_changed(self, key):
        """The tooltip of `key` came up (or went, None): draw the recognition site and the cut of a restriction enzyme."""
        from cut_sites import hover_site
        site = hover_site(self.editor, key) if isinstance(key, tuple) else None
        if site != self.hovered_site:
            self.hovered_site = site
            self.viewport().update()

    def hover_target(self, x, y):
        """What a tooltip is for under a viewport point: a primer (its key), a cut (a tuple, see hover.describe_cut) or a
        feature (its index), else None."""
        for found in (self.cut_at(x, y), self.site_label_at(x, y), self.primer_at(x, y)):   # an enzyme over a primer wins
            if found is not None:
                return found
        return self.feature_at(x, y)

    def cut_at(self, x, y):
        """Key ("cut", owner, i) of the cut mark (the strokes drawn through the bases where a digest or a guide cuts) under a
        viewport point, or None. The strokes are thin, so a few pixels on either side count."""
        ed = self.editor
        if not ed.cut_marks or not ed.length:
            return None
        pos, row, y_in_row = self.pos_at(x, y)
        row = self.real_row(row)
        _, seq_top, _, _ = self._parts(row)
        bpl = self.bases_per_line()
        base0 = row * bpl
        named = self._pinned_name_at(x, y_in_row, row, base0)   # the name of a pinned cut hovers like its strokes
        if named is not None:
            return named
        if not seq_top - 9 <= y_in_row <= seq_top + self.strand_h + 9:
            return None
        reach = 5
        mid = seq_top + self.lh
        for owner, marks in ed.cut_marks.items():
            for i, m in enumerate(marks):
                xt, xb = self.left + (m.top - base0) * self.cw, self.left + (m.bottom - base0) * self.cw
                if base0 <= m.top <= base0 + bpl and abs(x - xt) <= reach and y_in_row <= mid + 2:
                    return ("cut", owner, i)
                if base0 <= m.bottom <= base0 + bpl and abs(x - xb) <= reach and y_in_row >= mid - 2:
                    return ("cut", owner, i)
                if m.top != m.bottom and min(xt, xb) - reach <= x <= max(xt, xb) + reach and abs(y_in_row - mid) <= 4:
                    return ("cut", owner, i)
        return None

    def site_label_at(self, x, y):
        """Key ("site", enzyme name, cut position) of the enzyme name (View > Restriction Sites) or its tick under a point."""
        ed = self.editor
        if not ed.show_enzymes or not ed.length:
            return None
        _, row, y_in_row = self.pos_at(x, y)
        row = self.real_row(row)
        if row not in self.enzyme_rows:
            return None
        base0 = row * self.bases_per_line()
        for pos, name in ed.shown_sites():
            x0 = self.left + (pos - base0) * self.cw
            if (base0 <= pos <= base0 + self.bases_per_line() and x0 - 2 <= x <= x0 + (len(name) + 1) * self.cw
                    and 0 <= y_in_row - self._name_y(row, pos, len(name)) < self.lh):
                return ("site", name, pos)
        return None

    def row_lanes(self, row, feats=None):
        """Sorted feature lanes used in `row` (features drawn there are stacked in this
        order). `feats`: features already fetched for a range covering the row."""
        bpl = self.bases_per_line()
        a, b = row * bpl, min(self.editor.length, (row + 1) * bpl)
        if feats is None:
            feats = self.editor.visible_features(a, b)
        lanes = self.editor.lanes
        return sorted({lanes[f.index] for f in feats if f.index in lanes
                       and any(s < b and max(e, s + 1) > a for s, e in f.spans)})

    def offset(self, row):
        """y of the top of `row`; offset(row count) = total height."""
        return row * self.row_h + self.cum[bisect_left(self.special, row)]

    def row_at(self, y):
        """Row containing content y (scroll position included)."""
        lo, hi = 0, self.rows - 1
        while lo < hi:  # last row whose top is at or above y
            mid = (lo + hi + 1) // 2
            if self.offset(mid) <= y:
                lo = mid
            else:
                hi = mid - 1
        return lo

    def resizeEvent(self, event):
        super().resizeEvent(event)
        self._place_zoom_bar()
        self._relayout()

    def focusInEvent(self, event):
        super().focusInEvent(event)
        self.viewport().update()   # the cursor is drawn solid in the view that has the keyboard

    def focusOutEvent(self, event):
        super().focusOutEvent(event)
        self.viewport().update()

    def scrollContentsBy(self, dx, dy):
        if not self._wrap_scroll():  # (a wrap sets the value again, which comes back here)
            self.viewport().update()  # we repaint from the scrollbar value; nothing to move

    def ensure_visible(self, pos):
        """Scroll so sequence position `pos` shows, in whichever copy is nearest."""
        sb = self.verticalScrollBar()
        row = min(pos, max(0, self.editor.length - 1)) // self.bases_per_line()
        rows = [row + k * self.rpc for k in range(self.copies)]
        self._show_row(min(rows, key=lambda r: abs(self.offset(r) - sb.value())))

    def center_on_position(self, pos):
        """Another view (the alignment view) was clicked at `pos`: bring it to the middle of this one, unless it already is
        well inside (nothing moves then, so a click next to what you are reading does not make the page jump)."""
        sb, height = self.verticalScrollBar(), self.viewport().height()
        row = min(pos, max(0, self.editor.length - 1)) // self.bases_per_line()
        rows = [row + k * self.rpc for k in range(self.copies)]
        row = min(rows, key=lambda r: abs(self.offset(r) - sb.value()))
        top, bottom = self.offset(row), self.offset(row + 1)
        if top < sb.value() + 0.2 * height or bottom > sb.value() + 0.8 * height:
            sb.setValue(max(0, int((top + bottom) / 2 - height / 2)))

    def _show_row(self, row):
        sb = self.verticalScrollBar()
        top, bottom = self.offset(row), self.offset(row + 1)
        if top < sb.value():
            sb.setValue(top)
        elif bottom > sb.value() + self.viewport().height():
            sb.setValue(bottom - self.viewport().height())

    def _screen_row(self, v):
        """Screen row of screen position `v` (in the copy it belongs to)."""
        n, bpl = self.editor.length, self.bases_per_line()
        copy = min(self.copies - 1, max(0, v // n)) if n else 0
        return copy * self.rpc + min(max(0, v - copy * n), max(0, n - 1)) // bpl

    def _show_virtual(self, v):
        """Scroll to screen position `v` (in the copy it belongs to)."""
        self._show_row(self._screen_row(v))

    def pos_at(self, x, y):
        """(screen position between bases, screen row, y inside that row) under a viewport
        point. In the second copy of a circular sequence the position is length + the
        sequence position (Editor.select brings it back)."""
        bpl = self.bases_per_line()
        absolute = y + self.verticalScrollBar().value()
        row = self.row_at(absolute)
        copy, real = divmod(row, self.rpc)
        col = max(0, min(bpl, round((x - self.left) / self.cw)))
        return copy * self.editor.length + min(self.editor.length, real * bpl + col), row, absolute - self.offset(row)

    def feature_at(self, x, y):
        """Index of the feature drawn under a viewport point, or None."""
        pos, row, y_in_row = self.pos_at(x, y)
        row = self.real_row(row)
        lanes_y = self._parts(row)[3]
        used = self.row_lanes(row)
        slot = int((y_in_row - lanes_y) // (self.lane_h + 2)) if y_in_row >= lanes_y else -1
        if not 0 <= slot < len(used) or self.editor.length == 0:
            return None
        base = min(int((x - self.left) / self.cw) + row * self.bases_per_line(), self.editor.length - 1)
        for f in self.editor.visible_features(base, base + 1):
            if self.editor.lanes.get(f.index) == used[slot]:
                return f.index
        return None

    # ---- painting ------------------------------------------------------------------
    def paintEvent(self, event):
        ed = self.editor
        p = QPainter(self.viewport())
        p.setFont(self.font)
        pal = self.palette()
        p.fillRect(self.viewport().rect(), pal.color(QPalette.Base))
        text, dim = pal.color(QPalette.Text), pal.color(QPalette.PlaceholderText)
        sel_color = QColor(theme.selection_color())
        sel_color.setAlpha(80)

        bpl, y0 = self.bases_per_line(), self.verticalScrollBar().value()
        first, last = self.row_at(y0), self.row_at(y0 + self.viewport().height())
        # Screen rows in runs within one copy of the sequence: (first, last) screen row.
        runs = [(max(first, k * self.rpc), min(last, (k + 1) * self.rpc - 1))
                for k in range(first // self.rpc, last // self.rpc + 1)]
        window = (self.real_row(runs[0][0]) * bpl, min(ed.length, (self.real_row(runs[0][1]) + 1) * bpl))
        if window != self._window:  # the first copy on screen (the genome view frames it)
            self._window = window
            self.window_changed.emit(*window)
        selection = ed.selected_spans()
        with timed("sequence view: paint"):
            for a, b in runs:
                self._paint_run(p, a, b, y0, selection, text, dim, sel_color, pal)
        p.end()

    def _paint_run(self, p, first, last, y0, selection, text, dim, sel_color, pal):
        """Screen rows first..last, all in the same copy of the sequence."""
        ed, bpl = self.editor, self.bases_per_line()
        shift = first - self.real_row(first)  # screen row of real row 0 in this copy
        start, end = (first - shift) * bpl, min(ed.length, (last - shift + 1) * bpl)
        # Two extra bases each side so codons crossing a row edge can be translated.
        pad = min(2, start)
        seq = ed.doc.sequence(start - pad, end + 2)
        feats = ed.visible_features(start, end)
        primers = [site for site in ed.drawn_primers if any(s < end and e > start for s, e in site.extent)]
        sites = [(pos, name) for pos, name in ed.shown_sites() if start <= pos <= end]
        edited_color = theme.color("edited_fg")
        edited = [(a, b) for a, b in ed.doc.edited_ranges() if a < end and b > start]

        for screen_row in range(first, last + 1):
            row = screen_row - shift
            base0 = row * bpl
            i0 = base0 - start + pad
            line = seq[i0: i0 + bpl]
            n = len(line)
            _, seq_top, rev_top, lanes_y = self._parts(row)
            row_y = self.offset(screen_row) - y0
            y = row_y + seq_top  # top of the base letters

            for f in primers:
                used = (self.row_fwd if f.strand > 0 else self.row_rev).get(row, [])
                if ed.primer_lanes[f.key] in used:
                    top = row_y + self._primer_lane_top(row, f.strand, ed.primer_lanes[f.key])
                    self._draw_primer(p, f, base0, n, top, text, primers)
            p.setFont(self.font)
            if row in self.enzyme_rows:   # after the primers: a name over a primer is drawn on top of it
                self._draw_sites(p, [s for s in sites if base0 <= s[0] <= base0 + n and
                                     (s[0] < base0 + bpl or base0 + n == ed.length)], row, row_y, dim)
                self._draw_pinned_names(p, base0, n, row, row_y)
                p.setFont(self.font)

            # Bases edited since the last save: their letters are drawn blue further down (no background, so they
            # can never be mistaken for the yellow search hits); `marks` are their columns in this row.
            marks = [(max(s, base0) - base0, min(e, base0 + n) - base0) for s, e in edited if max(s, base0) < min(e, base0 + n)]
            for s, e in ed.hits_between(base0, base0 + n):   # everything the search bar found, on the strand it matched
                a, b = max(s, base0), min(e, base0 + n)
                strand = ed.search_strand.get((s, e), 0)
                x0, w = self.left + (a - base0) * self.cw, (b - a) * self.cw
                if strand >= 0:   # top strand (or both)
                    p.fillRect(QRectF(x0, y, w, self.lh if strand else self.strand_h), theme.hit_color())
                if strand < 0:    # bottom strand only (the one line there is, when only the top strand is shown)
                    p.fillRect(QRectF(x0, y + self.strand_h - self.lh, w, self.lh), theme.hit_color())
            # selection (two parts when it passes the origin)
            for s, e in selection:
                a, b = max(s, base0), min(e, base0 + n)
                if a < b:
                    p.fillRect(QRectF(self.left + (a - base0) * self.cw, y, (b - a) * self.cw, self.strand_h), sel_color)
            self._draw_hovered_site(p, base0, n, y)
            # [ ] around the selection: drag either one to move that end
            if ed.end > ed.start:
                first, last = selection[0][0], selection[-1][1]
                if base0 <= first < base0 + max(n, 1):
                    self._draw_bracket(p, self.left + (first - base0) * self.cw, y, 1, pal)
                if base0 < last <= base0 + n:
                    self._draw_bracket(p, self.left + (last - base0) * self.cw, y, -1, pal)

            # origin of a circular sequence: red bar after the last base (click: endless
            # scrolling on/off), and before base 1 of the next copy when looping
            if ed.circular:
                p.setPen(QPen(QColor(ORIGIN_COLOR), 2))
                if base0 + n == ed.length:
                    x = self.left + n * self.cw + 1
                    p.drawLine(QPointF(x, y - 2), QPointF(x, y + self.strand_h + 2))
                if row == 0 and screen_row > 0:
                    p.drawLine(QPointF(self.left - 2, y - 2), QPointF(self.left - 2, y + self.strand_h + 2))

            p.setPen(dim)
            p.setFont(self.number_font)   # right-aligned, ending just left of the first base
            number = f"{base0 + 1:,}"
            p.drawText(QPointF(self.base_left - 6 - self.number_metrics.horizontalAdvance(number), y + self.ascent), number)
            p.setFont(self.font)
            self._draw_letters(p, line, marks, y + self.ascent, text, edited_color)
            if ed.show_complement:
                self._draw_letters(p, complement(line), marks, y + self.lh + self.ascent, dim, edited_color)
            for a, b in marks:   # thin blue underline under the edited letter pairs
                p.fillRect(QRectF(self.left + a * self.cw, y + self.strand_h - 1, (b - a) * self.cw, 1), edited_color)
            self._draw_cuts(p, base0, n, y)

            if ed.show_translation:
                trans_y = row_y + rev_top + len(self.row_rev.get(row, ())) * self.primer_h
                self._draw_translation(p, seq, i0, base0, n, trans_y, dim)

            slots = {lane: i for i, lane in enumerate(self.row_lanes(row, feats))}
            for f in feats:
                slot = slots.get(ed.lanes.get(f.index))
                if slot is not None:
                    self._draw_feature(p, f, base0, n, row_y + lanes_y + slot * (self.lane_h + 2))

            # cursor last, so it is drawn over the primers and features: a line through the whole row, its
            # base number above. It belongs to this row unless it sits at the row end and another row follows.
            c = ed.start
            if ed.start == ed.end and base0 <= c <= base0 + n and (c < base0 + bpl or base0 + n == ed.length):
                x = self.left + (c - base0) * self.cw
                # Faint while another view has the keyboard: with several panes, the one where typing goes is the one
                # whose cursor is solid. (The selection is shared by every view of the document: it stays as it is.)
                caret = QColor(theme.selection_color())
                if not self.hasFocus():
                    caret.setAlpha(90)
                p.setPen(QPen(caret, 1.5))
                p.drawLine(QPointF(x, row_y), QPointF(x, self.offset(screen_row + 1) - y0))
                self._draw_bracket_number(p, x, y, 1, c + 1, pal, caret)
                p.setFont(self.font)

            # the selection's end numbers last of all: over every primer, feature and cut
            if ed.end > ed.start:
                first, last = selection[0][0], selection[-1][1]
                if base0 <= first < base0 + max(n, 1):
                    self._draw_bracket_number(p, self.left + (first - base0) * self.cw, y, 1, first % ed.length + 1, pal)
                if base0 < last <= base0 + n:
                    self._draw_bracket_number(p, self.left + (last - base0) * self.cw, y, -1, (last - 1) % ed.length + 1, pal)
                p.setFont(self.font)

    def _draw_letters(self, p, letters, marks, baseline, color, edited_color):
        """One line of base letters in `color`, except the columns in `marks` ((first, end) pairs), which are blue.
        The font is monospaced, so a blanked-out letter leaves its column free for the blue one."""
        gaps = list(letters)
        for a, b in marks:
            gaps[a:b] = " " * (b - a)
        p.setPen(color)
        p.drawText(QPointF(self.left, baseline), "".join(gaps))
        p.setPen(edited_color)
        for a, b in marks:
            p.drawText(QPointF(self.left + a * self.cw, baseline), letters[a:b])

    def _draw_primer(self, p, f, base0, n, top, text_color, others=()):
        """A primer as a thin line along the bases it anneals to, with a half arrow at the
        3' end; its 5' tail (bases that do not anneal) runs on a line further from the DNA.
        Forward primers sit above the bases, reverse primers below the complement. `others`: the primers drawn around it
        (its name keeps clear of those in its lane)."""
        ed = self.editor
        fwd = f.strand > 0
        line_y = top + self.primer_h - 3 if fwd else top + 3
        tail_y = line_y - 5 if fwd else line_y + 5
        away = -1 if fwd else 1  # direction away from the DNA
        selected = ed.current_primer == f.key
        color = theme.selection_color() if selected else (
            QColor(primer_color(f)))
        pen = QPen(color, 3 if selected else 2)
        if f.mismatches:
            pen.setStyle(Qt.DashLine)  # "almost binds"
        p.setPen(pen)
        length, bpl = ed.length, self.bases_per_line()
        last_row = base0 + n >= length and not ed.circular
        row_end = base0 + (bpl if last_row else n)   # (the last row of a linear sequence has room for a tail)
        low = base0

        def x(pos):
            return self.left + (pos - base0) * self.cw

        def pieces(a, b):
            """What of positions a..b (they may be before the start or past the end) is drawn in this row. Beyond the
            ends of the sequence a tail goes on in the previous / next row: across the origin of a circular
            sequence, and, for a linear one, in the first row's (last row's) own lane at the other end of the row."""
            if ed.circular:
                out = [(max(a, 0), min(b, length))]
                if a < 0:
                    out.append((a + length, min(b, 0) + length))
                if b > length:
                    out.append((max(a, length) - length, b - length))
                return out
            out = [(max(a, 0), b)]
            if a < 0 and base0 == 0:
                out.append((a + bpl, min(b, 0) + bpl))
            if b > base0 + bpl and last_row:
                out.append((max(a, base0 + bpl) - bpl, b - bpl))
            return out

        def segment(a, b, y):
            found = []
            for pa, pb in pieces(a, b):
                pa, pb = max(pa, low), min(pb, row_end)
                if pa < pb:
                    p.drawLine(QPointF(x(pa), y), QPointF(x(pb), y))
                    found.append((pa, pb))
            return found

        drawn = [d for a, b in f.spans for d in segment(a, b, line_y)]
        tail = f.tail
        if tail:
            five = f.spans[0][0] if fwd else f.spans[-1][1]  # where the tail leaves the DNA
            segment(five - tail, five, tail_y) if fwd else segment(five, five + tail, tail_y)
            if base0 <= five <= row_end:
                p.setRenderHint(QPainter.Antialiasing, True)   # the slanted strokes only: straight lines stay crisp
                p.drawLine(QPointF(x(five), line_y), QPointF(x(five) + (-4 if fwd else 4), tail_y))
                p.setRenderHint(QPainter.Antialiasing, False)
        three = f.spans[-1][1] if fwd else f.spans[0][0]
        if base0 <= three <= row_end:
            p.setRenderHint(QPainter.Antialiasing, True)
            p.drawLine(QPointF(x(three), line_y), QPointF(x(three) + (-6 if fwd else 6), line_y + 5 * away))
            p.setRenderHint(QPainter.Antialiasing, False)
        # The primer's own bases where they do not pair (mismatches, on its line) or hang off
        # the template (5' tail, on the tail line), in orange: what the primer brings in.
        primer = f.sequence
        positions = [q for s, e in f.spans for q in range(s, e)]
        annealed = primer[len(primer) - len(positions):]
        facing = annealed if fwd else annealed[::-1]  # primer base at each template base, left to right
        wrong = set(f.mismatch_at)
        letters = [(q, facing[k], line_y) for k, q in enumerate(positions) if q in wrong and k < len(facing)]
        if tail:
            hang = primer[:tail] if fwd else primer[:tail][::-1]
            first = five - tail if fwd else five
            letters += [(first + i, c, tail_y) for i, c in enumerate(hang)]
        p.setFont(self.labels)
        h = self.label_metrics.height()
        for pos, letter, ly in letters:
            for at, _ in pieces(pos, pos + 1):
                if low <= at < row_end:
                    box = QRectF(x(at), ly - h / 2, self.cw, h)
                    p.fillRect(box, self.palette().color(QPalette.Base))  # over the line, readable
                    p.setPen(QColor(ALMOST_COLOR))
                    p.drawText(box, Qt.AlignCenter, letter)
        p.setFont(self.font)
        if drawn:
            a, b = max(drawn, key=lambda s: s[1] - s[0])
            label_h = self.primer_h - 8
            label = f.label + (" (temporary)" if f.temporary else "") + (f" ({f.mismatches} mismatch)" if f.mismatches else "")
            # The name starts at the 5' end of its line and reads toward the 3' end, running on past it by up to 6 letters:
            # never toward the 5' tail, whose bases are written beside the line there (a reverse primer's name is
            # right-aligned for that), and never into the next primer of its lane.
            lane = ed.primer_lanes[f.key]
            near = [piece for o in others if o.key != f.key and o.strand == f.strand and ed.primer_lanes.get(o.key) == lane
                    for piece in o.extent]
            if fwd:
                x0 = x(a)
                x1 = min([x(b) + 6 * self.cw] + [x(s) - 4 for s, _ in near if s >= b])
            else:
                x1 = x(b)
                x0 = max([x(a) - 6 * self.cw, 2.0] + [x(e) + 4 for _, e in near if e <= a])
            rect = QRectF(x0, top if fwd else line_y + 4, max(x1 - x0, 0.0), label_h)
            name = self.label_metrics.elidedText(label, Qt.ElideRight, rect.width())
            p.setFont(self.labels)
            p.setPen(text_color)
            p.drawText(rect, (Qt.AlignLeft if fwd else Qt.AlignRight) | Qt.AlignVCenter, name)

    def _draw_cuts(self, p, base0, n, y):
        """Cuts (digest, Cas9) on both strands: a stroke through the top strand where it is cut and one through the bottom
        strand where that one is cut, each ending in a small cross (scissors) outside the bases, joined across the middle
        when they differ (a sticky end). The cross keeps a cut apart from the origin mark of a circular sequence, a plain bar."""
        ed = self.editor
        if not ed.cut_marks:
            return
        bpl = self.bases_per_line()
        # A gap at the end of a row belongs to the next row, unless this is the last one (same rule as the cursor).
        here = lambda q: base0 <= q <= base0 + n and (q < base0 + bpl or base0 + n == ed.length)
        mid = y + self.lh
        color = theme.color("cut")
        arm = 2.5   # half size of the cross
        for marks in ed.cut_marks.values():
            for m in marks:
                xt, xb = self.left + (m.top - base0) * self.cw, self.left + (m.bottom - base0) * self.cw
                p.setPen(QPen(color, 4 if m.pinned else 2))   # pinned: as if hovered, a heavier stroke
                if here(m.top):
                    p.drawLine(QPointF(xt, y - 1), QPointF(xt, mid))
                if here(m.bottom):
                    p.drawLine(QPointF(xb, mid), QPointF(xb, y + self.strand_h + 1))
                if m.top != m.bottom and here(m.top) and here(m.bottom):
                    p.setPen(QPen(color, 1.2, Qt.DotLine))
                    p.drawLine(QPointF(xt, mid), QPointF(xb, mid))
                p.setPen(QPen(color, 1.4))
                for q, cy in ((m.top, y - 2 - arm), (m.bottom, y + self.strand_h + 2 + arm)):
                    if here(q):
                        cx = self.left + (q - base0) * self.cw
                        p.drawLine(QPointF(cx - arm, cy - arm), QPointF(cx + arm, cy + arm))
                        p.drawLine(QPointF(cx - arm, cy + arm), QPointF(cx + arm, cy - arm))

    def _draw_hovered_site(self, p, base0, n, y):
        """The recognition site of the hovered restriction cut as a band over both strands, and its cut: a heavy stroke through
        each strand where it is cut, joined across the middle for a sticky end (the picture of the tooltip, on the bases)."""
        if self.hovered_site is None or not self.editor.length:
            return
        s, e, ut, ub = self.hovered_site
        L = self.editor.length
        band = QColor(theme.color("pill_site"))
        band.setAlpha(170)
        for shift in ((0, -L, L) if self.editor.circular else (0,)):   # a site (or cut) across the origin shows on both ends
            a, b = max(s + shift, base0), min(e + shift, base0 + n)
            if a < b:
                p.fillRect(QRectF(self.left + (a - base0) * self.cw, y, (b - a) * self.cw, self.strand_h), band)
            xt, xb = (self.left + (q + shift - base0) * self.cw for q in (ut, ub))
            mid = y + self.lh
            color = theme.color("cut")
            lo, hi = self.left - 1, self.left + n * self.cw + 1
            p.setPen(QPen(color, 4))
            if lo <= xt <= hi:
                p.drawLine(QPointF(xt, y - 1), QPointF(xt, mid))
            if lo <= xb <= hi:
                p.drawLine(QPointF(xb, mid), QPointF(xb, y + self.strand_h + 1))
            if ut != ub and lo <= xt <= hi and lo <= xb <= hi:
                p.setPen(QPen(color, 1.5, Qt.DotLine))
                p.drawLine(QPointF(xt, mid), QPointF(xb, mid))

    def _draw_pinned_names(self, p, base0, n, row, top_of_row):
        """The enzyme name of each pinned cut (see Editor.set_cut_marks) on the line above the bases, in the cut colour, on
        the base colour so it stays readable over the names of View > Restriction Sites."""
        ed = self.editor
        bpl = self.bases_per_line()
        bold = QFont(self.font)
        bold.setBold(True)
        fm = QFontMetricsF(bold)
        for marks in ed.cut_marks.values():
            for m in marks:
                if not (m.pinned and base0 <= m.top <= base0 + n and (m.top < base0 + bpl or base0 + n == ed.length)):
                    continue
                x = self.left + (m.top - base0) * self.cw
                row_y = top_of_row + self._name_y(row, m.top, len(m.label))
                box = QRectF(x - 1, row_y + 1, fm.horizontalAdvance(m.label) + 4, self.lh - 2)
                p.fillRect(box, self.palette().color(QPalette.Base))
                p.setFont(bold)
                p.setPen(theme.color("cut"))
                p.drawText(QPointF(x + 1, row_y + self.ascent - 2), m.label)
                p.setPen(QPen(theme.color("cut"), 2))
                p.drawLine(QPointF(x, row_y + self.lh - 3), QPointF(x, row_y + self.lh + 1))
        p.setFont(self.font)

    def _pinned_name_at(self, x, y_in_row, row, base0):
        """Key of the pinned cut whose name is drawn under (x, y_in_row) (the enzyme line of the row), or None."""
        if row not in self.enzyme_rows:
            return None
        for owner, marks in self.editor.cut_marks.items():
            for i, m in enumerate(marks):
                x0 = self.left + (m.top - base0) * self.cw
                if (m.pinned and x0 - 1 <= x <= x0 + (len(m.label) + 0.5) * self.cw and base0 <= m.top <= base0 + self.bases_per_line()
                        and 0 <= y_in_row - self._name_y(row, m.top, len(m.label)) < self.lh):
                    return ("cut", owner, i)
        return None

    def _draw_sites(self, p, sites, row, row_y, color):
        """Enzyme names on the line above the bases (raised over the primers when one is under the name), with a tick at the
        cut."""
        p.setPen(color)
        base0 = row * self.bases_per_line()
        free_x = {}   # per line: where the last name written there ends
        for pos, name in sites:
            x = self.left + (pos - base0) * self.cw
            top = row_y + self._name_y(row, pos, len(name))
            p.drawLine(QPointF(x, top + self.lh - 3), QPointF(x, top + self.lh + 1))
            if x >= free_x.get(top, -1.0):  # skip a name that would overlap the previous one
                p.drawText(QPointF(x + 1, top + self.ascent - 2), name)
                free_x[top] = x + (len(name) + 1) * self.cw

    def _draw_translation(self, p, seq, i0, base0, n, y, color):
        """Frames 1-3 of the top strand (frame 1 starts at base 1), one amino acid under
        the middle base of each codon (genetic code: the Tools > Settings choice)."""
        p.setPen(color)
        code = genetic_code.codon_map()
        for frame in range(3):
            # First codon whose middle base is in this row.
            c = base0 - 1 + (frame - (base0 - 1)) % 3
            chars = [" "] * n
            while c + 1 < base0 + n:
                j = c - base0 + i0
                if c >= 0 and j + 3 <= len(seq):
                    codon = seq[j:j + 3].upper().replace("U", "T")
                    chars[c + 1 - base0] = code.get(codon, "X")
                c += 3
            line = "".join(chars)
            p.setPen(color)
            p.drawText(QPointF(self.left, y + frame * self.lh + self.ascent), line.replace("*", " "))
            for i, aa in enumerate(line):
                if aa == "*":  # stops stand out: a red block over the codon
                    self._draw_stop(p, self.left + (i - 1) * self.cw, y + frame * self.lh + 1, self.lh - 2)

    def _draw_bracket(self, p, x, y, side, pal):
        """Selection bracket over both strands: "[" (side 1) at the start, "]" (-1) at the end."""
        p.setPen(QPen(theme.selection_color(), 2))
        top, bottom, arm = y - 2, y + self.strand_h + 2, side * min(5.0, self.cw * 0.6)
        p.drawPolyline(QPolygonF([QPointF(x + arm, top), QPointF(x, top), QPointF(x, bottom), QPointF(x + arm, bottom)]))

    def _draw_bracket_number(self, p, x, y, side, number, pal, color=None):
        """Tiny base number at a bracket: first selected base above the "[" (text starts at the bracket),
        last selected base below the "]" (text ends at the bracket), so the two never overlap."""
        small = QFont(self.labels)
        small.setPixelSize(max(7, int(self.lh * 0.5)))
        p.setFont(small)
        fm = QFontMetricsF(small)
        text = f"{number:,}"
        w = fm.horizontalAdvance(text)
        left = x if side > 0 else x - w
        top = y - 3 - fm.height() if side > 0 else y + self.strand_h + 3
        box = QRectF(left, top, w, fm.height())
        p.fillRect(box, pal.color(QPalette.Base))  # readable over whatever is under it
        p.setPen(color or theme.selection_color())
        p.drawText(QPointF(left, top + fm.ascent()), text)

    def origin_at(self, x, y):
        """Is (x, y) on the red origin bar of a circular sequence (end of the last row, or
        start of base 1's row in a later copy)?"""
        ed = self.editor
        if not ed.circular or not ed.length:
            return False
        _, row, y_in_row = self.pos_at(x, y)
        real = self.real_row(row)
        _, seq_top, _, _ = self._parts(real)
        if not seq_top - 3 <= y_in_row <= seq_top + self.strand_h + 3:
            return False
        if real == self.rpc - 1:  # the last row: bar right after its last base
            n_last = ed.length - real * self.bases_per_line()
            if abs(x - (self.left + n_last * self.cw + 1)) <= 4:
                return True
        return real == 0 and row > 0 and abs(x - (self.left - 2)) <= 4

    def toggle_loop(self):
        """Endless scrolling on / off (circular sequences)."""
        self.loop = not self.loop
        self._relayout()
        self.ensure_visible(self.editor.start)

    def bracket_at(self, x, y):
        """Screen position of the selection end whose bracket is under (x, y), with the
        other end's screen position: (end grabbed, other end), or None."""
        ed = self.editor
        if ed.end <= ed.start:
            return None
        pos, row, y_in_row = self.pos_at(x, y)
        _, seq_top, _, _ = self._parts(self.real_row(row))
        if not seq_top - 3 <= y_in_row <= seq_top + self.strand_h + 3:
            return None
        boundary_x = self.left + round((x - self.left) / self.cw) * self.cw  # x of the gap `pos`
        if abs(x - boundary_x) > self.cw * 0.4:
            return None  # over a base, not on the gap where a bracket is
        for here, other in ((ed.start, ed.end), (ed.end, ed.start)):
            v = self.near(here, pos)
            if v == pos:
                return v, v + (other - here)
        return None

    def _draw_stop(self, p, x, y, h):
        """Stop codon: red block 3 bases wide with a white "*", so early stops catch the eye."""
        p.fillRect(QRectF(x, y, 3 * self.cw, h), QColor(STOP_COLOR))
        font = QFont(self.labels)
        font.setPixelSize(max(8, int(h * 1.1)))
        font.setBold(True)
        p.setFont(font)
        p.setPen(QColor("white"))
        p.drawText(QRectF(x, y + h * 0.15, 3 * self.cw, h), Qt.AlignCenter, "*")
        p.setFont(self.font)

    def _draw_cas9(self, p, f, layout, base0, n, y):
        """A saved Cas9 / Cas12a guide site: the protospacer (orange box) and the PAM (yellow box) side by side, a blue bar
        under both with the guide's name (all text black), and a hollow red arrow pointing up at the cut (Cas9: 3 bp from the PAM,
        inside the protospacer)."""
        ed, h = self.editor, self.lane_h
        hb = h * 0.5   # the boxes take the top half of the lane, the name bar the bottom half
        x = lambda q: self.left + (q - base0) * self.cw
        clip = lambda parts: [(max(a, base0), min(b, base0 + n)) for a, b in parts if a < base0 + n and b > base0]
        proto, pam, label = theme.color("cas9_proto"), theme.color("cas9_pam"), theme.color("cas9_label")
        for parts, color in ((layout.proto, proto), (layout.pam, pam)):
            p.setPen(color.darker(150))
            p.setBrush(color)
            for a, b in clip(parts):
                p.drawRect(QRectF(x(a), y, x(b) - x(a), hb))
        bars = clip(runs(layout.cells))   # one bar, or two where the site crosses the origin
        p.setPen(label.darker(140))
        p.setBrush(label)
        for a, b in bars:
            p.drawRect(QRectF(x(a), y + hb, x(b) - x(a), h - hb))
        text = theme.color("cas9_text")   # black on the orange / yellow / blue fills, in both themes
        if bars:
            a, b = max(bars, key=lambda r: r[1] - r[0])
            p.setPen(text)
            p.setFont(self.half_font)
            name = QFontMetricsF(p.font()).elidedText(f.label, Qt.ElideRight, x(b) - x(a) - 4)
            p.drawText(QRectF(x(a), y + hb, x(b) - x(a), h - hb), Qt.AlignCenter, name)
            for a, b in clip(layout.pam):   # "PAM" in its box when it fits
                if x(b) - x(a) >= QFontMetricsF(p.font()).horizontalAdvance("PAM") + 4:
                    p.drawText(QRectF(x(a), y, x(b) - x(a), hb), Qt.AlignCenter, "PAM")
            p.setFont(self.font)
        # hollow red arrow on each cut (the gap between two bases), pointing UP at the cut: outline only, so the box stays
        # visible under it. The head is at the top of the lane, next to the bases; the shaft runs down through the lane.
        p.setPen(QPen(theme.color("cas9_arrow"), 1.6))
        p.setBrush(Qt.NoBrush)
        p.setRenderHint(QPainter.Antialiasing, True)   # its slanted head
        bpl = self.bases_per_line()
        shaft, head = max(2.0, self.cw * 0.3), max(5.0, self.cw * 0.9)
        for q in layout.gaps:
            if not (base0 <= q <= base0 + n and (q < base0 + bpl or base0 + n == ed.length)):
                continue   # (a gap at the very end of a circle, its origin, is drawn at the end of the last row)
            cx, neck, base = x(q), y + hb * 0.55, y + h - 1
            p.drawPolygon(QPolygonF([QPointF(cx, y), QPointF(cx + head / 2, neck), QPointF(cx + shaft / 2, neck),
                                     QPointF(cx + shaft / 2, base), QPointF(cx - shaft / 2, base), QPointF(cx - shaft / 2, neck),
                                     QPointF(cx - head / 2, neck)]))
        p.setRenderHint(QPainter.Antialiasing, False)

    def _draw_unit(self, p, f, g, base0, n, y):
        """A transcription unit (units.py): a thin line over its promoter, a thicker one from +1 on, a T bar where the
        terminator stops it (an arrow head when nothing does), its name above the line from +1."""
        h, color = self.lane_h, self.editor.colors[f.index].darker(130)
        x = lambda q: self.left + (q - base0) * self.cw
        ym = y + h * 0.72
        thin, thick = max(1.0, h * 0.07), max(3.0, h * 0.32)
        p.setBrush(Qt.NoBrush)
        for a, b, heavy in g.pieces:
            a, b = max(a, base0), min(b, base0 + n)
            if a < b:
                p.setPen(QPen(color, thick if heavy else thin, Qt.SolidLine, Qt.FlatCap))
                p.drawLine(QPointF(x(a), ym), QPointF(x(b), ym))
        if base0 <= g.tss <= base0 + n:   # +1: a short tick up where the line widens
            p.setPen(QPen(color, thin, Qt.SolidLine, Qt.FlatCap))
            p.drawLine(QPointF(x(g.tss), ym), QPointF(x(g.tss), y + h * 0.42))
        if base0 <= g.end <= base0 + n and not (g.end == base0 + n and g.strand < 0):
            ex = x(g.end)
            if g.stopped:   # the T of the terminator
                p.setPen(QPen(color, thick, Qt.SolidLine, Qt.FlatCap))
                p.drawLine(QPointF(ex, y + h * 0.35), QPointF(ex, y + h))
            else:           # nothing stops it: it runs on
                p.setPen(Qt.NoPen)
                p.setBrush(color)
                way = 1 if g.strand > 0 else -1
                p.setRenderHint(QPainter.Antialiasing, True)
                p.drawPolygon(QPolygonF([QPointF(ex + way * h * 0.3, ym), QPointF(ex - way * h * 0.2, y + h * 0.45),
                                         QPointF(ex - way * h * 0.2, y + h)]))
                p.setRenderHint(QPainter.Antialiasing, False)
        # The name, small, above the line: from +1 toward the 3' end, after a "+1" when the start base is known exactly.
        if base0 <= g.tss < base0 + n or g.strand < 0 and base0 < g.tss <= base0 + n:
            p.save()
            p.setFont(self.half_font)
            p.setPen(color)
            fm = QFontMetricsF(self.half_font)
            text = f"+1  {f.label}" if g.exact else f.label
            if g.strand > 0:
                p.drawText(QPointF(x(g.tss) + 2, y + fm.ascent()), text)
            else:
                text = f"{f.label}  +1" if g.exact else f.label
                p.drawText(QPointF(x(g.tss) - 2 - fm.horizontalAdvance(text), y + fm.ascent()), text)
            p.restore()

    def _draw_feature(self, p, f, base0, n, y):
        g = units.geometry(f, self.editor.length, self.editor.doc.circular)
        if g is not None:
            return self._draw_unit(p, f, g, base0, n, y)
        if any(k == CAS9_KEY for k, _ in f.qualifiers):
            layout = feature_layout(f, self.editor.length)
            if layout is not None:   # (None: the feature was resized by hand: it is drawn as any other)
                return self._draw_cas9(p, f, layout, base0, n, y)
        color = self.editor.colors[f.index]
        h = self.lane_h
        for i, (fs, fe) in enumerate(f.spans):
            a, b = max(fs, base0), min(fe, base0 + n)
            if a >= b:
                continue
            x0, x1 = self.left + (a - base0) * self.cw, self.left + (b - base0) * self.cw
            # Arrow head on the end of the feature that lies in this row.
            tip, way = min(h / 2, x1 - x0), arrow(f)
            if way > 0 and i == len(f.spans) - 1 and fe <= base0 + n:
                pts = [(x0, y), (x1 - tip, y), (x1, y + h / 2), (x1 - tip, y + h), (x0, y + h)]
            elif way < 0 and i == 0 and fs >= base0:
                pts = [(x0 + tip, y), (x1, y), (x1, y + h), (x0 + tip, y + h), (x0, y + h / 2)]
            else:
                pts = [(x0, y), (x1, y), (x1, y + h), (x0, y + h)]
            shape = QPolygonF([QPointF(*pt) for pt in pts])
            p.setPen(color.darker(140))
            p.setBrush(color)
            p.setRenderHint(QPainter.Antialiasing, len(pts) == 5)   # an arrow head is slanted: smooth; a box stays crisp
            p.drawPolygon(shape)
            p.setRenderHint(QPainter.Antialiasing, False)
            translated = self.editor.translated(f)
            if translated:  # amino acids in the top half (next to the DNA), name below
                # Clipped to the arrow, so a stop codon's red block takes the arrow head's
                # shape instead of hiding it; then the outline again on top.
                outline = QPainterPath()
                outline.addPolygon(shape)
                p.setClipPath(outline)
                self._draw_amino_acids(p, f, base0, a, b, y, color)
                p.setClipping(False)
                p.setPen(color.darker(140))
                p.setBrush(Qt.NoBrush)
                p.setRenderHint(QPainter.Antialiasing, len(pts) == 5)
                p.drawPolygon(shape)
                p.setRenderHint(QPainter.Antialiasing, False)
            if x1 - x0 > 3 * self.cw:
                p.setPen(text_color_on(color))
                p.setFont(self.half_font if translated else self.labels)
                top, height = (y + h / 2, h / 2) if translated else (y, h)
                label = QFontMetricsF(p.font()).elidedText(f.label, Qt.ElideRight, x1 - x0 - 2 * tip)
                p.drawText(QRectF(x0 + tip, top, x1 - x0 - 2 * tip, height), Qt.AlignCenter, label)
                p.setFont(self.font)

    def _draw_amino_acids(self, p, f, base0, a, b, y, color):
        """Top half of a translated CDS arrow, over its bases [a, b) in this row: each amino
        acid under the middle base of its codon, and small ticks along the top edge at the
        codon boundaries, so each letter's three bases can be seen."""
        p.setFont(self.half_font)
        p.setPen(text_color_on(color))
        h = self.lane_h / 2
        for mid, aa, lo, hi in self._codons(f):
            if hi <= a or lo >= b:
                continue
            if a <= mid < b:
                cx = self.left + (mid - base0 + 0.5) * self.cw
                if aa == "*":
                    self._draw_stop(p, cx - 1.5 * self.cw, y + 1, self.lane_h - 2)
                    p.setFont(self.half_font)
                    p.setPen(text_color_on(color))
                    continue
                p.drawText(QRectF(cx - 1.5 * self.cw, y + 1, 3 * self.cw, h), Qt.AlignCenter, aa)
            for edge in (lo, hi):  # ticks between codons
                if a <= edge <= b:
                    x = self.left + (edge - base0) * self.cw
                    p.drawLine(QPointF(x, y + 1), QPointF(x, y + 3))
        p.setFont(self.font)

    def _code_changed(self, tool):
        """The genetic code setting changed: amino acids are drawn again (the cache is keyed by the table, so it stays valid)."""
        if tool == genetic_code.TOOL:
            self.viewport().update()

    def _codons(self, f):
        """[(middle base, amino acid, first base, last base + 1)] of a CDS's codons, in
        reading order (/codon_start honoured, joins followed, minus strand read backwards).
        Positions are on the sequence."""
        key = (f.index, self.editor.version, genetic_code.table_id())
        if key not in self._codon_cache:
            if len(self._codon_cache) > 300:
                self._codon_cache.clear()
            self._codon_cache[key] = cds_codons(self.editor.doc, f)
        return self._codon_cache[key]

    # ---- mouse ---------------------------------------------------------------------
    def _menu_at(self, x, y, global_pos):
        """The right-click menu for viewport point (x, y): the primer's or the feature's under it, else that of the selected
        bases, else Cut / Copy / Paste / Select All."""
        site, feature = self.primer_at(x, y), self.feature_at(x, y)
        if site is not None:
            from feature_table import primer_menu  # here: feature_table imports the tools
            primer_menu(self, self.editor, site, global_pos)
        elif feature is not None:
            from feature_table import feature_menu
            feature_menu(self, self.editor, feature, global_pos)
        elif self.editor.end > self.editor.start:  # selected bases: primers, feature, gene, copy
            from feature_table import selection_menu
            selection_menu(self, self.editor, global_pos)
        else:
            edit_menu(self, self.editor, global_pos, self.editable, self._menu_action)

    def _menu_action(self, name):
        """Cut / Copy / Paste / Select All of the menu, sent as the key that does the same: a view that refuses the key
        (one that only shows a sequence) refuses the menu's action too."""
        key = {"cut": QKeySequence.Cut, "copy": QKeySequence.Copy, "paste": QKeySequence.Paste,
               "select_all": QKeySequence.SelectAll}[name]
        combo = QKeySequence(key)[0]
        self.keyPressEvent(QKeyEvent(QKeyEvent.Type.KeyPress, combo.key(), combo.keyboardModifiers()))

    def contextMenuEvent(self, event):
        """The Menu key or Shift+F10: the menu of the selection (or Cut / Copy / Paste / Select All), shown at the cursor. A
        right click already opened its menu when the button went down (mousePressEvent)."""
        if event.reason() == QContextMenuEvent.Mouse:
            event.accept()
            return
        self._keyboard_menu()

    def _keyboard_menu(self):
        """The menu of the Menu key: as a right click on the bases at the cursor (on the first row on screen when the cursor is
        not), where no primer or feature is drawn, so it is the selection's menu or Cut / Copy / Paste / Select All. Sent as
        that click, a view that refuses right clicks (one that only shows a sequence) refuses the key too."""
        ed, sb = self.editor, self.verticalScrollBar()

        def bases_line(row):   # viewport y of the middle of a row's bases (between the two strands)
            return self.offset(row) - sb.value() + self._parts(self.real_row(row))[1] + self.lh

        row = self._screen_row(self.near(ed.start, self._first_virtual()))
        x = self.left + max(0, ed.start - self.real_row(row) * self.bases_per_line()) * self.cw
        if not 0 <= bases_line(row) < self.viewport().height():   # the cursor is off screen
            row = self.row_at(sb.value())
            if bases_line(row) < 0:   # the top row's bases are scrolled off: the next row's
                row = min(row + 1, self.rows - 1)
            x = self.left
        at = QPointF(x, bases_line(row))
        self.mousePressEvent(QMouseEvent(QEvent.Type.MouseButtonPress, at, QPointF(self.viewport().mapToGlobal(at.toPoint())),
                                         Qt.RightButton, Qt.RightButton, Qt.NoModifier))

    def _first_virtual(self):
        """Screen position of the first base on screen (positions of the copy in view, for a circular sequence)."""
        row = self.row_at(self.verticalScrollBar().value())
        return (row // self.rpc) * self.editor.length + self.real_row(row) * self.bases_per_line()

    def mousePressEvent(self, event):
        x, y = event.position().x(), event.position().y()
        if event.button() == Qt.RightButton:
            self._menu_at(x, y, event.globalPosition().toPoint())
            return
        if event.button() != Qt.LeftButton:
            return
        site = self.primer_at(x, y)
        if site is not None:
            self.editor.select_primer(site, extend=bool(event.modifiers() & Qt.ShiftModifier))
            return
        feature = self.feature_at(x, y)
        if feature is not None:
            self.editor.select_feature(feature)
            return
        if self.origin_at(x, y):  # the red origin bar: endless scrolling on / off
            self.toggle_loop()
            return
        grabbed = self.bracket_at(x, y)
        if grabbed is not None:  # drag a bracket: the other end stays put
            self.anchor = grabbed[1]
            return
        pos, _, _ = self.pos_at(x, y)
        self._drag_at = event.position()
        if event.modifiers() & Qt.ShiftModifier:
            ed = self.editor
            start = self.near(ed.start, pos)  # the copy of the selection next to the click
            self.anchor = start if pos >= start else start + ed.end - ed.start
            ed.select(self.anchor, pos)
        else:
            self.anchor = pos
            self.editor.set_cursor(pos)

    def mouseDoubleClickEvent(self, event):
        x, y = event.position().x(), event.position().y()
        site = self.primer_at(x, y)
        if site is not None:
            self.editor.edit_primer_requested.emit(self.editor.site_by_key[site].label)
            return
        feature = self.feature_at(x, y)
        if feature is not None:
            self.editor.edit_feature_requested.emit(feature)

    def mouseMoveEvent(self, event):
        x, y = event.position().x(), event.position().y()
        if self.anchor is not None and event.buttons() & Qt.LeftButton:
            # Past the top or bottom edge the selection goes on to the edge row, and the timer scrolls while it stays out.
            self._drag_at = event.position()
            h = self.viewport().height()
            if y < 0 or y >= h:
                if not self.scroll_timer.isActive():
                    self.scroll_timer.start()
            else:
                self.scroll_timer.stop()
            pos, _, _ = self.pos_at(x, min(max(y, 0.0), h - 1.0))
            self.editor.select(self.anchor, pos)
        elif not event.buttons():  # hovering: show what a bracket / the origin bar does
            on_origin = self.origin_at(x, y)
            if self._origin_tip and not on_origin:   # its tooltip goes with the pointer, it does not linger for seconds
                QToolTip.hideText()
                self._origin_tip = False
            if on_origin:
                self.viewport().setCursor(Qt.PointingHandCursor)
                QToolTip.showText(event.globalPosition().toPoint(), "Circular sequence: the origin.\nClick to "
                                  + ("stop" if self.loop else "enable") + " endless scrolling.", self.viewport(),
                                  QRect(int(x) - 6, int(y) - 6, 12, 12))
                self._origin_tip = True
            elif self.bracket_at(x, y):
                self.viewport().setCursor(Qt.SizeHorCursor)
            elif self.hover_target(x, y) is not None:   # clickable: a feature, a primer, a cut
                self.viewport().setCursor(Qt.PointingHandCursor)
            else:
                self.viewport().setCursor(Qt.IBeamCursor)

    def mouseReleaseEvent(self, event):
        if event.button() == Qt.LeftButton:
            self._drag_at = None
            self.scroll_timer.stop()
        super().mouseReleaseEvent(event)

    def _autoscroll(self):
        """While a selection is dragged past the top or bottom edge: scroll (faster the further out the pointer is) and
        take the selection on to the row now at that edge."""
        at = self._drag_at
        h = self.viewport().height()
        out = 0 if at is None or self.anchor is None else at.y() if at.y() < 0 else at.y() - h + 1 if at.y() >= h else 0
        if not out:
            self.scroll_timer.stop()
            return
        sb = self.verticalScrollBar()
        sb.setValue(sb.value() + int(math.copysign(max(self.lh / 2, min(abs(out), 4 * self.lh)), out)))
        pos, _, _ = self.pos_at(at.x(), 0.0 if out < 0 else h - 1.0)
        self.editor.select(self.anchor, pos)

    # ---- keyboard ------------------------------------------------------------------
    def keyPressEvent(self, event):
        ed, key = self.editor, event.key()
        shift = bool(event.modifiers() & Qt.ShiftModifier)
        ctrl = bool(event.modifiers() & Qt.ControlModifier)
        if event.matches(QKeySequence.ZoomIn) or (ctrl and key in (Qt.Key_Equal, Qt.Key_Plus)):
            self.set_size(self.size + 1)
        elif event.matches(QKeySequence.ZoomOut) or (ctrl and key == Qt.Key_Minus):
            self.set_size(self.size - 1)
        elif ctrl and key == Qt.Key_0:
            self.set_size(DEFAULT_SIZE)
        elif event.matches(QKeySequence.Copy):
            QGuiApplication.clipboard().setMimeData(ed.copy_mime())
        elif event.matches(QKeySequence.Cut):
            if ed.end > ed.start:
                QGuiApplication.clipboard().setMimeData(ed.copy_mime())
                ed.backspace()
        elif event.matches(QKeySequence.Paste):
            ed.paste_mime(QGuiApplication.clipboard().mimeData())
        elif event.matches(QKeySequence.SelectAll):
            ed.select(0, ed.length)
        elif key == Qt.Key_Backspace:
            ed.backspace()
        elif key == Qt.Key_Delete:
            ed.delete_forward()
        elif key in (Qt.Key_Left, Qt.Key_Right, Qt.Key_Up, Qt.Key_Down, Qt.Key_Home, Qt.Key_End, Qt.Key_PageUp, Qt.Key_PageDown):
            self._move(key, shift, whole=ctrl)
        elif key == Qt.Key_Menu or (key == Qt.Key_F10 and shift):
            self._keyboard_menu()
        elif event.text() and not event.modifiers() & (Qt.ControlModifier | Qt.AltModifier):
            ed.type_bases(event.text())
        else:
            super().keyPressEvent(event)

    def _move(self, key, extend, whole=False):
        """Arrow keys, Home / End (of the row; with Ctrl, `whole`, of the sequence), Page Up / Down (a screenful of rows);
        with Shift, `extend`, the selection's moving end goes there."""
        ed, bpl, n = self.editor, self.bases_per_line(), self.editor.length
        # The selection has a fixed end (anchor) and a moving end (head). Positions here may
        # be screen positions of the second copy (circular): compare them modulo the length.
        same = (lambda a, b: (a - b) % n == 0) if self.copies > 1 and n else (lambda a, b: a == b)
        anchor = ed.end if self.anchor is not None and same(self.anchor, ed.end) and \
            not same(self.anchor, ed.start) else ed.start
        head = ed.end if anchor == ed.start else ed.start
        if self.copies > 1 and self.anchor is not None:  # stay in the copy being worked in
            anchor, head = self.near(anchor, self.anchor), self.near(head, self.anchor)
        if not extend and ed.start != ed.end and key in (Qt.Key_Left, Qt.Key_Right):
            pos = min(anchor, head) if key == Qt.Key_Left else max(anchor, head)  # collapse the selection
        else:
            real = head % n if self.copies > 1 else head  # position in the sequence
            in_row = real % bpl                           # column in its row
            page = max(1, int(self.viewport().height() // max(1, self.row_h)) - 1) * bpl   # the rows on screen, less one
            pos = {
                Qt.Key_Left: head - 1, Qt.Key_Right: head + 1,
                Qt.Key_Up: head - bpl, Qt.Key_Down: head + bpl,
                Qt.Key_Home: head - in_row, Qt.Key_End: head - in_row + min(bpl, n - real + in_row),
                Qt.Key_PageUp: head - page, Qt.Key_PageDown: head + page,
            }[key]
            if whole and key in (Qt.Key_Home, Qt.Key_End):   # Ctrl+Home / Ctrl+End: the first / last base (of this copy)
                pos = head - real + (n if key == Qt.Key_End else 0)
        if self.copies > 1:
            pos = max(0, min(self.copies * n, pos))  # the copies on screen
        else:
            pos = max(0, min(n, pos))
        if key in (Qt.Key_PageUp, Qt.Key_PageDown):   # the page scrolls with the cursor, which stays where it was on screen
            sb = self.verticalScrollBar()
            sb.setValue(sb.value() + self.offset(self._screen_row(pos)) - self.offset(self._screen_row(head)))
        if extend:
            self.anchor = anchor
            ed.select(anchor, pos)
        else:
            self.anchor = pos
            ed.set_cursor(pos)
        self._show_virtual(pos)
