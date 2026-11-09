"""The letters of the Protein view, with tracks under each row of amino acids and a legend beside them.

The letters stay plain. Above each row, small numbers mark the end of every group of ten (10, 20, ...); the row's first position is
in the number column. A selection gets the DNA view's brackets: "[" before its first residue with that position above it, "]" after
its last with that position below the tracks. Under the letters, one track per row, each with its label at the start of the row:
  Selection        a bar under the selected residues (always shown)
  Structure        the 3D model's secondary structure: a pink wave for a helix, a gold arrow for a strand, a thin line for the rest
  Confidence       the model's pLDDT as a bar (taller is surer) in AlphaFold DB's colours, for a prediction
  Hydrophobicity   a bar up (red) for a hydrophobic residue, down (blue) for a hydrophilic one, by the Kyte-Doolittle scale
  Size             a bar as tall as the side chain is large (its volume)
  Charge           "+" for K and R (pale for H, partly charged at pH 7), "−" for D and E
  Property         Zappo colour blocks: aliphatic, aromatic, positive, negative, polar, P and G, C
The legend (TrackLegend) shows each track with samples drawn by the same code (draw_cell), and a check box to switch it (remembered in
window.ini). Rows hold as many groups of ten as the width allows (wanted_per_line). Hover a cell for its value.

The text is a QTextEdit (a QPlainTextEdit ignores the room asked for above and under its lines): each line gets margins as tall as the
numbers and tracks need, and paintEvent draws them there.
"""
import math

from PySide6.QtCore import QEvent, QPoint, QPointF, QRectF, QSettings, Qt
from PySide6.QtGui import QColor, QFont, QFontMetricsF, QPainter, QPainterPath, QPalette, QPen, QPolygonF, QTextBlockFormat, QTextCursor
from PySide6.QtWidgets import (
    QCheckBox, QGridLayout, QLabel, QScrollArea, QSizePolicy, QTextEdit, QToolTip, QVBoxLayout, QWidget,
)

import structure
import theme
from primer_db import data_dir
from structure_view import KYTE_DOOLITTLE, PROPERTY, SS_COLOURS, VOLUME

TRACKS = {   # name -> (label at the start of its row, what it shows, needs a 3D model); in the order they are drawn
    "Selection": ("select", "The residues of the selection in the sequence", False),
    "Structure": ("struct", "Secondary structure of the 3D model; nothing where the model has no residue", True),
    "Confidence": ("pLDDT", "Confidence of the predicted model (pLDDT): the taller the bar, the surer, in AlphaFold DB's colours", True),
    "Hydrophobicity": ("hydro", "Kyte-Doolittle hydropathy: up and red is hydrophobic, down and blue hydrophilic", False),
    "Size": ("size", "Volume of the side chain (Zamyatnin): the taller the bar, the larger", False),
    "Charge": ("charge", "Charge at pH 7: + positive (K, R; pale for H), − negative (D, E)", False),
    "Property": ("prop", "Kind of side chain, in Jalview's Zappo colours", False),
}
SWITCHED = [name for name in TRACKS if name != "Selection"]   # the ones with a check box
DEFAULT_ON = ["Structure", "Hydrophobicity", "Charge"]
ROW, GAP, PAD = 10, 2, 3   # px: a track's height, between two tracks, between the letters and the first track
SMALL = 9                  # px: the numbers and the row labels
HYDROPHOBIC, HYDROPHILIC, LARGE = "#E0483F", "#3F7FD8", "#A0581C"
POSITIVE, PARTLY, NEGATIVE = "#2F63D8", "#8FAAF0", "#D83A3A"
GROUPS = [("ILVAM", "aliphatic"), ("FWY", "aromatic"), ("KRH", "positive"), ("DE", "negative"), ("STNQ", "polar"), ("PG", "P and G"),
          ("C", "cysteine")]
SS_NAMES = {"H": "helix", "E": "strand", "-": "neither helix nor strand"}


def _settings():
    return QSettings(str(data_dir() / "window.ini"), QSettings.IniFormat)


def saved_tracks():
    """The switched tracks that were on last time (DEFAULT_ON the first time)."""
    value = _settings().value("protein/tracks", None)
    if value is None:
        return list(DEFAULT_ON)
    return [n for n in str(value).split(",") if n in SWITCHED]   # (a comma list: QSettings loses an empty list)


def small_font(base):
    font = QFont(base)
    font.setPixelSize(SMALL)
    return font


def draw_cell(p, name, r, aa="", info=None, selected=False, strand_ends=False, grey=None):
    """Draw track `name` for one amino acid in the cell `r` (QRectF): `aa` its letter, `info` (ss, pLDDT) of its residue in the model
    (None: not in it), `selected`, `strand_ends` (the last residue of a strand: the arrow head). Shared by the tracks and the legend."""
    x, w, top, mid, bottom = r.x(), r.width(), r.top(), r.center().y(), r.bottom()
    grey = grey or QColor("#909090")
    if name == "Selection":
        if selected:
            p.fillRect(QRectF(x, mid - 2, w, 4), theme.selection_color())
    elif name == "Structure":
        if info is None:
            return
        ss = info[0]
        if ss == "H":   # a wave, one turn per residue (it joins up from cell to cell)
            path = QPainterPath(QPointF(x, mid))
            for k in range(1, 9):
                path.lineTo(x + w * k / 8, mid - (ROW * 0.38) * math.sin(2 * math.pi * k / 8))
            p.setPen(QPen(QColor(SS_COLOURS["H"]), 1.8))
            p.setBrush(Qt.NoBrush)
            p.drawPath(path)
        elif ss == "E":
            p.setPen(Qt.NoPen)
            p.setBrush(QColor(SS_COLOURS["E"]))
            if strand_ends:
                p.drawPolygon(QPolygonF([QPointF(x, top), QPointF(x + w, mid), QPointF(x, bottom)]))
            else:
                p.drawRect(QRectF(x, mid - ROW * 0.22, w + 0.5, ROW * 0.44))
        else:
            p.fillRect(QRectF(x, mid - 0.5, w, 1), grey)
    elif name == "Confidence":
        if info is None:
            return
        colour = next((c for lo, c, _ in structure.PLDDT_BANDS if info[1] >= lo), structure.PLDDT_BANDS[-1][1])
        h = max(1.5, ROW * min(info[1], 100) / 100)
        p.fillRect(QRectF(x, bottom - h, w - 0.5, h), QColor(colour))
    elif name == "Hydrophobicity":
        p.fillRect(QRectF(x, mid - 0.25, w, 0.5), grey)   # the zero line
        if aa in KYTE_DOOLITTLE:
            v = KYTE_DOOLITTLE[aa]
            h = max(1.0, abs(v) / 4.5 * ROW / 2)
            p.fillRect(QRectF(x, mid - h if v > 0 else mid, w - 0.5, h), QColor(HYDROPHOBIC if v > 0 else HYDROPHILIC))
    elif name == "Size":
        if aa in VOLUME:
            h = max(1.0, ROW * VOLUME[aa] / 227.8)
            p.fillRect(QRectF(x, bottom - h, w - 0.5, h), QColor(LARGE))
    elif name == "Charge":
        sign = {"K": ("+", POSITIVE), "R": ("+", POSITIVE), "H": ("+", PARTLY), "D": ("−", NEGATIVE), "E": ("−", NEGATIVE)}.get(aa)
        if sign:
            font = QFont(p.font())
            font.setPixelSize(ROW + 3)
            font.setBold(True)
            p.setFont(font)
            p.setPen(QColor(sign[1]))
            p.drawText(QRectF(x - 2, top - 2, w + 4, ROW + 4), Qt.AlignCenter, sign[0])
    elif name == "Property":
        if aa in PROPERTY:
            p.fillRect(QRectF(x, top + 1, w - 0.5, ROW - 2), QColor(PROPERTY[aa]))


def describe(name, aa, i, info=None, selected=False):
    """The value of track `name` for amino acid `aa` (number i + 1), for its tooltip."""
    head = f"{name} · {aa}{i + 1}: "
    if name == "Selection":
        return head + ("selected" if selected else "not selected")
    if name in ("Structure", "Confidence"):
        if info is None:
            return head + "not in the model"
        return head + (SS_NAMES[info[0]] if name == "Structure" else f"pLDDT {info[1]:.0f}")
    if name == "Hydrophobicity" and aa in KYTE_DOOLITTLE:
        return head + f"{KYTE_DOOLITTLE[aa]:+.1f} (Kyte-Doolittle; above 0 hydrophobic)"
    if name == "Size" and aa in VOLUME:
        return head + f"{VOLUME[aa]:.0f} Å³"
    if name == "Charge":
        return head + {"K": "positive", "R": "positive", "H": "partly positive at pH 7", "D": "negative", "E": "negative"}.get(aa, "neutral")
    if name == "Property":
        return head + next((g for letters, g in GROUPS if aa in letters), "unknown")
    return head + "–"


class TrackText(QTextEdit):
    """The numbered letters (protein.layout, `per_line` amino acids a row) with numbers, brackets and tracks. Set `protein`, `model`
    ({index: (ss, pLDDT)} of the residues in the 3D model) and `marked` (the residues selected), then call set_tracks. `relayout`
    is called when the width calls for another number of groups of ten a row (see wanted_per_line)."""

    def __init__(self, offset_of, index_at, relayout, per_line=50, block=10, margin=8):
        super().__init__(readOnly=True, lineWrapMode=QTextEdit.NoWrap, acceptRichText=False)
        self.document().setUndoRedoEnabled(False)
        self.offset_of, self.index_at, self.relayout = offset_of, index_at, relayout
        self.per_line, self.block, self.margin = per_line, block, margin
        self.protein, self.model, self.marked, self.tracks = "", {}, set(), ["Selection"]

    # ---- the rows ---------------------------------------------------------------------------------------------------------------------
    def wanted_per_line(self):
        """As many groups of ten as fit in the width (one at least). The width counts a vertical scroll bar even when there is none:
        else a row more or less could show or hide the bar and change the count again."""
        room = self.width() - 2 * self.frameWidth() - self.verticalScrollBar().sizeHint().width() - 2 * self.document().documentMargin() - 6
        chars = room / max(1.0, QFontMetricsF(self.font()).horizontalAdvance("M"))
        return max(1, int((chars - self.margin + 1) // (self.block + 1))) * self.block

    def resizeEvent(self, event):
        super().resizeEvent(event)
        if self.wanted_per_line() != self.per_line:
            self.relayout()

    def setPlainText(self, text):
        super().setPlainText(text)
        self._make_room()   # (a new text has plain lines)

    def set_tracks(self, names):
        """Draw the tracks `names` (in TRACKS' order), making room for them under each line."""
        names = [n for n in TRACKS if n in names]
        if names != self.tracks:
            self.tracks = names
            self._make_room()
        self.viewport().update()

    def _make_room(self):
        # All the room is under the lines: QTextDocument collapses a line's top margin into the one above (and ignores it on the first,
        # which gets the document margin instead). Under each line: the tracks, the last position of a selection, then the numbers of
        # the next line's groups.
        self.document().setDocumentMargin(SMALL + 6)
        fmt = QTextBlockFormat()
        fmt.setBottomMargin(PAD + len(self.tracks) * (ROW + GAP) + 2 * SMALL + 8)
        # the last line's bottom margin is dropped too: the document's frame keeps the room under it, so its tracks can be scrolled to
        frame = self.document().rootFrame()
        frame_format = frame.frameFormat()
        frame_format.setBottomMargin(fmt.bottomMargin() + SMALL + 6)
        frame.setFrameFormat(frame_format)
        cursor = QTextCursor(self.document())
        cursor.select(QTextCursor.Document)
        cursor.mergeBlockFormat(fmt)

    def _rows(self, block):
        """(x of each amino acid of the line `block`, the width of a letter, the y of the first track), in viewport coordinates."""
        line = block.layout().lineAt(0)
        origin = block.layout().position()
        dx, dy = -self.horizontalScrollBar().value(), -self.verticalScrollBar().value()
        width = QFontMetricsF(self.font()).horizontalAdvance("M")
        first = block.blockNumber() * self.per_line
        xs = {}
        for i in range(first, min(first + self.per_line, len(self.protein))):
            x = line.cursorToX(self.offset_of(i, self.per_line)[1])
            xs[i] = origin.x() + (x[0] if isinstance(x, tuple) else x) + dx
        return xs, width, origin.y() + line.y() + line.height() + PAD + dy

    def _text_top(self, block):
        line = block.layout().lineAt(0)
        return block.layout().position().y() + line.y() - self.verticalScrollBar().value()

    # ---- drawing ------------------------------------------------------------------------------------------------------------------------
    def paintEvent(self, event):
        super().paintEvent(event)
        if not self.protein:
            return
        p = QPainter(self.viewport())
        p.setRenderHint(QPainter.Antialiasing)
        pal = self.palette()
        grey, base = pal.color(QPalette.PlaceholderText), pal.color(QPalette.Base)
        small = small_font(self.font())
        fm = QFontMetricsF(small)
        first_marked, last_marked = (min(self.marked), max(self.marked)) if self.marked else (None, None)
        block = self.cursorForPosition(QPoint(0, 0)).block()
        bottom = self.viewport().height()
        while block.isValid():
            xs, w, y0 = self._rows(block)
            top = self._text_top(block)
            if top - SMALL > bottom:
                break
            if xs:
                tracks_end = y0 + len(self.tracks) * (ROW + GAP)
                left = min(xs.values())
                # the label of each track at the start of its row
                p.setFont(small)
                p.setPen(grey)
                for t, name in enumerate(self.tracks):
                    p.drawText(QRectF(0, y0 + t * (ROW + GAP) - 2, left - 6, ROW + 4), Qt.AlignRight | Qt.AlignVCenter, TRACKS[name][0])
                # the end of each group of ten, numbered above its last letter
                for i, x in xs.items():
                    if (i + 1) % self.block == 0:
                        text = f"{i + 1:,}"
                        p.drawText(QPointF(x + w - fm.horizontalAdvance(text), top - 3), text)
                # the tracks
                for t, name in enumerate(self.tracks):
                    y = y0 + t * (ROW + GAP)
                    for i, x in xs.items():
                        info = self.model.get(i)
                        ends = bool(info) and info[0] == "E" and (self.model.get(i + 1) or ("-",))[0] != "E"
                        draw_cell(p, name, QRectF(x, y, w, ROW), self.protein[i], info, i in self.marked, ends, grey)
                # the brackets of the selection, with its first and last positions
                p.setFont(small)
                if first_marked in xs:
                    self._bracket(p, xs[first_marked] - 1, top, tracks_end, 1)
                    self._number(p, fm, base, xs[first_marked] - 1, top - 3 - fm.height(), first_marked + 1, 1)
                if last_marked in xs:
                    self._bracket(p, xs[last_marked] + w + 1, top, tracks_end, -1)
                    self._number(p, fm, base, xs[last_marked] + w + 1, tracks_end + 1, last_marked + 1, -1)
            block = block.next()
        p.end()

    def _bracket(self, p, x, top, bottom, side):
        """'[' (side 1) or ']' (-1) from the letters down through the tracks, as in the DNA view."""
        p.setPen(QPen(theme.selection_color(), 2))
        arm = side * 4.0
        p.drawPolyline(QPolygonF([QPointF(x + arm, top), QPointF(x, top), QPointF(x, bottom), QPointF(x + arm, bottom)]))

    def _number(self, p, fm, base, x, top, number, side):
        """A position at a bracket: the text starts at a '[' and ends at a ']', on the background (readable over a group number)."""
        text = f"{number:,}"
        w = fm.horizontalAdvance(text)
        left = x if side > 0 else x - w
        p.fillRect(QRectF(left - 1, top, w + 2, fm.height()), base)
        p.setPen(theme.selection_color())
        p.drawText(QPointF(left, top + fm.ascent()), text)

    # ---- hovering ------------------------------------------------------------------------------------------------------------------------
    def track_at(self, pos):
        """(track name, amino acid index) under the viewport point `pos`, or None (over the letters, or past the tracks)."""
        block = self.cursorForPosition(pos).block()
        if not block.isValid() or not self.protein:
            return None
        xs, w, y0 = self._rows(block)
        t = int((pos.y() - y0) // (ROW + GAP)) if pos.y() >= y0 else -1
        hit = [i for i, x in xs.items() if x <= pos.x() < x + w]
        return (self.tracks[t], hit[0]) if 0 <= t < len(self.tracks) and hit else None

    def describe(self, name, i):
        return describe(name, self.protein[i], i, self.model.get(i), i in self.marked)

    def viewportEvent(self, event):
        if event.type() == QEvent.ToolTip:
            found = self.track_at(event.pos())
            if found:
                QToolTip.showText(event.globalPos(), self.describe(*found), self.viewport())
            else:
                QToolTip.hideText()
            return True
        return super().viewportEvent(event)


# ---- the legend -------------------------------------------------------------------------------------------------------------------------

KEYS = {   # track -> [(sample: [(aa, info, selected, strand_ends)] or "brackets", caption)]
    "Selection": [([("", None, True, False)] * 3, "selected"), ("brackets", "first and last position")],
    "Structure": [([("", ("H", 0), False, False)] * 3, "helix"), ([("", ("E", 0), False, False)] * 2 + [("", ("E", 0), False, True)], "strand"),
                  ([("", ("-", 0), False, False)] * 3, "loop, turn")],
    "Confidence": [([("", ("-", b), False, False)] * 2, text) for b, text in
                   ((95, "very high, > 90"), (80, "confident, 70–90"), (60, "low, 50–70"), (30, "very low, < 50"))],
    "Hydrophobicity": [([(a, None, False, False) for a in "IVL"], "hydrophobic (I V L)"),
                       ([(a, None, False, False) for a in "RKD"], "hydrophilic (R K D)")],
    "Size": [([(a, None, False, False) for a in "GAS"], "small (G A S)"), ([(a, None, False, False) for a in "FYW"], "large (F Y W)")],
    "Charge": [([("K", None, False, False)], "positive (K, R)"), ([("H", None, False, False)], "partly positive (H)"),
               ([("D", None, False, False)], "negative (D, E)")],
    "Property": [([(letters[0], None, False, False)], f"{name} ({' '.join(letters)})") for letters, name in GROUPS],
}


class Sample(QWidget):
    """A few cells of a track, drawn as in the rows (draw_cell): the colour and the symbol of one entry of the legend."""

    def __init__(self, track, cells):
        super().__init__()
        self.track, self.cells = track, cells
        n = 2 if cells == "brackets" else len(cells)
        self.setFixedSize(int(n * 8 + 4), ROW + 4)

    def paintEvent(self, event):
        p = QPainter(self)
        p.setRenderHint(QPainter.Antialiasing)
        p.setOpacity(1.0 if self.isEnabled() else 0.35)
        grey = self.palette().color(QPalette.PlaceholderText)
        if self.cells == "brackets":
            p.setPen(QPen(theme.selection_color(), 2))
            for x, side in ((2, 1), (self.width() - 2, -1)):
                p.drawPolyline(QPolygonF([QPointF(x + side * 4, 1), QPointF(x, 1), QPointF(x, self.height() - 1),
                                          QPointF(x + side * 4, self.height() - 1)]))
        else:
            for k, (aa, info, selected, ends) in enumerate(self.cells):
                draw_cell(p, self.track, QRectF(2 + k * 8, 2, 8, ROW), aa, info, selected, ends, grey)
        p.end()


class TrackLegend(QScrollArea):
    """What each track shows: its name (a check box switches it on and off, remembered; the selection is always shown), its label at
    the start of its row, and samples of its colours and symbols. `changed` is called when a track is switched."""

    def __init__(self, changed):
        super().__init__(widgetResizable=True, frameShape=QScrollArea.NoFrame)
        self.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Expanding)
        self.changed = changed
        inner = QWidget()
        column = QVBoxLayout(inner)
        column.setContentsMargins(4, 0, 4, 0)
        column.setSpacing(1)
        on = saved_tracks()
        self.buttons, self.keys = {}, {}
        for name, (label, about, _) in TRACKS.items():
            if name == "Selection":
                head = QLabel(f"<b>{name}</b>")
            else:
                head = QCheckBox(name, checked=name in on)
                head.toggled.connect(self._toggled)
                self.buttons[name] = head
            head.setToolTip(about)
            column.addWidget(head)
            key = QWidget()
            grid = QGridLayout(key)
            grid.setContentsMargins(18, 0, 0, 4)
            grid.setHorizontalSpacing(4)
            grid.setVerticalSpacing(0)
            caption = QLabel(f"row “{label}”")
            caption.setFont(small_font(caption.font()))
            grid.addWidget(caption, 0, 0, 1, 2)
            for n, (cells, text) in enumerate(KEYS[name], start=1):
                words = QLabel(text)
                words.setFont(small_font(words.font()))
                grid.addWidget(Sample(name, cells), n, 0)
                grid.addWidget(words, n, 1)
            grid.setColumnStretch(1, 1)
            key.setToolTip(about)
            column.addWidget(key)
            self.keys[name] = key
        column.addStretch()
        self.setWidget(inner)
        self.setFixedWidth(inner.sizeHint().width() + self.verticalScrollBar().sizeHint().width())
        self._grey()

    def on(self):
        return [n for n, b in self.buttons.items() if b.isChecked()]

    def set_model(self, has_model, plddt=False):
        """The tracks of the 3D model can be switched on only when there is one (and Confidence only for a prediction)."""
        for name, b in self.buttons.items():
            if TRACKS[name][2]:
                usable = has_model and (plddt or name != "Confidence")
                b.setEnabled(usable)
                why = "Get a 3D structure first (3D structure…)." if not has_model else "An experimental structure has no pLDDT."
                b.setToolTip(TRACKS[name][1] if usable else f"{TRACKS[name][1]}. {why}")
        self._grey()

    def _grey(self):
        for name, key in self.keys.items():
            b = self.buttons.get(name)
            key.setEnabled(b is None or (b.isChecked() and b.isEnabled()))

    def _toggled(self):
        _settings().setValue("protein/tracks", ",".join(self.on()))
        self._grey()
        self.changed()
