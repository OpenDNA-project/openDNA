"""The letters of the Protein view with lines drawn under each row of amino acids, one per track: the selection, the secondary structure
and the confidence of the 3D model, hydrophobicity, size, charge and property. The letters stay plain and readable; each track but the
selection is switched on or off with the buttons beside the letters (remembered in window.ini).

The text is a QTextEdit (a QPlainTextEdit ignores the room asked for under its lines): each line gets a bottom margin as tall as the
tracks shown, and paintEvent draws the tracks in it, one cell per amino acid under its letter. Hover a cell for its value.
"""
from PySide6.QtCore import QEvent, QPoint, QRectF, QSettings, Qt
from PySide6.QtGui import QColor, QPainter, QTextBlockFormat, QTextCursor
from PySide6.QtWidgets import QSizePolicy, QTextEdit, QToolButton, QToolTip, QVBoxLayout, QWidget

import structure
import theme
from primer_db import data_dir
from structure_view import AA_SCHEMES, KYTE_DOOLITTLE, SS_COLOURS, VOLUME, aa_colour

TRACKS = {   # name -> (what it shows, needs a 3D model); in the order they are drawn
    "Selection": ("The residues of the selection in the sequence", False),
    "Structure": ("Secondary structure of the 3D model: helices pink, strands gold, the rest a thin line; nothing where the model has no "
                  "residue", True),
    "Confidence": ("pLDDT of the predicted model, in AlphaFold DB's colours (blue reliable, orange not)", True),
    "Hydrophobicity": (AA_SCHEMES["Hydrophobicity"], False),
    "Size": (AA_SCHEMES["Size"] + " The taller the bar, the larger the side chain.", False),
    "Charge": (AA_SCHEMES["Charge"], False),
    "Property": (AA_SCHEMES["Property"], False),
}
SWITCHED = [name for name in TRACKS if name != "Selection"]   # the ones with a button
DEFAULT_ON = ["Structure", "Hydrophobicity", "Size"]
ROW, GAP, PAD = 5, 1, 2   # px: a track's height, between two tracks, between the letters and the first track
GROUPS = {"ILVAM": "aliphatic", "FWY": "aromatic", "KRH": "positive", "DE": "negative", "STNQ": "polar", "PG": "P or G", "C": "cysteine"}
SS_NAMES = {"H": "helix", "E": "strand", "-": "neither helix nor strand"}


def _settings():
    return QSettings(str(data_dir() / "window.ini"), QSettings.IniFormat)


def saved_tracks():
    """The switched tracks that were on last time (DEFAULT_ON the first time)."""
    value = _settings().value("protein/tracks", None)
    if value is None:
        return list(DEFAULT_ON)
    return [n for n in str(value).split(",") if n in SWITCHED]   # (a comma list: QSettings loses an empty list)


class TrackText(QTextEdit):
    """The numbered letters (see protein.layout) with the tracks under each line. Set `protein`, `model` ({index: (ss, pLDDT)} of the
    residues in the 3D model), `plddt` (the model is a prediction) and `marked` (the residues selected), then call set_tracks."""

    def __init__(self, per_line, offset_of, index_at):
        super().__init__(readOnly=True, lineWrapMode=QTextEdit.NoWrap, acceptRichText=False)
        self.document().setUndoRedoEnabled(False)
        self.per_line, self.offset_of, self.index_at = per_line, offset_of, index_at
        self.protein, self.model, self.plddt, self.marked, self.tracks = "", {}, False, set(), ["Selection"]

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
        fmt = QTextBlockFormat()
        fmt.setBottomMargin(PAD + len(self.tracks) * (ROW + GAP))
        cursor = QTextCursor(self.document())
        cursor.select(QTextCursor.Document)
        cursor.mergeBlockFormat(fmt)

    def _rows(self, block):
        """(x of each amino acid of the line `block` and its width, the y of the first track), in viewport coordinates."""
        line = block.layout().lineAt(0)
        origin = block.layout().position()
        dx, dy = -self.horizontalScrollBar().value(), -self.verticalScrollBar().value()
        width = self.fontMetrics().horizontalAdvance("M")
        first = block.blockNumber() * self.per_line
        xs = {}
        for i in range(first, min(first + self.per_line, len(self.protein))):
            x = line.cursorToX(self.offset_of(i)[1])
            xs[i] = origin.x() + (x[0] if isinstance(x, tuple) else x) + dx
        return xs, width, origin.y() + line.y() + line.height() + PAD + dy

    def paintEvent(self, event):
        super().paintEvent(event)
        if not self.protein:
            return
        p = QPainter(self.viewport())
        selection, grey = QColor(theme.selection_color()), self.palette().color(self.palette().ColorRole.PlaceholderText)
        block = self.cursorForPosition(QPoint(0, 0)).block()
        bottom = self.viewport().height()
        while block.isValid():
            xs, w, y0 = self._rows(block)
            if y0 - 40 > bottom:
                break
            for t, name in enumerate(self.tracks):
                y = y0 + t * (ROW + GAP)
                for i, x in xs.items():
                    aa, cell = self.protein[i], QRectF(x, y, w, ROW)
                    if name == "Selection":
                        colour = selection if i in self.marked else None
                    elif name in ("Structure", "Confidence"):
                        info = self.model.get(i)
                        if info is None:
                            continue
                        if name == "Confidence":
                            colour = next((c for lo, c, _ in structure.PLDDT_BANDS if info[1] >= lo), structure.PLDDT_BANDS[-1][1])
                        elif info[0] == "-":
                            p.fillRect(QRectF(x, y + ROW / 2 - 0.5, w, 1), grey)   # a loop: a thin line
                            continue
                        else:
                            colour = SS_COLOURS[info[0]]
                    elif name == "Size" and aa in VOLUME:
                        h = max(1.0, ROW * (VOLUME[aa] - 60.1) / (227.8 - 60.1))
                        cell = QRectF(x, y + ROW - h, w, h)
                        colour = aa_colour("Size", aa)
                    else:
                        colour = aa_colour(name, aa)
                    if colour:
                        p.fillRect(cell.adjusted(0, 0, -0.5, 0), QColor(colour))
            block = block.next()
        p.end()

    def track_at(self, pos):
        """(track name, amino acid index) under the viewport point `pos`, or None (over the letters, or past the tracks)."""
        cursor = self.cursorForPosition(pos)
        block = cursor.block()
        if not block.isValid() or not self.protein:
            return None
        xs, w, y0 = self._rows(block)
        t = int((pos.y() - y0) // (ROW + GAP)) if pos.y() >= y0 else -1
        hit = [i for i, x in xs.items() if x <= pos.x() < x + w]
        return (self.tracks[t], hit[0]) if 0 <= t < len(self.tracks) and hit else None

    def describe(self, name, i):
        """The value of track `name` at amino acid `i`, for its tooltip."""
        aa = self.protein[i]
        head = f"{name} · {aa}{i + 1}: "
        if name == "Selection":
            return head + ("selected" if i in self.marked else "not selected")
        if name in ("Structure", "Confidence"):
            info = self.model.get(i)
            if info is None:
                return head + "not in the model"
            return head + (SS_NAMES[info[0]] if name == "Structure" else f"pLDDT {info[1]:.0f}")
        if name == "Hydrophobicity" and aa in KYTE_DOOLITTLE:
            return head + f"{KYTE_DOOLITTLE[aa]:+.1f} (Kyte-Doolittle)"
        if name == "Size" and aa in VOLUME:
            return head + f"{VOLUME[aa]:.0f} Å³"
        if name == "Charge":
            return head + {"K": "positive", "R": "positive", "H": "positive (partly, at pH 7)", "D": "negative", "E": "negative"}.get(aa, "neutral")
        if name == "Property":
            return head + next((g for letters, g in GROUPS.items() if aa in letters), "unknown")
        return head + "–"

    def viewportEvent(self, event):
        if event.type() == QEvent.ToolTip:
            found = self.track_at(event.pos())
            if found:
                QToolTip.showText(event.globalPos(), self.describe(*found), self.viewport())
            else:
                QToolTip.hideText()
            return True
        return super().viewportEvent(event)


class TrackButtons(QWidget):
    """A column of buttons, one per switched track (SWITCHED): on or off. `changed` is called with the tracks that are on."""

    def __init__(self, changed):
        super().__init__()
        self.changed = changed
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(2)
        on = saved_tracks()
        self.buttons = {}
        for name in SWITCHED:
            b = QToolButton(text=name, checkable=True, checked=name in on, toolTip=TRACKS[name][0])
            b.setToolButtonStyle(Qt.ToolButtonTextOnly)
            b.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Fixed)
            b.toggled.connect(self._toggled)
            layout.addWidget(b)
            self.buttons[name] = b
        layout.addStretch()

    def on(self):
        return [n for n, b in self.buttons.items() if b.isChecked()]

    def set_model(self, has_model):
        """The tracks of the 3D model can be switched on only when there is one."""
        for name, b in self.buttons.items():
            if TRACKS[name][1]:
                b.setEnabled(has_model)
                b.setToolTip(TRACKS[name][0] if has_model else f"{TRACKS[name][0]}. Get a 3D structure first (3D structure…).")

    def _toggled(self):
        _settings().setValue("protein/tracks", ",".join(self.on()))
        self.changed()
