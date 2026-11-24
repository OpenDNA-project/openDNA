"""What to take from a GFF3 / GTF annotation track: which feature types (the dialog lists each with its number of lines), and the
extra fields the file has (a colour per line, a note shown when hovering). `read_color` reads the colour of a line in any of the
ways files write it; one it cannot read is drawn in the track's colour and counted in a warning, never dropped silently."""
import re

from PySide6.QtCore import Qt
from PySide6.QtGui import QColor
from PySide6.QtWidgets import QCheckBox, QDialog, QDialogButtonBox, QGroupBox, QLabel, QScrollArea, QVBoxLayout, QWidget

from uikit import plural

ALL = "all"   # the types (or fields) of a project saved before the choice existed: everything, nothing asked
FIELDS = {"color": "Colour", "note": "Note (shown when hovering)"}

# Artemis numbers its colours 0-17 (its manual, "colour" qualifier).
ARTEMIS = [(255, 255, 255), (100, 100, 100), (255, 0, 0), (0, 255, 0), (0, 0, 255), (0, 255, 255), (255, 0, 255), (255, 255, 0),
           (152, 251, 152), (135, 206, 250), (255, 165, 0), (200, 150, 100), (255, 200, 200), (170, 170, 170), (0, 0, 0),
           (255, 63, 63), (255, 127, 127), (255, 191, 191)]


def read_color(raw):
    """"#rrggbb" for a colour as a file writes it, or None when it cannot be read. Read: #RRGGBB, #RGB, RRGGBB, 0xRRGGBB,
    "r,g,b", "r g b" (Artemis), rgb(r, g, b), a colour name (red, darkgreen…), an Artemis colour number 0-17."""
    text = raw.strip().lower()
    numbers = re.fullmatch(r"(?:rgb\s*\()?\s*(\d{1,3})\s*[, ]\s*(\d{1,3})\s*[, ]\s*(\d{1,3})\s*\)?", text)
    if numbers:
        rgb = [int(n) for n in numbers.groups()]
        return "#%02x%02x%02x" % tuple(rgb) if all(n <= 255 for n in rgb) else None
    if text.isdigit():
        return "#%02x%02x%02x" % ARTEMIS[int(text)] if int(text) < len(ARTEMIS) else None
    hexa = re.fullmatch(r"(?:#|0x)?([0-9a-f]{6}|[0-9a-f]{3})", text)
    text = "#" + hexa.group(1) if hexa else text   # (Qt reads #RRGGBB, #RGB and the colour names)
    return QColor(text).name() if QColor.isValidColorName(text) else None


def colors_of(track):
    """({colour as written: "#rrggbb" or None}, a sentence about the ones not read or "") for a track's colours."""
    found = track.colors()
    read = {raw: read_color(raw) for raw, _ in found}
    bad = [(raw, n) for raw, n in found if read[raw] is None]
    if not bad:
        return read, ""
    lines = sum(n for _, n in bad)
    return read, (f"{plural(lines, 'line')} with a colour that cannot be read (e.g. '{bad[0][0]}') "
                  "are drawn in the track's colour.")


class TrackImportDialog(QDialog):
    """Import a GFF/GTF track: a box per feature type with its number of lines (ticked, except the types whose lines only cover
    their whole sequence, like a `region` line that names the chromosome), then the extra fields found."""

    def __init__(self, parent, name, types, fields):
        super().__init__(parent)
        self.setWindowTitle(f"Import {name}")
        layout = QVBoxLayout(self)
        layout.addWidget(QLabel("Feature types to import:"))
        body = QWidget()
        column = QVBoxLayout(body)
        column.setContentsMargins(0, 0, 0, 0)
        self.types = {}
        for kind, count, whole in types:
            box = QCheckBox(f"{kind} ({count:,})", checked=not whole)
            if whole:
                box.setToolTip("Every line of this type covers its whole sequence (it names the sequence): not ticked")
            self.types[kind] = box
            column.addWidget(box)
        column.addStretch()
        scroll = QScrollArea(widgetResizable=True)   # a genome annotation can have dozens of types
        scroll.setWidget(body)
        scroll.setMaximumHeight(320)
        layout.addWidget(scroll)
        self.fields = {}
        total = sum(count for _, count, _ in types)
        if fields:
            group = QGroupBox("Also import")
            inner = QVBoxLayout(group)
            for key, count in fields.items():
                self.fields[key] = QCheckBox(f"{FIELDS[key]} ({count:,} of {total:,} lines)", checked=True)
                inner.addWidget(self.fields[key])
            layout.addWidget(group)
        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel, Qt.Horizontal)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def choice(self):
        return ([k for k, b in self.types.items() if b.isChecked()], [k for k, b in self.fields.items() if b.isChecked()])


def detected_fields(track):
    """{field: number of lines that have it} for the extra fields the file has (none for a BAM)."""
    color, note = track.field_counts()
    return {k: n for k, n in (("color", color), ("note", note)) if n}
