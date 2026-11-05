"""Add / edit a feature: name, type, location, color, other qualifiers."""
import re

import validate

from PySide6.QtWidgets import (
    QColorDialog, QComboBox, QDialog, QDialogButtonBox, QFormLayout, QLineEdit,
    QMessageBox, QPlainTextEdit, QPushButton,
)

from editor import COLOR_KEYS, DIRECTION_KEY, DIRECTIONLESS_KINDS, KIND_COLORS, arrow, feature_color, site_location, type_color

DIRECTIONS = ["Forward →", "Reverse ←", "None (no arrow)"]


QUALIFIER_KEY = re.compile(r"[A-Za-z0-9_'-]+")   # GenBank qualifier names: letters, digits, _ ' -
MAX_NAME = 200
LOCATION_HELP = "GenBank syntax, 1-based: 10..200, complement(10..200), join(1..50,80..120)"
KINDS = sorted(set(KIND_COLORS) | {"misc_feature", "enhancer", "exon", "intron", "LTR", "sig_peptide"})


class FeatureDialog(QDialog):
    def __init__(self, parent, editor, index=None):
        super().__init__(parent)
        self.editor, self.index = editor, index
        f = editor.by_index.get(index) if index is not None else None
        self.setWindowTitle("Edit Feature" if f else "New Feature")
        self.setMinimumWidth(460)

        self.name = QLineEdit(f.label if f else "")
        self.name.setMaxLength(max(MAX_NAME, len(self.name.text())))   # never shorter than a label the file already has
        self.kind = QComboBox(editable=True)
        self.kind.addItems(KINDS)
        self.kind.setCurrentText(f.kind if f else "misc_feature")
        self.location = QLineEdit(f.location if f else site_location(editor.start, editor.end, 1, editor.length))
        self.location.setToolTip(LOCATION_HELP)
        self.direction = QComboBox()
        self.direction.addItems(DIRECTIONS)
        self.direction.setToolTip("None: drawn without an arrow (the default of misc_feature); saved as /directional")
        way = arrow(f) if f else (0 if self.kind.currentText() in DIRECTIONLESS_KINDS else 1)
        self.direction.setCurrentIndex({1: 0, -1: 1, 0: 2}[way])
        self.direction_set = False   # until chosen by hand, a new feature follows its type's default
        self.direction.activated.connect(self._direction_chosen)
        self.color = feature_color(f) if f else type_color(self.kind.currentText())
        # Write color qualifiers only if the feature had them or the user picks a color.
        self.color_set = f is not None and any(k in COLOR_KEYS for k, _ in f.qualifiers)
        self.color_button = QPushButton(clicked=self._pick_color)
        self._show_color()
        self.kind.currentTextChanged.connect(self._kind_changed)
        # Every other qualifier as "key=value" lines (label and colors have their own fields).
        hidden = {"label", DIRECTION_KEY, *COLOR_KEYS}
        rest = [(k, v) for k, v in (f.qualifiers if f else []) if k not in hidden]
        self.qualifiers = QPlainTextEdit("\n".join(k if v is None else f"{k}={v}" for k, v in rest))
        self.qualifiers.setPlaceholderText("one per line, e.g.\nnote=my note\ngene=lacZ")
        self.had_label = f is not None and any(k == "label" for k, _ in f.qualifiers)
        self.original_name = self.name.text()

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        form = QFormLayout(self)
        form.addRow("Name", self.name)
        form.addRow("Type", self.kind)
        form.addRow("Location", self.location)
        form.addRow("Direction", self.direction)
        form.addRow("Color", self.color_button)
        form.addRow("Qualifiers", self.qualifiers)
        form.addRow(buttons)

    def _direction_chosen(self, index):
        self.direction_set = True
        text = self.location.text().strip()
        reverse = index == 1
        if reverse and not text.startswith("complement("):
            self.location.setText(f"complement({text})")
        elif not reverse and text.startswith("complement(") and text.endswith(")"):
            self.location.setText(text[len("complement("):-1])

    def _kind_changed(self, kind):
        """Until a colour (or direction) is picked by hand, it follows the type's standard one."""
        if not self.color_set:
            self.color = type_color(kind)
            self._show_color()
        if self.index is None and not self.direction_set and self.direction.currentIndex() != 1:
            self.direction.setCurrentIndex(2 if kind in DIRECTIONLESS_KINDS else 0)

    def _pick_color(self):
        color = QColorDialog.getColor(self.color, self, "Feature color")
        if color.isValid():
            self.color, self.color_set = color, True
            self._show_color()

    def _show_color(self):
        self.color_button.setText(self.color.name())
        self.color_button.setStyleSheet(f"background: {self.color.name()}")

    def qualifier_list(self):
        out = []
        name = self.name.text().strip()
        # Only write /label when the user named it, so features named by /gene stay as they were.
        if name and (self.had_label or name != self.original_name):
            out.append(("label", name))
        for line in self.qualifiers.toPlainText().splitlines():
            line = line.strip().lstrip("/")
            if not line:
                continue
            key, sep, value = line.partition("=")
            key, value = key.strip(), value.strip()
            if not QUALIFIER_KEY.fullmatch(key):
                raise ValueError(f"'{line[:40]}': a qualifier is written key=value, and the key is one word (letters, digits, _), "
                                 "for example note=my note")
            if len(value) >= 2 and value[0] == value[-1] == '"':
                value = value[1:-1]
            if any(not c.isprintable() for c in value):
                raise ValueError(f"'{line[:40]}': the value has a tab or control character; use plain text (one qualifier per line).")
            out.append((key, value if sep else None))
        if self.color_set:
            out += [(key, self.color.name()) for key in COLOR_KEYS]
        # Only what differs from the type's default is written.
        directional = self.direction.currentIndex() != 2
        if directional == (self.kind.currentText().strip() in DIRECTIONLESS_KINDS):
            out.append((DIRECTION_KEY, "yes" if directional else "no"))
        return out

    def accept(self):
        problem = validate.text_problem(self.name.text().strip(), "name", max(MAX_NAME, len(self.original_name or "")))
        if problem:
            QMessageBox.warning(self, "Invalid feature", problem)
            self.name.setFocus()
            return
        location = self.location.text().strip()
        if not location:
            QMessageBox.warning(self, "Invalid feature", f"Give the location. {LOCATION_HELP}")
            self.location.setFocus()
            return
        if not re.search(r"\d", location) or not re.fullmatch(r"[\w\s().,<>^:-]+", location):   # no numbers: not a location at all
            QMessageBox.warning(self, "Invalid feature", f"Cannot read the location '{location}'. {LOCATION_HELP}")
            self.location.setFocus()
            return
        try:
            self.editor.set_feature(self.index, self.kind.currentText().strip(), location, self.qualifier_list())
        except (OSError, ValueError) as e:
            text = str(e)
            if "Syntax error" in text or text.startswith("bad location"):   # the parser's own words mean nothing to a user
                text = f"Cannot read the location '{location}'. {LOCATION_HELP}"
            QMessageBox.warning(self, "Invalid feature", text)
            return  # keep the dialog open so the user can fix it
        super().accept()


def edit_feature(parent, editor, index=None):
    """Open the dialog; `index` None = new feature over the current selection."""
    if index is None and editor.end <= editor.start:
        editor.error.emit("Select some bases first, then add a feature.")
        return False
    if index is not None and index not in editor.by_index:   # a tool's temporary feature (Editor.set_overlay)
        return False
    return FeatureDialog(parent, editor, index).exec() == QDialog.Accepted
