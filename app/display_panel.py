"""Display filter (side panel, View → Display Filter): what the sequence view, map and genome
view draw. One switch per feature type found in the sequence (with its count), primers,
unique cutters, CDS translations, the 3-frame translation and transcription units. Only the
drawing changes: features stay in the file, the feature table and Find. The broad switches are also toggles at the right
end of the tab row (QuickFilters).
"""

from types import SimpleNamespace

from PySide6.QtCore import QEvent, QPointF, QSize, Qt, Signal
from PySide6.QtGui import QColor, QIcon, QPainter, QPen
from PySide6.QtWidgets import (
    QApplication, QCheckBox, QComboBox, QFormLayout, QFrame, QHBoxLayout, QLabel, QProxyStyle, QPushButton, QScrollArea,
    QSpinBox, QStyle, QToolButton, QVBoxLayout, QWidget,
)

import genetic_code
import icons
from uikit import plural
import theme
import tool_settings
from enzymes import PRESETS, EnzymeFilter
from sequence_view import primer_color


class PrimerSwatch(QWidget):
    """A small sample of how a primer is drawn: the same colour (primer_color) and pen as the sites in the views, a line with a
    half arrow at its 3' end. `provisional`: the temporary style (designed here, or only in the primer database)."""

    def __init__(self, provisional):
        super().__init__()
        self.provisional = provisional
        self.setFixedSize(QSize(34, 14))

    @property
    def color(self):   # read when drawn: the primer colour has a dark-theme twin
        return QColor(primer_color(SimpleNamespace(provisional=self.provisional, temporary=self.provisional, mismatches=0)))

    def paintEvent(self, event):
        p = QPainter(self)
        p.setRenderHint(QPainter.Antialiasing)
        p.setPen(QPen(self.color, 2))
        y = self.height() / 2 + 1
        p.drawLine(QPointF(3, y), QPointF(self.width() - 3, y))
        p.drawLine(QPointF(self.width() - 3, y), QPointF(self.width() - 8, y - 4))   # the half arrow


class QuickFilters(QWidget):
    """The right end of the tab row: one toggle per broad category of the display filter (icons.category_icon), clear
    when shown, faded when hidden. Same switches as the side panel, which keeps the finer ones (each primer source, each
    annotation type); right-click opens it."""
    def __init__(self, window, editor):
        super().__init__()
        self.window, self.editor = window, editor
        self.kept_hidden = set()   # the annotation types hidden one by one before Annotations went off, given back when on
        ed = editor
        # (icon, tooltip, shown?, show/hide)
        self.categories = [
            ("annotations", "Annotations (features)", self._annotations_on, self._set_annotations),
            ("primers", "Primers", lambda: ed.show_primers, lambda on: ed.set_display(primers=on)),
            ("enzymes", "Restriction sites", lambda: ed.show_enzymes, lambda on: ed.set_display(enzymes=on)),
            ("cds", "CDS translations (amino acids)", lambda: ed.show_cds_translation,
             lambda on: ed.set_display(cds_translation=on)),
            ("units", "Transcription units", lambda: ed.show_units, lambda on: ed.set_display(units=on)),
        ]
        row = QHBoxLayout(self)
        row.setContentsMargins(0, 0, 4, 0)
        row.setSpacing(1)
        self.buttons = {}
        for name, tip, _, slot in self.categories:
            b = QToolButton(checkable=True, autoRaise=True, focusPolicy=Qt.NoFocus)
            b.setIconSize(QSize(18, 18))
            b.toggled.connect(slot)
            row.addWidget(b)
            self.buttons[name] = b
        self._paint()
        ed.changed.connect(self.refresh)            # annotation types may appear or vanish with edits
        ed.display_changed.connect(self.refresh)    # the side panel and the menus change the same switches
        self.refresh()

    def _annotations_on(self):
        kinds = set(self.editor.kind_counts)
        return bool(kinds - self.editor.hidden_kinds) if kinds else not self.editor.hidden_kinds

    def _set_annotations(self, on):
        ed = self.editor
        kinds = set(ed.kind_counts)
        if on:
            ed.set_display(hidden_kinds=self.kept_hidden if kinds - self.kept_hidden else ())
        else:
            self.kept_hidden = set(ed.hidden_kinds)
            ed.set_display(hidden_kinds=kinds | ed.hidden_kinds)

    def _paint(self):
        """Icons and backgrounds for the current theme."""
        self.setStyleSheet(f"QToolButton {{ border: none; border-radius: 3px; padding: 2px; }}"
                           f"QToolButton:hover {{ background: {theme.hex('toggle_hover_bg')}; }}"
                           f"QToolButton:checked {{ background: {theme.hex('highlight')}; }}")
        for name, *_ in self.categories:
            self.buttons[name].setIcon(icons.category_icon(name, self.palette(), theme.color("highlight_text")))   # white on the blue

    def refresh(self):
        for name, tip, shown, _ in self.categories:
            b, on = self.buttons[name], shown()
            b.blockSignals(True)   # mirror the state without sending it back
            b.setChecked(on)
            b.blockSignals(False)
            b.setToolTip(f"{tip}: {'shown' if on else 'hidden'} (click to {'hide' if on else 'show'}). "
                         "Right-click: the full display filter.")

    def contextMenuEvent(self, event):
        self.window.display_dock.show()

    def changeEvent(self, event):
        if event.type() == QEvent.PaletteChange:   # a theme switch: the icons follow (black on light, white on dark)
            self._paint()
        super().changeEvent(event)


_slim = []


def _slim_arrows():
    """The application's style with the up / down buttons of a spin box half as wide (one shared instance: a widget does not own
    its style)."""
    if not _slim:
        class Slim(QProxyStyle):
            def subControlRect(self, cc, opt, sc, widget=None):
                r = super().subControlRect(cc, opt, sc, widget)
                if cc != QStyle.CC_SpinBox or sc not in (QStyle.SC_SpinBoxUp, QStyle.SC_SpinBoxDown, QStyle.SC_SpinBoxEditField):
                    return r
                up = super().subControlRect(cc, opt, QStyle.SC_SpinBoxUp, widget)
                down = super().subControlRect(cc, opt, QStyle.SC_SpinBoxDown, widget)
                half = max(8, up.width() // 2)   # half the usual width: slim, still wide enough to see and hit
                if up.x() == down.x():   # stacked (Fusion): both at the right edge, half as wide
                    left = up.right() + 1 - half
                    lefts = {QStyle.SC_SpinBoxUp: left, QStyle.SC_SpinBoxDown: left}
                else:                    # side by side (Windows 11): the right-hand one at the edge, the other next to it
                    edge = max(up.right(), down.right()) + 1 - half
                    other = edge - half
                    lefts = ({QStyle.SC_SpinBoxUp: other, QStyle.SC_SpinBoxDown: edge} if up.x() < down.x()
                             else {QStyle.SC_SpinBoxUp: edge, QStyle.SC_SpinBoxDown: other})
                if sc == QStyle.SC_SpinBoxEditField:
                    r.setRight(min(lefts.values()) - 1)   # the text gets the room the buttons give up
                else:
                    r.setLeft(lefts[sc])
                    r.setWidth(half)
                return r
        _slim.append(Slim(QApplication.style().name()))   # the base style's name (the application's style is a PressedStyle proxy)
    return _slim[0]


FIELD = 64   # every field of the enzyme criteria (numbers and choices) is at least this wide: tidy, and it fits the panel


def field_width(widget):
    """The width of every criteria field: FIELD, or more when the font is larger (125-200 % text, a bigger system font), so that
    "999" or "any" and the slim arrows both fit (a fixed 64 px crushed the arrows to a sliver)."""
    return max(FIELD, widget.fontMetrics().horizontalAdvance("999") + 40)
ICON, GAP = 20, 6   # the icon of a category (at the right end of its first row, by the bar above it) and the gap before it


class EnzymeCriteria(QWidget):
    """The criteria of an enzyme filter (enzymes.EnzymeFilter), all of which must hold: how many sites, how long the recognition
    site, palindromic, Type IIS, the ends left, in the inventory; and the common choices in one click (PRESETS). Used by the
    Display panel (which sites are drawn) and the Digest (which enzymes are listed). `changed` sends the new filter."""
    changed = Signal(object)

    def __init__(self, presets=PRESETS):
        super().__init__()
        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        form = QFormLayout()
        form.setContentsMargins(0, 0, 0, 0)
        form.setHorizontalSpacing(6)
        form.setVerticalSpacing(2)

        def spin(low, tip, none_text=None):
            """A number from `low`; with `none_text`, one less than `low` stands for "no limit" and shows as that text."""
            box = QSpinBox(minimum=low - 1 if none_text else low, maximum=999, toolTip=tip)
            box.setFixedWidth(field_width(box))
            box.setStyle(_slim_arrows())
            if none_text:
                box.setSpecialValueText(none_text)   # 0 shows as "any"
            box.valueChanged.connect(self._changed)
            return box

        def pair(low, high, label):
            row = QHBoxLayout()
            row.setContentsMargins(0, 0, 0, 0)
            row.setSpacing(5)   # "to" with room on both sides, not touching the boxes
            for w in (low, QLabel("to"), high):
                row.addWidget(w)
            row.addStretch()
            form.addRow(label, row)

        # Sites can be 0: the enzymes that do not cut (0 to 0), e.g. to clone into a site the insert lacks.
        self.min_sites, self.max_sites = (spin(0, "At least this many sites in this sequence (0: also the enzymes that do not cut)"),
                                          spin(0, "At most this many (no limit: any; 0 to 0: the enzymes that do not cut)", "any"))
        pair(self.min_sites, self.max_sites, "Sites")
        self.min_size, self.max_size = spin(1, "Recognition site at least this many bases"), spin(1, "At most this many (no limit: any)", "any")
        pair(self.min_size, self.max_size, "Size")

        def choice(label, options):
            box = QComboBox()
            for text, value in options:
                box.addItem(text, value)
            box.setFixedWidth(field_width(box))
            box.currentIndexChanged.connect(self._changed)
            form.addRow(label, box)
            return box

        either = [("Any", None), ("Yes", True), ("No", False)]
        self.palindromic = choice("Palindromic", either)
        self.type_iis = choice("Type IIS", either)
        self.ends = choice("Ends", [("Any", None), ("5′", "5"), ("3′", "3"), ("Blunt", "blunt")])
        self.ends.setToolTip("The ends the enzyme leaves: a 5′ overhang, a 3′ overhang or blunt")
        self.inventory = choice("Owned", either)
        self.inventory.setToolTip("Only the enzymes tagged as owned in the lab (the Owned box of the Digest's enzyme list), or not")
        outer.addLayout(form)
        buttons = QHBoxLayout()   # the common choices in one click, side by side on one row
        buttons.setContentsMargins(0, 4, 0, 0)
        buttons.setSpacing(3)
        for i, (label, f, tip) in enumerate(presets):
            button = QPushButton(label, toolTip=tip, clicked=lambda _=False, f=f: self._preset(f))
            button.setStyleSheet("padding: 2px 4px;")   # slim buttons in the normal font (a 9 px font was hard to read)
            button.setFixedWidth(button.fontMetrics().horizontalAdvance(label) + 16)
            buttons.addWidget(button)
        buttons.addStretch()
        outer.addLayout(buttons)
        self.fields = [self.min_sites, self.max_sites, self.min_size, self.max_size, self.palindromic, self.type_iis, self.ends,
                       self.inventory]

    def filter(self):
        data = lambda c: c.currentData()
        limit = lambda box: None if box.value() == box.minimum() else box.value()   # the lowest value is "any"
        return EnzymeFilter(self.min_sites.value(), limit(self.max_sites), self.min_size.value(), limit(self.max_size),
                            data(self.palindromic), data(self.type_iis), data(self.ends), data(self.inventory))

    def set_filter(self, f):
        """Show `f` without sending it back."""
        for box in self.fields:
            box.blockSignals(True)
        for spin, value in ((self.min_sites, f.min_sites), (self.max_sites, self.max_sites.minimum() if f.max_sites is None else f.max_sites),
                            (self.min_size, f.min_size), (self.max_size, self.max_size.minimum() if f.max_size is None else f.max_size)):
            spin.setValue(value)
        for combo, value in ((self.palindromic, f.palindromic), (self.type_iis, f.type_iis), (self.ends, f.ends),
                             (self.inventory, f.inventory)):
            combo.setCurrentIndex(combo.findData(value))
        for box in self.fields:
            box.blockSignals(False)

    def _preset(self, f):
        self.set_filter(f)
        self.changed.emit(f)

    def _changed(self, *_):
        self.changed.emit(self.filter())


class DisplayPanel(QScrollArea):
    def __init__(self, editor):
        super().__init__()
        self.editor = editor
        self.setWidgetResizable(True)
        self.kinds = {}  # feature type -> its check box
        self.icons = {}  # category -> the label holding its icon (the same pictures as the quick bar's toggles)
        body = QWidget()
        self.box = QVBoxLayout(body)
        self.box.setAlignment(Qt.AlignTop)
        self.box.addWidget(QLabel("<b>Show</b>"))
        self.all_primers = self._switch("Primers", editor.show_primers, lambda on: editor.set_display(primers=on), icon="primers",
                                        tip="All primers at once; the two sources below can be switched one by one")
        self.file_primers, self.file_swatch = self._switch(   # two sources, stacked (not side by side) so the panel stays narrow
            "From this sequence", editor.show_file_primers, lambda on: editor.set_display(file_primers=on),
            swatch=PrimerSwatch(False),
            tip="Primers that are part of this file (stored in it or imported to it), drawn in the regular primer colour.")
        self.db_primers, self.db_swatch = self._switch(
            "From the database", editor.show_db_primers, lambda on: editor.set_display(db_primers=on),
            swatch=PrimerSwatch(True),
            tip="Primers that are only in the primer database (not in this file), drawn like temporary primers. "
                "Right-click one > Import to This Sequence to make it part of the file.")
        self._rule()
        self.expander = QToolButton(autoRaise=True, checkable=True, arrowType=Qt.RightArrow, focusPolicy=Qt.NoFocus,
                                    toolTip="Show / hide the criteria of the sites drawn (number of sites, size, ends, inventory...)")
        self.expander.setFixedWidth(max(8, self.expander.sizeHint().width() // 2))   # half the usual width: just an arrow
        self.enzymes = self._switch("Restriction sites", editor.show_enzymes, lambda on: editor.set_display(enzymes=on),
                                    icon="enzymes", trailing=self.expander)
        self._enzyme_criteria()   # folded away until the arrow opens it: the annotation types matter more
        self._rule()
        self.cds = self._switch("CDS translations", editor.show_cds_translation,
                                lambda on: editor.set_display(cds_translation=on), icon="cds")
        self.frames = self._switch("3-frame translation", editor.show_translation,
                                   lambda on: editor.set_display(translation=on))
        self.strands = self._switch("Complementary strand", editor.show_complement,
                                    lambda on: editor.set_display(complement=on),
                                    tip="The bottom strand under the top one. Off, the sequence view shows the top strand alone.")
        self._rule()
        self.units = self._switch("Transcription units", editor.show_units, lambda on: editor.set_display(units=on), icon="units")
        self._rule()
        self.hidden_row = QWidget()   # features hidden one by one (right-click > Hide): how many, and a way back
        hidden_layout = QHBoxLayout(self.hidden_row)
        hidden_layout.setContentsMargins(0, 0, 0, 0)
        self.hidden_label = QLabel()
        hidden_layout.addWidget(self.hidden_label)
        self.show_hidden = QPushButton("Show all", clicked=lambda: editor.set_features_hidden(editor.hidden_features, False))
        hidden_layout.addWidget(self.show_hidden)
        self.box.addWidget(self.hidden_row)
        self._heading("annotations", "<b>Annotation types</b>")

        buttons = QHBoxLayout()
        buttons.addWidget(QPushButton("All", clicked=lambda: editor.set_display(hidden_kinds=())))
        buttons.addWidget(QPushButton("None", clicked=lambda: editor.set_display(
            hidden_kinds=set(editor.kind_counts))))
        self.box.addLayout(buttons)
        self.kinds_box = QVBoxLayout()
        self.kinds_box.setContentsMargins(0, 0, 0, 0)
        self.box.addLayout(self.kinds_box)
        self.setWidget(body)
        editor.changed.connect(self.refresh)          # types may appear or vanish with edits
        editor.display_changed.connect(self.refresh)  # menus change the same switches
        tool_settings.notifier.changed.connect(self._code_changed)  # the genetic code is named on the translation switch
        self.refresh()

    def _enzyme_criteria(self):
        """Which restriction sites are drawn: independent criteria that all have to hold (enzymes.EnzymeFilter), the same form
        as the Digest's filter bar (EnzymeCriteria)."""
        self.criteria_box = QWidget()   # what the arrow next to the switch folds
        inner = QVBoxLayout(self.criteria_box)
        inner.setContentsMargins(GAP * 2, 0, 0, 0)
        self.criteria_form = EnzymeCriteria()
        self.criteria_form.changed.connect(lambda f: self.editor.set_display(enzyme_filter=f))
        inner.addWidget(self.criteria_form)
        for name in ("min_sites", "max_sites", "min_size", "max_size", "palindromic", "type_iis", "ends", "inventory"):
            setattr(self, name, getattr(self.criteria_form, name))
        self.box.addWidget(self.criteria_box)
        self.criteria_box.setVisible(False)
        self.expander.toggled.connect(self._fold_criteria)
        self.criteria = self.criteria_form.fields

    def _fold_criteria(self, open_):
        self.expander.setArrowType(Qt.DownArrow if open_ else Qt.RightArrow)
        self.criteria_box.setVisible(open_)

    def _rule(self):
        """A thin bar between two categories."""
        line = QFrame(frameShape=QFrame.HLine, frameShadow=QFrame.Plain)
        line.setStyleSheet(f"color: {theme.hex('pill_edge')};")
        line.theme_changed = lambda line=line: line.setStyleSheet(f"color: {theme.hex('pill_edge')};")   # theme.refresh_views
        line.setFixedHeight(1)
        self.box.addSpacing(3)
        self.box.addWidget(line)
        self.box.addSpacing(3)

    def _icon_label(self, category):
        label = QLabel()
        label.setFixedSize(ICON, ICON)
        self.icons[category] = label
        return label

    def _heading(self, category, text):
        """The title of a category: its icon (as on the quick bar) and `text`."""
        row = QHBoxLayout()
        row.setContentsMargins(0, 0, 0, 0)
        row.setSpacing(GAP)
        row.addWidget(QLabel(text), 1)
        row.addWidget(self._icon_label(category))
        self.box.addLayout(row)

    def _paint_icons(self):
        """Each category's icon, clear while it is shown and faded when not (as the quick bar's toggles)."""
        ed = self.editor
        kinds = set(ed.kind_counts)
        shown = {"primers": ed.show_primers, "enzymes": ed.show_enzymes, "cds": ed.show_cds_translation, "units": ed.show_units,
                 "annotations": bool(kinds - ed.hidden_kinds) if kinds else not ed.hidden_kinds}
        for name, label in self.icons.items():
            state = QIcon.On if shown[name] else QIcon.Off
            label.setPixmap(icons.category_icon(name, self.palette()).pixmap(18, 18, QIcon.Normal, state))

    def changeEvent(self, event):
        if event.type() == QEvent.PaletteChange and self.icons:   # a theme switch: the icons follow
            self._paint_icons()
        super().changeEvent(event)

    def _switch(self, text, on, slot, swatch=None, tip="", icon=None, trailing=None):
        """A check box in a row of the panel, always at the left edge. `icon`: the category's icon, at the right end of the row.
        `swatch` (a PrimerSwatch) sits right after the box, `trailing` (a widget) next to it. Returns the box, or (box, swatch)."""
        box = QCheckBox(text, checked=on)
        box.toggled.connect(slot)
        row = QHBoxLayout()
        row.setContentsMargins(0, 0, 0, 0)
        row.setSpacing(GAP)
        row.addWidget(box)
        if swatch is not None:
            row.addWidget(swatch)
        if trailing is not None:
            row.addWidget(trailing)
        row.addStretch()
        if icon is not None:
            row.addWidget(self._icon_label(icon))
        self.box.addLayout(row)
        for w in (box, swatch):
            if w is not None and tip:
                w.setToolTip(tip)
        return box if swatch is None else (box, swatch)

    def _code_changed(self, tool):
        if tool == genetic_code.TOOL:
            self.refresh()

    def refresh(self):
        ed = self.editor
        self._paint_icons()
        self.frames.setText("3-frame translation" + genetic_code.table_note())   # "· Table 11" when not the standard code
        for box, on in ((self.all_primers, ed.show_primers), (self.file_primers, ed.show_file_primers), (self.db_primers, ed.show_db_primers), (self.enzymes, ed.show_enzymes),
                        (self.cds, ed.show_cds_translation), (self.frames, ed.show_translation), (self.strands, ed.show_complement),
                        (self.units, ed.show_units)):
            box.blockSignals(True)  # mirror the state without sending it back
            box.setChecked(on)
            box.blockSignals(False)
        self.criteria_form.set_filter(ed.enzyme_filter)   # mirrored, not sent back
        self.criteria_form.setEnabled(ed.show_enzymes)
        # always there (appearing and vanishing, it pushed the list of types below up and down): greyed with nothing hidden
        self.show_hidden.setEnabled(bool(ed.hidden_features))
        self.hidden_label.setText(f"{plural(len(ed.hidden_features), 'feature')} hidden by hand" if ed.hidden_features
                                  else "No feature hidden by hand")
        counts = ed.kind_counts   # (units have their own switch)
        if sorted(counts) != sorted(self.kinds):  # rebuild the list when the types change
            while self.kinds_box.count():
                self.kinds_box.takeAt(0).widget().deleteLater()
            self.kinds = {}
            for kind in sorted(counts, key=str.lower):
                box = QCheckBox()
                box.toggled.connect(lambda on, k=kind: self._kind(k, on))
                self.kinds_box.addWidget(box)
                self.kinds[kind] = box
        for kind, box in self.kinds.items():
            box.blockSignals(True)
            box.setText(f"{kind}  ({counts[kind]})")
            box.setChecked(kind not in ed.hidden_kinds)
            box.blockSignals(False)

    def _kind(self, kind, on):
        hidden = set(self.editor.hidden_kinds)
        (hidden.discard if on else hidden.add)(kind)
        self.editor.set_display(hidden_kinds=hidden)

