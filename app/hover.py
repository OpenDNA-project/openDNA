"""Hovering a feature for a moment (Settings → Display, 0.3 s by default) shows its name,
type, location, length and description in a tooltip. Used by the sequence view, map and
genome view: each gives a function point -> feature index (or None).
"""
from html import escape

from PySide6.QtCore import QEvent, QObject, QPoint, QSize, Qt, QTimer
from PySide6.QtWidgets import QApplication, QLabel, QToolTip

import genetic_code
import tool_settings
from tool_settings import Param

TOOL = "Display"
tool_settings.register(TOOL, [
    Param("hover_delay", "Show a feature's name after hovering it for", 0.3, 0.1, 2.0, step=0.1, suffix=" s",
          tip="In the sequence view, map and genome view."),
])
DESCRIBE = ("gene", "product", "locus_tag", "note", "function")  # qualifiers worth showing


class Text(str):
    """A tooltip given as its own text (rich), for what is neither a feature, a primer nor a cut: a track's name in the
    alignment view. Two of them with the same text are the same tip (see HoverTips._show)."""


def describe_primer(editor, key):
    """Tooltip text (rich) for the primer site `key`: its bases (the tail in lowercase orange), where it binds, Tm and the checks."""
    from primer_panel import melting_temp, structure   # here: primer_panel imports the views
    site = editor.site_by_key.get(key)
    if site is None:
        return ""
    tail, core = site.sequence[:site.tail], site.sequence[site.tail:]
    where = f"{site.start + 1:,}–{(site.end - 1) % editor.length + 1:,} {'→' if site.strand > 0 else '←'}"
    t = melting_temp(core)
    s = structure(site.sequence) if len(site.sequence) >= 8 else None
    if site.temporary:
        kind = "temporary primer (not in the primer database yet)"
    elif site.provisional:
        kind = "from the primer database (not in this file)"
    else:
        kind = "primer in this file"
    lines = [f"<b>{escape(site.label)}</b> &nbsp;<i>{kind}</i>",
             f"5′-<span style='color:#bf360c'>{escape(tail.lower())}</span>{escape(core.upper())}-3′ &nbsp;·&nbsp; {len(site.sequence)} nt"
             + (f" ({site.tail} nt tail)" if site.tail else ""),
             f"binds {where}" + (f" with {site.mismatches} mismatch(es)" if site.mismatches else "") + (f" &nbsp;·&nbsp; Tm {t:.1f} °C" if t else "")]
    if s:
        lines.append(f"hairpin ΔG {min(s['hairpin'], 0):.1f}, self-dimer ΔG {s['self']:.1f} kcal/mol" + (f" &nbsp;⚠ {escape('; '.join(s['warnings']))}" if s["warnings"] else ""))
    if site.temporary:
        lines.append("<i>Save primers to the database to keep it; right-click: Import to This Sequence.</i>")
    elif site.provisional:
        lines.append("<i>Right-click: Import to This Sequence, to keep it in this file.</i>")
    return "<br>".join(lines)


def describe_cut(editor, key):
    """Tooltip text (rich) for a cut drawn on the sequence: key ("cut", owner, i) is a mark of a tool (a digest, the Spacer
    Finder: its `tip` is the text or a function that makes it), key ("site", enzyme name, cut position) an enzyme name drawn by
    View > Restriction Sites."""
    from cut_sites import unique_site_tip
    if key[0] == "site":
        return unique_site_tip(editor, key[1], key[2])
    marks = editor.cut_marks.get(key[1], [])
    tip = marks[key[2]].tip if key[2] < len(marks) else None
    return tip() if callable(tip) else tip or ""


def describe(editor, index):
    """Tooltip text (rich) for feature `index` (an int), primer site (its key, a str) or cut (a tuple, see describe_cut), or
    a Text (its own text)."""
    if isinstance(index, Text):
        return str(index)
    if isinstance(index, str):
        return describe_primer(editor, index)
    if isinstance(index, tuple):
        return describe_cut(editor, index)
    f = editor.feature_by_index(index)
    if f is None:
        return ""
    from cut_sites import cas9_tip, feature_layout
    layout = feature_layout(f, editor.length)   # a saved Cas9 / Cas12a guide site: say where it cuts
    if layout is not None:
        return cas9_tip(layout, editor.doc.sequence, editor.length, escape(f.label or f.kind))
    size = sum(e - s for s, e in f.spans)
    lines = [f"<b>{escape(f.label or f.kind)}</b> &nbsp;<i>{escape(f.kind)}</i>",
             f"{escape(f.location)} &nbsp;·&nbsp; {size:,} bp" + (f" ({size // 3:,} aa)" if f.kind == "CDS" else "")]
    for key, value in f.qualifiers:
        if key in DESCRIBE and value and value != f.label:
            lines.append(f"{key}: {escape(value[:200])}")
    if f.kind == "CDS":   # which genetic code it is read with, when not the standard one; and a translation we cannot follow
        note = genetic_code.cds_note(editor.doc, f)
        if note:
            lines.append(f"<i>{escape(note)}</i>")
    return "<br>".join(lines)


class HoverTips(QObject):
    """Watches mouse moves over `widget`; when the pointer rests on the same feature or primer for
    the set delay, shows its tooltip (hidden again when the pointer leaves it)."""

    def __init__(self, widget, editor, feature_at, on_change=None, anchor=None):
        """`on_change(key or None)`: told what is hovered once its tooltip is up, and None when it goes (the sequence view
        draws a restriction site's recognition site then). `anchor(pt, key)`: the points (viewport) above and below which a cut's tooltip is put."""
        super().__init__(widget)
        self.widget, self.editor, self.feature_at = widget, editor, feature_at
        self.on_change, self.anchor = on_change, anchor
        self.pos = QPoint()
        self.shown = None  # feature whose tooltip is up
        self.timer = QTimer(self, singleShot=True, timeout=self._show)
        widget.setMouseTracking(True)
        widget.installEventFilter(self)

    def eventFilter(self, obj, event):
        if event.type() == QEvent.MouseMove:
            self.pos = event.position().toPoint()
            if self.shown is not None and self.feature_at(self.pos) != self.shown:
                self._hide()
            self.timer.start(int(tool_settings.value(TOOL, "hover_delay") * 1000))
        elif event.type() in (QEvent.Leave, QEvent.MouseButtonPress, QEvent.Wheel):
            self.timer.stop()
            self._hide()
        elif event.type() == QEvent.ToolTip:
            return True  # our own tooltips only
        return False

    def _hide(self):
        QToolTip.hideText()
        if self.shown is not None and self.on_change:
            self.on_change(None)
        self.shown = None

    def _show(self):
        index = self.feature_at(self.pos)
        if index is not None and index == self.shown:
            return   # its tip is up already: it does not jump to the pointer, nor is its text (a ΔG, a pill) worked out again
        if index is not None and self.widget.isVisible():
            self.shown = index
            html = describe(self.editor, index)
            at, cut = self.pos, isinstance(index, tuple) and self.anchor is not None
            if cut:   # the tip goes above the enzyme name, off the sequence: worked out before it is shown, so it never jumps
                above, below = (self.widget.mapToGlobal(q) for q in self.anchor(self.pos, index))
                size = self._size_of(html)
                screen = self.widget.screen().availableGeometry()
                y = above.y() - size.height() - 6
                if y < screen.top():
                    y = below.y() + 6
                x = max(screen.left(), min(above.x(), screen.right() - size.width()))
                QToolTip.showText(QPoint(x - 2, y - 16), html, self.widget)   # Qt puts the tip at (+2, +16) from this point
                self._lift(above, below)   # settles the last pixels if the estimate of the size was off
            else:
                QToolTip.showText(self.widget.mapToGlobal(at), html, self.widget)
            if self.on_change:
                self.on_change(index)

    @staticmethod
    def _size_of(html):
        """About the size the tooltip of this rich text will have (a label in the tooltip font, plus its frame and margin)."""
        label = QLabel(html)
        label.setFont(QToolTip.font())
        label.setTextFormat(Qt.RichText)
        label.setWordWrap(True)   # as the tooltip does for rich text: it wraps a long line
        return label.sizeHint() + QSize(6, 6)

    @staticmethod
    def _lift(point, under):
        """Move the tooltip that is up so its bottom edge is just above `point` (Qt puts it below and right of the pointer);
        with no room above on the screen, just under `under` instead. Either way it stays off the sequence."""
        for w in QApplication.topLevelWidgets():
            if w.isVisible() and w.metaObject().className() == "QTipLabel":
                screen = w.screen().availableGeometry()
                y = point.y() - w.height() - 6
                if y < screen.top():
                    y = under.y() + 6
                x = min(point.x(), screen.right() - w.width())
                w.move(max(screen.left(), x), min(y, screen.bottom() - w.height()))
