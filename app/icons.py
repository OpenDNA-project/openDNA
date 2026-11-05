"""Small flat icons drawn here (no image files): one colour, the text colour of the palette, so they read as black on the
light theme and white on the dark one. A funnel for the display filter, an hourglass for the history, and one per category
of the quick filters in the tab row."""
from PySide6.QtCore import QPointF, QRectF, Qt
from PySide6.QtGui import QFont, QIcon, QPainter, QPalette, QPen, QPixmap, QPolygonF
from PySide6.QtWidgets import QApplication

import theme

SIZE = 64   # drawn big; Qt scales it down to the button


def _icon(draw, palette=None):
    """The icon in the text colour; pressed (QIcon.On, on the dark blue of a pressed button), in white."""
    icon = QIcon()
    for state, color in ((QIcon.Off, (palette or QApplication.palette()).color(QPalette.WindowText)),
                         (QIcon.On, theme.color("highlight_text"))):
        pix = QPixmap(SIZE, SIZE)
        pix.fill(Qt.transparent)
        p = QPainter(pix)
        p.setRenderHint(QPainter.Antialiasing)
        p.setPen(Qt.NoPen)
        p.setBrush(color)
        draw(p, SIZE)
        p.end()
        icon.addPixmap(pix, QIcon.Normal, state)
    return icon


def _funnel(p, s):
    p.drawPolygon(QPolygonF([QPointF(s * 0.08, s * 0.14), QPointF(s * 0.92, s * 0.14), QPointF(s * 0.58, s * 0.52),
                             QPointF(s * 0.58, s * 0.88), QPointF(s * 0.42, s * 0.78), QPointF(s * 0.42, s * 0.52)]))


def _hourglass(p, s):
    p.drawRect(QRectF(s * 0.18, s * 0.06, s * 0.64, s * 0.08))   # the two plates
    p.drawRect(QRectF(s * 0.18, s * 0.86, s * 0.64, s * 0.08))
    p.drawPolygon(QPolygonF([QPointF(s * 0.24, s * 0.16), QPointF(s * 0.76, s * 0.16), QPointF(s * 0.54, s * 0.5),
                             QPointF(s * 0.76, s * 0.84), QPointF(s * 0.24, s * 0.84), QPointF(s * 0.46, s * 0.5)]))


# ---- the quick filters of the tab row (display_panel.QuickFilters): one icon per broad category of the display filter ----
def _line(p, s, width=0.07):
    """Draw with lines from here on (same colour as the fills)."""
    p.setPen(QPen(p.brush().color(), s * width, Qt.SolidLine, Qt.RoundCap, Qt.RoundJoin))


def _feature(p, s):
    """An annotation: a feature arrow, as the map and sequence draw them."""
    p.drawPolygon(QPolygonF([QPointF(s * 0.06, s * 0.32), QPointF(s * 0.66, s * 0.32), QPointF(s * 0.94, s * 0.5),
                             QPointF(s * 0.66, s * 0.68), QPointF(s * 0.06, s * 0.68)]))


def _primer(p, s):
    """A primer: a line with a half arrow at its 3' end, over a strand (thinner)."""
    _line(p, s, 0.09)
    p.drawLine(QPointF(s * 0.08, s * 0.42), QPointF(s * 0.9, s * 0.42))
    p.drawLine(QPointF(s * 0.9, s * 0.42), QPointF(s * 0.7, s * 0.24))
    _line(p, s, 0.04)
    p.drawLine(QPointF(s * 0.04, s * 0.72), QPointF(s * 0.96, s * 0.72))


def _scissors(p, s):
    """A restriction enzyme: scissors."""
    _line(p, s, 0.08)
    p.setBrush(Qt.NoBrush)
    p.drawEllipse(QPointF(s * 0.24, s * 0.74), s * 0.13, s * 0.13)   # the two rings
    p.drawEllipse(QPointF(s * 0.76, s * 0.74), s * 0.13, s * 0.13)
    p.drawLine(QPointF(s * 0.33, s * 0.64), QPointF(s * 0.74, s * 0.08))   # the blades, crossing
    p.drawLine(QPointF(s * 0.67, s * 0.64), QPointF(s * 0.26, s * 0.08))


def _start_codon(p, s):
    """CDS translations: the start codon, ATG, over the bracket that marks a codon (big: it is read at 18 px)."""
    font = QFont()
    font.setBold(True)
    font.setPixelSize(round(s * 0.5))
    font.setStretch(QFont.Condensed)
    p.setFont(font)
    p.setPen(p.brush().color())
    p.drawText(QRectF(-s * 0.1, s * 0.08, s * 1.2, s * 0.6), Qt.AlignCenter, "ATG")
    _line(p, s, 0.08)
    p.drawPolyline(QPolygonF([QPointF(s * 0.06, s * 0.7), QPointF(s * 0.06, s * 0.84), QPointF(s * 0.94, s * 0.84),
                              QPointF(s * 0.94, s * 0.7)]))


def _unit(p, s):
    """A transcription unit: on a line, the bent arrow of a promoter and a hairpin terminator."""
    _line(p, s, 0.07)
    p.drawLine(QPointF(s * 0.02, s * 0.88), QPointF(s * 0.98, s * 0.88))   # the DNA
    p.drawLine(QPointF(s * 0.1, s * 0.88), QPointF(s * 0.1, s * 0.46))     # the promoter, bent to the right
    p.drawLine(QPointF(s * 0.1, s * 0.46), QPointF(s * 0.3, s * 0.46))
    p.setPen(Qt.NoPen)
    p.drawPolygon(QPolygonF([QPointF(s * 0.44, s * 0.46), QPointF(s * 0.28, s * 0.32), QPointF(s * 0.28, s * 0.6)]))
    _line(p, s, 0.07)
    p.setBrush(Qt.NoBrush)
    p.drawLine(QPointF(s * 0.62, s * 0.88), QPointF(s * 0.62, s * 0.5))   # the hairpin: a stem of two strands, apart...
    p.drawLine(QPointF(s * 0.86, s * 0.88), QPointF(s * 0.86, s * 0.5))
    p.drawEllipse(QRectF(s * 0.56, s * 0.06, s * 0.36, s * 0.5))           # ... closed by a loop


CATEGORY_DRAW = {"annotations": _feature, "primers": _primer, "enzymes": _scissors, "cds": _start_codon, "units": _unit}
OFF_OPACITY = 0.3   # a category switched off: the same icon, faded


def category_icon(name, palette=None, on_color=None):
    """The icon of a display filter category (CATEGORY_DRAW), clear when it is shown (QIcon.On) and faded when not (Off). `on_color`:
    the clear one in this colour (white, on the highlight of a pressed quick filter) instead of the text colour."""
    color = (palette or QApplication.palette()).color(QPalette.WindowText)
    icon = QIcon()
    for state, opacity in ((QIcon.On, 1.0), (QIcon.Off, OFF_OPACITY)):
        pix = QPixmap(SIZE, SIZE)
        pix.fill(Qt.transparent)
        p = QPainter(pix)
        p.setRenderHint(QPainter.Antialiasing)
        p.setOpacity(opacity)
        p.setPen(Qt.NoPen)
        p.setBrush(on_color if on_color is not None and state == QIcon.On else color)
        CATEGORY_DRAW[name](p, SIZE)
        p.end()
        icon.addPixmap(pix, QIcon.Normal, state)
    return icon


def filter_icon(palette=None):
    return _icon(_funnel, palette)


def history_icon(palette=None):
    return _icon(_hourglass, palette)


# ---- the theme button of the toolbar: a sun (Light), a dark sun (Dark), a sun with an A (System: follows the computer) ----
def _rays(p, s, inner, outer):
    import math
    _line(p, s, 0.075)
    for i in range(8):
        a = i * math.pi / 4
        p.drawLine(QPointF(s / 2 + inner * s * math.cos(a), s / 2 + inner * s * math.sin(a)),
                   QPointF(s / 2 + outer * s * math.cos(a), s / 2 + outer * s * math.sin(a)))


def _sun(p, s):
    """Light: an open sun, a ring and its rays."""
    _rays(p, s, 0.33, 0.47)
    p.setBrush(Qt.NoBrush)
    p.drawEllipse(QPointF(s / 2, s / 2), s * 0.2, s * 0.2)


def _dark_sun(p, s):
    """Dark: the same sun, filled in."""
    color = p.pen().color() if p.pen().style() != Qt.NoPen else p.brush().color()
    _rays(p, s, 0.33, 0.47)
    p.setBrush(color)
    p.drawEllipse(QPointF(s / 2, s / 2), s * 0.22, s * 0.22)


def _auto_sun(p, s):
    """System: the open sun, with an A inside (automatic)."""
    color = p.brush().color()
    _rays(p, s, 0.36, 0.49)
    p.setBrush(Qt.NoBrush)
    p.drawEllipse(QPointF(s / 2, s / 2), s * 0.26, s * 0.26)
    font = QFont()
    font.setBold(True)
    font.setPixelSize(round(s * 0.34))
    p.setFont(font)
    p.setPen(color)
    p.drawText(QRectF(0, 0, s, s), Qt.AlignCenter, "A")


THEME_DRAW = {"light": _sun, "dark": _dark_sun, "system": _auto_sun}


def theme_icon(name, palette=None):
    """The icon of a theme of theme.THEMES ("light", "dark", "system")."""
    return _icon(THEME_DRAW[name], palette)


# ---- "this sends data over the internet": a small globe (online.py puts it on such buttons) ----
def _globe(p, s):
    """A globe: a circle, its equator, two parallels and a meridian ellipse."""
    _line(p, s, 0.07)
    p.setBrush(Qt.NoBrush)
    c, r = s / 2, s * 0.42
    p.drawEllipse(QPointF(c, c), r, r)
    p.drawEllipse(QPointF(c, c), r * 0.42, r)                  # the meridians, seen from the side
    p.drawLine(QPointF(c - r, c), QPointF(c + r, c))           # the equator
    for dy in (-0.55, 0.55):                                   # two parallels
        half = r * (1 - dy * dy) ** 0.5
        p.drawLine(QPointF(c - half, c + dy * r), QPointF(c + half, c + dy * r))
    p.drawLine(QPointF(c, c - r), QPointF(c, c + r))


def globe_icon(palette=None):
    return _icon(_globe, palette)


def _genome(p, s):
    """A circular genome: a ring, thicker over a stretch (a gene on it)."""
    from PySide6.QtCore import QRectF
    _line(p, s, 0.08)
    p.setBrush(Qt.NoBrush)
    c, r = s / 2, s * 0.36
    p.drawEllipse(QPointF(c, c), r, r)
    _line(p, s, 0.2)
    p.drawArc(QRectF(c - r, c - r, 2 * r, 2 * r), 30 * 16, 75 * 16)


def genome_icon(palette=None):
    return _icon(_genome, palette)


# ---- the header of a pane (workspace.Pane): its options, split right / below, close. Drawn, so they look the same on every
# computer (the font symbols they were depended on the fonts installed: other sizes, other baselines, or empty boxes) ----
def _cog(p, s):
    import math
    from PySide6.QtGui import QPainterPath
    c, outer, inner, hole = s / 2, s * 0.44, s * 0.32, s * 0.13
    path = QPainterPath()
    teeth = 8
    for i in range(teeth * 4):
        a = (i / (teeth * 4)) * 2 * math.pi
        r = outer if i % 4 in (1, 2) else inner
        pt = QPointF(c + r * math.cos(a), c + r * math.sin(a))
        path.moveTo(pt) if i == 0 else path.lineTo(pt)
    path.closeSubpath()
    path.addEllipse(QPointF(c, c), hole, hole)   # (odd-even fill: the hole stays open)
    p.drawPath(path)


def _split(p, s, vertical):
    from PySide6.QtCore import QRectF
    _line(p, s, 0.08)
    p.setBrush(Qt.NoBrush)
    m = s * 0.14
    p.drawRoundedRect(QRectF(m, m, s - 2 * m, s - 2 * m), s * 0.06, s * 0.06)
    if vertical:
        p.drawLine(QPointF(s / 2, m), QPointF(s / 2, s - m))
    else:
        p.drawLine(QPointF(m, s / 2), QPointF(s - m, s / 2))


def _cross(p, s):
    _line(p, s, 0.09)
    m = s * 0.26
    p.drawLine(QPointF(m, m), QPointF(s - m, s - m))
    p.drawLine(QPointF(s - m, m), QPointF(m, s - m))


def pane_icon(name, palette=None):
    """'options' (a cog), 'right' / 'below' (split), 'close' (a cross)."""
    draw = {"options": _cog, "right": lambda p, s: _split(p, s, True), "below": lambda p, s: _split(p, s, False),
            "close": _cross}[name]
    return _icon(draw, palette)


# ---- the DNA lock of the toolbar: a padlock, shut (pressed: the bases cannot change) or open ----
def _padlock(p, s, shut):
    from PySide6.QtCore import QRectF
    body = QRectF(s * 0.2, s * 0.46, s * 0.6, s * 0.44)
    p.drawRoundedRect(body, s * 0.07, s * 0.07)
    _line(p, s, 0.09)
    p.setBrush(Qt.NoBrush)
    w = s * 0.36
    left = s * 0.32 if shut else s * 0.52   # open: the shackle swings to the right, off the body
    p.drawArc(QRectF(left, s * 0.12, w, w * 1.1), 0, 180 * 16)
    p.drawLine(QPointF(left + w, s * 0.12 + w * 0.55), QPointF(left + w, s * 0.46))
    if shut:
        p.drawLine(QPointF(left, s * 0.12 + w * 0.55), QPointF(left, s * 0.46))


def lock_icon(palette=None):
    """Off: an open padlock; On (pressed, white on the blue): shut."""
    icon = QIcon()
    for state, shut, color in ((QIcon.Off, False, (palette or QApplication.palette()).color(QPalette.WindowText)),
                               (QIcon.On, True, theme.color("highlight_text"))):
        pix = QPixmap(SIZE, SIZE)
        pix.fill(Qt.transparent)
        p = QPainter(pix)
        p.setRenderHint(QPainter.Antialiasing)
        p.setPen(Qt.NoPen)
        p.setBrush(color)
        _padlock(p, SIZE, shut)
        p.end()
        icon.addPixmap(pix, QIcon.Normal, state)
    return icon
