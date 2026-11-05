"""Draws the icon and the oDNA wordmark, and writes every file made from them.

Run:  python packaging/make_icon.py        (needs only PySide6, which the app already needs; no Pillow)

The outputs are committed, so a normal build does not run this. Re-run it only to change the picture:
    packaging/openDNA.ico          the exe icon, 7 sizes (16 to 256 px) on a light rounded tile
    app/logo/icon_<size>.png           the same sizes, for the window icon (main.py)
    app/logo/wordmark.png, ..._dark    the oDNA wordmark, for the welcome page (light and dark theme)
    docs/images/icon.svg, wordmark.svg, wordmark-dark.svg   the same as vector files, for the README

The mark: a plasmid drawn as an open green ring (the gap is the "open") with a dot in the middle, a promoter
(bent arrow, coral) where the ring starts and a terminator (T, ink) where it ends. Everything is built as filled
outlines (QPainterPath), text included, so the SVGs need no font and look the same everywhere.
The wordmark: the mark is the "o" of oDNA, followed by "DNA"; below them "open" (as wide as the ring) and
"DESIGN · NAVIGATION · ASSEMBLY" (as wide as "DNA"). The ring sits on the baseline of "DNA" and is a little
smaller than a real "o", so the promoter on top of it reaches about the top of the lowercase letters.
Alignment is done on the real outlines. The lettering is DM Sans (SIL Open Font License 1.1), loaded from
packaging/fonts/; it is needed only to run this script, the outputs carry the outlines.
Each icon size is drawn from scratch; below 32 px the promoter and terminator are left out (they turn to mush).
"""
import math
import struct
import sys
from pathlib import Path

from PySide6.QtCore import QBuffer, QByteArray, QIODevice, QPointF, QRectF, Qt
from PySide6.QtGui import (
    QColor, QFont, QFontDatabase, QGuiApplication, QImage, QPainter, QPainterPath, QPainterPathStroker, QTransform,
)

ROOT = Path(__file__).resolve().parent.parent
SIZES = (16, 24, 32, 48, 64, 128, 256)
ICO = ROOT / "packaging" / "openDNA.ico"
LOGO_DIR = ROOT / "app" / "logo"
DOCS = ROOT / "docs" / "images"
FONTS = ROOT / "packaging" / "fonts"
FONT_FAMILY = "DM Sans"
NAME = "oDNA"

# "muted" is the grey of the tagline under the name.
LIGHT = {"ring": "#1D9E75", "promoter": "#D85A30", "ink": "#2C2C2A", "muted": "#5F5E5A"}
DARK = {"ring": "#5DCAA5", "promoter": "#F0997B", "ink": "#F1EFE8", "muted": "#B4B2A9"}
TILE, TILE_EDGE = "#FBFAF7", "#D3D1C7"

# The mark, in a 100 x 100 box: ring centre (50, 50), radius R, stroke W. Angles in degrees, clockwise from 12 o'clock.
R, W = 33.0, 10.0
GAP = (304.0, 360.0)            # the ring runs from 0 to 304: the gap (the "Open") is top left
PROMOTER_AT, TERMINATOR_AT = 4.0, 300.0
LINE = 4.0                      # stroke of the promoter and terminator
REACH = 11.0                    # how far they stand out from the ring
CENTRE = QPointF(50, 50)        # the ring's centre


def pt(angle, radius):
    t = math.radians(angle - 90)
    return QPointF(50 + radius * math.cos(t), 50 + radius * math.sin(t))


def stroke(path, width):
    s = QPainterPathStroker()
    s.setWidth(width)
    s.setCapStyle(Qt.RoundCap)
    s.setJoinStyle(Qt.RoundJoin)
    return s.createStroke(path).simplified()


def circle_rect(radius):
    return QRectF(50 - radius, 50 - radius, 2 * radius, 2 * radius)


def mark(full=True):
    """[(outline, role)]: the ring and dot; with full=True also the promoter and the terminator."""
    ring = QPainterPath()
    ring.arcMoveTo(circle_rect(R), 90)                  # Qt angles: counter-clockwise from 3 o'clock
    ring.arcTo(circle_rect(R), 90, -GAP[0])
    body = stroke(ring, W)
    dot = QPainterPath()
    dot.addEllipse(QPointF(50, 50), W, W)
    parts = [(body.united(dot).simplified(), "ring")]
    if not full:
        return parts

    foot, top = R + W / 2 - 1, R + W / 2 + REACH        # the stems start just inside the ring's edge
    # Promoter: a stem out of the ring, then an arc along it, ending in an arrowhead pointing clockwise.
    span, head = 30.0, 6.0
    p = QPainterPath(pt(PROMOTER_AT, foot))
    p.lineTo(pt(PROMOTER_AT, top))
    p.arcTo(circle_rect(top), 90 - PROMOTER_AT, -span)
    end = PROMOTER_AT + span
    tip = QPainterPath(pt(end, top + 4))
    tip.lineTo(pt(end + math.degrees(head / top), top))
    tip.lineTo(pt(end, top - 4))
    tip.closeSubpath()
    tip = tip.united(stroke(tip, 1.5))                  # round the arrowhead's corners a little
    parts.append((stroke(p, LINE).united(tip).simplified(), "promoter"))

    # Terminator: a stem out of the ring with a bar across its end (a T).
    half = math.degrees(6.0 / top)
    t = QPainterPath(pt(TERMINATOR_AT, foot))
    t.lineTo(pt(TERMINATOR_AT, top))
    t.moveTo(pt(TERMINATOR_AT - half, top))
    t.lineTo(pt(TERMINATOR_AT + half, top))
    parts.append((stroke(t, LINE), "ink"))
    return parts


def bounds(parts):
    r = QRectF()
    for path, _ in parts:
        r = r.united(path.boundingRect())
    return r


def fit(parts, target, around=None):
    """Scale and move the parts so their real outline is centred in `target`, as large as fits.
    With `around` (a point of the parts), that point goes to the centre instead, still with everything inside."""
    b = bounds(parts)
    if around is not None:
        half = max(around.x() - b.left(), b.right() - around.x(), around.y() - b.top(), b.bottom() - around.y())
        b = QRectF(around.x() - half, around.y() - half, 2 * half, 2 * half)
    k = min(target.width() / b.width(), target.height() / b.height())
    tr = QTransform()
    tr.translate(target.center().x(), target.center().y())
    tr.scale(k, k)
    tr.translate(-b.center().x(), -b.center().y())
    return [(tr.map(path), role) for path, role in parts]


def load_fonts():
    """Registers the DM Sans files in packaging/fonts/ (the wordmark must not depend on the fonts installed here)."""
    files = sorted(FONTS.glob("DMSans-*.ttf"))
    if not files or any(QFontDatabase.addApplicationFont(str(f)) < 0 for f in files):
        sys.exit(f"make_icon: the DM Sans font files are missing or unreadable in {FONTS}")


def text_path(text, pixels, weight, spacing=0.0):
    """The outline of `text`, baseline at y = 0, starting near x = 0; `spacing` px are added after every character."""
    font = QFont(FONT_FAMILY)
    font.setPixelSize(round(pixels))
    font.setWeight(weight)
    font.setKerning(True)
    font.setLetterSpacing(QFont.AbsoluteSpacing, spacing)
    path = QPainterPath()
    path.addText(0, 0, font, text)
    if abs(pixels - round(pixels)) > 1e-6:               # fractional sizes: draw at the integer size, then scale
        path = QTransform.fromScale(pixels / round(pixels), pixels / round(pixels)).map(path)
    # The fonts are instances of a variable font, whose letters are built from overlapping contours: merge them,
    # or an even-odd fill punches holes where they overlap (the crossbar of the "A").
    path.setFillRule(Qt.WindingFill)
    return path.simplified()


def spread_text(text, pixels, weight, width):
    """`text` with the letter spacing chosen so its outline is exactly `width` wide (the size stays `pixels`)."""
    spacing = 0.0
    for _ in range(4):                                   # spacing in whole font pixels is not exact: converge on it
        path = text_path(text, pixels, weight, spacing)
        spacing += (width - path.boundingRect().width()) / (len(text) - 1) * round(pixels) / pixels
    return text_path(text, pixels, weight, spacing)


def place(path, left, baseline, scale=1.0):
    """`path` (baseline at 0) scaled by `scale`, its outline's left edge moved to `left` and its baseline to `baseline`."""
    tr = QTransform()
    tr.translate(left - path.boundingRect().left() * scale, baseline)
    tr.scale(scale, scale)
    return tr.map(path)


def wordmark():
    """[(outline, role)] of the oDNA wordmark: the mark as the "o", then "DNA"; "open" and the tagline below them.
    Baseline of "DNA" at y = 0."""
    size = 110.0                                          # "DNA"; everything else is in proportion to it
    k = 0.66                                              # the mark's scale: the ring a little smaller than a real "o"
    outer = (R + W / 2) * k                               # the ring's outer radius
    tr = QTransform()
    tr.translate(0, -outer)                               # the ring's outer edge on the baseline, its centre at x = 0
    tr.scale(k, k)
    tr.translate(-CENTRE.x(), -CENTRE.y())
    icon = [(tr.map(path), role) for path, role in mark()]

    dna = place(text_path("DNA", size, QFont.Bold), outer + 0.073 * size, 0)   # a letter's space after the ring
    db = dna.boundingRect()
    line = 0.255 * size                                   # the baseline of the line below
    word = text_path("open", 0.2 * size, QFont.Bold)
    word = place(word, -outer, line, 2 * outer / word.boundingRect().width())   # exactly as wide as the ring
    tag = place(spread_text("DESIGN · NAVIGATION · ASSEMBLY", 0.132 * size, QFont.Medium, db.width()), db.left(), line)
    return icon + [(dna, "ink"), (word, "ring"), (tag, "muted")]


def paint(painter, parts, colors):
    painter.setPen(Qt.NoPen)
    for path, role in parts:
        painter.setBrush(QColor(colors[role]))
        painter.drawPath(path)


def icon_image(size):
    """One square ARGB image: the tile, and the mark on it."""
    img = QImage(size, size, QImage.Format_ARGB32)
    img.fill(Qt.transparent)
    p = QPainter(img)
    p.setRenderHints(QPainter.Antialiasing)
    s = float(size)
    tile = QRectF(s * 0.03, s * 0.03, s * 0.94, s * 0.94)
    p.setBrush(QColor(TILE))
    p.setPen(QColor(TILE_EDGE) if size >= 32 else Qt.NoPen)
    p.drawRoundedRect(tile, s * 0.22, s * 0.22)
    full = size >= 32
    margin = s * (0.12 if full else 0.14)
    # The ring's centre is the tile's centre (the promoter and terminator stand out from it, inside the margin).
    paint(p, fit(mark(full), QRectF(margin, margin, s - 2 * margin, s - 2 * margin), CENTRE), LIGHT)
    p.end()
    return img


def wordmark_image(height, colors):
    parts = wordmark()
    b = bounds(parts)
    pad = 0.06 * b.height()
    k = height / (b.height() + 2 * pad)
    img = QImage(math.ceil((b.width() + 2 * pad) * k), height, QImage.Format_ARGB32)
    img.fill(Qt.transparent)
    p = QPainter(img)
    p.setRenderHints(QPainter.Antialiasing)
    p.scale(k, k)
    p.translate(pad - b.left(), pad - b.top())
    paint(p, parts, colors)
    p.end()
    return img


def svg_d(path):
    """SVG path data of a QPainterPath (move, line and cubic curve elements)."""
    out, i, n = [], 0, path.elementCount()
    f = lambda v: f"{v:.2f}".rstrip("0").rstrip(".")
    while i < n:
        e = path.elementAt(i)
        if e.type == QPainterPath.MoveToElement:
            out.append(f"M{f(e.x)} {f(e.y)}")
            i += 1
        elif e.type == QPainterPath.LineToElement:
            out.append(f"L{f(e.x)} {f(e.y)}")
            i += 1
        else:                                            # CurveTo, then its two CurveToData points
            c2, end = path.elementAt(i + 1), path.elementAt(i + 2)
            out.append(f"C{f(e.x)} {f(e.y)} {f(c2.x)} {f(c2.y)} {f(end.x)} {f(end.y)}")
            i += 3
    return "".join(out) + "Z"


def write_svg(target, parts, colors, tile=False, title=NAME):
    b = bounds(parts)
    pad = 0.06 * max(b.width(), b.height())
    x, y, w, h = (0, 0, 100, 100) if tile else (b.left() - pad, b.top() - pad, b.width() + 2 * pad, b.height() + 2 * pad)
    body = []
    if tile:
        body.append(f'<rect x="{x + w * 0.03:.2f}" y="{y + h * 0.03:.2f}" width="{w * 0.94:.2f}" height="{h * 0.94:.2f}" '
                    f'rx="{w * 0.22:.2f}" fill="{TILE}" stroke="{TILE_EDGE}"/>')
    for path, role in parts:
        body.append(f'<path fill="{colors[role]}" d="{svg_d(path)}"/>')
    target.write_text(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{x:.2f} {y:.2f} {w:.2f} {h:.2f}" role="img">'
                      f"<title>{title}</title>" + "".join(body) + "</svg>\n", encoding="utf-8")


def png_bytes(img):
    data = QByteArray()
    buf = QBuffer(data)
    buf.open(QIODevice.WriteOnly)
    img.save(buf, "PNG")
    return bytes(data)


def bmp_bytes(img):
    """ICO bitmap entry: a BITMAPINFOHEADER (height doubled), bottom-up BGRA pixels, then a 1-bit AND mask (all zero:
    the alpha channel does the transparency)."""
    w, h = img.width(), img.height()
    img = img.convertToFormat(QImage.Format_ARGB32)
    pixels = bytearray()
    for y in range(h - 1, -1, -1):
        for x in range(w):
            c = img.pixelColor(x, y)
            pixels += bytes((c.blue(), c.green(), c.red(), c.alpha()))
    mask_row = ((w + 31) // 32) * 4
    header = struct.pack("<IiiHHIIiiII", 40, w, h * 2, 1, 32, 0, len(pixels), 0, 0, 0, 0)
    return header + bytes(pixels) + bytes(mask_row * h)


def write_ico(path, images):
    """Sizes up to 64 px as classic bitmaps (the safest for old tools); 128 and 256 px as PNG (what Windows itself does)."""
    entries = [(size, png_bytes(img) if size >= 128 else bmp_bytes(img)) for size, img in images.items()]
    out = bytearray(struct.pack("<HHH", 0, 1, len(entries)))
    offset = 6 + 16 * len(entries)
    for size, data in entries:
        out += struct.pack("<BBBBHHII", size % 256, size % 256, 0, 0, 1, 32, len(data), offset)   # 0 means 256
        offset += len(data)
    for _, data in entries:
        out += data
    Path(path).write_bytes(bytes(out))


if __name__ == "__main__":
    app = QGuiApplication(sys.argv)                       # text outlines need the font system
    load_fonts()
    images = {size: icon_image(size) for size in SIZES}
    write_ico(ICO, images)
    LOGO_DIR.mkdir(exist_ok=True)
    for size, img in images.items():
        img.save(str(LOGO_DIR / f"icon_{size}.png"))
    wordmark_image(160, LIGHT).save(str(LOGO_DIR / "wordmark.png"))
    wordmark_image(160, DARK).save(str(LOGO_DIR / "wordmark_dark.png"))
    write_svg(DOCS / "icon.svg", fit(mark(), QRectF(12, 12, 76, 76), CENTRE), LIGHT, tile=True)
    write_svg(DOCS / "wordmark.svg", wordmark(), LIGHT)
    write_svg(DOCS / "wordmark-dark.svg", wordmark(), DARK)
    print(f"wrote {ICO.name} ({', '.join(map(str, SIZES))} px), app/logo/*.png, docs/images/*.svg")
