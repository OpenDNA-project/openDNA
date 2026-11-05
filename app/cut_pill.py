"""The "pill" of a restriction enzyme: a small picture, shown in the hover tooltip of its cut, of the recognition site as a
double-stranded piece of DNA that has really been cut. Both strands are written out with the base pairs between them; the
two pieces are pulled apart at the cut and the strands end where the enzyme cuts them, so a 5' overhang, a 3' overhang or
a blunt end can be seen at a glance:

    EcoRI                         5'-…G     AATTC…-3'
    G^AATTC                       3'-…CTTAA     G…-5'

Nothing is drawn on the sequence: this is tooltip content only (an image inside the rich text of the tooltip).

Two parts: `pill_layout` (no drawing: the columns of the picture and which piece each strand's letter belongs to, which the
tests read) and `pill_pixmap` / `pill_html` (QPainter on a QPixmap, in the theme's colours, at the screen's pixel density and
at least twice the logical one, so it is crisp). `enzyme_rich_tip` is the whole tooltip: the pill and one line of words.
"""
import math
from collections import namedtuple
from types import SimpleNamespace

from PySide6.QtCore import QBuffer, QByteArray, QIODevice, QRectF, Qt, QPointF
from PySide6.QtGui import QFont, QFontMetricsF, QGuiApplication, QPainter, QPen, QPixmap

import theme
from cut_sites import _after, _distance, locate_site
from editor import mono_font
from enzymes import notation

# One column of the picture: the base on the top strand and its complement, the piece (0 = left of the cut, 1 = right) each
# strand's letter belongs to, whether it is part of the recognition site, whether the two letters are still paired in
# their piece (a tick between them), and whether it is a real base (False: a "…" marker or a blank past the end of a
# linear sequence).
Column = namedtuple("Column", "top bottom tpiece bpiece site paired real")

_COMP = str.maketrans("ACGTURYSWKMBDHVN", "TGCAAYRSWMKVHDBN")
FLANK_GAP = 2      # at most this many hidden bases are written out; more are replaced by a "…" (Type IIS: site far from its cut)
DPR = 2            # the picture is drawn at least at this many pixels per logical pixel (more on a denser screen: see density)
_cache = {}


def pill_layout(enzyme, seq, circular, top, bottom):
    """The columns of the picture of `enzyme` cutting `seq` at (`top`, `bottom`) (gaps between bases, as find_cuts gives
    them). Returns a namespace: `columns`, `kind` ("blunt", "5′" or "3′" overhang), `overhang` (its bases, read on the top
    strand), `site` ((start, end, strand) of the recognition site or None), `ut` / `ub` (the cuts, counted on from the site so
    that a cut across the origin of a circle stays right), `pieces` (the letters of the four pieces: left top, left bottom,
    right top, right bottom, markers left out)."""
    n = len(seq)
    found = locate_site(enzyme, seq, circular, top, bottom)
    if found:
        s, e, strand, ut, ub = found
        site = (s, e, strand)
    else:
        site, ut, ub, s, e = None, top, bottom, 0, 0
    lo, hi = min(ut, ub), max(ut, ub)
    # The bases drawn: the site and the cut with one base on each side of it. If the cut is far from the site (Type IIS)
    # the bases between are left out and marked "…".
    wanted = set(range(lo - 1, hi + 1)) | (set(range(s, e)) if site else set(range(lo - 3, hi + 3)))
    idx = sorted(wanted)

    def letter(i):
        if circular and n:
            return seq[i % n].upper()
        return seq[i].upper() if 0 <= i < n else None

    def base(i):
        c = letter(i)
        in_site = site is not None and s <= i < e
        tp, bp = int(i >= ut), int(i >= ub)
        if c is None:
            return Column(" ", " ", tp, bp, in_site, False, False)
        return Column(c, c.translate(_COMP), tp, bp, in_site, tp == bp, True)

    def marker(i):   # a "…" standing for the bases around index i
        return Column("…", "…", int(i >= ut), int(i >= ub), False, False, False)

    cols = [marker(idx[0] - 1)]
    prev = None
    for i in idx:
        if prev is not None and i - prev > 1:
            if i - prev - 1 <= FLANK_GAP:
                cols += [base(j) for j in range(prev + 1, i)]
            else:
                cols.append(marker(prev + 1))
        cols.append(base(i))
        prev = i
    cols.append(marker(idx[-1] + 1))
    kind = "blunt" if ut == ub else ("5′" if ut < ub else "3′")
    overhang = "".join(letter(i) or "" for i in range(lo, hi)) if ut != ub else ""
    real = [c for c in cols if c.real]
    pieces = tuple("".join(c.top if which == "top" else c.bottom for c in real if (c.tpiece if which == "top" else c.bpiece) == piece)
                   for piece, which in ((0, "top"), (0, "bottom"), (1, "top"), (1, "bottom")))
    return SimpleNamespace(columns=cols, kind=kind, overhang=overhang, site=site, ut=ut, ub=ub, pieces=pieces,
                           name=enzyme.name, notation=notation(enzyme))


# ---- drawing -------------------------------------------------------------------------------------------------------------


def density():
    """Pixels per logical pixel to draw the picture at: the screen's (rounded up: a 250 % screen gets 3), at least DPR. A fixed 2
    was blurred on a 300 % screen."""
    screen = QGuiApplication.primaryScreen()
    return max(DPR, math.ceil(screen.devicePixelRatio())) if screen is not None else DPR


def pill_pixmap(layout, dpr=None):
    """The picture of `layout` as a QPixmap (logical size = pixmap size / devicePixelRatio), in the current theme's colours,
    drawn at `dpr` pixels per logical pixel (default: density())."""
    dpr = dpr or density()
    bg, fg, edge = theme.color("pill_bg"), theme.color("pill_fg"), theme.color("pill_edge")
    band, tick, cut = theme.color("pill_site"), theme.color("pill_tick"), theme.color("cut")
    font = mono_font(11)
    bold = QFont(font)
    bold.setBold(True)
    name_font = QFont(QGuiApplication.font())
    name_font.setBold(True)
    name_font.setPointSizeF(10)
    small = QFont(QGuiApplication.font())
    small.setPointSizeF(9)
    fm, nfm, sfm = QFontMetricsF(font), QFontMetricsF(name_font), QFontMetricsF(small)
    cw, lh = fm.horizontalAdvance("A") + 4, fm.height() + 2     # a cell: one letter of one strand
    gap, tick_h, pad = cw * 1.1, 7.0, 11.0                      # how far the pieces are pulled apart; the pairing row
    cols = layout.columns
    end_l, end_r = fm.horizontalAdvance("5′-") + 3, fm.horizontalAdvance("-3′") + 3
    body = end_l + len(cols) * cw + gap + end_r
    head = nfm.horizontalAdvance(layout.name) + 12 + sfm.horizontalAdvance(layout.notation)
    width = max(body, head) + 2 * pad
    x0 = (width - body) / 2 + end_l                              # x of the first column
    name_y = pad - 2
    top_y = name_y + nfm.height() + 5
    bot_y = top_y + lh + tick_h
    height = bot_y + lh + pad

    pm = QPixmap(int(width * dpr + 1), int(height * dpr + 1))
    pm.setDevicePixelRatio(dpr)
    pm.fill(Qt.transparent)
    p = QPainter(pm)
    p.setRenderHint(QPainter.Antialiasing)
    p.setRenderHint(QPainter.TextAntialiasing)
    p.setPen(QPen(edge, 1.2))
    p.setBrush(bg)
    p.drawRoundedRect(QRectF(0.7, 0.7, width - 1.4, height - 1.4), 15, 15)   # the pill

    p.setFont(name_font)
    p.setPen(fg)
    p.drawText(QPointF(pad, name_y + nfm.ascent()), layout.name)
    p.setFont(small)
    p.setPen(tick)
    p.drawText(QPointF(pad + nfm.horizontalAdvance(layout.name) + 12, name_y + nfm.ascent()), layout.notation)

    x = lambda c, piece: x0 + c * cw + gap * piece
    for strand, y in (("top", top_y), ("bottom", bot_y)):
        piece_of = (lambda col: col.tpiece) if strand == "top" else (lambda col: col.bpiece)
        for c, col in enumerate(cols):
            cell = QRectF(x(c, piece_of(col)), y, cw, lh)
            if col.site and col.real:
                p.fillRect(cell, band)
            overhang = col.real and col.tpiece != col.bpiece
            p.setFont(bold if overhang or col.site else font)
            p.setPen(cut if overhang else (tick if not col.real else fg))
            p.drawText(cell, Qt.AlignCenter, col.top if strand == "top" else col.bottom)
        # the 5' / 3' ends of the two pieces of this strand
        first_right = next((c for c, col in enumerate(cols) if piece_of(col) == 1), len(cols))
        p.setFont(small)
        p.setPen(fg)
        left, right = ("5′-", "-3′") if strand == "top" else ("3′-", "-5′")
        p.drawText(QRectF(x0 - end_l, y, end_l - 1, lh), Qt.AlignRight | Qt.AlignVCenter, left)
        p.drawText(QRectF(x(len(cols), 1), y, end_r, lh), Qt.AlignLeft | Qt.AlignVCenter, right)
        # the break: where this strand is cut, a red stroke in the gap between the two pieces
        bx = x0 + first_right * cw + gap / 2
        p.setPen(QPen(cut, 2.2, Qt.SolidLine, Qt.RoundCap))
        p.drawLine(QPointF(bx, y - 1), QPointF(bx, y + lh + 1))
    # the pairing ticks between the letters that are still paired
    p.setPen(QPen(tick, 1.2))
    for c, col in enumerate(cols):
        if col.paired:
            cx = x(c, col.tpiece) + cw / 2
            p.drawLine(QPointF(cx, top_y + lh), QPointF(cx, top_y + lh + tick_h))
    p.end()
    return pm


def pill_html(layout):
    """The picture as an <img> for rich text (a data: URI, so a tooltip needs no file)."""
    dpr = density()
    key = (tuple(layout.columns), layout.name, layout.notation, theme.is_dark(), dpr)
    if key not in _cache:
        pm = pill_pixmap(layout, dpr)
        data = QByteArray()
        buf = QBuffer(data)
        buf.open(QIODevice.WriteOnly)
        pm.save(buf, "PNG")
        b64 = bytes(data.toBase64()).decode("ascii")
        if len(_cache) > 64:
            _cache.clear()
        _cache[key] = (f'<img src="data:image/png;base64,{b64}" width="{round(pm.width() / dpr)}" '
                       f'height="{round(pm.height() / dpr)}">')
    return _cache[key]


def summary(enzyme, layout, n, circular, top, bottom):
    """One line under the picture: the site with its cut, the ends, where each strand is cut."""
    parts = [notation(enzyme)]
    if layout.kind == "blunt":
        parts.append("blunt ends")
    else:
        parts.append(f"{layout.kind} overhang {abs(layout.ub - layout.ut)} nt ({layout.overhang})")
    if layout.site is not None:
        s, e, _ = layout.site
        if min(layout.ut, layout.ub) >= e or max(layout.ut, layout.ub) <= s:
            parts.append(f"Type IIS, cuts {_distance(layout.ut, s, e)} (top strand)")
        parts.append(f"site at bases {s % n + 1:,}–{(e - 1) % n + 1:,}")
    parts.append(f"top strand cut {_after(top, n, circular)}, bottom strand {_after(bottom, n, circular)}")
    return " · ".join(parts)


def enzyme_rich_tip(enzyme, seq, circular, top, bottom):
    """The tooltip (rich text) of a restriction cut: the pill and one line of words."""
    layout = pill_layout(enzyme, seq, circular, top, bottom)
    return pill_html(layout) + "<br>" + summary(enzyme, layout, len(seq), circular, top, bottom)
