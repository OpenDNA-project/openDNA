"""The 3D view of the Protein panel: the structure of the protein shown, drawn here (no browser), with the residues of the selection
highlighted and a click on a residue selecting its codons.

The backbone is a smooth tube through the alpha carbons (thick in helices, medium in strands: see structure.secondary), shaded
across so it looks round; the highlighted residues are thicker and show their atoms while the rest fades, ligands are small balls. QPainter draws it, far parts first and fogged, so it needs no
OpenGL (it runs in the offscreen tests and on any graphics driver).
shortcut: QPainter redraws everything on each frame; past ~10 000 residues turning gets slow, move to QOpenGLWidget if that matters.

Mouse: drag to turn, right or middle drag to move, wheel to zoom, double-click to fit, click a residue to select it (Shift: up to
it), hover for its name. Keys: arrows turn, + and - zoom, Home fits.
"""
import math
from html import escape

from PySide6.QtCore import QPointF, Qt, Signal
from PySide6.QtGui import QColor, QIcon, QLinearGradient, QPainter, QPalette, QPen, QPixmap
from PySide6.QtWidgets import (
    QComboBox, QGridLayout, QHBoxLayout, QLabel, QListWidget, QListWidgetItem, QSizePolicy, QToolButton, QToolTip, QVBoxLayout,
    QWidget,
)

import structure

STEPS = 4                                   # spline points per residue
WIDTH = {"H": 1.7, "E": 1.3, "-": 0.55}     # tube width (Å) in helices, strands and the rest
SS_COLOURS = {"H": "#D6457A", "E": "#E8B62C", "-": "#9AA0A6"}
CHAIN_COLOURS = ["#4E79A7", "#F28E2B", "#59A14F", "#B07AA1", "#76B7B2", "#EDC948", "#FF9DA7", "#9C755F"]
ELEMENT_COLOURS = {"N": "#3050F8", "O": "#E03030", "S": "#D8B020", "P": "#FF8000", "FE": "#C06030"}
# Colours by the kind of amino acid: the tracks under the letters of the Protein panel use them too (protein_tracks.py).
KYTE_DOOLITTLE = {"I": 4.5, "V": 4.2, "L": 3.8, "F": 2.8, "C": 2.5, "M": 1.9, "A": 1.8, "G": -0.4, "T": -0.7, "S": -0.8, "W": -0.9,
                  "Y": -1.3, "P": -1.6, "H": -3.2, "E": -3.5, "Q": -3.5, "D": -3.5, "N": -3.5, "K": -3.9, "R": -4.5}
VOLUME = {"G": 60.1, "A": 88.6, "S": 89.0, "C": 108.5, "D": 111.1, "P": 112.7, "N": 114.1, "T": 116.1, "E": 138.4, "V": 140.0,
          "Q": 143.8, "H": 153.2, "M": 162.9, "I": 166.7, "L": 166.7, "K": 168.6, "R": 173.4, "F": 189.9, "Y": 193.6, "W": 227.8}   # Å³, Zamyatnin 1972
PROPERTY = {**dict.fromkeys("ILVAM", "#FFAFAF"), **dict.fromkeys("FWY", "#FFC800"), **dict.fromkeys("KRH", "#7B7BFF"),
            **dict.fromkeys("DE", "#FF5050"), **dict.fromkeys("STNQ", "#50E050"), **dict.fromkeys("PG", "#FF70FF"), "C": "#FFFF00"}   # Zappo
CHARGE = {"K": "#6C8EF5", "R": "#6C8EF5", "H": "#B4C7F9", "D": "#F07070", "E": "#F07070"}
AA_SCHEMES = {   # name -> what the colours mean (the tooltips)
    "Property": "Zappo colours: aliphatic (ILVAM) salmon, aromatic (FWY) orange, positive (KRH) blue, negative (DE) red, "
                "polar (STNQ) green, P and G pink, C yellow.",
    "Hydrophobicity": "Kyte-Doolittle scale: hydrophobic red (I, V, L), hydrophilic blue (R, K, D, E), in between white.",
    "Size": "Volume of the side chain (Zamyatnin): small pale (G, A, S), large brown (W, Y, F).",
    "Charge": "At pH 7: positive blue (K, R; H pale), negative red (D, E), the others uncoloured.",
}
MODES = ["Confidence", "Rainbow", "Chain", "Secondary structure", *AA_SCHEMES]


def _mix(a, b, t):
    a, b = QColor(a), QColor(b)
    return QColor.fromRgbF(*(x + (y - x) * t for x, y in zip(a.getRgbF()[:3], b.getRgbF()[:3]))).name()


def aa_colour(scheme, aa):
    """The colour (#rrggbb) of the amino acid `aa` in `scheme` (a key of AA_SCHEMES), or None (no colour: X, a stop, an uncharged one)."""
    if scheme == "Property":
        return PROPERTY.get(aa)
    if scheme == "Charge":
        return CHARGE.get(aa)
    if scheme == "Hydrophobicity" and aa in KYTE_DOOLITTLE:
        t = (KYTE_DOOLITTLE[aa] + 4.5) / 9
        return _mix("#5B9BE0", "#FFFFFF", t * 2) if t < 0.5 else _mix("#FFFFFF", "#E8544E", t * 2 - 1)
    if scheme == "Size" and aa in VOLUME:
        return _mix("#FFF6D5", "#A0581C", (VOLUME[aa] - 60.1) / (227.8 - 60.1))
    return None
PICK_RADIUS = 10                            # px


def _matmul(a, b):
    return [[sum(a[i][k] * b[k][j] for k in range(3)) for j in range(3)] for i in range(3)]


def _orthonormal(m):
    """The rotation `m` with the rounding drift of many small turns taken out (Gram-Schmidt on its rows)."""
    x = m[0]
    n = math.sqrt(sum(c * c for c in x))
    x = [c / n for c in x]
    y = m[1]
    d = sum(a * b for a, b in zip(x, y))
    y = [b - d * a for a, b in zip(x, y)]
    n = math.sqrt(sum(c * c for c in y))
    y = [c / n for c in y]
    return [x, y, [x[1] * y[2] - x[2] * y[1], x[2] * y[0] - x[0] * y[2], x[0] * y[1] - x[1] * y[0]]]


def _spline(points):
    """Catmull-Rom points through `points` (STEPS per gap) and, for each, the index in `points` of the nearest one."""
    if len(points) == 1:
        return list(points), [0]
    out, owner = [], []
    for i in range(len(points) - 1):
        p0, p1, p2 = points[max(i - 1, 0)], points[i], points[i + 1]
        p3 = points[min(i + 2, len(points) - 1)]
        for s in range(STEPS):
            t = s / STEPS
            t2, t3 = t * t, t * t * t
            out.append(tuple(0.5 * (2 * b + (c - a) * t + (2 * a - 5 * b + 4 * c - d) * t2 + (3 * b - a - 3 * c + d) * t3)
                             for a, b, c, d in zip(p0, p1, p2, p3)))
            owner.append(i if t < 0.5 else i + 1)
    out.append(points[-1])
    owner.append(len(points) - 1)
    return out, owner


class StructureCanvas(QWidget):
    picked = Signal(int, bool)        # index in the protein, Shift held

    def __init__(self):
        super().__init__()
        self.setMinimumHeight(160)
        self.setMouseTracking(True)
        self.setFocusPolicy(Qt.StrongFocus)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.residues, self.ligands, self.to_model, self.to_protein = [], [], {}, {}
        self.highlight, self.mode, self.plddt = set(), "Rainbow", False
        self.hidden = set()   # chains ("A") and ligand groups ("HEM A900") not drawn
        self.points, self.owner, self.links, self.cas, self.atoms = [], [], [], [], {}
        self.reset_view()

    # ---- the model ----------------------------------------------------------------------------------------------------------------
    def set_structure(self, found, protein):
        """Show `found` (a structure.Structure) for the amino acids `protein`."""
        self.residues, self.ligands = structure.model_of(found.pdb)
        self.plddt = found.plddt
        self.to_model = structure.residue_map(protein, self.residues)
        self.to_protein = {k: i for i, ks in self.to_model.items() for k in ks}
        centre = [sum(r.ca[n] for r in self.residues) / max(len(self.residues), 1) for n in range(3)]
        self.points, self.owner, self.links = [], [], []
        for run in structure.segments(self.residues):
            pts, owner = _spline([self.residues[k].ca for k in run])
            first = len(self.points)
            self.points += [tuple(c - m for c, m in zip(p, centre)) for p in pts]
            self.owner += [run[o] for o in owner]
            self.links += [(j, j + 1) for j in range(first, len(self.points) - 1)]
        self.cas = [tuple(c - m for c, m in zip(r.ca, centre)) for r in self.residues]
        self.atoms = {k: [(name, el, tuple(c - m for c, m in zip(xyz, centre))) for name, el, xyz in r.atoms]
                      for k, r in enumerate(self.residues)}
        self.ligands = [(el, tuple(c - m for c, m in zip(xyz, centre)), group) for el, xyz, group in self.ligands]
        self.hidden = set()
        # fit the folded part: a floppy tail (a signal peptide in a prediction, a tenth of the chain) would leave it tiny, so the
        # radius is that of 85 % of the backbone, plus a little
        far = sorted(math.sqrt(sum(c * c for c in p)) for p in self.points)
        self.radius = (1.15 * far[int(0.85 * (len(far) - 1))] if far else 1.0) or 1.0
        self.mode = "Confidence" if found.plddt else "Rainbow"
        self.highlight = set()
        self.reset_view()

    def set_highlight(self, indices):
        """Highlight the residues of the protein at `indices` (and show their atoms)."""
        wanted = {k for i in indices for k in self.to_model.get(i, ())}
        if len(wanted) == len(self.to_protein):   # the whole protein (its CDS is selected): nothing to point out
            wanted = set()
        if wanted != self.highlight:
            self.highlight = wanted
            self.update()

    def set_shown(self, key, shown):
        """Show or hide the chain or ligand group `key`."""
        (self.hidden.discard if shown else self.hidden.add)(key)
        self.update()

    def set_mode(self, mode):
        self.mode = mode
        self.update()

    def colour(self, k):
        r = self.residues[k]
        if self.mode == "Confidence" and self.plddt:
            return next((c for lo, c, _ in structure.PLDDT_BANDS if r.b >= lo), structure.PLDDT_BANDS[-1][1])
        if self.mode == "Chain":
            chains = sorted({x.chain for x in self.residues})
            return CHAIN_COLOURS[chains.index(r.chain) % len(CHAIN_COLOURS)]
        if self.mode == "Secondary structure":
            return SS_COLOURS[r.ss]
        if self.mode in AA_SCHEMES:
            return aa_colour(self.mode, r.aa) or "#C8C8C8"
        return QColor.fromHsvF(0.7 * (1 - k / max(len(self.residues) - 1, 1)), 0.7, 0.92).name()   # blue (N) to red (C)

    # ---- the view -----------------------------------------------------------------------------------------------------------------
    def reset_view(self):
        self.rotation, self.zoom, self.pan = [[1, 0, 0], [0, 1, 0], [0, 0, 1]], 1.0, QPointF(0, 0)
        self.update()

    def turn(self, ax, ay):
        """Turn about the screen's horizontal axis by `ax` and its vertical axis by `ay` (radians)."""
        ca, sa, cb, sb = math.cos(ax), math.sin(ax), math.cos(ay), math.sin(ay)
        rx = [[1, 0, 0], [0, ca, -sa], [0, sa, ca]]
        ry = [[cb, 0, sb], [0, 1, 0], [-sb, 0, cb]]
        self.rotation = _orthonormal(_matmul(ry, _matmul(rx, self.rotation)))
        self.update()

    def scale(self):
        return 0.45 * min(self.width(), self.height()) / self.radius * self.zoom if self.points else 1.0

    def project(self, p):
        """Screen x, y and depth (larger is nearer) of the model point `p`."""
        (a, b, c), (d, e, f), (g, h, i) = self.rotation
        s = self.scale()
        x, y, z = p
        return (self.width() / 2 + self.pan.x() + s * (a * x + b * y + c * z),
                self.height() / 2 + self.pan.y() - s * (d * x + e * y + f * z), g * x + h * y + i * z)

    def residue_at(self, pos):
        """Index in `residues` of the alpha carbon under `pos` (the nearest to the viewer), or None."""
        best = None
        for k, ca in enumerate(self.cas):
            if self.residues[k].chain in self.hidden:
                continue
            x, y, z = self.project(ca)
            if (x - pos.x()) ** 2 + (y - pos.y()) ** 2 <= PICK_RADIUS ** 2 and (best is None or z > best[0]):
                best = (z, k)
        return best and best[1]

    def paintEvent(self, event):
        p = QPainter(self)
        background = self.palette().color(QPalette.Base)
        p.fillRect(self.rect(), background)
        if not self.points:
            p.setPen(self.palette().color(QPalette.PlaceholderText))
            p.drawText(self.rect(), Qt.AlignCenter, "No structure")
            return
        p.setRenderHint(QPainter.Antialiasing)
        s, far = self.scale(), self.radius
        bg = (background.redF(), background.greenF(), background.blueF())
        fogged = {}

        def fog(hex_colour, depth, faded=False):
            """The colour mixed with the background as it gets further (8 steps, cached); `faded`: much more (not highlighted)."""
            step = max(0, min(7, int((depth + far) / (2 * far) * 8)))
            key = (hex_colour, step, faded)
            if key not in fogged:
                c = QColor(hex_colour)
                f = (0.35 + 0.65 * step / 7) * (0.45 if faded else 1)   # the farthest parts are 65 % background
                c.setRgbF(*(f * a + (1 - f) * b for a, b in zip((c.redF(), c.greenF(), c.blueF()), bg)))
                fogged[key] = (c, c.lighter(150), c.darker(160))
            return fogged[key]

        screen = [self.project(q) for q in self.points]
        colours = {k: self.colour(k) for k in set(self.owner)}
        dim = bool(self.highlight)
        items = []   # (depth, kind, ...): drawn far to near
        for a, b in self.links:
            k = self.owner[a]
            if self.residues[k].chain in self.hidden:
                continue
            w = WIDTH[self.residues[k].ss] * (1.4 if k in self.highlight else 1.0)
            items.append(((screen[a][2] + screen[b][2]) / 2, 0, a, b, colours[k], w * s, dim and k not in self.highlight))
        for k in self.highlight:
            if self.residues[k].chain in self.hidden:
                continue
            atoms = self.atoms[k]
            seen = [self.project(xyz) for _, _, xyz in atoms]
            for i in range(len(atoms)):
                for j in range(i + 1, len(atoms)):
                    if math.dist(atoms[i][2], atoms[j][2]) <= 1.9:   # bonded
                        items.append(((seen[i][2] + seen[j][2]) / 2, 1, seen[i], seen[j], "#808080", 0.3 * s))
            for (_, el, _), pa in zip(atoms, seen):
                items.append((pa[2] + 0.01, 2, pa, ELEMENT_COLOURS.get(el, colours.get(k, "#909090")), 0.35 * s))
        for el, xyz, group in self.ligands:
            if group in self.hidden:
                continue
            pa = self.project(xyz)
            items.append((pa[2], 2, pa, ELEMENT_COLOURS.get(el, "#5FAE57"), 0.45 * s, dim))
        items.sort(key=lambda item: item[0])
        for item in items:
            depth, kind = item[0], item[1]
            if kind == 0:
                _, _, a, b, colour, w, faded = item
                (x1, y1, _), (x2, y2, _) = screen[a], screen[b]
                mid, light, dark = fog(colour, depth, faded)
                # shaded across the tube, lit from the top left: dark edge, light band, colour, dark edge
                length = math.hypot(x2 - x1, y2 - y1) or 1.0
                nx, ny = (y1 - y2) / length * w / 2, (x2 - x1) / length * w / 2
                if nx + ny > 0:
                    nx, ny = -nx, -ny
                shade = QLinearGradient(x1 + nx, y1 + ny, x1 - nx, y1 - ny)
                shade.setColorAt(0, mid)
                shade.setColorAt(0.3, light)
                shade.setColorAt(0.65, mid)
                shade.setColorAt(1, dark)
                p.setPen(QPen(shade, max(w, 1.0), Qt.SolidLine, Qt.RoundCap))
                p.drawLine(QPointF(x1, y1), QPointF(x2, y2))
            elif kind == 1:
                _, _, pa, pb, colour, w = item
                p.setPen(QPen(fog(colour, depth)[0], max(w, 1.0), Qt.SolidLine, Qt.RoundCap))
                p.drawLine(QPointF(pa[0], pa[1]), QPointF(pb[0], pb[1]))
            else:
                _, _, pa, colour, r, *faded = item
                p.setPen(Qt.NoPen)
                p.setBrush(fog(colour, depth, bool(faded and faded[0]))[0])
                p.drawEllipse(QPointF(pa[0], pa[1]), r, r)
        p.end()

    # ---- mouse and keys -----------------------------------------------------------------------------------------------------------
    def mousePressEvent(self, event):
        self.press = self.last = event.position()
        self.moved = False

    def mouseMoveEvent(self, event):
        pos = event.position()
        if event.buttons() == Qt.NoButton:
            k = self.residue_at(pos) if self.residues else None
            if k is None:
                QToolTip.hideText()
            else:
                r, i = self.residues[k], self.to_protein.get(k)
                text = (f"{r.aa}{r.number} · chain {r.chain}" + (f" · pLDDT {r.b:.0f}" if self.plddt else "")
                        + (f" · aa {i + 1} of the protein" if i is not None else " · not part of this protein"))
                QToolTip.showText(event.globalPosition().toPoint(), text, self)
            return
        d = pos - self.last
        self.last = pos
        if (pos - self.press).manhattanLength() > 3:
            self.moved = True
        if event.buttons() & Qt.LeftButton:
            self.turn(d.y() * 0.01, d.x() * 0.01)
        else:
            self.pan += d
            self.update()

    def mouseReleaseEvent(self, event):
        if event.button() == Qt.LeftButton and not self.moved and self.residues:
            k = self.residue_at(event.position())
            if k is not None and k in self.to_protein:
                self.picked.emit(self.to_protein[k], bool(event.modifiers() & Qt.ShiftModifier))

    def mouseDoubleClickEvent(self, event):
        self.reset_view()

    def wheelEvent(self, event):
        self.zoom = max(0.2, min(20.0, self.zoom * 1.0015 ** event.angleDelta().y()))
        self.update()

    def keyPressEvent(self, event):
        key, step = event.key(), 0.12
        turns = {Qt.Key_Up: (-step, 0), Qt.Key_Down: (step, 0), Qt.Key_Left: (0, -step), Qt.Key_Right: (0, step)}
        if key in turns:
            self.turn(*turns[key])
        elif key in (Qt.Key_Plus, Qt.Key_Equal, Qt.Key_Minus):
            self.zoom = max(0.2, min(20.0, self.zoom * (1.2 if key != Qt.Key_Minus else 1 / 1.2)))
            self.update()
        elif key == Qt.Key_Home:
            self.reset_view()
        else:
            super().keyPressEvent(event)


class StructureView(QWidget):
    """The canvas under what the database says about the structure (its name, then its fields two by two: source, entry with a link to
    its page, method, resolution, organism...), with the colouring and a button that hides it."""
    picked = Signal(int, bool)
    closed = Signal()

    def __init__(self):
        super().__init__()
        self.canvas = StructureCanvas()
        self.canvas.picked.connect(self.picked)
        self.title = QLabel(wordWrap=True, textInteractionFlags=Qt.TextSelectableByMouse)
        self.title.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Preferred)
        font = self.title.font()
        font.setBold(True)
        self.title.setFont(font)
        self.mode = QComboBox(toolTip="Colour the structure by", sizeAdjustPolicy=QComboBox.AdjustToContents)
        self.mode.currentTextChanged.connect(lambda m: self.mode.setToolTip(AA_SCHEMES.get(m, "Colour the structure by")))
        self.mode.currentTextChanged.connect(self.canvas.set_mode)
        close = QToolButton(text="✕", toolTip="Hide the 3D view", autoRaise=True, clicked=self._close)
        top = QHBoxLayout()
        top.setContentsMargins(0, 0, 0, 0)
        top.addWidget(self.title, 1)
        top.addWidget(self.mode, 0, Qt.AlignTop)
        top.addWidget(close, 0, Qt.AlignTop)
        self.fields = QGridLayout()
        self.fields.setContentsMargins(0, 0, 0, 4)
        self.fields.setHorizontalSpacing(8)
        self.fields.setVerticalSpacing(1)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addLayout(top)
        layout.addLayout(self.fields)
        # what the model holds, a click on a row shows or hides it: the chains, then the other molecules (cofactors, ions...)
        self.parts = QListWidget(toolTip="Click to show or hide", selectionMode=QListWidget.NoSelection)
        self.parts.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Expanding)
        self.parts.setMaximumWidth(150)
        self.parts.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)   # long names are cut short (whole in the tooltip)
        self.parts.itemClicked.connect(self._toggle)
        middle = QHBoxLayout()
        middle.setContentsMargins(0, 0, 0, 0)
        middle.addWidget(self.canvas, 1)
        middle.addWidget(self.parts)
        layout.addLayout(middle, 1)
        self.protein = None   # the amino acids shown
        self.shown_fields = {}   # field -> text, as shown

    def set_structure(self, found, protein):
        self.protein = protein
        self.canvas.set_structure(found, protein)
        self.mode.blockSignals(True)
        self.mode.clear()
        self.mode.addItems(MODES if found.plddt else MODES[1:])
        self.mode.setCurrentText(self.canvas.mode)
        self.mode.blockSignals(False)
        info = dict(found.info)
        link = info.pop("Link", "")
        self.title.setText(info.pop("Protein", "") or (found.detail if not info else found.source))
        if not info:   # a structure cached before the fields existed: its one line says it all
            info = {"Source": found.source}
        mapped = len(self.canvas.to_model)
        info["In the model"] = f"{mapped:,} of {len(protein):,} aa" + (f" ({mapped / len(protein):.0%})" if mapped < len(protein) else "")
        self.shown_fields = info
        while self.fields.count():
            self.fields.takeAt(0).widget().deleteLater()
        muted = self.palette().color(QPalette.PlaceholderText).name()
        for n, (key, value) in enumerate(info.items()):
            name = QLabel(f"<span style='color:{muted}'>{escape(key)}</span>")
            text = escape(value)
            if key == "Entry" and link.startswith("https://"):
                text = f"<a href='{escape(link)}'>{text}</a>"   # opens the database's page in the browser
            shown = QLabel(text, openExternalLinks=True, toolTip=value,   # one line (the whole text in the tooltip if it is cut)
                           textInteractionFlags=Qt.TextBrowserInteraction | Qt.TextSelectableByMouse)
            shown.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Preferred)
            row, column = divmod(n, 2)
            self.fields.addWidget(name, row, column * 2, Qt.AlignRight | Qt.AlignTop)
            self.fields.addWidget(shown, row, column * 2 + 1, Qt.AlignLeft | Qt.AlignTop)   # (top: in line with its name)
        self.fields.setColumnStretch(1, 1)
        self.fields.setColumnStretch(3, 1)
        self.parts.clear()
        residues = self.canvas.residues
        chains = sorted({r.chain for r in residues})
        for n, c in enumerate(chains):
            swatch = QPixmap(10, 10)
            swatch.fill(QColor(CHAIN_COLOURS[n % len(CHAIN_COLOURS)]))   # the colours of the "Chain" colouring
            self._add_part(f"Chain {c} · {sum(r.chain == c for r in residues):,} aa", c, QIcon(swatch))
        for group in dict.fromkeys(g for _, _, g in self.canvas.ligands):
            self._add_part(group, group)
        self.parts.setVisible(self.parts.count() > 1)   # one chain alone: nothing to choose

    def _add_part(self, text, key, icon=None):
        item = QListWidgetItem(text, self.parts)
        item.setToolTip(f"{text}: click to show or hide")
        if icon:
            item.setIcon(icon)
        item.setData(Qt.UserRole, key)
        item.setFlags(Qt.ItemIsEnabled)   # not user-checkable: the click anywhere on the row toggles it (_toggle), once
        item.setCheckState(Qt.Checked)

    def _toggle(self, item):
        shown = item.checkState() != Qt.Checked
        item.setCheckState(Qt.Checked if shown else Qt.Unchecked)
        self.canvas.set_shown(item.data(Qt.UserRole), shown)

    def set_highlight(self, indices):
        self.canvas.set_highlight(indices)

    def _close(self):
        self.hide()
        self.closed.emit()
