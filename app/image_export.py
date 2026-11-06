"""Export > Image of the Current Tab: every pane of the tab shown, as they are laid out on screen,
written as PNG or TIFF (raster, at the resolution and DPI you choose) or SVG (vector).

The views paint themselves into the file (QWidget.render with a scaled painter), so text and
lines are redrawn at the export resolution rather than blown up from a screenshot. The views
that draw from a cached screen-sized pixmap (the map, the alignment view) check `ACTIVE` and
draw directly while an export runs, so they stay sharp (and vector in an SVG) too.

A figure is for paper: by default it is drawn in light colours even when the program is dark
(the "Light background" choice, remembered), and without the controls that only make sense on
screen: each view lists them in `export_hidden()` (zoom bars, scroll bars, the alignment view's
header row) and they are hidden while it paints.
"""
import re
from contextlib import contextmanager
from pathlib import Path

from PySide6.QtCore import QCoreApplication, QEvent, QPoint, QRect, QSettings, QSize, Qt
from PySide6.QtGui import QColor, QImage, QImageWriter, QPainter, QPalette, QRegion
from PySide6.QtSvg import QSvgGenerator
from PySide6.QtWidgets import (
    QApplication, QCheckBox, QComboBox, QDialog, QDialogButtonBox, QFormLayout, QLabel, QSpinBox, QWidget,
)

from primer_db import data_dir
import theme

ACTIVE = False   # True while an export paints: views with cached pixmaps draw directly
FORMATS = {"PNG": "png", "TIFF": "tiff", "SVG": "svg"}
SCREEN_DPI = 96  # what a "100 %" export assumes the screen is
INCH_CM = 2.54


def settings():
    return QSettings(str(data_dir() / "image_export.ini"), QSettings.IniFormat)


def light_wanted():
    """The "Light background" choice of the dialog, as last made (on until it is turned off)."""
    return settings().value("light", True, type=bool)


def light_palette():
    """Paper colours: white everywhere (the window colour too, so no grey shows around the views), black text, the faint
    text derived from it."""
    return QPalette(QColor(Qt.white))


@contextmanager
def exporting(root=None, light=False):
    """While an export paints: ACTIVE is set (the views draw directly), the controls that each view under `root` lists in
    export_hidden() are hidden (layouts follow, so the alignment view's tracks move up into its header's place), and with
    `light` the whole program has a light palette (the views read their colours from it, and theme.is_dark() turns False,
    so the theme's colours are the light ones). All of it is put back afterwards, before anything is drawn on screen."""
    global ACTIVE
    app = QApplication.instance()
    palette = QPalette(app.palette()) if light and app is not None else None
    hidden, owners = [], []
    ACTIVE = True
    try:
        if palette is not None:
            app.setPalette(light_palette())
            # (the widgets hear of a new palette through posted events: delivered now, before they paint)
            QCoreApplication.sendPostedEvents(None, QEvent.ApplicationPaletteChange)
        if root is not None:
            for w in [root] + root.findChildren(QWidget):
                parts = getattr(w, "export_hidden", None)
                if callable(parts) and w.isVisible():
                    owners.append(w)
                    for c in parts():
                        if c is not None and c.isVisible():
                            c.hide()
                            hidden.append(c)
            for w in owners:
                if w.layout() is not None:
                    w.layout().activate()
        yield
    finally:
        for c in hidden:
            c.show()
        for w in owners:
            if w.layout() is not None:
                w.layout().activate()
        if palette is not None:
            app.setPalette(palette)
            QCoreApplication.sendPostedEvents(None, QEvent.ApplicationPaletteChange)
        ACTIVE = False


def parts_of(root, views):
    """[(widget, offset in root)]: the root itself (pane headers and all), or only the given views at their place."""
    if views is None:
        return [(root, QPoint(0, 0))]
    return [(v, v.mapTo(root, QPoint(0, 0))) for v in views]


MAX_PIXELS = 100_000_000   # a raster of this many pixels takes 400 MB of memory; more would stop the program


def output_size(root, scale):
    return QSize(max(1, round(root.width() * scale)), max(1, round(root.height() * scale)))


def _paint(device, root, parts, scale, light):
    p = QPainter(device)
    p.setRenderHints(QPainter.Antialiasing | QPainter.TextAntialiasing | QPainter.SmoothPixmapTransform)
    p.scale(scale, scale)
    with exporting(root, light):
        p.fillRect(QRect(0, 0, root.width(), root.height()), root.palette().color(QPalette.Base))
        for w, offset in parts:
            w.render(p, offset, QRegion(), QWidget.DrawWindowBackground | QWidget.DrawChildren)
    p.end()


def _fix_svg(path, size, dpi):
    """Two touches to what QSvgGenerator writes:
    - width / height (in mm) become the size it prints at: its pixels at `dpi`;
    - text is filled "nonzero". Qt puts fill-rule="evenodd" on every group (right for its paths) and viewers apply it to
      text too: where the outlines of a glyph overlap themselves (variable fonts such as Cascadia Mono are drawn that
      way) even-odd filling leaves a hole, and the letters look blotched with the background colour."""
    mm = lambda px: f"{px / dpi * INCH_CM * 10:.3f}mm"
    path = Path(path)
    text = path.read_text(encoding="utf-8")
    text = re.sub(r'(<svg\b[^>]*?)width="[^"]*"', rf'\1width="{mm(size.width())}"', text, count=1)
    text = re.sub(r'(<svg\b[^>]*?)height="[^"]*"', rf'\1height="{mm(size.height())}"', text, count=1)
    text = text.replace("<text ", '<text fill-rule="nonzero" ')
    path.write_text(text, encoding="utf-8")


def write(root, path, fmt, scale=1.0, dpi=SCREEN_DPI, views=None, light=None):
    """Write `root` (a tab) to `path`. `fmt`: png / tiff / svg; `scale`: output pixels per screen pixel; `dpi`: written
    in the file, so it prints at size = pixels / dpi; `views`: only these widgets (the views without their pane headers),
    or None for the whole tab; `light`: in light colours whatever the theme (None: the dialog's last choice). Raises OSError
    when the file cannot be written."""
    light = light_wanted() if light is None else light
    parts = parts_of(root, views)
    size = output_size(root, scale)
    if fmt == "svg":
        gen = QSvgGenerator()
        gen.setFileName(str(path))
        gen.setSize(size)
        gen.setViewBox(QRect(QPoint(0, 0), size))
        # Font point sizes are turned into SVG units with the generator's resolution: it must be the screen's, as for
        # the views' own drawing (the painter's scale already enlarges everything). At the export DPI the text came out
        # dpi / 96 times too big. The printed size asked for is written into the width / height afterwards.
        gen.setResolution(root.logicalDpiY())
        gen.setTitle(root.window().windowTitle() or "openDNA")
        gen.setDescription("Exported from openDNA")
        _paint(gen, root, parts, scale, light)
        _fix_svg(path, size, dpi)
        return
    img = QImage(size, QImage.Format_ARGB32)
    img.fill(light_palette().color(QPalette.Base) if light else root.palette().color(QPalette.Base))
    dots_per_meter = round(dpi / INCH_CM * 100)
    img.setDotsPerMeterX(dots_per_meter)
    img.setDotsPerMeterY(dots_per_meter)
    _paint(img, root, parts, scale, light)
    writer = QImageWriter(str(path), fmt.encode())
    if fmt == "tiff":
        writer.setCompression(1)  # LZW: lossless, and much smaller than raw
    if not writer.write(img):
        raise OSError(writer.errorString())


class ImageExportDialog(QDialog):
    """Format, DPI and size of the image; the size in pixels follows the DPI unless you change it (100 % = the size on
    screen printed at screen resolution, so a 300 DPI export of a pane prints as big as it looks)."""

    def __init__(self, parent, root):
        super().__init__(parent)
        self.root = root
        self.setWindowTitle("Export image of the current tab")
        self.format = QComboBox()
        self.format.addItems(list(FORMATS))
        self.dpi = QSpinBox(minimum=72, maximum=2400, value=300, singleStep=50, suffix=" DPI")
        self.scale = QSpinBox(minimum=25, maximum=1000, value=round(300 / SCREEN_DPI * 100), singleStep=25, suffix=" %",
                              toolTip="Output pixels per screen pixel")
        self.headers = QCheckBox("Include the pane headers (view pickers, buttons)")
        self.light = QCheckBox("Light background (white paper, dark text), also when the program is dark",
                               checked=light_wanted())
        # Two lines kept for the warning even while there is none: the dialog no longer grows (and its buttons jump) when
        # it appears.
        self.warning = QLabel(styleSheet=f"color: {theme.hex('status_bad')}", wordWrap=True)
        # a fixed height (three lines, room for the warning wrapped in a wider font): it appears and goes without the dialog
        # changing size (with a minimum only, a font that wraps it on more lines made the dialog grow)
        self.warning.setFixedHeight(3 * self.warning.fontMetrics().lineSpacing() + 4)
        self.info = QLabel()
        self.dpi.valueChanged.connect(self._dpi_changed)
        for w in (self.format, self.scale):
            (w.currentIndexChanged if w is self.format else w.valueChanged).connect(self._update)
        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        self.ok = buttons.button(QDialogButtonBox.Ok)
        self.ok.setText("Export…")
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        form = QFormLayout(self)
        form.addRow("Format", self.format)
        form.addRow("Resolution", self.dpi)
        form.addRow("Size", self.scale)
        form.addRow("", self.headers)
        form.addRow("", self.light)
        form.addRow("", self.info)
        form.addRow("", self.warning)
        form.addRow(buttons)
        self._update()

    def _dpi_changed(self, dpi):
        self.scale.setValue(round(dpi / SCREEN_DPI * 100))  # same printed size as on screen
        self._update()

    def fmt(self):
        return FORMATS[self.format.currentText()]

    def _update(self):
        size = output_size(self.root, self.scale.value() / 100)
        dpi = self.dpi.value()
        cm = f"{size.width() / dpi * INCH_CM:.1f} × {size.height() / dpi * INCH_CM:.1f} cm"
        kind = "vector (lines and text stay sharp at any size)" if self.fmt() == "svg" else "pixels"
        self.info.setText(f"{size.width():,} × {size.height():,} {'px' if self.fmt() != 'svg' else 'units'}, "
                          f"prints at {cm} at {dpi} DPI\n{kind}")
        too_big = self.fmt() != "svg" and size.width() * size.height() > MAX_PIXELS
        self.warning.setText(f"Too large: {size.width() * size.height() / 1e6:,.0f} million pixels (at most "
                             f"{MAX_PIXELS / 1e6:,.0f}). Lower the resolution or the size." if too_big else "")
        self.ok.setEnabled(not too_big)

    def accept(self):
        settings().setValue("light", self.light.isChecked())   # write() reads it (light=None), and the next export starts so
        super().accept()

    def options(self):
        """(format, scale, dpi, headers). The light background is remembered by accept(), see light_wanted()."""
        return self.fmt(), self.scale.value() / 100, self.dpi.value(), self.headers.isChecked()
