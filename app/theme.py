"""Light and dark themes (View > Theme): System follows the operating system, Light and Dark force one.
Qt (6.8 and later) restyles every widget, the title bar included, when the colour scheme is set, and the
views that draw themselves read their colours from the palette, so nothing else needs to know."""
from PySide6.QtCore import QRectF, QSettings, Qt
from PySide6.QtGui import QColor, QGuiApplication, QPainter, QPalette
from PySide6.QtWidgets import QProxyStyle, QStyle, QStyleFactory

from primer_db import data_dir

THEMES = {"system": ("System", Qt.ColorScheme.Unknown), "light": ("Light", Qt.ColorScheme.Light), "dark": ("Dark", Qt.ColorScheme.Dark)}


def _settings():
    return QSettings(str(data_dir() / "theme.ini"), QSettings.IniFormat)


def saved():
    """The theme chosen last time ("system" if none)."""
    name = _settings().value("theme", "system")
    return name if name in THEMES else "system"


def apply(name, remember=True):
    """Switch the whole application to a theme and remember it for the next start."""
    name = name if name in THEMES else "system"
    hints = QGuiApplication.styleHints()
    hints.setColorScheme(THEMES[name][1])
    global _linked
    if not _linked:   # the palette is made again whenever the scheme changes (a theme switch, or the system's own)
        hints.colorSchemeChanged.connect(lambda *_: apply_highlight())
        _linked = True
    apply_highlight()
    if remember:
        _settings().setValue("theme", name)
    return name


class PressedStyle(QProxyStyle):
    """Every pressed (checked) tool button or push button of the interface: the dark blue of the quick filters under white (the
    style's own pressed look is its accent, a light blue in the dark theme, which white icons cannot be read on).

    Drawn by the style itself, not by a style sheet: a style sheet rule for the checked state only gives a checked button another
    box (border, padding) than an unchecked one, so a button changed size when it was pressed and pushed its neighbours along.
    Here only the colours change; every size and margin stays the base style's."""

    def drawPrimitive(self, element, option, painter, widget=None):
        if element in (QStyle.PE_PanelButtonTool, QStyle.PE_PanelButtonCommand) and option.state & QStyle.State_On:
            painter.save()
            painter.setRenderHint(QPainter.Antialiasing)
            painter.setPen(Qt.NoPen)
            painter.setBrush(color("highlight"))
            painter.drawRoundedRect(QRectF(option.rect).adjusted(0.5, 0.5, -0.5, -0.5), 4, 4)
            painter.restore()
            return
        super().drawPrimitive(element, option, painter, widget)

    def drawControl(self, element, option, painter, widget=None):
        if element in (QStyle.CE_ToolButtonLabel, QStyle.CE_PushButtonLabel) and option.state & QStyle.State_On:
            option = type(option)(option)   # a copy: the label in white on the blue
            for role in (QPalette.ButtonText, QPalette.WindowText, QPalette.Text):
                option.palette.setColor(role, color("highlight_text"))
        super().drawControl(element, option, painter, widget)


def install_style(app):
    """The application's own style (Fusion, Windows 11...) wrapped in PressedStyle, once."""
    if not isinstance(app.style(), PressedStyle):
        app.setStyle(PressedStyle(QStyleFactory.create(app.style().name())))


_linked = False
_system_highlight = None


def selection_color():
    """What the sequence, map and genome views draw their selection, its brackets and the search marks with: the colour that
    was theirs before the interface took the blue highlight (the palette's own is for the bars, lists and text fields)."""
    return QColor(_system_highlight) if _system_highlight is not None else QColor(QGuiApplication.palette().color(QPalette.ColorRole.Highlight))


def apply_highlight():
    """Put the one highlight colour of ROLES into the application palette: what selects (rows of lists and tables, the selected
    text), and what the style draws in its accent (ticked boxes, the bar beside the current row), which would otherwise be the
    system's accent colour (a yellow with white text on it, for one). Menus and every widget read the palette: set here only."""
    QGuiApplication.setPalette(QPalette())   # back to the system's palette for this scheme: its highlight is kept for the views
    global _system_highlight
    _system_highlight = QColor(QGuiApplication.palette().color(QPalette.ColorRole.Highlight))
    p = QPalette(QGuiApplication.palette())
    for role, key in ((QPalette.ColorRole.Highlight, "highlight"), (QPalette.ColorRole.HighlightedText, "highlight_text"),
                      (QPalette.ColorRole.Accent, "accent")):
        p.setColor(role, color(key))
    QGuiApplication.setPalette(p)
    app = QGuiApplication.instance()
    if app is not None and hasattr(app, "setStyle"):
        install_style(app)
    refresh_views()


def refresh_views():
    """Every widget of every window draws again in the new colours, now: the views that keep a drawn image (the map, the
    alignment view) drop it through their theme_changed(), the others are just repainted (a palette change alone does not
    reach a widget whose own picture is cached, nor a colour that was put in a style sheet when the widget was made)."""
    app = QGuiApplication.instance()
    if app is None or not hasattr(app, "allWidgets"):
        return
    for w in app.allWidgets():
        changed = getattr(w, "theme_changed", None)
        if callable(changed):
            changed()
        w.update()


# Colours that belong to the theme, so every widget that needs them reads the same values (light, dark).
# "group_*": a list grouped under headings, like the view picker: a dark band for a heading, a lighter one for its items.
ROLES = {
    "group_header_bg": ("#cfd6df", "#1a1d21"),
    "group_header_fg": ("#1f2933", "#e8eaed"),
    "group_item_bg": ("#f3f5f8", "#2d3136"),
    "group_item_hover_bg": ("#dfe6ee", "#3b4148"),
    "group_item_fg": ("#1f2933", "#e8eaed"),
    "group_border": ("#aab4c0", "#4a5057"),
    # bases typed, pasted or replaced since the last save (sequence view): blue letters with a thin blue underline,
    # never a background (yellow is the search hit). Blue on white / on the dark base colour, far from the text colour.
    "edited_fg": ("#0b5fff", "#5aa9ff"),
    # a cut on the sequence (restriction digest, Cas9): a stroke through each cut strand ending in a small cross (scissors)
    "cut": ("#e53935", "#ff6f66"),
    # a Cas9 / Cas12a guide site: the protospacer box, the PAM box beside it and the bar that names the guide
    "cas9_proto": ("#ffa726", "#d98200"),
    "cas9_pam": ("#ffee58", "#b3a000"),
    "cas9_label": ("#42a5f5", "#3b8fd9"),
    # the arrow on a guide's cut (hollow, pointing up at the cut): red against the orange / yellow boxes in both themes
    "cas9_arrow": ("#d50000", "#ff1744"),
    # text on those boxes ("PAM", the guide's name): black in both themes, the box fills are light enough for it
    "cas9_text": ("#000000", "#000000"),
    # the restriction enzyme pill of the hover tooltip (see cut_pill.py): the rounded body, its edge, the recognition site's
    # band, the faint parts (pairing ticks, the notation); letters use pill_fg and the cut stroke the "cut" role
    "pill_bg": ("#f4f7fb", "#2b3138"),
    "pill_fg": ("#1f2933", "#e8eaed"),
    "pill_edge": ("#8a97a8", "#6b7480"),
    "pill_site": ("#cfe0f7", "#3a4f6e"),
    "pill_tick": ("#7a8696", "#8b95a1"),
    # THE highlight of the interface, set here and nowhere else (apply_highlight puts it in the palette): the dark blue of a
    # pressed quick filter of the tab row, for what is selected or switched on, with white on it (icons and text); and "accent",
    # what the style fills ticked boxes with (it draws a white tick in the light theme, a dark one in the dark theme, so the dark
    # theme's accent is a lighter blue of the same hue); pressed buttons do not use it: see PressedStyle
    "highlight": ("#3a4a5e", "#3a4a5e"),
    "highlight_text": ("#ffffff", "#ffffff"),
    "accent": ("#3a4a5e", "#8db4e8"),
    # under the mouse, on the quick filters
    "toggle_hover_bg": ("#e3e8ee", "#33383f"),
    # what a check says (primer checks, Spacer Finder verdicts, assembly cards, Plasmid Checkup, database status), as text:
    # each pair readable on its own theme's background (contrast 4.5 or more; the yellow warning was 2.0 on white, the green
    # and red under 3.2 on a dark background)
    "status_ok": ("#1e7b34", "#6fcf7f"),
    "status_warn": ("#a15c00", "#f2b84b"),
    "status_bad": ("#c62828", "#ff7a70"),
    "status_info": ("#5b6573", "#a3adb8"),
    "status_likely": ("#6b7a00", "#c5d34f"),   # Plasmid Checkup's "Likely"
    # a primer's exact binding site in the views (the deep purple was 2.0 on a dark background)
    "primer": ("#5e35b1", "#b39ddb"),
}

# The search hit band over the bases (with its transparency): the light theme's yellow; in the dark theme a deeper amber, so the
# light letters on it stay readable (the yellow made them 3.4 to 1)
HIT = {False: (255, 214, 0, 130), True: (138, 109, 0, 217)}


def hit_color():
    return QColor(*HIT[is_dark()])


class Live:
    """A table of theme colours read when used, so a theme switch is followed: {key: role} gives the role's "#rrggbb", and
    {key: (mark, role)} gives (mark, "#rrggbb"). Reads like the dict it replaces ([key], get, in, iteration over the keys)."""

    def __init__(self, spec):
        self.spec = dict(spec)

    def _value(self, v):
        return (v[0], hex(v[1])) if isinstance(v, tuple) else hex(v)

    def __getitem__(self, key):
        return self._value(self.spec[key])

    def get(self, key, default=None):
        return self._value(self.spec[key]) if key in self.spec else default

    def __contains__(self, key):
        return key in self.spec

    def __iter__(self):
        return iter(self.spec)

    def __len__(self):
        return len(self.spec)


def is_dark():
    """True when the application is drawn with a dark palette (whatever chose it: the theme or the system)."""
    return QGuiApplication.palette().color(QPalette.ColorRole.Window).lightness() < 128


def hex(role):
    """The colour of a role in the current theme, as "#rrggbb"."""
    light, dark = ROLES[role]
    return dark if is_dark() else light


def color(role):
    return QColor(hex(role))
