"""Work in progress: the orange "WIP" tag on the features that are not finished yet.

FEATURES BELOW IS THE SWITCH. Every place that shows a tag names its feature by key (wip.Badge("protein"),
wip.mark_action(action, "mapping"), ...) and shows the tag only while that feature is marked True there:
- a feature is finished: set it to False (or delete its line once its tags are removed from the code): every tag of it goes;
- a new unfinished feature: add a line, then tag the places it shows with the helpers here (an unknown key raises, so a typo
  cannot hide a tag without anyone noticing);
- SHOW_TAGS = False hides every tag at once (a release where nothing is labelled).

One look everywhere (a small rounded orange pill, white bold letters), drawn by `paint`:
- `Badge(key)`: a widget, beside a title, a view's name, a choice;
- `mark_action(action, key)`: a menu item (the pill at its right end, drawn by theme.PressedStyle through `menu_pill`);
- `Delegate(key_of)`: the rows of a list or of a combo box's popup.
"""
from PySide6.QtCore import QRectF, QSize, Qt
from PySide6.QtGui import QFont, QFontMetricsF, QPainter
from PySide6.QtWidgets import QStyle, QStyledItemDelegate, QToolTip, QWidget

import theme

SHOW_TAGS = True
# key: (what it is, said in the tag's tooltip; still a work in progress?)
FEATURES = {
    "golden_gate": ("Golden Gate assembly", True),
    "sidewinder": ("Sidewinder assembly", True),
    "protein": ("The protein interface (the Protein view and 3D structures)", True),
    "mapping": ("Mapping reads (bwa / minimap2, through WSL on Windows, bash on Linux and macOS)", True),
}
# Views of a pane's picker (workspace.VIEWS) that are a feature above: view name -> key.
VIEWS = {"Protein": "protein"}

TEXT = "WIP"
PAD_X, GAP = 5, 6   # inner room either side of the letters; the room between the tag and what it marks


def is_wip(key):
    """Is the feature `key` (a key of FEATURES) still tagged?"""
    if key not in FEATURES:
        raise KeyError(f"unknown WIP feature {key!r}: add it to wip.FEATURES")
    return SHOW_TAGS and FEATURES[key][1]


def tip(key):
    return f"Work in progress: {FEATURES[key][0]} is not finished yet. It works, but check its results; it may still change."


def _font(base):
    f = QFont(base)
    f.setBold(True)
    f.setPointSizeF(max(6.5, (base.pointSizeF() if base.pointSizeF() > 0 else 9) * 0.78))
    f.setLetterSpacing(QFont.AbsoluteSpacing, 0.6)
    return f


def size(base_font):
    """The tag's size next to text in `base_font`."""
    fm = QFontMetricsF(_font(base_font))
    return QSize(round(fm.horizontalAdvance(TEXT) + 2 * PAD_X), round(fm.height() + 2))


def paint(painter, rect, base_font):
    """The tag, centred in `rect` (QRect / QRectF) at its own size."""
    s = size(base_font)
    r = QRectF(rect)
    box = QRectF(r.left() + (r.width() - s.width()) / 2, r.top() + (r.height() - s.height()) / 2, s.width(), s.height())
    painter.save()
    painter.setRenderHint(QPainter.Antialiasing)
    painter.setPen(Qt.NoPen)
    painter.setBrush(theme.color("wip"))
    painter.drawRoundedRect(box, box.height() / 2, box.height() / 2)
    painter.setPen(theme.color("wip_text"))
    painter.setFont(_font(base_font))
    painter.drawText(box, Qt.AlignCenter, TEXT)
    painter.restore()


class Badge(QWidget):
    """The tag as a widget, shown only while `key` is WIP (hidden, and taking no room, once it is finished). `mouse=False`: clicks
    go through to what is under it (a tag drawn on a button)."""

    def __init__(self, key, parent=None, mouse=True):
        super().__init__(parent)
        self.key = None
        if not mouse:
            self.setAttribute(Qt.WA_TransparentForMouseEvents)
        self.set_key(key)

    def set_key(self, key):
        """Tag another feature (or none: None), e.g. the method chosen in a list."""
        self.key = key
        on = key is not None and is_wip(key)
        self.setToolTip(tip(key) if on else "")
        self.setFixedSize(size(self.font()))
        if not on:
            self.hide()
        elif self.parentWidget() is not None:   # (show() on a widget with no parent yet would open it as a window of its own)
            self.show()

    def paintEvent(self, event):
        paint(QPainter(self), self.rect(), self.font())


class Delegate(QStyledItemDelegate):
    """Rows of a list or of a combo box's popup, the WIP ones with the tag at their right end. `key_of(index)` gives a row's
    feature key, or None."""

    def __init__(self, key_of, parent=None):
        super().__init__(parent)
        self.key_of = key_of

    def _wip(self, index):
        key = self.key_of(index)
        return key is not None and is_wip(key)

    def sizeHint(self, option, index):
        s = super().sizeHint(option, index)
        if self._wip(index):
            s.setWidth(s.width() + size(option.font).width() + 2 * GAP)
        return s

    def paint(self, painter, option, index):
        super().paint(painter, option, index)
        if self._wip(index):
            s = size(option.font)
            r = option.rect
            paint(painter, QRectF(r.right() - GAP - s.width(), r.top(), s.width(), r.height()), option.font)

    def helpEvent(self, event, view, option, index):
        if self._wip(index):
            QToolTip.showText(event.globalPos(), tip(self.key_of(index)), view)
            return True
        return super().helpEvent(event, view, option, index)


# ---- menu items -------------------------------------------------------------------------------------------------------------
PROPERTY = "wip_key"


def mark_action(action, key):
    """A menu item of a WIP feature: the tag at its right end (while `key` is WIP), and its tooltip says so. Returns `action`."""
    action.setProperty(PROPERTY, key)
    if is_wip(key):
        action.setToolTip((action.toolTip() + "\n\n" if action.toolTip() != action.text().replace("&", "") else "") + tip(key))
    return action


def _menu_action(option, menu):
    """The WIP action a menu item's style option is drawn for, or None. (Found by its text: QMenu measures its items before
    they have a place, so the position cannot tell.)"""
    if menu is None or not hasattr(menu, "actions"):
        return None
    text = (option.text or "").split("\t")[0]
    for a in menu.actions():
        key = a.property(PROPERTY)
        if key and a.text() == text and is_wip(key):
            return a
    return None


def menu_extra(option, menu):
    """Width a menu item needs for its tag (0 for the others): theme.PressedStyle.sizeFromContents."""
    return size(option.font).width() + 2 * GAP if _menu_action(option, menu) is not None else 0


def menu_pill(option, painter, menu, style):
    """Draw the tag of a WIP menu item at its right end, left of a shortcut or a submenu arrow: theme.PressedStyle.drawControl."""
    if _menu_action(option, menu) is None:
        return
    s, r = size(option.font), option.rect
    right = r.right() - 2 * GAP
    parts = (option.text or "").split("\t")
    if len(parts) > 1 and parts[1]:
        right -= QFontMetricsF(option.font).horizontalAdvance(parts[1]) + 2 * GAP
    if option.menuItemType == option.MenuItemType.SubMenu:
        right -= style.pixelMetric(QStyle.PM_MenuHMargin) + 12
    paint(painter, QRectF(right - s.width(), r.top(), s.width(), r.height()), option.font)
