"""The view picker of a pane: a button that opens a popup with one column per category (Views, Annotations, Tools).
Each column has a dark heading band and, under it, the views on a lighter background, lined up with the heading. The colours come
from theme.py (roles "group_*"), so the popup follows the light / dark theme and the same roles can be used for any
other grouped list."""
from PySide6.QtCore import QPoint, Qt, Signal
from PySide6.QtWidgets import QFrame, QHBoxLayout, QLabel, QPushButton, QVBoxLayout, QWidget

import theme
import wip

INDENT = 16   # px: left padding of an item: the same as its heading's, so each column lines up (the columns already group them)


def popup_style():
    """The stylesheet of the popup, from the current theme's group colours."""
    c = theme.hex
    return (f"QFrame#viewPicker {{ background: {c('group_border')}; border: 1px solid {c('group_border')}; }}"
            f"QWidget#column {{ background: {c('group_item_bg')}; }}"
            f"QLabel#heading {{ background: {c('group_header_bg')}; color: {c('group_header_fg')}; font-weight: 700; padding: 7px 16px; }}"
            f"QPushButton#item {{ background: {c('group_item_bg')}; color: {c('group_item_fg')}; text-align: left; border: none;"
            f" padding: 5px 16px 5px {INDENT}px; }}"
            f"QPushButton#item:hover {{ background: {c('group_item_hover_bg')}; }}"
            f"QPushButton#item[current=\"true\"] {{ font-weight: 700; }}"
            f"QFrame#divider {{ background: {c('group_border')}; border: none; margin: 0 {INDENT}px; }}")


class ViewPicker(QPushButton):
    """Shows the current view's name; a click opens the popup. `groups()` returns [(category, [view names])] when it opens."""
    chosen = Signal(str)

    def __init__(self, groups):
        super().__init__()
        self.groups, self.name, self.popup = groups, "", None
        # As wide as the longest view name it can show: picking "Plasmid Checkup" after "Map" no longer widens the button and
        # pushes the options and help buttons beside it along.
        names = [n for _, group in groups() for n in group if n]
        widest = max((self.fontMetrics().horizontalAdvance(f"{n}  ▾") for n in names), default=0)
        self.setFixedWidth(max(140, widest + 32))
        self.clicked.connect(self.open_popup)

    def current(self):
        return self.name

    def set_current(self, name):
        self.name = name
        self.setText(f"{name}  ▾")

    def open_popup(self, anchor=None):
        """Open the popup under the button, or under `anchor` (any widget: lets another button use this picker), a column
        per category, and return it."""
        popup = QFrame(self, Qt.Popup)
        popup.setObjectName("viewPicker")
        popup.setStyleSheet(popup_style())
        row = QHBoxLayout(popup)
        row.setContentsMargins(1, 1, 1, 1)
        row.setSpacing(1)
        for category, names in self.groups():
            column = QWidget()
            column.setObjectName("column")
            column.setAttribute(Qt.WA_StyledBackground, True)
            box = QVBoxLayout(column)
            box.setContentsMargins(0, 0, 0, 6)
            box.setSpacing(0)
            heading = QLabel(category)
            heading.setObjectName("heading")
            box.addWidget(heading)
            for name in names:
                if name is None:   # a thin line between two sets of views
                    line = QFrame(frameShape=QFrame.HLine)
                    line.setObjectName("divider")
                    line.setFixedHeight(1)
                    box.addSpacing(4)
                    box.addWidget(line)
                    box.addSpacing(4)
                    continue
                item = QPushButton(name)
                item.setObjectName("item")
                item.setFlat(True)
                item.setCursor(Qt.PointingHandCursor)
                item.setProperty("current", name == self.name)
                item.clicked.connect(lambda _checked=False, name=name: self._choose(name))
                if name in wip.VIEWS and wip.is_wip(wip.VIEWS[name]):   # an unfinished view: its tag after its name
                    tag = wip.Badge(wip.VIEWS[name], mouse=False)
                    tagged = QHBoxLayout(item)
                    tagged.setContentsMargins(0, 0, 16, 0)
                    tagged.addStretch()
                    tagged.addWidget(tag)
                    item.setMinimumWidth(item.sizeHint().width() + tag.width() + 2 * wip.GAP)
                    item.setToolTip(wip.tip(wip.VIEWS[name]))
                box.addWidget(item)
            box.addStretch()
            row.addWidget(column)
        anchor = anchor or self
        popup.move(anchor.mapToGlobal(QPoint(0, anchor.height())))
        popup.show()
        self.popup = popup
        return popup

    def _choose(self, name):
        if self.popup is not None:
            self.popup.close()
        self.chosen.emit(name)
