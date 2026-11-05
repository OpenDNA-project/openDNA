"""Buttons and menu entries that send something over the internet carry a small globe (icons.globe_icon) and say, in their tooltip,
what is sent and to whom. `mark` does both, so every such control looks and reads the same."""
from PySide6.QtWidgets import QAbstractButton

import icons

NOTE = "🌐 Sends data over the internet"


def mark(control, what):
    """Put the globe on a button (QAbstractButton) or menu entry (QAction) and add "Sends data over the internet: `what`" to its
    tooltip. Returns the control."""
    control.setIcon(icons.globe_icon())
    tip = control.toolTip()
    note = f"{NOTE}: {what}"
    control.setToolTip(f"{tip}\n\n{note}" if tip and tip != control_text(control) else note)
    control.setProperty("sends_online", True)   # how the tests find every marked control
    return control


def control_text(control):
    return control.text().replace("&", "")


def is_marked(control):
    return bool(control.property("sends_online"))


def marked(widget):
    """Every marked button under `widget` (for the tests)."""
    return [b for b in widget.findChildren(QAbstractButton) if is_marked(b)]
