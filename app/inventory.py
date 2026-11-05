"""The restriction enzymes the user has in the freezer (the Owned box of the Digest's enzyme list). Kept in the data
folder (enzyme_inventory.json), the same for every sequence; the display filter can show only these (enzymes.EnzymeFilter)."""
import json

from PySide6.QtCore import QObject, Signal

from primer_db import data_dir


class _Notifier(QObject):
    changed = Signal()


notifier = _Notifier()
_owned = None


def _file():
    return data_dir() / "enzyme_inventory.json"


def owned():
    """The names of the enzymes in the inventory."""
    global _owned
    if _owned is None:
        try:
            _owned = set(json.loads(_file().read_text(encoding="utf-8")))
        except (OSError, ValueError, TypeError):
            _owned = set()
    return _owned


def set_owned(name, on):
    """Put an enzyme in / take it out of the inventory."""
    if on == (name in owned()):
        return
    (owned().add if on else owned().discard)(name)
    try:
        _file().write_text(json.dumps(sorted(owned())), encoding="utf-8")
    except OSError:
        pass   # shortcut: not saved this time; the tag still holds until the app closes
    notifier.changed.emit()
