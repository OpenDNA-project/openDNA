"""Help → Keyboard Shortcuts (F1): one tab for the menus, read from the menu bar itself (so it
never goes stale), and one tab per panel for the keys and mouse gestures the panels handle
themselves. When a panel learns a new key or gesture, add it to PANELS here.
"""
from PySide6.QtGui import QKeySequence
from PySide6.QtWidgets import (
    QDialog, QDialogButtonBox, QHeaderView, QTabWidget, QTableWidget, QTableWidgetItem, QVBoxLayout,
)

PANELS = {
    "Sequence view": [
        ("Click / drag", "Place the cursor / select bases"),
        ("Shift+click", "Extend the selection"),
        ("Drag the blue [ or ]", "Move that end of the selection"),
        ("Arrows, Home, End", "Move the cursor (Shift: extend the selection)"),
        ("A C G T N", "Type bases (replace the selection)"),
        ("Backspace / Delete", "Delete the selection or one base"),
        ("Ctrl+C / Ctrl+X / Ctrl+V", "Copy / cut / paste (features come along)"),
        ("Ctrl+A", "Select the whole sequence"),
        ("Click a feature or primer", "Select it"),
        ("Shift+click a primer", "With a primer on the other strand selected: select the product"),
        ("Double-click a feature or primer", "Edit it"),
        ("Right-click a feature", "Edit / delete / show its translation"),
        ("Right-click selected bases", "Add primer F / R, add feature, add gene, copy, copy reverse complement"),
        ("Ctrl+wheel, Ctrl+= / Ctrl+-", "Bigger / smaller letters"),
        ("Ctrl+0", "Default letter size"),
        ("Scroll past the end", "Circular sequences continue at base 1 (and back)"),
    ],
    "Map": [
        ("Click a feature", "Select it"),
        ("Click / drag the backbone", "Place the cursor / select bases"),
        ("Double-click a feature", "Edit it"),
        ("Right-click a feature", "Edit / delete / show its translation"),
    ],
    "Alignment view": [
        ("Wheel", "Zoom at the pointer"),
        ("Drag", "Pan"),
        ("Shift+drag", "Select bases"),
        ("Double-click", "Zoom in"),
        ("Click the overview bar", "Jump there"),
        ("Click a feature / primer / read", "Select it (a read: its details in the header)"),
        ("Shift+click a primer", "With a primer on the other strand selected: select the product"),
        ("Right-click a feature", "Edit / delete / show its translation"),
        ("Click a Query hits bar", "Select the hit (its details in the header)"),
        ("Click a track name", "Hide / show the track"),
        ("Right-click a track", "Expanded / condensed / hidden, scale"),
        ("Drag the line under a track", "Resize it"),
        ("Locus box + Enter", "Go to 12,000-15,000, a position or a name"),
    ],
    "Features table": [
        ("Click a row", "Select the feature"),
        ("Delete", "Delete the selected feature(s)"),
        ("Click ▸ / ▾", "Show / hide the parts of a gene"),
        ("Right-click", "Delete, hide/show (all selected rows), group as gene, add to the Features Database…"),
    ],
    "History panel": [
        ("Select a change + button (or double-click)", "Go back to just before it (undoes it and all later ones)"),
        ("Select a greyed change + button", "Redo up to it"),
        ("Drag ▶ Current sequence up / down", "Go back to before the changes it passes / redo them"),
        ("Opened … + button", "Back to the file as opened"),
    ],
    "Tabs and panes": [
        ("+ after the last tab", "New tab"),
        ("Double-click a tab", "Rename it"),
        ("Right-click a tab", "Rename / duplicate / close"),
        ("Middle-click a tab", "Duplicate it"),
        ("Drag a tab", "Reorder tabs"),
        ("Pane header ◫ / ⬓ / ✕", "Split right / split below / close the pane"),
    ],
}


def menu_shortcuts(window):
    """(menu path, shortcut) for every menu action that has a shortcut."""
    rows = []

    def walk(menu, path):
        for action in menu.actions():
            if action.menu():
                walk(action.menu(), path + [action.text()])
            elif not action.shortcut().isEmpty():
                name = " → ".join(path + [action.text()]).replace("&", "")
                rows.append((action.shortcut().toString(QKeySequence.NativeText), name))

    for top in window.menuBar().actions():
        if top.menu():
            walk(top.menu(), [top.text()])
    return rows


class ShortcutsDialog(QDialog):
    def __init__(self, window):
        super().__init__(window)
        self.setWindowTitle("Keyboard Shortcuts")
        self.resize(640, 520)
        self.tabs = QTabWidget()
        self.tabs.addTab(self._table(menu_shortcuts(window)), "Menus (anywhere)")
        for name, rows in PANELS.items():
            self.tabs.addTab(self._table(rows), name)
        buttons = QDialogButtonBox(QDialogButtonBox.Close)
        buttons.rejected.connect(self.reject)
        layout = QVBoxLayout(self)
        layout.addWidget(self.tabs)
        layout.addWidget(buttons)

    @staticmethod
    def _table(rows):
        table = QTableWidget(len(rows), 2)
        table.setHorizontalHeaderLabels(["Key / gesture", "Does"])
        table.verticalHeader().hide()
        table.setEditTriggers(QTableWidget.NoEditTriggers)
        table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeToContents)
        table.horizontalHeader().setStretchLastSection(True)
        for i, (key, text) in enumerate(rows):
            table.setItem(i, 0, QTableWidgetItem(key))
            table.setItem(i, 1, QTableWidgetItem(text))
        return table
