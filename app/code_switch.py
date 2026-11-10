"""The genetic-code quick switch of the main toolbar: a button showing the table in use ("Code 11") that opens the list of tables.
It reads and writes the same setting as Tools > Settings > General > Genetic code, so the two never disagree; that page has the
checkbox that hides it."""
from PySide6.QtCore import Signal
from PySide6.QtGui import QAction, QActionGroup
from PySide6.QtWidgets import QMenu, QPushButton

import genetic_code
import tool_settings


class CodeSwitch(QPushButton):
    manage = Signal()   # "Manage codon tables…" was chosen: the owner opens the editor

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setToolTip("Genetic code used to translate (Tools > Settings > General > Genetic code)")
        self.menu = QMenu(self)
        self.setMenu(self.menu)
        self.actions_by_label = {}
        self.rebuild()
        tool_settings.notifier.changed.connect(self._changed)

    def rebuild(self):
        """The NCBI tables, a divider, then the recoded and your own codes, then the editor."""
        self.menu.clear()
        self.group = QActionGroup(self.menu)
        self.actions_by_label = {}
        divider = False
        for label in genetic_code.LABELS:
            if not divider and genetic_code.is_custom(int(label.split()[0])):
                self.menu.addSeparator()
                divider = True
            action = QAction(label, self.menu, checkable=True)
            action.triggered.connect(lambda _=False, label=label: self.choose(label))
            self.group.addAction(action)
            self.menu.addAction(action)
            self.actions_by_label[label] = action
        self.menu.addSeparator()
        self.menu.addAction("Manage codon tables…", self.manage.emit)
        self.update_label()

    def choose(self, label):
        """Use the table `label` (as the Settings page would: saved, and the views told)."""
        if label != tool_settings.value(genetic_code.TOOL, "table"):
            tool_settings.set_value(genetic_code.TOOL, "table", label)
            tool_settings.notifier.changed.emit(genetic_code.TOOL)

    def update_label(self):
        label = tool_settings.value(genetic_code.TOOL, "table")
        self.setText(f"Code {genetic_code.table_id()}")
        self.actions_by_label.get(label, self.actions_by_label[genetic_code.DEFAULT_LABEL]).setChecked(True)

    def _changed(self, tool):
        if tool == genetic_code.TOOL:
            if len(self.actions_by_label) != len(genetic_code.LABELS) or any(l not in self.actions_by_label for l in genetic_code.LABELS):
                self.rebuild()   # a code was added, renamed or deleted
            else:
                self.update_label()


def shown():
    return bool(tool_settings.value(genetic_code.TOOL, "show_switch"))
