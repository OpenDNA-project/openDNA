"""The cog and "?" next to each pane's view picker: every view has a guide and a cog entry, the cog opens the matching options
page (or is disabled with a tooltip), and the guide popup opens."""
import os
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ["GENOMEEDIT_DATA_DIR"] = tempfile.mkdtemp(prefix="genomeedit_test_")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "app"))

from PySide6.QtWidgets import QApplication

import tool_settings
import view_help
import workspace
from test_cas9 import close, dna, open_doc
from genomeedit_core import Document

app = QApplication.instance() or QApplication([])


def picker_views():
    return [n for _, names in workspace.grouped_views() for n in names if n is not None]


def test_every_picker_view_has_a_guide_and_an_options_entry():
    """Every view has a guide and an options entry, and every options target is a real settings tab."""
    # --- guide and options entry per view ---
    for name in picker_views():
        lines = view_help.HELP[name]
        assert 4 <= len(lines) <= 8, name
        assert sum(len(x.split()) for x in lines) <= 130, f"{name}: the guide is too long"
        assert name in view_help.OPTIONS, name
    assert set(view_help.HELP) == set(view_help.OPTIONS) == set(picker_views())
    # --- options targets exist ---
    for name, target in view_help.OPTIONS.items():
        if target and target[0] == "settings":
            assert target[1] in tool_settings.TOOLS, f"{name}: no settings tab {target[1]}"


def test_cog_and_help_buttons_in_every_pane(monkeypatch):
    """The cog opens the matching options page (or is disabled with a tooltip) in every pane, and the "?" opens the guide popup."""
    # --- cog in every view ---
    doc = Document.new_sequence('t', dna(2000), False)
    w = open_doc(doc)
    try:
        opened = []
        monkeypatch.setattr(tool_settings.SettingsDialog, "exec",
                            lambda self: opened.append(self.tabs.tabText(self.tabs.currentIndex())))
        monkeypatch.setattr(w, "primer_options", lambda: opened.append("Primer options"))
        pane = w.workspace.panes(w.workspace.currentWidget())[0]
        for name in workspace.VIEWS:
            pane.set_view(name)
            target = view_help.OPTIONS[name]
            assert pane.cog.isEnabled() == (target is not None), name
            assert pane.cog.toolTip() == view_help.options_tip(name)
            if target is None:
                assert pane.cog.toolTip() == "No options for this view"
                continue
            opened.clear()
            pane.cog.click()
            expected = tool_settings.LABELS.get(target[1], target[1]) if target[0] == "settings" else "Primer options"
            assert opened == [expected], name
    finally:
        close(w)

    # --- the guide popup opens ---
    shown = []
    monkeypatch.setattr(view_help.ViewHelpDialog, "exec", lambda self: shown.append(self))
    w = open_doc(Document.new_sequence('t', dna(2000), False))
    try:
        pane = w.workspace.panes(w.workspace.currentWidget())[0]
        pane.set_view("Digest")
        pane.help.click()
        assert len(shown) == 1 and "restriction enzymes" in shown[0].label.text()
        assert pane.help.toolTip() and pane.cog.focusPolicy() and pane.help.focusPolicy()
    finally:
        close(w)
