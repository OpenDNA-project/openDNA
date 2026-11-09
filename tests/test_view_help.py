"""The cog and "?" next to each pane's view picker: every view has a guide and a cog entry, the cog opens the matching options
page (or is disabled with a tooltip), and the guide popup opens."""
import os
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ["OPENDNA_DATA_DIR"] = tempfile.mkdtemp(prefix="opendna_test_")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "app"))

from PySide6.QtWidgets import QApplication

import tool_settings
import view_help
import workspace
from test_cas9 import close, dna, open_doc
from opendna_core import Document

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
    # --- every view has a page in Tools > Settings ---
    pages = {p for _, names in tool_settings.page_layout() for p in names}
    for name, page in view_help.OPTIONS.items():
        assert page in pages, f"{name}: no settings page {page}"


def test_cog_and_help_buttons_in_every_pane(monkeypatch):
    """The cog opens the matching options page (or is disabled with a tooltip) in every pane, and the "?" opens the guide popup."""
    # --- cog in every view ---
    doc = Document.new_sequence('t', dna(2000), False)
    w = open_doc(doc)
    try:
        opened = []
        monkeypatch.setattr(tool_settings.SettingsDialog, "exec",
                            lambda self: opened.append(self.tabs.currentPage()))
        pane = w.workspace.panes(w.workspace.currentWidget())[0]
        for name in workspace.VIEWS:
            pane.set_view(name)
            assert pane.cog.isEnabled() and pane.cog.toolTip() == view_help.options_tip(name)
            opened.clear()
            pane.cog.click()
            assert opened == [name], name   # the cog opens the page of its own view
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
