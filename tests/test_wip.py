"""The orange "WIP" tag of the unfinished features (wip.py): shown in each place a feature appears while it is marked unfinished
in wip.FEATURES, and gone from all of them when it is switched off there. Headless (Qt offscreen)."""
import os
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("OPENDNA_DATA_DIR", tempfile.mkdtemp(prefix="opendna_test_"))
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "app"))

import pytest
from PySide6.QtWidgets import QApplication, QMenu, QStyle, QStyleOptionMenuItem

import theme
import wip

app = QApplication.instance() or QApplication([])
theme.install_style(app)
DEMO = ROOT / "demo" / "pFG051_Chlor.gb"


def _menu_actions(win):
    return [a for m in win.menuBar().findChildren(QMenu) for a in m.actions()]


def _tagged(win):
    return {a.text().replace("&", "") for a in _menu_actions(win) if a.property(wip.PROPERTY) and wip.is_wip(a.property(wip.PROPERTY))}


def _extra(menu, action):
    option = QStyleOptionMenuItem()
    option.text, option.font = action.text(), menu.font()
    return wip.menu_extra(option, menu)


def test_an_unknown_feature_is_an_error_not_a_missing_tag():
    with pytest.raises(KeyError):
        wip.is_wip("golden_gat")
    assert set(wip.VIEWS.values()) <= set(wip.FEATURES)


def test_each_unfinished_feature_is_tagged_where_it_shows():
    from main import MainWindow
    from structure_dialog import StructureDialog
    import runner
    win = MainWindow()
    win.show()
    win.load(str(DEMO))
    app.processEvents()
    # menus: the protein interface and mapping reads
    assert _tagged(win) >= {"Protein FASTA (3D Structure)…", "From NCBI (Protein)…", "Map Reads (FASTQ) to This Sequence…",
                            "Map Reads to This Sequence…", "Check Mapping Tools…", "Install Mapping Tools (conda)…"}
    assert "Assembly" not in _tagged(win) and "Save" not in _tagged(win)
    tools = next(m for m in win.menuBar().findChildren(QMenu) if m.title() == "T&ools")
    map_reads = next(a for a in tools.actions() if a.text() == "&Map Reads to This Sequence…")
    digest = next(a for a in tools.actions() if a.text().startswith("Restriction"))
    assert _extra(tools, map_reads) > 0 and _extra(tools, digest) == 0, "room for the tag on the tagged items only"
    assert "Work in progress" in map_reads.toolTip()
    size = app.style().sizeFromContents(QStyle.CT_MenuItem, _option(tools, map_reads), tools.sizeHint(), tools)
    assert size.width() > tools.sizeHint().width()
    tools.popup(win.pos())
    app.processEvents()
    tools.grab()   # drawn with the tag (no error)
    tools.hide()
    # Assembly: Golden Gate and Sidewinder, not Gibson
    win.show_assembly()
    app.processEvents()
    panel = win.assembly
    assert panel.method_wip.isHidden(), "Gibson is finished"
    for m, key in ((1, "golden_gate"), (2, "sidewinder")):
        panel.method.setCurrentIndex(m)
        assert not panel.method_wip.isHidden() and panel.method_wip.key == key
    panel.method.setCurrentIndex(0)
    assert panel.method_wip.isHidden()
    # the Protein view: beside the view's name in the pane header
    from workspace import Pane
    pane = Pane(win.workspace, "Protein")
    assert not pane.wip.isHidden() and pane.wip.toolTip()
    popup = pane.picker.open_popup()   # the picker's list: the tag after "Protein", the three columns intact
    tags = popup.findChildren(wip.Badge)
    assert [t.key for t in tags] == ["protein"] and tags[0].parentWidget().text() == "Protein"
    assert len(popup.layout().children()) == 0 and popup.layout().count() == 3
    popup.close()
    pane.set_view("Map")
    assert pane.wip.isHidden()
    # the 3D structure window and the Map Reads window
    for dialog, key in ((StructureDialog(win, "p", "MKTAYIAKQR"), "protein"), (runner.MapReadsDialog(win, "s"), "mapping")):
        badges = dialog.findChildren(wip.Badge)
        assert [b.key for b in badges] == [key] and not badges[0].isHidden()


def _option(menu, action):
    option = QStyleOptionMenuItem()
    option.initFrom(menu)
    option.text, option.font = action.text(), menu.font()
    option.menuItemType = QStyleOptionMenuItem.Normal
    return option


def test_a_finished_feature_loses_every_tag_with_one_switch(monkeypatch):
    """Setting a feature to False in wip.FEATURES (or SHOW_TAGS to False) is all it takes."""
    from main import MainWindow
    from structure_dialog import StructureDialog
    from workspace import Pane
    monkeypatch.setitem(wip.FEATURES, "protein", (wip.FEATURES["protein"][0], False))
    win = MainWindow()
    win.show()
    win.load(str(DEMO))
    app.processEvents()
    tagged = _tagged(win)
    assert "Protein FASTA (3D Structure)…" not in tagged and "Map Reads to This Sequence…" in tagged
    proteins = next(m for m in win.menuBar().findChildren(QMenu) if m.title() == "&Proteins")
    assert all(_extra(proteins, a) == 0 for a in proteins.actions())
    assert "Work in progress" not in proteins.actions()[0].toolTip()
    assert Pane(win.workspace, "Protein").wip.isHidden()
    assert all(b.isHidden() for b in StructureDialog(win, "p", "MKTAYIAKQR").findChildren(wip.Badge))
    monkeypatch.setattr(wip, "SHOW_TAGS", False)
    assert not any(wip.is_wip(k) for k in wip.FEATURES)
    win.show_assembly()
    app.processEvents()
    win.assembly.method.setCurrentIndex(1)
    assert win.assembly.method_wip.isHidden()
