"""Where a primer comes from decides how it is drawn: in this file (regular colour), only in the primer database or designed and
unsaved (the temporary colour). A primer in both is the file's. Import to This Sequence makes a database primer part of the file."""
import os
import random
import sys
import tempfile
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ["OPENDNA_DATA_DIR"] = tempfile.mkdtemp(prefix="opendna_test_")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "app"))

from PySide6.QtGui import QColor
from PySide6.QtWidgets import QApplication, QMenu

from opendna_core import Document

import feature_table
import primer_db
from display_panel import DisplayPanel, PrimerSwatch
from editor import Editor
from history_panel import group
from hover import describe_primer
from map_view import MapView
from primer_db import add_primer, db
from sequence_view import PRIMER_COLOR, TEMP_COLOR, SequenceView, primer_color

app = QApplication.instance() or QApplication([])
rng = random.Random(21)
SEQ = "".join(rng.choice("ACGT") for _ in range(600))
DBONLY, MINE, BOTH, TEMP = SEQ[50:74], SEQ[150:174], SEQ[250:274], SEQ[350:374]


@pytest.fixture(autouse=True)
def private_primer_db(tmp_path, monkeypatch):
    monkeypatch.setattr(primer_db, "_db", primer_db.PrimerDB(tmp_path / "primers.tsv"))


def make(tmp_path, name="o.gb"):
    """An editor whose file stores MINE and BOTH; the database has DBONLY and BOTH; TEMP is designed here, not saved."""
    doc = Document.new_sequence("o", SEQ, True)
    doc.add_features([("primer_bind", "151..174", [("label", "mine"), ("note", "sequence: " + MINE)]),
                      ("primer_bind", "251..274", [("label", "both"), ("note", "sequence: " + BOTH)])])
    db().add("dbonly", DBONLY)
    db().add("both", BOTH)
    ed = Editor(doc, tmp_path / name)
    ed.set_extra([], [("temp", TEMP)])
    return ed


def sites(ed):
    return {s.label: s for s in ed.primer_sites}



def test_a_primer_is_drawn_by_where_it_comes_from(tmp_path):
    """The four origins (database only, file, both, temporary) get their flags and colour; the views draw each origin in its
    colour (sequence view and map); the display filter swatches show the two styles."""
    ed = make(tmp_path)
    # --- the flags and colour of each origin ---
    s = sites(ed)
    assert set(s) == {"dbonly", "mine", "both", "temp"}
    assert (s["dbonly"].from_db, s["dbonly"].from_file, s["dbonly"].provisional) == (True, False, True)
    assert (s["mine"].from_db, s["mine"].from_file, s["mine"].provisional) == (False, True, False)
    assert (s["both"].from_db, s["both"].from_file, s["both"].provisional) == (True, True, False), "in both: the file's style wins"
    assert (s["temp"].temporary, s["temp"].provisional) == (True, True)
    assert primer_color(s["dbonly"]) == primer_color(s["temp"]) == TEMP_COLOR, "the database primer looks like a temporary one"
    assert primer_color(s["mine"]) == primer_color(s["both"]) == PRIMER_COLOR

    # --- the views draw each origin in its colour ---
    ed.set_display(primers=True)

    def pixels(widget, color):
        image = widget.grab().toImage()
        want = QColor(color)
        return sum(1 for y in range(image.height()) for x in range(0, image.width())
                   if (lambda c: abs(c.red() - want.red()) + abs(c.green() - want.green()) + abs(c.blue() - want.blue()) < 12)(
                       image.pixelColor(x, y)))
    view = SequenceView(ed)
    view.resize(1000, 700)
    view.show()
    app.processEvents()
    assert pixels(view, TEMP_COLOR) > 0 and pixels(view, PRIMER_COLOR) > 0
    ed.set_display(file_primers=False)    # only the database's and the designed one are left: no regular colour
    app.processEvents()
    assert pixels(view, PRIMER_COLOR) == 0 and pixels(view, TEMP_COLOR) > 0
    ed.set_display(file_primers=True, db_primers=False)   # the designed one is drawn with "this sequence"
    app.processEvents()
    assert {x.label for x in ed.drawn_primers} == {"mine", "both", "temp"}
    ed.set_display(db_primers=True)
    MapView(ed).grab()   # the map draws them with the same function (no error, same colours by primer_color)

    # --- the display filter swatches show the two styles ---
    panel = DisplayPanel(ed)
    assert isinstance(panel.file_swatch, PrimerSwatch) and isinstance(panel.db_swatch, PrimerSwatch)
    assert panel.file_swatch.color == QColor(PRIMER_COLOR) and panel.db_swatch.color == QColor(TEMP_COLOR)
    assert panel.file_primers.toolTip() and panel.db_primers.toolTip() and panel.db_swatch.toolTip() == panel.db_primers.toolTip()
    panel.grab()
    panel.db_primers.setChecked(False)
    assert {x.label for x in ed.drawn_primers} == {"mine", "both", "temp"}


def test_the_primers_tab_lists_the_sources_the_display_filter_shows(tmp_path):
    """Issue #7: the Primers tab follows the Display filter's two primer switches (From this sequence / From the database), as
    the views do, and says how many primers the filter keeps out of the list."""
    from primer_panel import PrimerPanel
    ed = make(tmp_path)
    panel = PrimerPanel(ed)

    def listed():
        return {panel.tree.topLevelItem(i).data(0, 0x0100) for i in range(panel.tree.topLevelItemCount())}   # (Qt.UserRole)
    assert listed() == {"dbonly", "mine", "both", "temp"}
    assert "not shown" not in panel.count.text()
    ed.set_display(db_primers=False)
    assert listed() == {"mine", "both", "temp"}, "a primer in the file and the database is the file's; a designed one too"
    assert "1 more not shown" in panel.count.text()
    ed.set_display(file_primers=False, db_primers=True)
    assert listed() == {"dbonly"}
    assert "3 more not shown" in panel.count.text()
    ed.set_display(primers=False)
    assert listed() == set()
    ed.set_display(primers=True)
    assert listed() == {"dbonly", "mine", "both", "temp"}
    ed.set_primer_hidden("mine", True)   # hidden by hand: still listed (greyed) while its source is shown
    assert "mine" in listed()


def test_import_to_this_sequence_menu_entries_and_one_undo_step(tmp_path):
    """The menu entries for each kind of primer; Import to This Sequence adds the primer to the file as one undo step."""
    ed = make(tmp_path)
    # --- menu entries for each kind of primer ---
    def entries(name):
        menu = QMenu()
        feature_table.primer_actions(menu, ed, name)
        return {a.text(): a.isEnabled() for a in menu.actions() if not a.isSeparator()}
    assert entries("dbonly")["Import to This Sequence"] and entries("dbonly")["Edit…"]
    assert entries("temp")["Import to This Sequence"] and not entries("temp")["Edit…"], "temporary: import, but no edit/delete"
    assert not entries("mine")["Import to This Sequence"] and not entries("both")["Import to This Sequence"]
    assert feature_table.primers_to_import(ed, ["dbonly", "temp", "both"]) == [("dbonly", DBONLY), ("temp", TEMP)]

    # --- import is one undo step ---
    steps = len(ed.doc.history()[0])
    menu = QMenu()
    feature_table.primer_actions(menu, ed, "dbonly")
    actions = {a.text(): a for a in menu.actions() if not a.isSeparator()}
    assert actions["Import to This Sequence"].isEnabled()
    assert actions["Edit…"].isEnabled() and actions["Delete from Primer Database"].isEnabled() and actions["Hide in This File"].isEnabled()
    actions["Import to This Sequence"].trigger()
    assert len(ed.doc.history()[0]) == steps + 1, "one undo step"
    assert [t for _, _, t in group(ed.doc.history()[0])][-1] == "Imported primer dbonly to this sequence"
    site = sites(ed)["dbonly"]
    assert site.from_file and not site.provisional and primer_color(site) == PRIMER_COLOR, "now a regular primer annotation"
    assert "in this file" in describe_primer(ed, site.key)
    again = QMenu()
    feature_table.primer_actions(again, ed, "dbonly")
    assert not next(a for a in again.actions() if a.text() == "Import to This Sequence").isEnabled()
    assert ed.undo()
    assert sites(ed)["dbonly"].provisional and "not in this file" in describe_primer(ed, sites(ed)["dbonly"].key)
    assert ed.redo() and not sites(ed)["dbonly"].provisional


def test_imported_and_made_primers_survive_save_and_reload_as_on_sequence_primers(tmp_path):
    """An imported primer (and a temporary one) is written and reloaded as a regular on-sequence primer; a primer made here
    is the file's too."""
    # --- import, save, reload ---
    ed = make(tmp_path)
    ed.import_primers([("dbonly", DBONLY), ("temp", TEMP)])
    assert not sites(ed)["temp"].provisional, "a temporary primer becomes the file's too"
    ed.save(tmp_path / "o.gb")
    stored = {f.label for f in Document(str(tmp_path / "o.gb")).features() if f.kind == "primer_bind"}
    assert stored == {"mine", "both", "dbonly", "temp"}, "a primer only in the database would not be written"
    db().remove(["dbonly"])            # not even in the database any more: the file still has it
    back = Editor(Document(str(tmp_path / "o.gb")), tmp_path / "o.gb")
    s = sites(back)
    assert all(s[n].from_file and not s[n].provisional for n in ("dbonly", "temp", "mine", "both")), "all regular after reload"
    assert not [x for x in sites(make(tmp_path, "p.gb")).values() if x.label == "temp" and not x.provisional], \
        "a fresh editor on the unsaved document still draws the designed primer as temporary"

    # --- a primer made here is the file's ---
    made = make(tmp_path, "m.gb")
    add_primer(made, "madehere", SEQ[450:474])
    assert not sites(made)["madehere"].provisional
    made.save(tmp_path / "m.gb")
    assert "madehere" in {f.label for f in Document(str(tmp_path / "m.gb")).features() if f.kind == "primer_bind"}
