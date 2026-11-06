"""The History: annotations and primer changes are steps in it, and the primer ones can be undone with the document's Undo; the
primer menu (right-click Edit / Delete) on a primer in the sequence view, the map and the primer panel."""
import itertools
import os
import random
import sys
import tempfile
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ["OPENDNA_DATA_DIR"] = tempfile.mkdtemp(prefix="opendna_test_")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "app"))

from PySide6.QtCore import QPoint, QPointF, Qt
from PySide6.QtGui import QMouseEvent
from PySide6.QtWidgets import QApplication, QMenu

from opendna_core import Document

import feature_table
import primer_db
import primer_panel
from assembly_panel import FragmentEditor
from history_panel import HistoryPanel, group
from primer_db import add_primer, add_primers, db, remove_primers, update_primer
from sequence_view import SequenceView

app = QApplication.instance() or QApplication([])
rng = random.Random(8)
SEQ = "".join(rng.choice("ACGT") for _ in range(600))
FWD, REV = SEQ[100:124], SEQ[300:324]


@pytest.fixture(autouse=True)
def private_primer_db(tmp_path, monkeypatch):
    monkeypatch.setattr(primer_db, "_db", primer_db.PrimerDB(tmp_path / "primers.tsv"))


def editor():
    return FragmentEditor(Document.new_sequence("t", SEQ, False), Path("t.gb"))


def texts(ed):
    return [text for _, _, text in group(ed.doc.history()[0])]


def names():
    return [r["name"] for r in db().rows]


_databases = itertools.count()


def fresh_db(monkeypatch, tmp_path):
    """Start the next section of a merged test with an empty primer database, as a test of its own would have."""
    monkeypatch.setattr(primer_db, "_db", primer_db.PrimerDB(tmp_path / f"primers_{next(_databases)}.tsv"))


class Recorder(QMenu):
    """A menu that records what it would show instead of waiting for a click."""
    shown = []

    def exec(self, *args):
        Recorder.shown.append(self)


# ---- annotations -----------------------------------------------------------------------------------------------------------

def test_history_steps(monkeypatch, tmp_path):
    """Annotations and primer changes are History steps: adding / editing / deleting an annotation; primer changes undo and redo
    in order with the edits; editing / deleting / importing primers undo (a primer changed since is left alone); the History
    panel lists them and goes back through them; the primer dialog records new and edited primers."""
    # --- annotations ---
    ed = editor()
    ed.add_features([("misc_feature", "11..40", [("label", "alpha")])])
    index = ed.features[0].index
    ed.set_feature(index, "misc_feature", "11..60", [("label", "alpha")])
    ed.set_feature(index, "misc_feature", "11..60", [("label", "beta")])
    ed.delete_feature(index)
    assert texts(ed) == ["Added feature alpha", "Changed feature alpha", "Added feature beta; Deleted feature alpha", "Deleted feature beta"]
    # --- primer changes undo and redo with the edits ---
    fresh_db(monkeypatch, tmp_path)
    ed = editor()
    ed.doc.save_genbank(tmp_path / "t.gb")   # saved: not modified
    assert not ed.doc.modified
    assert add_primer(ed, "pFwd", FWD) == "added" and names() == ["pFwd"]
    assert not ed.doc.modified, "a primer database change does not make the document modified"
    assert [s.label for s in ed.primer_sites] == ["pFwd"]
    ed.select(10, 10)
    ed.type_bases("A")   # a base edit between the two primer steps
    assert add_primer(ed, "pRev", REV) == "added"
    assert texts(ed) == ["Added primer pFwd", "Inserted 1 bp at 11: A", "Added primer pRev"]
    assert "pRev" in names() and ed.doc.modified
    # Undo: the primer goes, then the base, then the other primer.
    assert ed.undo() and names() == ["pFwd"] and ed.doc.modified
    assert ed.undo() and ed.doc.sequence() == SEQ and not ed.doc.modified and names() == ["pFwd"]
    assert ed.undo() and names() == [] and not ed.primer_sites, "the primer is no longer matched on the sequence"
    assert not ed.undo()
    # Redo brings them back, in order.
    assert ed.redo() and names() == ["pFwd"] and ed.redo() and ed.doc.modified and ed.redo() and names() == ["pFwd", "pRev"]
    assert {s.label for s in ed.primer_sites} == {"pFwd", "pRev"}
    # A new edit drops what could have been redone, primer steps included.
    assert ed.undo() and ed.undo()
    ed.select(0, 0)
    ed.type_bases("C")
    assert not ed.redo() and names() == ["pFwd"]
    # --- editing and deleting a primer ---
    fresh_db(monkeypatch, tmp_path)
    ed = editor()
    add_primer(ed, "p1", FWD, "first")
    add_primer(ed, "p2", REV)
    row_before = dict(db().find("p1"))
    assert update_primer(ed, "p1", "p1b", FWD, "renamed") and names() == ["p1b", "p2"]
    assert texts(ed)[-1] == "Renamed primer p1 to p1b"
    assert update_primer(ed, "p1b", "p1b", FWD, "edited") and texts(ed)[-1] == "Edited primer p1b"
    assert not update_primer(ed, "p1b", "p2", FWD, ""), "a name that belongs to another primer is refused"
    assert len(texts(ed)) == 4, "and is not a history step"
    remove_primers(ed, ["p1b", "p2"])
    assert names() == [] and texts(ed)[-1] == "Deleted primers p1b, p2"
    assert ed.undo() and names() == ["p1b", "p2"] and db().find("p1b")["note"] == "edited"
    assert ed.undo() and db().find("p1b")["note"] == "renamed"
    assert ed.undo() and names() == ["p1", "p2"] and db().find("p1") == row_before, "back to the original, date included"
    assert ed.redo() and ed.redo() and ed.redo() and names() == []
    # Several steps at once, as the History panel does.
    assert ed.undo_steps(3) == 3 and names() == ["p1", "p2"] and ed.undo_steps(-3) == 3 and names() == []
    # --- importing primers ---
    fresh_db(monkeypatch, tmp_path)
    ed = editor()
    report = add_primers(ed, [("a1", FWD), ("a2", REV), ("a3", SEQ[400:424])], "Imported")
    assert report["added"] == 3 and texts(ed) == ["Imported 3 primers: a1, a2, a3"]
    assert add_primers(ed, [("a1", FWD)])["exists"] == 1 and len(texts(ed)) == 1, "nothing added, nothing recorded"
    errors = []
    ed.error.connect(errors.append)
    update_primer(None, "a2", "a2", REV, "changed elsewhere")   # another document changes a2 meanwhile
    assert ed.undo() and names() == ["a2"] and db().find("a2")["note"] == "changed elsewhere", "a1 and a3 go, a2 stays as it is"
    assert errors and "a2" in errors[-1] and "changed since" in errors[-1]
    assert ed.redo() and names() == ["a1", "a2", "a3"] and "a2" in errors[-1], "redo puts a1 and a3 back in their places; a2 is there already"
    # Undoing a delete whose name was taken by another primer meanwhile does not overwrite it.
    remove_primers(ed, ["a1"])
    add_primer(None, "a1", SEQ[500:524])
    assert ed.undo() and db().find("a1")["sequence"] == SEQ[500:524] and "a1" in errors[-1] and "used now" in errors[-1]
    # --- the History panel ---
    fresh_db(monkeypatch, tmp_path)
    ed = editor()
    panel = HistoryPanel(ed)
    add_primer(ed, "pA", FWD)
    ed.select(0, 0)
    ed.type_bases("G")
    remove_primers(ed, ["pA"])
    items = [panel.list.item(i).text() for i in range(panel.list.count())]
    assert items[1:4] == ["Added primer pA", "Inserted 1 bp at 1: G", "Deleted primer pA"], items
    panel.list.setCurrentRow(1)   # back to before "Added primer pA": everything after it is undone too
    panel.go()
    assert names() == [] and ed.doc.sequence() == SEQ
    panel.list.setCurrentRow(panel.list.count() - 1)   # the last (greyed) row: redo it all
    panel.go()
    assert names() == [] and ed.doc.sequence()[0] == "G" and len(ed.doc.history()[0]) == 3
    # --- the primer dialog ---
    fresh_db(monkeypatch, tmp_path)
    ed = editor()
    dialog = primer_panel.PrimerDialog(None, ed)
    dialog.name.setText("viaDialog")
    dialog.sequence.setText(FWD)
    dialog.accept()
    assert names() == ["viaDialog"] and texts(ed) == ["Added primer viaDialog"]
    dialog = primer_panel.PrimerDialog(None, ed, "viaDialog")
    dialog.note.setText("a note")
    dialog.accept()
    assert texts(ed)[-1] == "Edited primer viaDialog" and db().find("viaDialog")["note"] == "a note"
    ed.undo()
    ed.undo()
    assert names() == []


# ---- primers in the history ----------------------------------------------------------------------------------------------------


# ---- the primer menu --------------------------------------------------------------------------------------------------------------

def test_primer_menu(monkeypatch, tmp_path):
    """The primer menu (Edit / Delete / Import / Hide): its actions on a database primer and a temporary one, and the right-click
    on a primer in the sequence view, the map (linear and circular) and the primer panel."""
    # --- primer actions ---
    fresh_db(monkeypatch, tmp_path)
    ed = editor()
    add_primer(ed, "menuP", FWD)
    asked = []
    ed.edit_primer_requested.connect(asked.append)
    errors = []
    ed.error.connect(errors.append)
    menu = QMenu()
    feature_table.primer_actions(menu, ed, "menuP")
    edit, delete, imported, hide = [a for a in menu.actions() if not a.isSeparator()]
    assert imported.text() == "Import to This Sequence" and not imported.isEnabled(), "made from this sequence: already in the file"
    assert edit.text().startswith("Edit") and delete.text().startswith("Delete") and edit.isEnabled() and delete.isEnabled()
    assert hide.text() == "Hide in This File" and hide.isEnabled()
    edit.trigger()
    assert asked == ["menuP"]
    delete.trigger()
    assert names() == [] and "Undo brings it back" in errors[-1] and texts(ed)[-1] == "Deleted primer menuP"
    ed.undo()
    assert names() == ["menuP"]
    other = QMenu()
    feature_table.primer_actions(other, ed, "temp1")   # not in the database (a temporary primer)
    acts = [a for a in other.actions() if not a.isSeparator()]
    assert not acts[0].isEnabled() and not acts[1].isEnabled(), "Edit and Delete need a database primer"
    assert acts[3].text() == "Hide in This File" and acts[3].isEnabled(), "any named primer can be hidden in this file"
    # --- sequence view ---
    fresh_db(monkeypatch, tmp_path)
    ed = editor()
    add_primer(ed, "viewP", FWD)
    view = SequenceView(ed)
    view.resize(900, 500)
    view.show()
    app.processEvents()
    view.viewport().repaint()   # lays the primers out
    site = ed.primer_sites[0]
    row, col = divmod(site.start + 3, view.bases_per_line())   # the row that holds the primer's 4th base
    y = view.offset(row) - view.verticalScrollBar().value() + view._primer_lane_top(row, site.strand, ed.primer_lanes[site.key]) + view.primer_h / 2
    x = view.left + col * view.cw
    calls = []
    with monkeypatch.context() as m:
        m.setattr(feature_table, "primer_menu", lambda parent, editor, key, pos: calls.append(key))
        view.mousePressEvent(QMouseEvent(QMouseEvent.Type.MouseButtonPress, QPointF(x, y), QPointF(view.mapToGlobal(QPoint(int(x), int(y)))),
                                         Qt.RightButton, Qt.RightButton, Qt.NoModifier))
    assert calls == [site.key]
    # The real menu: the section is the primer's name, with Edit and Delete.
    monkeypatch.setattr(feature_table, "QMenu", Recorder)
    Recorder.shown.clear()
    feature_table.primer_menu(view, ed, site.key, QPoint(0, 0))
    menu = Recorder.shown[-1]
    assert ed.current_primer == site.key and [a.text() for a in menu.actions() if not a.isSeparator()][0].startswith("Edit")
    assert [a.text() for a in menu.actions() if a.isSeparator() and a.text()] == ["viewP"]
    # --- map and primer panel ---
    fresh_db(monkeypatch, tmp_path)
    from map_view import MapView
    for circular in (False, True):
        ed = FragmentEditor(Document.new_sequence("t", SEQ, circular), Path("t.gb"))
        add_primer(ed, "mapP", FWD)
        site = ed.drawn_primers[0]
        view = MapView(ed)
        view.resize(500, 500)
        view.show()
        app.processEvents()
        mid = (site.spans[0][0] + site.spans[0][1]) / 2
        if circular:
            import math
            c, r, _ = view._circle()
            angle = math.radians(view._angle(mid))
            point = QPointF(c.x() + (r + 5) * math.cos(angle), c.y() - (r + 5) * math.sin(angle))
        else:
            x0, w, y, _ = view._linear()
            point = QPointF(x0 + mid / ed.length * w, y - 3)
        assert view._primer_at(point) == site.key
        assert view._primer_at(QPointF(point.x(), point.y() + 40)) is None
        view.mouseMoveEvent(QMouseEvent(QMouseEvent.Type.MouseMove, point, QPointF(point), Qt.NoButton, Qt.NoButton, Qt.NoModifier))
        assert view.cursor().shape() == Qt.PointingHandCursor, "a hand over a primer mark"
        view.mousePressEvent(QMouseEvent(QMouseEvent.Type.MouseButtonPress, point, QPointF(point), Qt.LeftButton, Qt.LeftButton, Qt.NoModifier))
        assert ed.current_primer == site.key, "a click selects the primer"
        db().remove(["mapP"])
    # The primer panel: right-click on a row.
    ed = editor()
    add_primer(ed, "rowP", FWD)
    panel = primer_panel.PrimerPanel(ed)
    panel.resize(600, 300)
    panel.show()
    app.processEvents()
    Recorder.shown.clear()
    monkeypatch.setattr(primer_panel, "QMenu", Recorder)
    item = panel.tree.topLevelItem(0)
    panel._menu(panel.tree.visualItemRect(item).center())
    menu = Recorder.shown[-1]
    assert [a.text() for a in menu.actions() if a.isSeparator() and a.text()] == ["rowP"]
    delete = [a for a in menu.actions() if a.text().startswith("Delete")][0]
    delete.trigger()
    assert names() == []


