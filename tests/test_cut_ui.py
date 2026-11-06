"""The cut-site interface: the restriction enzyme pill in the hover tooltip, the Cas9 guide look (red upward arrow, black text,
lane heights from font metrics), the Spacer Finder's automatic specificity check and right-click menu, and the Digest's
fragment click (selection plus the pinned end cuts)."""
import os
import random
import re
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ["OPENDNA_DATA_DIR"] = tempfile.mkdtemp(prefix="opendna_test_")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "app"))

from PySide6.QtCore import QPoint, QRectF, Qt
from PySide6.QtGui import QColor
from PySide6.QtWidgets import QApplication, QLabel

from opendna_core import Document

import cas9
import cut_pill
import cut_sites as cs
import enzymes
import hover
import theme
from assembly_panel import FragmentEditor
from reference import DocumentReference
from sequence_view import SequenceView
from test_cas9 import GUIDE, close, dna, open_doc, wait_for
from test_cut_sites import count, guide_doc, shown
from tools import DigestDialog

app = QApplication.instance() or QApplication([])
ENZ = {e.name: e for e in enzymes.COMMON}


def layout_of(name, seq, circular=False):
    top, bottom = enzymes.find_cuts(seq, circular, [ENZ[name]])[name][0]
    return cut_pill.pill_layout(ENZ[name], seq, circular, top, bottom), top, bottom


# ---- the pill ---------------------------------------------------------------------------------------------------------------


def test_pill_layouts():
    """The pill's pieces for each kind of end (5', 3', blunt, Type IIS), a site on the bottom strand, a site across the
    origin of a circle, and a cut far from its site marked with an ellipsis."""
    # --- pieces for each kind of end ---
    lay, *_ = layout_of("EcoRI", "TTGAATTCAA")
    assert lay.kind == "5′" and lay.pieces == ("G", "CTTAA", "AATTC", "G"), "the staggered pieces of the example"
    lay, *_ = layout_of("KpnI", "TTGGTACCAA")
    assert lay.kind == "3′" and lay.pieces == ("GGTAC", "C", "C", "CATGG")
    lay, *_ = layout_of("EcoRV", "TTGATATCAA")
    assert lay.kind == "blunt" and lay.pieces == ("GAT", "CTA", "ATC", "TAG")
    lay, *_ = layout_of("BsaI", "TTGGTCTCAAGGAATTTT")   # Type IIS: site then cut one / five bases on
    assert lay.kind == "5′" and lay.overhang == "AGGA" and lay.site[2] == 1

    # --- a site on the bottom strand and across the origin ---
    lay, *_ = layout_of("BsaI", "TTTTAAAAGGAAGAGACCTT")
    assert lay.site[2] == -1 and lay.overhang == "AGGA", "the cut is before a site that reads on the bottom strand"
    lay, *_ = layout_of("EcoRI", "ATTCAAAAAAGA", circular=True)   # GA at the end + ATTC at the start
    assert lay.pieces == ("G", "CTTAA", "AATTC", "G")

    # --- a cut far from its site is marked with an ellipsis ---
    far = enzymes.Enzyme("FarI", "GGTCTC", ((20, 24),), False)   # a made-up Type IIS enzyme, cuts 14 / 18 bases past its site
    seq = "AA" + "GGTCTC" + "ACGT" * 8
    lay = cut_pill.pill_layout(far, seq, False, 2 + 20, 2 + 24)
    inner = [c for c in lay.columns[1:-1] if not c.real]
    assert len(inner) == 1 and inner[0].top == "…", "the gap between the site and the cut is one marker"
    assert lay.columns[0].top == "…" and lay.columns[-1].top == "…"
    assert lay.kind == "5′" and len(lay.overhang) == 4


def test_pill_tooltips_and_pixmaps(monkeypatch):
    """The tooltip of a cut holds the pill and a line of words (Qt draws it), reaches unique cutters and the Digest's cut
    marks; the pill pixmaps are drawn in both themes and differ by kind of end."""
    # --- the tooltip of a cut holds the pill and a line of words ---
    seq = "TTGAATTCAA"
    top, bottom = enzymes.find_cuts(seq, False, [ENZ["EcoRI"]])["EcoRI"][0]
    html = cs.enzyme_mark_tip(FragmentEditor(Document.new_sequence("e", seq, False), Path("e.gb")), "EcoRI", top, bottom)
    assert re.search(r'<img src="data:image/png;base64,[A-Za-z0-9+/=]{200,}"', html)
    words = re.sub(r"<[^>]+>", "", html)
    assert "G^AATTC" in words and "5′ overhang 4 nt (AATT)" in words and "after base 3" in words
    # Qt draws it: a label with that rich text has the pill's body colour in it
    label = QLabel(html)
    label.adjustSize()
    img = label.grab().toImage()
    pill = theme.color("pill_bg").rgb()
    assert sum(1 for x in range(img.width()) for y in range(img.height()) if QColor(img.pixel(x, y)).rgb() == pill) > 500

    # --- the pill tooltips reach unique cutters and digest cut marks ---
    ed = FragmentEditor(Document.new_sequence("e", "TTGAATTCAA" + "T" * 30, False), Path("e.gb"))
    ed.set_display(enzymes=True)
    (pos, name), = [s for s in ed.shown_sites() if s[1] == "EcoRI"]
    assert "data:image/png" in hover.describe(ed, ("site", name, pos))
    dlg = DigestDialog(None, ed)
    dlg.ticked.add("EcoRI")
    dlg.reload()
    assert "data:image/png" in hover.describe(ed, ("cut", "digest", 0))

    # --- the pixmaps are drawn and differ by kind, in the light and the dark theme ---
    for dark in (False, True):
        monkeypatch.setattr(theme, "is_dark", lambda: dark)
        pms = {}
        for name, seq in (("EcoRI", "TTGAATTCAA"), ("KpnI", "TTGGTACCAA"), ("EcoRV", "TTGATATCAA"), ("BsaI", "TTGGTCTCAAGGAATTTT")):
            lay, *_ = layout_of(name, seq)
            pm = cut_pill.pill_pixmap(lay)
            assert not pm.isNull() and pm.width() > 100 and pm.height() > 50
            img = pm.toImage()
            pms[name] = img
            body = {QColor(img.pixel(x, y)).rgb() for x in range(img.width()) for y in range(img.height()) if QColor(img.pixel(x, y)).alpha() == 255}
            assert theme.color("pill_bg").rgb() in body, f"the rounded body is in the theme colour (dark={dark})"
            assert theme.color("cut").rgb() in body, f"the cut strokes are in the cut colour (dark={dark})"
        names = list(pms)
        for i, a in enumerate(names):
            for b in names[i + 1:]:
                assert pms[a] != pms[b], f"{a} and {b} must look different (dark={dark})"


# ---- the guide look ----------------------------------------------------------------------------------------------------------


def test_guide_look_and_interface_highlight(monkeypatch):
    """The guide arrow is red, hollow and points up; the guide text is black in both themes and fits its lane (a feature
    lane is one letter tall); the interface highlight is set in one place and follows the theme."""
    # --- the arrow is red, hollow and points up, in both themes ---
    for dark in (False, True):
        monkeypatch.setattr(theme, "is_dark", lambda: dark)
        doc, lead, n = guide_doc(1)
        view = shown(FragmentEditor(doc, Path("g.gb")), (1800, 300))
        red = theme.color("cas9_arrow")
        assert red.red() > 190 and red.green() < 60 and red.blue() < 80
        img = view.viewport().grab().toImage()
        lane = view._parts(0)[3] + 2
        h = view.lane_h
        cx = view.left + (lead + 20 - 3) * view.cw
        widths = []   # red pixels per row of the lane, near the cut
        for dy in range(int(h)):
            xs = [x for x in range(int(cx - view.cw * 2), int(cx + view.cw * 2)) if QColor(img.pixel(x, int(lane + dy))).rgb() == red.rgb()]
            widths.append(len(xs))
        assert sum(widths) > 10, f"the arrow is painted in its theme colour (dark={dark})"
        first = next(i for i, w in enumerate(widths) if w)
        widest = max(range(len(widths)), key=lambda i: widths[i])
        assert first <= 2 and widest <= int(h * 0.7) and widths[first] <= widths[widest], \
            f"the head (the wide part) is at the top: it points up (dark={dark})"
        assert widths[first] < widths[widest] or widest == first
        view.close()

    # --- the text is black in both themes and fits its lane ---
    assert theme.ROLES["cas9_text"] == ("#000000", "#000000")
    for dark in (False, True):
        monkeypatch.setattr(theme, "is_dark", lambda d=dark: d)
        doc, lead, n = guide_doc(1)
        ed = FragmentEditor(doc, Path("g.gb"))
        view = shown(ed, (1800, 300))
        assert theme.color("cas9_text") == QColor("black")
        # nothing in the guide's boxes is white: the text, where fonts exist here, is black
        img = view.viewport().grab().toImage()
        lane = int(view._parts(0)[3] + 2)
        x0, x1 = int(view.left + lead * view.cw), int(view.left + (lead + 23) * view.cw)
        whites = sum(1 for x in range(x0 + 2, x1 - 2) for y in range(lane + 2, lane + int(view.lane_h) - 2)
                     if QColor(img.pixel(x, y)).rgb() == QColor("white").rgb())
        assert whites == 0, f"dark={dark}"
        view.close()
    # a feature lane is one letter tall (it was never made taller for the guide text)
    ed = FragmentEditor(Document.new_sequence("a", "ACGT" * 50, False), Path("a.gb"))
    view = SequenceView(ed)
    for size in (6, 8, 11, 16, 24, 32):
        view.set_size(size)
        assert view.lane_h == view.lh
    view.close()

    # --- the interface highlight is set in one place and follows the theme ---
    monkeypatch.undo()   # the real is_dark again
    from PySide6.QtGui import QPalette
    for name in ("dark", "light"):
        theme.apply(name, remember=False)
        app.processEvents()
        pal = app.palette()
        for role, key in ((QPalette.Highlight, "highlight"), (QPalette.HighlightedText, "highlight_text"), (QPalette.Accent, "accent")):
            assert pal.color(role) == theme.color(key), (name, key)
    theme.apply("system", remember=False)


# ---- Spacer Finder: automatic check and right-click ---------------------------------------------------------------------------


def many_guides_doc():
    rng = random.Random(3)
    return Document.new_sequence("many", "".join(rng.choice("ACGT") for _ in range(2500)), False)


def test_spacer_finder_automatic_check(monkeypatch):
    """The first 100 guides are checked automatically against the reference (and sorting keeps working); without a
    reference nothing is checked and the header says so; the automatic check can be cancelled."""
    # --- the first 100 guides are checked against the reference ---
    doc = many_guides_doc()
    monkeypatch.setattr(cas9, "genome", lambda: DocumentReference(Document.new_sequence("ref", dna(300) + GUIDE + "TGG" + dna(300), False)))
    monkeypatch.setattr(cas9, "CHUNK", 7)
    win = open_doc(doc)
    tool = cas9.Cas9Dialog(win, win.editor)
    total = len(tool.guides)
    assert total > 100
    wait_for(lambda: tool.check.isEnabled() and not tool.busy.isVisible(), 20000)
    assert all(s is not None for s in tool.specs[:100]) and all(s is None for s in tool.specs[100:])
    assert f"100 of {total} checked" in tool.status.text()
    assert tool.table.item(0, 6).text() != "–" and tool.table.item(0, 7).text() != "–", "here and genome columns are filled"
    assert tool.table.item(100, 6).text() == "–"
    # sorting keeps working: by the Score column, a row still maps to its guide
    tool.table.sortItems(5, cas9.Qt.AscendingOrder)
    scores = [int(tool.table.item(r, 5).text()) for r in range(tool.table.rowCount())]
    assert scores == sorted(scores)
    r = 3
    i = tool._index(r)
    assert tool.table.item(r, 5).text() == str(tool.guides[i].score)
    close(win)

    # --- without a reference nothing is checked and the header says so ---
    monkeypatch.undo()
    monkeypatch.setattr(cas9, "genome", lambda: None)
    win = open_doc(many_guides_doc())
    tool = cas9.Cas9Dialog(win, win.editor)
    assert "no reference" in tool.table.horizontalHeaderItem(7).text().lower()
    assert "No background genome" in tool.table.horizontalHeaderItem(7).toolTip()
    assert not tool.busy.isVisible() and all(s is None for s in tool.specs)
    assert "not checked automatically" in tool.status.text()
    close(win)

    # --- the automatic check can be cancelled ---
    monkeypatch.undo()
    monkeypatch.setattr(cas9, "genome", lambda: DocumentReference(Document.new_sequence("ref", dna(300), False)))
    monkeypatch.setattr(cas9, "CHUNK", 1)
    real = cas9.scan

    def slow(*a, **k):
        import time
        time.sleep(0.05)
        return real(*a, **k)

    monkeypatch.setattr(cas9, "scan", slow)
    win = open_doc(many_guides_doc())
    tool = cas9.Cas9Dialog(win, win.editor)
    assert tool.busy.isVisibleTo(tool) and not tool.check.isEnabled()
    tool.stop_auto()
    wait_for(lambda: tool.check.isEnabled(), 10000)
    done = sum(s is not None for s in tool.specs)
    assert done < 100 and "(cancelled)" in tool.status.text()
    close(win)


def test_right_click_adds_the_same_feature_as_the_button(monkeypatch):
    monkeypatch.setattr(cas9, "genome", lambda: None)
    seq = dna(80) + GUIDE + "TGG" + dna(60) + GUIDE + "AGG" + dna(60)
    win = open_doc(Document.new_sequence("t", seq, False))
    ed = win.editor
    tool = cas9.Cas9Dialog(win, ed)
    tool.resize(900, 400)
    tool.show()
    app.processEvents()
    tool.table.setCurrentCell(0, 0)
    tool.add()
    via_button = ed.doc.features()[-1]
    # a fresh window: the same row through the menu
    win2 = open_doc(Document.new_sequence("t", seq, False))
    tool2 = cas9.Cas9Dialog(win2, win2.editor)
    tool2.resize(900, 400)
    tool2.show()
    app.processEvents()
    pos = tool2.table.visualItemRect(tool2.table.item(0, 2)).center()
    menu = tool2.build_menu(pos)
    (action,) = menu.actions()
    assert action.text() == "Add as feature"
    action.trigger()
    via_menu = win2.editor.doc.features()[-1]
    assert (via_menu.kind, via_menu.location, via_menu.qualifiers) == (via_button.kind, via_button.location, via_button.qualifiers)
    assert any(k == cs.CAS9_KEY for k, _ in via_menu.qualifiers)
    # several rows selected: one entry adding them all
    tool2.table.selectAll()
    menu = tool2.build_menu(tool2.table.visualItemRect(tool2.table.item(1, 2)).center())
    (action,) = menu.actions()
    assert "(2)" in action.text()
    before = len(win2.editor.doc.features())
    action.trigger()
    assert len(win2.editor.doc.features()) == before + 2
    assert tool2.build_menu(tool2.table.viewport().rect().bottomRight()) is None, "no menu below the last row"
    close(win)
    close(win2)


# ---- Digest: a fragment clicked ------------------------------------------------------------------------------------------------


def digest_of(seq, circular, names):
    ed = FragmentEditor(Document.new_sequence("d", seq, circular), Path("d.gb"))
    dlg = DigestDialog(None, ed)
    for n in names:
        dlg.ticked.add(n)
    dlg.reload()
    return ed, dlg


def test_clicking_a_digest_fragment():
    """A clicked linear fragment is selected and its two end cuts pinned (drawn and hoverable, cleared by another
    selection or a new listing); a fragment across the origin of a circle is one selection passing the origin."""
    # --- a linear fragment selects it and pins its two end cuts ---
    seq = "T" * 20 + "GAATTC" + "T" * 30 + "GGATCC" + "T" * 20
    ed, dlg = digest_of(seq, False, ["EcoRI", "BamHI"])
    rows = {tuple(dlg.table.item(r, c).text() for c in (1, 2)): r for r in range(dlg.table.rowCount())}
    mid = rows[("22", "57")]   # between the two cuts: EcoRI after base 21, BamHI after base 57
    dlg.table.setCurrentCell(mid, 0)
    assert (ed.start, ed.end) == (21, 57)
    pinned = ed.cut_marks["digest_ends"]
    assert sorted(m.label for m in pinned) == ["BamHI", "EcoRI"] and all(m.pinned for m in pinned)
    assert "data:image/png" in hover.describe(ed, ("cut", "digest_ends", 0)), "the same pill tooltip as a hover"
    # the sequence view draws them and its hover reaches them
    view = shown(ed, (1200, 300))
    assert any(r in view.enzyme_rows for r in (0,)), "their names take a line above the bases"
    _, seq_top, _, _ = view._parts(0)
    x = view.left + pinned[0].top * view.cw
    assert view.cut_at(x, view.offset(0) + seq_top + 2)[1] in ("digest_ends", "digest")
    colour = theme.color("cut")
    assert count(view, colour, (x - 3, view.offset(0) + seq_top, x + 3, view.offset(0) + seq_top + 8)) > 0
    view.close()
    # another selection clears them and the table row
    ed.select(0, 3)
    assert "digest_ends" not in ed.cut_marks and not dlg.table.selectedIndexes()
    # an end of the sequence has no enzyme: only the one cut is pinned
    first = rows[("1", "21")]
    dlg.table.setCurrentCell(first, 0)
    assert (ed.start, ed.end) == (0, 21) and [m.label for m in ed.cut_marks["digest_ends"]] == ["EcoRI"]
    dlg.tabs.setCurrentIndex(0)
    dlg.reload()
    assert "digest_ends" not in ed.cut_marks, "a new listing drops the pinned cuts"
    dlg.destroyed.emit()
    assert "digest" not in ed.cut_marks

    # --- a fragment across the origin of a circle ---
    seq = "T" * 10 + "GAATTC" + "T" * 40 + "GGATCC" + "T" * 10   # 72 bp; cuts after base 11 and after base 57
    ed, dlg = digest_of(seq, True, ["EcoRI", "BamHI"])
    rows = {tuple(dlg.table.item(r, c).text() for c in (1, 2)): r for r in range(dlg.table.rowCount())}
    wrap = rows[("58", "11")]   # from the BamHI cut, through the origin, to the EcoRI cut
    dlg.table.setCurrentCell(wrap, 0)
    assert ed.wraps and (ed.start, ed.end) == (57, 72 + 11), "one selection that passes the origin"
    assert ed.selected_spans() == [(57, 72), (0, 11)]
    assert sorted(m.label for m in ed.cut_marks["digest_ends"]) == ["BamHI", "EcoRI"]
    # the whole circle cut once: the fragment is the whole sequence, both ends the same cut
    ed2, dlg2 = digest_of(seq, True, ["EcoRI"])
    dlg2.table.setCurrentCell(0, 0)
    assert ed2.end - ed2.start == 72 and [m.label for m in ed2.cut_marks["digest_ends"]] == ["EcoRI", "EcoRI"]


# ---- the recognition site drawn on the sequence while a cut is hovered --------------------------------------------------------


def test_hovered_cut_site_enzyme_names_and_map_click(monkeypatch, tmp_path):
    """Hovering a cut draws its recognition site and puts the tip above; clicking an enzyme on the map selects its site
    (circle and line); enzyme names stay on the bases whether or not a primer is there."""
    from PySide6.QtTest import QTest
    from map_view import MapView
    # --- hovering a cut draws its site and puts the tip above ---
    seq = "TTTTTTTTTTGAATTCAAAAAAAAAA"
    ed = FragmentEditor(Document.new_sequence("t", seq, False), Path("t.gb"))
    view = shown(ed, (900, 200))
    top, bottom = enzymes.find_cuts(seq, False, [ENZ["EcoRI"]])["EcoRI"][0]
    ed.set_cut_marks("digest", [(top, bottom, "EcoRI", None)])
    assert cs.hover_site(ed, ("cut", "digest", 0)) == (10, 16, 11, 15), "site GAATTC, cuts G^AATTC"
    assert cs.hover_site(ed, ("cut", "digest", 5)) is None
    before = view.viewport().grab().toImage()
    view._hover_changed(("cut", "digest", 0))
    after = view.viewport().grab().toImage()
    assert view.hovered_site == (10, 16, 11, 15) and after != before, "the site is drawn on the bases"
    view._hover_changed(None)
    assert view.hovered_site is None
    ed.set_cut_marks("digest", [(1, 1, "Cas9", None)])
    assert cs.hover_site(ed, ("cut", "digest", 0)) is None, "a guide has no recognition site"
    view.close()

    # --- clicking an enzyme on the map selects its site ---
    rng = random.Random(7)
    seq = "".join(rng.choice("ACGT") for _ in range(300))
    seq = seq[:120] + "GAATTC" + seq[126:]
    for circular in (True, False):
        ed = FragmentEditor(Document.new_sequence("m", seq, circular), Path("m.gb"))
        ed.set_display(enzymes=True)
        view = MapView(ed)
        view.resize(500, 500)
        view.show()
        app.processEvents()
        view.grab()   # paints, which records where each enzyme is drawn
        hits = [(shape, key) for shape, key in view.site_hits if key[1] == "EcoRI"]
        assert hits, "EcoRI is a unique cutter of this sequence"
        shape, key = hits[0]
        pt = shape.center() if isinstance(shape, QRectF) else (shape[0] + shape[1]) / 2
        assert view._site_at(pt) == key
        revealed = []
        ed.reveal.connect(revealed.append)
        QTest.mouseClick(view, Qt.LeftButton, Qt.NoModifier, pt.toPoint())
        assert (ed.start, ed.end) == (120, 126) and revealed == [120], "the recognition site is selected and shown"
        view.close()

    # --- enzyme names stay on the bases whether or not a primer is there ---
    import primer_db
    monkeypatch.setattr(primer_db, "_db", primer_db.PrimerDB(tmp_path / "primers.tsv"))
    rng = random.Random(5)
    seq = "TTGAATTCAA" + "".join(rng.choice("ACGT") for _ in range(30)) + "GGATCC" + "".join(rng.choice("ACGT") for _ in range(60))
    ed = FragmentEditor(Document.new_sequence("e", seq, False), Path("e.gb"))
    ed.set_display(enzymes=True)
    view = shown(ed)
    assert view._name_y(0, 3, 5) == view._parts(0)[1] - view.lh == 0, "no primer: the name line is on top of the bases"
    primer_db.add_primer(ed, "pF", seq[0:20])   # one forward primer under the first name
    app.processEvents()
    _, seq_top, _, _ = view._parts(0)
    assert seq_top == view.primer_h, "the primer lane is next to the bases, with no empty line above it"
    for pos in (3, 60):
        assert view._name_y(0, pos, 5) == max(0, seq_top - view.lh), "the name stays on the bases, over the primer"
    above, below = view._tip_anchor(QPoint(int(view.left + 3 * view.cw), int(view.offset(0) + seq_top + 3)), ("site", "EcoRI", 3))
    assert above.x() == int(view.left + 3 * view.cw) and below.y() > above.y()
    view.close()


# ---- which restriction sites are shown; hiding features one by one ---------------------------------------------------------------


def filter_doc():
    rng = random.Random(11)
    # EcoRI once, AluI (AGCT, 4 bp) many times, BamHI twice
    body = "".join(rng.choice("CT") for _ in range(40))
    seq = body + "GAATTC" + body + "GGATCC" + body + "GGATCC" + body + "AGCT" + body + "AGCT" + body + "AGCT" + body
    return Document.new_sequence("f", seq, False)


def test_the_enzyme_filter():
    """The enzyme filter combines independent criteria (counts, size, palindromic, Type IIS, ends - all must hold); the Digest
    has a filter bar of its own."""
    from enzymes import EnzymeFilter, wanted
    # --- the filter combines independent criteria ---
    ed = FragmentEditor(filter_doc(), Path("f.gb"))
    ed.set_display(enzymes=True)
    names = lambda: {n for _, n in ed.shown_sites()}
    assert ed.enzyme_filter == EnzymeFilter() and "EcoRI" in names() and "BamHI" not in names() and "AluI" not in names()
    ed.set_display(enzyme_filter=EnzymeFilter(min_sites=1, max_sites=None, min_size=6))
    assert {"EcoRI", "BamHI"} <= names() and "AluI" not in names(), "sites of 6+ bp: every cut of those enzymes"
    assert sum(1 for _, n in ed.shown_sites() if n == "BamHI") == 2
    ed.set_display(enzyme_filter=EnzymeFilter(1, None, 1, None))
    assert {"EcoRI", "BamHI", "AluI"} <= names()
    assert [p for p, _ in ed.shown_sites()] == sorted(p for p, _ in ed.shown_sites())
    ed.set_display(enzyme_filter=EnzymeFilter(2, 3, 1, 4))   # cuts 2 to 3 times, sites of up to 4 bp: AluI (3), Sau3AI (2)
    assert names() == {"AluI", "Sau3AI"}

    # --- criteria: palindromic, Type IIS, ends ---
    e = ENZ
    anything = dict(min_sites=1, max_sites=None)
    assert wanted(EnzymeFilter(palindromic=True, **anything), e["EcoRI"], 1)
    assert not wanted(EnzymeFilter(palindromic=True, **anything), e["BsaI"], 1), "GGTCTC is not its own reverse complement"
    assert wanted(EnzymeFilter(type_iis=True, **anything), e["BsaI"], 1) and not wanted(EnzymeFilter(type_iis=True, **anything), e["EcoRI"], 1)
    assert wanted(EnzymeFilter(type_iis=False, **anything), e["EcoRI"], 1)
    assert wanted(EnzymeFilter(ends="5", **anything), e["EcoRI"], 1) and wanted(EnzymeFilter(ends="3", **anything), e["KpnI"], 1)
    assert wanted(EnzymeFilter(ends="blunt", **anything), e["EcoRV"], 1) and not wanted(EnzymeFilter(ends="blunt", **anything), e["EcoRI"], 1)
    assert wanted(EnzymeFilter(min_size=8, **anything), e["NotI"], 1) and not wanted(EnzymeFilter(min_size=8, **anything), e["EcoRI"], 1)
    both = EnzymeFilter(palindromic=True, ends="5", min_size=6, **anything)
    assert wanted(both, e["EcoRI"], 1) and not wanted(both, e["KpnI"], 1) and not wanted(both, e["EcoRV"], 1), "all criteria must hold"

    # --- the Digest's filter bar is its own: it lists, it does not change what the views draw ---
    ed = FragmentEditor(filter_doc(), Path("f.gb"))
    dlg = DigestDialog(None, ed)
    assert dlg.enzymes.names() == ["AluI", "BamHI", "EcoRI", "Sau3AI"], "every enzyme that cuts, by name"
    dlg.criteria.max_sites.setValue(1)
    assert dlg.enzymes.names() == ["EcoRI"] and ed.enzyme_filter == EnzymeFilter(), "the display filter is untouched"
    # Sites can be 0: 0 to 0 lists the enzymes that do not cut, 0 to any every common enzyme.
    from enzymes import COMMON
    dlg.criteria.min_sites.setValue(0)
    dlg.criteria.max_sites.setValue(0)
    assert dlg.criteria.filter().max_sites == 0 and "PstI" in dlg.enzymes.names() and "EcoRI" not in dlg.enzymes.names()
    assert dlg.enzymes.item_of("PstI").text(1) == "0"
    dlg.criteria.max_sites.setValue(-1)
    assert len(dlg.enzymes.names()) == len(COMMON)
    # The presets side by side, on one row.
    from PySide6.QtWidgets import QPushButton
    dlg.show()   # (laid out once shown)
    QApplication.processEvents()
    buttons = dlg.criteria.findChildren(QPushButton)
    tops = {b.mapTo(dlg, b.rect().topLeft()).y() for b in buttons}
    assert len({b.mapTo(dlg, b.rect().topLeft()).x() for b in buttons}) == 3, "side by side"
    assert len(dlg.criteria.findChildren(QPushButton)) == 3 and len(tops) == 1
    dlg.close()


def feature_key_of(f):
    from editor import feature_key
    return feature_key(f)


def test_display_panel_and_hiding_features(monkeypatch):
    """The display panel edits and mirrors the enzyme filter; the default filter, its presets and the inventory square;
    aligned squares, one primer master square and folded enzyme criteria; a feature can be hidden and shown again."""
    import inventory
    from enzymes import PRESETS, EnzymeFilter, wanted
    from display_panel import DisplayPanel
    from feature_table import FeatureModel, feature_actions
    from PySide6.QtWidgets import QMenu, QPushButton
    # --- the panel edits and mirrors the enzyme filter ---
    ed = FragmentEditor(filter_doc(), Path("f.gb"))
    panel = DisplayPanel(ed)
    panel.max_sites.setValue(3)
    panel.palindromic.setCurrentIndex(panel.palindromic.findData(True))
    panel.ends.setCurrentIndex(panel.ends.findData("5"))
    assert ed.enzyme_filter == EnzymeFilter(1, 3, 6, None, True, None, "5"), "each control changes only its own criterion"
    ed.set_display(enzyme_filter=EnzymeFilter(2, None, 6, 8, False, True, "3"))   # e.g. a loaded project
    assert (panel.min_sites.value(), panel.max_sites.value(), panel.max_size.value()) == (2, -1, 8), "no limit: the lowest value, shown as any"
    assert panel.max_sites.text() == "any"
    assert (panel.palindromic.currentData(), panel.type_iis.currentData(), panel.ends.currentData()) == (False, True, "3")
    panel.deleteLater()

    # --- default filter, presets and the inventory square ---
    monkeypatch.setattr(inventory, "_owned", set())
    monkeypatch.setattr(inventory, "_file", lambda: Path(tempfile.mkdtemp()) / "inv.json")
    assert EnzymeFilter() == EnzymeFilter(1, 1, 6), "the default is the unique cutters of 6+ bp"
    assert not wanted(EnzymeFilter(), ENZ["AluI"], 1) and wanted(EnzymeFilter(), ENZ["EcoRI"], 1)
    ed = FragmentEditor(filter_doc(), Path("f.gb"))
    ed.set_display(enzymes=True)
    panel = DisplayPanel(ed)
    for label, f, tip in PRESETS:
        button = next(b for b in panel.findChildren(QPushButton) if b.text() == label)
        button.click()
        assert ed.enzyme_filter == f, label
    # the inventory: a criterion of the filter, set by the square of the Digest's list
    ed.set_display(enzyme_filter=EnzymeFilter(1, None, 1, None, None, None, None, True))
    assert ed.shown_sites() == [], "nothing is in the inventory yet"
    dlg = DigestDialog(None, ed)
    dlg.enzymes.item_of("EcoRI").setCheckState(4, Qt.Checked)   # the Owned box of the Digest's list
    QApplication.processEvents()   # (the list is refilled after the click)
    assert "EcoRI" in inventory.owned() and {n for _, n in ed.shown_sites()} == {"EcoRI"}, "tagged: the filter shows it"
    assert dlg.enzymes.item_of("EcoRI").checkState(0) == Qt.Unchecked, "Owned is not the digest tick"
    dlg.criteria.inventory.setCurrentIndex(dlg.criteria.inventory.findData(True))   # the bar: only the owned enzymes
    assert dlg.enzymes.names() == ["EcoRI"]
    dlg.enzymes.item_of("EcoRI").setCheckState(4, Qt.Unchecked)
    QApplication.processEvents()
    assert inventory.owned() == set() and ed.shown_sites() == [] and dlg.enzymes.names() == []
    dlg.close()

    # --- aligned squares, one primer master square and folded enzyme criteria ---
    ed = FragmentEditor(filter_doc(), Path("f.gb"))
    ed.add_features([("CDS", "10..30", [("label", "c")]), ("promoter", "40..60", [("label", "p")])])
    panel = DisplayPanel(ed)
    panel.resize(320, 900)
    panel.show()
    app.processEvents()
    assert not panel.criteria_box.isVisible(), "the enzyme criteria start folded"
    xs = {box.mapTo(panel.widget(), box.rect().topLeft()).x() for box in
          [panel.all_primers, panel.file_primers, panel.db_primers, panel.enzymes, panel.cds, panel.frames, panel.units,
           *panel.kinds.values()]}
    assert len(xs) == 1, f"every square of the panel starts at the same x: {xs}"
    panel.expander.setChecked(True)
    assert panel.criteria_box.isVisible()
    panel.all_primers.setChecked(False)
    assert not ed.show_file_primers and not ed.show_db_primers and not panel.file_primers.isChecked(), "one square for all primers"
    panel.all_primers.setChecked(True)
    assert ed.show_file_primers and ed.show_db_primers
    panel.close()

    # --- a feature can be hidden and shown again ---
    ed = FragmentEditor(Document.new_sequence("h", "ACGT" * 100, False), Path("h.gb"))
    ed.add_features([("misc_feature", "10..40", [("label", "mine")])])
    f = ed.features[0]
    assert ed.visible(f) and ed.visible_features(0, 400)
    menu = QMenu()
    feature_actions(menu, None, ed, f.index)
    hide = next(a for a in menu.actions() if a.text() == "Hide")
    hide.trigger()
    assert not ed.visible(f) and not ed.visible_features(0, 400) and ed.shown == [], "hidden: no view draws it"
    assert f in ed.features, "it is still in the file"
    model = FeatureModel(ed)
    assert model.data(model.index(0, 0), Qt.ForegroundRole) is not None, "greyed in the feature table"
    menu = QMenu()
    feature_actions(menu, None, ed, f.index)
    next(a for a in menu.actions() if a.text() == "Show on Sequence").trigger()
    assert ed.visible(f) and ed.shown == [f]
    ed.set_features_hidden([feature_key_of(f)], True)
    ed.set_features_hidden(ed.hidden_features, False)
    assert ed.visible(f)
