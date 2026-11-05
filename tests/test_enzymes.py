"""Restriction enzymes: the cut finder on degenerate sites and on enzymes that cut on both sides, and the Digest tool."""
import os
import random
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ["OPENDNA_DATA_DIR"] = tempfile.mkdtemp(prefix="opendna_test_")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "app"))

from PySide6.QtCore import Qt
from PySide6.QtWidgets import QApplication

from opendna_core import Document

import enzymes
from assembly_panel import FragmentEditor
from enzymes import Enzyme
from tools import DigestDialog

app = QApplication.instance() or QApplication([])
rng = random.Random(11)
dna = lambda n: "".join(rng.choice("ACGT") for _ in range(n))

# Enzymes written the way a longer list would hold them: site (IUPAC) and the cuts as offsets of the gap from the site's start.
ACCI = Enzyme("AccI", "GTMKAC", ((2, 4),), False)               # GT^MKAC, M = A/C, K = G/T
BGLI = Enzyme("BglI", "GCCNNNNNGGC", ((7, 4),), False)          # GCCNNNN^NGGC: a 3' overhang
BCGI = Enzyme("BcgI", "CGANNNNNNTGC", ((-10, -12), (24, 22)), False)   # cuts on both sides of its site
BSAI = Enzyme("BsaI", "GGTCTC", ((7, 11),), False)
ECORI = Enzyme("EcoRI", "GAATTC", ((1, 5),), False)


def test_cut_finder_and_the_common_list():
    """Degenerate sites, two-sided cutters, both strands and the origin; the common list still cuts where it always did."""
    # --- AccI GT^MKAC: M is A or C, K is G or T.
    s = "AAA" + "GTAGAC" + "TTT" + "GTCTAC" + "AAA" + "GTGGAC"
    assert enzymes.find_cuts(s, False, [ACCI])["AccI"] == [(5, 7), (14, 16)], "GTAGAC and GTCTAC match; GTGGAC does not (M is A or C)"
    # --- BglI leaves a 3' overhang: the top cut is after the bottom one on the top strand's numbering.
    assert enzymes.find_cuts("TT" + "GCCAAGTTGGC" + "TT", False, [BGLI])["BglI"] == [(9, 6)]
    # --- BcgI cuts on both sides of CGANNNNNNTGC: two cut pairs per site, 2 nt 3' overhangs.
    s = dna(40) + "CGAGGGGGGTGC" + dna(40)
    assert enzymes.find_cuts(s, False, [BCGI])["BcgI"] == [(30, 28), (64, 62)]
    # --- a site that is not its own reverse complement is found on the other strand, with the cuts mirrored.
    fwd = enzymes.find_cuts("AA" + "GGTCTC" + "A" * 20, False, [BSAI])["BsaI"]
    rev = enzymes.find_cuts("A" * 20 + "GAGACC" + "AA", False, [BSAI])["BsaI"]
    assert fwd == [(9, 13)] and rev == [(15, 19)] and rev[0][1] - rev[0][0] == 4
    # --- across the origin of a circle.
    assert enzymes.find_cuts("ATTCAAAAAAGA", True, [ECORI])["EcoRI"] == [(11, 3)]
    # --- the common list
    common = {e.name: e for e in enzymes.COMMON}
    assert common["EcoRI"].cuts == ((1, 5),) and common["KpnI"].cuts == ((5, 1),) and common["EcoRV"].blunt
    assert common["BsaI"].cuts == ((7, 11),) and common["SapI"].cuts == ((8, 11),)
    assert enzymes.find_sites("TTGAATTCAA" + "GAGACCTTTTT" + "TCAA", False) == {"EcoRI": [3], "BsaI": [5]}
    assert enzymes.describe(common["EcoRI"]) == "GAATTC: G^AATTC, 5′ overhang 4 nt"
    assert enzymes.describe(common["KpnI"]).endswith("3′ overhang 4 nt") and "blunt" in enzymes.describe(common["EcoRV"])
    assert "more place" in enzymes.describe(BCGI)
    # A Type IIS enzyme cuts outside its site: written the way suppliers do, GGTCTC(1/5), so it is clear that the
    # overhang's sequence is whatever the user puts between the site and the cut (the point of Golden Gate).
    assert enzymes.describe(common["BsaI"]) == "GGTCTC: GGTCTC(1/5), 5′ overhang 4 nt"
    assert enzymes.describe(common["BbsI"]).startswith("GAAGAC: GAAGAC(2/6),")
    assert enzymes.describe(common["SapI"]) == "GCTCTTC: GCTCTTC(1/4), 5′ overhang 3 nt"


# ---- the Digest tool ---------------------------------------------------------------------------------------------------


def digest(seq, circular=False):
    doc = Document.new_sequence("t", seq, circular)
    return DigestDialog(None, FragmentEditor(doc, Path("t.gb")))


def listed(dlg):
    return dlg.enzymes.names()


def tick(dlg, name, state=Qt.Checked):
    dlg.enzymes.item_of(name).setCheckState(0, state)


def test_the_digest_dialog_filter_gel_ladders_and_lanes():
    """The Digest tool: the filter bar (name / site, then the criteria), ticks kept under a filter, the enzyme table's columns
    and Owned box; the gel's five ladders chosen by clicking it; lanes side by side, each with its own enzymes and table, named on
    the gel by their enzymes."""
    from enzymes import EnzymeFilter
    # --- the filter
    seq = dna(1500) + "GAATTC" + dna(1500)
    while seq.count("GAATTC") != 1:   # random bases can hold a second EcoRI site: keep the one we planted
        seq = dna(1500) + "GAATTC" + dna(1500)
    dlg = digest(seq)
    assert "EcoRI" in listed(dlg) and "common enzymes" in dlg.note.text()
    assert dlg.criteria.filter() == EnzymeFilter(1, None, 1), "every enzyme that cuts, to start with"
    assert set(listed(dlg)) == set(dlg.cuts) and listed(dlg) == sorted(listed(dlg), key=str.lower), "by name"
    # The columns: the tick, the sites here, the site recognised, the ends, owned.
    row = dlg.enzymes.item_of("EcoRI")
    assert [row.text(c) for c in range(1, 4)] == ["1", "GAATTC", "5′ 4 nt"] and row.checkState(4) == Qt.Unchecked
    assert dlg.enzymes.item_of("EcoRV") is None or dlg.enzymes.item_of("EcoRV").text(3) == "blunt"
    # The list as narrow as its columns, each column as wide as its widest cell or its title (not the style's own measure of a
    # header, which gives every column the room of a sort arrow: 100 px more on Windows).
    table = dlg.enzymes
    widths = [table.columnWidth(c) for c in range(5)]
    assert widths == [table.column_need(c) for c in range(5)]
    assert all(w <= max(table.sizeHintForColumn(c) + 4, table.header().fontMetrics().horizontalAdvance(table.headerItem().text(c)) + 40)
               for c, w in enumerate(widths))
    assert table.maximumWidth() == sum(widths) + table.verticalScrollBar().sizeHint().width() + 2 * table.frameWidth() + 4
    # Sorting by sites: numbers, not text.
    dlg.enzymes.sortByColumn(1, Qt.DescendingOrder)
    counts = [int(dlg.enzymes.topLevelItem(i).text(1)) for i in range(dlg.enzymes.topLevelItemCount())]
    assert counts == sorted(counts, reverse=True)
    dlg.enzymes.sortByColumn(0, Qt.AscendingOrder)
    # The filter looks at the name and at the site.
    dlg.filter.setText("ecor")
    assert "EcoRI" in listed(dlg) and all("ECOR" in n.upper() for n in listed(dlg))
    dlg.filter.setText("GAATTC")
    assert "EcoRI" in listed(dlg) and all(dlg.info[n].site == "GAATTC" for n in listed(dlg))
    # Unique cutters: one site only (the criteria of the bar); then 5' overhangs only.
    dlg.filter.setText("")
    dlg.criteria.max_sites.setValue(1)
    assert "EcoRI" in listed(dlg) and all(len(dlg.cuts[n]) == 1 for n in listed(dlg))
    dlg.criteria.max_sites.setValue(-1)   # back to any
    dlg.criteria.ends.setCurrentIndex(dlg.criteria.ends.findData("5"))
    assert "EcoRI" in listed(dlg) and all(dlg.enzymes.item_of(n).text(3).startswith("5′") for n in listed(dlg))
    dlg.criteria.ends.setCurrentIndex(0)
    # A tick stays when the filter hides the enzyme, and its cut stays on the sequence.
    tick(dlg, "EcoRI")
    dlg.filter.setText("zzzz")
    assert listed(dlg) == [] and dlg.chosen() == ["EcoRI"]
    assert [(m.top, m.bottom) for m in dlg.editor.cut_marks["digest"]] == [(1501, 1505)]
    dlg.filter.setText("")
    assert dlg.enzymes.item_of("EcoRI").checkState(0) == Qt.Checked
    # The tooltip says what the enzyme leaves.
    tip = dlg.enzymes.item_of("EcoRI").toolTip(0)
    assert "G^AATTC" in tip and "5′ overhang 4 nt" in tip
    # --- the gel has five ladders named in brackets, chosen by clicking it
    from PySide6.QtCore import QPointF
    from PySide6.QtGui import QMouseEvent
    from tools import LADDERS
    assert {n: (min(b), max(b)) for n, (b, _) in LADDERS.items()} == {
        "1 kb": (250, 10000), "1 kb Plus": (250, 25000), "100 bp": (100, 1500), "100 bp Plus": (100, 3000), "50 bp": (50, 500)}
    assert all(bright in bands for bands, bright in LADDERS.values())
    dlg = digest(dna(1500) + "GAATTC" + dna(1500))
    gel = dlg.gel
    assert gel.lanes[0][0] == "Ladder (1 kb)"
    gel.resize(300, 400)
    gel.show()
    assert gel.over_ladder(40) and not gel.over_ladder(gel.width() - 5)
    chosen = []
    gel.ladder_changed.connect(chosen.append)
    menu = gel.menu()
    assert [a.text().split("   ")[0] for a in menu.actions()] == list(LADDERS) and [a.isChecked() for a in menu.actions()][0]
    menu.actions()[2].trigger()   # "100 bp"
    assert gel.ladder == "100 bp" and chosen == ["100 bp"] and gel.lanes[0][0] == "Ladder (100 bp)" and gel.lanes[0][1] == LADDERS["100 bp"][0]
    for name in LADDERS:   # every ladder paints, with or without a digest
        gel.set_ladder(name)
        gel.grab()
    # The pointer tells the ladder can be clicked.
    gel.mouseMoveEvent(QMouseEvent(QMouseEvent.Type.MouseMove, QPointF(40, 100), QPointF(40, 100), Qt.NoButton, Qt.NoButton, Qt.NoModifier))
    assert gel.cursor().shape() == Qt.PointingHandCursor
    # --- lanes are side-by-side digests, each with its own enzymes and table
    seq = dna(1500) + "GAATTC" + dna(1500) + "GGATCC" + dna(500)
    while seq.count("GAATTC") != 1 or seq.count("GGATCC") != 1:
        seq = dna(1500) + "GAATTC" + dna(1500) + "GGATCC" + dna(500)
    dlg = digest(seq)
    assert len(dlg.lanes) == 1 and dlg.tabs.count() == 1 and not dlg.tabs.tabsClosable() and len(dlg.gel.lanes) == 2
    tick(dlg, "EcoRI")
    assert len(dlg.frags) == 2 and dlg.tabs.tabText(0) == "EcoRI"
    dlg.add_lane()
    assert dlg.tabs.count() == 2 and dlg.tabs.currentIndex() == 1 and dlg.tabs.tabsClosable()
    assert dlg.ticked == set() and len(dlg.frags) == 1, "a new lane starts uncut, with its own table"
    assert not dlg.editor.cut_marks, "the cuts drawn on the sequence are the selected digest's"
    tick(dlg, "BamHI")
    tick(dlg, "EcoRI")
    assert len(dlg.frags) == 3 and sum(f[0] for f in dlg.frags) == len(seq) and dlg.tabs.tabText(1) == "BamHI+EcoRI"
    assert [n for n, _ in dlg.gel.lanes] == ["Ladder (1 kb)", "EcoRI", "BamHI+EcoRI"]   # every digest named by its enzymes
    assert dlg.gel.lane_lines("BamHI+EcoRI") == ["BamHI", "EcoRI"], "one enzyme per line over the lane"
    assert dlg.gel.lane_lines("A+B+C+D+E") == ["A", "B", "+3 more"]
    assert dlg.table.columnWidth(0) < 120, "the fragment table's columns fit what they hold"
    assert not dlg.table.horizontalHeader().stretchLastSection() and dlg.table.columnWidth(3) < 200, "Ends too: not stretched"
    assert [len(b) for _, b in dlg.gel.lanes[1:]] == [2, 3], "one gel lane per digest, next to the ladder"
    dlg.gel.resize(400, 400)
    dlg.gel.grab()
    dlg.tabs.setCurrentIndex(0)   # back to the first digest: its ticks, its table, its cuts
    assert dlg.chosen() == ["EcoRI"] and len(dlg.frags) == 2 and [m.label for m in dlg.editor.cut_marks["digest"]] == ["EcoRI"]
    assert dlg.enzymes.item_of("EcoRI").checkState(0) == Qt.Checked
    assert dlg.enzymes.item_of("BamHI").checkState(0) == Qt.Unchecked
    dlg.tabs.tabCloseRequested.emit(0)   # ✕ on the first: the second remains, alone, and cannot be closed
    assert len(dlg.lanes) == 1 and dlg.chosen() == ["BamHI", "EcoRI"] and not dlg.tabs.tabsClosable() and len(dlg.gel.lanes) == 2
    dlg.remove_lane(0)
    assert len(dlg.lanes) == 1, "the last digest stays"
