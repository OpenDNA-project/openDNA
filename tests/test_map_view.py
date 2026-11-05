"""The plasmid map on a crowded record: it keeps the largest features that fit its size (the others are left out, a note says
how many), clicks and hovers find only what is drawn, and enlarging the map brings features back."""
import os
import random
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ["GENOMEEDIT_DATA_DIR"] = tempfile.mkdtemp(prefix="genomeedit_test_")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "app"))

from PySide6.QtCore import QPointF
from PySide6.QtWidgets import QApplication

from genomeedit_core import Document

from editor import Editor
from map_view import MIN_BAND, MapView

app = QApplication.instance() or QApplication([])


def crowded(circular, count=300, seed=1):
    rng = random.Random(seed)
    seq = "".join(rng.choice("ACGT") for _ in range(6000))
    ed = Editor(Document.new_sequence("t", seq, circular), Path("t.gb"))
    feats = [("CDS", "101..1600", [("label", "big")])]
    for i in range(count):
        s = rng.randrange(5800)
        feats.append(("misc_feature", f"{s + 1}..{s + rng.choice([6, 12, 20, 40, 150])}", [("label", f"site{i}")]))
    ed.add_features(feats)
    return ed


def size(f):
    return sum(e - s for s, e in f.spans)


def test_map_layout_crowding_clicks_scale_and_geometry():
    """A crowded map keeps its largest features (a bigger map holds more), an uncrowded one draws everything, clicks find only drawn
    features (a left-out feature's place is backbone), the position scale steps from an oligo to a genome, and the circle takes
    80% of the room while a linear map is centred."""
    from PySide6.QtGui import QFont, QFontMetricsF
    from map_view import tick_label, ticks
    import math
    # --- a crowded map keeps its largest features ---
    for circular in (True, False):
        ed = crowded(circular)
        view = MapView(ed)
        view.resize(500, 500)
        drawn, lanes, count, left_out = view._layout()
        assert left_out and len(drawn) + left_out == len(ed.map_features)
        assert any(f.label == "big" for f in drawn)
        assert min(map(size, drawn)) >= max(size(f) for f in ed.map_features if f not in drawn), "the largest are kept"
        band = view._circle()[2] if circular else view._linear()[3]
        assert band >= MIN_BAND and set(lanes) == {f.index for f in drawn}
        view.grab()   # paints, with the note
        # A bigger map holds more.
        view.resize(1100, 1100)
        assert view._layout()[3] < left_out
    # --- an uncrowded map draws everything ---
    ed = crowded(True, count=5)
    view = MapView(ed)
    view.resize(500, 500)
    drawn, lanes, count, left_out = view._layout()
    assert left_out == 0 and drawn is ed.map_features and lanes is ed.map_lanes
    # --- clicks find only drawn features: a point over a left-out feature's place is backbone, not that feature ---
    ed = crowded(True)
    view = MapView(ed)
    view.resize(500, 500)
    drawn, lanes, _, _ = view._layout()
    c, r, band = view._circle()
    kept = {f.index for f in drawn}
    hits = 0
    for f in ed.map_features[:80]:
        s, e = f.spans[0]
        a = math.radians(view._angle((s + e) / 2))
        for lane in range(4):
            outer, inner = view._lane_radii(lane, r, band)
            radius = (outer + inner) / 2
            hit = view._hit(QPointF(c.x() + radius * math.cos(a), c.y() - radius * math.sin(a)))[1]
            assert hit is None or hit in kept
            hits += hit is not None
    assert hits > 40, f"only {hits} of 320 clicks on drawn features hit one: the click test would pass with a _hit that finds nothing"
    # --- position scale steps from an oligo to a genome: ticks every 1, 2 or 5 x 10^k bp, far enough apart for their labels ---
    fm = QFontMetricsF(QFont())
    for length in (100, 999, 5_936, 48_502, 4_641_652, 10_000_000):
        for px in (300.0, 1400.0):
            step, minor = ticks(length, px, fm)
            assert str(step).strip("0") in ("1", "2", "5")
            assert 2 <= length // step <= 40, (length, px, step)
            assert step / length * px >= fm.horizontalAdvance(tick_label(length // step * step, step)) + 8
            assert minor is None or (step % minor == 0 and minor / length * px >= 4)
    assert all(ticks(n, 1400.0, fm)[0] <= ticks(n, 300.0, fm)[0] for n in (100, 10_000_000)), "a bigger map, finer ticks"
    assert [tick_label(p, s) for p, s in ((750, 250), (2500, 500), (15_000, 5_000), (2_500_000, 500_000), (1, 1_000_000))] == \
        ["750", "2,500", "15 kb", "2.5 Mb", "1"]
    # Painted on both shapes, at any size, without error.
    for circular in (True, False):
        ed = crowded(circular, count=10)
        view = MapView(ed)
        for px in (250, 900):
            view.resize(px, px)
            view.grab()

    # --- the circle takes 80% of the room and a linear map is centred ---
    view = MapView(crowded(True, count=5))
    for w, h in ((500, 500), (900, 600), (400, 800)):
        view.resize(w, h)
        assert abs(2 * view._radius() / min(w, h) - 0.8) < 0.01
    view = MapView(crowded(False, count=5))
    view.resize(800, 600)
    _, _, y, band = view._linear()
    assert 150 < y < 300, "the backbone sits in the middle of the pane, not at its top edge"
    view.grab()


# ---- linear maps, and a FASTA file of several sequences (one line per chromosome) ---------------------------------------
HUMAN = Path(os.environ.get("GENOMEEDIT_HUMAN_FASTA", "V:/HM_Softwares/TestData/human_genome/hg38.fa"))   # outside git; tests skip if absent


def multi_fasta_editor(lengths, current=0, folder=None):
    from genomeedit_core import Fasta
    folder = Path(folder or tempfile.mkdtemp(prefix="genomeedit_multi_"))
    rng = random.Random(3)
    path = folder / "multi.fa"
    path.write_text("".join(f">chr{i + 1}\n" + "".join(rng.choice("ACGT") for _ in range(n)) + "\n" for i, n in enumerate(lengths)))
    records = [(n, int(l)) for n, l in Fasta(str(path)).sequences()]
    ed = Editor(Document.from_fasta(str(path), records[current][0]), path)
    ed.fasta_record, ed.fasta_records = records[current][0], records
    return ed, records


def test_a_multi_fasta_lists_one_line_per_chromosome_and_a_click_opens_it():
    """A multi-FASTA lists one line per chromosome and a click opens it; a long list scrolls and shows the open one."""
    from PySide6.QtCore import QPointF, Qt
    from PySide6.QtGui import QMouseEvent
    from PySide6.QtCore import QEvent
    ed, records = multi_fasta_editor([300, 200, 100, 50], current=1)
    view = MapView(ed)
    view.resize(600, 500)
    assert view._rows() == 4 and view._top() > 0 and view._room() == 500 - view._top()
    assert [view._record_at(QPointF(10, 2 + i * 17 + 3)) for i in range(4)] == [n for n, _ in records]
    assert view._first_row() == 0
    got = []
    ed.record_requested.connect(got.append)
    press = QMouseEvent(QEvent.MouseButtonPress, QPointF(20, 2 + 2 * 17 + 3), Qt.LeftButton, Qt.LeftButton, Qt.NoModifier)
    view.mousePressEvent(press)
    app.processEvents()   # the request is sent after the click
    assert got == ["chr3"]
    view.grab()
    # One sequence: no list.
    ed1, _ = multi_fasta_editor([100])
    assert MapView(ed1)._top() == 0
    # --- a long chromosome list scrolls and shows the open one ---
    ed, records = multi_fasta_editor([60 + i for i in range(80)], current=70)
    view = MapView(ed)
    view.resize(600, 400)
    rows = view._rows()
    assert rows < 80 and view._top() <= 0.45 * 400
    first = view._first_row()
    assert first <= 70 < first + rows, "opens with the current sequence in view"
    view.scroll = 0
    assert view._record_at(QPointF(10, 5)) == records[0][0]
    view.scroll = 10_000
    assert view._first_row() == 80 - rows, "scrolling stops at the end"
    view.grab()


def test_human_genome_chromosomes_listed_drawn_and_zoomed():
    """The human chromosomes are listed and a linear chromosome is drawn, and a chromosome zooms to 50 bp."""
    import pytest
    if not HUMAN.exists():
        pytest.skip(f"no human genome at {HUMAN}")
    from genomeedit_core import Fasta
    records = [(n, int(l)) for n, l in Fasta(str(HUMAN)).sequences()]
    assert len(records) > 24 and dict(records)["chr1"] == 248_956_422
    chr21 = dict(records)["chr21"]
    ed = Editor(Document.from_fasta(str(HUMAN), "chr21"), HUMAN)
    ed.fasta_record, ed.fasta_records = "chr21", records
    assert ed.length == chr21 and not ed.doc.circular
    view = MapView(ed)
    view.resize(1000, 700)
    view.grab()
    assert view._first_row() <= [n for n, _ in records].index("chr21") < view._first_row() + view._rows()
    assert view._linear()[2] > view._top()
    # --- a zoomed human chromosome ---
    ed = Editor(Document.from_fasta(str(HUMAN), "chr21"), HUMAN)
    view = MapView(ed)
    view.resize(1000, 600)
    x0, w, y, _ = view._linear()
    for _ in range(80):
        wheel(view, x0 + w / 2, y)
    a, b = view._win()
    assert round(b - a) == 50
    view.grab()
    assert 0 <= a < b <= ed.length


# ---- zoom on a linear map ------------------------------------------------------------------------------------------------
def wheel(view, x, y, dy=120, dx=0, mods=0):
    from PySide6.QtCore import QPoint, QPointF, Qt
    from PySide6.QtGui import QWheelEvent
    pt = QPointF(x, y)
    view.wheelEvent(QWheelEvent(pt, view.mapToGlobal(pt), QPoint(), QPoint(dx, dy), Qt.NoButton, Qt.KeyboardModifier(mods),
                                Qt.NoScrollPhase, False))


def test_the_wheel_zooms_a_linear_map_at_the_pointer_and_shows_small_features():
    """The wheel zooms a linear map at the pointer and shows small features, the selection and hits are found at the zoomed scale,
    zooming out ends on the whole map; the pan slider and Shift+wheel pan, and the map's zoom is independent of the sequence view."""
    # --- wheel zoom ---
    ed = crowded(False, count=300)
    view = MapView(ed)
    view.resize(700, 500)
    assert view.view is None and view._layout()[3] > 0
    x0, w, y, _ = view._linear()
    x = x0 + w * 0.25
    under = view._pos_at(x)
    wheel(view, x, y)
    assert view.view is not None
    a, b = view._win()
    assert b - a < ed.length and abs(view._pos_at(x) - under) < 1, "the base under the pointer stays put"
    for _ in range(30):
        wheel(view, x, y)
    a, b = view._win()
    assert round(b - a) == 50, "deepest zoom"
    assert view._layout()[3] == 0, "few features in view: none left out"
    assert all(any(e > a and s < b for s, e in f.spans) for f in view._layout()[0])
    view.grab()
    # The selection and hits are found at the zoomed scale.
    wheel(view, x, y, dy=-120)
    assert view._hit(QPointF(x, y))[0] == int(view._pos_at(x))
    for _ in range(60):
        wheel(view, x, y, dy=-120)
    assert view.view is None and view._win() == (0.0, float(ed.length)), "zooming out ends on the whole map"
    # --- pan, slider and independence from the sequence view ---
    ed = crowded(False, count=10)
    view = MapView(ed)
    view.resize(700, 500)
    x0, w, y, _ = view._linear()
    wheel(view, x0 + 10, y, dy=120)
    wheel(view, x0 + 10, y, dy=120)
    a, b = view._win()
    wheel(view, x0 + 100, y, dy=120, mods=0x02000000)   # Shift+wheel pans (up = left)
    a2, b2 = view._win()
    assert abs((b2 - a2) - (b - a)) < 1e-6 and a2 != a
    view.zoom_slider.setValue(100)
    a, b = view._win()
    assert round(b - a) == 50
    view.zoom_slider.setValue(0)
    assert view.view is None
    # Another view going to a base outside the zoomed map brings it in.
    view.zoom_slider.setValue(80)
    far = 5500 if view._win()[1] < 5000 else 100
    ed.navigate.emit(far)
    assert view._win()[0] <= far < view._win()[1]
    # The map's zoom is its own: a second map of the same sequence (and the sequence view) are not moved.
    other = MapView(ed)
    assert other.view is None
    from sequence_view import SequenceView
    sv = SequenceView(ed)
    size = sv.size
    view.zoom_slider.setValue(100)
    assert sv.size == size
    # A circular map has no zoom.
    circ = MapView(crowded(True, count=3))
    circ.resize(500, 500)
    wheel(circ, 250, 250)
    assert circ.view is None and circ.zoom_bar.isHidden()

