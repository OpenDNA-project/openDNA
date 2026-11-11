"""Repeats and the dot plot: the engine finds planted repeats at the right places and nothing in random DNA; Local Align's paired mode lists the
stretches two sequences share and draws their dot plot in both themes."""
import os
import random
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ["OPENDNA_DATA_DIR"] = tempfile.mkdtemp(prefix="opendna_test_")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "app"))

from PySide6.QtCore import QPoint, Qt
from PySide6.QtGui import QPalette
from PySide6.QtTest import QTest
from PySide6.QtWidgets import QApplication

from opendna_core import Document

import repeats
import theme
from editor import reverse_complement as rc

app = QApplication.instance() or QApplication([])
rng = random.Random(11)


def dna(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def kinds(result, kind):
    return [r for r in result.repeats if r.kind == kind]


def test_planted_direct_inverted_and_tandem_repeats_are_found_at_their_places():
    """A direct repeat of 3 copies, an inverted pair and a tandem array come back with every copy at its place."""
    d, h, unit = dna(30), dna(16), "GATTACAGCT"
    seq = dna(100) + d + dna(150) + d + dna(80) + d + dna(60) + h + dna(90) + rc(h) + dna(70) + unit * 5 + dna(100)
    found = repeats.find_repeats(seq)
    (direct,) = kinds(found, "direct")
    # (a random neighbouring base may match by chance and lengthen a repeat by a base: positions are checked within that)
    near = lambda got, want: len(got) == len(want) and all(abs(g - w) <= 2 for g, w in zip(got, want))
    assert 30 <= direct.length <= 32 and near([c[0] for c in direct.copies], [seq.find(d), seq.find(d, 140), seq.rfind(d)])
    assert all(d in seq[s - 2:e + 2] for s, e, _ in direct.copies)
    (inv,) = kinds(found, "inverted")
    assert 16 <= inv.length <= 18 and near([c[0] for c in inv.copies], [seq.find(h), seq.find(rc(h))])
    assert [c[2] for c in inv.copies] == [1, -1]
    (tandem,) = kinds(found, "tandem")
    assert tandem.unit == 10 and len(tandem.copies) == 5 and abs(tandem.copies[0][0] - seq.find(unit)) <= 2
    assert not found.capped


def test_random_dna_has_no_repeat_at_the_default_thresholds():
    assert not repeats.find_repeats(dna(20_000)).repeats
    assert not repeats.find_repeats(dna(5_000), circular=True).repeats
    assert len(repeats.find_repeats(dna(20_000), max_loop=None).repeats) > 0, "a 12-base stem anywhere in 20 kb is chance, hence the loop limit"


def test_thresholds_and_mismatches_decide_what_is_listed():
    """A 14-base direct repeat is below the default 20 and found at 12; one difference needs a mismatch allowance."""
    d = dna(14)
    seq = dna(80) + d + dna(80) + d + dna(80)
    assert not kinds(repeats.find_repeats(seq), "direct")
    assert len(kinds(repeats.find_repeats(seq, min_direct=12), "direct")) == 1
    long_one = dna(30)
    bent = long_one[:15] + ("A" if long_one[15] != "A" else "C") + long_one[16:]
    seq = dna(80) + long_one + dna(80) + bent + dna(80)
    assert not [r for r in kinds(repeats.find_repeats(seq), "direct") if r.length >= 30]
    (r,) = kinds(repeats.find_repeats(seq, mismatches=1), "direct")
    assert r.length >= 30 and r.mismatches >= 1


def test_a_repeat_over_the_origin_of_a_circle_is_found():
    d = dna(25)
    near = lambda got, want: len(got) == len(want) and all(abs(g - w) <= 2 for g, w in zip(got, want))
    seq = d[10:] + dna(100) + d + dna(100) + d[:10]
    (r,) = kinds(repeats.find_repeats(seq, circular=True), "direct")
    assert 25 <= r.length <= 28 and near(sorted(c[0] for c in r.copies), [115, 240])


def test_a_long_input_is_capped_and_says_so():
    result = repeats.find_repeats("ACGT" * 200_000)   # nothing but a microsatellite: every k-mer is skipped
    assert result.capped and "skipped" in result.capped


def test_pair_repeats_find_the_shared_stretches_in_both_orientations_and_a_sequence_against_itself():
    d, e = dna(40), dna(30)
    a = dna(50) + d + dna(60) + e + dna(50)
    b = dna(30) + rc(e) + dna(40) + d + dna(20)
    found = {r[0]: r for r in repeats.pair_repeats(a, b, 20).repeats}
    near = lambda got, want: all(abs(g - w) <= 2 for g, w in zip(got, want))   # (a random flank base may extend a stretch)
    assert near(found["direct"][1:], (50, 100, 40)) and near(found["inverted"][1:], (150, 30, 30))
    assert not repeats.pair_repeats(dna(3000), dna(3000), 20).repeats
    seq = dna(100) + d + dna(100) + d + dna(100)
    (r,) = repeats.pair_repeats(seq, seq, 30).repeats   # itself: the pair once, not the main diagonal
    assert r[0] == "direct" and near(r[1:4], (100, 240, 40))
    direct, inverted, _ = repeats.pair_dots(a, b)
    assert direct and inverted


def test_local_align_pair_mode_aligns_the_two_boxes_and_the_dot_plot_draws_in_both_themes():
    import main as M
    from local_align import LocalAlignDialog
    d = dna(60)
    a, b = dna(100) + d + dna(100), dna(80) + rc(d) + dna(120)
    win = M.MainWindow()
    win.show()
    win.show_document(Document.new_sequence("t", dna(500), True), Path(tempfile.mkdtemp()) / "t.gb", untitled=True)
    dlg = LocalAlignDialog(win, win.editor)
    dlg.query.setPlainText(a)
    dlg.paired.setChecked(True)
    assert not dlg.query_b.isHidden() and not dlg.tabs.isTabVisible(0) and dlg.tabs.isTabVisible(3)
    dlg.run()   # the second box is empty: said, nothing drawn
    assert "empty" in dlg.status.text() and dlg.canvas.alignment is None
    dlg.query_b.setPlainText(b)
    dlg.run()
    assert dlg.canvas.alignment.names == ["A", "B"] and "identity" in dlg.status.text()
    pair, plot = dlg.pair, dlg.pair.plot
    assert [pair.tree.topLevelItem(i).text(0) for i in range(pair.tree.topLevelItemCount())] == ["Inverted"]
    assert plot.inverted and (plot.length, plot.rows) == (len(a), len(b))
    plot.resize(400, 400)
    before = plot.grab().toImage()
    pal = QApplication.palette()
    dark = QPalette(pal)
    dark.setColor(QPalette.Window, Qt.black)
    QApplication.setPalette(dark)
    theme.refresh_views()
    assert theme.is_dark() and plot.grab().toImage() != before
    QApplication.setPalette(pal)
    theme.refresh_views()
    rect = plot.plot_rect()
    got = []
    plot.picked.connect(lambda x, y: got.append((x, y)))
    QTest.mouseClick(plot, Qt.LeftButton, Qt.NoModifier, QPoint(int(rect.left() + rect.width() * 0.5), int(rect.top() + rect.height() * 0.5)))
    assert got and abs(got[0][0] - len(a) // 2) <= len(a) // 100 and abs(got[0][1] - len(b) // 2) <= len(b) // 100
    dlg.paired.setChecked(False)
    assert dlg.tabs.isTabVisible(0) and not dlg.tabs.isTabVisible(3)
    win.editor.doc.modified = False
    win.reset_layout()
    win.close()
