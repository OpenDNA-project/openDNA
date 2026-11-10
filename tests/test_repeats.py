"""Repeats and the dot plot: the engine finds planted repeats at the right places and nothing in random DNA; the tool lists them,
selects a copy on a click, and the plot draws in both themes."""
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


def test_the_tool_lists_copies_selects_on_a_click_and_the_plot_draws_in_both_themes():
    import main as M
    d, e = dna(30), dna(24)
    seq = dna(100) + d + dna(200) + rc(d) + dna(100) + e + dna(100) + e + dna(100)
    win = M.MainWindow()
    win.show()
    win.show_document(Document.new_sequence("t", seq, True), Path(tempfile.mkdtemp()) / "t.gb", untitled=True)
    tool = repeats.RepeatsDialog(win, win.editor)
    assert tool.tree.topLevelItemCount() == 2 and "2 repeats" in tool.status.text()
    top = [tool.tree.topLevelItem(i) for i in range(2)]
    top[0].child(1).setSelected(True)
    tool.tree.setCurrentItem(top[0].child(1))
    assert win.editor.end - win.editor.start >= 30 and seq[win.editor.start:win.editor.end].find(rc(d)) >= 0
    # the plot: a dot on the anti-diagonal for the inverted repeat, picture changes with the theme
    plot = tool.plot
    plot.resize(400, 400)
    assert plot.inverted and plot.length == len(seq)
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
    at = QPoint(int(rect.left() + rect.width() * 0.5), int(rect.top() + rect.height() * 0.5))
    QTest.mouseClick(plot, Qt.LeftButton, Qt.NoModifier, at)
    assert abs(win.editor.start - len(seq) // 2) <= len(seq) // 100
    win.editor.doc.modified = False
    win.reset_layout()
    win.close()
