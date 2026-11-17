"""File tabs without a window: the bookkeeping of the files of one window, and that two open files share nothing."""
import os
import sys
import tempfile
import json
import subprocess
import time
from pathlib import Path
from types import SimpleNamespace

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ["OPENDNA_DATA_DIR"] = tempfile.mkdtemp(prefix="opendna_test_")
APP = Path(__file__).resolve().parent.parent / "app"
sys.path.insert(0, str(APP))

from opendna_core import Document

from editor import Editor
from file_tabs import FileSet
from single_instance import Listener, send_to_running


def page(path=None, unsaved=False):
    return SimpleNamespace(file_path=lambda: path, has_unsaved=lambda: unsaved)


def test_file_set_keeps_order_dedupes_and_picks_the_next_page_on_close(tmp_path):
    a, b, c = (tmp_path / n for n in ("a.gb", "b.gb", "c.gb"))
    for f in (a, b, c):
        f.write_text("")
    pa, pb, pc, untitled = page(a), page(b, unsaved=True), page(c), page(None, unsaved=True)
    files = FileSet()
    for p in (pa, pb, pc, untitled):
        files.add(p)
    assert files.current is untitled                                  # the last added is shown
    assert files.find(b) is pb and files.find(tmp_path / "x.gb") is None
    assert files.find(tmp_path / "sub" / ".." / "a.gb") is pa         # the same file, spelled another way
    assert files.find(str(c).upper()) is (pc if os.name == "nt" else None)   # case only matters on Windows
    assert files.unsaved() == [pb, untitled]
    files.move(0, 2)                                                  # drag a.gb to the third place
    assert files.pages == [pb, pc, pa, untitled]
    assert files.step(1) is pb and files.step(-1) is pa               # from the shown last page: round the end, and back one
    files.move(2, 0)
    assert files.close(untitled) is pc                                # the shown (last) page closes: the one to its left
    files.current = pb
    assert files.close(pb) is pc                                      # the shown page closes: the one that took its place
    assert files.close(pa) is pc                                      # another page closes: the shown one stays
    assert files.close(pc) is None and files.pages == []


def test_two_open_files_share_no_editor_state():
    one = Editor(Document.new_sequence("one", "ACGT" * 50, False), Path("one.gb"))
    two = Editor(Document.new_sequence("two", "TTTT" * 50, False), Path("two.gb"))
    before = two.version
    one.set_cursor(0)
    one.type_bases("GGGG")
    assert one.version > before and two.version == before
    assert one.length == 204 and two.length == 200
    assert (one.start, one.end) == (4, 4) and (two.start, two.end) == (0, 0)
    assert two.doc.sequence() == "TTTT" * 50


def test_a_second_launch_hands_its_files_to_the_running_program(tmp_path):
    from PySide6.QtWidgets import QApplication
    app = QApplication.instance() or QApplication([])
    assert send_to_running(["x.gb"], timeout=200) is False, "nothing is running: the launch goes on to start the program"
    listener = Listener()
    assert listener.listening
    got = []
    listener.files.connect(got.append)
    files = [str(tmp_path / "a.gb"), str(tmp_path / "b c.gb"), "ünï.gb"]
    # the second launch is a process of its own, as in life (it waits for the answer, which this process's loop gives)
    code = ("import json, sys; sys.path.insert(0, sys.argv[1]); from PySide6.QtWidgets import QApplication; app = QApplication([]); "
            "import single_instance; sys.exit(0 if single_instance.send_to_running(json.loads(sys.argv[2]), 2000, 8000) else 1)")
    sender = subprocess.Popen([sys.executable, "-c", code, str(APP), json.dumps(files)])
    for _ in range(1000):
        app.processEvents()
        if sender.poll() is not None:
            break
        time.sleep(0.01)
    assert sender.wait(10) == 0
    assert got == [files], "every file arrives, in order, whatever the characters in the names"
    listener.server.close()
