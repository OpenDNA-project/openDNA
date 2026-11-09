"""fileio.replace: the rename over a file that another thread still reads waits for it (Windows refuses it at once)."""
import os
import sys
import threading
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "app"))

import pytest

import fileio


def test_replace_waits_for_a_reader_to_let_go(tmp_path):
    target, tmp = tmp_path / "a.gb", tmp_path / "a.gb.tmp"
    target.write_text("old")
    tmp.write_text("new")
    reader = open(target, "rb")   # on Windows this blocks the rename (opened without delete sharing)
    threading.Timer(0.5, reader.close).start()
    start = time.time()
    fileio.replace(tmp, target)
    assert target.read_text() == "new" and not tmp.exists()
    assert os.name != "nt" or time.time() - start >= 0.4   # it did wait on Windows


def test_replace_does_not_wait_when_waiting_cannot_help(tmp_path):
    tmp = tmp_path / "x.tmp"
    tmp.write_text("new")
    (tmp_path / "folder").mkdir()
    start = time.time()
    with pytest.raises(OSError):
        fileio.replace(tmp, tmp_path / "folder")   # a folder: an error at once, not after 5 s
    assert time.time() - start < 2
