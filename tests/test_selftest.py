"""Runs every in-app self-test check (Tools → Settings → Self-test) under pytest."""
import os
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ["GENOMEEDIT_DATA_DIR"] = tempfile.mkdtemp(prefix="genomeedit_test_")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "app"))

from PySide6.QtWidgets import QApplication

import selftest

app = QApplication.instance() or QApplication([])


def test_selftest_checks():
    """Every in-app self-test check passes (all failures are listed together). The panel that runs them in the app is tested
    with fake checks in test_app_logic.py::test_self_test_panel_counts_passes_and_failures_and_shows_the_message."""
    assert len(selftest.CHECKS) >= 9, "the 9 checks (November 2026) are registered"
    failed = []
    for area, what, fn in selftest.CHECKS:
        ok, _, message = selftest.run_check(fn)
        if not ok:
            failed.append(f"{area}: {what}: {message}")
    assert not failed, f"{len(failed)} self-test check(s) failed:\n" + "\n".join(failed)
