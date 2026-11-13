"""Runs every in-app self-test check (Help > Troubleshooting > Self-Test) under pytest."""
import os
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ["OPENDNA_DATA_DIR"] = tempfile.mkdtemp(prefix="opendna_test_")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "app"))

import selftest


def test_every_selftest_check_passes():
    """Every in-app self-test check is registered and passes (all failures are listed together)."""
    assert len(selftest.CHECKS) >= 9, "the 9 checks (November 2026) are registered"
    failed = []
    for area, what, fn in selftest.CHECKS:
        ok, _, message = selftest.run_check(fn)
        if not ok:
            failed.append(f"{area}: {what}: {message}")
    assert not failed, f"{len(failed)} self-test check(s) failed:\n" + "\n".join(failed)
