"""The long runs of the editing / annotation-drift property tests (docs/TESTING.md sections 1
and 2). The machinery lives in app/selftest.py, which runs short versions in the app."""
import os
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ["GENOMEEDIT_DATA_DIR"] = tempfile.mkdtemp(prefix="genomeedit_test_")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "app"))

import selftest


def test_500_random_edits_keep_bases_and_features_then_undo_redo_exactly():
    """500 random edits (insert, delete, replace, overwrite, paste with features, feature
    add / delete, reverse complement; rotate and edits across the origin on plasmids): after
    each one the bases equal a plain string and every feature covers exactly the bases the
    per-base identity model says. Then undo everything and redo everything, snapshot by snapshot.
    Three seeds, linear and circular (the app's selftest runs more)."""
    for circular in (False, True):
        for seed in range(3):
            steps = selftest.run_random_edits(seed, circular, 500)
            assert steps == 500, (seed, circular)
