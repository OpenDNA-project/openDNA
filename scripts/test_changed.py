"""Run only the tests a change can affect, in parallel (one worker per test file), instead of the whole suite.

    python scripts/test_changed.py              # what differs from main, plus uncommitted files
    python scripts/test_changed.py app/cas9.py  # or name the files yourself
    python scripts/test_changed.py --all        # everything, still in parallel

A changed tests/test_x.py runs itself; a changed app/m.py runs every test file that imports the module m (and tests/test_m.py);
a change in core/, conftest.py or pytest.ini runs everything. It is a shortcut: CI still runs the whole suite."""
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TESTS = sorted((ROOT / "tests").glob("test_*.py"))


def git(*args):
    return subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True).stdout.split()


def pick(changed):
    chosen = set()
    for path in changed:
        p = Path(path)
        if path.startswith("core/") or p.name in ("conftest.py", "pytest.ini", "fake_reads.py", "sample_files.py"):
            return TESTS
        if path.startswith("tests/test_") and (ROOT / path).exists():
            chosen.add(ROOT / path)
        elif path.startswith("app/") and p.suffix == ".py":
            uses = re.compile(rf"^\s*(import {p.stem}\b|from {p.stem} import)", re.M)
            chosen |= {t for t in TESTS if t == ROOT / "tests" / f"test_{p.stem}.py" or uses.search(t.read_text(encoding="utf-8"))}
    return sorted(chosen)


if __name__ == "__main__":
    args = sys.argv[1:]
    changed = [a for a in args if a != "--all"] or git("diff", "--name-only", "main") + git("ls-files", "--others", "--exclude-standard")
    files = TESTS if "--all" in args else pick(changed)
    if not files:
        sys.exit("No test file is affected by: " + ", ".join(changed))
    print(f"{len(files)} of {len(TESTS)} test files: " + " ".join(f.name for f in files))
    sys.exit(subprocess.run([sys.executable, "-m", "pytest", "-q", "-n", "8", "--dist", "loadfile", *map(str, files)], cwd=ROOT).returncode)
