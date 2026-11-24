"""CHANGELOG.md: what changed in each version, newest first. The start page shows the running version's section and a window with the whole file.

Each version is "## x.y.z (date)" with the bold labels Features, Bugfixes and UI/UX (tests/test_project.py checks the shape); the top section is
"## Unreleased", where every change adds its line (CLAUDE.md) and which is never shown to users.
"""
import re
from pathlib import Path

# next to this file in the bundled exe (packaging/openDNA.spec), at the repo root from source
FILE = next((p for p in (Path(__file__).with_name("CHANGELOG.md"), Path(__file__).parent.parent / "CHANGELOG.md") if p.is_file()), None)
CATEGORIES = ("Features", "Bugfixes", "UI/UX")


def text():
    """The whole file ("" if it is missing: the page then shows nothing)."""
    return FILE.read_text(encoding="utf-8") if FILE else ""


def section(version, source=None):
    """The body of the "## version (date)" section (markdown), or None when there is none."""
    m = re.search(rf"^## {re.escape(version)}\b[^\n]*\n(.*?)(?=^## |\Z)", text() if source is None else source, re.S | re.M)
    return m.group(1).strip() if m else None
