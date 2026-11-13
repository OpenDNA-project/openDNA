"""The Features Database: sources (one file each, your own ones the source "Local", which can be removed again), the strength a
source measured, and switching a source off; and the data-quality rules every shipped part must meet."""
import os
import sys
import tempfile
from pathlib import Path

import pytest

os.environ["OPENDNA_DATA_DIR"] = tempfile.mkdtemp(prefix="opendna_test_")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "app"))

import parts_library as pl
import tool_settings

ROOT = Path(__file__).resolve().parent.parent


@pytest.fixture(autouse=True)
def own_library(tmp_path, monkeypatch):
    """A collection folder of this test's own, and no user library, so the bundled files are not touched."""
    monkeypatch.setattr(pl, "BUNDLED_DIR", tmp_path / "parts")
    (tmp_path / "parts").mkdir()
    monkeypatch.setattr(pl, "user_library", lambda: tmp_path / "mine.tsv")
    # The data folder is shared by the whole suite: put the switches back as they were.
    before = pl.off_ids()
    for cid in list(before):
        pl.set_on(cid, True)
    yield
    tool_settings.store().setValue(f"{pl.SETTING[0]}/{pl.SETTING[1]}", ",".join(sorted(before)))


def write(name, body):
    (pl.BUNDLED_DIR / f"{name}.tsv").write_text(body, encoding="utf-8")


MEASURED = """#title: Terminators (Someone 2013)
#reference: Someone A, et al. J Test 1, 1-2 (2013).
#doi: 10.0000/test
#note: Termination efficiency measured in E. coli.
name\ttype\tsequence\torganism\tstrength\tunits\tsource
Tstrong\tterminator\tGCCGCTTAAGCGGCTTTTTT\tE. coli\t88\tTE\tpFAB1
Tweak\tterminator\tGCCGCTTAAGCGGCATTTTA\tE. coli\t3.5\tTE\tpFAB2
"""
PLAIN = """#title: Cloning vectors (NCBI records)
name\ttype\tsequence\torganism\tstrength\tunits\tsource
AmpR\tCDS\tATGAGTATTCAACATTTCCG\t\t\t\tpBR322 (J01749)
"""


def clear_library():
    """Start the next section of a merged test from an empty library (no collection files, no own parts, all switches on)."""
    for f in pl.BUNDLED_DIR.glob("*.tsv"):
        f.unlink()
    if pl.user_library().exists():
        pl.user_library().unlink()
    pl.set_on(pl.YOURS, True)


def test_collections_sources_and_switches():
    """A collection carries its source and measurements; a source switched off is used by nothing; your own parts are a
    collection too (the source "Local", keeping the file they came from), and only they can be removed again, the file keeping
    its header; a broken row is skipped, not fatal."""
    # --- source and measurements ---
    write("measured", MEASURED)
    write("vectors", PLAIN)
    by_id = {c.id: c for c in pl.collections()}
    assert set(by_id) == {"measured", "vectors"}
    c = by_id["measured"]
    assert c.title == "Terminators (Someone 2013)" and c.doi == "10.0000/test" and "Someone A" in c.reference
    assert c.counts == {"terminator": 2} and c.measured == 2
    assert by_id["vectors"].counts == {"CDS": 1} and by_id["vectors"].measured == 0
    parts = {p.name: p for p in pl.load_parts()}
    assert parts["Tstrong"].strength == 88 and parts["Tstrong"].units == "TE"
    assert parts["Tstrong"].organism == "E. coli" and parts["Tstrong"].collection == "measured"
    assert pl.describe(parts["Tstrong"]) == "TE 88 (E. coli)"
    assert pl.describe(parts["AmpR"]) == "", "nothing measured, nothing claimed"
    # --- a source switched off ---
    write("measured", MEASURED)
    write("vectors", PLAIN)
    assert len(pl.load_parts()) == 3
    pl.set_on("measured", False)
    assert [p.name for p in pl.load_parts()] == ["AmpR"] and not pl.is_on("measured")
    assert [(c.id, c.on) for c in pl.collections()] == [("measured", False), ("vectors", True)], "the dialog's switches"
    pl.set_on("measured", True)
    assert len(pl.load_parts()) == 3
    # --- your own parts ---
    clear_library()
    write("vectors", PLAIN)
    assert pl.add_parts([("mine", "promoter", "ttgacaattaatcatcggct", "my.gb"),
                         ("again", "promoter", "TTGACAATTAATCATCGGCT", "my.gb"),   # the same sequence, skipped
                         ("AmpR copy", "CDS", "ATGAGTATTCAACATTTCCG", "my.gb")]) == 1
    mine = next(c for c in pl.collections() if c.id == pl.YOURS)
    assert mine.title == pl.LOCAL and mine.counts == {"promoter": 1} and mine.measured == 0
    assert {p.name for p in pl.load_parts()} == {"AmpR", "mine"}
    pl.set_on(pl.YOURS, False)
    assert {p.name for p in pl.load_parts()} == {"AmpR"}
    # --- Local features keep their file and only they can be removed ---
    clear_library()
    write("vectors", PLAIN)
    pl.add_parts([("mine", "promoter", "TTGACAATTAATCATCGGCT", "my.gb"), ("other", "RBS", "AAGGAGGTAAAAAATG", "b.gb")])
    by_name = {p.name: p for p in pl.all_parts()}
    assert (by_name["mine"].source_name, by_name["mine"].source) == (pl.LOCAL, "my.gb")
    assert by_name["AmpR"].source_name == "Cloning vectors (NCBI records)"
    assert pl.remove_parts([("mine", "TTGACAATTAATCATCGGCT"), ("AmpR", "ATGAGTATTCAACATTTCCG")]) == 1, "only Local"
    assert {p.name for p in pl.all_parts()} == {"AmpR", "other"}
    assert pl.user_library().read_text(encoding="utf-8").startswith(f"#title: {pl.LOCAL}\n"), "the header is kept"
    # --- a broken row ---
    clear_library()
    write("odd", "#title: Odd\nname\ttype\tsequence\torganism\tstrength\tunits\tsource\n"
                 "no sequence\tpromoter\t\t\t\t\t\n"
                 "protein\tCDS\tMKPG\t\t\t\t\n"
                 "bad number\tterminator\tACGTACGTACGTACGTACGT\t\tvery strong\tTE\t\n")
    parts = pl.load_parts()
    assert [p.name for p in parts] == ["bad number"], "only the row with real DNA"
    assert parts[0].strength is None and pl.describe(parts[0]) == ""


# ---- the collections that actually ship ---------------------------------------------------------
# A GenBank submitter annotates their own construct, not a parts catalogue, so the build scripts curate
# what they take (scripts/build_parts.py, docs/PARTS_LIBRARY_METHOD.md). These tests hold that line:
# the rules are enforced here, not merely written down.

sys.path.insert(0, str(ROOT / "scripts"))
import build_parts  # noqa: E402


@pytest.fixture
def shipped(monkeypatch):
    monkeypatch.setattr(pl, "BUNDLED_DIR", ROOT / "app" / "parts")
    monkeypatch.setattr(pl, "user_library", lambda: ROOT / "nonexistent.tsv")
    return pl.load_parts()


def test_shipped_parts_meet_the_data_quality_rules(shipped):
    """The collections actually shipped: every file parses into a titled, non-empty source; no part is junk (primer sites,
    cloning sites and slices of a part are not parts: the primer database draws primers, the digest view finds cutters, and five
    overlapping pieces of one FRT site light up five times; nor a type you cannot name); every part clears the length floor of
    its type (a phage promoter really is 19 bp, an 18 bp operator matches by chance); no part twice by name or by sequence,
    inside a collection or across two (six KanR rows and three names for one origin were the complaint); and every part cites
    the record it came from (attribution is the licence condition for the whole library)."""
    # --- every shipped file parses ---
    found = pl.collections()
    assert found and all(c.title and sum(c.counts.values()) for c in found)
    assert len(shipped) == sum(sum(c.counts.values()) for c in found)
    for p in shipped:
        # --- not junk ---
        assert not build_parts.NEVER.search(p.name), f"{p.name!r} ({p.collection}) is not a reusable part"
        assert p.kind not in ("misc_feature", "regulatory", "sig_peptide"),             f"{p.name!r} has type {p.kind}: a part whose type you cannot name is one you have not understood"
        # --- the floor for its type ---
        floor = build_parts.FLOORS.get(p.kind)
        assert floor, f"{p.kind!r} ({p.name}) has no agreed minimum length"
        assert len(p.sequence) >= floor, f"{p.name!r} is {len(p.sequence)} bp, under the {floor} bp floor for {p.kind}"
        # --- says where it came from ---
        assert p.source.strip(), f"{p.name!r} cites no record"
    # --- one part, one row (identical sequences collapse at load time, so the raw rows are checked too) ---
    names = [p.name for p in shipped]
    assert len(names) == len(set(names)), f"repeated names: {sorted({n for n in names if names.count(n) > 1})}"
    raw = {}
    for cid, path in pl.files():
        for row in pl._meta_and_rows(path)[1]:
            seq = (row.get("sequence") or "").upper()
            if seq:
                raw.setdefault(seq, []).append(f"{row['name']} ({cid})")
    both = {s: who for s, who in raw.items() if len(who) > 1}
    assert not both, f"the same sequence is stored more than once: {list(both.values())[:5]}"
