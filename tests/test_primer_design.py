"""Primers: the checks and their limits, where a primer comes from (the file, the database, or designed here) and what
that changes, designing a pair across the origin of a circle, and the files written for ordering."""
import csv
import os
import random
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ["OPENDNA_DATA_DIR"] = tempfile.mkdtemp(prefix="opendna_test_")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "app"))

import pytest

from opendna_core import Document

import feature_table
import primer_db
import primer_panel as PP
import tool_settings
from editor import Editor, reverse_complement
from history_panel import group
from hover import describe_primer
from primer_db import add_primer, db
from sequence_view import PRIMER_COLOR, TEMP_COLOR, primer_color

DEFAULT_LIMITS = {"tm_low": 52.0, "tm_high": 65.0, "tm_margin": 7.0, "gc_low": 40, "gc_high": 60, "gc_margin": 10,
                  "clamp_min": 1, "clamp_max": 3, "hairpin_bad": -6.0, "dimer3_bad": -8.0, "dimer_bad": -12.0}

rng = random.Random(21)
SEQ = "".join(rng.choice("ACGT") for _ in range(600))
DBONLY, MINE, BOTH, TEMP = SEQ[50:74], SEQ[150:174], SEQ[250:274], SEQ[350:374]


@pytest.fixture(autouse=True)
def private_primer_db(tmp_path, monkeypatch):
    """A primer database of this file's own: these tests add primers to it."""
    monkeypatch.setattr(primer_db, "_db", primer_db.PrimerDB(tmp_path / "primers.tsv"))


def random_dna(n, seed):
    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))


def editor_of(seq, circular=False, path="t.gb"):
    return Editor(Document.new_sequence("t", seq, circular), Path(path))


@pytest.fixture
def checks(monkeypatch):
    """primer_checks with the documented default limits, the default salt conditions, no genome and a Tm we choose."""
    limits = dict(DEFAULT_LIMITS)
    tm = {"value": 60.0}
    monkeypatch.setattr(PP, "limit", lambda key: limits[key])
    monkeypatch.setattr(PP, "conditions", lambda: {key: default for key, _, default, _ in PP.CONDITIONS})
    monkeypatch.setattr(PP, "melting_temp", lambda seq: tm["value"])
    monkeypatch.setattr(PP, "genome", lambda: None)
    ed = editor_of(random_dna(400, 9))

    def run(seq, name=None):
        result = PP.primer_checks(ed, seq)
        return {n: lv for n, lv, _ in result} if name is None else next(lv for n, lv, _ in result if n == name)
    run.limits, run.tm = limits, tm
    return run


def gc_primer(n_gc, n=20):
    body = list("G" * n_gc + "A" * (n - n_gc))
    random.Random(3).shuffle(body)
    return "".join(body)


def test_primer_check_levels_follow_their_limits(checks):
    """The documented defaults are what the settings declare and live in the data folder; the Tm, GC, 3' clamp, hairpin
    and self-dimer levels change exactly at their thresholds and move with the limits; the verdicts roll up and a short
    primer has none."""
    # --- the limits are the declared defaults, and the settings live in the data folder, not the registry ---
    from PySide6.QtCore import QSettings
    s = primer_db.settings()
    assert s.format() == QSettings.IniFormat and Path(s.fileName()) == primer_db.data_dir() / "primers.ini"
    assert PP.settings is primer_db.settings
    assert {p.key: p.default for p in tool_settings.TOOLS[PP.CHECKS]} == DEFAULT_LIMITS
    assert [PP.band(v, 52, 65, 7) for v in (51.9, 52, 65, 65.1, 72, 72.1, 45, 44.9)] == \
        ["warn", "ok", "ok", "warn", "warn", "bad", "warn", "bad"]

    # --- Tm ---
    seq = "ACGTAGCTAGCTAGGATCCA"
    for tm, level in [(51.9, "warn"), (52.0, "ok"), (58.0, "ok"), (65.0, "ok"), (65.1, "warn"), (72.0, "warn"), (72.1, "bad"),
                      (45.0, "warn"), (44.9, "bad"), (20.0, "bad")]:
        checks.tm["value"] = tm
        assert checks(seq, "Tm") == level, tm
    checks.limits["tm_low"], checks.limits["tm_margin"] = 55.0, 2.0      # the limits are the settings: moving them moves the verdict
    checks.tm["value"] = 53.0
    assert checks(seq, "Tm") == "warn"
    checks.tm["value"] = 52.9
    assert checks(seq, "Tm") == "bad"

    # --- GC ---
    checks.limits.update(DEFAULT_LIMITS)
    got = {k: checks(gc_primer(k), "GC") for k in (4, 5, 6, 7, 8, 12, 13, 14, 15, 20)}   # of 20 bases: 5 % per base
    assert got == {4: "bad", 5: "bad", 6: "warn", 7: "warn", 8: "ok", 12: "ok", 13: "warn", 14: "warn", 15: "bad", 20: "bad"}
    checks.limits["gc_low"], checks.limits["gc_margin"] = 45, 5
    assert [checks(gc_primer(k), "GC") for k in (8, 9, 10)] == ["warn", "ok", "ok"]
    assert checks("gcgcgcgcATATATAT", "GC") == "ok"   # lower case counts as the same bases (50 %)

    # --- 3' clamp ---
    checks.limits.update(DEFAULT_LIMITS)
    ends = {"ATATA": 0, "ATATG": 1, "AGCTA": 2, "AGCTG": 3, "CGCTG": 4, "GCGCG": 5}
    assert {n: checks("ACTGATGCATGACAT" + end, "3′ clamp") for end, n in ends.items()} == \
        {0: "warn", 1: "ok", 2: "ok", 3: "ok", 4: "warn", 5: "warn"}
    checks.limits["clamp_max"] = 4
    assert checks("ACTGATGCATGACAT" + "CGCTG", "3′ clamp") == "ok"

    # --- hairpin and self-dimer, and the exact edge where the limit crosses the primer's own ΔG ---
    checks.limits.update(DEFAULT_LIMITS)
    clean, loop, strong = "ATGCCGATTAGCCTAGCAAGTCCGA", "GGATCCAAAAAGGATCC", "GCGCGCGCAAAAGCGCGCGC"
    assert checks(clean, "Hairpin") == "ok" and checks(clean, "Self-dimer") == "ok"
    assert checks(loop, "Hairpin") == "warn" and checks(strong, "Hairpin") == "bad"     # usual limit (3' end < -2 or any < -3), then < -6
    assert checks("CAGTCGGAGTCCGACTG", "Self-dimer") == "warn" and checks("ACGCGTACGCGTACGCGT", "Self-dimer") == "bad"
    structure = PP.structure(strong)
    checks.limits["hairpin_bad"] = structure["hairpin"] - 0.01
    assert checks(strong, "Hairpin") == "warn"
    checks.limits["hairpin_bad"] = structure["hairpin"] + 0.01
    assert checks(strong, "Hairpin") == "bad"
    d = PP.structure("ACGCGTACGCGTACGCGT")
    checks.limits.update(dimer3_bad=d["self3"] - 0.01, dimer_bad=d["self"] - 0.01)
    assert checks("ACGCGTACGCGTACGCGT", "Self-dimer") == "warn"
    checks.limits["dimer3_bad"] = d["self3"] + 0.01
    assert checks("ACGCGTACGCGTACGCGT", "Self-dimer") == "bad"
    checks.limits.update(dimer3_bad=-50.0, dimer_bad=d["self"] + 0.01)     # only the "anywhere" limit
    assert checks("ACGCGTACGCGTACGCGT", "Self-dimer") == "bad"

    # --- the verdicts roll up and short primers have none ---
    assert PP.primer_checks(editor_of("ACGT" * 50), "ACGT") == [("Primer", "info", "Type at least 8 bases.")]
    assert PP.worst([("a", "ok", ""), ("b", "warn", ""), ("c", "info", "")]) == "warn"
    assert PP.worst([("a", "ok", ""), ("b", "bad", ""), ("c", "warn", "")]) == "bad"
    assert PP.worst([("a", "info", "")]) == "info" and PP.worst([]) == "info"
    html = PP.checks_html([("Tm", "ok", "60.0 °C"), ("GC", "bad", "10 %")])
    assert PP.MARKS["ok"][0] in html and PP.MARKS["bad"][0] in html and "60.0 °C" in html


def make(tmp_path, name="o.gb"):
    """An editor whose file stores MINE and BOTH; the database has DBONLY and BOTH; TEMP is designed here, not saved."""
    doc = Document.new_sequence("o", SEQ, True)
    doc.add_features([("primer_bind", "151..174", [("label", "mine"), ("note", "sequence: " + MINE)]),
                      ("primer_bind", "251..274", [("label", "both"), ("note", "sequence: " + BOTH)])])
    db().add("dbonly", DBONLY)
    db().add("both", BOTH)
    ed = Editor(doc, tmp_path / name)
    ed.set_extra([], [("temp", TEMP)])
    return ed


def sites(ed):
    return {s.label: s for s in ed.primer_sites}


def test_where_a_primer_comes_from_decides_how_it_is_kept(tmp_path):
    """The four origins (database only, file, both, designed here) get their flags, their colour and their place in the
    display filter; Import to This Sequence is offered for the ones that are not in the file and is one undo step; an
    imported, a designed and a newly made primer are all written to the file and come back as the file's own."""
    ed = make(tmp_path)
    # --- the flags and colour of each origin ---
    s = sites(ed)
    assert set(s) == {"dbonly", "mine", "both", "temp"}
    assert (s["dbonly"].from_db, s["dbonly"].from_file, s["dbonly"].provisional) == (True, False, True)
    assert (s["mine"].from_db, s["mine"].from_file, s["mine"].provisional) == (False, True, False)
    assert (s["both"].from_db, s["both"].from_file, s["both"].provisional) == (True, True, False), "in both: the file's style wins"
    assert (s["temp"].temporary, s["temp"].provisional) == (True, True)
    assert primer_color(s["dbonly"]) == primer_color(s["temp"]) == TEMP_COLOR, "the database primer looks like a temporary one"
    assert primer_color(s["mine"]) == primer_color(s["both"]) == PRIMER_COLOR

    # --- the display filter keeps each source apart ---
    ed.set_display(primers=True, db_primers=True)
    assert {x.label for x in ed.drawn_primers} == {"dbonly", "mine", "both", "temp"}
    ed.set_display(file_primers=False)
    assert {x.label for x in ed.drawn_primers} == {"dbonly"}
    ed.set_display(file_primers=True, db_primers=False)   # a designed primer goes with "this sequence"
    assert {x.label for x in ed.drawn_primers} == {"mine", "both", "temp"}
    ed.set_display(primers=False)
    assert ed.drawn_primers == []
    ed.set_display(primers=True, db_primers=True)

    # --- which primers can be imported, and importing one is a single undo step ---
    assert feature_table.primers_to_import(ed, ["dbonly", "temp", "both"]) == [("dbonly", DBONLY), ("temp", TEMP)]
    steps = len(ed.doc.history()[0])
    ed.import_primers([("dbonly", DBONLY)])
    assert len(ed.doc.history()[0]) == steps + 1, "one undo step"
    assert [t for _, _, t in group(ed.doc.history()[0])][-1] == "Imported primer dbonly to this sequence"
    site = sites(ed)["dbonly"]
    assert site.from_file and not site.provisional and primer_color(site) == PRIMER_COLOR, "now a regular primer annotation"
    assert "in this file" in describe_primer(ed, site.key)
    assert feature_table.primers_to_import(ed, ["dbonly"]) == [], "it is in the file now"
    assert ed.undo()
    assert sites(ed)["dbonly"].provisional and "not in this file" in describe_primer(ed, sites(ed)["dbonly"].key)
    assert ed.redo() and not sites(ed)["dbonly"].provisional

    # --- what is written to the file, and what comes back ---
    ed.import_primers([("temp", TEMP)])
    assert not sites(ed)["temp"].provisional, "a temporary primer becomes the file's too"
    ed.save(tmp_path / "o.gb")
    stored = {f.label for f in Document(str(tmp_path / "o.gb")).features() if f.kind == "primer_bind"}
    assert stored == {"mine", "both", "dbonly", "temp"}, "a primer only in the database would not be written"
    db().remove(["dbonly"])            # not even in the database any more: the file still has it
    back = Editor(Document(str(tmp_path / "o.gb")), tmp_path / "o.gb")
    s = sites(back)
    assert all(s[n].from_file and not s[n].provisional for n in ("dbonly", "temp", "mine", "both")), "all regular after reload"
    assert [x for x in sites(make(tmp_path, "p.gb")).values() if x.label == "temp" and x.provisional], \
        "a fresh editor on the unsaved document still has the designed primer as temporary"
    # a primer made here is the file's
    made = make(tmp_path, "m.gb")
    add_primer(made, "madehere", SEQ[450:474])
    assert not sites(made)["madehere"].provisional
    made.save(tmp_path / "m.gb")
    assert "madehere" in {f.label for f in Document(str(tmp_path / "m.gb")).features() if f.kind == "primer_bind"}


def test_primers_are_designed_for_a_target_at_or_across_the_origin_of_a_circle(tmp_path):
    """On a circular sequence the designer looks through the origin: a target at base 1, at the last base, or one that
    passes the origin gets pairs whose primers are the template's own bases and whose product holds the whole target. A
    linear sequence has nothing before base 1."""
    from primer_design import design
    seq = random_dna(3000, 5)
    n = len(seq)
    ed = Editor(Document.new_sequence("c", seq, True), tmp_path / "c.gb")
    for start, end in [(0, 200), (n - 200, n), (n - 100, n + 100)]:
        pairs = design(ed, start, end)
        assert pairs, (start, end)
        for p in pairs:
            f, r = p["fwd"], p["rev"]
            assert (seq * 2)[f["start"]:f["end"]] == f["seq"] and (seq * 2)[r["start"]:r["end"]] == reverse_complement(r["seq"])
            before, after = (start - f["end"]) % n, (r["start"] - end) % n   # template between each primer and the target
            assert before <= 150 and after <= 150
            assert p["size"] == len(f["seq"]) + before + (end - start) + after + len(r["seq"])
    linear = Editor(Document.new_sequence("l", seq, False), tmp_path / "l.gb")
    assert design(linear, 0, 200) == [] and design(linear, 500, 700)


def test_the_files_written_for_ordering(tmp_path):
    """The IDT plate file fills down its columns and starts a new plate after 96 oligos, carries the scale and the
    purification on every row, and the plain name/sequence file writes the names as they are; nothing is written for an
    empty order."""
    read = lambda path: list(csv.reader(open(path, encoding="utf-8", newline="")))
    rows = lambda n: [{"name": f"P{i}", "sequence": "ACGTACGTAC", "note": "", "added": ""} for i in range(n)]
    # --- the wells fill down the columns before moving right ---
    [out] = primer_db.export_order(rows(10), tmp_path / "o.csv")
    table = read(out)
    assert table[0] == ["Well Position", "Name", "Sequence", "Scale", "Purification"]
    assert [r[0] for r in table[1:]] == ["A1", "B1", "C1", "D1", "E1", "F1", "G1", "H1", "A2", "B2"]
    # --- more than 96 oligos continue on a second plate, starting again at A1 ---
    written = primer_db.export_order(rows(100), tmp_path / "o.csv")
    assert [w.name for w in written] == ["o_plate1.csv", "o_plate2.csv"]
    first, second = read(written[0]), read(written[1])
    assert len(first) - 1 == 96 and first[-1][0] == "H12"
    assert len(second) - 1 == 4 and second[1][0] == "A1" and second[1][1] == "P96"
    # --- the scale and the purification are on every row ---
    [out] = primer_db.export_order(rows(2), tmp_path / "s.csv", scale="100nm", purification="HPLC")
    assert {tuple(r[3:]) for r in read(out)[1:]} == {("100nm", "HPLC")}
    # --- odd names and degenerate bases ---
    odd = [{"name": 'fwd, "A" primer', "sequence": "ACGTACGT", "note": "", "added": ""}]
    [plate] = primer_db.export_order(odd, tmp_path / "p.csv")
    [plain] = primer_db.export_order(odd, tmp_path / "q.csv", primer_db.ORDER_FORMATS[1])
    assert read(plate)[1][1] == 'fwd, "A" primer' and read(plain) == [['fwd, "A" primer', "ACGTACGT"]]
    broken = [{"name": "a\tb\nc", "sequence": "ACGTACGT", "note": "", "added": ""}]
    [out] = primer_db.export_order(broken, tmp_path / "r.csv", primer_db.ORDER_FORMATS[1])
    assert read(out) == [["a b c", "ACGTACGT"]], "a tab or a line break in a name becomes a space"
    [out] = primer_db.export_order([{"name": "deg", "sequence": "acgtNNRYacgt", "note": "", "added": ""}],
                                   tmp_path / "d.csv", primer_db.ORDER_FORMATS[1])
    assert read(out) == [["deg", "ACGTNNRYACGT"]], "degenerate bases are exported as they are"
    # --- an empty order is refused and writes nothing ---
    empty = tmp_path / "empty"
    empty.mkdir()
    with pytest.raises(ValueError):
        primer_db.export_order([], empty / "o.csv")
    assert not os.listdir(empty)
