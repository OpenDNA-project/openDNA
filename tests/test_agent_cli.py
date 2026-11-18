"""The command line for agents (app/cli.py and the recipes in app/agent/): what an agent reads from a plasmid must be what the
GenBank file says, a mistake must be one `ERROR:` line, and the example in every recipe must still run."""
import json
import os
import re
import shlex
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ["OPENDNA_DATA_DIR"] = tempfile.mkdtemp(prefix="opendna_test_")
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "app"))

import cli
import mutation
from opendna_core import Document

DEMO = str(ROOT / "demo" / "pBAD30_eilR.gb")


def run(capsys, *argv):
    code = cli.main([str(a) for a in argv])
    out = capsys.readouterr()
    return code, out.out, out.err


def test_agent_reads_the_plasmid(capsys):
    doc = Document(DEMO)
    code, out, _ = run(capsys, "inspect", DEMO)
    info = json.loads(out)
    assert code == 0 and (info["length"], info["circular"]) == (doc.length, doc.circular)
    assert len(info["features"]) == len(doc.features())
    # a reverse-strand CDS (araC, 3655..4533) is read from its own start codon: the reverse complement of the file's bases
    code, out, _ = run(capsys, "feature", DEMO, "araC")
    f = json.loads(out)
    assert (f["start"], f["end"], f["strand"]) == (3655, 4533, -1)
    assert f["sequence"] == mutation.reverse_complement(doc.sequence(3654, 4533)) and f["sequence"].startswith("ATG")
    assert f["protein"]["sequence"].startswith("M") and f["protein"]["sequence"].endswith("*") and f["protein"]["length"] == (4533 - 3655 + 1) // 3 - 1
    # region is 1-based inclusive, like the file
    code, out, _ = run(capsys, "region", DEMO, 1, 12)
    assert json.loads(out)["sequence"] == doc.sequence(0, 12)


def test_agent_mistakes_are_one_error_line(capsys):
    for argv, words in [(("inspect", "nofile.gb"), "File not found"), (("feature", DEMO, "nope"), "No feature called"),
                        (("region", DEMO, 0, 12), "Positions must satisfy"), (("region", DEMO, 5, 99999), "Positions must satisfy"),
                        (("skills", "show", "nope"), "No recipe called")]:
        code, out, err = run(capsys, *argv)
        assert code == 1 and out == "" and err.startswith("ERROR: ") and words in err and err.count("\n") == 1


def test_every_recipe_example_still_runs(capsys):
    """The `## Example` block of each recipe is a list of commands, one per line: all of them must succeed (a recipe that drifted from
    the tool is caught here)."""
    code, out, _ = run(capsys, "skills", "list")
    names = [line.split(":")[0] for line in out.splitlines()]
    assert code == 0 and {"index", "inspect", "kld", "pcr"} <= set(names)
    for name in names:
        code, text, _ = run(capsys, "skills", "show", name)
        assert code == 0
        for block in re.findall(r"## Example\s+```\n(.*?)```", text, re.S):
            for line in block.splitlines():
                code, _, err = run(capsys, *shlex.split(line.replace("demo/", (ROOT / "demo").as_posix() + "/").replace("app/agent/", (ROOT / "app" / "agent").as_posix() + "/")))
                assert code == 0, f"{name}: {line}: {err}"


# ---- mutations in HGVS and the KLD batch ---------------------------------------------------------------------------------------

import csv
import random

import kld
import pytest
from user_error import UserError

CDS = "ATGGCTGAAAAATTTGGTTAA"   # M A E K F G *


def dna(n, seed):
    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))


def plasmid():
    """521 bp circle: gx = CDS at 101..121 (forward), gy = the same CDS on the reverse strand at 401..421."""
    seq = dna(100, 1) + CDS + dna(279, 2) + mutation.reverse_complement(CDS) + dna(100, 3)
    doc = Document.new_sequence("t", seq, True)
    doc.add_features([("CDS", "101..121", [("label", "gx")]), ("CDS", "complement(401..421)", [("label", "gy")])])
    return doc, seq


def change(doc, spec):
    r = mutation.resolve(doc, spec)
    return r.start, r.end, r.text, r.hgvs_dna, r.hgvs_protein, r.effect


def test_hgvs_mutations_become_the_change_they_name():
    doc, seq = plasmid()
    # hand-computed: gx is bases 101..121, codon 2 (GCT) at 104..106, codon 3 (GAA) at 107..109, c.4 = base 104
    assert change(doc, "gx:c.4G>A") == (103, 104, "A", "o.104G>A", "p.(Ala2Thr)", "missense")
    assert change(doc, "gx:p.Ala2Thr")[:4] == (103, 104, "A", "o.104G>A")       # the protein form picks the codon with the fewest changes
    assert change(doc, "gx:A2T")[:3] == (103, 104, "A") and change(doc, "o.104G>A")[:3] == (103, 104, "A")   # one letter; plasmid numbers
    assert change(doc, "gx:p.Glu3del") == (106, 109, "", "o.107_109del", "p.(Glu3del)", "in-frame deletion")
    assert change(doc, "gx:c.9_10insCCC") == (109, 109, "CCC", "o.109_110insCCC", "p.(Glu3_Lys4insPro)", "in-frame insertion")
    assert change(doc, "gx:c.9_10insC")[5] == "frameshift"
    assert change(doc, "gx:c.7_9dup") == (109, 109, "GAA", "o.109_110insGAA", "p.(Glu3dup)", "in-frame duplication")
    assert change(doc, "gx:p.Glu3_Lys4delinsGln")[:3] == (106, 110, "C") and change(doc, "gx:c.7_12delinsCAA")[:3] == (106, 110, "C")   # (GAAAAA -> CAA: only GAAA -> C differs)
    assert change(doc, "gx:p.Ala2Ter")[5] == "nonsense"
    # the reverse-strand copy: c. numbers and bases are read on the gene's own strand, the plasmid change is on the forward one
    assert change(doc, "gy:c.4G>A") == (417, 418, "T", "o.418C>T", "p.(Ala2Thr)", "missense")
    assert change(doc, "gy:p.Ala2Thr")[:3] == (417, 418, "T")
    assert change(doc, "gy:p.Glu3del")[:3] == (412, 415, "")
    assert change(doc, "gy:c.9_10insCCC")[:3] == (412, 412, "GGG")


def test_hgvs_mistakes_say_what_to_fix_and_ambiguities_only_warn():
    doc, _ = plasmid()
    for spec, words in [("gx:c.4C>A", "not C"), ("gx:p.Lys2Thr", "is Ala (GCT), not Lys"), ("gx:p.Ala2Ala", "already"), ("o.600A>T", "outside the plasmid"),
                        ("gx:c.50G>A", "past the end"), ("nope:c.1A>T", "No feature called"), ("gx:c.9_11insA", "neighbouring"),
                        ("gx:c.0A>T", "no position 0"), ("c.4G>A", "name the feature"), ("lacZ del 5", "Cannot read"), ("gx:c.*3A>T", "stop codon")]:
        with pytest.raises(UserError, match=re.escape(words)):
            mutation.resolve(doc, spec)
    r = mutation.resolve(doc, "gx:p.Ala2Gly")   # GCT -> GGT is one change; GGC, GGA, GGG are two: no tie, no warning
    assert (r.start, r.end, r.text) == (104, 105, "G") and not any("chosen" in w for w in r.warnings)
    r = mutation.resolve(doc, "gx:p.Ala2Ser")   # GCT -> TCT or AGT (one change each)... TCT, TCC, TCA, TCG, AGT, AGC: TCT is nearest
    assert (r.start, r.end, r.text) == (103, 104, "T")
    r = mutation.resolve(doc, "gx:p.Glu3Asp")   # GAA -> GAT or GAC: two codons one change apart: a default (alphabetical without a usage table) and a warning naming the other
    assert r.text == "C" and any("GAC chosen" in w and "GAT" in w for w in r.warnings)   # (no usage table: alphabetical)


def test_kld_batch_gives_the_files_and_primers_the_panel_gives(tmp_path, capsys):
    doc, seq = plasmid()
    gb = tmp_path / "t.gb"
    doc.write_genbank(str(gb), None)
    table = tmp_path / "m.csv"
    table.write_text("id,mutation\nsub,gx:p.Ala2Thr\nrev,gy:c.4G>A\ndel,gx:p.Glu3del\nbad,gx:c.4C>A\n", encoding="utf-8")
    out = tmp_path / "out"
    code, text, _ = run(capsys, "kld", "run", gb, table, "--out", out)
    res = json.loads(text)
    assert code == 1 and res["summary"] == {"rows": 3 + 1, "ok": 3, "errors": 1, "with_warnings": res["summary"]["with_warnings"]}
    assert [r["status"] for r in res["rows"]] == ["ok", "ok", "ok", "error"] and "not C" in res["rows"][3]["error"]
    assert sorted(p.name for p in out.iterdir()) == ["del.gb", "primers.csv", "report.jsonl", "rev.gb", "settings.json", "sub.gb"]
    for name, edit in (("sub", (103, 104, "A")), ("rev", (417, 418, "T")), ("del", (106, 109, ""))):
        s, e, t = edit
        made = Document(str(out / f"{name}.gb"))
        assert made.sequence(0, None).upper() == seq[:s] + t + seq[e:]                # the edited plasmid, nothing else changed
        assert {f.label for f in made.features() if f.kind == "primer_bind"} == {"KLD-F", "KLD-R"}               # with the two primers placed on it, as the panel saves it
        direct = kld.design_selected(seq, [edit])                                      # the panel's own function gives the same primers
        row = next(r for r in res["rows"] if r["id"] == name)
        assert [p["sequence"] for p in row["primers"]] == [seq_ for _, seq_, _ in direct.rows]
    rows = list(csv.DictReader(open(out / "primers.csv", encoding="utf-8")))
    assert [r["id"] for r in rows] == ["sub", "rev", "del"] and rows[0]["forward_primer"] == res["rows"][0]["primers"][0]["sequence"]
    assert res["settings"]["kld"]["tm"] == {"value": 60.0, "source": "default"} and "pcr_conditions" in res["settings"]
    # nothing is overwritten, and a bad table is refused whole
    code, text, err = run(capsys, "kld", "run", gb, table, "--out", out)
    assert code == 1 and text == "" and "already exists" in err
    table.write_text("id,mutation,note\nx,gx:A2T,hello\n", encoding="utf-8")
    assert "exactly id,mutation" in run(capsys, "kld", "check", gb, table)[2]
    table.write_text("id,mutation\na b,gx:A2T\nc,gx:A2T\nc,gx:A2T\n", encoding="utf-8")
    assert "used twice" in run(capsys, "kld", "check", gb, table)[2]


# ---- PCR primers ---------------------------------------------------------------------------------------------------------------

import primer_design


def test_pcr_batch_gives_the_pairs_the_panel_gives(tmp_path, capsys):
    doc, seq = plasmid()
    gb = tmp_path / "t.gb"
    doc.write_genbank(str(gb), None)
    table = tmp_path / "t.csv"
    table.write_text("id,target\na,150..400\nb,200..450\nshort,150..161\nbad,nope\n", encoding="utf-8")
    code, text, _ = run(capsys, "pcr", "check", gb, table)   # check writes nothing
    res = json.loads(text)
    assert code == 1 and not (tmp_path / "t_pcr").exists() and res["summary"]["ok"] == 2 and res["summary"]["errors"] == 2
    assert [r["status"] for r in res["rows"]] == ["ok", "ok", "error", "error"]
    assert "No pair found" in res["rows"][2]["error"] and "No feature called" in res["rows"][3]["error"]   # (12 bp: shorter than a primer)
    for row, (a, b) in zip(res["rows"][:2], ((150, 400), (200, 450))):
        direct = primer_design.design_region(doc, a - 1, b)      # the panel's own call
        assert [(p["fwd"]["sequence"], p["rev"]["sequence"], p["product_size"]) for p in row["pairs"]] == \
               [(d["fwd"]["seq"], d["rev"]["seq"], d["size"]) for d in direct]
        best = row["pairs"][0]   # on the ends of the target: the product is the target
        assert best["product_size"] == b - a + 1 and best["fwd"]["start"] == a and best["rev"]["end"] == b and best["mark"] in ("ok", "warn", "bad")
        assert best["fwd"]["sequence"] == seq[a - 1:a - 1 + len(best["fwd"]["sequence"])]
    assert res["settings"]["pcr_design"]["tm"] == {"value": 60.0, "source": "default"} and "pcr_conditions" in res["settings"]
    out = tmp_path / "out"
    code, text, _ = run(capsys, "pcr", "run", gb, table, "--out", out, "--tm", 62)
    res = json.loads(text)
    assert code == 1 and res["settings"]["pcr_design"]["tm"] == {"value": 62.0, "source": "override (command line)"}
    assert sorted(p.name for p in out.iterdir()) == ["primers.csv", "report.jsonl", "settings.json"]
    rows = list(csv.DictReader(open(out / "primers.csv", encoding="utf-8")))
    assert len(rows) == sum(len(r.get("pairs", [])) for r in res["rows"]) and {r["id"] for r in rows} == {"a", "b"}
    assert "already exists" in run(capsys, "pcr", "run", gb, table, "--out", out)[2]
