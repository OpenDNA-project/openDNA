"""The command line for agents (app/cli.py and the recipes in app/agent/): what an agent reads from a plasmid must be what the
GenBank file says, a mistake must be one `ERROR:` line, and the example in every recipe must still run."""
import csv
import json
import os
import random
import re
import shlex
import sys
import tempfile
from types import SimpleNamespace
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ["OPENDNA_DATA_DIR"] = tempfile.mkdtemp(prefix="opendna_test_")
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "app"))

import pytest

import cas9
import cli
import cli_tools
import enzymes
import kld
import local_align
import msa
import mutation
import primer_design
import recombineering
import repeats
import search
import synthesis_salis
import tools
from cut_sites import NUCLEASES
from opendna_core import Document
from user_error import UserError

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


def test_every_recipe_example_still_runs(capsys, tmp_path, monkeypatch):
    """The `## Example` block of each recipe is a list of commands, one per line: all of them must succeed (a recipe that drifted from
    the tool is caught here)."""
    monkeypatch.chdir(tmp_path)   # (an example that writes a file writes it here)
    fake_ncbi(monkeypatch)        # (and none of them reaches the network)
    code, out, _ = run(capsys, "skills", "list")
    names = [line.split(":")[0] for line in out.splitlines()]
    assert code == 0 and {"index", "inspect", "kld", "pcr", "recombineering", "sites", "guides", "repeats", "search", "align", "msa", "codon", "export", "order", "variants", "checkup", "annotate", "gibson", "ncbi", "synthesis"} <= set(names)
    for name in names:
        if name == "synthesis" and not synthesis_salis.have_vienna():
            continue   # (needs ViennaRNA, which is optional)
        code, text, _ = run(capsys, "skills", "show", name)
        assert code == 0
        for block in re.findall(r"## Example.*?```\n(.*?)```", text, re.S):
            for line in block.splitlines():
                code, _, err = run(capsys, *shlex.split(line.replace("demo/", (ROOT / "demo").as_posix() + "/").replace("app/agent/", (ROOT / "app" / "agent").as_posix() + "/")))
                assert code == 0, f"{name}: {line}: {err}"


# ---- mutations in HGVS and the KLD batch ---------------------------------------------------------------------------------------
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
    # the sender's own design is applied without comment: a lost stop codon is named for what it is, a changed start codon is not warned about
    r = mutation.resolve(doc, "gx:p.Ter7Gln")
    assert (r.effect, r.hgvs_protein, r.warnings) == ("stop lost", "p.(Ter7Glnext*?)", [])
    assert mutation.resolve(doc, "gx:p.Met1Val").warnings == []
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


def test_a_change_through_the_origin_is_made_and_the_file_keeps_its_origin(tmp_path, capsys):
    doc, seq = plasmid()   # 521 bp: bases 520, 521, 1, 2, 3 are neighbours
    gb = tmp_path / "t.gb"
    doc.write_genbank(str(gb), None)
    table = tmp_path / "m.csv"
    table.write_text("id,mutation\ndel,o.520_3del\nins,o.521_1insATG\nrep,o.520_2delinsGG\ndup,o.520_2dup\nlow,o.3_1del\n", encoding="utf-8")
    code, text, _ = run(capsys, "kld", "run", gb, table, "--out", tmp_path / "out")
    res = json.loads(text)
    assert code == 1 and [r["status"] for r in res["rows"]] == ["ok"] * 4 + ["error"] and "o.520_3del" == res["rows"][0]["hgvs_dna"]
    # hand-computed on the string: the origin stays at base 1 unless the change replaced it (then the first base after the change)
    expected = {"del": seq[3:519], "ins": seq + "ATG", "rep": seq[2:519] + "GG", "dup": seq[:2] + seq[519:] + seq[:2] + seq[2:]}
    for name, want in expected.items():
        assert Document(str(tmp_path / "out" / f"{name}.gb")).sequence(0, None).upper() == want, name


def test_kld_batch_gives_the_files_and_primers_the_panel_gives(tmp_path, capsys):
    doc, seq = plasmid()
    gb = tmp_path / "t.gb"
    doc.write_genbank(str(gb), None)
    table = tmp_path / "m.csv"
    table.write_text("id,mutation\nsub,gx:p.Ala2Thr\nrev,gy:c.4G>A\ndel,gx:p.Glu3del\nbad,gx:c.4C>A\n", encoding="utf-8")
    out = tmp_path / "out"
    code, text, _ = run(capsys, "kld", "run", gb, table, "--out", out)
    res = json.loads(text)
    assert code == 1 and res["summary"] == {"rows": 3 + 1, "ok": 3, "errors": 1, "with_warnings": 3, "files_written": 3}
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
    assert res["settings"]["kld"]["tm"] == {"value": 60.0, "source": "default", "min": 45.0, "max": 72.0} and "pcr_conditions" in res["settings"]
    # nothing is overwritten, and a bad table is refused whole
    code, text, err = run(capsys, "kld", "run", gb, table, "--out", out)
    assert code == 1 and text == "" and "already exists" in err
    table.write_text("id,mutation,note\nx,gx:A2T,hello\n", encoding="utf-8")
    assert "exactly id,mutation" in run(capsys, "kld", "check", gb, table)[2]


# ---- PCR primers ---------------------------------------------------------------------------------------------------------------



def test_pcr_batch_gives_the_pairs_the_panel_gives(tmp_path, capsys, monkeypatch):
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
    assert res["settings"]["pcr_design"]["tm"] == {"value": 60.0, "source": "default", "min": 45.0, "max": 75.0} and "pcr_conditions" in res["settings"]
    out = tmp_path / "out"
    code, text, _ = run(capsys, "pcr", "run", gb, table, "--out", out, "--tm", 62)
    res = json.loads(text)
    assert code == 1 and res["settings"]["pcr_design"]["tm"] == {"value": 62.0, "source": "override (command line)", "min": 45.0, "max": 75.0}
    assert sorted(p.name for p in out.iterdir()) == ["primers.csv", "report.jsonl", "settings.json"]
    rows = list(csv.DictReader(open(out / "primers.csv", encoding="utf-8")))
    assert len(rows) == sum(len(r.get("pairs", [])) for r in res["rows"]) and {r["id"] for r in rows} == {"a", "b"}
    assert "already exists" in run(capsys, "pcr", "run", gb, table, "--out", out)[2]
    # by default each primer is a summary (the checks that are not fine); the file and --full keep every check
    first = res["rows"][0]["pairs"][0]["fwd"]
    assert "checks" not in first and first["checks_passed"] + len(first["issues"]) == 8 and all(c["level"] != "ok" for c in first["issues"])
    saved = json.loads((out / "report.jsonl").read_text(encoding="utf-8").splitlines()[0])["pairs"][0]["fwd"]
    full = json.loads(run(capsys, "pcr", "check", gb, table, "--tm", 62, "--full")[1])["rows"][0]["pairs"][0]["fwd"]
    assert len(saved["checks"]) == len(full["checks"]) == 8 and "issues" not in full and [c for c in saved["checks"] if c["level"] != "ok"] == first["issues"]
    # every pair bad: the row warns and offers a wider window, and the window itself stays at 0 (the user's choice)
    import primer_panel
    monkeypatch.setattr(primer_panel, "worst", lambda checks: "bad")
    table.write_text("id,target\na,150..400\n", encoding="utf-8")
    row = json.loads(run(capsys, "pcr", "check", gb, table)[1])
    warns = row["rows"][0]["warnings"]
    assert any("--window 0" in w and "offer" in w for w in warns) and row["settings"]["pcr_design"]["window"]["value"] == 0
    assert any("Only the best 5 pairs are listed" in w for w in warns)   # (a cut list is a warning, never silent)
    # on a circle a window wider than half of what is left would make the two sides meet: it is narrowed, and said
    wide = json.loads(run(capsys, "pcr", "check", gb, table, "--window", 150)[1])["rows"][0]["warnings"]
    assert any("narrowed to 135 bases" in w for w in wide)


# ---- recombineering ------------------------------------------------------------------------------------------------------------



def test_recombineering_batch_gives_the_designs_the_panel_gives(tmp_path, capsys):
    genome, seq = plasmid()
    cas = dna(120, 9)
    genome.write_genbank(str(tmp_path / "g.gb"), None)
    cassette = Document.new_sequence("cas", cas, False)
    cassette.write_genbank(str(tmp_path / "cas.gb"), None)
    table = tmp_path / "r.csv"
    table.write_text("id,region,cassette,orientation\nko,gx,cas.gb,\nflip,150..160,cas.gb,reverse\ndel,gy,none,\nbad,nope,cas.gb,\nori,150..160,cas.gb,up\n", encoding="utf-8")
    out = tmp_path / "out"
    code, text, _ = run(capsys, "recombineering", "run", tmp_path / "g.gb", table, "--out", out)
    res = json.loads(text)
    assert code == 1 and [r["status"] for r in res["rows"]] == ["ok", "ok", "ok", "error", "error"]
    assert "No feature called" in res["rows"][3]["error"] and "forward or reverse" in res["rows"][4]["error"]
    assert sorted(p.name for p in out.iterdir()) == ["del.gb", "flip.gb", "ko.gb", "primers.csv", "report.jsonl", "settings.json"]
    # hand-computed: gx is bases 101..121 (0-based 100..121); the cassette replaces them; 50 bp of the genome either side are the homologies
    ko = res["rows"][0]
    assert Document(str(out / "ko.gb")).sequence(0, None).upper() == seq[:100] + cas + seq[121:]
    f, r = ko["primers"]
    assert f["sequence"].startswith(seq[50:100]) and f["sequence"][50:] == cas[:f["annealing_length"]]
    assert r["sequence"].startswith(mutation.reverse_complement(seq[121:171])) and r["sequence"][50:] == mutation.reverse_complement(cas)[:r["annealing_length"]]
    direct = recombineering.design(genome, cassette, 100, 121)   # the panel's own function gives the same rows
    assert [(p["name"], p["sequence"]) for p in ko["primers"]] == [(n, s) for n, s, _, _ in direct.rows]
    assert ko["warnings"] == ["The region removes the CDS 'gx' (101..121)."]
    assert {x.label for x in Document(str(out / "ko.gb")).features() if x.kind == "primer_bind"} == {"Rec-F", "Rec-R"}
    assert Document(str(out / "flip.gb")).sequence(0, None).upper() == seq[:149] + mutation.reverse_complement(cas) + seq[160:]
    assert Document(str(out / "del.gb")).sequence(0, None).upper() == seq[:400] + seq[421:] and res["rows"][2]["primers"] == []   # a deletion has nothing to amplify
    assert len(list(csv.DictReader(open(out / "primers.csv", encoding="utf-8")))) == 3 and "already exists" in run(capsys, "recombineering", "run", tmp_path / "g.gb", table, "--out", out)[2]
    code, text, _ = run(capsys, "recombineering", "run", tmp_path / "g.gb", table, "--out", tmp_path / "o2", "--primers-only", "--homology", 600)
    assert "no room" in json.loads(text)["rows"][0]["error"] and not list((tmp_path / "o2").glob("*.gb"))


def test_batches_keep_going_past_bad_rows_and_never_collide_on_file_names(tmp_path, capsys, monkeypatch):
    doc, seq = plasmid()
    gb = tmp_path / "t.gb"
    doc.write_genbank(str(gb), None)
    table = tmp_path / "m.csv"
    # cp1252, as Excel saves it: ids that differ only by case, with unsafe characters, a Windows device name, a ragged row, an empty cell,
    # a mutation that changes nothing (base 104 is already G), a failure nobody foresaw
    table.write_bytes("id,mutation\nA,gx:p.Ala2Thr\na,gx:p.Ala2Thr\nbad id!,gx:p.Ala2Thr\nCON,gx:p.Ala2Thr\nragged,gx:p.Ala2Thr,extra\nempty,\n"
                      "noop,o.104G>G\ncaf\xe9,gx:p.Glu3del\nboom,gx:c.7_9dup\n".encode("cp1252"))
    real = mutation.resolve
    monkeypatch.setattr(mutation, "resolve", lambda d, spec, usage=None: (_ for _ in ()).throw(RuntimeError("surprise")) if spec == "gx:c.7_9dup" else real(d, spec, usage))
    out = tmp_path / "out"
    code, text, err = run(capsys, "kld", "run", gb, table, "--out", out)
    res = json.loads(text)
    rows = {r["id"]: r for r in res["rows"]}
    assert code == 1 and [r["status"] for r in res["rows"]] == ["ok", "ok", "ok", "ok", "error", "error", "ok", "ok", "error"]
    files = sorted(p.name for p in out.glob("*.gb"))
    assert files == ["A.gb", "_CON.gb", "a_2.gb", "bad_id_.gb", "caf_.gb", "noop.gb"] and len({f.lower() for f in files}) == len(files)   # every ok row has its own file
    assert all((out / r["file"]).is_file() for r in res["rows"] if r["status"] == "ok") and res["summary"]["files_written"] == 6
    assert "a_2.gb" in rows["a"]["warnings"][-1] and not any("file name" in w for w in rows["A"]["warnings"])   # (only the renamed one is told)
    assert "3 cells but the header has 2" in rows["ragged"]["error"] and "mutation is empty" in rows["empty"]["error"] and "skipped" in rows["empty"]["error"]
    assert rows["boom"]["error"].startswith("Unexpected RuntimeError: surprise") and "Traceback" in rows["boom"]["detail"]
    assert rows["noop"]["effect"] == "no change" and rows["noop"]["primers"] == [] and "no change" in rows["noop"]["warnings"][0] and rows["noop"]["hgvs_dna"] == "o.104="
    assert Document(str(out / "noop.gb")).sequence(0, None).upper() == seq
    assert "ERROR [row ragged]" in err and "ERROR [row boom]" in err and "Traceback" in err    # loud, on stderr too
    assert len(list(csv.DictReader(open(out / "primers.csv", encoding="utf-8")))) == 6
    # a value beyond what the windows accept is used, with a warning that it is the sender's responsibility
    table.write_text("id,mutation\nm,gx:p.Ala2Thr\n", encoding="utf-8")
    code, text, err = run(capsys, "kld", "check", gb, table, "--tm", 80)
    res = json.loads(text)
    assert code == 0 and res["settings"]["kld"]["tm"]["max"] == 72.0 and "responsibility" in res["warnings"][0] and "WARNING: --tm = 80.0" in err
    # the other batch commands keep the same promises: a ragged row or an empty cell is that row's error, and ids differing by case get two files
    ptable = tmp_path / "p.csv"
    ptable.write_text("id,target\nok,150..400\nragged,150..400,x\nempty,\n", encoding="utf-8")
    code, text, err = run(capsys, "pcr", "check", gb, ptable)
    prow = json.loads(text)["rows"]
    assert code == 1 and [r["status"] for r in prow] == ["ok", "error", "error"] and "3 cells but the header has 2" in prow[1]["error"] and "target is empty" in prow[2]["error"]
    assert "ERROR [row ragged]" in err and "ERROR [row empty]" in err
    rtable = tmp_path / "r.csv"
    rtable.write_text("id,region,cassette\nA,gx,none\na,gy,none\nragged,gx,none,extra\nempty,gx,\n", encoding="utf-8")
    code, text, err = run(capsys, "recombineering", "run", gb, rtable, "--out", tmp_path / "rout")
    rrows = json.loads(text)["rows"]
    assert code == 1 and [r["status"] for r in rrows] == ["ok", "ok", "error", "error"] and sorted(p.name for p in (tmp_path / "rout").glob("*.gb")) == ["A.gb", "a_2.gb"]
    assert "cassette is empty" in rrows[3]["error"] and "a_2.gb" in rrows[1]["warnings"][-1]


# ---- the analysis commands (sites, guides, repeats, search): read-only, one JSON answer -------------------------------------------
def planted(tmp_path, *plants, n=400, circular=True):
    """A random n bp sequence with each (0-based position, bases) written over it, saved as a GenBank file: (path, sequence)."""
    seq = list(dna(n, 5))
    for at, bases in plants:
        seq[at:at + len(bases)] = bases
    seq = "".join(seq)
    gb = tmp_path / "p.gb"
    Document.new_sequence("p", seq, circular).write_genbank(str(gb), None)
    return gb, seq


def answer(capsys, *argv):
    code, text, err = run(capsys, *argv)
    assert code == 0, err
    return json.loads(text)


def test_sites_command_gives_the_cuts_and_fragments_of_the_windows(tmp_path, capsys):
    gb, seq = planted(tmp_path, (200, "GAATTC"), (300, "GGATCC"))
    assert seq.count("GAATTC") == seq.count("GGATCC") == 1
    res = answer(capsys, "sites", gb, "--enzymes", "EcoRI,BamHI,NotI", "--digest", "EcoRI,BamHI")
    eco, bam, notl = res["enzymes"]
    # G^AATTC at bases 201..206 cuts the top strand after base 201 and the bottom after 205; GGATCC at 301..306 likewise; NotI is absent
    assert (eco["notation"], eco["sites"], eco["cuts"]) == ("G^AATTC", [{"start": 201, "end": 206, "strand": 1}], [{"top_after_base": 201, "bottom_after_base": 205}])
    assert bam["cuts"][0]["top_after_base"] == 301 and notl["count"] == 0
    assert [(f["size"], f["start"], f["end"], f["left"], f["right"]) for f in res["digest"]["fragments"]] == [(300, 302, 201, ["BamHI"], ["EcoRI"]), (100, 202, 301, ["EcoRI"], ["BamHI"])]
    # the default list is the window's: the enzymes that cut once and recognise 6 bases or more, and it is what the module says
    want = {n for n, c in enzymes.find_cuts(seq, True).items() if len(c) == 1 and len(next(e for e in enzymes.COMMON if e.name == n).site) >= 6}
    assert {e["enzyme"] for e in answer(capsys, "sites", gb)["enzymes"]} == want >= {"EcoRI", "BamHI"}
    assert sorted(f[0] for f in tools.fragments(400, True, {201: ["EcoRI"], 301: ["BamHI"]})) == [100, 300]
    assert "No enzyme called" in run(capsys, "sites", gb, "--enzymes", "Nope")[2]
    both = answer(capsys, "sites", gb, "--digest", "EcoRI,BamHI")    # a digest alone keeps the list of single cutters and adds the fragments
    assert {e["enzyme"] for e in both["enzymes"]} == want and sorted(f["size"] for f in both["digest"]["fragments"]) == [100, 300]


def test_guides_command_gives_the_guides_and_marks_of_the_spacer_finder(tmp_path, capsys):
    guide = "GCGCATGCAAGTCCGATGAC"
    gb, seq = planted(tmp_path, (100, guide + "AGG"))
    res = answer(capsys, "guides", gb, "--no-genome", "--top", "400")
    # the planted site is there: bases 101..123 (guide and PAM) on the top strand, cut 3 bp from the PAM, before base 118
    mine = next(g for g in res["guides"] if g["guide"] == guide)
    assert (mine["start"], mine["end"], mine["strand"], mine["pam"]) == (101, 123, 1, "AGG") and mine["cut"].endswith("before base 118")
    # the window's own three functions give the same list and marks
    nuc = NUCLEASES["SpCas9 (NGG, 20 nt)"]
    found = cas9.candidates(seq, True, nuc, 0, True)
    assert [g["guide"] for g in res["guides"]] == [g.guide for g in found] and res["guides_found"] == len(found)
    specs = cas9.scan(seq, True, "p", found, nuc, 3, None)
    assert [g["mark"] for g in res["guides"]] == [cas9.verdict(s, g)[0] for s, g in zip(specs, found)]
    assert res["settings"]["cas9"]["max_mismatches"] == {"value": 3, "source": "default", "min": 1, "max": 4}
    # a region keeps to its bases, and the limits warn without stopping
    sub = answer(capsys, "guides", gb, "--no-genome", "--region", "90..130", "--max-mm", "9")
    assert all(90 <= g["start"] and g["end"] <= 130 for g in sub["guides"]) and "responsibility" in sub["warnings"][0]
    assert "No nuclease called" in run(capsys, "guides", gb, "--nuclease", "Cas99")[2]
    assert "--top must be at least 1" in run(capsys, "guides", gb, "--no-genome", "--top", "-5")[2]
    assert len(answer(capsys, "guides", gb, "--no-genome", "--top", "2")["guides"]) == 2


def test_repeats_command_gives_the_stretches_the_compare_view_gives(tmp_path, capsys):
    gb, seq = planted(tmp_path, (20, dna(30, 9)), (220, dna(30, 9)))
    res = answer(capsys, "repeats", gb, "--min-len", "20")
    assert [(r["kind"], r["start_in_a"], r["start_in_b"], r["length"], r["differences"]) for r in res["repeats"]] == [("direct", 21, 221, 30, 0)]
    assert [(k, i + 1, j + 1, n) for k, i, j, n, _ in repeats.pair_repeats(seq, seq, 20).repeats] == [("direct", 21, 221, 30)]
    other = tmp_path / "q.gb"
    Document.new_sequence("q", dna(100, 7) + dna(30, 9) + dna(100, 8), False).write_genbank(str(other), None)
    two = answer(capsys, "repeats", gb, other, "--min-len", "25")
    assert {(r["start_in_a"], r["start_in_b"]) for r in two["repeats"] if r["kind"] == "direct"} == {(21, 101), (221, 101)}
    twin = tmp_path / "twin.gb"      # another file, the same bases: they share their whole length (it is only one file with itself that leaves that out)
    Document.new_sequence("twin", seq, True).write_genbank(str(twin), None)
    same = answer(capsys, "repeats", gb, twin, "--min-len", "20")
    first = same["repeats"][0]
    assert same["compared"] == {"a": gb.name, "b": "twin.gb", "same_file": False, "identical_sequences": True} and (first["kind"], first["length"], first["start_in_a"], first["start_in_b"]) == ("direct", len(seq), 1, 1)
    assert "identical" in same["note"] and not answer(capsys, "repeats", gb, "--min-len", "20")["compared"]["identical_sequences"]


@pytest.mark.skipif(not synthesis_salis.have_vienna(), reason="needs ViennaRNA (pip install ViennaRNA)")
def test_synthesis_command_gives_what_the_synthesis_module_gives(tmp_path, capsys):
    gb, seq = planted(tmp_path, n=260)
    res = answer(capsys, "synthesis", gb, "--region", "1..200")
    direct = synthesis_salis.evaluate(seq[:200])
    assert (res["probability"], res["features"], res["predicted_quick_synthesis"]) == (direct["probability"], direct["features"], direct["outcome"])
    assert res["region"] == {"start": 1, "end": 200} and res["bases"] == 200 and "Halper" in res["citation"]
    assert any("straight piece" in w for w in answer(capsys, "synthesis", gb)["warnings"])   # a circular file: the model reads it straight
    through = answer(capsys, "synthesis", gb, "--region", "250..10")   # through the origin of the circle: the piece is joined, read straight
    assert through["bases"] == 21 and through["region"] == {"start": 250, "end": 10, "crosses_origin": True} and through["features"] == synthesis_salis.evaluate(seq[249:] + seq[:10])["features"]
    linear = tmp_path / "linear.gb"
    Document.new_sequence("l", seq, False).write_genbank(str(linear), None)
    assert run(capsys, "synthesis", linear, "--region", "200..100")[2].startswith("ERROR: ")


def test_search_command_finds_dna_on_both_strands_and_names(tmp_path, capsys):
    site = "TTGACAGCTAGCTCAGTCC"
    gb, seq = planted(tmp_path, (50, site))
    rc = mutation.reverse_complement(site)
    both = answer(capsys, "search", gb, rc)["hits"]
    assert [(h["kind"], h["start"], h["end"], h["strand"]) for h in both] == [("dna", 51, 69, -1)]
    assert answer(capsys, "search", gb, rc, "--single-strand")["hits"] == []
    assert [(h["start"], h["strand"]) for h in answer(capsys, "search", gb, site[:-1] + "A", "--fuzzy", "1")["hits"]] == [(51, 1)]
    assert [h["label"] for h in answer(capsys, "search", DEMO, "araC")["hits"]][0] == "araC"
    ed = SimpleNamespace(doc=Document(DEMO), features=Document(DEMO).features(), tracks=[], primer_sites=[], start=0, end=0)
    direct, _ = search.search_sequence(ed, "araC", search.Options(primers=False))
    assert [(h.start + 1, h.end) for h in direct] == [(h["start"], h["end"]) for h in answer(capsys, "search", DEMO, "araC")["hits"]]


def test_every_command_has_a_recipe_example_and_a_test_that_compares_it_with_the_window(capsys):
    """A command added without a recipe whose example runs, or without a `test_<command>_...` here, fails this test."""
    commands = set(cli.parser()._subparsers._group_actions[0].choices) - {"skills", "inspect", "feature", "region"}   # (the three that read a file are tested above, with skills)
    examples = chr(10).join(re.findall(r"## Example.*?```\n(.*?)```", "".join(p.read_text(encoding="utf-8") for p in cli.recipes().values()), re.S))
    mine = [n for n in globals() if n.startswith("test_")]
    for name in commands:
        assert re.search(rf"^{name} ", examples, re.M), f"{name}: no recipe has a runnable example for it"
        assert any(n.startswith(f"test_{name}_") for n in mine), f"{name}: no test_{name}_... compares it with the module the window calls"


def test_align_command_finds_pieces_and_structural_changes_like_local_align(tmp_path, capsys):
    gb, seq = planted(tmp_path)
    one = answer(capsys, "align", gb, seq[100:160])
    assert [(p["query_start"], p["query_end"], p["start"], p["end"], p["strand"], p["identity"]) for p in one["pieces"]] == [(1, 60, 101, 160, 1, 1.0)]
    rev = answer(capsys, "align", gb, mutation.reverse_complement(seq[100:160]))
    assert [(p["start"], p["end"], p["strand"]) for p in rev["pieces"]] == [(101, 160, -1)]
    gap = seq[100:130] + seq[160:190]   # 30 bp of the template skipped between the two pieces
    two = answer(capsys, "align", gb, gap)
    p0, p1 = two["pieces"]   # (the second piece may start a base or two early: the bases next to a junction can match by chance)
    assert (p0["start"], p0["end"], p1["end"]) == (101, 130, 190) and 28 <= p1["start"] - p0["end"] - 1 <= 30 and p1["before"].startswith("deletion from the query")
    base = "A" if seq[130].upper() != "A" else "C"      # one base of the query changed: it is named, with its place, on both strands
    changed = seq[100:130] + base + seq[131:160]
    for q, strand in ((changed, 1), (mutation.reverse_complement(changed), -1)):
        piece = answer(capsys, "align", gb, q)["pieces"][0]
        assert piece["strand"] == strand and [(m["position"], m["plasmid_base"].upper(), m["query_base"].upper()) for m in piece["mismatch_list"]] == [(131, seq[130].upper(), base)]
    cut = answer(capsys, "align", gb, seq[100:130] + seq[135:165])["pieces"][0]     # 5 bases of the plasmid missing from the query (found as one gapped piece or two)
    assert cut["deletion_list"] == [{"start": 131, "end": 135}] or len(answer(capsys, "align", gb, seq[100:130] + seq[135:165])["pieces"]) == 2
    direct, notes = local_align.find_hits(seq, gap, True)    # the window's own function
    assert notes == [] and two["warnings"] == []
    # a changed base just past the origin of the circle is named at its own place (it used to crash the command)
    other = "A" if seq[1].upper() != "A" else "C"
    over = answer(capsys, "align", gb, seq[-30:] + seq[0] + other + seq[2:30])["pieces"][0]
    assert [(m["position"], m["query_base"].upper()) for m in over["mismatch_list"]] == [(2, other)]
    # a query longer than the circle: matched once round, and the answer says so
    assert any("longer than the circular template" in w for w in answer(capsys, "align", gb, seq + seq[:100])["warnings"])
    # a time limit already over: no place aligned, the warning says so and how to search everything; --threads is shown, never saved
    cut = answer(capsys, "align", gb, seq[100:160], "--time-limit", "0.000000001", "--threads", "2")
    assert cut["pieces"] == [] and any(w.startswith("The time limit") for w in cut["warnings"]) and any("--time-limit 0" in w for w in cut["warnings"])
    assert (cut["settings"]["align"]["time_limit"]["value"], cut["settings"]["align"]["cores"]["value"]) == (1e-09, 2)
    assert [(h.start + 1, h.end) for h in direct] == [(p["start"], p["end"]) for p in two["pieces"]]
    err = run(capsys, "align", gb, "MKLVPQ")[2]
    assert err.startswith("ERROR: ") and "protein" in err


def test_msa_command_aligns_the_files_and_never_overwrites(tmp_path, capsys):
    a = dna(200, 11)
    b = a[:80] + a[90:]   # 10 bases deleted
    for name, s in (("a", a), ("b", b)):
        Document.new_sequence(name, s, False).write_genbank(str(tmp_path / f"{name}.gb"), None)
    out = tmp_path / "aln.fa"
    res = answer(capsys, "msa", tmp_path / "a.gb", tmp_path / "b.gb", "--out", out)
    ra, rb = (s["aligned"] for s in res["sequences"])
    assert ra.replace("-", "") == a and rb.replace("-", "") == b and rb.count("-") == 10 and res["width"] == 200
    assert [s["file"] for s in res["sequences"]] == [str(tmp_path / "a.gb"), str(tmp_path / "b.gb")]
    assert [s["identity_to_first"] for s in res["sequences"]] == [100.0, 95.0]    # 190 identical columns of the 200 where either has a base
    assert msa.parse_aligned_fasta(out.read_text())[1] == [ra, rb]
    code, text, err = run(capsys, "msa", tmp_path / "a.gb", tmp_path / "b.gb", "--out", out)
    assert code == 1 and "already exists" in err
    assert "at least two" in run(capsys, "msa", tmp_path / "a.gb")[2]
    # a multi-FASTA is aligned whole (no splitting), the same as its records given one by one
    fa = tmp_path / "both.fasta"
    fa.write_text(f">a\n{a}\n>b\n{b}\n", encoding="utf-8")
    whole = answer(capsys, "msa", fa)
    assert [(s["name"], s["file"], s["aligned"]) for s in whole["sequences"]] == [("a", f"{fa}#a", ra), ("b", f"{fa}#b", rb)]
    assert answer(capsys, "msa", f"{fa}#a", f"{fa}#b")["sequences"][1]["aligned"] == rb


def test_codon_batch_recodes_in_place_keeping_the_protein(tmp_path, capsys):
    doc, seq = plasmid()
    gb = tmp_path / "t.gb"
    doc.write_genbank(str(gb), None)
    table = ROOT / "app" / "agent" / "examples" / "codon_usage_ecoli.csv"
    rows = tmp_path / "r.csv"
    rows.write_text("id,feature\nfwd,gx\nrev,gy\nbad,nope\n", encoding="utf-8")
    out = tmp_path / "out"
    code, text, _ = run(capsys, "codon", "run", gb, rows, "--usage", table, "--method", "Highest frequency", "--out", out)
    res = json.loads(text)
    assert code == 1 and [r["status"] for r in res["rows"]] == ["ok", "ok", "error"] and "No feature called" in res["rows"][2]["error"]
    # hand-computed from the table: the most used synonym of each codon (M ATG, A GCG, E GAA, K AAA, F TTT, G GGC, stop TAA)
    best = "ATGGCGGAAAAATTTGGCTAA"
    assert res["rows"][0]["sequence"] == res["rows"][1]["sequence"] == best and res["rows"][0]["codons_changed"] == 2
    assert Document(str(out / "fwd.gb")).sequence(0, None).upper() == seq[:100] + best + seq[121:]
    assert Document(str(out / "rev.gb")).sequence(0, None).upper() == seq[:400] + mutation.reverse_complement(best) + seq[421:]   # the reverse CDS is written on the other strand
    # the window's own function, with the same seed, gives the same recoding
    import codon_tools
    import codon_usage
    import genetic_code
    s = codon_tools.feature_source(doc, doc.features()[0])
    want = codon_tools.recode(s.seq, genetic_code.codon_map(s.table), codon_usage.read_csv(str(table)), gc_range=codon_tools.DEFAULT_GC, left=s.left, right=s.right, seed=1)
    matched = json.loads(run(capsys, "codon", "check", gb, rows, "--usage", table)[1])
    assert matched["rows"][0]["sequence"] == want.seq and codon_tools.translate(want.seq, genetic_code.codon_map(1)) == codon_tools.translate(CDS, genetic_code.codon_map(1))
    assert "--gc-low is above" in run(capsys, "codon", "check", gb, rows, "--usage", table, "--gc-low", "70")[2]
    one = json.loads(run(capsys, "codon", "check", gb, rows, "--usage", table, "--gc-low", "0", "--gc-high", "20", "--tries", "3")[1])
    assert one["settings"]["codon_usage_table"]["source"] == "command line" and one["settings"]["codon_usage_table"]["value"] == one["settings"]["codon"]["usage_table"]["value"]
    assert one["rows"][0]["seed_used"] in (1, 2, 3) and "above 20 %" in one["rows"][0]["warnings"][0]     # a GC rule that cannot be met says which limit it breaks
    s = codon_tools.feature_source(doc, doc.features()[0])
    best = codon_tools.recode_tries(s, codon_usage.read_csv(str(table)), 4, 1, gc_range=(0, 20))
    assert best.seed in (1, 2, 3, 4) and len(best.unmet) == min(len(codon_tools.recode_source(s, codon_usage.read_csv(str(table)), seed=k, gc_range=(0, 20)).unmet) for k in (1, 2, 3, 4))
    assert "--tries is between" in run(capsys, "codon", "check", gb, rows, "--usage", table, "--tries", "99")[2]
    nogc = json.loads(run(capsys, "codon", "check", gb, rows, "--usage", table, "--no-gc")[1])["settings"]["codon"]   # the settings say what was applied
    assert nogc["gc_limits"]["value"] == "none" and "gc_low" not in nogc and "gc_window" not in nogc


def test_export_command_writes_the_formats_and_never_overwrites(tmp_path, capsys):
    import exports
    doc, seq = plasmid()
    gb = tmp_path / "t.gb"
    doc.write_genbank(str(gb), None)
    ed = SimpleNamespace(doc=doc, features=doc.features(), length=doc.length, circular=True, fasta_record=None)
    for fmt, want in (("fasta", exports.fasta_text("t", doc.sequence(0, None), "[topology=circular]")), ("gff3", exports.gff3_text(ed)),
                      ("bed", exports.bed_text(ed)), ("tbl", exports.tbl_text(ed)), ("proteins", exports.proteins_fasta(ed))):
        out = tmp_path / f"o.{fmt}"
        assert answer(capsys, "export", gb, fmt, out)["file"] == str(out) and out.read_text(encoding="utf-8") == want, fmt
    assert ">gx" in (tmp_path / "o.proteins").read_text() and "MAEKFG" in (tmp_path / "o.proteins").read_text()
    answer(capsys, "export", gb, "genbank", tmp_path / "o.gb")
    assert Document(str(tmp_path / "o.gb")).sequence(0, None).upper() == seq.upper()
    assert "already exists" in run(capsys, "export", gb, "fasta", tmp_path / "o.fasta")[2]


def test_order_command_writes_the_supplier_file_or_nothing(tmp_path, capsys):
    table = tmp_path / "p.csv"
    table.write_text("id,sequence\nF1,atgcgtacgttagcagct\nR1,TTAGCGGCATCGATCGAA\n", encoding="utf-8")
    res = answer(capsys, "order", table, tmp_path / "o.csv", "--scale", "100nm", "--purification", "PAGE")
    assert res["primers"] == 2 and list(csv.reader(open(tmp_path / "o.csv", encoding="utf-8"))) == [
        ["Well Position", "Name", "Sequence", "Scale", "Purification"], ["A1", "F1", "ATGCGTACGTTAGCAGCT", "100nm", "PAGE"], ["B1", "R1", "TTAGCGGCATCGATCGAA", "100nm", "PAGE"]]
    answer(capsys, "order", table, tmp_path / "plain.csv", "--format", "plain")
    assert list(csv.reader(open(tmp_path / "plain.csv", encoding="utf-8"))) == [["F1", "ATGCGTACGTTAGCAGCT"], ["R1", "TTAGCGGCATCGATCGAA"]]
    assert "already exists" in run(capsys, "order", table, tmp_path / "o.csv")[2]
    table.write_text("id,sequence\nF1,ATGCGTACGTTAGCAGCT\nshort,ACGT\nbad,ATGCXTACGTTAG\n", encoding="utf-8")
    code, _, err = run(capsys, "order", table, tmp_path / "never.csv")
    assert code == 1 and "short" in err and "bad" in err and not (tmp_path / "never.csv").exists()


READS = ROOT / "app" / "agent" / "examples" / "pBAD30_eilR_clone_reads.bam"


def test_variants_command_finds_what_the_simulated_reads_were_given(tmp_path, capsys):
    """The reads were simulated from the plasmid with the changes listed in tests/pBAD30_eilR_reads/truth.tsv (the 100 percent ones are consensus changes,
    the 30 percent one is a mixed colony)."""
    res = answer(capsys, "variants", DEMO, READS)
    found = {(v["position"], v["kind"], v["reference"], v["alternative"]): v for v in res["variants"]}
    assert found[(4000, "snp", "C", "A")]["level"] == "consensus" and found[(5300, "snp", "T", "A")]["level"] == "consensus"
    assert found[(4200, "snp", "A", "T")]["level"] == "variant" and 0.2 <= found[(4200, "snp", "A", "T")]["fraction"] < 0.8      # truth.tsv: 30 % of the plasmids
    assert found[(4700, "ins", "", "A")]["level"] == "consensus" and found[(5400, "del", "GGA", "")]["can_apply"]
    assert any(v["kind"] == "gap" and v["position"] == 2000 for v in res["variants"])     # the 100-bp deletion: no read covers it
    # the Variants pane's own scan gives the same list
    import alignments
    from opendna_core import Bam
    direct = alignments.scan_variants(Bam(str(READS)), "pBAD30_eilR", Document(DEMO).sequence(0, None), (5, 20, 80, 13, 0), True)
    assert [(v.pos + 1, v.kind, v.count) for v in direct] == [(v["position"], v["kind"], v["reads"]) for v in res["variants"]]
    # --apply writes the consensus changes into a copy: 1 base inserted, 3 + 1 deleted, the substitution at 4000 made
    out = tmp_path / "fixed.gb"
    applied = answer(capsys, "variants", DEMO, READS, "--apply", out)["applied"]
    fixed = Document(str(out))
    assert fixed.length == 5936 + 1 - 3 - 1 and fixed.sequence(3999, 4000).upper() == "A" and Document(DEMO).sequence(3999, 4000).upper() == "C"
    assert applied["changes"] == sum(1 for v in direct if v.level == "consensus" and v.kind in ("snp", "del", "ins"))
    code, _, err = run(capsys, "variants", DEMO, READS, "--apply", out)
    assert code == 1 and "already exists" in err
    assert "Cannot read" in run(capsys, "variants", DEMO, DEMO)[2] and "responsibility" in answer(capsys, "variants", DEMO, READS, "--min-depth", "0")["warnings"][0]


def test_tool_settings_override_restores_the_outer_value_and_is_never_saved():
    import checkup
    import tool_settings
    key = (checkup.TOOL, "min_promoter")
    with tool_settings.override({key: 41}):
        with tool_settings.override({key: 42}):
            assert tool_settings.value(*key) == 42
        assert tool_settings.value(*key) == 41
    assert key not in tool_settings._overrides and tool_settings.value(*key) != 41


def test_checkup_command_gives_the_verdicts_of_the_checkup_window(tmp_path, capsys):
    doc, seq = plasmid()
    gb = tmp_path / "t.gb"
    doc.write_genbank(str(gb), None)
    bare = answer(capsys, "checkup", gb, "--no-predict", "--no-library")
    assert [(c["cds"], c["verdict"]) for c in bare["cds"]] == [("gx", "Not expressed"), ("gy", "Not expressed")]   # a CDS with nothing upstream
    assert "no bacterial promoter upstream on this strand" in bare["cds"][0]["reasons"]
    # annotated promoter, RBS and terminator around gx (forward, 101..121): the design is complete
    doc.add_features([("promoter", "40..70", [("label", "Pgx")]), ("RBS", "88..99", [("label", "rbs_gx")]), ("terminator", "130..160", [("label", "Tgx")])])
    doc.write_genbank(str(tmp_path / "full.gb"), None)
    full = answer(capsys, "checkup", tmp_path / "full.gb", "--no-predict", "--no-library")
    gx = full["cds"][0]
    assert gx["verdict"] == "Expressed" and gx["hosts"]["Bacteria"]["promoter"] == "Pgx" and gx["hosts"]["Bacteria"]["terminator"] == "Tgx" and full["cds"][1]["verdict"] == "Not expressed"
    assert {(e["kind"], e["origin"]) for e in full["elements"]} == {("promoter", "annotated"), ("RBS", "annotated"), ("terminator", "annotated")}
    # the window's own function gives the same verdicts
    import checkup
    import tool_settings
    with tool_settings.override({(checkup.TOOL, "predict"): False, (checkup.TOOL, "library"): False}):
        direct = checkup.run(cli_tools.doc_editor(Document(str(tmp_path / "full.gb"))))
    assert [c.verdict for c in direct.cds] == [c["verdict"] for c in full["cds"]]
    assert "responsibility" in answer(capsys, "checkup", gb, "--min-promoter", "99")["warnings"][0]
    doc.add_features([("rep_origin", "200..400", [("label", "p15A")])])   # an origin comes out as fields, not as a Python object's text
    doc.write_genbank(str(tmp_path / "ori.gb"), None)
    origin = answer(capsys, "checkup", tmp_path / "ori.gb", "--no-predict", "--no-library")["origins"][0]
    assert (origin["label"], origin["start"], origin["end"], origin["name"], origin["incompatibility_group"]) == ("p15A", 200, 400, "p15A", "p15A")


def test_annotate_command_finds_library_parts_and_adds_only_new_ones(tmp_path, capsys):
    import annotation
    part = max((p for p in annotation.load_parts() if p.kind == "terminator" and set(p.sequence.upper()) <= set("ACGT")), key=lambda p: len(p.sequence))
    gb, seq = planted(tmp_path, (100, part.sequence.upper()), n=600)
    res = answer(capsys, "annotate", gb, "--kinds", "terminator")
    mine = [h for h in res["hits"] if h["start"] == 101]
    assert mine and mine[0]["end"] == 100 + len(part.sequence) and mine[0]["strand"] == 1 and mine[0]["identity"] == 1.0 and not mine[0]["already_annotated"]
    direct = annotation.detect(cli_tools.doc_editor(Document(str(gb))), {"terminator"})
    assert [(h.start + 1, h.end, h.name) for h in direct] == [(h["start"], h["end"], h["name"]) for h in res["hits"]]
    out = tmp_path / "annotated.gb"
    added = answer(capsys, "annotate", gb, "--kinds", "terminator", "--add", out)["added"]
    made = Document(str(out))
    assert added["features"] == len(res["hits"]) and any(f.kind == "terminator" and f.spans[0] == (100, 100 + len(part.sequence)) for f in made.features())
    again = tmp_path / "again.gb"      # a second pass over the annotated copy finds the same parts, already annotated, and adds nothing
    assert answer(capsys, "annotate", out, "--kinds", "terminator", "--add", again)["added"]["features"] == 0
    assert "responsibility" in answer(capsys, "annotate", gb, "--min-dna", "10")["warnings"][0]


@pytest.fixture(autouse=True)
def no_wsl(monkeypatch):
    """No test in this file reaches the real WSL or a real pLannotate (the answer would depend on the machine): every call into WSL fails
    as if it were not installed. A test of the merge puts canned rows in place of the search itself."""
    import plannotate_tool

    def gone(*args, **kwargs):
        raise plannotate_tool.PlannotateError(plannotate_tool.WSL_HELP)
    monkeypatch.setattr(plannotate_tool, "_run", gone)
    plannotate_tool.forget()
    yield
    plannotate_tool.forget()


def test_annotate_command_lists_plannotate_hits_beside_the_library_and_says_who_found_each(tmp_path, capsys, monkeypatch):
    import annotation
    import plannotate_tool
    part = max((p for p in annotation.load_parts() if p.kind == "terminator" and set(p.sequence.upper()) <= set("ACGT")), key=lambda p: len(p.sequence))
    gb, seq = planted(tmp_path, (100, part.sequence.upper()), n=600)
    n = len(part.sequence)
    rows = [{"name": "rrnB T1", "kind": "terminator", "start": 100, "end": 100 + n, "strand": 1, "identity": 0.99, "match": 1.0, "database": "snapgene", "fragment": False, "note": ""},
            {"name": "AmpR", "kind": "CDS", "start": 300, "end": 420, "strand": -1, "identity": 0.9, "match": 0.8, "database": "swissprot", "fragment": True, "note": "beta-lactamase"}]
    monkeypatch.setattr(plannotate_tool, "_search", lambda sequence, circular: rows)
    res = answer(capsys, "annotate", gb, "--kinds", "terminator,CDS")
    by = {h["found_by"]: h for h in res["hits"]}
    assert "doi:10.1093/nar/gkab374" in res["plannotate_citation"]
    assert res["plannotate"] == "ran" and len(res["hits"]) == 2 and set(by) == {"openDNA + pLannotate", "pLannotate"}   # nothing dropped, nothing doubled
    both, theirs = by["openDNA + pLannotate"], by["pLannotate"]
    assert (both["name"], both["start"], both["end"]) == (part.name, 101, 100 + n) and both["plannotate"] == {"name": "rrnB T1", "database": "snapgene", "identity": 0.99}
    assert "pLannotate: snapgene" in both["source"] and part.source in both["source"]
    assert (theirs["name"], theirs["kind"], theirs["start"], theirs["end"], theirs["strand"], theirs["protein_match"], theirs["source"]) == (
        "AmpR (fragment)", "CDS", 301, 420, -1, True, "pLannotate: swissprot")
    off = answer(capsys, "annotate", gb, "--kinds", "terminator,CDS", "--no-plannotate")
    assert off["plannotate"] == "off" and {h["found_by"] for h in off["hits"]} == {"openDNA"} and off["settings"]["detect_features"]["use_plannotate"]["value"] is False
    out = tmp_path / "annotated.gb"
    answer(capsys, "annotate", gb, "--kinds", "terminator,CDS", "--add", out)
    notes = {f.label: dict(f.qualifiers)["note"] for f in Document(str(out)).features() if f.label}
    assert notes["AmpR (fragment)"].startswith("detected by pLannotate (database swissprot): 90 % identity") and "beta-lactamase" in notes["AmpR (fragment)"]
    assert "also found by pLannotate as rrnB T1 (database snapgene, 99 % identity)" in notes[part.name]
    monkeypatch.setattr(plannotate_tool, "_search", lambda sequence, circular: (_ for _ in ()).throw(plannotate_tool.NotInstalled("pLannotate is not installed.")))
    plannotate_tool.forget()
    missing = answer(capsys, "annotate", gb, "--kinds", "terminator,CDS")
    assert missing["plannotate"] == "missing" and "wsltools install" in missing["warnings"][-1] and {h["found_by"] for h in missing["hits"]} == {"openDNA"}


def test_detect_features_merges_only_the_same_thing_in_the_same_place():
    import annotation
    hit = lambda kind, s, e, strand, origin: SimpleNamespace(name=kind, kind=kind, start=s, end=e, strand=strand, origin=origin, also=None, source=origin)
    own = [hit("terminator", 100, 200, 1, "openDNA"), hit("CDS", 300, 600, 1, "openDNA"), hit("promoter", 700, 760, 1, "openDNA")]
    theirs = [hit("terminator", 105, 200, -1, "pLannotate"),   # same bases, other strand: a terminator is annotated on either
              hit("CDS", 300, 600, -1, "pLannotate"),          # same bases but a coding sequence on the other strand: another gene
              hit("misc_feature", 700, 760, 1, "pLannotate"),  # only says misc_feature: the same thing
              hit("terminator", 1000, 1050, 1, "pLannotate")]
    merged = annotation.merge(own, theirs)
    assert [(h.kind, h.start, h.origin) for h in merged] == [("terminator", 100, annotation.BOTH), ("CDS", 300, "openDNA"), ("CDS", 300, "pLannotate"),
                                                              ("promoter", 700, annotation.BOTH), ("terminator", 1000, "pLannotate")]
    assert merged[0].also is theirs[0] and merged[0].source == "openDNA + pLannotate"


def test_plannotate_rows_become_positions_the_way_the_windows_draw_them():
    import plannotate_run
    rows = plannotate_run.rows_from([
        {"name": "AmpR", "type": "CDS", "qstart": 520, "qend": 1381, "sframe": 1, "pident": 100.0, "abs percmatch": 100.0, "db": "snapgene", "fragment": False, "blurb": None},
        {"name": "f1 ori", "type": "origin of replication", "qstart": 5900, "qend": 40, "sframe": -1, "pident": 99.34, "abs percmatch": 96.67, "db": "snapgene", "fragment": True, "blurb": "x"}], 5936)
    assert [(r["start"], r["end"], r["strand"], r["kind"], r["identity"], r["fragment"]) for r in rows] == [
        (520, 1381, 1, "CDS", 1.0, False), (5900, 5976, -1, "rep_origin", 0.9934, True)]   # across the origin: the end runs past the length, like a library hit


def test_wsltools_command_says_what_is_missing_and_never_installs_wsl(capsys, monkeypatch):
    import plannotate_tool
    import runner
    status = answer(capsys, "wsltools", "status")
    assert status["wsl"] is False and status["plannotate"] == {"state": "no_wsl"} and set(status["terminal_tools"]) == set(runner.terminal_tools())
    code, out, err = run(capsys, "wsltools", "install")
    assert code == 1 and err.startswith("ERROR: ") and "wsl --install" in err
    scripts = []   # WSL works: one start answers both, the status line then the installed tools after the mark
    monkeypatch.setattr(plannotate_tool, "_run", lambda script, timeout=None: scripts.append(script) or
                        SimpleNamespace(stdout=f"ready 2.0.0\n{plannotate_tool.TOOLS_MARK}\nfastp\nsamtools\n", returncode=0))
    status = answer(capsys, "wsltools", "status")
    assert len(scripts) == 1 and status["wsl"] is True and status["plannotate"] == {"state": "ready", "version": "2.0.0"}
    assert status["terminal_tools"] == {t: t in ("fastp", "samtools") for t in runner.terminal_tools()}
    both = plannotate_tool.install_script()
    assert "opendna-plannotate" in both and 'plannotate" setupdb' in both and all(t in both for t in runner.terminal_tools()) and both.endswith("exit $r")
    only = plannotate_tool.install_script(plannotate=True, terminal=False)
    assert "fastp" not in only and "opendna-tools" not in only


def test_detect_features_calls_a_part_known_only_when_a_feature_covers_exactly_its_bases():
    import annotation
    s, e = 100, 160
    having = lambda start, end, strand=1: SimpleNamespace(length=1000, visible_features=lambda a, b, overlays=True: [SimpleNamespace(kind="terminator", strand=strand, spans=[(start, end)])])
    assert annotation.already_there(having(s, e), "terminator", s, e, 1) and annotation.already_there(having(s, e, -1), "terminator", s, e, 1)   # (a terminator on either strand)
    assert not annotation.already_there(having(s + 5, e), "terminator", s, e, 1) and not annotation.already_there(having(s - 20, e + 20), "terminator", s, e, 1)   # overlapping is not the same bases
    assert annotation.already_there(having(s + 5, e), "terminator", s, e, 1, exact=False)    # (what the Plasmid Checkup still treats as covered)
    assert not annotation.already_there(having(s, e, -1), "CDS", s, e, 1)    # a CDS: its strand counts


def test_gibson_command_adds_overlaps_and_builds_the_product_the_panel_builds(tmp_path, capsys):
    a, b = dna(300, 21), dna(300, 22)
    for name, s in (("a", a), ("b", b)):
        Document.new_sequence(name, s, False).write_genbank(str(tmp_path / f"{name}.gb"), None)
    table = tmp_path / "f.csv"
    table.write_text("id,source,region,orientation\nfa,a.gb,,\nfb,b.gb,11..290,forward\n", encoding="utf-8")
    out = tmp_path / "out"
    res = answer(capsys, "gibson", "run", table, "--out", out)
    # hand-computed: two linear fragments with nothing in common are joined end to end into a circle (the window starts it at the first
    # fragment's own overlap, a few bases before the fragment: the same circle, turned)
    product = Document(str(out / "assembly.gb"))
    circle = product.sequence(0, None).upper()
    assert res["verdict"] in ("ok", "warn") and product.circular and len(circle) == 580 and (a + b[10:290]) in circle + circle
    prim = {p["name"]: p for p in res["primers"]}
    tail = lambda p: p["sequence"][:len(p["sequence"]) - p["annealing_length"]].upper()   # the 5' part of a primer that does not anneal
    assert a.endswith(tail(prim["fb-F"])) and b[10:290].startswith(mutation.reverse_complement(tail(prim["fa-R"])))   # fb-F carries the end of a; fa-R the start of b, reverse complemented
    assert 20 <= len(tail(prim["fb-F"])) + len(tail(prim["fa-R"])) <= 40 and len(res["junctions"]) == 2
    assert all(p["tail_nt"] == len(tail(p)) for p in res["primers"])    # (what the primer adds at its 5' end)
    # the Assembly panel's own functions give the same primers, checks and product
    import assembly
    from assembly import Frag, snapshot
    frags = [Frag("fa", snapshot(Document(str(tmp_path / "a.gb")))), Frag("fb", snapshot(Document(str(tmp_path / "b.gb")), "fb"), 10, 290)]
    assembly.auto_primers(frags, 0, True)
    d = assembly.analyse(frags, 0, True)
    assert [p["sequence"] for p in res["primers"]] == [p[1] for f in frags for p in (f.fwd, f.rev)]
    assert assembly.simulate(frags, d).sequence(0, None).upper() == product.sequence(0, None).upper() and [c[1] for c in d.checks] == [c["level"] for c in res["checks"]]
    # an assembly needs every fragment: a bad row stops the command, and the message names it
    table.write_text("id,source,region,orientation\nfa,a.gb,,\nfb,b.gb,0..5,\n", encoding="utf-8")
    code, _, err = run(capsys, "gibson", "check", table)
    assert code == 1 and err.startswith("ERROR: ") and "positions must satisfy" in err
    table.write_text("id,source\nfa,a.gb\nfa,b.gb\n", encoding="utf-8")
    assert "must be different" in run(capsys, "gibson", "check", table)[2]
    table.write_text("id,source\nfa,a.gb\nfb,b.gb\n", encoding="utf-8")
    assert "responsibility" in answer(capsys, "gibson", "check", table, "--min-overlap", "5")["warnings"][0]
    assert "already exists" in run(capsys, "gibson", "run", table, "--out", out)[2]


def test_ncbi_command_searches_and_fetches_through_the_import_dialogs_functions(tmp_path, capsys, monkeypatch):
    seen = fake_ncbi(monkeypatch)
    res = answer(capsys, "ncbi", "search", "pUC19", "--limit", "3")
    assert [(r["accession"], r["length"], r["topology"]) for r in res["records"]] == [("L09137.2", 60, "circular")] and "pUC19" in res["settings"]["ncbi"]["sends"]["value"]
    assert "retmax=3" in seen[0] and all("eutils.ncbi.nlm.nih.gov" in u for u in seen) and not any("api_key" in u for u in seen)    # only the words are sent
    out = tmp_path / "puc.gb"
    got = answer(capsys, "ncbi", "get", "L09137", "--out", out)    # an accession without its version finds the current one
    assert got["record"]["accession"] == "L09137.2" and Document(str(out)).length == 60 and out.read_text().startswith("LOCUS")
    import ncbi
    direct = ncbi.search("nuccore", "pUC19", "", 3)
    assert [(r.accession, r.length) for r in direct] == [("L09137.2", 60)]
    for argv, words in ((("get", "L09137", "--out", out), "already exists"), (("get", "L09137"), "needs --out"), (("get", "XYZ", "--out", tmp_path / "x.gb"), "is not the accession"),
                        (("search", "  "), "Type an accession")):
        code, text, err = run(capsys, "ncbi", *argv)
        assert code == 1 and text == "" and words in err and err.startswith("ERROR: ")
    asked = len(seen)
    for argv, words in ((("search", "pUC19", "--limit", "0"), "--limit must be at least 1"), (("get", "L09137", "--out", out), "already exists")):
        assert words in run(capsys, "ncbi", *argv)[2]
    assert len(seen) == asked    # (a mistake that can be seen here never reaches the network)


def fake_ncbi(monkeypatch):
    """Replace the one network call of ncbi.py with a canned NCBI: a search finds pUC19-like record L09137.2 (60 bp); the requests made are listed."""
    import ncbi
    gb = Path(tempfile.mkdtemp()) / "r.gb"
    Document.new_sequence("L09137", dna(60, 3), True).write_genbank(str(gb), None)
    asked = []

    def fetch(url, timeout=60):
        asked.append(url)
        if "esearch" in url:
            return json.dumps({"esearchresult": {"idlist": ["111"]}})
        if "esummary" in url:
            return json.dumps({"result": {"111": {"accessionversion": "L09137.2", "title": "Cloning vector pUC19", "slen": 60, "organism": "synthetic construct",
                                                 "updatedate": "2026/01/01", "topology": "circular", "sourcedb": "genbank"}}})
        return gb.read_text(encoding="utf-8")
    monkeypatch.setattr(ncbi, "fetch", fetch)
    monkeypatch.setattr(ncbi.time, "sleep", lambda s: None)
    monkeypatch.setattr(ncbi, "path_of", lambda record: Path(tempfile.mkdtemp()) / "cache.gb")   # (never a file left by an earlier test)
    return asked


def test_split_command_writes_one_file_per_record_and_file_hash_name_reads_one(tmp_path, capsys):
    a, b = dna(120, 21), dna(90, 22)
    for name, seq, circ in (("p one", a, True), ("p/two", b, False)):
        Document.new_sequence(name.replace(" ", "_").replace("/", "_"), seq, circ).write_genbank(str(tmp_path / f"{name[-3:]}.gb"), None)
    import multigb
    multi = tmp_path / "multi.gb"
    multi.write_text("".join((tmp_path / f"{n}.gb").read_text(encoding="utf-8").rstrip("\n") + "\n" for n in ("one", "two")), encoding="utf-8")
    assert [n for n, *_ in multigb.records(multi)] == ["p_one", "p_two"]
    err = run(capsys, "inspect", multi)[2]   # a bare name is refused with the names, and the way out
    assert err.startswith("ERROR: ") and "p_one, p_two" in err and "split" in err
    assert answer(capsys, "region", f"{multi}#p_two", 1, 10)["sequence"] == b[:10]
    assert "no record called" in run(capsys, "region", f"{multi}#nope", 1, 10)[2]
    res = answer(capsys, "split", multi, "--out", tmp_path / "parts")
    assert [(f["record"], f["file"], f["length"], f["circular"]) for f in res["files"]] == [("p_one", "p_one.gb", 120, True), ("p_two", "p_two.gb", 90, False)]
    assert Document(str(tmp_path / "parts" / "p_two.gb")).sequence(0, 90) == b
    assert "not empty" in run(capsys, "split", multi, "--out", tmp_path / "parts")[2] and "one record" in run(capsys, "split", tmp_path / "one.gb")[2]
    # a multi-FASTA gives FASTA files; two records with one name are told apart, with a warning
    fa = tmp_path / "m.fasta"
    fa.write_text(f">x\n{a}\n>X\n{b}\n", encoding="utf-8")
    res = answer(capsys, "split", fa, "--out", tmp_path / "fparts")
    assert [f["file"] for f in res["files"]] == ["x.fasta", "X_2.fasta"] and res["warnings"] and Document.from_fasta(str(tmp_path / "fparts" / "X_2.fasta"), "X").sequence(0, 90) == b


def test_merge_command_joins_files_and_records_and_keeps_what_the_format_keeps(tmp_path, capsys):
    a, b = dna(120, 31), dna(90, 32)
    Document.new_sequence("pa", a, True).write_genbank(str(tmp_path / "a.gb"), None)
    Document.new_sequence("pb", b, False).write_genbank(str(tmp_path / "b.gb"), None)
    out = tmp_path / "all.gb"
    res = answer(capsys, "merge", tmp_path / "a.gb", tmp_path / "b.gb", tmp_path / "a.gb", "--out", out)
    assert [(r["name"], r["length"], r["circular"]) for r in res["records"]] == [("pa", 120, True), ("pb", 90, False), ("a", 120, True)] and res["warnings"]
    assert answer(capsys, "region", f"{out}#pb", 1, 90)["sequence"] == b and answer(capsys, "region", f"{out}#a", 1, 120)["sequence"] == a   # round trip
    assert "already exists" in run(capsys, "merge", tmp_path / "a.gb", tmp_path / "b.gb", "--out", out)[2]
    # two records with one name: the second takes its file's name; --names sets them all (one each)
    Document.new_sequence("pa", b, False).write_genbank(str(tmp_path / "other.gb"), None)
    clash = answer(capsys, "merge", tmp_path / "a.gb", tmp_path / "other.gb", "--out", tmp_path / "clash.gb")
    assert [r["name"] for r in clash["records"]] == ["pa", "other"] and clash["warnings"]
    named = answer(capsys, "merge", tmp_path / "a.gb", tmp_path / "other.gb", "--names", "x,y", "--out", tmp_path / "named.gb")
    assert [r["name"] for r in named["records"]] == ["x", "y"] and named["warnings"] == []
    assert "one name each" in run(capsys, "merge", tmp_path / "a.gb", tmp_path / "other.gb", "--names", "x", "--out", tmp_path / "bad.gb")[2]
    assert "two or more" in run(capsys, "merge", tmp_path / "a.gb", "--out", tmp_path / "x.gb")[2]
    # a record picked with FILE#NAME; a .fasta name writes only the bases and says so
    fa = answer(capsys, "merge", f"{out}#pb", tmp_path / "a.gb", "--out", tmp_path / "two.fasta")
    assert [r["name"] for r in fa["records"]] == ["pb", "pa"] and any("FASTA keeps only the bases" in w for w in fa["warnings"])
    assert (tmp_path / "two.fasta").read_text(encoding="utf-8").split(">")[1:] and answer(capsys, "region", f"{tmp_path / 'two.fasta'}#pa", 1, 5)["sequence"] == a[:5]


def test_circular_tools_give_the_same_answers_on_a_plasmid_turned_round_the_origin(tmp_path, capsys):
    """A circle turned by k bases is the same molecule. With a restriction site, a guide, a repeat and a CDS each made to cross the origin in
    turn, sites, guides, repeats, DNA and protein search, feature, region, pcr and kld find the same things as on the unturned plasmid."""
    guide = "GACTCGTACGTAGCATCGAT"
    rep, cds = dna(30, 41), "ATGGCTAAAGACGCTCATCATCACTTTGAAGCGGGTAATTGGCATTTAACCTAGGGTTAA"[:57] + "TAA"
    n = 900
    gb, seq = planted(tmp_path, (299, "GAATTC"), (99, guide + "TGG"), (49, rep), (599, rep), (699, cds), n=n)
    doc = Document(str(gb))
    doc.add_features([("CDS", "700..759", [("label", "wrapCDS")])])
    doc.write_genbank(str(gb), None)
    pcr_table, kld_table, codon_table = tmp_path / "p.csv", tmp_path / "k.csv", tmp_path / "c.csv"
    codon_table.write_text("id,feature" + chr(10) + "c,wrapCDS" + chr(10), encoding="utf-8")
    usage = ROOT / "app" / "agent" / "examples" / "codon_usage_ecoli.csv"
    pcr_table.write_text("id,target" + chr(10) + "w,wrapCDS" + chr(10), encoding="utf-8")
    kld_table.write_text("id,mutation" + chr(10) + "m,wrapCDS:p.Ala2Gly" + chr(10), encoding="utf-8")

    def signature(path, k):
        back = lambda p: (p - 1 + k) % n + 1   # (a position on the turned plasmid, on the original)
        sites = sorted(back(h["start"]) for e in answer(capsys, "sites", path, "--enzymes", "EcoRI")["enzymes"] for h in e["sites"])
        guides = [(g["guide"], g["strand"]) for g in answer(capsys, "guides", path, "--no-genome", "--top", 5000)["guides"] if g["guide"] == guide]
        reps = [(r["kind"], r["length"]) for r in answer(capsys, "repeats", path, "--min-len", 25)["repeats"]]
        found = sorted(back(h["start"]) for h in answer(capsys, "search", path, "GAATTC")["hits"])
        prot = [back(h["start"]) for h in answer(capsys, "search", path, "MAKDAHHHF", "--protein")["hits"] if h["kind"] == "protein"]
        feat = answer(capsys, "feature", path, "wrapCDS")
        row = json.loads(run(capsys, "pcr", "check", path, pcr_table)[1])["rows"][0]["pairs"][0]
        kld = json.loads(run(capsys, "kld", "check", path, kld_table)[1])["rows"][0]
        codon = json.loads(run(capsys, "codon", "check", path, codon_table, "--usage", usage)[1])["rows"][0]
        return (sites, guides, reps, found, prot, feat["length"], feat["protein"]["sequence"], row["product_size"], row["fwd"]["sequence"], kld["hgvs_protein"], kld["effect"],
                codon["codons_changed"], codon["sequence"])

    base = signature(gb, 0)
    assert 300 in base[0] and base[1] == [(guide, 1)] and [r[0] for r in base[2]] == ["direct"] and 30 <= base[2][0][1] <= 32 and base[4] == [700] and base[9] == "p.(Ala2Gly)"
    for k in (302, 110, 64, 614, 720):   # the site, the guide, each copy of the repeat, the CDS: the origin falls inside it
        turned = Document(str(gb))
        turned.rotate(k)
        path = tmp_path / f"turned{k}.gb"
        turned.write_genbank(str(path), None)
        assert signature(path, k) == base, k
    # a region through the origin, in both commands that take one
    turned = tmp_path / "turned720.gb"
    start = (700 - 1 - 720) % n + 1
    assert answer(capsys, "region", turned, start, start + 59 - n)["sequence"] == cds
    assert answer(capsys, "inspect", turned)["features"][0]["crosses_origin"] is True
    wrapped = [c for c in answer(capsys, "checkup", turned, "--no-predict", "--no-library")["cds"] if c["cds"] == "wrapCDS"][0]   # (the checkup says where a CDS across the origin is)
    assert (wrapped["start"], wrapped["end"], wrapped["crosses_origin"]) == (start, start + 59 - n, True)
    # a residue that is split by the origin cannot be one block of bases: said, not turned into a change over the whole plasmid
    split = Document(str(gb))
    split.rotate(721)
    split.write_genbank(str(tmp_path / "split.gb"), None)
    kld_table.write_text("id,mutation" + chr(10) + "m,wrapCDS:p.His8Ala" + chr(10), encoding="utf-8")
    assert "cross the origin" in json.loads(run(capsys, "kld", "check", tmp_path / "split.gb", kld_table)[1])["rows"][0]["error"]
    # the recoded plasmid keeps the protein, the gene's place and every base outside it
    answer(capsys, "codon", "run", turned, codon_table, "--usage", usage, "--out", tmp_path / "recoded")
    before, after = Document(str(turned)), Document(str(tmp_path / "recoded" / "c.gb"))
    cds_places = set(range(880, 901)) | set(range(1, 40))
    changed = [i + 1 for i in range(n) if before.sequence(i, i + 1) != after.sequence(i, i + 1)]
    assert changed and set(changed) <= cds_places
    assert answer(capsys, "feature", tmp_path / "recoded" / "c.gb", "wrapCDS")["protein"]["sequence"] == answer(capsys, "feature", turned, "wrapCDS")["protein"]["sequence"]


def test_a_limit_that_is_reached_is_always_a_warning(tmp_path, capsys, monkeypatch):
    """A list cut short, a search stopped or a query not read in full is never silent: each says so in `warnings`, in every command."""
    import ncbi
    import search
    gb, seq = planted(tmp_path, (20, "GAATTC"), (60, "GAATTC"), (120, dna(30, 9)), (220, dna(30, 9)), (160, dna(30, 10)), (260, dna(30, 10)))
    assert any("Only the best 1 of the" in w for w in answer(capsys, "guides", gb, "--no-genome", "--top", 1)["warnings"])
    (tmp_path / "n").mkdir()
    holes, _ = planted(tmp_path / "n", (150, "NNNN"), n=400)   # (bases that are N: no guide over them)
    assert any("4 base(s) are N" in w for w in answer(capsys, "guides", holes, "--no-genome")["warnings"])
    monkeypatch.setattr(search, "MAX_HITS", 1)
    assert any("stopped at" in w for w in answer(capsys, "search", gb, "GAATTC")["warnings"])
    monkeypatch.setattr(repeats, "MAX_REPEATS", 1)
    assert any("only the longest 1 of" in w for w in answer(capsys, "repeats", gb, "--min-len", "20")["warnings"])
    monkeypatch.setattr(local_align, "MAX_HITS", 1)
    assert any("Only 1 pieces are listed" in w for w in answer(capsys, "align", gb, seq[10:40] + seq[100:130])["warnings"])
    asked = fake_ncbi(monkeypatch)
    plain = ncbi.fetch
    monkeypatch.setattr(ncbi, "fetch", lambda url, timeout=60: plain(url).replace('"idlist"', '"count": "40", "idlist"') if "esearch" in url else plain(url))
    words = answer(capsys, "ncbi", "search", "pUC19", "--limit", 1)["warnings"]
    assert any("40 matches; only the first 1 are shown" in w for w in words) and asked


def test_variants_command_sees_changes_on_both_sides_of_the_origin_and_apply_writes_them(tmp_path, capsys):
    """Reads of a circular plasmid that run over the origin are mapped as two pieces (one clipped at each end of the reference): SNPs at bases 1, 5 and
    the last base are found with their full depth and written by --apply, the clips at the ends are not reported as breaks, and a deletion across the origin
    is called as one deletion through the origin."""
    from test_core_formats import read, write_bam
    n = 600
    gb, ref = planted(tmp_path, n=n)
    other = lambda p: "A" if ref[p - 1] != "A" else "C"

    def reads_of(sample, first=0):
        """Reads of 100 bases every 10 along the circle `sample` (whose base 0 is reference base `first`); one over its end is two pieces."""
        m, out = len(sample), []
        for i, s in enumerate(range(0, len(sample), 10)):
            seq, a = (sample + sample)[s:s + 100], min(100, m - s)
            if s + 100 <= m:
                out.append(read(f"r{i}", first + s, [(100, "M")], seq))
            else:
                out += [read(f"r{i}", first + s, [(a, "M"), (100 - a, "S")], seq), read(f"r{i}", first, [(a, "S"), (100 - a, "M")], seq, flag=2048)]
        return out

    snps = {1: other(1), 5: other(5), n: other(n)}
    sample = "".join(snps.get(p, c) for p, c in enumerate(ref, 1))
    bam = tmp_path / "s.bam"
    write_bam(bam, [("p", n)], reads_of(sample))
    fixed = tmp_path / "fixed.gb"
    res = answer(capsys, "variants", gb, bam, "--min-depth", 3, "--apply", fixed)
    assert [(v["position"], v["kind"], v["alternative"], v["reads"], v["depth"], v["level"]) for v in res["variants"]] == \
           [(p, "snp", snps[p], 10, 10, "consensus") for p in (1, 5, n)]
    corrected = Document(str(fixed))
    assert corrected.sequence(0, None) == sample and corrected.circular
    # 80 bases deleted across the origin (reference bases 561..600 and 1..40): one deletion, not two breaks and a gap, and --apply writes it
    deleted = ref[40:560]
    write_bam(bam.with_name("d.bam"), [("p", n)], reads_of(deleted, first=40))
    cut = tmp_path / "cut.gb"
    found = answer(capsys, "variants", gb, bam.with_name("d.bam"), "--min-depth", 3, "--apply", cut)["variants"]
    assert [(v["position"], v["kind"], len(v["reference"]), v["level"]) for v in found] == [(561, "del", 80, "consensus")] and "through the origin" in found[0]["change"]
    assert Document(str(cut)).sequence(0, None) == deleted
    # reads hard-clipped at the origin hold no bases to compare with the circle: reported at the start, with that said (never dropped)
    hard = [r if not r["flag"] else read(r["name"], 0, [(r["cigar"][0][0], "H"), r["cigar"][1]], r["seq"][r["cigar"][0][0]:], flag=2048) for r in reads_of(ref)]
    write_bam(bam.with_name("h.bam"), [("p", n)], hard)
    notes = answer(capsys, "variants", gb, bam.with_name("h.bam"), "--min-depth", 3)["variants"]
    assert [(v["position"], v["kind"]) for v in notes] == [(1, "clip")] and "could not be checked" in notes[0]["change"]
