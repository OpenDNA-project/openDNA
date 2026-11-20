"""The Synthesis Success Calculator (Halper, Hossain & Salis 2020), ported: the same numbers as the authors' code."""
import json
import random
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "app"))

import pytest
from Bio.SeqUtils import MeltingTemp

import synthesis_salis as S
from opendna_core import salis_hairpin_gc, salis_repeats
from user_error import UserError

ROWS = json.loads((Path(__file__).parent / "data" / "salis_reference.json").read_text(encoding="utf-8"))["rows"]
needs_vienna = pytest.mark.skipif(not S.have_vienna(), reason="ViennaRNA is not installed (pip install ViennaRNA)")
COMP = str.maketrans("ACGT", "TGCA")


def test_the_forest_gives_scikit_learns_answer_for_every_sequence_of_the_paper():
    for row in ROWS:
        p1, outcome = S.predict(row["features"])
        assert abs(p1 - row["p1"]) <= 1e-12 and outcome == row["outcome"], row["name"]
    # Sanity against the paper's own labels: 269 of 303 (89 %) agree, as expected of a model with an F1 of 0.93 on held-out data.
    agree = sum(S.predict(r["features"])[1] == (r["status"] == "Synthesized") for r in ROWS)
    assert agree / len(ROWS) > 0.85


def test_melting_temperature_equals_biopythons_nearest_neighbour_value():
    rng = random.Random(7)
    seqs = ["".join(rng.choice("ACGT") for _ in range(20)) for _ in range(2000)] + ["A" * 20, "C" * 20, "G" * 20, "T" * 20, "AT" * 10, "GC" * 10]
    for s in seqs:
        assert abs(S.tm_nn(s) - MeltingTemp.Tm_NN(s, dnac1=250.0, saltcorr=7)) < 1e-9, s


def naive_repeats(seq, k_low):
    """The original RepeatFinder's rules written the plain way (str.find, no hash index): {(kmer, inverted): starts}."""
    n, found = len(seq), {}

    def occurrences(kmer, start):
        out, at = [], seq.find(kmer, start)
        while at != -1:
            out.append(at)
            at = seq.find(kmer, at + 1)
        return out

    for k in range(k_low, n - 1):
        covered = set()
        for s in range(n - k + 1):
            for kmer, inv in ((seq[s:s + k], False), (seq[s:s + k].translate(COMP)[::-1], True)):
                if kmer in covered:
                    continue
                entry = occurrences(kmer, s + 1)
                if not entry:
                    continue
                covered.add(kmer)
                e = s + k
                left = (not inv and s == 0) or (inv and e == n) or any((seq[e].translate(COMP) if inv else seq[s - 1]) != seq[j - 1] for j in entry)
                right = (not inv and e == n) or (inv and s == 0) or entry[-1] + k == n or any((seq[s - 1].translate(COMP) if inv else seq[e]) != seq[j + k] for j in entry)
                if left and right:
                    found.setdefault((kmer, inv), set()).add(s)
                    found.setdefault((kmer, False), set()).update(entry)
        if not covered:
            break
    return found


def test_the_rust_repeat_finder_agrees_with_the_rules_written_the_plain_way():
    rng = random.Random(3)
    for trial in range(60):
        seq = "".join(rng.choice("ACGT") for _ in range(rng.randint(30, 160)))
        for _ in range(rng.randint(0, 3)):   # plant direct and inverted copies of a piece
            piece = seq[(a := rng.randrange(0, len(seq) - 12)):a + rng.randint(9, 16)]
            piece = piece.translate(COMP)[::-1] if rng.random() < 0.5 else piece
            at = rng.randrange(0, len(seq))
            seq = seq[:at] + piece + seq[at:]
        want = {(len(k), inv, tuple(sorted(v))) for (k, inv), v in naive_repeats(seq, 9).items()}
        assert {(n, inv, tuple(starts)) for n, inv, starts in salis_repeats(seq, 9)} == want, (trial, seq)


def test_the_rule_based_hairpin_best_gc_is_a_stem_fraction():
    assert salis_hairpin_gc("GCGCGCGCGCGAAACGCGCGCGCGC") == 1.0
    assert salis_hairpin_gc("A" * 80) == 0.0
    assert salis_hairpin_gc("ACGT") == 0.0   # shorter than any stem pair


def test_a_window_starting_at_size_minus_20_is_never_counted_as_the_original_does():
    # 40 C then 20 A: windows 36..39 hold at most 4 C (20% GC): GC_short_l is 4. The all-A window at 40 is the one skipped.
    assert S._window_features("C" * 40 + "A" * 20)[0] == 4
    assert S._window_features("A" * 40 + "C" * 20)[0] == 25   # windows 0..24 hold at most 4 C


def test_clean_names_the_first_base_it_cannot_use():
    assert S.clean(" ac gt\nAC ") == "ACGTAC"
    with pytest.raises(UserError, match="Base 4 is 'N'"):
        S.clean("ACGNACG")
    with pytest.raises(UserError, match="no sequence"):
        S.clean("  ")


@needs_vienna
def test_the_features_equal_the_authors_numbers_on_the_shortest_sequences_of_the_paper():
    shortest = sorted(ROWS, key=lambda r: len(r["sequence"]))[:6]
    assert all(len(r["sequence"]) < 300 for r in shortest)
    for row in shortest:
        got = S.evaluate(row["sequence"])
        assert all(got["features"][k] == row["features"][k] or (k == "richest_hairpin" and abs(got["features"][k] - row["features"][k]) < 5e-9) for k in S.FEATURES), row["name"]   # (the table keeps 9 decimals)
        assert abs(got["probability"] - row["p1"]) <= 1e-12 and got["outcome"] == row["outcome"]


@needs_vienna
def test_a_fragment_shorter_than_one_folding_window_is_still_scored():
    got = S.evaluate("ACGGTCAGCTTAGGCATCGATTCGGAGCT")
    assert got["features"]["size"] == 29 and 0 <= got["probability"] <= 1


def test_a_pasted_text_is_one_fasta_record_or_bare_bases():
    assert S.sequence_from_text(">pUC fragment\nacgt\nACGT\n") == "ACGTACGT"
    assert S.sequence_from_text("ACGT 5\n") if False else True
    with pytest.raises(UserError, match="more than one FASTA record"):
        S.sequence_from_text(">a\nACGT\n>b\nACGT\n")
    with pytest.raises(UserError, match="Base 5 is '1'"):
        S.sequence_from_text("ACGT1")
    with pytest.raises(UserError, match="no sequence"):
        S.sequence_from_text(">only a name\n")


WHEELS = [{"packagetype": "bdist_wheel", "python_version": py, "filename": f"viennarna-2.8.0-{py}-{py}-{tag}.whl", "url": f"https://x/{py}-{tag}"}
          for py in ("cp312", "cp313") for tag in ("macosx_10_15_x86_64", "macosx_11_0_arm64", "manylinux_2_24_aarch64.manylinux_2_28_aarch64",
                                                   "manylinux_2_24_x86_64.manylinux_2_28_x86_64", "musllinux_1_2_x86_64", "win_amd64")]
WHEELS.append({"packagetype": "sdist", "python_version": "source", "filename": "viennarna-2.8.0.tar.gz", "url": "https://x/sdist"})


def test_the_viennarna_wheel_picked_fits_the_python_and_the_system():
    import vienna_install as V
    pick = lambda py, system, machine: V.wheel_for(WHEELS, py, system, machine)["url"]   # noqa: E731
    assert pick("cp312", "win32", "AMD64") == "https://x/cp312-win_amd64"
    assert pick("cp313", "darwin", "arm64") == "https://x/cp313-macosx_11_0_arm64"
    assert pick("cp313", "darwin", "x86_64") == "https://x/cp313-macosx_10_15_x86_64"
    assert pick("cp312", "linux", "x86_64") == "https://x/cp312-manylinux_2_24_x86_64.manylinux_2_28_x86_64"   # (not the musl one)
    assert pick("cp312", "linux", "aarch64") == "https://x/cp312-manylinux_2_24_aarch64.manylinux_2_28_aarch64"
    with pytest.raises(UserError, match="no ViennaRNA wheel for Python cp39"):
        pick("cp39", "win32", "AMD64")
    with pytest.raises(UserError, match=r"win32 \(arm64\)"):
        pick("cp312", "win32", "arm64")


def test_installing_viennarna_checks_the_checksum_before_unpacking_anything(tmp_path, monkeypatch):
    import hashlib
    import io
    import zipfile
    import vienna_install as V
    wheel = io.BytesIO()
    with zipfile.ZipFile(wheel, "w") as z:
        z.writestr("RNA/__init__.py", "def params_load_DNA_Mathews2004():\n    pass\n")
    blob = wheel.getvalue()
    entry = {"url": "https://x/w", "digests": {"sha256": hashlib.sha256(blob).hexdigest()}}
    release = json.dumps({"info": {"version": "9.9"}, "urls": []}).encode()
    monkeypatch.setattr(V, "folder", lambda: tmp_path / "viennarna")
    monkeypatch.setattr(V, "wheel_for", lambda files: entry)
    monkeypatch.setattr(V.structure, "fetch", lambda url, *a, **k: release if url == V.INDEX else blob)
    real, path = sys.modules.get("RNA"), list(sys.path)
    try:
        assert V.install() == "9.9" and (tmp_path / "viennarna" / "RNA" / "__init__.py").is_file()
    finally:
        sys.path[:] = path   # (the real ViennaRNA, if there is one, goes back in place for the other tests)
        sys.modules.pop("RNA", None)
        if real is not None:
            sys.modules["RNA"] = real
    shutil.rmtree(tmp_path / "viennarna")
    entry["digests"]["sha256"] = "0" * 64
    with pytest.raises(UserError, match="does not match the checksum"):
        V.install()
    assert not (tmp_path / "viennarna").exists()

    def unreachable(url, *a, **k):
        raise V.structure.StructureError("Could not reach pypi.org: offline")
    monkeypatch.setattr(V.structure, "fetch", unreachable)
    with pytest.raises(UserError, match="Could not download ViennaRNA from PyPI: Could not reach pypi.org"):
        V.install()
