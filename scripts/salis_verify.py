"""Check, on every sequence of the paper, that openDNA's Synthesis Success Calculator gives exactly the original's numbers.

    .venv/Scripts/python.exe scripts/salis_verify.py [--workers 2] [--only NAME ...]

Reads tests/data/salis_reference.json (303 sequences of Halper, Hossain & Salis 2020 with the nine features of the
authors' own table, the repeat metric from their code, and the forest's answer from scikit-learn) and recomputes all of
it. About 80 CPU-minutes for the whole file (folding is the slow part), so it is a script, not a test; the test suite
runs the shortest sequences. It exits with 1 if any number differs. Needs ViennaRNA. Keep --workers small: each one
uses a full core.
"""
import argparse
import json
import multiprocessing as mp
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "app"))
import synthesis_salis as S  # noqa: E402


def check(row):
    got = S.evaluate(row["sequence"])
    # (richest_hairpin: the authors' table keeps 9 decimals; every other number must be identical)
    bad = [f"{k}: got {got['features'][k]!r}, original {row['features'][k]!r}" for k in S.FEATURES
           if got["features"][k] != row["features"][k] and not (k == "richest_hairpin" and abs(got["features"][k] - row["features"][k]) < 5e-9)]
    if abs(got["probability"] - row["p1"]) > 1e-12 or got["outcome"] != row["outcome"]:
        bad.append(f"probability: got {got['probability']!r}, scikit-learn {row['p1']!r}")
    return row["name"], len(row["sequence"]), bad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--only", nargs="*", help="names of the sequences to check (default: all)")
    args = ap.parse_args()
    rows = json.load(open(ROOT / "tests" / "data" / "salis_reference.json", encoding="utf-8"))["rows"]
    if args.only:
        rows = [r for r in rows if r["name"] in args.only]
    rows.sort(key=lambda r: -len(r["sequence"]))   # (longest first: the workers finish together)
    failed = 0
    with mp.Pool(max(1, args.workers)) as pool:
        for done, (name, size, bad) in enumerate(pool.imap_unordered(check, rows), 1):
            failed += bool(bad)
            print(f"{done}/{len(rows)} {name} ({size} bp): {'DIFFERS' if bad else 'same'}", flush=True)
            for line in bad:
                print("    " + line, flush=True)
    print(f"{len(rows) - failed} of {len(rows)} identical to the original" + ("" if not failed else f"; {failed} differ"))
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
