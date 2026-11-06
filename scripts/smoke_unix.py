"""Smoke test of a macOS or Linux build: start the program the way a user would (with a file on the command line), check that it
stays alive, shows the file in its window, finds every data file it bundles, and has printed no Python error. The macOS / Linux
counterpart of scripts/smoke_exe.ps1, run by scripts/build_unix.sh (and so by the release workflow).

  python scripts/smoke_unix.py dist/openDNA/openDNA [--file demo/demo_vector.gb] [--seconds 8]
  python scripts/smoke_unix.py dist/openDNA.app/Contents/MacOS/openDNA

It runs headless (QT_QPA_PLATFORM=offscreen, a fresh data folder): no window appears on the build machine. The program writes what
it shows into a report file (OPENDNA_SMOKE_REPORT, see main.write_smoke_report). Exit code 0 = passed.
"""
import argparse
import json
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("program", type=Path)
    parser.add_argument("--file", type=Path, default=ROOT / "demo" / "demo_vector.gb")
    parser.add_argument("--seconds", type=float, default=8)
    args = parser.parse_args()
    program, file = args.program.resolve(), args.file.resolve()
    if not program.is_file() or not os.access(program, os.X_OK):
        sys.exit(f"FAIL: {program} is not an executable file")
    with tempfile.TemporaryDirectory(prefix="opendna_smoke_") as tmp:
        report = Path(tmp) / "report.json"
        env = {k: v for k, v in os.environ.items() if not k.startswith(("PYTHON", "QT_", "OPENDNA_"))}   # nothing from the build
        env.update(QT_QPA_PLATFORM="offscreen", OPENDNA_DATA_DIR=str(Path(tmp) / "data"), OPENDNA_SMOKE_REPORT=str(report),
                   HOME=str(Path(tmp) / "home"))
        Path(env["HOME"]).mkdir()
        log = open(Path(tmp) / "output.txt", "w+b")
        started = time.time()
        proc = subprocess.Popen([str(program), str(file)], env=env, cwd=tmp, stdout=log, stderr=subprocess.STDOUT)
        try:
            while time.time() - started < args.seconds and proc.poll() is None:
                time.sleep(0.25)
            alive = proc.poll() is None
        finally:
            if proc.poll() is None:
                proc.terminate()
                try:
                    proc.wait(10)
                except subprocess.TimeoutExpired:
                    proc.kill()
        log.seek(0)
        output = log.read().decode("utf-8", "replace")
        problems = []
        if not alive:
            problems.append(f"the program ended after {time.time() - started:.1f} s (exit code {proc.returncode})")
        if "Traceback (most recent call last)" in output:
            problems.append("it printed a Python error")
        shown = json.loads(report.read_text(encoding="utf-8")) if report.exists() else None
        if shown is None:
            problems.append("it wrote no report (it did not get as far as showing its window)")
        else:
            if not shown["frozen"]:
                problems.append("this is not a built program (sys.frozen is not set)")
            if not shown["file"] or Path(shown["file"]).resolve() != file:
                problems.append(f"the window does not show {file.name} (it shows {shown['file']})")
            if not shown["length"]:
                problems.append("the sequence is empty")
            if shown["missing"]:
                problems.append("bundled files not found: " + ", ".join(shown["missing"]))
        print(f"program  {program}")
        if shown:
            print(f"window   {shown['title']!r}: {shown['sequence']}, {shown['length']:,} bp (Qt {shown['qt']}, {shown['platform']})")
        if problems:
            print("output:\n" + output[-4000:])
            sys.exit("FAIL: " + "; ".join(problems))
        print("smoke test passed")


if __name__ == "__main__":
    main()
