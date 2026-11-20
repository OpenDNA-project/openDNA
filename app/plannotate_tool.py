"""pLannotate, the plasmid annotator of the Barrick lab, used beside openDNA's own Features Database. It is a separate tool, not
part of openDNA: it is installed by itself (conda, from bioconda, with its own databases) in its own environment, and openDNA calls
pLannotate's own function in it (plannotate_run.py), so what comes back is exactly what pLannotate finds.

    McGuffie MJ, Barrick JE. pLannotate: engineered plasmid annotation. Nucleic Acids Res 49(W1), W516-W522 (2021).
    doi:10.1093/nar/gkab374, https://github.com/barricklab/pLannotate (GPL-3.0)

Where it runs is the same as the Terminal's tools (runner.py): inside WSL on Windows, bash elsewhere, in the conda environment
`opendna-plannotate`. Installing WSL itself is left to the user (it needs an administrator and a restart).
"""
import hashlib
import json
import os
import shlex
import subprocess
import sys
import tempfile
from pathlib import Path
from types import SimpleNamespace

import runner

ENV = "opendna-plannotate"
NAME = "pLannotate"
CITATION = "McGuffie MJ, Barrick JE. pLannotate: engineered plasmid annotation. Nucleic Acids Res 49(W1), W516-W522 (2021)."
DOI = "10.1093/nar/gkab374"
MAX_LENGTH = 50_000    # pLannotate's own limit ("are you sure this is an engineered plasmid?"); a longer sequence is not sent
WSL_HELP = "WSL is not available. Install it once, in an administrator terminal: wsl --install (then restart), and try again."
HELPER = Path(__file__).with_name("plannotate_run.py")
_cache = {}    # (sequence hash, circular) -> rows: the window asks again at every change of its settings
_missing = []  # [True] once the search found pLannotate not installed


class PlannotateError(Exception):
    pass


class NotInstalled(PlannotateError):
    pass


def _run(script, timeout=None, **kwargs):
    """Run bash text (in WSL on Windows). Returns the CompletedProcess; no console window flashes up behind the application."""
    flags = getattr(subprocess, "CREATE_NO_WINDOW", 0) if os.name == "nt" else 0
    try:
        return subprocess.run(runner.command_line(["-lc", script]), capture_output=True, text=True, encoding="utf-8", errors="replace",
                              timeout=timeout, creationflags=flags, **kwargs)
    except OSError as e:   # no wsl.exe
        raise PlannotateError(f"{WSL_HELP} ({e})") from None
    except subprocess.TimeoutExpired as e:   # it hangs
        raise PlannotateError(f"Timed out: {e}") from None


def env_lookup():
    """bash that sets E to the folder of the pLannotate environment and puts its programs first on PATH (blastn, diamond and
    cmscan are found by name); E stays empty when there is none."""
    bases = " ".join(runner.CONDA_BASES)
    return (f'E=; for b in {bases}; do [ -x "$HOME/$b/envs/{ENV}/bin/python" ] && E="$HOME/$b/envs/{ENV}" && break; done; '
            '[ -n "$E" ] && export PATH="$E/bin:$PATH"; ')


def status():
    """SimpleNamespace(state, version, text). state: "ready"; "no_database" (installed, `plannotate setupdb` not run);
    "missing" (WSL works, pLannotate is not installed); "no_wsl". Takes a second or two (it starts WSL)."""
    try:
        done = _run(env_lookup() + f'[ -n "$E" ] && python {shlex.quote(runner.shell_path(HELPER))} --status || echo missing', timeout=120)
    except PlannotateError:
        return SimpleNamespace(state="no_wsl", version="", text=WSL_HELP)
    out = done.stdout.strip().splitlines()[-1:] if done.returncode == 0 else []
    word, _, version = (out[0] if out else "").partition(" ")
    if word == "ready":
        return SimpleNamespace(state="ready", version=version, text=f"{NAME} {version}")
    if word == "nodb":
        return SimpleNamespace(state="no_database", version=version, text=f"{NAME} {version} is installed but its databases are not: run the install again.")
    if word == "missing":
        return SimpleNamespace(state="missing", version="", text=f"{NAME} is not installed.")
    return SimpleNamespace(state="no_wsl", version="", text=WSL_HELP)


def tools_state():
    """(status(), {tool: installed}) of pLannotate and the Terminal's tools (oDNA-cli wsltools, Help > Run Setup Again). Without WSL
    every tool is False and WSL is not started a second time to ask."""
    state, tools = status(), runner.terminal_tools()
    have = set(runner.installed_tools(tools)) if state.state != "no_wsl" else set()
    return state, {t: t in have for t in tools}


def install_script(plannotate=True, terminal=True):
    """bash that installs pLannotate (conda package, then its databases with `plannotate setupdb`) and/or the Terminal's tools.
    Each part is tried even when the other failed; the exit status is that of the worst."""
    parts = []
    if plannotate:
        parts.append(runner.install(["plannotate"], env=ENV, after='"$E/bin/plannotate" setupdb'))
    if terminal:
        parts.append(runner.install(runner.terminal_tools()))
    return "r=0; " + "".join(f"{p}; r=$((r|$?)); " for p in parts) + "exit $r"


def install(plannotate=True, terminal=True, out=None):
    """Run the install in the foreground, its text going to the file `out` (default: standard error). Returns the exit code."""
    return subprocess.run(runner.command_line(["-lc", install_script(plannotate, terminal)]), stdout=out or sys.stderr).returncode


def annotate(sequence, circular, run=None):
    """pLannotate's annotations of `sequence`: [{name, kind, start, end, strand, identity, match, database, fragment, note}] (see
    plannotate_run.rows_from), the same call returning the same rows from a cache. `run(work)`: how the search is run (in a thread
    behind a busy dialog); by default, here. Raises PlannotateError (a message for the user)."""
    if len(sequence) > MAX_LENGTH:
        raise PlannotateError(f"{NAME} annotates engineered plasmids of up to {MAX_LENGTH:,} bp; this sequence has {len(sequence):,}.")
    if _missing:    # asked once, and not again until `forget`: every change of a setting would start WSL for nothing
        raise NotInstalled(f"{NAME} is not installed.")
    key = (hashlib.sha1(sequence.encode()).hexdigest(), bool(circular))
    if key not in _cache:
        try:
            _cache[key] = (run or (lambda work: work()))(lambda: _search(sequence, circular))
        except NotInstalled:
            _missing.append(True)
            raise
    return _cache[key]


def _search(sequence, circular):
    with tempfile.TemporaryDirectory(prefix="opendna_plannotate_") as folder:
        job, result = Path(folder) / "job.json", Path(folder) / "result.json"
        job.write_text(json.dumps({"sequence": sequence, "circular": bool(circular)}), encoding="utf-8")
        script = env_lookup() + (f'[ -n "$E" ] || {{ echo "{NAME} is not installed" >&2; exit 3; }}; python {shlex.quote(runner.shell_path(HELPER))} '
                                 f'{shlex.quote(runner.shell_path(job))} {shlex.quote(runner.shell_path(result))}')
        done = _run(script, timeout=1800)
        if done.returncode == 3:
            raise NotInstalled(f"{NAME} is not installed.")
        if not result.exists():
            raise PlannotateError((done.stderr or done.stdout or f"exit code {done.returncode}").strip()[-400:])
        answer = json.loads(result.read_text(encoding="utf-8"))
    if "error" in answer:
        raise PlannotateError(answer["error"])
    return answer["rows"]


def forget():
    """After an install: ask again, and search again."""
    _cache.clear()
    _missing.clear()
