"""The logic behind the Terminal tab (terminal.py): pipelines read from tools.toml turned into bash,
paths between Windows and WSL, what `open` does with a file, and the QC summary of a run.

On Windows the session runs in WSL (`wsl -e bash`), with Windows paths turned into /mnt/<drive>/…
paths; elsewhere in bash directly. Tools are looked up first in a dedicated conda env
(opendna-tools, see install()), then on the normal PATH.

Each pipeline run gets a folder next to its output (<out>_run/) with its inputs, the exact command
(command.txt), the log (log.txt) and the tools' own reports, so it can be re-run by hand.
"""
import json
import re
import shlex
import sys
import tomllib
from datetime import datetime
from pathlib import Path

RECIPES = Path(__file__).with_name("tools.toml")
ENV = "opendna-tools"
# Where the tools env usually is (Miniforge, Miniconda, Mambaforge); searched before PATH.
TOOL_DIRS = [f"$HOME/{base}/envs/{ENV}/bin" for base in ("miniforge3", "miniconda3", "mambaforge")]
FASTQ_FILTER = "Reads (*.fastq *.fq *.fastq.gz *.fq.gz *.fasta *.fa *.fa.gz);;All files (*)"
SEQ_FILTER = "Sequences (*.gb *.gbk *.gbff *.fa *.fasta *.fna *.dna);;All files (*)"
TRACK_EXT = (".bam", ".vcf", ".vcf.gz", ".gff", ".gff3", ".gtf", ".bed")
SEQ_EXT = (".fa", ".fasta", ".fna", ".fas", ".gb", ".gbk", ".gbff", ".dna")
FASTA_EXT = (".fa", ".fasta", ".fna", ".fas")
MSA_SNIFF_LIMIT = 50_000_000   # bytes: a bigger FASTA is not read to look for an alignment
PLACEHOLDER = re.compile(r"\{(\w+)\}")


def pipelines():
    with open(RECIPES, "rb") as fh:
        return tomllib.load(fh)


def use_wsl():
    return sys.platform == "win32"


def shell_path(path):
    """A path as the shell sees it: on Windows V:\\data\\x.fq -> /mnt/v/data/x.fq (WSL)."""
    path = str(Path(path).resolve())
    if not use_wsl():
        return path
    if len(path) < 2 or path[1] != ":":
        raise ValueError(f"{path}: only paths on a drive letter can be reached from WSL")
    return f"/mnt/{path[0].lower()}{path[2:].replace(chr(92), '/')}"


def windows_path(path):
    """The reverse of shell_path for a /mnt/<drive>/… path; any other path is returned as it is."""
    m = re.match(r"^/mnt/([a-z])(/.*)?$", path)
    if not use_wsl() or not m:
        return path
    return f"{m.group(1).upper()}:{(m.group(2) or '/').replace('/', chr(92))}"


def command_line(args):
    """argv that starts bash (in WSL on Windows) with `args` (e.g. ["-lc", text], or ["-ls"] for a session)."""
    return (["wsl", "-e", "bash"] if use_wsl() else ["bash"]) + list(args)


def fill(template, values):
    """Put values in for {name} placeholders; unknown names (bash's ${R2:+...}) are left alone."""
    return PLACEHOLDER.sub(lambda m: str(values[m.group(1)]) if m.group(1) in values else m.group(0), template)


def options_for(pipeline, preset):
    """{option: value} of a preset: its own values, the options' defaults for the rest."""
    chosen = pipeline.get("presets", {}).get(preset, {})
    return {key: chosen.get(key, opt["default"]) for key, opt in pipeline.get("options", {}).items()}


def inputs_for(pipeline, preset):
    return pipeline.get("presets", {}).get(preset, {}).get("inputs", pipeline["inputs"])


INPUT_NAMES = {"reads1": "Reads", "reads2": "Mate reads", "long": "Long reads", "reference": "Reference", "sequences": "Sequences"}
ASSIGN_ONLY = re.compile(r"^\s*(\w+=\S*\s*;?\s*)+$")


def diagram(pipeline, options, inputs):
    """The flow chart of a pipeline as the Terminal draws it: [(kind, text)] in order, kind "input", "step",
    "skipped" (a choice that does nothing, e.g. cleanup off) or "output". Read from tools.toml alone, so a new
    pipeline gets its diagram with no code: a step that is one {choice} becomes "<option label>: <choice>",
    each tool of `needs` in a step becomes "<tool> <subcommand>", the other choices of the step are added in
    brackets, `open` is the output, and steps that only set variables are left out."""
    opts = pipeline.get("options", {})
    nodes = [("input", INPUT_NAMES.get(name, name)) for name in inputs]
    for step in pipeline["steps"]:
        if ASSIGN_ONLY.match(step):
            continue
        if step.lstrip().startswith("open "):
            nodes.append(("output", "Open in openDNA"))
            continue
        used = [k for k in PLACEHOLDER.findall(step) if "choices" in opts.get(k, {})]
        if step.strip() in (f"{{{k}}}" for k in used):
            key = used[0]
            nodes.append(("skipped" if opts[key]["choices"][options[key]].strip() == "true" else "step",
                          f"{opts[key]['label']}: {options[key]}"))
            continue
        for key in used:
            if step.lstrip().startswith(f"{{{key}}}"):   # {aligner} | samtools sort …: the choice, then the rest
                nodes.append(("step", f"{opts[key]['label']}: {options[key]}"))
        rest = PLACEHOLDER.sub(" ", step)
        extra = ", ".join(f"{opts[k]['label']}: {options[k]}" for k in used if not step.lstrip().startswith(f"{{{k}}}"))
        for m in re.finditer(r"(?:^|[;|&(])\s*(" +"|".join(map(re.escape, pipeline.get("needs", []))) + r")(?: ([a-z]\w*))?", rest):
            nodes.append(("step", " ".join(filter(None, m.groups())) + (f" ({extra})" if extra else "")))
    return nodes


def build(pipeline, options, files, threads, folder):
    """(bash command, tools it uses) for a run: steps joined with &&, run in `folder` in a subshell
    (so the session's own directory and variables stay as they were). `files`: {placeholder: path or ""}."""
    values = {key: shlex.quote(shell_path(v)) if v else "''" for key, v in files.items()}
    values["threads"] = int(threads)
    for key, opt in pipeline.get("options", {}).items():
        value = options[key]
        values[key] = opt["choices"][value] if "choices" in opt else int(value) if isinstance(opt["default"], int) else value
    # two passes: a choice's text holds placeholders of its own ({threads}, {min_qual}…)
    steps = [fill(fill(step, values), values) for step in pipeline["steps"]]
    body = " && ".join(steps)
    used = [t for t in pipeline.get("needs", []) if re.search(rf"(^|[\s;|&(]){re.escape(t)}\s", body)]
    need = f"need {' '.join(used)} && " if used else ""
    return f"( {need}cd {shlex.quote(shell_path(folder))} && {body} )", used


def install(tools):
    """bash that installs `tools` from bioconda into the opendna-tools env (created if missing), installing
    Miniforge into ~/miniforge3 first when no conda is found."""
    names = " ".join(shlex.quote(t) for t in tools)
    miniforge = "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-$(uname -m).sh"
    return ('( c=; for x in "$HOME"/miniforge3/bin/conda "$HOME"/miniconda3/bin/conda "$HOME"/mambaforge/bin/conda; do '
            '[ -x "$x" ] && { c=$x; break; }; done; '
            'if [ -z "$c" ]; then echo "No conda found: installing Miniforge into ~/miniforge3"; '
            f'curl -fsSL -o /tmp/miniforge.sh "{miniforge}" && bash /tmp/miniforge.sh -b -p "$HOME/miniforge3" '
            '&& rm -f /tmp/miniforge.sh && c="$HOME/miniforge3/bin/conda" || exit 1; fi; '   # exit leaves the ( ) only, not the session
            f'if [ -d "${{c%/bin/conda}}/envs/{ENV}" ]; then "$c" install -y -n {ENV} -c conda-forge -c bioconda {names}; '
            f'else "$c" create -y -n {ENV} -c conda-forge -c bioconda {names}; fi )')


def session_setup(folder):
    """What the session runs first: its PID (for Stop), tools on PATH, the `need` and `open` commands, the start folder."""
    return "\n".join([
        # output through a pipe: in WSL it is a socket otherwise, and tools that open /dev/stderr (MAFFT) fail on it
        'exec > >(cat) 2>&1',
        'echo "__PID__ $$ $!"',
        f'export PATH="{":".join(TOOL_DIRS)}:$PATH"; set -o pipefail',
        'need() { local m=; for t; do command -v "$t" >/dev/null || m="$m $t"; done; '
        '[ -z "$m" ] || { echo "Missing tools:$m. Press Install in the pipeline form."; return 1; }; }',
        # open: the app catches the marker; a path outside /mnt goes through wslpath so Windows can read it
        'open() { local tag=__OPEN__; [ "$1" = --msa ] && { tag=__OPENMSA__; shift; }; local f p; for f; do '
        'p=$(realpath -- "$f") || return 1; case "$p" in /mnt/*) ;; *) command -v wslpath >/dev/null && p=$(wslpath -w "$p");; esac; '
        'echo "$tag $p"; done; }',
        f"cd {shlex.quote(shell_path(folder))} 2>/dev/null",
        "",
    ])


def wrap(command):
    """A user's line as the session runs it: eval (a syntax error must not end the shell), no stdin (a
    program waiting for input would read the next commands), then the exit status after a marker."""
    return f"eval {shlex.quote(command)} </dev/null; printf '\\n__DONE__ %s\\n' $?\n"


def stop_children(pid, keep):
    """bash that stops every program the session started (all descendants of its shell), not the shell
    and not `keep` (the cat that carries the session's output)."""
    pid, keep = int(pid), int(keep)
    return (f"d() {{ for c in $(pgrep -P $1); do [ $c = {keep} ] && continue; echo $c; d $c; done; }}; "
            f"kill -TERM $(d {pid}) 2>/dev/null; true")


def open_kind(path, msa=False):
    """What `open` does with a file: "track" (genome view), "msa" (Multiple Alignment view), "document"
    (a new window) or "browser" (the system's default program: HTML reports and anything else)."""
    name = str(path).lower()
    if msa:
        return "msa"
    if name.endswith(TRACK_EXT):
        return "track"
    if name.endswith(FASTA_EXT) and Path(path).is_file() and Path(path).stat().st_size <= MSA_SNIFF_LIMIT:
        if is_alignment(Path(path).read_text(encoding="utf-8", errors="replace")):
            return "msa"
    return "document" if name.endswith(SEQ_EXT) else "browser"


def is_alignment(fasta_text):
    """An aligned FASTA: at least two records, all of one length, with gaps."""
    names, rows = [], []
    for line in fasta_text.splitlines():
        if line.startswith(">"):
            names.append(line)
            rows.append([])
        elif rows:
            rows[-1].append(line.strip())
    rows = ["".join(r) for r in rows]
    return len(rows) >= 2 and len(set(map(len, rows))) == 1 and any("-" in r for r in rows)


def summary(folder):
    """The QC summary lines of a run, read from what its tools left in `folder`."""
    folder, lines = Path(folder), []
    fastp = folder / "fastp.json"
    if fastp.is_file():
        s = json.loads(fastp.read_text(encoding="utf-8"))["summary"]
        before, after = s["before_filtering"]["total_reads"], s["after_filtering"]["total_reads"]
        lines.append(f"fastp: {before:,} reads in, {after:,} kept ({after / max(before, 1):.0%}), "
                     f"Q30 {s['after_filtering']['q30_rate']:.0%}")
    flagstat = folder / "flagstat.txt"
    if flagstat.is_file():
        m = re.search(r"^(\d+) \+ \d+ mapped \(([\d.]+)%", flagstat.read_text(encoding="utf-8"), re.M)
        if m:
            lines.append(f"samtools flagstat: {int(m.group(1)):,} reads mapped ({m.group(2)}%)")
    quast = folder / "quast" / "report.tsv"
    if quast.is_file():
        rows = dict(line.split("\t")[:2] for line in quast.read_text(encoding="utf-8").splitlines() if "\t" in line)
        lines.append(f"QUAST: {rows.get('# contigs', '?')} contigs, {int(rows.get('Total length', 0)):,} bp, "
                     f"N50 {int(rows.get('N50', 0)):,} bp")
    return lines


def reports(folder):
    """The tools' own HTML reports in a run folder."""
    return sorted(str(p) for p in Path(folder).rglob("*.html") if "icarus" not in p.parts[-2] + p.name)


def default_output(pipeline, doc, reads1, doc_path):
    """The output path a run suggests: the pipeline's `output` name beside the reads (or the document)."""
    stem = Path(reads1).name.split(".")[0] if reads1 else ""
    for tag in ("_R1", "_1"):  # sample_R1.fastq.gz -> sample
        stem = stem[: -len(tag)] if stem.endswith(tag) else stem
    name = pipeline["output"].format(doc=safe_name(doc), reads=stem or "reads")
    base = Path(reads1).parent if reads1 else Path(doc_path).parent if doc_path else Path.home()
    return str(base / name)


def safe_name(name):
    """A sequence name usable in a FASTA header and a file name."""
    return "".join(c if c.isalnum() or c in "._-" else "_" for c in name) or "sequence"


def run_folder(out):
    path = Path(out)
    folder = path.with_name((path.stem if path.suffix else path.name) + "_run")
    if folder.exists() and any(folder.iterdir()):
        folder = folder.with_name(folder.name + datetime.now().strftime("_%Y%m%d_%H%M%S"))
    return folder
