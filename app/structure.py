"""Protein structures: look one up or predict one from the amino acids, keep it in a cache, and make the page that shows it in 3D.

Sources (`SOURCES`; "Auto" tries them in this order, stopping at the first that works):
- cache: a structure fetched before for the same sequence (and accession) is read from `<data folder>/structures`;
- AlphaFold DB: a predicted structure for a UniProt accession (alphafold.ebi.ac.uk; only the accession is sent);
- PDB: RCSB's sequence search for an experimental (or modelled) structure of 90 % identity or more, then its PDB file;
- ESMFold: a prediction from the sequence alone (the ESM Metagenomic Atlas API; up to 400 amino acids);
- Local command: a command of the user's (AlphaFold, ColabFold, Boltz... installed here) with {fasta} and {pdb} in it.
Every source gives PDB text; AlphaFold DB and ESMFold put the pLDDT confidence (0-100) in the B-factor column, which the page colours by.
The network calls are the functions `fetch` (replaced in the tests) and run in a worker thread of the caller.
"""
import hashlib
import json
import re
import shlex
import subprocess
import tempfile
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace

from primer_db import data_dir

SOURCES = ["Auto", "PDB", "AlphaFold DB", "ESMFold", "Local command"]
NOTICE = {
    "Auto": "Tries the cache, then AlphaFold DB (if there is a UniProt accession), then PDB, then ESMFold: so the sequence may be "
            "sent to RCSB PDB and to the ESM Atlas (Meta).",
    "PDB": "The sequence is sent to RCSB PDB (search.rcsb.org) to look for a matching structure.",
    "AlphaFold DB": "Only the UniProt accession is sent (to alphafold.ebi.ac.uk).",
    "ESMFold": "The sequence is sent to api.esmatlas.com (the ESM Metagenomic Atlas of Meta) to be predicted.",
    "Local command": "Nothing leaves this computer: your command runs here.",
}
ESM_MAX = 400          # amino acids the ESM Atlas accepts
PDB_MIN = 25           # amino acids RCSB's sequence search accepts
TIMEOUT, ESM_TIMEOUT, LOCAL_TIMEOUT = 30, 300, 3600   # seconds
USER_AGENT = "OpenGenEdit"
AMINO = set("ACDEFGHIKLMNPQRSTVWY")
VENDOR = Path(__file__).with_name("vendor") / "3Dmol-min.js"
PLDDT_BANDS = [(90, "#0053D6", "very high (> 90)"), (70, "#65CBF3", "confident (70-90)"), (50, "#FFDB13", "low (50-70)"),
               (0, "#FF7D45", "very low (< 50)")]   # AlphaFold DB's colours


class StructureError(Exception):
    """Something the user can read: why no structure came."""


@dataclass
class Structure:
    pdb: str
    source: str          # which source gave it
    detail: str          # "AlphaFold DB AF-P69905-F1, mean pLDDT 98.1"
    plddt: bool          # the B-factor column holds pLDDT (0-100)
    cached: bool = False


# ---- the network ---------------------------------------------------------------------------------------------------------------------


def fetch(url, data=None, headers=None, timeout=TIMEOUT):
    """The bytes at `url` (a POST when `data` is given). StructureError when it cannot be reached or says no."""
    request = urllib.request.Request(url, data=data, headers={"User-Agent": USER_AGENT, **(headers or {})})
    host = url.split("/")[2]
    try:
        with urllib.request.urlopen(request, timeout=timeout) as reply:
            return reply.read()
    except urllib.error.HTTPError as e:
        raise StructureError(f"{host} answered {e.code} {e.reason}.") from None
    except (urllib.error.URLError, TimeoutError, OSError) as e:
        raise StructureError(f"Could not reach {host}: {getattr(e, 'reason', e)}") from None


def clean(sequence):
    """The amino acids of `sequence` as capital letters; StructureError when it has others or is empty."""
    seq = "".join(sequence.split()).upper().rstrip("*")
    if not seq:
        raise StructureError("There are no amino acids.")
    bad = sorted(set(seq) - AMINO)
    if bad:
        raise StructureError(f"Letters that are not amino acids (or are unknown): {' '.join(bad)}.")
    return seq


def accession_of(qualifiers):
    """The UniProt accession in a feature's db_xref (UniProtKB/Swiss-Prot:P69905), or ''."""
    for key, value in qualifiers:
        if key == "db_xref" and value.lower().startswith("uniprotkb"):
            return value.split(":", 1)[-1].strip()
    return ""


def scale_bfactors(pdb, factor):
    """The PDB text with every B-factor multiplied by `factor` (columns 61-66)."""
    out = []
    for line in pdb.splitlines():
        if line.startswith(("ATOM", "HETATM")) and len(line) >= 66:
            try:
                line = line[:60] + f"{float(line[60:66]) * factor:6.2f}" + line[66:]
            except ValueError:
                pass
        out.append(line)
    return "\n".join(out) + "\n"


def mean_bfactor(pdb):
    """Mean B-factor of the alpha carbons (the mean pLDDT of a prediction), or None."""
    values = [float(line[60:66]) for line in pdb.splitlines() if line.startswith("ATOM") and line[12:16] == " CA "]
    return sum(values) / len(values) if values else None


def looks_like_pdb(text):
    return any(line.startswith("ATOM") for line in text.splitlines()[:5000])


# ---- the sources ---------------------------------------------------------------------------------------------------------------------


def from_alphafold(accession):
    acc = accession.strip().upper()
    if not re.fullmatch(r"[A-Z0-9]{6,10}", acc):
        raise StructureError("AlphaFold DB needs a UniProt accession (such as P69905).")
    try:
        entries = json.loads(fetch(f"https://alphafold.ebi.ac.uk/api/prediction/{acc}"))
        entry = entries[0]
        url = entry["pdbUrl"]
        if not (isinstance(url, str) and url.startswith("https://alphafold.ebi.ac.uk/")):
            raise StructureError("AlphaFold DB's answer points somewhere else, so it was not fetched.")   # never follow a reply to another host or file://
        pdb = fetch(url).decode("utf-8", "replace")
        if not looks_like_pdb(pdb):
            raise StructureError("AlphaFold DB did not return a structure.")
    except StructureError as e:
        raise StructureError(f"AlphaFold DB has no model for {acc} ({e})") from None
    except (ValueError, KeyError, IndexError, TypeError):
        raise StructureError(f"AlphaFold DB's answer for {acc} could not be read.") from None
    return Structure(pdb, "AlphaFold DB", f"AlphaFold DB {entry.get('modelEntityId', acc)}, mean pLDDT {entry.get('globalMetricValue', '?')}", True)


def from_pdb(sequence):
    if len(sequence) < PDB_MIN:
        raise StructureError(f"RCSB PDB looks up sequences of {PDB_MIN} amino acids or more.")
    query = {"query": {"type": "terminal", "service": "sequence",
                       "parameters": {"evalue_cutoff": 1, "identity_cutoff": 0.9, "sequence_type": "protein", "value": sequence}},
             "return_type": "polymer_entity", "request_options": {"paginate": {"start": 0, "rows": 5}, "scoring_strategy": "sequence"}}
    try:
        answer = fetch("https://search.rcsb.org/rcsbsearch/v2/query", json.dumps(query).encode(), {"Content-Type": "application/json"})
        found = json.loads(answer).get("result_set", []) if answer else []   # no answer at all: nothing matches
    except (StructureError, ValueError) as e:
        raise StructureError(f"The PDB search failed ({e})") from None
    for entry in found:
        code = entry["identifier"].split("_")[0]
        try:
            pdb = fetch(f"https://files.rcsb.org/download/{code}.pdb").decode("utf-8", "replace")
        except StructureError:
            continue   # large structures have no PDB-format file: try the next
        if looks_like_pdb(pdb):
            return Structure(pdb, "PDB", f"PDB {code} (sequence identity 90 % or more; the whole entry, with its ligands and other chains)", False)
    raise StructureError("No structure in the PDB matches this sequence (90 % identity or more).")


def from_esmfold(sequence):
    if len(sequence) > ESM_MAX:
        raise StructureError(f"ESMFold (the ESM Atlas) predicts up to {ESM_MAX} amino acids; this has {len(sequence):,}.")
    pdb = fetch("https://api.esmatlas.com/foldSequence/v1/pdb/", sequence.encode(), {"Content-Type": "text/plain"}, ESM_TIMEOUT).decode("utf-8", "replace")
    if not looks_like_pdb(pdb):
        raise StructureError("The ESM Atlas did not return a structure.")
    if (mean_bfactor(pdb) or 0) <= 1.0:
        pdb = scale_bfactors(pdb, 100)   # it gives pLDDT as 0-1; the page and AlphaFold use 0-100
    return Structure(pdb, "ESMFold", f"ESMFold prediction, mean pLDDT {mean_bfactor(pdb) or 0:.1f}", True)


def from_local(sequence, command, job=None):
    """Run the user's `command` (words with {fasta} and {pdb}) and read the PDB file it writes. `job` (optional): a namespace
    whose `proc` is set to the running process, so another thread can stop it (a prediction can run for an hour); it sets
    `cancelled` first, and then this raises StructureError("Cancelled.")."""
    if "{fasta}" not in command or "{pdb}" not in command:
        raise StructureError("The local command needs {fasta} (the input file) and {pdb} (the file it writes) in it.")
    with tempfile.TemporaryDirectory(prefix="genomeedit_fold_") as folder:
        fasta, out = Path(folder) / "protein.fasta", Path(folder) / "model.pdb"
        fasta.write_text(f">protein\n{sequence}\n")
        # posix=False keeps Windows paths whole, and keeps the quotes round a word: take them off
        words = [w.strip('"').replace("{fasta}", str(fasta)).replace("{pdb}", str(out)) for w in shlex.split(command, posix=False)]
        try:
            proc = subprocess.Popen(words, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, cwd=folder)
        except OSError as e:
            raise StructureError(f"The local command did not run: {e}") from None
        if job is not None:
            job.proc = proc
            if job.cancelled:   # cancelled while it was starting
                proc.terminate()
        try:
            _, stderr = proc.communicate(timeout=LOCAL_TIMEOUT)
        except subprocess.TimeoutExpired as e:
            proc.kill()
            proc.communicate()
            raise StructureError(f"The local command did not run: {e}") from None
        if job is not None and job.cancelled:
            raise StructureError("Cancelled.")
        done = SimpleNamespace(returncode=proc.returncode, stderr=stderr or "")
        if not out.exists() or not looks_like_pdb(out.read_text(errors="replace")):
            raise StructureError(f"The local command wrote no PDB file (exit code {done.returncode}). {done.stderr.strip()[-300:]}")
        return Structure(out.read_text(errors="replace"), "Local command", f"Local command ({words[0]})", False)


# ---- the cache ------------------------------------------------------------------------------------------------------------------------


def cache_dir():
    path = data_dir() / "structures"
    path.mkdir(exist_ok=True)
    return path


def cache_key(sequence, accession, source):
    digest = hashlib.sha256(f"{sequence}|{accession.strip().upper() if source == 'AlphaFold DB' else ''}".encode()).hexdigest()[:20]
    return cache_dir() / f"{digest}_{re.sub('[^A-Za-z]', '', source)}"


def cached(sequence, accession, source):
    base = cache_key(sequence, accession, source)
    try:
        meta = json.loads(base.with_suffix(".json").read_text())
        return Structure(base.with_suffix(".pdb").read_text(), meta["source"], meta["detail"], meta["plddt"], True)
    except (OSError, ValueError, KeyError):
        return None


def store(sequence, accession, structure):
    base = cache_key(sequence, accession, structure.source)
    base.with_suffix(".pdb").write_text(structure.pdb)
    base.with_suffix(".json").write_text(json.dumps({"source": structure.source, "detail": structure.detail, "plddt": structure.plddt}))


# ---- asking ---------------------------------------------------------------------------------------------------------------------------


def get(sequence, source="Auto", accession="", command="", use_cache=True, job=None):
    """A Structure for the amino acids `sequence` from `source` (see SOURCES); StructureError says why not. A local command's
    results are not cached (the command may change); `accession` is a UniProt accession for AlphaFold DB; `job`: see
    from_local (a local command can be stopped through it)."""
    seq = clean(sequence)
    order = ["AlphaFold DB", "PDB", "ESMFold"] if source == "Auto" else [source]
    errors = []
    for name in order:
        if use_cache and name != "Local command" and (hit := cached(seq, accession, name)):
            return hit
        try:
            if name == "AlphaFold DB":
                result = from_alphafold(accession)
            elif name == "PDB":
                result = from_pdb(seq)
            elif name == "ESMFold":
                result = from_esmfold(seq)
            else:
                return from_local(seq, command, job)
        except StructureError as e:
            if source != "Auto":
                raise
            if not (name == "AlphaFold DB" and not accession.strip()):   # "no accession" is not worth a line in the report
                errors.append(f"{name}: {e}")
            continue
        store(seq, accession, result)
        return result
    raise StructureError("No structure found.\n" + "\n".join(errors))


# ---- the 3D page ---------------------------------------------------------------------------------------------------------------------


def viewer_html(structure, title):
    """A self-contained page (3Dmol.js inlined) that shows the structure: cartoon coloured by pLDDT / chain / rainbow, with buttons."""
    from html import escape
    library = VENDOR.read_text(encoding="utf-8").replace("</script", "<\\/script")
    bands = json.dumps([[lo, colour] for lo, colour, _ in PLDDT_BANDS])
    # The PDB text sits inside a <script>: no "<" may stay in it ("</script>" or "<!--" in a hostile file or a remark line would change
    # how the page is read); inside a JSON string the six characters backslash, u, 0, 0, 3, c mean "<".
    pdb_json = json.dumps(structure.pdb).replace("<", chr(92) + "u003c")
    legend = "".join(f"<span><i style='background:{c}'></i>{escape(label)}</span>" for _, c, label in PLDDT_BANDS)
    return f"""<!doctype html><html><head><meta charset="utf-8"><title>{escape(title)}</title>
<style>body{{margin:0;font:14px system-ui,sans-serif}}#bar{{padding:6px 10px;background:#f2f2f2;border-bottom:1px solid #ccc}}
#bar button{{margin-right:4px}}#v{{position:absolute;top:64px;bottom:0;left:0;right:0}}
#legend span{{margin-right:12px;font-size:12px}}#legend i{{display:inline-block;width:10px;height:10px;margin-right:4px}}</style></head><body>
<div id="bar"><b>{escape(title)}</b> &nbsp; {escape(structure.detail)}<br>
<button onclick="paint('plddt')">Colour by confidence</button><button onclick="paint('spectrum')">Rainbow</button>
<button onclick="paint('chain')">By chain</button><button onclick="shape('cartoon')">Cartoon</button>
<button onclick="shape('stick')">Sticks</button><button onclick="shape('surface')">Surface</button>
<span id="legend" style="display:{'inline' if structure.plddt else 'none'}">{legend}</span></div><div id="v"></div>
<script>{library}</script>
<script>
const pdb = {pdb_json}, bands = {bands}, hasConfidence = {json.dumps(structure.plddt)};
const viewer = $3Dmol.createViewer("v", {{backgroundColor: "white"}});
viewer.addModel(pdb, "pdb");
let mode = hasConfidence ? "plddt" : "spectrum", form = "cartoon";
function band(b) {{ for (const [lo, c] of bands) if (b >= lo) return c; return bands[bands.length - 1][1]; }}
function colour() {{ return mode === "plddt" ? {{colorfunc: a => band(a.b)}} : {{colorscheme: mode === "chain" ? "chain" : "spectrum"}}; }}
function draw() {{
  viewer.setStyle({{}}, form === "stick" ? {{stick: colour()}} : {{cartoon: colour()}});
  viewer.removeAllSurfaces();
  if (form === "surface") {{ viewer.setStyle({{}}, {{cartoon: colour()}}); viewer.addSurface($3Dmol.SurfaceType.VDW, {{opacity: 0.8, colorfunc: mode === "plddt" ? (a => band(a.b)) : undefined, colorscheme: mode === "plddt" ? undefined : "spectrum"}}); }}
  viewer.render();
}}
function paint(m) {{ mode = m; draw(); }}
function shape(f) {{ form = f; draw(); }}
draw(); viewer.zoomTo(); viewer.render();
</script></body></html>"""
