"""Protein structures: look one up or predict one from the amino acids, keep it in a cache, and make the page that shows it in 3D.

Sources (`SOURCES`; "Auto" tries them in this order, stopping at the first that works):
- cache: a structure fetched before for the same sequence (and accession) is read from `<data folder>/structures`;
- AlphaFold DB: a predicted structure for a UniProt accession (alphafold.ebi.ac.uk; only the accession is sent);
- PDB: RCSB's sequence search for an experimental (or modelled) structure of 90 % identity or more, then its PDB file (the window
  lists the entries found, down to 50 % identity, to pick one: search_pdb and pdb_entry);
- Linked entries (the window only): the PDB entries the feature links (db_xref "PDB:1BTL", as in ENA records) and those the PDB maps
  to its UniProt accession, listed to pick one (linked_pdb); only the codes and the accession are sent;
- ESMFold: a prediction from the sequence alone (the ESM Metagenomic Atlas API; up to 400 amino acids);
- Local command: a command of the user's (AlphaFold, ColabFold, Boltz... installed here) with {fasta} and {pdb} in it.
Every source gives PDB text; AlphaFold DB and ESMFold put the pLDDT confidence (0-100) in the B-factor column, which the page colours by.
The network calls are the functions `fetch` (replaced in the tests) and run in a worker thread of the caller.
"""
import hashlib
import json
import math
import re
import shlex
import subprocess
import tempfile
import urllib.error
import urllib.request
from dataclasses import dataclass, field
from difflib import SequenceMatcher
from pathlib import Path
from types import SimpleNamespace

from primer_db import data_dir
import project

SOURCES = ["Auto", "PDB", "Linked entries", "AlphaFold DB", "ESMFold", "Local command"]
NOTICE = {
    "Auto": "Tries the cache, then AlphaFold DB (if there is a UniProt accession), then PDB, then ESMFold: so the sequence may be "
            "sent to RCSB PDB and to the ESM Atlas (Meta).",
    "PDB": "The sequence is sent to RCSB PDB (search.rcsb.org) to look for a matching structure.",
    "Linked entries": "Lists the PDB entries the feature links (PDB:…) and those the PDB maps to the UniProt accession: only the codes "
                      "and the accession are sent to RCSB, the sequence stays here.",
    "AlphaFold DB": "Only the UniProt accession is sent (to alphafold.ebi.ac.uk).",
    "ESMFold": "The sequence is sent to api.esmatlas.com (the ESM Metagenomic Atlas of Meta) to be predicted.",
    "Local command": "Nothing leaves this computer: your command runs here.",
}
ESM_MAX = 400          # amino acids the ESM Atlas accepts
PDB_MIN = 25           # amino acids RCSB's sequence search accepts
TIMEOUT, ESM_TIMEOUT, LOCAL_TIMEOUT = 30, 300, 3600   # seconds
USER_AGENT = f"{project.NAME}/{project.VERSION}"   # AlphaFold DB answers 403 to a bare "openDNA"
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
    info: dict = field(default_factory=dict)   # what the database says, field -> text, in order ("Link": its page); see structure_view


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
    found = Structure(pdb, "AlphaFold DB", f"AlphaFold DB {entry.get('modelEntityId', acc)}, mean pLDDT {entry.get('globalMetricValue', '?')}", True)
    start, end = entry.get("uniprotStart"), entry.get("uniprotEnd")
    found.info = _fields(Protein=entry.get("uniprotDescription"), Source="AlphaFold DB (predicted)", Entry=entry.get("modelEntityId"),
                         UniProt=acc, Gene=entry.get("gene"), Organism=entry.get("organismScientificName"),
                         Residues=f"{start}–{end} of the UniProt sequence" if start and end else None,
                         **{"Mean pLDDT": entry.get("globalMetricValue"), "Model date": str(entry.get("modelCreatedDate") or "")[:10]},
                         Link=f"https://alphafold.ebi.ac.uk/entry/{acc}")
    return found


@dataclass
class PdbHit:
    """An entry of the PDB that holds the protein (or one like it), as search_pdb lists it."""
    entry: str                # "1BTL"
    identity: float           # 0-1, over the aligned part
    covered: float            # 0-1: the part of the protein the matching chain covers
    title: str = ""
    method: str = ""          # "X-RAY DIFFRACTION"...
    resolution: float = None  # Å
    organism: str = ""
    chains: str = ""          # the chains that are this protein ("A, C")


def search_pdb(sequence, identity=0.5, rows=25, details=True):
    """[PdbHit] of the PDB entries with a chain like the amino acids `sequence` (RCSB's sequence search, best first). `details`: also
    ask RCSB's data API for their titles, methods, resolutions, organisms and chains."""
    if len(sequence) < PDB_MIN:
        raise StructureError(f"RCSB PDB looks up sequences of {PDB_MIN} amino acids or more.")
    query = {"query": {"type": "terminal", "service": "sequence",
                       "parameters": {"evalue_cutoff": 1, "identity_cutoff": identity, "sequence_type": "protein", "value": sequence}},
             "return_type": "polymer_entity", "request_options": {"paginate": {"start": 0, "rows": rows}, "scoring_strategy": "sequence",
                                                                  "results_verbosity": "verbose"}}
    try:
        answer = fetch("https://search.rcsb.org/rcsbsearch/v2/query", json.dumps(query).encode(), {"Content-Type": "application/json"})
        found = json.loads(answer).get("result_set", []) if answer else []   # no answer at all: nothing matches
    except (StructureError, ValueError) as e:
        raise StructureError(f"The PDB search failed ({e})") from None
    hits, entities = [], []
    for item in found:
        try:
            match = item["services"][0]["nodes"][0]["match_context"][0]
            same, covered = float(match["sequence_identity"]), (match["query_end"] - match["query_beg"] + 1) / len(sequence)
        except (KeyError, IndexError, TypeError, ValueError):
            same, covered = float(item.get("score", 0)), 1.0
        code = str(item.get("identifier", "")).split("_")[0].upper()
        if re.fullmatch(r"[0-9A-Z]{4,12}", code) and code not in (h.entry for h in hits):   # (an entry with two matching entities: once)
            hits.append(PdbHit(code, same, min(covered, 1.0)))
            entities.append(item["identifier"])
    if details and hits:
        try:
            _pdb_details(hits, entities)
        except (StructureError, ValueError, KeyError, TypeError, AttributeError):
            pass   # the list without its titles is still worth showing
    return hits


def _pdb_details(hits, entities):
    """Fill in the titles, methods, resolutions, organisms and chains of `hits` (RCSB's data API, in one request)."""
    graphql = ("{ entries(entry_ids: %s) { rcsb_id struct { title } exptl { method } rcsb_entry_info { resolution_combined } } "
               "polymer_entities(entity_ids: %s) { rcsb_id rcsb_entity_source_organism { ncbi_scientific_name } "
               "rcsb_polymer_entity_container_identifiers { auth_asym_ids } } }") % (json.dumps([h.entry for h in hits]), json.dumps(entities))
    data = json.loads(fetch("https://data.rcsb.org/graphql", json.dumps({"query": graphql}).encode(), {"Content-Type": "application/json"}))["data"]
    by_code = {h.entry: h for h in hits}
    for e in data.get("entries") or []:
        h = by_code.get(e["rcsb_id"])
        if h:
            h.title = ((e.get("struct") or {}).get("title") or "").strip()
            h.method = ", ".join(x["method"] for x in e.get("exptl") or [] if x.get("method"))
            h.resolution = ((e.get("rcsb_entry_info") or {}).get("resolution_combined") or [None])[0]
    for e in data.get("polymer_entities") or []:
        h = by_code.get(e["rcsb_id"].split("_")[0])
        if h:
            h.organism = ", ".join(sorted({o["ncbi_scientific_name"] for o in e.get("rcsb_entity_source_organism") or []
                                           if o.get("ncbi_scientific_name")}))
            h.chains = ", ".join((e.get("rcsb_polymer_entity_container_identifiers") or {}).get("auth_asym_ids") or [])


def pdb_links(qualifiers):
    """The PDB codes a feature links in its db_xref ("PDB:1BTL"; ENA records list them), in order, each once."""
    codes = []
    for key, value in qualifiers:
        if key == "db_xref" and value.upper().startswith("PDB:"):
            code = re.match(r"[0-9][0-9A-Z]{3}", value[4:].strip().upper())
            if code and code.group() not in codes:
                codes.append(code.group())
    return codes


def similarity(protein, other):
    """(identity, covered) of the amino acids `other` against `protein`: identity over the stretch of the protein they line up with
    (the same letters / its length), covered = that stretch / the protein."""
    same, ends = 0, []
    for op, a1, a2, b1, b2 in SequenceMatcher(None, protein, other, autojunk=False).get_opcodes():
        if op == "equal" or (op == "replace" and a2 - a1 == b2 - b1):
            same += (a2 - a1) if op == "equal" else 0
            ends += [a1, a2]
    if not ends:
        return 0.0, 0.0
    span = max(ends) - min(ends)
    return same / span, span / len(protein)


def linked_pdb(sequence, codes=(), accession="", most=100):
    """[PdbHit] of the PDB entries linked to the protein: `codes` (the feature's PDB links) and the entries whose chains the PDB maps
    to the UniProt `accession` (curated by UniProt and the PDB: SIFTS), best match first. Each is compared with `sequence` here, on its
    chain that matches best; only the codes and the accession are sent.
    shortcut: the first `most` entries only (a well studied protein has hundreds); raise it if one is missed."""
    codes = list(codes)
    acc = accession.strip().upper()
    if acc:
        if not re.fullmatch(r"[A-Z0-9]{6,10}", acc):
            raise StructureError("A UniProt accession is 6 to 10 letters and digits (such as P62593).")
        query = {"query": {"type": "group", "logical_operator": "and", "nodes": [
                     {"type": "terminal", "service": "text", "parameters": {
                         "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                         "operator": "exact_match", "value": acc}},
                     {"type": "terminal", "service": "text", "parameters": {
                         "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_name",
                         "operator": "exact_match", "value": "UniProt"}}]},
                 "return_type": "entry", "request_options": {"paginate": {"start": 0, "rows": most}}}
        try:
            answer = fetch("https://search.rcsb.org/rcsbsearch/v2/query", json.dumps(query).encode(), {"Content-Type": "application/json"})
            found = json.loads(answer).get("result_set", []) if answer else []
        except (StructureError, ValueError) as e:
            raise StructureError(f"The PDB could not be asked about {acc} ({e})") from None
        codes += [str(item.get("identifier", "")).upper() for item in found]
    codes = [c for i, c in enumerate(codes) if re.fullmatch(r"[0-9A-Z]{4,12}", c) and c not in codes[:i]][:most]
    if not codes:
        return []
    graphql = ("{ entries(entry_ids: %s) { rcsb_id struct { title } exptl { method } rcsb_entry_info { resolution_combined } "
               "polymer_entities { entity_poly { pdbx_seq_one_letter_code_can } rcsb_entity_source_organism { ncbi_scientific_name } "
               "rcsb_polymer_entity_container_identifiers { auth_asym_ids } } } }") % json.dumps(codes)
    try:
        data = json.loads(fetch("https://data.rcsb.org/graphql", json.dumps({"query": graphql}).encode(),
                                {"Content-Type": "application/json"}))["data"]
        hits = []
        for e in data.get("entries") or []:
            best = None
            for entity in e.get("polymer_entities") or []:
                chain_seq = "".join(((entity.get("entity_poly") or {}).get("pdbx_seq_one_letter_code_can") or "").split())
                score = similarity(sequence, chain_seq)
                if best is None or score > best[0]:
                    best = (score, entity)
            if best is None:
                continue   # no protein chain (a nucleic acid only)
            (identity, covered), entity = best
            hits.append(PdbHit(
                e["rcsb_id"], identity, covered, ((e.get("struct") or {}).get("title") or "").strip(),
                ", ".join(x["method"] for x in e.get("exptl") or [] if x.get("method")),
                ((e.get("rcsb_entry_info") or {}).get("resolution_combined") or [None])[0],
                ", ".join(sorted({o["ncbi_scientific_name"] for o in entity.get("rcsb_entity_source_organism") or []
                                  if o.get("ncbi_scientific_name")})),
                ", ".join((entity.get("rcsb_polymer_entity_container_identifiers") or {}).get("auth_asym_ids") or [])))
    except (StructureError, ValueError, KeyError, TypeError, AttributeError) as e:
        raise StructureError(f"The linked PDB entries could not be described ({e})") from None
    hits.sort(key=lambda h: (-round(h.identity, 2), -round(h.covered, 2), h.resolution or 99))
    return hits


def pdb_entry(hit):
    """The Structure of the PDB entry of `hit` (a PdbHit): its PDB file, the whole entry with its ligands and other chains."""
    pdb = fetch(f"https://files.rcsb.org/download/{hit.entry}.pdb").decode("utf-8", "replace")
    if not looks_like_pdb(pdb):
        raise StructureError(f"PDB {hit.entry} has no PDB-format file (very large entries only come as mmCIF): pick another.")
    about = (f"{hit.identity:.0%} identity" + (f", {hit.method.lower()}" if hit.method else "")
             + (f", {hit.resolution:.2f} Å" if hit.resolution else "") + (f"; chain{'s' if ',' in hit.chains else ''} {hit.chains}" if hit.chains else ""))
    found = Structure(pdb, "PDB", f"PDB {hit.entry}" + (f": {hit.title}" if hit.title else "") + f" ({about})", False)
    found.info = _fields(Protein=hit.title, Source="RCSB PDB (experimental)", Entry=hit.entry, Method=hit.method.capitalize(),
                         Resolution=f"{hit.resolution:.2f} Å" if hit.resolution else None, Organism=hit.organism,
                         Identity=f"{hit.identity:.0%}", Chain_covers=f"{hit.covered:.0%} of the protein", Chains=hit.chains,
                         Link=f"https://www.rcsb.org/structure/{hit.entry}")
    return found


def _fields(**values):
    """The info of a Structure: the fields that have a value, in the order given (keyword names, "_" for a space)."""
    return {key.replace("_", " "): str(value) for key, value in values.items() if value not in (None, "", [])}


def from_pdb(sequence):
    for hit in search_pdb(sequence, 0.9, 5):
        try:
            return pdb_entry(hit)
        except StructureError:
            continue   # large structures have no PDB-format file: try the next
    raise StructureError("No structure in the PDB matches this sequence (90 % identity or more).")


def from_esmfold(sequence):
    if len(sequence) > ESM_MAX:
        raise StructureError(f"ESMFold (the ESM Atlas) predicts up to {ESM_MAX} amino acids; this has {len(sequence):,}.")
    pdb = fetch("https://api.esmatlas.com/foldSequence/v1/pdb/", sequence.encode(), {"Content-Type": "text/plain"}, ESM_TIMEOUT).decode("utf-8", "replace")
    if not looks_like_pdb(pdb):
        raise StructureError("The ESM Atlas did not return a structure.")
    if (mean_bfactor(pdb) or 0) <= 1.0:
        pdb = scale_bfactors(pdb, 100)   # it gives pLDDT as 0-1; the page and AlphaFold use 0-100
    found = Structure(pdb, "ESMFold", f"ESMFold prediction, mean pLDDT {mean_bfactor(pdb) or 0:.1f}", True)
    found.info = _fields(Source="ESMFold (predicted, ESM Atlas of Meta)", **{"Mean pLDDT": f"{mean_bfactor(pdb) or 0:.1f}"})
    return found


def from_local(sequence, command, job=None):
    """Run the user's `command` (words with {fasta} and {pdb}) and read the PDB file it writes. `job` (optional): a namespace
    whose `proc` is set to the running process, so another thread can stop it (a prediction can run for an hour); it sets
    `cancelled` first, and then this raises StructureError("Cancelled.")."""
    if "{fasta}" not in command or "{pdb}" not in command:
        raise StructureError("The local command needs {fasta} (the input file) and {pdb} (the file it writes) in it.")
    with tempfile.TemporaryDirectory(prefix="opendna_fold_") as folder:
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
        found = Structure(out.read_text(errors="replace"), "Local command", f"Local command ({words[0]})", False)
        found.info = _fields(Source="Local command (predicted here)", Command=command)
        return found


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
        return Structure(base.with_suffix(".pdb").read_text(), meta["source"], meta["detail"], meta["plddt"], True, meta.get("info", {}))
    except (OSError, ValueError, KeyError):
        return None


def store(sequence, accession, structure):
    base = cache_key(sequence, accession, structure.source)
    base.with_suffix(".pdb").write_text(structure.pdb)
    base.with_suffix(".json").write_text(json.dumps({"source": structure.source, "detail": structure.detail, "plddt": structure.plddt,
                                                     "info": structure.info}))


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


# ---- the model, for the 3D view of the Protein panel (structure_view.py) --------------------------------------------------------------

THREE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
         "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
         "MSE": "M"}   # selenomethionine: a methionine in many crystal structures
WATER = {"HOH", "WAT", "DOD"}
CHAIN_BREAK = 4.2      # Å: two alpha carbons further apart than this are not bonded (missing residues between them)


@dataclass
class Residue:
    chain: str
    number: str          # residue number and insertion code, as in the file ("42", "42A")
    aa: str              # one letter
    b: float = 0.0       # B-factor of the alpha carbon: pLDDT in a prediction
    ca: tuple = None     # (x, y, z) of the alpha carbon
    atoms: list = field(default_factory=list)   # [(name, element, (x, y, z))]
    ss: str = "-"        # H helix, E strand, - neither (see secondary)


def model_of(pdb):
    """(residues, ligand atoms) of the first model of the PDB text: the amino acids that have an alpha carbon, in file order (first
    alternate location only), and [(element, (x, y, z))] of the other HETATM records, water left out."""
    residues, ligands, by_key = [], [], {}
    for line in pdb.splitlines():
        if line.startswith("ENDMDL"):
            break
        if not line.startswith(("ATOM  ", "HETATM")) or len(line) < 54 or line[16] not in " A":
            continue
        try:
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except ValueError:
            continue
        name, kind = line[12:16].strip(), line[17:20].strip()
        element = (line[76:78].strip() if len(line) >= 78 else "") or name.lstrip("0123456789")[:1]
        if kind in THREE:
            key = (line[21], line[22:27].strip())
            r = by_key.get(key)
            if r is None:
                r = by_key[key] = Residue(line[21], key[1], THREE[kind])
                residues.append(r)
            r.atoms.append((name, element.upper(), xyz))
            if name == "CA":
                r.ca = xyz
                try:
                    r.b = float(line[60:66])
                except ValueError:
                    pass
        elif line.startswith("HETATM") and kind not in WATER:
            ligands.append((element.upper(), xyz))
    residues = [r for r in residues if r.ca]
    secondary(residues)
    return residues, ligands


def segments(residues):
    """[[indices]] of the runs of residues bonded one to the next (same chain, alpha carbons close enough)."""
    runs = []
    for k, r in enumerate(residues):
        if runs and residues[k - 1].chain == r.chain and math.dist(residues[k - 1].ca, r.ca) <= CHAIN_BREAK:
            runs[-1].append(k)
        else:
            runs.append([k])
    return runs


def secondary(residues):
    """Set the `ss` of each residue from the alpha carbons alone (the distances from i to i+2, i+3, i+4 of P-SEA, Labesse et al. 1997):
    a prediction has no HELIX/SHEET records, and one rule for every source keeps the pictures alike. A strand also needs a
    neighbouring strand (an extended loop alone is not a sheet); helices of fewer than 4 and strands of fewer than 3 are dropped."""
    helix, strand = ((5.5, 0.5), (5.3, 0.5), (6.4, 0.6)), ((6.7, 0.6), (9.9, 0.9), (12.4, 1.1))
    def fits(run, i, ranges):
        return i + 4 < len(run) and all(abs(math.dist(residues[run[i]].ca, residues[run[i + n]].ca) - mid) <= tol
                                        for n, (mid, tol) in zip((2, 3, 4), ranges))
    for run in segments(residues):
        for i in range(len(run) - 4):
            if fits(run, i, helix) and fits(run, i + 1, helix):   # two windows in a row: one alone spills into the next strand
                for k in run[i:i + 6]:
                    residues[k].ss = "H"
        for i in range(len(run) - 4):
            if fits(run, i, strand):
                for k in run[i:i + 5]:
                    if residues[k].ss != "H":
                        residues[k].ss = "E"
    # a strand pairs with another one: an alpha carbon of another strand within 5.5 Å (not its own neighbours); a grid finds them
    cell, grid = 5.5, {}
    strands = [k for k, r in enumerate(residues) if r.ss == "E"]
    for k in strands:
        grid.setdefault(tuple(int(c // cell) for c in residues[k].ca), []).append(k)
    for k in strands:
        x, y, z = (int(c // cell) for c in residues[k].ca)
        near = (j for dx in (-1, 0, 1) for dy in (-1, 0, 1) for dz in (-1, 0, 1) for j in grid.get((x + dx, y + dy, z + dz), ()))
        if not any((abs(j - k) > 2 or residues[j].chain != residues[k].chain) and math.dist(residues[j].ca, residues[k].ca) <= cell
                   for j in near):
            residues[k].ss = "-"
    for ss, shortest in (("H", 4), ("E", 3)):
        k = 0
        while k < len(residues):
            end = k
            while end < len(residues) and residues[end].ss == ss and residues[end].chain == residues[k].chain:
                end += 1
            if end - k < shortest:
                for j in range(k, end):
                    residues[j].ss = "-"
            k = max(end, k + 1)


def residue_map(protein, residues):
    """{index in `protein`: [indices in `residues`]}: each amino acid of the protein -> the residues of the model that stand for it.
    The chains are aligned to the protein, so a PDB entry with its own numbering, missing residues or a few substitutions still
    lines up; every chain that matches about as well as the best one counts (each copy of a homodimer)."""
    chains = {}
    for k, r in enumerate(residues):
        chains.setdefault(r.chain, []).append(k)
    found = []
    for ks in chains.values():
        pairs = []
        for op, a1, a2, b1, b2 in SequenceMatcher(None, protein, "".join(residues[k].aa for k in ks), autojunk=False).get_opcodes():
            if op == "equal" or (op == "replace" and a2 - a1 == b2 - b1):   # a substitution keeps its place
                pairs += [(a1 + d, ks[b1 + d]) for d in range(a2 - a1)]
        found.append(pairs)
    best = max(map(len, found), default=0)
    out = {}
    for pairs in found:
        if pairs and len(pairs) >= 0.9 * best:
            for i, k in pairs:
                out.setdefault(i, []).append(k)
    return out
