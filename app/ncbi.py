"""NCBI import: search the Nucleotide (GenBank, RefSeq, genomes) and Protein databases and fetch a record, through NCBI's E-utilities
(https://www.ncbi.nlm.nih.gov/books/NBK25497/): esearch finds the record numbers, esummary describes them, efetch gives the record
(a full GenBank file with its annotation for nucleotides, FASTA for proteins).

NCBI asks for at most 3 requests a second (10 with an API key, free in an NCBI account): `call` waits between requests, and waits
longer and asks again (twice) if NCBI still answers "too many requests". Fetched records are kept in `<data folder>/ncbi`, so the
same record is not downloaded twice. Only what you type (the search words, an accession, your optional API key) is sent.
The network call is `fetch` (replaced in the tests); it runs in a worker thread of the caller.

Reference genomes (File > Import > Reference Genome from NCBI) come from NCBI Datasets
(https://www.ncbi.nlm.nih.gov/datasets/docs/v2/api/): `search_genomes` lists the assemblies of an organism (or one
GCF_/GCA_ accession), `download_genome` fetches its sequences (FASTA) and annotation (GFF3) as one zip, streamed to
disk, and keeps them in `<data folder>/genomes/<accession>/`. The download is `open_url` (replaced in the tests).
"""
import json
import re
import threading
import time
import urllib.error
import urllib.request
import zipfile
from dataclasses import dataclass
from urllib.parse import quote, urlencode

import structure
from primer_db import data_dir

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
DATABASES = {"Nucleotide / genome (GenBank)": "nuccore", "Protein": "protein"}
TOOL = "openDNA"
MAX_RESULTS = 25
WARN_BASES = 20_000_000     # a record this long takes minutes: ask first
MAX_BASES = 100_000_000     # longer: not fetched here (NCBI Datasets is made for those)
DATASETS = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/"
GENOME_WARN_BASES = 200_000_000   # a bigger genome is a download of hundreds of MB: ask first
ASSEMBLY = re.compile(r"GC[AF]_\d{9}(\.\d+)?", re.I)
_gap = {"last": 0.0}
_lock = threading.Lock()


class NcbiError(Exception):
    """Something the user can read."""


@dataclass
class Record:
    uid: str
    accession: str        # with version: NC_000913.3
    title: str
    length: int           # bases or amino acids
    organism: str
    updated: str
    db: str               # nuccore | protein
    topology: str = "linear"
    source: str = ""      # sourcedb: refseq, genbank, swiss_prot...

    @property
    def uniprot(self):
        """The UniProt accession of a Swiss-Prot / TrEMBL record (the AlphaFold DB key), else ''."""
        return self.accession.split(".")[0] if self.source in ("swiss_prot", "uniprot", "tpa") and self.db == "protein" else ""


def fetch(url, timeout=60):
    """The text at `url`; NcbiError when NCBI cannot be reached or refuses."""
    try:
        return structure.fetch(url, None, None, timeout).decode("utf-8", "replace")
    except structure.StructureError as e:
        raise NcbiError(str(e)) from None


def call(endpoint, params, key="", timeout=60):
    """One E-utilities request (waiting first so the rate limit is kept); the answer text."""
    params = {**params, "tool": TOOL, **({"api_key": key} if key else {})}
    for attempt in range(3):
        with _lock:
            wait = _gap["last"] + (0.11 if key else 0.4) - time.monotonic()
            if wait > 0:
                time.sleep(wait)
            _gap["last"] = time.monotonic()
        try:
            return fetch(f"{BASE}{endpoint}?{urlencode(params)}", timeout)
        except NcbiError as e:
            if "429" not in str(e) or attempt == 2:
                raise
            time.sleep(1.5 * (attempt + 1))   # "too many requests": NCBI counts by address, so wait longer and ask again


def search(db, term, key="", limit=MAX_RESULTS):
    """The records of database `db` ("nuccore" / "protein") that match `term` (an accession, or words: E. coli K-12 complete genome)."""
    term = term.strip()
    if not term:
        raise NcbiError("Type an accession or some words to search for.")
    try:
        found = json.loads(call("esearch.fcgi", {"db": db, "term": term, "retmax": limit, "retmode": "json"}, key))
        ids = found["esearchresult"]["idlist"]
        if not ids:
            return []
        summary = json.loads(call("esummary.fcgi", {"db": db, "id": ",".join(ids), "retmode": "json"}, key))["result"]
    except (ValueError, KeyError) as e:
        raise NcbiError(f"NCBI's answer could not be read ({e!r}).") from None
    return [Record(uid=uid, accession=s.get("accessionversion") or s.get("caption", uid), title=s.get("title", ""), length=int(s.get("slen") or 0),
                   organism=s.get("organism", ""), updated=s.get("updatedate", ""), db=db, topology=s.get("topology", "linear"),
                   source=s.get("sourcedb", "")) for uid in ids if (s := summary.get(uid))]


def path_of(record):
    return data_dir() / "ncbi" / (re.sub(r"[^A-Za-z0-9_.-]", "_", record.accession) + (".gb" if record.db == "nuccore" else ".fasta"))


def get(record, key=""):
    """The file of `record`, fetched if it is not in the cache: (path, was_cached). Nucleotides: GenBank with annotation
    (efetch gbwithparts); proteins: FASTA."""
    path = path_of(record)
    if path.exists() and path.stat().st_size:
        return path, True
    if record.db == "nuccore" and record.length > MAX_BASES:
        raise NcbiError(f"{record.accession} has {record.length:,} bases: too long to fetch here (use NCBI Datasets).")
    rettype = "gbwithparts" if record.db == "nuccore" else "fasta"
    text = call("efetch.fcgi", {"db": record.db, "id": record.uid, "rettype": rettype, "retmode": "text"}, key, timeout=600)
    if not text.strip() or text.lstrip().startswith(("<?xml", "Error")) or ("LOCUS" not in text[:2000] if record.db == "nuccore" else not text.startswith(">")):
        raise NcbiError(f"NCBI did not return {record.accession} ({text.strip()[:120] or 'empty answer'}).")
    path.parent.mkdir(exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path, False


def protein_of(path):
    """(name, amino acids) of the first record of a FASTA file."""
    lines = path.read_text(encoding="utf-8").splitlines()
    body = []
    for line in lines[1:]:
        if line.startswith(">"):
            break
        body.append(line.strip())
    return lines[0][1:].split(" ")[0], "".join(body).upper()


# ---- reference genomes (NCBI Datasets) ---------------------------------------------------------------------------------
@dataclass
class Assembly:
    accession: str        # GCF_000005845.2 (RefSeq) or GCA_... (GenBank)
    organism: str
    name: str             # ASM584v2
    level: str            # Complete Genome, Chromosome, Scaffold, Contig
    length: int           # bases
    sequences: int        # chromosomes, plasmids, contigs...
    released: str
    category: str = ""    # "reference genome" when NCBI picked it for its species


def _datasets_json(url, key=""):
    try:
        return json.loads(fetch(url + (("&" if "?" in url else "?") + urlencode({"api_key": key}) if key else "")))
    except ValueError as e:
        raise NcbiError(f"NCBI's answer could not be read ({e!r}).") from None
    except NcbiError as e:
        if "404" in str(e) or "400" in str(e):
            return {}   # unknown organism / accession: nothing found
        raise


def search_genomes(term, reference_only=True, key="", limit=MAX_RESULTS):
    """Assemblies of an organism (name or taxonomy number: "Escherichia coli", 562), or the one with a GCF_/GCA_
    accession. `reference_only`: only NCBI's reference genome of each species."""
    term = term.strip()
    if not term:
        raise NcbiError("Type an organism (Escherichia coli) or an assembly accession (GCF_000005845.2).")
    if ASSEMBLY.fullmatch(term):
        url = f"{DATASETS}accession/{quote(term.upper())}/dataset_report"
    else:
        url = f"{DATASETS}taxon/{quote(term)}/dataset_report?" + urlencode(
            {"filters.reference_only": "true" if reference_only else "false", "page_size": limit})
    reports = _datasets_json(url, key).get("reports", [])
    out = []
    for r in reports:
        info, stats = r.get("assembly_info", {}), r.get("assembly_stats", {})
        out.append(Assembly(accession=r.get("accession", ""), organism=r.get("organism", {}).get("organism_name", ""),
                            name=info.get("assembly_name", ""), level=info.get("assembly_level", ""),
                            length=int(stats.get("total_sequence_length") or 0),
                            sequences=int(stats.get("number_of_component_sequences") or stats.get("number_of_scaffolds") or 0),
                            released=info.get("release_date", ""), category=info.get("refseq_category", "")))
    return out


def genome_folder(assembly):
    return data_dir() / "genomes" / re.sub(r"[^A-Za-z0-9_.-]", "_", assembly.accession)


def open_url(url, headers, timeout=120):
    """A file-like answer for `url` (streamed, not read at once); NcbiError when it cannot be reached."""
    request = urllib.request.Request(url, headers={"User-Agent": structure.USER_AGENT, **headers})
    try:
        return urllib.request.urlopen(request, timeout=timeout)
    except urllib.error.HTTPError as e:
        raise NcbiError(f"NCBI answered {e.code} {e.reason}.") from None
    except (urllib.error.URLError, TimeoutError, OSError) as e:
        raise NcbiError(f"Could not reach NCBI: {getattr(e, 'reason', e)}") from None


def download_genome(assembly, key="", progress=None):
    """(FASTA path, GFF3 path or None, was_cached) of `assembly`, downloaded once into its genome folder.
    `progress(bytes so far)` is called while the zip streams in."""
    folder = genome_folder(assembly)
    fasta, gff = folder / f"{folder.name}.fna", folder / f"{folder.name}.gff"
    if fasta.exists() and fasta.stat().st_size:
        return fasta, (gff if gff.exists() else None), True
    folder.mkdir(parents=True, exist_ok=True)
    url = f"{DATASETS}accession/{quote(assembly.accession)}/download?" + urlencode(
        [("include_annotation_type", "GENOME_FASTA"), ("include_annotation_type", "GENOME_GFF")])
    archive = folder / "download.zip.part"
    got = 0
    try:
        with open_url(url, {"api-key": key} if key else {}, timeout=600) as reply, open(archive, "wb") as out:
            while chunk := reply.read(1 << 20):
                out.write(chunk)
                got += len(chunk)
                if progress:
                    progress(got)
        # Only the sequence and annotation files, written under our own names (never a path taken from the zip).
        with zipfile.ZipFile(archive) as z:
            members = z.namelist()
            fna = [m for m in members if m.endswith((".fna", ".fa", ".fasta"))]
            gffs = [m for m in members if m.endswith((".gff", ".gff3"))]
            if not fna:
                raise NcbiError(f"NCBI sent no sequence for {assembly.accession}.")
            with open(fasta, "wb") as dst:   # several FASTA files (one per chromosome): one file
                for m in fna:
                    with z.open(m) as src:
                        last = b"\n"
                        while block := src.read(1 << 20):
                            dst.write(block)
                            last = block[-1:]
                        if last != b"\n":
                            dst.write(b"\n")   # a member without a final newline would glue the next header to its last line
            if gffs:
                with z.open(gffs[0]) as src, open(gff, "wb") as dst:
                    while block := src.read(1 << 20):
                        dst.write(block)
    except zipfile.BadZipFile:
        fasta.unlink(missing_ok=True)
        raise NcbiError(f"The download of {assembly.accession} is damaged: try again.") from None
    except BaseException:
        fasta.unlink(missing_ok=True)   # no half file that would look downloaded next time
        raise
    finally:
        archive.unlink(missing_ok=True)
    return fasta, (gff if gff.exists() else None), False
