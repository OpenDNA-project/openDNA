"""Build one collection of the bundled parts library (app/parts/ncbi_vectors.tsv) from public NCBI GenBank records of
well-known vectors (NCBI says most of its database content is free to use, but a submitter may claim rights in a record: see https://www.ncbi.nlm.nih.gov/home/about/policies/ and the records listed below). Run once; the output is committed:

    python scripts/build_parts.py

Every annotated feature of a useful kind becomes a part: name (label > gene > product >
short note), type, sequence in the part's own orientation, source accession. Identical
sequences are kept once. A record says what a part is, never how strong it is, so the
organism / strength / units columns of the collection format stay empty here.
"""
import csv
import re
import sys
import tempfile
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "app"))
sys.stdout.reconfigure(encoding="utf-8", errors="replace")   # part names carry Greek letters (lacZα)

from opendna_core import Document  # noqa: E402

from editor import reverse_complement  # noqa: E402

ACCESSIONS = {
    "L09137": "pUC19", "L08752": "pUC18", "J01749": "pBR322", "X06403": "pACYC184", "X01654": "pSC101",
    "AY048742": "pKD3", "AY048743": "pKD4", "AY048744": "pKD13", "AY048746": "pKD46",
    "LC760188": "pBAD33", "X81838": "pBAD18", "M22744": "pTrc99A", "U23751": "pBBR1MCS-2",
    "JX560327": "pSEVA221", "CP035489": "pZE21", "MT127274": "pLacI", "MW147276": "pmCherry",
    "AY678264": "mRFP", "KJ541673": "pCellFree_G09 (sfGFP)", "KM099231": "pCas9", "FJ437239": "pK18mobsacB",
    "AF047518": "pEX18Gm", "U55762": "pEGFP-N1", "X65316": "pGEM-luc", "AF234296": "pCAMBIA-1300",
}
KINDS = {"CDS", "promoter", "terminator", "rep_origin", "RBS", "protein_bind", "oriT", "regulatory",
         "misc_feature", "polyA_signal", "enhancer", "LTR", "misc_RNA", "ncRNA", "sig_peptide"}
NAME_KEYS = ("label", "gene", "product", "standard_name", "note")
MIN_LEN, MAX_LEN = 15, 6000
# Names that describe where a stretch came from, not what it is.
JUNK = ("derived from", "fragment", "part of", "from ", "region", "sequence", "hypothetical", "priming site",
        "paraslipin", "dyad", "orf", "_yqr")
# The short names people use for common markers (as on most plasmid maps).
ALIASES = {
    "bla": "AmpR", "beta-lactamase": "AmpR", "cat": "CmR", "chloramphenicol acetyl transferase": "CmR",
    "neo": "KanR", "nptii": "KanR", "tn5 neomycin phosphotransferase": "KanR", "neomycin phosphotransferase": "KanR",
    "aminoglycoside phosphotransferase": "KanR", "aacc1": "GmR", "hygromycin phosphotransferase": "HygR",
    "tet": "TetR (tetA)", "egfp": "EGFP", "mcherry": "mCherry", "monomeric red fluorescent protein": "mRFP1",
    "luc": "luciferase (luc)", "araC protein": "araC", "rop protein": "rop", "cas9 endonuclease": "Cas9",
    "lacZ": "lacZα", "lacZalpha": "lacZα", "ccdb": "ccdB",
}


def kind_of(kind, name):
    """GenBank 'regulatory' / 'misc_feature' -> promoter / terminator when the name says so."""
    low = name.lower()
    if "promoter" in low or low in ("parab", "pc"):
        return "promoter"
    if "terminator" in low or low in ("t0", "t1", "t2", "tl3", "rrnb"):
        return "terminator"
    if "origin" in low or low.endswith("_ori") or low == "orit":
        return "rep_origin" if "orit" not in low else "oriT"
    return kind
HEADER = [
    "#title: Cloning vectors (NCBI records)",
    "#reference: Annotated features of public NCBI GenBank records of well-known cloning vectors "
    "(scripts/build_parts.py lists the records).",
    "#note: Names, short feature sequences and their accession numbers; no strength was measured. "
    "NCBI says most of its database content is free to use, but a submitter may claim rights in a record.",
]
URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gbwithparts&retmode=text&id="

# ---- curation (docs/PARTS_LIBRARY_METHOD.md) ---------------------------------------------------
# A GenBank submitter annotates their own construct, not a parts catalogue: they label primer sites,
# polylinkers, and whichever slice of an FRT site survived their cloning. Everything below decides
# what of that is a reusable part. Rebuilding runs these rules again, so the junk cannot come back.

# Names that are never a part, whatever the record calls them.
NEVER = re.compile(r"primer|polylinker|multiple cloning|\bMCS\b|\bdistal\b|-nt of|effector site"
                   r"|cleavage site|^SITS$|^operator |stuffer|\bspacer\b|\bfragment\b|\bpartial\b"
                   r"|\bregion\b|derived from", re.IGNORECASE)
NEVER_TYPES = {"sig_peptide"}   # a piece of a CDS the library already holds whole

# `misc_feature` and `regulatory` say nothing about what a part is, so a row of either is dropped
# unless it is named here, which gives it a real type and the name people use for it.
RETYPE = {
    "lac operator": ("protein_bind", "lac operator (lacO)"),
    "natural FRT site": ("misc_recomb", "FRT"),
    "CaMV 3'UTR (polyA signal)": ("polyA_signal", "CaMV polyA signal"),
    "yeast centromere": ("centromere", "CEN (yeast centromere)"),
}

# One name per biological part: these all mean the same thing, and only the longest row survives.
CANONICAL = {
    "p15A Origin": "p15A ori", "p15A origin of replication": "p15A ori", "p15A_ori": "p15A ori",
    "pBR322 origin": "ColE1 / pMB1 ori", "pBR322 origin of replication": "ColE1 / pMB1 ori",
    "pMB1 replication origin": "ColE1 / pMB1 ori",
    "R6K gamma replication origin": "R6K gamma ori",
    "f1_ori": "f1 ori", "M13 origin": "M13 ori", "oriR101": "pSC101 oriR",
    "rrnB": "rrnB T1T2 terminator", "rrnB transcription terminator": "rrnB T1T2 terminator",
    "5S RNA, T1, T2, rrnB": "rrnB T1T2 terminator",
    "lambda terminator": "lambda tL3 terminator", "tL3": "lambda tL3 terminator",
    "Lac promoter": "lac promoter", "E. coli lac operon promoter": "lac promoter",
    "PlacZ; lacZ promoter": "lac promoter",
    "Pc": "Pc promoter", "ParaB": "araBAD promoter", "Bla promoter": "bla (AmpR) promoter",
    "oriT": "RP4 oriT", "rep": "pBBR1 rep", "repA101 protein": "pSC101 repA101",
    "trfA": "RK2 trfA", "blasticidin S deaminase": "BsdR (blasticidin)",
}

# Vector labels that name a position in one plasmid rather than a part. Which terminator each of
# these actually is cannot be told from the record, and a part of the right name with the wrong
# sequence is worse than no part, so they go.
AMBIGUOUS = {"T0", "T1", "T2"}

# Coding parts that app/parts/engineering_proteins.tsv now holds, verified against UniProt with an
# organism and a length gate. That collection is the authority for a CDS; these rows are its
# worse-documented twins (defect D5) and are dropped here rather than left to collide at load time.
SUPERSEDED = {
    "AmpR", "Cas9", "CmR", "GmR", "HygR", "KanR", "TetR (tetA)", "araC", "bet protein", "ccdB",
    "exonuclease", "gam protein", "luciferase (luc)", "rop", "sacB",
    "neomycin phosphoribosyl transferase",   # a KanR variant: the protein matcher finds it from nptII
}

# The shortest length worth searching for, per type. One global minimum cannot work: a phage
# promoter really is 19 bp and highly specific, while an 18 bp "operator" matches by chance.
FLOORS = {"CDS": 150, "promoter": 19, "terminator": 30, "rep_origin": 200, "oriT": 50,
          "protein_bind": 25, "misc_recomb": 25, "polyA_signal": 50, "centromere": 100}


def curate(parts):
    """(name, kind, sequence, source) rows in, curated rows out: junk dropped, types and names made
    canonical, one row per part (the longest), nothing under its type's floor."""
    kept = {}
    for name, kind, seq, source in parts:
        if kind in NEVER_TYPES or NEVER.search(name) or name in SUPERSEDED or name in AMBIGUOUS:
            continue
        if kind in ("misc_feature", "regulatory"):
            if name not in RETYPE:
                continue
            kind, name = RETYPE[name]
        name = CANONICAL.get(name, name)
        if len(seq) < FLOORS.get(kind, 50):
            continue
        # One part, one row — but which row? The matcher needs 95 % identity across the whole part with
        # no gaps, so a long part is the harder one to find: a 426 bp rrnB does not match a plasmid
        # carrying 260 bp of it. For a CDS the longest row is the complete gene and the right answer;
        # for everything else the records annotate different extents of one element, and the shortest
        # surviving row is its conserved core, which is found in more constructs.
        if name in kept:
            better = len(seq) > len(kept[name][2]) if kind == "CDS" else len(seq) < len(kept[name][2])
            if not better:
                continue
        kept[name] = (name, kind, seq, source)
    return sorted(kept.values(), key=lambda r: (r[1], r[0].lower()))


def name_of(f):
    for key in NAME_KEYS:
        v = next((v for k, v in f.qualifiers if k == key and v), None)
        if v and (key != "note" or len(v) <= 40):  # long notes are prose, not names
            return v.strip()
    return None


OUT = ROOT / "app" / "parts" / "ncbi_vectors.tsv"


def write(parts):
    OUT.parent.mkdir(exist_ok=True)
    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        for line in HEADER:
            fh.write(line + "\n")
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["name", "type", "sequence", "organism", "strength", "units", "source"])
        w.writerows((name, kind, seq, "", "", "", source) for name, kind, seq, source in parts)
    print(f"{len(parts)} parts written to {OUT}")


def recurate():
    """Run the curation rules over the collection already on disk, without fetching anything.
    The result is what a full rebuild would write, so the two cannot drift apart."""
    rows = [r for r in csv.DictReader([line for line in open(OUT, encoding="utf-8")
                                       if not line.startswith("#")], delimiter="\t")]
    before = [(r["name"], r["type"], r["sequence"], r["source"]) for r in rows]
    after = curate(before)
    # A renamed or merged part is still there under its canonical name: report it apart from a real drop.
    kept_seqs = {r[2]: r[0] for r in after}
    for name, _, seq, _ in sorted(before, key=lambda r: r[0].lower()):
        if seq in kept_seqs:
            if kept_seqs[seq] != name:
                print(f"  renamed  {name}  ->  {kept_seqs[seq]}")
        else:
            print(f"  DROPPED  {name}")
    print(f"{len(before)} -> {len(after)} parts")
    write(after)


def main():
    if "--recurate" in sys.argv:
        recurate()
        return
    parts, seen = [], set()
    with tempfile.TemporaryDirectory() as tmp:
        for acc, vector in ACCESSIONS.items():
            path = Path(tmp) / f"{acc}.gb"
            path.write_bytes(urllib.request.urlopen(URL + acc, timeout=60).read())
            time.sleep(0.4)  # NCBI asks for at most 3 requests per second
            doc = Document(str(path))
            for f in doc.features():
                name = name_of(f)
                if f.kind not in KINDS or not name or len(f.spans) != 1 or any(j in name.lower() for j in JUNK):
                    continue
                s, e = f.spans[0]
                if not MIN_LEN <= e - s <= (MAX_LEN if f.kind == "CDS" else 1500):
                    continue
                name = ALIASES.get(name.lower(), ALIASES.get(name, name))
                kind = kind_of(f.kind, name)
                seq = doc.sequence(s, e).upper()
                if f.strand < 0:
                    seq = reverse_complement(seq)
                if seq in seen or set(seq) - set("ACGT"):
                    continue
                seen.add(seq)
                parts.append((name, kind, seq, f"{vector} ({acc})"))
            print(f"{acc} {vector}: {len(parts)} parts so far")
    print(f"{len(parts)} features annotated in the records")
    write(curate(parts))


if __name__ == "__main__":
    main()
