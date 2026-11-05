"""Build app/parts/engineering_proteins.tsv: the coding parts a plasmid designer reaches for
(markers, counter-selection, recombinases, CRISPR nucleases, polymerases, reporters, regulators).

Each part is found in UniProtKB/Swiss-Prot (reviewed entries only, CC BY 4.0), which names the
canonical protein and cross-references the nucleotide record holding its coding sequence; the DNA
itself is then taken from that record at NCBI (public domain). So every row traces to a UniProt
accession and a GenBank accession, and nothing is reconstructed by hand.

    python scripts/build_engineering_parts.py            # writes the collection
    python scripts/build_engineering_parts.py --report   # only prints what it would write

Swiss-Prot holds natural proteins: engineered variants (sfGFP, mScarlet, dCas9, SB100X, piggyBac)
are deliberately absent. The library finds them anyway, because the detector also matches by protein
identity, and a variant of a part stays close to the parent it was made from.
"""
import csv
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
UNIPROT = "https://rest.uniprot.org/uniprotkb"
E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
PAUSE = 0.4   # NCBI asks for at most 3 requests a second

HEADER = [
    "#title: Engineering parts (UniProt + GenBank)",
    "#reference: UniProt Consortium. UniProt: the Universal Protein Knowledgebase in 2025. "
    "Nucleic Acids Res 53, D609-D617 (2025). Coding sequences from the GenBank records UniProt cites.",
    "#doi: 10.1093/nar/gkae1010",
    "#note: Swiss-Prot (reviewed) entries, CC BY 4.0; the DNA is the coding sequence of the GenBank "
    "record each entry cites (NCBI public domain). No strength was measured for these parts.",
]

# (name on the map, feature type, Swiss-Prot query). The query names the protein and its organism,
# so the reviewed entry that comes back is the canonical one, not a homologue in another species.
WANTED = [
    # ---- antibiotic resistance markers ------------------------------------------------------
    ("AmpR (bla, TEM-1)", "CDS", 'protein_name:"Beta-lactamase TEM" AND organism_id:562'),
    ("KanR (nptII, Tn5)", "CDS", 'gene:nptII AND protein_name:"aminoglycoside"'),
    ("TetA (tetracycline efflux)", "CDS", 'protein_name:"Tetracycline resistance protein, class B"'),
    ("TetR (tet repressor)", "CDS", 'protein_name:"Tetracycline repressor protein class B"'),
    ("HygR (hph)", "CDS", 'protein_name:"Hygromycin-B 4-O-kinase"'),
    ("SpecR / StrR (aadA)", "CDS", 'protein_name:"Streptomycin 3\'\'-adenylyltransferase"'),
    ("GmR (aacC1)", "CDS", 'protein_name:"Gentamicin 3-N-acetyltransferase"'),
    ("BleoR / ZeoR (ble)", "CDS", 'protein_name:"Bleomycin resistance protein"'),
    ("PuroR (pac)", "CDS", 'protein_name:"Puromycin N-acetyltransferase"'),
    ("ErmC (MLS resistance)", "CDS", 'gene:ermC AND protein_name:"rRNA adenine N-6-methyltransferase"'),
    # ---- counter-selection markers ----------------------------------------------------------
    ("SacB (levansucrase)", "CDS", 'gene:sacB AND organism_name:"Bacillus subtilis"'),
    ("CcdB (toxin)", "CDS", 'gene:ccdB AND organism_id:562'),
    ("CcdA (antitoxin)", "CDS", 'gene:ccdA AND organism_id:562'),
    ("RpsL (streptomycin sensitivity)", "CDS", 'gene:rpsL AND organism_id:83333'),
    ("URA3", "CDS", 'gene:URA3 AND organism_name:"Saccharomyces cerevisiae"'),
    ("ThyA", "CDS", 'gene:thyA AND organism_id:83333'),
    ("GalK", "CDS", 'gene:galK AND organism_id:83333'),
    ("MazF (toxin)", "CDS", 'gene:mazF AND organism_id:83333'),
    ("PheS", "CDS", 'gene:pheS AND organism_id:83333'),
    # ---- recombinases, integrases, transposases ---------------------------------------------
    ("Cre recombinase", "CDS", 'protein_name:"Recombinase cre"'),
    ("Flp recombinase", "CDS", 'protein_name:"Site-specific recombinase Flp"'),
    ("phiC31 integrase", "CDS", 'gene:int AND organism_name:"Streptomyces phage phiC31"'),
    ("Bxb1 integrase", "CDS", 'protein_name:"Serine integrase" AND organism_name:"Mycobacterium phage Bxb1"'),
    ("Tn3 transposase (tnpA)", "CDS", 'accession:P03008'),
    ("Tn5 transposase", "CDS", 'protein_name:"Transposase" AND organism_name:"Transposon Tn5"'),
    ("Tn7 transposase (tnsB)", "CDS", 'gene:tnsB AND organism_name:"Transposon Tn7"'),
    ("Lambda Red exo", "CDS", 'gene:exo AND organism_name:"Escherichia phage lambda"'),
    ("Lambda Red beta", "CDS", 'gene:bet AND organism_name:"Escherichia phage lambda"'),
    ("Lambda Red gam", "CDS", 'gene:gam AND organism_name:"Escherichia phage lambda"'),
    ("Lambda integrase", "CDS", 'gene:int AND organism_name:"Escherichia phage lambda"'),
    # ---- CRISPR ------------------------------------------------------------------------------
    ("SpCas9", "CDS", 'accession:Q99ZW2'),
    ("SaCas9", "CDS", 'accession:J7RUA5'),
    ("NmCas9", "CDS", 'accession:C9X1G5'),
    ("StCas9", "CDS", 'accession:G3ECR1'),
    ("Cas12a (Cpf1, FnCas12a)", "CDS", 'accession:A0Q7Q2'),
    ("Cas3", "CDS", 'gene:cas3 AND organism_id:83333'),
    # ---- RNA polymerases ----------------------------------------------------------------------
    ("T7 RNA polymerase", "CDS", 'protein_name:"DNA-directed RNA polymerase" AND organism_name:"Escherichia phage T7"'),
    ("T3 RNA polymerase", "CDS", 'protein_name:"DNA-directed RNA polymerase" AND organism_name:"Escherichia phage T3"'),
    ("SP6 RNA polymerase (alt)", "CDS", 'protein_name:"DNA-directed RNA polymerase" AND organism_name:"Enterobacteria phage SP6"'),
    ("SP6 RNA polymerase", "CDS", 'protein_name:"DNA-directed RNA polymerase" AND organism_name:"Salmonella phage SP6"'),
    # ---- reporters and tags --------------------------------------------------------------------
    ("GFP (avGFP)", "CDS", 'accession:P42212'),
    ("DsRed", "CDS", 'accession:P83689'),
    ("Firefly luciferase (luc)", "CDS", 'accession:P08659'),
    ("Renilla luciferase", "CDS", 'accession:P27652'),
    ("LacZ (beta-galactosidase)", "CDS", 'gene:lacZ AND organism_id:83333'),
    ("GusA (beta-glucuronidase)", "CDS", 'gene:uidA AND organism_id:83333'),
    ("CmR (cat, Tn9)", "CDS", 'protein_name:"Chloramphenicol acetyltransferase" AND organism_name:"Escherichia coli"'),
    # ---- regulators used as switches -------------------------------------------------------------
    ("LacI", "CDS", 'gene:lacI AND organism_id:83333'),
    ("AraC", "CDS", 'gene:araC AND organism_id:83333'),
    ("RhaS", "CDS", 'gene:rhaS AND organism_id:83333'),
    ("CI repressor (lambda)", "CDS", 'gene:cI AND organism_name:"Escherichia phage lambda"'),
    ("LuxR", "CDS", 'gene:luxR AND organism_name:"Aliivibrio fischeri"'),
    ("LuxI", "CDS", 'gene:luxI AND organism_name:"Aliivibrio fischeri"'),
    # ---- assembly and maturation ------------------------------------------------------------------
    ("Rop / Rom", "CDS", 'gene:rop AND organism_id:562'),
    ("MBP (malE)", "CDS", 'gene:malE AND organism_id:83333'),
    ("SUMO (Smt3)", "CDS", 'gene:SMT3 AND organism_name:"Saccharomyces cerevisiae"'),
    ("TEV protease", "CDS", 'protein_name:"Nuclear inclusion protein A" AND organism_name:"Tobacco etch virus"'),
]


def get(url, tries=3):
    for attempt in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=90) as fh:
                return fh.read().decode("utf-8", "replace")
        except Exception as e:
            if attempt == tries - 1:
                raise
            print(f"    retry after {e}")
            time.sleep(2)


def best_entry(query):
    """The best hit: (accession, protein name, organism, [(embl acc, protein id)]). Reviewed
    (Swiss-Prot) entries are preferred; a few phage and transposon proteins are only in TrEMBL,
    which carries the same CC BY 4.0 terms. The JSON form is used because only it names the
    protein id of each EMBL cross-reference, which is what points at the coding sequence."""
    for clause in (f"({query}) AND (reviewed:true)", query):
        url = f"{UNIPROT}/search?query={urllib.parse.quote(clause)}&format=json&size=1"
        results = json.loads(get(url)).get("results", [])
        if not results:
            continue
        r = results[0]
        xrefs = []
        for x in r.get("uniProtKBCrossReferences", []):
            if x.get("database") != "EMBL":
                continue
            props = {p["key"]: p["value"] for p in x.get("properties", [])}
            pid, kind = props.get("ProteinId"), props.get("MoleculeType", "")
            if pid and pid != "-" and "RNA" not in kind:   # genomic DNA, not an mRNA-only entry
                xrefs.append((x["id"], pid))
        name = r.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", "")
        return (r["primaryAccession"], name, r.get("organism", {}).get("scientificName", ""), xrefs,
                r.get("sequence", {}).get("length", 0))
    return None


def coding_dna(protein_id):
    """The coding sequence behind a GenBank protein accession, through its /coded_by."""
    gp = get(f"{E}/efetch.fcgi?db=protein&id={protein_id}&rettype=gp&retmode=text")
    time.sleep(PAUSE)
    m = re.search(r'/coded_by="(complement\()?([A-Z_0-9.]+):(<?\d+)\.\.(>?\d+)', gp)
    if not m:
        return None, None
    rev, acc, start, stop = m.group(1), m.group(2), m.group(3).lstrip("<"), m.group(4).lstrip(">")
    fa = get(f"{E}/efetch.fcgi?db=nuccore&id={acc}&rettype=fasta&retmode=text"
             f"&seq_start={start}&seq_stop={stop}&strand={2 if rev else 1}")
    time.sleep(PAUSE)
    seq = "".join(line.strip() for line in fa.splitlines() if not line.startswith(">")).upper()
    return (seq if seq and not set(seq) - set("ACGTN") else None), acc


def main():
    report_only = "--report" in sys.argv
    rows, missing = [], []
    for name, kind, query in WANTED:
        try:
            found = best_entry(query)
        except Exception as e:
            missing.append((name, f"UniProt error: {e}"))
            continue
        if not found:
            missing.append((name, "no reviewed entry"))
            continue
        acc, protein, organism, xrefs, aa = found
        # The length gate: a coding sequence must be exactly the protein plus a stop codon. Cross-references
        # to a partial CDS or to the wrong feature of a record come back short or long, and are refused —
        # a part of the right name with the wrong sequence would be worse than no part at all.
        seq = source = None
        for _, protein_id in xrefs[:6]:
            try:
                dna, record = coding_dna(protein_id)
            except Exception as e:
                print(f"    {name}: {e}")
                continue
            if dna and aa and len(dna) in (3 * aa + 3, 3 * aa):
                seq, source = dna, f"{record} ({protein_id}), UniProt {acc}"
                break
        if not seq:
            missing.append((name, f"{acc}: no cross-reference gave the whole {aa} aa coding sequence"))
            continue
        if len(seq) < 150:          # see FLOORS in build_parts.py: a coding part shorter than this is a tag
            missing.append((name, f"{len(seq)} bp is under the 150 bp floor for a CDS"))
            continue
        organism = organism.split(" (")[0]
        rows.append((name, kind, seq, organism, "", "", source))
        print(f"  {name:36s} {len(seq):>6,} bp  {acc}  {protein[:44]}")
    print(f"\n{len(rows)} parts, {len(missing)} not found")
    for name, why in missing:
        print(f"  MISSING  {name}: {why}")
    if report_only:
        return
    out = ROOT / "app" / "parts" / "engineering_proteins.tsv"
    out.parent.mkdir(exist_ok=True)
    with open(out, "w", newline="", encoding="utf-8") as fh:
        for line in HEADER:
            fh.write(line + "\n")
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["name", "type", "sequence", "organism", "strength", "units", "source"])
        w.writerows(rows)
    print(f"written to {out}")


if __name__ == "__main__":
    main()
