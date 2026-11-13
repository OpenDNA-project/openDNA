"""Online sources, with the network replaced by fakes: the NCBI client (search, records, its cache, its rate limit, errors) and
the 3D structure sources (AlphaFold DB, PDB, ESMFold, a local command; the cache and Auto's order), the viewer page's script
safety, PDB entry lists and the model read from a PDB text."""
import json
import math
import os
import sys
import tempfile
import time
from pathlib import Path
from urllib.parse import parse_qs, urlparse

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ["OPENDNA_DATA_DIR"] = tempfile.mkdtemp(prefix="opendna_test_")
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "app"))

import pytest

import ncbi
import structure
from ncbi import NcbiError
from structure import StructureError

GENBANK = (ROOT / "demo" / "demo_vector.gb").read_text(encoding="utf-8")
SUMMARY = {"uids": ["11"], "11": {"uid": "11", "caption": "U1", "title": "A made-up plasmid", "slen": 1834, "organism": "synthetic construct",
                                   "updatedate": "2026/01/02", "accessionversion": "U1.1", "topology": "circular", "sourcedb": "genbank"}}
PROTEIN_SUMMARY = {"uids": ["22"], "22": {"uid": "22", "caption": "P12345", "title": "RecName: Full=Tiny; AltName: Full=x", "slen": 6,
                                          "organism": "Nowhere", "updatedate": "2026/01/02", "accessionversion": "P12345.1",
                                          "moltype": "aa", "sourcedb": "swiss_prot"}}


PROTEIN = "TTCCPSIVARSNFNVCRLPGTPEALCATYTGCIIIPGATCPGDYAN"


def atom(i, b):
    return f"ATOM  {i:5d}  CA  THR A{i:4d}      -8.418  -7.381   3.408  1.00{b:6.2f}           C  "


def line(serial, res, chain, num, xyz, b=50.0, name="CA", record="ATOM  ", element="C"):
    """One atom record, in the PDB's columns."""
    return (f"{record}{serial:5d}  {name:<3s} {res} {chain}{num:4d}    {xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00{b:6.2f}"
            f"          {element:>2s}")


def pdb(b):
    return "\n".join(atom(i, b) for i in range(1, 4)) + "\n"


class Net:
    """A fake of structure.fetch: answers by address, and remembers what was asked."""

    def __init__(self, esm=0.5, pdb_hits=True):
        self.calls, self.esm, self.pdb_hits = [], esm, pdb_hits

    def __call__(self, url, data=None, headers=None, timeout=30):
        self.calls.append(url)
        if "alphafold.ebi.ac.uk/api" in url:
            if url.endswith("/BADACC1"):
                raise StructureError("alphafold.ebi.ac.uk answered 404 Not Found.")
            return json.dumps([{"modelEntityId": "AF-P1-F1", "globalMetricValue": 91.5, "pdbUrl": "https://alphafold.ebi.ac.uk/files/x.pdb"}]).encode()
        if url.endswith("x.pdb"):
            return pdb(91.5).encode()
        if "search.rcsb.org" in url and json.loads(data)["query"]["type"] == "group":   # the entries mapped to a UniProt accession
            nodes = json.loads(data)["query"]["nodes"]
            assert nodes[0]["parameters"]["value"] == "P12345" and nodes[1]["parameters"]["value"] == "UniProt"
            return json.dumps({"result_set": [{"identifier": "3DEF"}, {"identifier": "1ABC"}]}).encode()
        if "data.rcsb.org/graphql" in url and "entity_poly" in json.loads(data)["query"]:   # the linked entries, with their chains
            asked = json.loads(data)["query"]
            assert "1ABC" in asked and "3DEF" in asked and asked.count("1ABC") == 1, "each entry asked about once"

            def entity(seq, chains, organism):
                return {"entity_poly": {"pdbx_seq_one_letter_code_can": seq}, "rcsb_entity_source_organism": [{"ncbi_scientific_name": organism}],
                        "rcsb_polymer_entity_container_identifiers": {"auth_asym_ids": chains}}
            return json.dumps({"data": {"entries": [
                {"rcsb_id": "1ABC", "struct": {"title": "Unrelated"}, "exptl": [{"method": "X-RAY DIFFRACTION"}],
                 "rcsb_entry_info": {"resolution_combined": [2.0]}, "polymer_entities": [entity("GGGGGGGG", ["A"], "Mus musculus")]},
                {"rcsb_id": "3DEF", "struct": {"title": "The protein in a complex"}, "exptl": [{"method": "ELECTRON MICROSCOPY"}],
                 "rcsb_entry_info": {"resolution_combined": [3.1]},
                 "polymer_entities": [entity("WWWWWWWWWW", ["X"], "Homo sapiens"), entity(PROTEIN[10:], ["B", "D"], "Crambe hispanica")]},
            ]}}).encode()
        if "search.rcsb.org" in url:
            body = json.loads(data)
            assert body["query"]["parameters"]["value"] == PROTEIN and body["query"]["parameters"]["identity_cutoff"] in (0.9, 0.5)
            match = {"sequence_identity": 0.95, "query_beg": 1, "query_end": len(PROTEIN) // 2}
            found = [{"identifier": "1ABC_1", "score": 1.0, "services": [{"nodes": [{"match_context": [match]}]}]},
                     {"identifier": "2XYZ_2", "score": 0.8},             # no match context: the score stands for the identity
                     {"identifier": "1ABC_2", "score": 0.7},             # the same entry again: listed once
                     {"identifier": "../etc_1", "score": 0.6}]           # not a PDB code: left out
            return json.dumps({"result_set": found}).encode() if self.pdb_hits else b""
        if "data.rcsb.org/graphql" in url:
            assert "1ABC" in json.loads(data)["query"]
            return json.dumps({"data": {
                "entries": [{"rcsb_id": "1ABC", "struct": {"title": "Crambin"}, "exptl": [{"method": "X-RAY DIFFRACTION"}],
                             "rcsb_entry_info": {"resolution_combined": [1.5]}}],
                "polymer_entities": [{"rcsb_id": "1ABC_1", "rcsb_entity_source_organism": [{"ncbi_scientific_name": "Crambe hispanica"}],
                                      "rcsb_polymer_entity_container_identifiers": {"auth_asym_ids": ["A", "B"]}}]}}).encode()
        if "files.rcsb.org" in url:
            return pdb(12.0).encode()
        if "esmatlas" in url:
            assert data.decode() == PROTEIN
            return pdb(self.esm).encode()
        raise AssertionError(url)


def clean_cache():
    for f in structure.cache_dir().iterdir():
        f.unlink()


@pytest.fixture
def net(monkeypatch):
    fake = Net()
    monkeypatch.setattr(structure, "fetch", fake)
    clean_cache()   # a clean cache for each test
    return fake


class NcbiFake:
    """A fake of ncbi.fetch: answers by endpoint and remembers the requests."""

    def __init__(self):
        self.requests = []

    def __call__(self, url, timeout=60):
        parts = urlparse(url)
        q = {k: v[0] for k, v in parse_qs(parts.query).items()}
        self.requests.append((parts.path.split("/")[-1], q))
        endpoint = parts.path.split("/")[-1]
        if endpoint == "esearch.fcgi":
            if q["term"] == "nothing":
                return json.dumps({"esearchresult": {"idlist": []}})
            if q["term"] == "broken":
              raise NcbiError("eutils.ncbi.nlm.nih.gov answered 503 Service Unavailable.")
            return json.dumps({"esearchresult": {"idlist": ["22" if q["db"] == "protein" else "11"]}})
        if endpoint == "esummary.fcgi":
            return json.dumps({"result": PROTEIN_SUMMARY if q["db"] == "protein" else SUMMARY})
        if endpoint == "efetch.fcgi":
            return ">P12345.1 Tiny [Nowhere]\nMKWHAA\n" if q["db"] == "protein" else GENBANK
        raise AssertionError(url)


@pytest.fixture
def ncbi_net(monkeypatch):
    f = NcbiFake()
    monkeypatch.setattr(ncbi, "fetch", f)
    monkeypatch.setattr(ncbi, "_gap", {"last": 0.0})
    for p in (ncbi.data_dir() / "ncbi").glob("*"):
        p.unlink()
    return f


def test_ncbi_client_search_cache_rate_limit_and_errors(ncbi_net, monkeypatch):
    """NCBI E-utilities: search and records (the tool name and key sent, no e-mail), records cached and size-checked, at most 3
    requests a second without a key, a bad answer not kept in the cache, and a 429 retried."""
    # --- search and records ---
    records = ncbi.search("nuccore", "U1.1", key="K")
    assert [(r.accession, r.length, r.topology, r.title) for r in records] == [("U1.1", 1834, "circular", "A made-up plasmid")]
    first = ncbi_net.requests[0][1]
    assert first["tool"] == "openDNA" and first["api_key"] == "K" and "email" not in first and first["db"] == "nuccore"
    assert ncbi.search("nuccore", "nothing") == []
    with pytest.raises(NcbiError, match="Type an accession"):
        ncbi.search("nuccore", "  ")
    with pytest.raises(NcbiError, match="503"):
        ncbi.search("nuccore", "broken")
    protein = ncbi.search("protein", "P12345")[0]
    assert protein.uniprot == "P12345" and records[0].uniprot == ""
    # --- get is cached and checked ---
    record = ncbi.search("nuccore", "U1.1")[0]
    path, cached = ncbi.get(record)
    assert path.suffix == ".gb" and not cached and path.read_text(encoding="utf-8").startswith("LOCUS")
    ncbi_net.requests.clear()
    assert ncbi.get(record) == (path, True) and not ncbi_net.requests, "the second time comes from the cache"
    path.unlink()
    big = ncbi.Record("1", "NC_BIG.1", "huge", ncbi.MAX_BASES + 1, "x", "", "nuccore")
    with pytest.raises(NcbiError, match="too long"):
        ncbi.get(big)
    protein_path, _ = ncbi.get(ncbi.search("protein", "P12345")[0])
    assert ncbi.protein_of(protein_path) == ("P12345.1", "MKWHAA")
    # --- rate limit ---
    t = time.monotonic()
    for _ in range(3):
        ncbi.call("esearch.fcgi", {"db": "nuccore", "term": "x"})
    assert time.monotonic() - t >= 0.6, "at most 3 requests a second without a key"
    t = time.monotonic()
    for _ in range(3):
        ncbi.call("esearch.fcgi", {"db": "nuccore", "term": "x"}, key="K")
    assert time.monotonic() - t < 0.6
    # --- bad answer is not kept ---
    record = ncbi.Record("5", "X1.1", "t", 10, "x", "", "nuccore")
    monkeypatch.setattr(ncbi, "fetch", lambda url, timeout=60: "Error: nothing here")
    with pytest.raises(NcbiError, match="did not return"):
        ncbi.get(record)
    assert not ncbi.path_of(record).exists()
    # --- too many requests is retried ---
    answers = [NcbiError("eutils.ncbi.nlm.nih.gov answered 429 Too Many Requests."), "fine"]

    def flaky(url, timeout=60):
        answer = answers.pop(0)
        if isinstance(answer, Exception):
            raise answer
        return answer
    monkeypatch.setattr(ncbi, "fetch", flaky)
    monkeypatch.setattr(ncbi.time, "sleep", lambda seconds: None)
    assert ncbi.call("esearch.fcgi", {"db": "nuccore", "term": "x"}) == "fine"


def test_structure_sources_cache_and_auto_order(net, monkeypatch, tmp_path):
    """Each structure source and its input checks, the cache, Auto's order (AlphaFold DB first; PDB then ESMFold without an
    accession), Auto reporting every failure, and the local command."""
    # --- each source ---
    af = structure.get(PROTEIN, "AlphaFold DB", "p12345")
    assert af.plddt and af.source == "AlphaFold DB" and "AF-P1-F1" in af.detail and not af.cached
    entry = structure.get(PROTEIN, "PDB")
    assert not entry.plddt and "1ABC" in entry.detail
    esm = structure.get(PROTEIN, "ESMFold")
    assert esm.plddt and structure.mean_bfactor(esm.pdb) == pytest.approx(50.0), "the 0-1 values are given as 0-100"
    with pytest.raises(StructureError, match="UniProt accession"):
        structure.get(PROTEIN, "AlphaFold DB", "")
    with pytest.raises(StructureError, match="25 amino acids"):
        structure.get("MKWV", "PDB")
    with pytest.raises(StructureError, match="up to 400"):
        structure.get("A" * 401, "ESMFold")
    with pytest.raises(StructureError, match="not amino acids"):
        structure.get("MKWV1*B", "ESMFold")
    with pytest.raises(StructureError, match="404"):
        structure.get(PROTEIN, "AlphaFold DB", "BADACC1")
    # --- cache and Auto's order ---
    clean_cache()
    net.calls.clear()
    first = structure.get(PROTEIN, "Auto", "P12345")          # AlphaFold DB comes first
    assert first.source == "AlphaFold DB" and net.calls[0].endswith("/P12345")
    net.calls.clear()
    again = structure.get(PROTEIN, "Auto", "P12345")
    assert again.cached and again.source == "AlphaFold DB" and not net.calls, "the cache answers without the network"
    assert not structure.get(PROTEIN, "Auto", "P12345", use_cache=False).cached and net.calls
    # --- Auto without an accession falls through ---
    clean_cache()
    net.calls.clear()
    net.pdb_hits = False
    result = structure.get(PROTEIN, "Auto")
    assert result.source == "ESMFold" and any("search.rcsb.org" in u for u in net.calls) and not any("alphafold" in u for u in net.calls)
    # --- Auto reports every failure ---
    def down(url, data=None, headers=None, timeout=30):
        raise StructureError("no network")
    monkeypatch.setattr(structure, "fetch", down)
    clean_cache()
    with pytest.raises(StructureError) as error:
        structure.get(PROTEIN, "Auto")
    assert "No structure found" in str(error.value) and "PDB:" in str(error.value) and "ESMFold:" in str(error.value)
    # --- local command ---
    script = tmp_path / "fold.py"
    script.write_text("import sys\nopen(sys.argv[2], 'w').write(open(sys.argv[1]).read().count('A') * 'x' and "
                      f"{pdb(80.0)!r})\n")
    command = f'"{sys.executable}" "{script}" {{fasta}} {{pdb}}'
    result = structure.get(PROTEIN, "Local command", command=command)
    assert result.source == "Local command" and "ATOM" in result.pdb
    with pytest.raises(StructureError, match="needs {fasta}"):
        structure.get(PROTEIN, "Local command", command="echo hi")
    with pytest.raises(StructureError, match="wrote no PDB"):
        structure.get(PROTEIN, "Local command", command=f'"{sys.executable}" -c "pass" {{fasta}} {{pdb}}')


def test_structure_page_script_safety_and_alphafold_reply_not_followed(monkeypatch):
    """The 3D viewer page: names escaped, a stray </script> in the details or the PDB cannot end the page script; an AlphaFold
    reply pointing at another address is not followed. Also the UniProt accession read from a feature."""
    # --- accession and page ---
    assert structure.accession_of([("label", "x"), ("db_xref", "UniProtKB/Swiss-Prot:P69905")]) == "P69905"
    assert structure.accession_of([("db_xref", "GeneID:1")]) == ""
    s = structure.Structure(pdb(90.0), "AlphaFold DB", "detail </script> here", True)
    html = structure.viewer_html(s, "My <protein>")
    assert "$3Dmol" in html and "ATOM" in html and "My &lt;protein&gt;" in html
    assert html.count("</script>") == 2, "a stray </script> in the data cannot end a script early"
    # --- an AlphaFold reply pointing elsewhere is not followed ---
    # (another scheme, another host, a look-alike host, AlphaFold only in the query string, plain http)
    for elsewhere in ("file:///C:/Windows/win.ini", "https://evil.example.com/x.pdb", "http://alphafold.ebi.ac.uk/files/x.pdb",
                      "https://alphafold.ebi.ac.uk.evil.example.com/x.pdb", "https://evil.example.com/?https://alphafold.ebi.ac.uk/"):
        calls = []

        def fake(url, data=None, headers=None, timeout=30, elsewhere=elsewhere):
            calls.append(url)
            return json.dumps([{"modelEntityId": "AF-X", "pdbUrl": elsewhere}]).encode()
        monkeypatch.setattr(structure, "fetch", fake)
        with pytest.raises(StructureError, match="somewhere else"):
            structure.from_alphafold("P12345")
        assert len(calls) == 1, f"only the API address was fetched, never {elsewhere}"
    # --- a closing script tag in the PDB cannot end the page script ---
    s = structure.Structure(pdb(90.0) + "REMARK </script><script>alert(1)</script>\n", "ESMFold", "d", True)
    html = structure.viewer_html(s, "t")
    assert "<script>alert(1)" not in html and html.count("</script>") == 2


def test_pdb_entries_and_the_model_read_from_a_pdb_text(net):
    """PDB search lists each entry once with its identity, coverage and details and gets the one picked (requests say program and
    version); linked entries (feature links and UniProt mapping) compared on their best chain, closest first; the model of a PDB
    text (first model, ligands without water, helices and sheets from alpha carbons) and the protein's residues found on chains
    with their own numbering, a missing residue and a substitution."""
    # --- PDB search
    assert structure.USER_AGENT.startswith("openDNA/") and structure.USER_AGENT != "openDNA"
    hits = structure.search_pdb(PROTEIN)
    assert [h.entry for h in hits] == ["1ABC", "2XYZ"]
    first = hits[0]
    assert (first.identity, round(first.covered, 2), first.title, first.organism, first.chains, first.resolution) == \
        (0.95, 0.5, "Crambin", "Crambe hispanica", "A, B", 1.5)
    assert hits[1].identity == 0.8 and hits[1].title == "", "an entry the data API did not describe is still listed"
    found = structure.pdb_entry(first)
    assert found.source == "PDB" and "PDB 1ABC: Crambin" in found.detail and "95% identity" in found.detail and "chains A, B" in found.detail
    assert found.info == {"Protein": "Crambin", "Source": "RCSB PDB (experimental)", "Entry": "1ABC", "Method": "X-ray diffraction",
                          "Resolution": "1.50 Å", "Organism": "Crambe hispanica", "Identity": "95%", "Chain covers": "50% of the protein",
                          "Chains": "A, B", "Link": "https://www.rcsb.org/structure/1ABC"}
    with pytest.raises(StructureError, match="25 amino acids"):
        structure.search_pdb("MKWV")

    # --- linked entries
    quals = [("db_xref", "PDB:1abc"), ("db_xref", "PDB:1ABC"), ("db_xref", "UniProtKB/Swiss-Prot:P12345"), ("db_xref", "PDB:2xyz_A"),
             ("db_xref", "PDB:bad"), ("note", "PDB:9ZZZ")]
    assert structure.pdb_links(quals) == ["1ABC", "2XYZ"]
    assert structure.similarity("MKWHLLEGRS", "KWHLLEGRS") == (1.0, 0.9)
    assert structure.similarity("MKWHLLEGRS", "") == (0.0, 0.0)
    hits = structure.linked_pdb(PROTEIN, ["1ABC"], "P12345")
    assert [h.entry for h in hits] == ["3DEF", "1ABC"], "the entry that holds the protein first"
    best = hits[0]
    assert (best.identity, round(best.covered, 2), best.chains, best.organism, best.title) == \
        (1.0, round((len(PROTEIN) - 10) / len(PROTEIN), 2), "B, D", "Crambe hispanica", "The protein in a complex"), "its best chain counts"
    assert hits[1].identity < 0.5
    assert structure.linked_pdb(PROTEIN) == []
    with pytest.raises(StructureError, match="6 to 10"):
        structure.linked_pdb(PROTEIN, accession="P1")

    # --- the model of a PDB text
    helix = [(2.3 * math.cos(math.radians(100 * i)), 2.3 * math.sin(math.radians(100 * i)), 1.5 * i) for i in range(12)]   # ideal alpha helix

    def strand(z):   # an extended chain, zigzag along x
        return [(3.3 * i, 0.9 if i % 2 else -0.9, z) for i in range(6)]
    lines = [line(i + 1, "ALA", "A", i + 1, p) for i, p in enumerate(helix)]
    lines += [line(100 + i, "VAL", "B", i + 1, p) for i, p in enumerate(strand(30))]
    lines += [line(200 + i, "VAL", "C", i + 1, p) for i, p in enumerate(strand(34.8))]   # beside B: a sheet
    lines += [line(300 + i, "VAL", "D", i + 1, p) for i, p in enumerate(strand(80))]     # alone: a loop
    lines += [line(400, "HEM", "A", 900, (0, 0, 0), name="FE", record="HETATM", element="FE"),
              line(401, "HOH", "A", 901, (1, 1, 1), name="O", record="HETATM", element="O"),
              "ENDMDL", line(500, "ALA", "E", 1, (9, 9, 9))]
    residues, ligands = structure.model_of("\n".join(lines))
    assert ligands == [("FE", (0.0, 0.0, 0.0), "HEM A900")],"water and the second model are left out"
    ss = {c: "".join(r.ss for r in residues if r.chain == c) for c in "ABCDE"}
    assert ss == {"A": "H" * 12, "B": "EEEEEE", "C": "EEEEEE", "D": "------", "E": ""}
    # --- the protein on three chains: A numbered from 10 without its W and with an A for an L, B whole, C unrelated
    protein = "MKWHLLEGRS"
    three = {one: name for name, one in structure.THREE.items() if name != "MSE"}

    def chain(c, seq, start):
        return [line(i + 1, three[a], c, start + i, (3.8 * i, 0, 0)) for i, a in enumerate(seq)]
    model, _ = structure.model_of("\n".join(chain("A", "MKHLAEGRS", 10) + chain("B", protein, 1) + chain("C", "GGGG", 1)))
    found = structure.residue_map(protein, model)
    names = {i: [model[k].chain + model[k].number for k in ks] for i, ks in found.items()}
    assert len(names) == 10 and names[0] == ["A10", "B1"]
    assert names[2] == ["B3"], "W is missing from chain A"
    assert names[5] == ["A14", "B6"], "a substitution keeps its place"
    assert not any(n.startswith("C") for ns in names.values() for n in ns), "a chain that does not match is left out"
