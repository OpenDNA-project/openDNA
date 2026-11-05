"""3D structures: the sources (with the network replaced by fakes), the cache, Auto's order, the local command, the page, the window."""
import json
import os
import sys
import tempfile
import time
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ["GENOMEEDIT_DATA_DIR"] = tempfile.mkdtemp(prefix="genomeedit_test_")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "app"))

import pytest
from PySide6.QtTest import QTest
from PySide6.QtWidgets import QApplication

import structure
from structure import StructureError

app = QApplication.instance() or QApplication([])
PROTEIN = "TTCCPSIVARSNFNVCRLPGTPEALCATYTGCIIIPGATCPGDYAN"


def atom(i, b):
    return f"ATOM  {i:5d}  CA  THR A{i:4d}      -8.418  -7.381   3.408  1.00{b:6.2f}           C  "


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
        if "search.rcsb.org" in url:
            body = json.loads(data)
            assert body["query"]["parameters"]["value"] == PROTEIN and body["query"]["parameters"]["identity_cutoff"] == 0.9
            return json.dumps({"result_set": [{"identifier": "1ABC_1", "score": 1.0}]}).encode() if self.pdb_hits else b""
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


def test_sources_cache_auto_order_and_local_command(net, monkeypatch, tmp_path):
    """Each source and its input checks, the cache and Auto's order, Auto without an accession, Auto reporting every failure, the local command."""
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


def test_accession_page_and_script_tag_safety(monkeypatch):
    """Accession lookup and the viewer page, a stray </script> in the data or the PDB cannot end the page script, and an AlphaFold reply pointing elsewhere is not followed."""
    # --- accession and page ---
    assert structure.accession_of([("label", "x"), ("db_xref", "UniProtKB/Swiss-Prot:P69905")]) == "P69905"
    assert structure.accession_of([("db_xref", "GeneID:1")]) == ""
    s = structure.Structure(pdb(90.0), "AlphaFold DB", "detail </script> here", True)
    html = structure.viewer_html(s, "My <protein>")
    assert "$3Dmol" in html and "ATOM" in html and "My &lt;protein&gt;" in html
    assert html.count("</script>") == 2, "a stray </script> in the data cannot end a script early"
    # --- an AlphaFold reply pointing elsewhere is not followed ---
    calls = []

    def fake(url, data=None, headers=None, timeout=30):
        calls.append(url)
        return json.dumps([{"modelEntityId": "AF-X", "pdbUrl": "file:///C:/Windows/win.ini"}]).encode()
    monkeypatch.setattr(structure, "fetch", fake)
    with pytest.raises(StructureError, match="somewhere else"):
        structure.from_alphafold("P12345")
    assert len(calls) == 1, "only the API address was fetched, never the file the reply named"
    # --- a closing script tag in the PDB cannot end the page script ---
    s = structure.Structure(pdb(90.0) + "REMARK </script><script>alert(1)</script>\n", "ESMFold", "d", True)
    html = structure.viewer_html(s, "t")
    assert "<script>alert(1)" not in html and html.count("</script>") == 2


def test_window(net, monkeypatch):
    from structure_dialog import StructureDialog
    d = StructureDialog(None, "tiny", PROTEIN, "P12345")
    assert not d.open.isEnabled() and "ESM" in structure.NOTICE["ESMFold"]
    d.source.setCurrentText("ESMFold")
    assert not d.accession.isEnabled() and "sent to api.esmatlas.com" in d.notice.text()
    d.start()
    t = time.perf_counter()
    while not d.get.isEnabled() and time.perf_counter() - t < 5:
        QTest.qWait(20)
    assert d.open.isEnabled() and d.save.isEnabled() and "ESMFold" in d.status.text()
    page = d.view_file()
    assert page.exists() and "3Dmol" in page.read_text(encoding="utf-8")
    d.source.setCurrentText("AlphaFold DB")
    d.accession.setText("BADACC1")
    d.start()
    t = time.perf_counter()
    while not d.get.isEnabled() and time.perf_counter() - t < 5:
        QTest.qWait(20)
    assert "404" in d.status.text()
    d.settings.clear()
