"""NCBI import: search, records, the cache, the rate limit, errors, the dialog and File > Import from NCBI (the network replaced by fakes)."""
import json
import os
import sys
import tempfile
import time
from pathlib import Path
from urllib.parse import parse_qs, urlparse

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ["GENOMEEDIT_DATA_DIR"] = tempfile.mkdtemp(prefix="genomeedit_test_")
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "app"))

import pytest
from PySide6.QtTest import QTest
from PySide6.QtWidgets import QApplication

import ncbi
import structure
from ncbi import NcbiError

app = QApplication.instance() or QApplication([])
GENBANK = (ROOT / "demo" / "demo_vector.gb").read_text(encoding="utf-8")
SUMMARY = {"uids": ["11"], "11": {"uid": "11", "caption": "U1", "title": "A made-up plasmid", "slen": 1834, "organism": "synthetic construct",
                                   "updatedate": "2026/01/02", "accessionversion": "U1.1", "topology": "circular", "sourcedb": "genbank"}}
PROTEIN_SUMMARY = {"uids": ["22"], "22": {"uid": "22", "caption": "P12345", "title": "RecName: Full=Tiny; AltName: Full=x", "slen": 6,
                                          "organism": "Nowhere", "updatedate": "2026/01/02", "accessionversion": "P12345.1",
                                          "moltype": "aa", "sourcedb": "swiss_prot"}}


class Fake:
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
def fake(monkeypatch):
    f = Fake()
    monkeypatch.setattr(ncbi, "fetch", f)
    monkeypatch.setattr(ncbi, "_gap", {"last": 0.0})
    for p in (ncbi.data_dir() / "ncbi").glob("*"):
        p.unlink()
    return f


def test_ncbi_calls(fake, monkeypatch):
    """Search and records, the cache and its checks, the rate limit, a bad answer that is not kept, and a 429 that is retried."""
    # --- search and records ---
    records = ncbi.search("nuccore", "U1.1", key="K")
    assert [(r.accession, r.length, r.topology, r.title) for r in records] == [("U1.1", 1834, "circular", "A made-up plasmid")]
    first = fake.requests[0][1]
    assert first["tool"] == "OpenGenEdit" and first["api_key"] == "K" and "email" not in first and first["db"] == "nuccore"
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
    fake.requests.clear()
    assert ncbi.get(record) == (path, True) and not fake.requests, "the second time comes from the cache"
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


def wait(dialog, what=lambda d: d.find.isEnabled()):
    t = time.perf_counter()
    QTest.qWait(20)
    while not what(dialog) and time.perf_counter() - t < 5:
        QTest.qWait(20)


def test_dialog(fake, monkeypatch):
    """The NCBI dialog opens a record in a window (the next one in a new window), opens a protein as a 3D structure request,
    and says when nothing is found or the service fails."""
    # --- opens a record in a window ---
    import main
    win = main.MainWindow()
    win.show()
    d = main.NcbiDialog(win)
    d.query.setText("U1.1")
    d.search()
    wait(d)
    assert d.table.rowCount() == 1 and d.table.item(0, 0).text() == "U1.1" and d.table.item(0, 2).text() == "1,834 bp"
    assert d.open_button.isEnabled() and d.open_button.text() == "Open"
    d.open()
    wait(d, lambda d: d.opened is not None)
    assert d.opened is win and win.editor.untitled and win.editor.length == 1834 and win.editor.features
    # the next record opens in a new window: nothing open is touched
    d2 = main.NcbiDialog(win)
    d2.query.setText("U1.1")
    d2.search()
    wait(d2)
    d2.open()
    wait(d2, lambda d: d.opened is not None)
    other = d2.opened
    assert other is not win and other.editor.length == 1834
    for w in (win, other):
        w.editor.doc.modified = False
        w.close()
    d.settings.clear()
    # --- protein and errors ---
    import main
    win = main.MainWindow()
    d = main.NcbiDialog(win)
    d.database.setCurrentIndex(1)
    d.query.setText("P12345")
    d.search()
    wait(d)
    assert d.open_button.text() == "3D structure…"
    shown = []
    monkeypatch.setattr(structure, "get", lambda *a, **k: shown.append(a))
    d.open()
    wait(d, lambda d: d.opened is not None)
    assert d.opened.protein == "MKWHAA" and d.opened.accession.text() == "P12345" and d.opened.name == "RecName: Full=Tiny"
    d.opened.close()
    d.query.setText("nothing")
    d.search()
    wait(d)
    assert "Nothing found" in d.status.text() and not d.open_button.isEnabled()
    d.query.setText("broken")
    d.search()
    wait(d)
    assert "503" in d.status.text() and d.find.isEnabled()
    d.settings.clear()
    win.close()


