"""Self-update: versions, the release JSON, the checksum, the allowed hosts, install and rollback on fake files, the once-a-day throttle,
skip-this-version, the disabled flag and the dialog (the network and the relaunch replaced by fakes; nothing is downloaded or started)."""
import hashlib
import json
import os
import sys
import tempfile
import threading
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ["OPENDNA_DATA_DIR"] = tempfile.mkdtemp(prefix="opendna_test_")
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "app"))

import pytest
from PySide6.QtTest import QTest
from PySide6.QtWidgets import QApplication, QWidget

import project
import tool_settings
import updater
from primer_db import data_dir
from updater import Update, UpdateError

app = QApplication.instance() or QApplication([])

EXE_BYTES = b"MZ" + b"fake exe " * 1000
DIGEST = hashlib.sha256(EXE_BYTES).hexdigest()
BASE = f"https://github.com/{project.REPO}/releases/download/v0.2.0"


def release(tag="v0.2.0", assets=True, **extra):
    """The JSON GitHub gives for releases/latest."""
    data = {"tag_name": tag, "html_url": f"{project.REPO_URL}/releases/tag/{tag}", "body": "Fixes things.", "draft": False, "prerelease": False,
            "assets": [{"name": updater.ASSET_NAME, "browser_download_url": f"{BASE}/openDNA.exe"},
                       {"name": updater.CHECKSUM_NAME, "browser_download_url": f"{BASE}/openDNA.exe.sha256"},
                       {"name": "other.zip", "browser_download_url": f"{BASE}/other.zip"}] if assets else []}
    data.update(extra)
    return json.dumps(data)


class Fake:
    """A fake of updater.fetch: the release JSON, the checksum and the exe by URL; remembers the requests."""

    def __init__(self, latest=None, checksum=None, exe=EXE_BYTES, error=None):
        self.latest = release() if latest is None else latest
        self.checksum = f"{DIGEST}  {updater.ASSET_NAME}\n" if checksum is None else checksum
        self.exe, self.error, self.requests = exe, error, []

    def __call__(self, url, timeout=30, sink=None, progress=None, cancel=None):
        self.requests.append(url)
        if self.error:
            raise self.error
        if url == project.LATEST_RELEASE_API:
            return self.latest
        if url.endswith(".sha256"):
            return self.checksum
        if sink is None:
            raise AssertionError("the exe must be streamed to a file")
        half = len(self.exe) // 2
        sink.write(self.exe[:half])
        if progress:
            progress(half, len(self.exe))
        if cancel is not None and cancel.is_set():
            raise UpdateError("Download cancelled.")
        sink.write(self.exe[half:])
        if progress:
            progress(len(self.exe), len(self.exe))


def forget_state():
    """The remembered update state and the setting emptied (merged tests call this between their sections)."""
    (data_dir() / "updates.ini").unlink(missing_ok=True)
    store = tool_settings.store()
    store.remove("Updates/check")
    store.sync()
    tool_settings._cache.clear()


@pytest.fixture(autouse=True)
def sandbox(monkeypatch):
    """No waiting between retries; the remembered state and the setting empty before and after; updates switched on by default."""
    monkeypatch.setattr(updater, "RETRY_DELAY", 0)
    monkeypatch.setattr(project, "UPDATES_ENABLED", True)
    forget_state()
    yield
    forget_state()


@pytest.fixture
def messages(monkeypatch):
    """Message boxes recorded instead of shown."""
    shown = []
    monkeypatch.setattr(updater, "show_message", lambda parent, title, text: shown.append(text))
    return shown


def wait_for(condition, ms=5000):
    for _ in range(ms // 10):
        if condition():
            return True
        QTest.qWait(10)
    return condition()


def update(**kw):
    values = dict(version="0.2.0", tag="v0.2.0", asset_url=f"{BASE}/openDNA.exe", sha256_url=f"{BASE}/openDNA.exe.sha256",
                  notes="Notes", page_url=f"{project.REPO_URL}/releases/tag/v0.2.0")
    values.update(kw)
    return Update(**values)



def fresh_dir(parent, name):
    """A new empty folder inside `parent`: merged tests that need an empty folder per section use one each."""
    folder = parent / name
    folder.mkdir()
    return folder


# ---- versions ------------------------------------------------------------------------------------------------------------------------

def test_versions():
    """Version comparison (v prefix, longer numbers, pre-releases, junk, None) and parsing."""
    # --- comparison
    for candidate, current, newer in [
        ("v0.10.0", "0.9.0", True), ("0.2.0", "0.1.0", True), ("v1.0.0", "0.99.99", True), ("0.1.1", "0.1.0", True),
        ("1.0", "1.0.0", False), ("1.0.0.1", "1.0", True), ("0.1.0", "0.1.0", False), ("v0.1.0", "0.1.0", False),
        ("0.9.0", "0.10.0", False), ("0.0.9", "0.1.0", False),
        ("0.2.0-beta1", "0.1.0", False), ("latest", "0.1.0", False), ("", "0.1.0", False), (None, "0.1.0", False), ("0.2.0", "dev", False),
    ]:
        assert updater.is_newer(candidate, current) is newer, f"{candidate!r} vs {current!r}"

    # --- parsing tolerates a leading v and refuses prereleases
    assert updater.parse_version("v0.10.0") == (0, 10, 0) and updater.parse_version(" V2.1 ") == (2, 1)
    assert updater.parse_version("1.2.3-rc1") is None and updater.parse_version("1.2.3+build") is None


# ---- the release JSON ----------------------------------------------------------------------------------------------------------------

def test_release_json(monkeypatch):
    """A newer release gives the update (through the default API address and the module-level fetch too); up to date, older, draft
    and prerelease offer nothing; a tag that is not a version is reported."""
    # --- a newer release gives the update
    result = updater.check("0.1.0", Fake())
    u = result.update
    assert result.error is None and u.version == "0.2.0" and u.tag == "v0.2.0" and u.notes == "Fixes things."
    assert u.asset_url == f"{BASE}/openDNA.exe" and u.sha256_url == f"{BASE}/openDNA.exe.sha256"
    assert u.page_url == f"{project.REPO_URL}/releases/tag/v0.2.0"
    assert updater.check_for_update("0.1.0", Fake()) == u

    # --- it asks the releases API with the project version by default
    fake = Fake(release(tag=f"v{project.VERSION}"))
    assert updater.check(fetch=fake).update is None and fake.requests == [project.LATEST_RELEASE_API]

    # --- the module-level fetch is used when none is given
    monkeypatch.setattr(updater, "fetch", Fake())
    assert updater.check_for_update("0.1.0").version == "0.2.0"

    # --- up to date, older, draft and prerelease offer nothing
    for tag, extra in (("v0.1.0", {}), ("v0.0.9", {}), ("v0.3.0", {"draft": True}), ("v0.3.0", {"prerelease": True}),
                       ("v0.3.0-beta1", {})):
        result = updater.check("0.1.0", Fake(release(tag, **extra)))
        assert result.update is None, tag
        # the beta tag name is explained, nothing is offered; the others are simply not newer
        assert (result.error is not None and "not understood" in result.error) == (tag == "v0.3.0-beta1"), (tag, result.error)

    # --- a tag that is not a version is reported
    result = updater.check("0.1.0", Fake(release("nightly")))
    assert result.update is None and "nightly" in result.error


def test_release_assets_and_errors():
    """A newer release without the exe or its checksum cannot be installed; assets on other hosts or repositories are refused; a
    page link that is not GitHub falls back to the releases page; network errors and garbage never raise."""
    # --- without the exe or its checksum
    only_exe = json.loads(release())
    only_exe["assets"] = only_exe["assets"][:1]
    only_sum = json.loads(release())
    only_sum["assets"] = only_sum["assets"][1:2]
    for data, missing in ((json.loads(release(assets=False)), updater.ASSET_NAME), (only_exe, updater.CHECKSUM_NAME),
                          (only_sum, updater.ASSET_NAME)):
        result = updater.check("0.1.0", Fake(json.dumps(data)))
        assert result.update is None and missing in result.error and result.page_url, missing

    # --- assets on other hosts or other repositories are refused
    for url in ("https://evil.example.com/openDNA.exe", "http://github.com/x/openDNA.exe",
                "https://github.com/someone/else/releases/download/v9/openDNA.exe"):
        data = json.loads(release())
        data["assets"][0]["browser_download_url"] = url
        result = updater.check("0.1.0", Fake(json.dumps(data)))
        assert result.update is None and "Refused" in result.error, url

    # --- a page link that is not GitHub falls back to the releases page
    result = updater.check("0.1.0", Fake(release(html_url="https://evil.example.com/x")))
    assert result.update.page_url == project.RELEASES_URL

    # --- network errors and garbage never raise
    assert updater.check("0.1.0", Fake(error=UpdateError("Could not reach api.github.com: offline"))).error.endswith("offline")
    assert "failed" in updater.check("0.1.0", Fake(error=RuntimeError("boom"))).error
    for junk in ("not json", "[1, 2]", "", b"\xff\xfe\x00"):
        result = updater.check("0.1.0", Fake(junk))
        assert result.update is None and "understood" in result.error
    assert updater.check_for_update("0.1.0", Fake(error=UpdateError("x"))) is None


# ---- the allowed addresses -----------------------------------------------------------------------------------------------------------

def test_allowed_addresses(monkeypatch):
    """Host validation (https GitHub hosts only), the real fetch refusing a bad address before any network, and a redirect to
    another host refused."""
    # --- host validation
    for url, ok in [
        ("https://api.github.com/repos/a/b/releases/latest", True), ("https://github.com/a/b/releases/download/v1/x.exe", True),
        ("https://objects.githubusercontent.com/github-production-release-asset/1", True), ("https://GitHub.com/a", True),
        ("https://release-assets.githubusercontent.com/x", True),
        ("http://github.com/a", False), ("ftp://github.com/a", False), ("https://evil.example.com/github.com", False),
        ("https://github.com.evil.example.com/a", False), ("https://evilgithub.com/a", False), ("https://user:pw@github.com/a", False),
        ("https://github.com@evil.example.com/a", False), ("https://github.com:8443/a", False), ("https://gist.github.com/a", False),
        ("", False), (None, False), ("file:///C:/x.exe", False),
    ]:
        if ok:
            assert updater.check_url(url) == url, f"{url!r}"
        else:
            with pytest.raises(UpdateError):
                updater.check_url(url)

    # --- the real fetch refuses a bad address before any network
    def no_network(*a, **k):
        raise AssertionError("must not connect")
    monkeypatch.setattr(updater.urllib.request, "build_opener", no_network)
    with pytest.raises(UpdateError):
        updater.fetch("http://github.com/a/b")
    with pytest.raises(UpdateError):
        updater.fetch("https://evil.example.com/a")

    # --- a redirect to another host is refused
    handler = updater._SafeRedirects()
    request = updater.urllib.request.Request("https://github.com/a/b/releases/download/v1/openDNA.exe")
    with pytest.raises(UpdateError):
        handler.redirect_request(request, None, 302, "Found", {}, "https://evil.example.com/openDNA.exe")
    assert handler.redirect_request(request, None, 302, "Found", {}, "https://objects.githubusercontent.com/x") is not None


# ---- the checksum and the download ---------------------------------------------------------------------------------------------------

def test_checksum_and_download(tmp_path):
    """Checksum file forms; the download verifies and returns the new exe; it refuses a mismatch, a missing or unreadable
    checksum and addresses off GitHub; a failed or cancelled download leaves no partial file; a missing folder is a readable error."""
    # --- checksum file forms
    upper = DIGEST.upper()
    assert updater.parse_checksum(f"{DIGEST}  {updater.ASSET_NAME}\n") == DIGEST
    assert updater.parse_checksum(f"{DIGEST} *{updater.ASSET_NAME}") == DIGEST
    assert updater.parse_checksum(DIGEST + "\n") == DIGEST
    assert updater.parse_checksum(f"{upper}  {updater.ASSET_NAME}") == DIGEST
    assert updater.parse_checksum(f"{'a' * 64}  other.zip\n{DIGEST}  {updater.ASSET_NAME}\n") == DIGEST
    assert updater.parse_checksum((DIGEST + "\r\n").encode("utf-16")) == DIGEST          # PowerShell's > writes UTF-16 with a BOM
    assert updater.parse_checksum(b"\xef\xbb\xbf" + DIGEST.encode()) == DIGEST
    for bad in ("", "hello", "abc  openDNA.exe", f"{DIGEST}  other.zip", DIGEST[:-1], DIGEST + "0"):
        assert updater.parse_checksum(bad) is None, bad

    # --- verifies and returns the new exe
    folder = fresh_dir(tmp_path, "ok")
    seen = []
    path = updater.download_update(update(), folder, progress=lambda d, t: seen.append((d, t)), fetch=Fake())
    assert path == folder / f"{updater.ASSET_NAME}.new" and path.read_bytes() == EXE_BYTES
    assert not (folder / f"{updater.ASSET_NAME}.part").exists()
    assert seen[-1] == (len(EXE_BYTES), len(EXE_BYTES))

    # --- refuses a checksum mismatch and deletes the file
    folder = fresh_dir(tmp_path, "mismatch")
    with pytest.raises(UpdateError, match="does not match"):
        updater.download_update(update(), folder, fetch=Fake(exe=EXE_BYTES + b"tampered"))
    assert list(folder.iterdir()) == []

    # --- refuses a missing or unreadable checksum
    folder = fresh_dir(tmp_path, "nosum")
    fake = Fake()
    with pytest.raises(UpdateError, match="no checksum"):
        updater.download_update(update(sha256_url=""), folder, fetch=fake)
    assert not any(u.endswith(".exe") for u in fake.requests)                # nothing was even downloaded
    with pytest.raises(UpdateError, match="expected form"):
        updater.download_update(update(), folder, fetch=Fake(checksum="<html>404</html>"))
    assert list(folder.iterdir()) == []

    # --- refuses addresses off github
    folder = fresh_dir(tmp_path, "offgithub")
    for bad in (update(asset_url="https://evil.example.com/openDNA.exe"), update(sha256_url="http://github.com/x.sha256")):
        fake = Fake()
        with pytest.raises(UpdateError, match="Refused"):
            updater.download_update(bad, folder, fetch=fake)
        assert fake.requests == []

    # --- a failed or cancelled download leaves no partial file
    folder = fresh_dir(tmp_path, "cancelled")
    fake = Fake()
    stop = threading.Event()
    stop.set()
    with pytest.raises(UpdateError, match="cancelled"):
        updater.download_update(update(), folder, fetch=fake, cancel=stop)
    assert list(folder.iterdir()) == []
    with pytest.raises(UpdateError, match="Could not reach"):
        def broken(url, timeout=30, sink=None, progress=None, cancel=None):
            if sink is None:
                return DIGEST
            sink.write(b"half")
            raise UpdateError("Could not reach github.com: reset")
        updater.download_update(update(), folder, fetch=broken)
    assert list(folder.iterdir()) == []

    # --- a missing folder is a readable error
    with pytest.raises(UpdateError, match="Could not write"):
        updater.download_update(update(), tmp_path / "nowhere", fetch=Fake())


# ---- installing ----------------------------------------------------------------------------------------------------------------------

def installed_in(folder):
    """A 'running exe' and a downloaded new one in a new folder."""
    folder.mkdir()
    current, new = folder / "openDNA.exe", folder / "openDNA.exe.new"
    current.write_bytes(b"OLD")
    new.write_bytes(b"NEW")
    return current, new


def fail_replace(monkeypatch, when):
    """os.replace in updater failing with PermissionError for the calls `when(src, dst)` says."""
    real, calls = os.replace, []

    def replace(src, dst):
        calls.append((Path(src).name, Path(dst).name))
        if when(Path(src), Path(dst)):
            raise PermissionError(13, "Access is denied")
        return real(src, dst)
    monkeypatch.setattr(updater.os, "replace", replace)
    return calls


def test_install_swaps_and_refuses(tmp_path, monkeypatch):
    """install_update swaps the files and relaunches, replaces a stale .old file, uses the running exe when frozen, and refuses
    from source, as the folder build, with missing files or another folder."""
    # --- swaps the files and relaunches
    current, new = installed_in(tmp_path / "swap")
    started = []
    updater.install_update(new, current, relaunch=started.append)
    assert current.read_bytes() == b"NEW" and (current.parent / "openDNA.exe.old").read_bytes() == b"OLD"
    assert not new.exists() and started == [current]

    # --- replaces a stale old file
    current, new = installed_in(tmp_path / "stale")
    (current.parent / "openDNA.exe.old").write_bytes(b"ANCIENT")
    updater.install_update(new, current, relaunch=lambda p: None)
    assert (current.parent / "openDNA.exe.old").read_bytes() == b"OLD"

    # --- refuses when not frozen
    assert not updater.is_frozen()
    _, new = installed_in(tmp_path / "source")
    with pytest.raises(UpdateError, match="source code"):
        updater.install_update(new, relaunch=lambda p: pytest.fail("must not start anything"))
    assert new.read_bytes() == b"NEW"
    assert "source code" in updater.install_blocker() and updater.install_blocker(Path("x.exe")) is None

    # --- the folder build does not replace itself
    current, new = installed_in(tmp_path / "folderbuild")
    with monkeypatch.context() as frozen:
        frozen.setattr(sys, "frozen", True, raising=False)
        frozen.setattr(sys, "executable", str(current))
        frozen.setattr(sys, "_MEIPASS", str(tmp_path / "somewhere" / "_MEI1234"), raising=False)   # one-file: a temp folder
        for system, download in (("darwin", "macos-arm64.zip"), ("linux", "linux-x86_64.tar.gz")):   # never the Windows exe there
            frozen.setattr(sys, "platform", system)
            assert "does not update itself" in updater.install_blocker() and download in updater.install_blocker()
        frozen.setattr(sys, "platform", "win32")   # (the one-file exe that updates itself is the Windows one)
        assert not updater.is_folder_build() and updater.install_blocker() is None
        for bundle in (current.parent / "_internal", current.parent):     # PyInstaller 6 and 5 one-folder layouts
            frozen.setattr(sys, "_MEIPASS", str(bundle), raising=False)
            assert updater.is_folder_build()
            assert "openDNA-folder.zip" in updater.install_blocker()
            with pytest.raises(UpdateError, match="folder version"):
                updater.install_update(new, relaunch=lambda p: pytest.fail("must not start anything"))
    assert current.read_bytes() == b"OLD" and new.read_bytes() == b"NEW"

    # --- uses the running exe when frozen
    current, new = installed_in(tmp_path / "running")
    with monkeypatch.context() as frozen:
        frozen.setattr(sys, "frozen", True, raising=False)
        frozen.setattr(sys, "platform", "win32")
        frozen.setattr(sys, "executable", str(current))
        updater.install_update(new, relaunch=lambda p: None)
    assert current.read_bytes() == b"NEW"

    # --- refuses missing files and another folder
    current, new = installed_in(tmp_path / "refuse")
    with pytest.raises(UpdateError, match="could not be found"):
        updater.install_update(tmp_path / "nothing.exe", current, relaunch=lambda p: None)
    other = tmp_path / "sub"
    other.mkdir()
    (other / "new.exe").write_bytes(b"NEW")
    with pytest.raises(UpdateError, match="same folder"):
        updater.install_update(other / "new.exe", current, relaunch=lambda p: None)
    assert current.read_bytes() == b"OLD"


def test_install_rolls_back(tmp_path, monkeypatch):
    """install_update rolls back when the new file cannot be moved in, points to the releases page when the folder cannot be
    written, retries while an antivirus holds the file, and rolls back when the new program does not start."""
    # --- the new file cannot be moved in
    current, new = installed_in(tmp_path / "movein")
    with monkeypatch.context() as patch:
        fail_replace(patch, lambda s, d: s == new)
        with pytest.raises(UpdateError, match="administrator"):
            updater.install_update(new, current, relaunch=lambda p: pytest.fail("must not start"))
    assert current.read_bytes() == b"OLD" and new.read_bytes() == b"NEW" and not (current.parent / "openDNA.exe.old").exists()

    # --- a folder that cannot be written points to the releases page
    current, new = installed_in(tmp_path / "readonly")
    with monkeypatch.context() as patch:
        calls = fail_replace(patch, lambda s, d: s == current)      # even the rename of the running exe is refused
        with pytest.raises(UpdateError) as e:
            updater.install_update(new, current, relaunch=lambda p: None)
    assert project.RELEASES_URL in str(e.value) and "administrator" in str(e.value)
    assert len(calls) == updater.RETRIES                              # it retried before giving up
    assert current.read_bytes() == b"OLD" and new.read_bytes() == b"NEW"

    # --- retries while an antivirus holds the file
    current, new = installed_in(tmp_path / "antivirus")
    refused = []

    def when(s, d):
        if s == current and len(refused) < 2:
            refused.append(1)
            return True
        return False
    with monkeypatch.context() as patch:
        fail_replace(patch, when)
        updater.install_update(new, current, relaunch=lambda p: None)
    assert len(refused) == 2 and current.read_bytes() == b"NEW"

    # --- the new program does not start
    current, new = installed_in(tmp_path / "nostart")

    def broken(path):
        raise OSError("not a valid Win32 application")
    with pytest.raises(UpdateError, match="old one was kept"):
        updater.install_update(new, current, relaunch=broken)
    assert current.read_bytes() == b"OLD" and new.read_bytes() == b"NEW" and not (current.parent / "openDNA.exe.old").exists()


def test_relaunch_and_cleanup(tmp_path, monkeypatch):
    """The relaunched child gets a clean PyInstaller environment; cleanup_old deletes the leftovers only, does nothing from source
    and waits for a held file."""
    # --- relaunch hands the child a clean PyInstaller environment
    seen = {}
    with monkeypatch.context() as patch:
        patch.setattr(updater.subprocess, "Popen", lambda args, **kw: seen.update(args=args, **kw))
        patch.setenv("_MEIPASS2", "C:/temp/_MEI1")
        patch.setenv("_PYI_APPLICATION_HOME_DIR", "C:/temp/_MEI1")
        updater.relaunch(tmp_path / "openDNA.exe")
    env = seen["env"]
    assert seen["args"] == [str(tmp_path / "openDNA.exe")] and env["PYINSTALLER_RESET_ENVIRONMENT"] == "1"
    assert "_MEIPASS2" not in env and not any(k.startswith("_PYI") for k in env)

    # --- cleanup_old deletes the leftovers only
    current, _ = installed_in(tmp_path / "leftovers")
    for suffix in (".old", ".part"):
        (current.parent / ("openDNA.exe" + suffix)).write_bytes(b"x")
    assert updater.cleanup_old(current) is True
    assert [p.name for p in current.parent.iterdir()] == ["openDNA.exe"] and current.read_bytes() == b"OLD"
    assert updater.cleanup_old(current) is True                          # nothing to do is fine

    # --- does nothing from source and waits for a held file
    current, _ = installed_in(tmp_path / "held")
    old = current.parent / "openDNA.exe.old"
    old.write_bytes(b"x")
    assert updater.cleanup_old() is True and old.exists()                 # not frozen, no exe given: not touched
    real = Path.unlink
    tries = []

    def held(self, missing_ok=False):
        tries.append(1)
        if len(tries) < 3:
            raise PermissionError(13, "in use")
        return real(self, missing_ok=missing_ok)
    with monkeypatch.context() as patch:
        patch.setattr(Path, "unlink", held)
        patch.setattr(updater.time, "sleep", lambda s: None)
        assert updater.cleanup_old(current, wait=60) is True and not old.exists()


# ---- what is remembered: throttle, skip, the switches --------------------------------------------------------------------------------

def test_remembered_state(monkeypatch):
    """The automatic check obeys the setting and the project flag and runs at most once a day (also when the clock goes back);
    a skipped version is remembered."""
    # --- obeys the setting and the project flag
    assert updater.check_enabled() and updater.auto_check_due()
    tool_settings.set_value(updater.TOOL, "check", False)
    assert not updater.check_enabled() and not updater.auto_check_due()
    tool_settings.set_value(updater.TOOL, "check", True)
    assert updater.auto_check_due()
    with monkeypatch.context() as patch:
        patch.setattr(project, "UPDATES_ENABLED", False)
        assert not updater.auto_check_due()

    # --- at most once a day
    forget_state()
    now = 1_000_000.0
    assert updater.auto_check_due(now)                                    # never checked
    updater.mark_checked(now)
    assert not updater.auto_check_due(now + 60) and not updater.auto_check_due(now + 23 * 3600)
    assert updater.auto_check_due(now + 24 * 3600)
    assert updater.auto_check_due(now - 3600)                             # the clock went back: do not wait for years

    # --- a skipped version is remembered
    forget_state()
    assert updater.skipped_version() == ""
    updater.skip_version("0.2.0")
    assert updater.skipped_version() == "0.2.0"


# ---- start-up and the menu action ----------------------------------------------------------------------------------------------------

def test_startup_check(monkeypatch, messages):
    """The start-up check opens the dialog for a new version and only once a day, is silent when up to date, skipped or offline,
    and does nothing at all while updates are disabled."""
    monkeypatch.setattr(project, "VERSION", "0.1.0")
    # --- opens the dialog for a new version, only once a day
    window, fake = QWidget(), Fake()
    assert updater.startup_check(window, fetch=fake, now=1e6, folder=Path(tempfile.gettempdir())) is True
    assert wait_for(lambda: getattr(window, "_update_dialog", None) is not None)
    assert window._update_dialog.update.version == "0.2.0" and window._update_dialog.isVisible()
    window._update_dialog.reject()
    assert updater.startup_check(window, fetch=fake, now=1e6 + 100) is False        # asked a minute ago
    assert fake.requests == [project.LATEST_RELEASE_API]

    # --- silent when up to date, skipped or offline
    forget_state()
    for i, fake in enumerate((Fake(release("v0.1.0")), Fake(error=UpdateError("offline")))):
        window = QWidget()
        assert updater.startup_check(window, fetch=fake, now=1e6 * (i + 1))
        assert wait_for(lambda: not updater._alive)
        assert not hasattr(window, "_update_dialog")
    updater.skip_version("0.2.0")
    window = QWidget()
    assert updater.startup_check(window, fetch=Fake(), now=1e7)
    assert wait_for(lambda: not updater._alive)
    assert not hasattr(window, "_update_dialog") and messages == []

    # --- nothing at all while updates are disabled
    forget_state()
    with monkeypatch.context() as patch:
        patch.setattr(project, "UPDATES_ENABLED", False)
        fake = Fake()
        assert updater.startup_check(QWidget(), fetch=fake) is False and fake.requests == []


def test_menu_action(monkeypatch, messages):
    """The menu action says 'not available yet' when disabled without touching the network, reports up to date and errors, and
    opens the dialog (also for a skipped version)."""
    # --- not available yet while disabled
    with monkeypatch.context() as patch:
        patch.setattr(project, "UPDATES_ENABLED", False)
        fake = Fake()
        assert updater.check_now(QWidget(), fetch=fake) is False
        assert fake.requests == [] and "not available yet" in messages[0] and project.RELEASES_URL in messages[0]

    # --- up to date, errors, and the dialog
    messages.clear()
    monkeypatch.setattr(project, "VERSION", "0.1.0")
    window = QWidget()
    assert updater.check_now(window, fetch=Fake(release("v0.1.0")))
    assert wait_for(lambda: messages) and "up to date" in messages[-1]
    assert updater.check_now(window, fetch=Fake(error=UpdateError("Could not reach api.github.com: offline")))
    assert wait_for(lambda: len(messages) == 2) and "offline" in messages[-1] and project.RELEASES_URL in messages[-1]
    updater.skip_version("0.2.0")                                                   # a manual check shows a skipped version again
    assert updater.check_now(window, fetch=Fake())
    assert wait_for(lambda: getattr(window, "_update_dialog", None) is not None)
    window._update_dialog.reject()


# ---- the dialog ----------------------------------------------------------------------------------------------------------------------

def make_dialog(tmp_path, fetch=None, **kw):
    log = {"installed": [], "quit": 0}
    dialog = updater.UpdateDialog(update(), None, folder=tmp_path, fetch=fetch or Fake(),
                                  install=lambda p: log["installed"].append(Path(p)),
                                  quit=lambda: log.__setitem__("quit", log["quit"] + 1), current_version="0.1.0", **kw)
    dialog.show()
    return dialog, log


def test_job_gc_handling():
    """A job pauses the cycle collector while its thread runs and restores it (back on only after the last of several jobs),
    and leaves it off if it was off before."""
    import gc
    # --- pauses and restores the collector
    assert wait_for(lambda: updater._gc_jobs == 0)           # jobs of earlier tests have delivered their answers
    gate1, gate2, got = threading.Event(), threading.Event(), []
    assert gc.isenabled()
    first = updater._Job(lambda progress: gate1.wait(5), got.append)
    second = updater._Job(lambda progress: gate2.wait(5), got.append)      # two at once: back on only after the last one
    try:
        assert not gc.isenabled()
        gate1.set()                                                           # the first one ends: the second still runs
        assert wait_for(lambda: len(got) == 1)
        assert not gc.isenabled() and updater._gc_jobs == 1, "still off while the second job runs"
    finally:
        gate1.set()
        gate2.set()
    assert wait_for(lambda: len(got) == 2) and wait_for(gc.isenabled)
    assert first is not second and updater._gc_jobs == 0

    # --- leaves the collector off if it was off before
    assert wait_for(lambda: updater._gc_jobs == 0)
    got = []
    gc.disable()
    try:
        updater._Job(lambda progress: 1, got.append)
        assert wait_for(lambda: got) and not gc.isenabled()
    finally:
        gc.enable()


def test_dialog_install_flow(tmp_path):
    """The dialog shows the notes and installs after a verified download, does not install a download that fails the checksum,
    lets the confirm stop the install, and remembers a skipped version while Later closes without remembering."""
    # --- notes, and install after a verified download
    folder = fresh_dir(tmp_path, "ok")
    dialog, log = make_dialog(folder)
    assert "0.2.0" in dialog.findChildren(type(dialog.status))[0].text() and dialog.notes.toPlainText() == "Notes"
    dialog.install_button.click()
    assert wait_for(lambda: dialog.installed)
    assert log["installed"] == [folder / f"{updater.ASSET_NAME}.new"] and log["quit"] == 1
    assert (folder / f"{updater.ASSET_NAME}.new").read_bytes() == EXE_BYTES

    # --- not installed when the checksum fails
    folder = fresh_dir(tmp_path, "tampered")
    dialog, log = make_dialog(folder, fetch=Fake(exe=b"MZ-tampered"))
    dialog.install_button.click()
    assert wait_for(lambda: not dialog.busy)
    assert log == {"installed": [], "quit": 0} and not dialog.installed
    assert "does not match" in dialog.status.text() and dialog.install_button.isEnabled() and dialog.install_button.text() == "Try Again"
    assert list(folder.iterdir()) == []

    # --- confirm can stop the install
    folder = fresh_dir(tmp_path, "confirm")
    dialog, log = make_dialog(folder, confirm=lambda: False)
    dialog.install_button.click()
    assert wait_for(lambda: not dialog.busy)
    assert log["installed"] == [] and list(folder.iterdir()) == []

    # --- skip remembers the version and closes
    dialog, log = make_dialog(tmp_path)
    dialog.skip.click()
    assert updater.skipped_version() == "0.2.0" and not dialog.isVisible() and log["installed"] == []

    # --- later closes without remembering
    forget_state()
    dialog, _ = make_dialog(tmp_path)
    dialog.later.click()
    assert updater.skipped_version() == "" and not dialog.isVisible()


def test_dialog_cancel_and_from_source(tmp_path):
    """Closing the dialog during the download cancels it and cleans up; from source the dialog cannot install and points to the
    releases page."""
    # --- closing during the download cancels it and cleans up
    gate = threading.Event()

    class Slow(Fake):
        def __call__(self, url, timeout=30, sink=None, progress=None, cancel=None):
            if sink is not None:
                sink.write(b"part")
                gate.wait(5)
            return super().__call__(url, timeout, sink, progress, cancel)
    folder = fresh_dir(tmp_path, "slow")
    dialog, log = make_dialog(folder, fetch=Slow())
    dialog.install_button.click()
    assert dialog.busy and dialog.later.text() == "Cancel"
    dialog.later.click()
    gate.set()
    assert wait_for(lambda: not updater._alive)
    assert log == {"installed": [], "quit": 0} and list(folder.iterdir()) == []

    # --- from source: cannot install, points to the releases page
    dialog = updater.UpdateDialog(update(), None, folder=tmp_path)           # no injected install: the real rule applies
    assert not dialog.install_button.isEnabled() and "source code" in dialog.status.text() and project.RELEASES_URL in dialog.status.text()
    assert dialog.skip.isEnabled() and dialog.later.isEnabled()
