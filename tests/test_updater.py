"""Self-update (app/updater.py): which releases are offered, where they may come from, the checksum, install and rollback on fake
files, and the once-a-day throttle. The network and the relaunch are replaced by fakes: nothing is downloaded or started."""
import hashlib
import json
import os
import platform
import sys
import tempfile
import threading
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ["OPENDNA_DATA_DIR"] = tempfile.mkdtemp(prefix="opendna_test_")
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "app"))

import pytest

import project
import tool_settings
import updater
from primer_db import data_dir
from updater import Update, UpdateError

EXE_BYTES = b"MZ" + b"fake exe " * 1000
DIGEST = hashlib.sha256(EXE_BYTES).hexdigest()
BASE = f"https://github.com/{project.REPO}/releases/download/v0.2.0"


def release(tag="v0.2.0", assets=True, **extra):
    """The JSON GitHub gives for releases/latest."""
    data = {"tag_name": tag, "html_url": f"{project.REPO_URL}/releases/tag/{tag}", "body": "Fixes things.", "draft": False, "prerelease": False,
            "assets": [{"name": updater.release_asset(tag), "browser_download_url": f"{BASE}/openDNA.exe"},
                       {"name": updater.release_asset(tag) + ".sha256", "browser_download_url": f"{BASE}/openDNA.exe.sha256"},
                       {"name": "other.zip", "browser_download_url": f"{BASE}/other.zip"}] if assets else []}
    data.update(extra)
    return json.dumps(data)


class Fake:
    """A fake of updater.fetch: the release JSON, the checksum and the exe by URL; remembers the requests."""

    def __init__(self, latest=None, checksum=None, exe=EXE_BYTES, error=None):
        self.latest = release() if latest is None else latest
        self.checksum = f"{DIGEST}  openDNA.exe\n" if checksum is None else checksum
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
    """The remembered update state and the setting emptied."""
    (data_dir() / "updates.ini").unlink(missing_ok=True)
    store = tool_settings.store()
    store.remove("Updates/check")
    store.sync()
    tool_settings._cache.clear()


@pytest.fixture(autouse=True)
def sandbox(monkeypatch):
    """No waiting between retries; the remembered state and the setting empty before and after; updates switched on."""
    monkeypatch.setattr(updater, "RETRY_DELAY", 0)
    monkeypatch.setattr(project, "UPDATES_ENABLED", True)
    forget_state()
    yield
    forget_state()


def update(**kw):
    values = dict(version="0.2.0", tag="v0.2.0", asset_url=f"{BASE}/openDNA.exe", sha256_url=f"{BASE}/openDNA.exe.sha256",
                  notes="Notes", page_url=f"{project.REPO_URL}/releases/tag/v0.2.0")
    values.update(kw)
    return Update(**values)


def fresh_dir(parent, name):
    folder = parent / name
    folder.mkdir()
    return folder


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


def test_release_check_offers_only_a_newer_complete_release(monkeypatch):
    """Only a newer, final release with both its exe and checksum on our GitHub repository is offered; versions compare
    numerically; anything else (older, draft, prerelease, missing asset, asset elsewhere, garbage, network error) offers nothing
    and never raises."""
    # --- version comparison and parsing
    for candidate, current, newer in [
        ("v0.10.0", "0.9.0", True), ("0.2.0", "0.1.0", True), ("v1.0.0", "0.99.99", True), ("0.1.1", "0.1.0", True),
        ("1.0", "1.0.0", False), ("1.0.0.1", "1.0", True), ("0.1.0", "0.1.0", False), ("v0.1.0", "0.1.0", False),
        ("0.9.0", "0.10.0", False), ("0.0.9", "0.1.0", False),
        ("0.2.0-beta1", "0.1.0", False), ("latest", "0.1.0", False), ("", "0.1.0", False), (None, "0.1.0", False), ("0.2.0", "dev", False),
    ]:
        assert updater.is_newer(candidate, current) is newer, f"{candidate!r} vs {current!r}"
    assert updater.parse_version("v0.10.0") == (0, 10, 0) and updater.parse_version(" V2.1 ") == (2, 1)
    assert updater.parse_version("1.2.3-rc1") is None and updater.parse_version("1.2.3+build") is None

    # --- a newer release gives the update
    result = updater.check("0.1.0", Fake())
    u = result.update
    assert result.error is None and u.version == "0.2.0" and u.tag == "v0.2.0" and u.notes == "Fixes things."
    assert u.asset_url == f"{BASE}/openDNA.exe" and u.sha256_url == f"{BASE}/openDNA.exe.sha256"
    assert u.page_url == f"{project.REPO_URL}/releases/tag/v0.2.0"
    assert updater.check_for_update("0.1.0", Fake()) == u

    # --- the releases API is asked with the project version by default; the module-level fetch is used when none is given
    fake = Fake(release(tag=f"v{project.VERSION}"))
    assert updater.check(fetch=fake).update is None and fake.requests == [project.LATEST_RELEASE_API]
    monkeypatch.setattr(updater, "fetch", Fake())
    assert updater.check_for_update("0.1.0").version == "0.2.0"

    # --- up to date, older, draft and prerelease offer nothing; a tag that is not a version is reported
    for tag, extra in (("v0.1.0", {}), ("v0.0.9", {}), ("v0.3.0", {"draft": True}), ("v0.3.0", {"prerelease": True}),
                       ("v0.3.0-beta1", {})):
        result = updater.check("0.1.0", Fake(release(tag, **extra)))
        assert result.update is None, tag
        assert (result.error is not None and "not understood" in result.error) == (tag == "v0.3.0-beta1"), (tag, result.error)
    result = updater.check("0.1.0", Fake(release("nightly")))
    assert result.update is None and "nightly" in result.error

    # --- without the exe or its checksum it cannot be installed
    only_exe = json.loads(release())
    only_exe["assets"] = only_exe["assets"][:1]
    only_sum = json.loads(release())
    only_sum["assets"] = only_sum["assets"][1:2]
    for data, missing in ((json.loads(release(assets=False)), updater.release_asset("v0.2.0")),
                          (only_exe, updater.release_asset("v0.2.0") + ".sha256"), (only_sum, updater.release_asset("v0.2.0"))):
        result = updater.check("0.1.0", Fake(json.dumps(data)))
        assert result.update is None and missing in result.error and result.page_url, missing

    # --- assets on other hosts or other repositories are refused; a page link off GitHub falls back to the releases page
    for url in ("https://evil.example.com/openDNA.exe", "http://github.com/x/openDNA.exe",
                "https://github.com/someone/else/releases/download/v9/openDNA.exe"):
        data = json.loads(release())
        data["assets"][0]["browser_download_url"] = url
        result = updater.check("0.1.0", Fake(json.dumps(data)))
        assert result.update is None and "Refused" in result.error, url
    # the checksum file is held to the same rule: a checksum from elsewhere could vouch for a tampered exe
    for url in ("https://evil.example.com/x.sha256", "http://github.com/x.sha256",
                "https://github.com/someone/else/releases/download/v9/x.sha256"):
        data = json.loads(release())
        data["assets"][1]["browser_download_url"] = url
        result = updater.check("0.1.0", Fake(json.dumps(data)))
        assert result.update is None and "Refused" in result.error, url
    assert updater.check("0.1.0", Fake(release(html_url="https://evil.example.com/x"))).update.page_url == project.RELEASES_URL

    # --- network errors and garbage never raise
    assert updater.check("0.1.0", Fake(error=UpdateError("Could not reach api.github.com: offline"))).error.endswith("offline")
    assert "failed" in updater.check("0.1.0", Fake(error=RuntimeError("boom"))).error
    for junk in ("not json", "[1, 2]", "", b"\xff\xfe\x00"):
        result = updater.check("0.1.0", Fake(junk))
        assert result.update is None and "understood" in result.error
    assert updater.check_for_update("0.1.0", Fake(error=UpdateError("x"))) is None


def test_download_comes_only_from_github_and_is_kept_only_when_it_matches_its_checksum(tmp_path, monkeypatch):
    """The host allow-list (https GitHub hosts only: no look-alikes, credentials, ports or other schemes), enforced by the real
    fetch before any network and on redirects; checksum file forms read or refused; the download verifies and returns the new
    exe, refuses a mismatch, a missing or unreadable checksum and addresses off GitHub; a failed or cancelled download leaves no
    partial file."""
    # --- the host allow-list
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
    # ... and fetch() really uses that guard: the opener it builds carries _SafeRedirects (no network: the opener is recorded, then stopped)
    seen = []

    def recorder(*handlers):
        seen.extend(handlers)
        raise updater.urllib.error.URLError("stopped by the test")
    with monkeypatch.context() as patch:
        patch.setattr(updater.urllib.request, "build_opener", recorder)
        with pytest.raises(UpdateError):
            updater.fetch("https://github.com/a/b/releases/download/v1/openDNA.exe")
    assert any(h is updater._SafeRedirects or isinstance(h, updater._SafeRedirects) for h in seen), seen

    # --- checksum file forms
    assert updater.parse_checksum(f"{DIGEST}  {updater.ASSET_NAME}\n") == DIGEST
    assert updater.parse_checksum(f"{DIGEST} *{updater.ASSET_NAME}") == DIGEST
    assert updater.parse_checksum(DIGEST + "\n") == DIGEST
    assert updater.parse_checksum(f"{DIGEST.upper()}  {updater.ASSET_NAME}") == DIGEST
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

    # --- a checksum mismatch is refused and the file deleted
    folder = fresh_dir(tmp_path, "mismatch")
    with pytest.raises(UpdateError, match="does not match"):
        updater.download_update(update(), folder, fetch=Fake(exe=EXE_BYTES + b"tampered"))
    assert list(folder.iterdir()) == []

    # --- a missing or unreadable checksum is refused (nothing downloaded without one)
    folder = fresh_dir(tmp_path, "nosum")
    fake = Fake()
    with pytest.raises(UpdateError, match="no checksum"):
        updater.download_update(update(sha256_url=""), folder, fetch=fake)
    assert not any(u.endswith(".exe") for u in fake.requests)
    with pytest.raises(UpdateError, match="expected form"):
        updater.download_update(update(), folder, fetch=Fake(checksum="<html>404</html>"))
    assert list(folder.iterdir()) == []

    # --- addresses off GitHub are refused before any request
    folder = fresh_dir(tmp_path, "offgithub")
    for bad in (update(asset_url="https://evil.example.com/openDNA.exe"), update(sha256_url="http://github.com/x.sha256")):
        fake = Fake()
        with pytest.raises(UpdateError, match="Refused"):
            updater.download_update(bad, folder, fetch=fake)
        assert fake.requests == []

    # --- a cancelled or failed download leaves no partial file
    folder = fresh_dir(tmp_path, "cancelled")
    stop = threading.Event()
    stop.set()
    with pytest.raises(UpdateError, match="cancelled"):
        updater.download_update(update(), folder, fetch=Fake(), cancel=stop)
    assert list(folder.iterdir()) == []

    def broken(url, timeout=30, sink=None, progress=None, cancel=None):
        if sink is None:
            return DIGEST
        sink.write(b"half")
        raise UpdateError("Could not reach github.com: reset")
    with pytest.raises(UpdateError, match="Could not reach"):
        updater.download_update(update(), folder, fetch=broken)
    assert list(folder.iterdir()) == []

    # --- a missing folder is a readable error
    with pytest.raises(UpdateError, match="Could not write"):
        updater.download_update(update(), tmp_path / "nowhere", fetch=Fake())


def test_install_swaps_relaunches_and_cleans_up(tmp_path, monkeypatch):
    """A successful install: the new exe replaces the running one (the old kept as .old, a stale .old replaced), the running exe
    is used when frozen, the relaunched child gets a clean PyInstaller environment, and the leftovers are removed later (waiting
    for a held file, never from source)."""
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

    # --- uses the running exe when frozen
    current, new = installed_in(tmp_path / "running")
    with monkeypatch.context() as frozen:
        frozen.setattr(sys, "frozen", True, raising=False)
        frozen.setattr(sys, "platform", "win32")
        frozen.setattr(sys, "executable", str(current))
        updater.install_update(new, relaunch=lambda p: None)
    assert current.read_bytes() == b"NEW"

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

    # --- does nothing from source, and waits for a held file
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


def test_install_refuses_or_rolls_back_leaving_the_old_program(tmp_path, monkeypatch):
    """Install refuses from source, as the folder build (pointing to the right download elsewhere), with missing files or a file
    in another folder; when the swap fails (no write access, antivirus that never lets go, new program does not start) the old
    program is kept intact; a briefly held file is retried."""
    # --- refuses from source
    assert not updater.is_frozen()
    _, new = installed_in(tmp_path / "source")
    with pytest.raises(UpdateError, match="source code"):
        updater.install_update(new, relaunch=lambda p: pytest.fail("must not start anything"))
    assert new.read_bytes() == b"NEW"
    assert "source code" in updater.install_blocker() and updater.install_blocker(Path("x.exe")) is None

    # --- the folder build and other systems do not replace themselves
    current, new = installed_in(tmp_path / "folderbuild")
    with monkeypatch.context() as frozen:
        frozen.setattr(sys, "frozen", True, raising=False)
        frozen.setattr(sys, "executable", str(current))
        frozen.setattr(sys, "_MEIPASS", str(tmp_path / "somewhere" / "_MEI1234"), raising=False)   # one-file: a temp folder
        for system, machine, download in (("darwin", "arm64", "macos-arm64.zip"), ("darwin", "x86_64", "macos-x86_64.zip"),
                                          ("linux", "x86_64", "linux-x86_64.tar.gz")):   # never the Windows exe there
            frozen.setattr(sys, "platform", system)
            frozen.setattr(platform, "machine", lambda m=machine: m)
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

    # --- the new file cannot be moved in: rolled back
    current, new = installed_in(tmp_path / "movein")
    with monkeypatch.context() as patch:
        fail_replace(patch, lambda s, d: s == new)
        with pytest.raises(UpdateError, match="administrator"):
            updater.install_update(new, current, relaunch=lambda p: pytest.fail("must not start"))
    assert current.read_bytes() == b"OLD" and new.read_bytes() == b"NEW" and not (current.parent / "openDNA.exe.old").exists()

    # --- a folder that cannot be written: retried, then points to the releases page
    current, new = installed_in(tmp_path / "readonly")
    with monkeypatch.context() as patch:
        calls = fail_replace(patch, lambda s, d: s == current)
        with pytest.raises(UpdateError) as e:
            updater.install_update(new, current, relaunch=lambda p: None)
    assert project.RELEASES_URL in str(e.value) and "administrator" in str(e.value)
    assert len(calls) == updater.RETRIES
    assert current.read_bytes() == b"OLD" and new.read_bytes() == b"NEW"

    # --- an antivirus holding the file briefly: retried until it lets go
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

    # --- the new program does not start: rolled back
    current, new = installed_in(tmp_path / "nostart")

    def broken(path):
        raise OSError("not a valid Win32 application")
    with pytest.raises(UpdateError, match="old one was kept"):
        updater.install_update(new, current, relaunch=broken)
    assert current.read_bytes() == b"OLD" and new.read_bytes() == b"NEW" and not (current.parent / "openDNA.exe.old").exists()


def test_automatic_check_obeys_the_switches_runs_once_a_day_and_remembers_a_skip(monkeypatch):
    """The automatic check obeys the setting and the project flag, runs at most once a day (also when the clock goes back), and a
    skipped version is remembered."""
    # --- the setting and the project flag
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
    assert updater.auto_check_due(now)
    updater.mark_checked(now)
    assert not updater.auto_check_due(now + 60) and not updater.auto_check_due(now + 23 * 3600)
    assert updater.auto_check_due(now + 24 * 3600)
    assert updater.auto_check_due(now - 3600)                             # the clock went back: do not wait for years

    # --- a skipped version is remembered
    forget_state()
    assert updater.skipped_version() == ""
    updater.skip_version("0.2.0")
    assert updater.skipped_version() == "0.2.0"

    # --- the start-up check: offers a new version once a day, never a skipped one (run synchronously, the dialog recorded)
    shown = []
    monkeypatch.setattr(updater, "cleanup_old_in_background", lambda: None)
    monkeypatch.setattr(updater, "check_in_background", lambda got, fetch=None, **kw: got(updater.check("0.1.0", fetch)))
    monkeypatch.setattr(updater, "open_dialog", lambda window, update, **kw: shown.append(update.version))
    now = 2_000_000.0
    forget_state()
    assert updater.startup_check(None, fetch=Fake(release()), now=now) is True and shown == ["0.2.0"]
    assert updater.startup_check(None, fetch=Fake(release()), now=now + 60) is False and shown == ["0.2.0"], "once a day"
    forget_state()
    updater.skip_version("0.2.0")
    assert updater.startup_check(None, fetch=Fake(release()), now=now) is True and shown == ["0.2.0"], "a skipped version is not offered"
    forget_state()
    updater.skip_version("0.1.5")                               # skipping an older version does not hide a newer one
    assert updater.startup_check(None, fetch=Fake(release()), now=now) is True and shown == ["0.2.0", "0.2.0"]
