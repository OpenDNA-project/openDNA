"""ViennaRNA on demand. The Synthesis tab folds DNA hairpins with ViennaRNA (Lorenz et al., Algorithms Mol Biol 2011), whose licence does not
let openDNA redistribute it for a fee and asks commercial products to ask its authors, so it is not bundled. This installs it on the user's
request: the official wheel for this Python and this system is fetched from PyPI (https://pypi.org/project/ViennaRNA/), checked against the
SHA-256 PyPI lists for it and unpacked into the openDNA data folder (vendor/viennarna). Nothing of the user's sequences is sent. A copy
installed with `pip install ViennaRNA` is used just the same."""
import hashlib
import importlib
import io
import json
import platform
import shutil
import sys
import zipfile

import structure
from primer_db import data_dir
from user_error import UserError

INDEX = "https://pypi.org/pypi/ViennaRNA/json"


def folder():
    return data_dir() / "vendor" / "viennarna"


def activate():
    """Make the copy in the data folder importable (a no-op when there is none)."""
    path = folder()
    if path.is_dir() and str(path) not in sys.path:
        sys.path.insert(0, str(path))
        importlib.invalidate_caches()


def wheel_for(files, python=None, system=None, machine=None):
    """The entry of PyPI's `urls` list (a wheel) that fits this Python and system; UserError, with what was looked for, when there is none."""
    python = python or f"cp{sys.version_info.major}{sys.version_info.minor}"
    system = system or sys.platform
    machine = (machine or platform.machine()).lower()
    arm = machine in ("arm64", "aarch64")
    for f in files:
        if f.get("packagetype") != "bdist_wheel" or f.get("python_version") != python:
            continue
        tag = f["filename"][:-4].split("-")[-1]
        if system == "win32" and tag == ("win_arm64" if arm else "win_amd64"):
            return f
        if system == "darwin" and "macosx" in tag and tag.endswith("arm64" if arm else "x86_64"):
            return f
        if system.startswith("linux") and "manylinux" in tag and ("aarch64" in tag if arm else "x86_64" in tag):
            return f
    raise UserError(f"PyPI has no ViennaRNA wheel for Python {python} on {system} ({machine}). Install it yourself: pip install ViennaRNA")


def install():
    """Download and unpack ViennaRNA; returns its version. UserError when it cannot be reached, does not match its checksum or does not load."""
    try:
        release = json.loads(structure.fetch(INDEX, timeout=60))
        wheel = wheel_for(release["urls"])
        blob = structure.fetch(wheel["url"], timeout=300)
    except (structure.StructureError, KeyError, ValueError) as e:
        raise UserError(f"Could not download ViennaRNA from PyPI: {e}") from None
    if hashlib.sha256(blob).hexdigest() != wheel["digests"]["sha256"]:
        raise UserError("The ViennaRNA download does not match the checksum PyPI lists for it, so it was not installed.")
    target, staging = folder(), folder().with_name("viennarna.new")
    shutil.rmtree(staging, ignore_errors=True)
    with zipfile.ZipFile(io.BytesIO(blob)) as z:
        z.extractall(staging)   # (extractall drops any leading ../ or / of a name, so nothing lands outside the folder)
    shutil.rmtree(target, ignore_errors=True)
    staging.rename(target)
    activate()
    sys.modules.pop("RNA", None)
    try:
        import RNA
        if not hasattr(RNA, "params_load_DNA_Mathews2004"):   # (the call the Synthesis tab needs)
            raise ImportError("this ViennaRNA has no DNA parameters")
    except ImportError as e:
        shutil.rmtree(target, ignore_errors=True)
        raise UserError(f"ViennaRNA was downloaded but cannot be loaded here: {e}") from None
    return release["info"]["version"]
