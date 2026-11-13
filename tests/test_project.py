"""app/project.py (name, version, links) agrees with the build files, and the exe recipe bundles every data file the app reads.

A release bumps the version in three places (project.py, core/Cargo.toml, core/pyproject.toml); the first test fails if one is forgotten.
Also: the data folder of the program's old name is copied over once.
"""
import os
import re
import struct
import sys
import tempfile
import tomllib
from pathlib import Path

os.environ["OPENDNA_DATA_DIR"] = tempfile.mkdtemp(prefix="opendna_test_")
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "app"))

import project

SPEC = (ROOT / "packaging" / "openDNA.spec").read_text(encoding="utf-8")


def test_release_metadata_and_exe_recipe_are_consistent():
    """project.py agrees with the build files: the version is the same in all three places and is major.minor.patch, the links
    are built from the repository name, the name and author are filled in, the updates switch is a bool, every file type
    extension starts with a dot. The exe recipe bundles every data file the app reads (Path(__file__).with_name("x") only works
    in the frozen exe when the spec ships x next to the modules), the files it names exist, and the icon is a multi-size ico."""
    # --- version is the same in all three places ---
    cargo = tomllib.loads((ROOT / "core" / "Cargo.toml").read_text(encoding="utf-8"))["package"]["version"]
    pyproject = tomllib.loads((ROOT / "core" / "pyproject.toml").read_text(encoding="utf-8"))["project"]["version"]
    assert project.VERSION == cargo == pyproject, (project.VERSION, cargo, pyproject)
    # --- version looks like a release number ---
    # The exe's file version (packaging/openDNA.spec) and the release tag (v + VERSION) both rely on major.minor.patch.
    assert re.fullmatch(r"\d+\.\d+\.\d+", project.VERSION), project.VERSION
    # --- links are built from the repository name ---
    owner, name = project.REPO.split("/")
    assert owner and name and "/" not in name
    assert project.REPO_URL == f"https://github.com/{project.REPO}"
    assert project.ISSUES_URL == project.REPO_URL + "/issues"
    assert project.RELEASES_URL == project.REPO_URL + "/releases"
    assert project.LATEST_RELEASE_API == f"https://api.github.com/repos/{project.REPO}/releases/latest"
    for url in (project.REPO_URL, project.ISSUES_URL, project.RELEASES_URL, project.LATEST_RELEASE_API):
        assert url.startswith("https://") and " " not in url
    # --- name and author are filled in ---
    assert project.NAME == "openDNA"      # the exe, the settings keys and the release asset are all called this
    assert project.AUTHOR.strip()
    # --- updates switch is a bool ---
    assert isinstance(project.UPDATES_ENABLED, bool)
    # --- every file type extension starts with a dot ---
    assert project.FILE_TYPES
    for ext, label in project.FILE_TYPES.items():
        assert ext.startswith(".") and len(ext) > 1 and ext == ext.lower() and " " not in ext, ext
        assert label.strip(), ext

    # --- every data file the app reads is bundled ---
    wanted = set()
    for py in (ROOT / "app").glob("*.py"):
        text = py.read_text(encoding="utf-8")
        wanted.update(re.findall(r'Path\(__file__\)\.with_name\("([^"]+)"\)', text))
    assert {"parts", "sponsors.json", "tools.toml", "vendor"} <= wanted, wanted   # the scan itself still finds them
    for name in wanted:
        assert f'"{name}"' in SPEC, f"app reads {name} but packaging/openDNA.spec does not bundle it"
    # Run the spec's data section (everything before the version info, which needs PyInstaller) and read its real `datas`:
    # a name only in a comment or an excludes list is no longer enough. Each file the app reads must land beside the modules
    # (its own name, folder ".") or be the folder a list entry is copied into.
    head, marker, _ = SPEC.partition('# Windows "Details"')
    assert marker, "the spec's data section is no longer followed by the version info: update this split"
    ns = {"SPECPATH": str(ROOT / "packaging"), "__name__": "spec_data"}
    exec(compile(head, "openDNA.spec", "exec"), ns)
    for name in wanted:
        assert any((Path(src).name == name and dest == ".") or dest == name for src, dest in ns["datas"]), \
            f"app reads {name} but packaging/openDNA.spec's datas does not bundle it"
    # The protein viewer's library is reached through the vendor folder.
    assert 'APP / "vendor" / "3Dmol-min.js"' in SPEC and 'APP / "vendor" / "3Dmol.LICENSE"' in SPEC

    # --- files named in the exe recipe exist ---
    for rel in re.findall(r'\(APP / "([^"]+)"', SPEC) + re.findall(r'\(ROOT / "([^"]+)"', SPEC):
        assert (ROOT / "app" / rel).exists() or (ROOT / rel).exists(), rel
    for rel in ("app/parts/ncbi_vectors.tsv", "app/tools.toml", "app/sponsors.json", "app/vendor/3Dmol-min.js", "app/vendor/3Dmol.LICENSE", "LICENSE"):
        assert (ROOT / rel).exists(), rel

    # --- icon is a multi-size ico ---
    data = (ROOT / "packaging" / "openDNA.ico").read_bytes()
    reserved, kind, count = struct.unpack_from("<HHH", data)
    assert (reserved, kind) == (0, 1)
    sizes = set()
    for i in range(count):
        w, h, _, _, _, bpp, length, offset = struct.unpack_from("<BBBBHHII", data, 6 + 16 * i)
        sizes.add(w or 256)
        assert w == h and bpp == 32 and offset + length <= len(data)
    assert {16, 32, 48, 256} <= sizes, sizes


def test_old_data_folder_is_copied_once(tmp_path):
    """The settings and primers of the program's old name come along to the new data folder, once, and the old folder stays."""
    from primer_db import OLD_DATA_NAME, copy_old_data
    old, new = tmp_path / OLD_DATA_NAME, tmp_path / "openDNA"
    (old / "logs").mkdir(parents=True)
    (old / "primers.tsv").write_text("name\tsequence\n", encoding="utf-8")
    assert copy_old_data(new) is True
    assert (new / "primers.tsv").read_text(encoding="utf-8") == "name\tsequence\n" and (new / "logs").is_dir()
    assert (old / "primers.tsv").exists(), "the old folder is left as it was"
    (new / "primers.tsv").write_text("changed", encoding="utf-8")
    assert copy_old_data(new) is False and (new / "primers.tsv").read_text(encoding="utf-8") == "changed", "never copies over the new folder"
    assert copy_old_data(tmp_path / "elsewhere" / "openDNA") is False, "nothing to copy when there is no old folder"
