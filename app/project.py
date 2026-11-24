"""What the program says about itself, in one place (the welcome screen, About, the update check, the issue links all read it).

Change a value here, not in the places that show it. VERSION must equal the version in core/Cargo.toml and core/pyproject.toml
(tests/test_project.py checks it).
"""
NAME = "openDNA"
# The file name of the released program (openDNA.exe is renamed to this in the release, and the installer is EXE_NAME-setup.exe).
# Only the downloaded file is renamed: NAME (the settings, the window titles, the data folder) stays.
EXE_NAME = "oDNA_beta"
VERSION = "0.8.5"
AUTHOR = "Antoine Champie"
CONTACT_EMAIL = ""   # a project-only address, to be filled in by the author before the first release (empty: not shown)
REPO = "OpenDNA-project/openDNA"   # owner/name on GitHub
REPO_URL = f"https://github.com/{REPO}"
RELEASES_URL = f"{REPO_URL}/releases"
ISSUES_URL = f"{REPO_URL}/issues"
LATEST_RELEASE_API = f"https://api.github.com/repos/{REPO}/releases/latest"
# The program looks for newer versions and offers to install them only when this is True. On since 0.7.0, when the repository became
# public (a private repository's releases cannot be read without a login).
UPDATES_ENABLED = True
# File extensions the program opens (File > Open, "Open with", the file-type association).
FILE_TYPES = {
    ".gb": "GenBank file", ".gbk": "GenBank file", ".gbff": "GenBank file", ".genbank": "GenBank file", ".ape": "GenBank file", ".embl": "EMBL file",
    ".dna": ".dna file",
    ".fa": "FASTA file", ".fasta": "FASTA file", ".fna": "FASTA file", ".fas": "FASTA file",
}
