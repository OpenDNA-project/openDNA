"""What the program says about itself, in one place (the welcome screen, About, the update check, the issue links all read it).

Change a value here, not in the places that show it. VERSION must equal the version in core/Cargo.toml and core/pyproject.toml
(tests/test_project.py checks it).
"""
NAME = "openDNA"
# The file name of the released program (openDNA.exe is renamed to this in the release, and the folder zip is EXE_NAME-folder.zip).
# Only the downloaded file is renamed: NAME (the settings, the window titles, the data folder) stays.
EXE_NAME = "oDNA_beta"
VERSION = "0.6.3"
AUTHOR = "Antoine Champie"
CONTACT_EMAIL = ""   # a project-only address, to be filled in by the author before the first release (empty: not shown)
REPO = "OpenDNA-project/openDNA"   # owner/name on GitHub
REPO_URL = f"https://github.com/{REPO}"
RELEASES_URL = f"{REPO_URL}/releases"
ISSUES_URL = f"{REPO_URL}/issues"
LATEST_RELEASE_API = f"https://api.github.com/repos/{REPO}/releases/latest"
# The program looks for newer versions and offers to install them only when this is True. The repository is private for now (a private
# repository's releases cannot be read without a login), so it stays False until the repository is made public: then switch it on.
UPDATES_ENABLED = False
# What changed in each version, newest first: the start page shows the lines of VERSION once, at the first start after an update.
WHATS_NEW = {
    "0.6.3": ["New forward / reverse primer from the selected bases with Ctrl+R / Ctrl+Shift+R.",
              "Local Align laid out as Query, Search in and results; the buttons that run a tool share one blue."],
    "0.6.2": ["The updater finds releases whose files carry the version in their name."],
    "0.6.1": ["Database primers are hidden by default when a sequence opens (Display > Primers > From the database).",
              "Released files carry the version in their name; the old openDNA file-association entries are cleaned up and follow the program when it moves."],
    "0.6.0": ["Restriction and ligation cloning, codon tools, degenerate primers and oligo order export.",
              "Multiple sequence alignment view, repeats dot plot, library and folder search.",
              "Safer file writes and a settings dialog."],
    "0.5.0": ["A cleaner, steadier interface: buttons and toolbars no longer move or change size, tables fit what they hold.",
              "Menus in the usual order (File, Edit, View...), Cut / Copy / Paste in Edit, Undo and Redo say what they undo.",
              "The window remembers its size, the View switches and the column widths you set.",
              "One range box in the toolbar (101..200), the background genome box, View > Theme, Help > Troubleshooting.",
              "Status colours readable in light and dark; images export with a light background by default."],
    "0.1.0": ["First test version."],
}
# File extensions the program opens (File > Open, "Open with", the file-type association).
FILE_TYPES = {
    ".gb": "GenBank file", ".gbk": "GenBank file", ".gbff": "GenBank file", ".genbank": "GenBank file",
    ".dna": ".dna file",
    ".fa": "FASTA file", ".fasta": "FASTA file", ".fna": "FASTA file", ".fas": "FASTA file",
}
