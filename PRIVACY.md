# Privacy

OpenGenEdit is an offline program. **Nothing leaves your computer unless you use one of the network features listed below, and only
what is listed.** There is no telemetry, no analytics, no crash reporting that sends anything by itself, no account and no advertising network.
The one thing that runs by itself is an optional once-a-day check for a new version (see the table; you can turn it off). Your sequence files are never uploaded.

This page describes version 0.1.0 (an early test version), including the Help-menu feedback items, the update check and the file-association question that go with it. If you find the program doing something that is not on this page, please
report it as a bug (see [CONTRIBUTING.md](CONTRIBUTING.md)); if it is a security matter, see [SECURITY.md](SECURITY.md).

## What contacts the internet, and what is sent

| You do this | It connects to | What is sent |
|---|---|---|
| **File > Import from NCBI**, press Search or Open | `eutils.ncbi.nlm.nih.gov` (NCBI) | The words or accession you typed, the name `OpenGenEdit` (as the tool name and the User-Agent), and, only if you typed it in that window, your NCBI API key. Your IP address is visible to NCBI, as for any website. |
| **3D structure > Get structure**, source *AlphaFold DB* | `alphafold.ebi.ac.uk` (EMBL-EBI) | The UniProt accession only. |
| 3D structure, source *PDB* | `search.rcsb.org` and `files.rcsb.org` (RCSB PDB) | **The protein sequence** (to search for a matching structure), then the PDB code of a hit. |
| 3D structure, source *ESMFold* | `api.esmatlas.com` (the ESM Metagenomic Atlas, run by Meta) | **The protein sequence** (up to 400 amino acids), to be predicted. |
| 3D structure, source *Auto* | In this order: your local cache, AlphaFold DB (only if there is an accession), PDB, ESMFold, stopping at the first that works | As each row above, so the sequence may go to RCSB PDB and to the ESM Atlas. |
| 3D structure, source *Local command* | Nothing, from OpenGenEdit | The command you typed runs on your own computer with a temporary FASTA file. What that command does on its own (some predictors download data) is up to it. |
| **Open 3D view in the browser** | Nothing | It opens a local HTML file in your default browser. The page contains the 3Dmol.js viewer and your structure inside itself and loads nothing from the internet. |
| **Tools > Install Mapping Tools** | conda, inside your WSL, downloads bwa, samtools and minimap2 from conda-forge / bioconda | Whatever conda sends. OpenGenEdit sends nothing itself. Only when you run that command. |
| Click a link in Help > About, Help > Citations or the welcome window | Your web browser opens it | Nothing from OpenGenEdit. |
| **Help > Report a Problem / Request a Feature** | Your web browser opens a GitHub "new issue" page | A link that pre-fills the form with the OpenGenEdit version, your Windows / Python / Qt versions and the last lines of the program's log. **Nothing is sent until you press "Submit" on GitHub**, and you can edit or delete the text first. Before it goes into the link or into the log file, DNA and protein sequences, file names, folder paths and your Windows user name are replaced by placeholders (`<sequence>`, `<file>`, `<path>`, `<user>`). GitHub then knows your account and IP address, as for any issue. |
| **Help > Copy Diagnostics**, **Open Log Folder** | Nothing | The same text goes to your clipboard; the log folder opens in Explorer. |
| **Check for Updates** (Tools menu), and once a day at start-up when *Tools > Settings > Updates* is ticked | `api.github.com` (GitHub) | A request for "the latest release of this repository" with the User-Agent `OpenGenEdit/<version>`. Your IP address is visible to GitHub; nothing about you or your files is sent. **Off in test versions while the repository is private** (the program then says "updates are not available yet"). A new exe is downloaded only if you press the button in the update window; it is fetched over HTTPS from GitHub only, checked against the SHA-256 published with the release, and refused if it does not match. |
| The question "Open .gb / .dna / FASTA files with OpenGenEdit?" (Windows exe only) | Nothing | Only if you answer Yes, it writes the registry entries listed below. |

The 3D structure window tells you, before you press *Get structure*, which service your sequence will go to. A sequence is sent only after you press the
button. If a sequence is confidential (an unpublished protein, for example), use *Local command*, or do not use the PDB / ESMFold / Auto sources: those
services' own privacy terms apply to what you send them.

Everything else works without a connection: editing, features, primers, PCR, Cas9 / Cas12a, digest, assembly, local alignment, protein view, search, tracks,
read mapping with local tools.

## What is stored on your computer

OpenGenEdit keeps its settings and caches in its **data folder**, and a few things in the Windows registry. The program's own folder is left alone (the one exception: when you press Install in the update window, the exe replaces itself).

**Data folder.** `%APPDATA%\OpenGenEdit\` (usually `C:\Users\<you>\AppData\Roaming\OpenGenEdit\`), whatever the program file is called
(the downloaded `oDNA_beta.exe`, or `python app/main.py` from source: the program names the folder itself). It can be moved by setting the
environment variable `GENOMEEDIT_DATA_DIR` to another folder before starting the program.

What is in it:

| File or folder | What it holds |
|---|---|
| `primers.tsv` | Your primer database (name, sequence, note, date added). |
| `parts.tsv` | Your own part library, if you imported parts (the program's built-in library is inside the program). |
| `settings.ini`, `theme.ini`, `workspace.ini`, `window.ini`, `image_export.ini`, `recent.ini`, `search.ini` | Tool settings, light / dark theme, window and tab layout, window size and place, side-panel width, View switches and column widths, image export choices, "snooze the welcome window" date, recently opened file paths, search options. |
| `ncbi.ini` | The NCBI window's options: your **NCBI API key if you typed it, stored as plain text** (the field is masked on screen only). |
| `structure.ini` | The 3D structure window's options, including the text of the *Local command*, if you typed one. |
| `primers.ini` | PCR conditions (salt, Mg2+, dNTP, primer concentration) and the primer options, including the path of the reference-genome file you chose. |
| `updates.ini`, `file_assoc.ini` | When the last update check ran and which version you chose to skip; your answer to the file-association question. |
| `logs\` | `genomeedit.log` (3 files of at most 500 KB): start-up lines with versions, warnings and errors, with sequences, file names, paths and user name replaced by placeholders. It stays on your computer; you may attach it to an issue. |
| `ncbi\` | Records downloaded from NCBI, so asking again does not download again. |
| `structures\` | Cached 3D structures (named by a hash of the sequence) and the last 3D page written. |

**Windows registry.** The program's settings are not in the registry: the PCR conditions (salt, Mg2+, dNTP, primer concentration) and the primer options
(mismatches allowed, homology seed, and the path of the reference-genome file you chose for non-specific priming checks) are in `primers.ini` in the data
folder, like the other `.ini` files. **Only if you answer Yes** to the file-association question (or ask for it again later), the exe writes, under
`HKEY_CURRENT_USER\Software\Classes` and `HKEY_CURRENT_USER\Software\RegisteredApplications`, the entries that make "OpenGenEdit" appear in Windows' "Open with"
list and in Settings > Default apps (no administrator rights; Windows still lets only you choose the default program).

**Next to your own files.** While you edit, a recovery copy is written every minute beside the file as `<file>.autosave.gb`; it is deleted when you save or close
normally (for a `.dna` file it sits beside it, and the `.dna` file itself is never written). *Map Reads* writes a `<output>_run\` folder
beside the output you chose (the reference, the exact commands and a log). A *Local command* uses a temporary folder named `genomeedit_fold_...` in your
temp folder and removes it afterwards.

## How to delete it

Close OpenGenEdit first. In the **uninstall** sense there is nothing else: no service, no scheduled task, no start-up entry.

- **Everything in the data folder** (primers, parts, settings, caches, the stored NCBI key): delete the folder. In PowerShell, for the downloaded exe:

  ```powershell
  Remove-Item -Recurse -Force "$env:APPDATA\OpenGenEdit"
  ```

  (from source: `"$env:APPDATA\python"` is shared by every program that runs under that name, so delete only the files and folders listed above, not the whole folder).
  **This deletes your primer database**: export it first (Primers > Primer Database > Export) if you want to keep it.
- **Only the caches**: delete the `ncbi` and `structures` folders inside it. Delete `ncbi.ini` to forget the NCBI key.
- **The registry entries**:

  ```powershell
  Remove-Item -Recurse "HKCU:\Software\OpenGenEdit"
  Remove-Item -Recurse "HKCU:\Software\Classes\OpenGenEdit.SequenceFile" -ErrorAction SilentlyContinue
  Remove-ItemProperty "HKCU:\Software\RegisteredApplications" -Name OpenGenEdit -ErrorAction SilentlyContinue
  ```

  The file-association entries also sit under `HKCU:\Software\Classes\Applications\oDNA_beta.exe` and, one value each, in `...\Classes\<.gb, .dna, ...>\OpenWithProgids`; the
  complete list is at the top of `app/file_assoc.py`. They are harmless if left (they only name the program in the "Open with" list).
- **Recovery copies** (`*.autosave.gb`) and `*_run` folders: delete them where they are, next to your files.
- The program itself is a single file with no installer: delete `oDNA_beta.exe` (and a leftover `oDNA_beta.exe.old` beside it after an update).

## Your files

OpenGenEdit reads and writes only the files you open, save or export, the autosave copies and the run folders described above. It never modifies a `.dna`
file (edits are saved as GenBank). Keep copies of your originals anyway: this is an early test version.

## Questions

Open an issue on GitHub (https://github.com/Ant-Champ/GenomeEdit/issues/new/choose), or write to `<project email>`.
