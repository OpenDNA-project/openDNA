# Release 0.8.5: changelog and what to try by hand

State when this was written (2026-11-20, `main` at 961fd53, not pushed): 159 Python + 71 Rust tests pass, pyflakes clean, test report rebuilt
(230 tests, all described). **No window was opened for any change below**: every GUI line was checked by reading the code only (no widget tests exist).
Open problems are in `docs/KNOWN_ISSUES.md`.

## Changelog (since 0.8.4)

### The window never freezes
- No analysis runs on the window's own thread any more. After an edit on a sequence over 50 kb, enzyme sites, primer sites, the Plasmid Checkup
  and the Enzymes/Digest table are computed in the background; until they arrive the sites are not drawn and the tool is greyed (nothing stale is shown).
- The status bar says "Running in the background: ..." for any job that takes over a second.
- Restriction sites are found in the Rust core, all enzymes in one pass: a 4.6 Mb genome takes 1.2 s instead of 7.8 s, and the 250 kb limit is gone
  (sites are shown on whole genomes). An N in the sequence now fits an N of a site, as in the search bar.
- The search bar searches long sequences in the background; autosave writes in the background; saving a large GenBank no longer blocks edits.

### Local Align and Folder Search
- Every candidate place is aligned, on several CPU cores. Before, only the 200 places with most seeds were, which lost real hits on a genome
  (a 5 kb query: 15 hits instead of 3). A 5 kb query on a 4.6 Mb genome takes 0.7 s on 8 cores.
- Optional time limit (Tools > Settings > Local Align, default 0 = none). When it is reached, a warning box says how many places were left out
  and where to change the limit.
- Every other search limit is said in the status line instead of silently cutting the answer.
- A difference in the last few bases of a query is now listed as a mismatch instead of being trimmed off the hit.
- Folder Search says which files had more than 5 hits; Copy results lists every note.

### CPU cores and setup
- New setting Tools > Settings > General > Performance: the CPU cores openDNA may use (default: all but one, at most 8). The installer asks for it.
- Help > Run Setup Again: the installer's choices for an existing install (openDNA in "Open with", CPU cores, pLannotate and the Terminal's tools
  in WSL), shown with their current state; nothing changes until OK.

### New tools (from the pre-release work)
- **Synthesis tab** (Tools, or right-click "Test Synthesis..."): the Synthesis Success Calculator (Salis lab model; needs ViennaRNA, with an install
  helper) and an IDT complexity check (login kept in the system's secure store only if you tick Remember).
- **pLannotate** in Detect Common Features: a "Found by" column, an install dialog, a credit line, and the "Also search with pLannotate" setting.
- **Tables** (Detect, Checkup, Local Align, Folder Search and others): sort by any column, filter by value, right-click a cell to show only / hide.
- Welcome page: up to 4 slots that can expire, a heading, a tip in the status bar every 5 minutes.

### Circular sequences: nothing near the origin is ignored
- Variants: clips that continue the circle are reads over the origin; a deletion through the origin is one deletion; nothing near an end is dropped.
- Golden Gate closing junction on a circle, recombineering homology through the origin, DNA and protein find through the origin, repeats across it.

### Command line for agents (`oDNA-cli`)
- New commands: `split`, `merge`, `synthesis`, `wsltools`. `align` has `--time-limit` and `--threads` and no longer crashes on a mismatch just past
  the origin. Every limit reached is a warning in every command. `msa` takes a multi-record file; `annotate` has `--no-plannotate`.

## Text for the release (in-app "What's new" is `WHATS_NEW["0.8.5"]` in `app/project.py`)

For the GitHub draft, above the **Downloads.** line:

> - The window no longer freezes: analyses run in the background, and the status bar says what is running. Restriction sites now show on whole genomes.
> - Local Align aligns every candidate place on several CPU cores and finds hits it used to miss; an optional time limit warns when it cuts the search.
> - New setting for the CPU cores openDNA may use (the installer asks), and Help > Run Setup Again for the installer's choices.
> - New Synthesis tab (synthesis success check, IDT complexity), pLannotate in Detect Common Features, sortable and filterable tables.
> - Sequences that cross the origin of a plasmid are handled everywhere (variants, assembly, find, repeats).
> - `oDNA-cli`: new `split`, `merge`, `synthesis` and `wsltools` commands; every limit reached is reported.

## Try by hand before publishing

Setup: rebuild the core and start from source (`D:\GitHub_Projects\openDNA-test.bat` does both, on `demo/BW25113_2026.gb`). Also have a plasmid with a
feature across base 1 (`demo/pFG051_Chlor.gb`, then Edit > Set Origin to Selection Start in the middle of a feature).

### The window never freezes (genome open)
1. Type a few bases in the genome. The window stays responsive; enzyme and primer sites disappear for a moment, then come back. The status bar shows
   "Running in the background: ..." if it takes over a second.
2. Show enzyme sites (Display) on the genome: sites appear on the whole 4.6 Mb (before: none above 250 kb). Open the Enzymes/Digest table: greyed while it computes, then filled.
3. Open the Plasmid Checkup on the genome, then edit: it greys and refreshes without freezing.
4. Search bar: search a 10-base motif in the genome; the window stays usable while it searches. Edit during a search: the old answer is dropped.
5. Save the genome (Ctrl+S) and keep scrolling while it writes. Close the file right after an edit: no error, no stray recovery copy at the next start.

### Local Align and Folder Search
6. Local Align on the genome with a 1-2 kb piece of pFG051 as query: hits listed, no note. Then a 5 kb query: more hits than in 0.8.4, under about a second.
7. Settings > Local Align > Time limit = 0.1 s, align the 5 kb query again: the **warning box** says how many places were left out and names the setting.
   Its text should be clear to you; set it back to 0.
8. Query with one base changed among its last 3 bases: the hit still covers the whole query and the mismatch is listed in Differences.
9. Folder Search on `demo/` with a short repeated query: the status line names files with more than 5 hits; Copy results ends with "Note: ..." lines.

### Settings and setup
10. Settings > General > Performance shows "CPU cores openDNA may use" with your core count; change it, close, reopen: kept.
11. Help > Run Setup Again: the dialog shows the current state of each choice; Cancel changes nothing; tick "Open with" and press OK (from the installed
    copy only: from source it says it does not apply).
12. **Installer** (after the build, step 2 of the release): install over 0.8.4. The new "CPU cores" page appears with min(8, cores - 1), or your saved value;
    a number outside 1..cores is refused. After the install, Settings shows that value. Let the 0.8.4 copy update itself to 0.8.5 instead: no page, your value kept.

### Earlier pre-release work (unchanged from the last checklist, still untried)
13. Start-up: no import error; status bar tip; welcome page slots.
14. Tables in Detect / Checkup / Local Align / Folder Search: sort a column, then click a row: the right feature is selected; ticked boxes after a sort add the right features.
15. Detect Common Features with and without pLannotate (install button, "Found by", credit line).
16. Synthesis tab on a 500 bp and a 2 kb piece (window stays responsive); without ViennaRNA a clear message. IDT with fake credentials: a readable error.
17. Circular origin: Find a query through base 1; protein motif through the origin; Variants pane with a deletion through the origin is one row.
18. Command line: run each recipe example once (`oDNA-cli skills show index`, then the examples); `align` with `--time-limit 0.001` gives a warning.

## Publish (when you give the go)
Version bump to 0.8.5 in `app/project.py`, `core/Cargo.toml`, `core/pyproject.toml` (+ `Cargo.lock`), tests, commit, push `main`, wait for CI, tag
`v0.8.5`, push the tag, wait for the release workflow, put the text above in the draft, check the files, publish (docs/PUBLISHING.md section 8).
