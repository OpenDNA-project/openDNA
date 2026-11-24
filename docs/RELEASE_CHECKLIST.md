# Release 0.8.5: what to try by hand

State when this was written (2026-11-20, `main`, not pushed; updated after your first round of feedback): 159 Python + 71 Rust tests pass, pyflakes clean,
test report rebuilt (230 tests, all described). **No window was opened for any change below**: every GUI line was checked by reading the code only (no widget tests exist).
Open problems are in `docs/KNOWN_ISSUES.md`.

What changed is in `CHANGELOG.md` (rename `## Unreleased` to the version at release, see CLAUDE.md); this file is only what to try by hand.

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

### Your first round of feedback
- Digest on the genome: tick AciI. The warning line appears at once, no cut is drawn, the table lists 1,000 fragments; untick it: back to normal.
- Local Align: the "Aligned (bp)" column matches Query end - start + 1.
- Synthesis: open it from a pane's picker (Tools > Synthesis; it is no longer in the Tools menu). Paste a FASTA record, run the calculator, point at
  each row of its table: the stretches light up in the sequence view; click a row: they stay and the view scrolls there; edit a base: they go.
- Detect Common Features: the window is half sequence, half list; click a feature: the left view jumps to it and selects it.
  On pFG051_Chlor.gb: Tn5 ME (misc_feature) shows "annotated the other way", unticked. Annotate a promoter backwards (Edit a feature, flip it), reopen
  Detect: it is ticked as "annotated the other way"; Add turns it round, Undo puts it back.
- Help > Run Setup Again: a ViennaRNA row says whether it is installed; ticking it downloads it (needs the internet).
- Display panel on the genome: open the criteria arrow, set every enzyme and every number of sites: the orange warning appears and no site is drawn;
  set them back: sites return. The Synthesis table: three columns, no scrollbar, tooltips say what each counts.

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
`v0.8.5`, push the tag, wait for the release workflow, paste the version's section of `CHANGELOG.md` (without its `##` line) in the draft above **Downloads.**, check the files, publish (docs/PUBLISHING.md section 8).
