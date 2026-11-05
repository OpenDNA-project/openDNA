# Feedback and roadmap (for the owner)

How testers send problems and ideas, how you sort them, and how everyone can see what happens next. Everything lives in the
GitHub repository: issues are the single list, and **labels** say what each one is and where it stands. There is no project board
on purpose: for one maintainer, a filtered list of issues is enough.

## How testers report

In the program, **Help** has four items (they are in `app/feedback.py`):

| Help item | What it does |
|---|---|
| Report a Problem... | Opens GitHub's *Report a problem* form in the browser, already filled with the version, the system (Windows, Python, Qt) and the last 30 lines of the log. The tester only describes what happened. |
| Request a Feature... | Opens the *Request a feature* form (the problem, the wanted behaviour, who it helps, alternatives). |
| Open Log Folder | Opens the folder with `opendna.log` (3 files of at most 500 KB, `logs` inside the program's data folder). |
| Copy Diagnostics | Copies versions plus the last log lines to the clipboard, to paste into an email or a chat. |

When something unexpected happens, the program does not close: it logs the full error and shows one small window,
*Something went wrong. Report it?* with **Report on GitHub**, **Copy details** and **Close**. The same error is shown only once
per run.

**Privacy.** Nothing is sent by the program itself; a report is a link that the tester opens and submits on GitHub. Before
anything is put in a link or written in the log, sequences (DNA and protein), file names, file paths and the Windows user name
are replaced by `<sequence>`, `<file>`, `<path>` and `<user>`. The forms remind testers not to attach private sequences.

A tester without the program can use *New issue* on GitHub: the forms are the same (blank issues are off, so every issue has
the same fields and labels). The testers need a free GitHub account **and access to the repository while it is private**:
add them under Settings > Collaborators (read access is enough to open issues; they cannot see the code without it).

## The labels (the public roadmap)

Defined in `.github/labels.yml`, already created in the repository. Each issue has one `type:` and one `status:` label; the
forms add `type: ...` and `status: new` by themselves.

| Label | Meaning |
|---|---|
| `type: bug`, `type: feature`, `type: question` | What it is |
| `status: new` | Nobody has looked yet (every issue starts here) |
| `status: under review` | You are reading it; you may ask the tester for details |
| `status: planned` | Accepted, not started |
| `status: in progress` | Being worked on |
| `status: done (next release)` | Finished on `main`, will be in the next release |
| `status: released` | Part of a published version |
| `status: declined` | Will not be done; always say why in a comment |
| `good first issue`, `help wanted` | For future contributors |

To re-create or change the labels later: edit `labels.yml`, then run, for each entry (Git Bash or PowerShell):

```
gh label create "status: planned" --color 0e8a16 --description "Accepted; will be done, not started" --force
```

### The roadmap: saved searches

A link to the issues with one status label is the roadmap (replace the label in the address):

```
https://github.com/OpenDNA-project/openDNA/issues?q=is%3Aissue+label%3A%22status%3A+planned%22
```

Useful ones: `status: planned`, `status: in progress`, `status: done (next release)`, `status: released`.

## Triage workflow

1. **New.** Look at new issues about once a week (`gh issue list --label "status: new"`). Duplicate: comment with the link to
   the other issue, close as duplicate. Not enough information: ask, set *under review*.
2. **Under review.** Reproduce a bug, or decide whether a feature fits. Move on within a few days or say what you wait for.
3. **Planned** (accepted) or **declined** (close the issue with a short reason; people accept a "no" that is explained).
4. **In progress.** Set it when you start; mention the issue in the commits and pull request with `Refs #12`.
5. **Done (next release).** Merged to `main`. Keep the issue **open**: it closes at the release, so the list of open issues always
   shows what is coming.
6. **Released.** At release time (below) change the label to `status: released` and close the issue with a comment naming
   the version.

Change a status from the command line (the issue number is the example):

```
gh issue edit 12 --remove-label "status: new" --add-label "status: planned"
```

Do not write `Fixes #12` in a commit message if you want step 5 to stay open: GitHub would close the issue when it reaches
`main`.

### Release-notes habit

- While working, when you finish an issue, set `status: done (next release)`.
- At release: list them, put them in the release notes (grouped as *New*, *Fixed*), then flip the labels and close:

```
gh issue list --label "status: done (next release)" --state open --json number,title
gh issue list --label "status: done (next release)" --state open --json number -q ".[].number"
# for each number N:
gh issue edit N --remove-label "status: done (next release)" --add-label "status: released"
gh issue close N --comment "Released in v0.2.0."
gh release create v0.2.0 --notes-file RELEASE_NOTES.md      # or --generate-notes, then edit
```

- Thank the testers by their GitHub names in the release notes when they reported something that was fixed.

## Discussions (optional, not enabled)

GitHub Discussions gives a forum for questions and ideas that are not yet issues. It is off. To enable it: repository >
Settings > General > Features > tick *Discussions*. Then uncomment the *Ask the community* entry in
`.github/ISSUE_TEMPLATE/config.yml`. For a small team of testers, issues alone are enough; consider it when there are more
questions than you can answer in issues.

## What to put in the first test announcement

Short, in the message that goes to the first testers:

- **What it is** in one line (an offline DNA editor for Windows), the version, and where to download it (the release link).
- **What you want tested:** open your usual files (GenBank, FASTA, .dna), edit and save, the tools you use in your work. Say
  that it is a first version and rough edges are expected.
- **Back up** the files they open before editing them (the program autosaves and offers recovery, but a first version deserves
  care).
- **How to report:** *Help > Report a Problem...* for anything wrong or any error window; *Help > Request a Feature...* for
  ideas. One problem per report. Say that no reply is needed to be helpful and that every report is read.
- **Privacy:** reports contain no sequences or file names; please do not attach private files; a made-up example that shows the
  problem is perfect.
- **What happens next:** the labels on each issue show its status (the list of issues labelled `status: planned` is the roadmap);
  releases list what was fixed and who reported it.
- **Access:** while the repository is private, they need a GitHub account that you added as a collaborator (give their GitHub
  names first), or you can ask them to send *Copy Diagnostics* text through your usual channel and file it yourself.
- **Time frame:** for example "tell me by the end of the month; I answer within a week".
