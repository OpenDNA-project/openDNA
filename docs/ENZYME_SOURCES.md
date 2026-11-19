# Enzyme data sources

## Update 2026-11-19: about 200 more enzymes, from Biopython's data

The list grew from 47 to about 250 enzymes (`app/enzyme_table.py`). Requested sources: REBASE, ApE, UGENE. REBASE's own files stay unused (see below,
"All rights reserved"); ApE's and UGENE's lists are their projects' files, not read. Instead `scripts/build_enzymes.py` takes only the facts (site and
cut offsets of each commercially sold enzyme) from `Bio.Restriction.Restriction_Dictionary` in the installed Biopython, checks that it reproduces the 47
hand-written enzymes exactly, drops enzymes that need methylated DNA or have fewer than 4 specific bases (their cuts are not what Biopython itself
reports), and lists isoschizomers once. The residual risk is the one in the text below: the facts are free, but the data file is derived from REBASE, so
the owner should decide whether that is acceptable. To list every isoschizomer, drop the `by_kind` de-duplication in the script.

Isoschizomers are kept as `aliases` of one listed enzyme (views show one name; search, the CLI, the Digest filter and the inventory know all of them). Enzymes that cut only methylated DNA (DpnI, GlaI, BisI, BlsI) carry a `methyl` flag and the filters hide them unless "Needs methylation" is set. The CNNR-type enzymes (MspJI, FspEI, LpnPI, AbaSI) are still left out: their sites have fewer than 4 specific bases.

To regenerate: run the script with the repo's `.venv` python (Biopython is a dev-only dependency; the app does not import it).

## Earlier: NEB list evaluated, not imported (2026-11-02)

Request: replace the hand-written 47-enzyme table (`app/enzymes.py`) with the enzymes of NEB's "Alphabetized List of Recognition
Specificities" (https://www.neb.com/en-ca/tools-and-resources/selection-charts/alphabetized-list-of-recognition-specificities).

## Finding: not importable on what we can verify

- The page and NEB's terms of use could not be read by program: neb.com answers automated requests (curl, WebFetch) with a Cloudflare
  "enable JavaScript and cookies" challenge (HTTP 403). We do not bypass bot checks, so nothing was copied and the page's exact
  content and any page-specific notice were not seen.
- The terms of use found through search (https://www.neb.com/en-us/terms-of-use and sibling locale pages), as reported by the search tool, say:
  content may be used or copied only "for non-commercial individual reference purposes" with notices kept; all site information is NEB
  copyright unless identified otherwise; no licence is implied. Every neb.com page footer reads "All Rights Reserved".
- That is not "free to use, no strings attached". A non-commercial, individual-reference permission does not cover redistributing the
  data inside a program. This repository is AGPL-3.0 (`LICENSE`), which lets anyone redistribute and use it commercially, so NEB's
  terms conflict with shipping their compiled list in it.
- REBASE, the original source, marks its files "All rights reserved" with no licence (see `docs/TODO.md`), and is not used here.

## What is free regardless

- Facts (an enzyme's recognition sequence and cut position) are not copyrightable. `app/enzymes.py` already holds 47 such facts written by
  hand, with no REBASE or NEB file behind them (`THIRD_PARTY_NOTICES.md` says so).
- NEB's compilation (selection and arrangement of the list), descriptions, and wording are theirs; we must not copy those.

## Options to continue

1. Ask NEB (legal@neb.com) for written permission to bundle the list, with attribution. Then import, adding the NEB URL, access date and
   licence to Help > Citations and `THIRD_PARTY_NOTICES.md`.
2. Build a larger table from primary literature or from a source with an open licence (for example Biopython's `Bio.Restriction`, which
   is BSD-licensed code and ships its data derived from REBASE, so check its notice before copying; or enzymes' own publications).
3. Hand-add more enzymes from textbook facts, in our own wording, one at a time.

No code, citation, or notice was changed. `docs/TODO.md` item "NEB enzyme catalogue" still describes the blocked state.
