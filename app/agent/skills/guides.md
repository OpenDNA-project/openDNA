# guides: Cas9 and Cas12a guide RNAs

**Use when:** the user wants guides (spacers) to cut a plasmid or a stretch of it, ranked, with off-target counts. Same as the Spacer Finder window. Read-only.

## Command

```
guides PLASMID.gb [--nuclease SpCas9|SaCas9|Cas12a] [--region START..END] [--max-mm N] [--top N] [--no-genome]
```

- `--region` (1-based, inclusive; START > END on a circular plasmid is a region through the origin): look only there (for example inside a gene from `inspect`). Without it the whole sequence is searched, including sites across the origin of a circular plasmid.
- `--top N`: report and check the best N guides (default 300, the window's list). Checking every guide against a genome takes longer for bigger N.
- `--max-mm`: other sites with up to this many mismatches are counted (default 3; the window accepts 1 to 4).
- The **background genome** is the one in openDNA's Primers > Options (shown in `settings.background_genome`; the source checkout falls back to the demo E. coli genome). `--no-genome` skips it. If none is set, only this sequence is checked: tell the user.

JSON: `guides` best first, each with `start`, `end` (1-based), `strand`, `guide` (the 20 nt written 5' to 3'), `pam`, `gc`, `score` (0-100 checklist, not a measured activity), `cut` (where the nuclease cuts), `other_sites_here` and `other_sites_genome` (lists of counts of the *other* places the guide matches, by number of mismatches: the item at index 0 is exact copies, index 1 one mismatch, up to `--max-mm`; `[0, 0, 2, 5]` = no copy, none with 1 mismatch, 2 with 2, 5 with 3), `places_next_to_n` (how many places could be one more site but have an N in the DNA where the PAM is: not counted, and the mark is `warn` when it is not 0), `mark` (`ok`, `warn`, `bad`) and `verdict`.

## What to tell the user

The best guides with their mark and the reason when it is not ok (`verdict`). That the score is a rule-of-thumb checklist, and the genome used. The guide's `start`..`end` includes the PAM.

## Errors you may see

| Message | Fix |
|---|---|
| `No nuclease called ...` | SpCas9, SaCas9 or Cas12a |
| `positions must satisfy` | `START..END` inside the plasmid |
| `select a region first` | the sequence is over 2 Mb: pass `--region` |

## Example

```
guides demo/pBAD30_eilR.gb --region 1..300 --top 3 --no-genome
```
