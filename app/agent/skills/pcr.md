# pcr: primer pairs to amplify a list of regions

**Use when:** the user wants PCR primers (forward and reverse) for one or several regions of a plasmid: a feature (a gene, a part) or a stretch of bases. Result: the best pairs per region, each with the Primers panel's checks. The results are the same as the PCR panel's "Design".

## 1. Look at the plasmid first

```
inspect PLASMID.gb
```

Take the feature names from it, and check each start and end the user gave against the features.

## 2. Build the input: a CSV with exactly two columns

```
id,target
araC,araC
ori_region,2839..3384
```

- `id`: a name for the row (no file is named after it); choose one that suits the user.
- `target`: a feature name (or `#INDEX` when several share a name), or `START..END`, 1-based and inclusive, along the plasmid (START smaller than END). A region across the origin is not accepted: ask for it in two parts, or ask the user.
- Primers sit exactly on the two ends of the target (`--window 0`, the default), so the product is the target. A primer needs a product at least as long as itself, so a target of about 25 bases or fewer gives no pair: widen it or say so.
- Anything else in the user's table you translate yourself when the meaning is clear (for example a start and an end in two columns become `START..END`), and ask when it is not. Other columns are an error.

## 2b. Options (all optional)

`--tm` (target Tm in deg C, default 60), `--min-len` (18) and `--max-len` (25) for the primer length, `--window` (default 0: how many bases before the start and after the end a primer may begin; above 0 the product may be longer than the target and only primers within 5 deg C of the target Tm are used).

## 3. Check

```
pcr check PLASMID.gb TARGETS.csv
```

Nothing is written. JSON: `settings`, `summary` and one record per row: `status` (ok or error), `error`, and up to five `pairs` (rank 1 is the best): both primers with `sequence`, `tm`, `gc`, position, the `checks` of each (Tm, GC, 3' clamp, hairpin, self-dimer, where it binds, non-specific priming in this sequence and in the background genome), the product size, the penalty (lower is better), `notes` and a `mark` (`ok`, `warn`, `bad`). Exit code 1: at least one row failed (for example `No pair found`).

## 4. Show the settings, then wait for the user's OK

Show the `settings` block: the Tm and length options, the PCR conditions (salt, dNTP, primer concentrations: `saved in openDNA` or default), and the **background genome** that primers are checked against (the `background_genome` of the settings; the user sets it in openDNA under Primers > Options). If none is set the genome check says `no genome set` and primers are only checked against the plasmid: tell the user. A source checkout falls back to the demo E. coli genome, so a primer for an E. coli gene is then marked `bad` ("the whole primer binds" in the genome) for the right reason. Every setting shows the range the windows accept (`min`, `max`): tell the user; a value beyond it is used as given with a warning that the result is their responsibility. **Do not run before the user agrees.**

## 5. Run

```
pcr run PLASMID.gb TARGETS.csv --out NEW_FOLDER
```

`--out` must not exist or be empty. It writes `primers.csv` (one line per id and rank: both primers, Tm, product size, penalty, mark, notes), `report.jsonl` (the full records with the checks) and `settings.json`. Nothing is ever overwritten.

## What to tell the user

- Per target: the best pair (rank 1) with its Tm, the product size and the `mark`. If the mark is `warn` or `bad`, say which check says so (read the `checks` of each primer) and offer the next pair.
- Targets with no pair, and why.
- That the checks use the background genome named in the settings.

## Errors you may see

| Message | Fix |
|---|---|
| `the columns must be exactly id,target` | remove or rename the extra columns |
| `No feature called ...` / `N features are called ...` | take the name from `inspect`; use `#INDEX` |
| `positions must satisfy` | `START..END` inside the plasmid, START smaller than END |
| `No pair found` | the target is too short, or widen `--window`, the length range, or the Tm |
| `already exists and is not empty` | choose a new `--out` |

## Example

This is the sample in the repository (`demo/pBAD30_eilR.gb`); the commands run from the repository root.

```
pcr check demo/pBAD30_eilR.gb app/agent/examples/pcr_targets.csv
```
