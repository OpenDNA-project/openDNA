# codon: recode a gene for another organism's codon usage

**Use when:** the user wants one or several CDS features rewritten with the codons an organism prefers, keeping the protein exactly, optionally with no sites of given enzymes, with GC and run-length limits. Same as the Recode window. A recoded plasmid is written per row.

## 1. Look first

`inspect PLASMID.gb` gives the CDS names; a CDS that is a join, crosses the origin or has `/transl_except` is refused with a message.

## 2. Build the input: a CSV with exactly two columns

```
id,feature
araC,araC
```

`id` names the output file (`araC.gb`); `feature` is a CDS name or `#INDEX`.

## 3. Options

`--usage`: the codon usage table: the name of one in openDNA, or **a CSV file** (`codon,aa,per_thousand`; the window's Import reads the same). Without it the table chosen in openDNA is used; if there is none, ask the user for one. `--method` `"Frequency matched"` (default: each codon drawn in proportion to the table, no long run of the best one) or `"Highest frequency"`. `--avoid EcoRI,BsaI`: sites to keep out. `--gc-low`, `--gc-high` (default 35 to 65 percent in every `--window` of 50 bases; `--no-gc` switches it off). `--max-run` (default 5 identical bases), `--no-rbs` (allow an internal AGGAGG). `--seed` (default 1): a frequency-matched recoding is a random draw, so the **same seed gives the same sequence**; try another seed for another draw.

## 4. Check, show the settings, run

```
codon check PLASMID.gb TARGETS.csv --usage TABLE.csv
codon run PLASMID.gb TARGETS.csv --usage TABLE.csv --out NEW_FOLDER
```

Show `settings` (the table, method and limits, each with the range the window accepts) and wait for the user's OK before `run`. Per row: `codons`, `codons_changed`, `gc_before`/`gc_after`, the recoded `sequence` (read on the coding strand) and `warnings`: a rule that could not be met (`Could not be met: GC content at bases ...`). Relay them: the protein is identical but the rule is not met. `run` writes one plasmid per row (the CDS replaced in place, same length, features untouched) plus `primers.csv` (a summary here), `report.jsonl` and `settings.json`.

## What to tell the user

How many codons changed, the GC before and after, any rule not met, and that the protein is unchanged (the tool refuses to write a recoding that changes it). Recoding changes the DNA only: the effect on expression is not predicted.

## Errors you may see

| Message | Fix |
|---|---|
| `No codon usage table is set` | give `--usage` |
| `is a gene, not a CDS` | name a CDS feature |
| `made of several parts` | not supported for a join or an origin-crossing CDS |
| `--gc-low is above --gc-high` | swap them |

## Example

```
codon check demo/pBAD30_eilR.gb app/agent/examples/codon_targets.csv --usage app/agent/examples/codon_usage_ecoli.csv
```
