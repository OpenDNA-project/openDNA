# kld: primers and edited plasmids for a table of mutations (kinase-ligase-DpnI)

**Use when:** the user has one circular plasmid (GenBank) and a list of mutations (substitutions, deletions, insertions, duplications) to make by inverse PCR + KLD. Result: one edited `.gb` per mutation, and one primer pair (forward and reverse) per mutation to order. The results are the same as the Assembly panel's KLD method.

## 1. Look at the plasmid first

```
inspect PLASMID.gb
```

Check that the features the user names exist, which are CDS, their strand, and (for protein mutations such as `D120A`) that residue 120 of that protein is really the amino acid the user says (`protein.sequence`, 1-based). Use `feature` or `region` when you need the bases.

## 2. Build the input: a CSV with exactly two columns

```
id,mutation
m1,lacZ:p.Asp120Ala
m2,lacZ:c.360G>C
m3,o.1234_1240del
```

- `id`: names the output file `<id>.gb`. Use the user's own name for the mutation if it has one, else a convention that suits them (`m1`, `m2`... or `<gene>_<change>`). Ids that would collide (also by letter case) or hold unsafe characters are renamed by openDNA with a warning, and no file is ever lost; but a clear unique id is better.
- `mutation`: HGVS, as in the table of `skills show index` (feature name before the colon; `o.` positions along the plasmid need no feature). One mutation per row. Anything else in the user's table (notes, positions in other formats) you translate yourself into one of the accepted forms when the meaning is clear, and ask when it is not. Other columns are an error.
- A mutation always carries what it expects to find (the reference base, or the amino acid before the number). Keep it: it is the guard against a wrong coordinate.

## 3. Check

```
kld check PLASMID.gb MUTATIONS.csv
```

Nothing is written. The answer is JSON: `settings`, `summary` (rows, ok, errors, with_warnings) and one record per row (`status` ok or error, `error` text, `hgvs_dna`, `hgvs_protein`, `effect`, `warnings`, `primers`). Every row is checked, so fix all the errors you see at once, then check again until `errors` is 0. Exit code 1 means at least one row failed.

## 4. Show the settings, then wait for the user's OK

Show the user the `settings` block: Tm and the other KLD values (default or your override), the PCR conditions (`saved in openDNA` or default), the genetic code table, the codon usage table. Tell them which matter (Tm and the codon usage table change the primers and the codon chosen). Options if they want other values: `--tm` (annealing Tm, default 60), `--group` (changes closer than this share a primer pair, default 30), `--warn` (warn for insertions above this many bp, default 120). Every setting shows the range the windows accept (`min`, `max`): tell the user. A value beyond it is used as given, with a warning that the result is their responsibility. **Do not run before the user agrees.**

## 5. Run

```
kld run PLASMID.gb MUTATIONS.csv --out NEW_FOLDER
```

`--out` must be a folder that does not exist or is empty: nothing is ever overwritten. It writes `<id>.gb` (the edited plasmid with both primers drawn on it), `primers.csv` (id, mutation, hgvs_dna, hgvs_protein, forward and reverse primer with Tm, warnings), `report.jsonl` (one record per row) and `settings.json`. Rows that failed get no file and are listed in the answer (and on stderr); exit code 1; the other rows are done. `summary.files_written` counts the files, and the run checks that each is there.

## What to tell the user

- Where the folder is, how many rows succeeded, and each failed row with its error in plain words.
- **Every warning**, per mutation. Warnings never stop a row: openDNA picks a default and says so. The ones you will see: several codons for the new amino acid (the one used and the alternatives), a change in a repeat that could be placed elsewhere (same plasmid, slightly different primers), a frameshift, a stop codon appearing or lost, an insertion above `--warn`, a mutation that changes nothing ("no change": the plasmid already has those bases, no primers, the file is the plasmid as it is), a file renamed because its id collided or was unsafe, a primer with a low GC or a long primer.
- The `hgvs_protein` and `effect` of each mutation: they say what the user's edit does to the protein.

## Errors you may see

| Message | Fix |
|---|---|
| `the columns must be exactly id,mutation` | remove or rename the extra columns (stops the run) |
| `line N ... has X cells but the header has 2` / `... is empty` | that row only: fix the line (a comma inside a cell needs double quotes) and re-run it |
| `Unexpected ...` | a bug in openDNA: give the user the row and the `detail` |
| `Cannot read the mutation` | the message lists the accepted forms: rewrite it |
| `residue N of X is Y (CODON), not Z` | the amino acid the user gave is not what the plasmid has there: ask the user (a different plasmid or numbering?) |
| `the plasmid has X at that place, not Y` | the expected base is wrong: same question |
| `No feature called ...` / `N features are called ...` | take the name from `inspect`; use `#INDEX` when several |
| `already exists and is not empty` | choose a new `--out` |
| `The sequence is too short for KLD` / `too many` | the change is outside what KLD can do: tell the user |
| `KLD needs a circular sequence` | the file is linear: tell the user |

## Example

This is the sample in the repository (`demo/pBAD30_eilR.gb`); the commands run from the repository root.

```
kld check demo/pBAD30_eilR.gb app/agent/examples/kld_mutations.csv
```
