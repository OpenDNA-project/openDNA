# recombineering: lambda Red primers and recombinant genomes for cassettes that replace regions

**Use when:** the user wants to replace a region of a genome (or a big plasmid) with a cassette (a resistance gene, a tag...) by lambda Red recombineering, or to delete a region. Result per row: the two primers that amplify the cassette (the genome's homologies are their 5' tails) and the genome as it will be after recombination. The results are the same as the Recombineering tool's design.

## 1. Look at the genome and the cassette first

```
inspect GENOME.gb
inspect CASSETTE.gb
```

The genome must be a GenBank file (circular for a bacterial chromosome). The cassette is a GenBank file whose **whole sequence** is what goes in: with the promoter, gene and terminator it needs, and without extra bases. Check the region the user names against the genome's features.

## 2. Build the input: a CSV

```
id,region,cassette,orientation
ko_eilR,eilR,cassettes/kan.gb,forward
del_araC,araC,none,
```

- `id`: names the output genome `<id>.gb`; choose a convention that suits the user. Colliding or unsafe ids are renamed with a warning (see the index), never lost.
- `region`: what is replaced: a feature name (or `#INDEX`), or `START..END`, 1-based and inclusive, along the genome. On a circular genome START > END is a region through the origin (`890..10`); a feature that crosses it is used as it is.
- `cassette`: the cassette's GenBank file (relative to the CSV's folder, or an absolute path), or `none` to delete the region (no primers: there is nothing to amplify).
- `orientation`: optional column, `forward` (default) or `reverse` (the cassette goes in reverse-complemented). Say to the user which one you used; for a gene knock-out it usually matches the direction of the gene it replaces.
- Insertions between two bases are not offered here: replace a small region instead.
- Anything else in the user's table you translate yourself when the meaning is clear, and ask when it is not. Other columns are an error.

## 2b. Options (all optional)

`--homology` (bp of homology on each side of the cassette, default 50), `--tm` (annealing Tm of the primers, default 60), `--primers-only` (with `run`: write only the primers, not the recombinant genomes. Each recombinant is a full copy of the genome: for a bacterial genome that is several megabytes per row, so 200 rows are about a gigabyte of files and a lot of time. Use it when the user only wants the primers to order; leave it off when they want the edited genomes). `check` never builds genomes.

## 3. Check

```
recombineering check GENOME.gb TABLE.csv
```

Nothing is written. JSON: `settings`, `summary` and one record per row: `status`, `error`, the region, the cassette length, the `left_homology` and `right_homology` sequences, both `primers` (name, full sequence with the homology tail, annealing length and Tm) and `warnings`. Exit code 1: at least one row failed.

## 4. Show the settings, then wait for the user's OK

Show the `settings` block (homology length, annealing Tm, whether genomes are written, the PCR conditions) and **the warnings of every row**: they say which features the region removes and which it cuts in two. Check with the user that those are meant. Every setting shows the range the windows accept (`min`, `max`): tell the user; a value beyond it is used as given with a warning that the result is their responsibility. **Do not run before the user agrees.**

## 5. Run

```
recombineering run GENOME.gb TABLE.csv --out NEW_FOLDER
```

`--out` must not exist or be empty; nothing is ever overwritten. It writes `<id>.gb` (the recombinant genome, with the cassette's features and the two primers drawn on it; not with `--primers-only`), `primers.csv`, `report.jsonl` and `settings.json`.

## What to tell the user

- Per row: the region replaced, the cassette used and its orientation, the two primers (to order) and the warnings.
- Rows that failed, with the error in plain words.
- That the primers are for amplifying the cassette: their 5' tails are the genome's homologies. Homology on a bacterial genome is checked against the genome itself only; the tool does not know the user's strain.

## Errors you may see

| Message | Fix |
|---|---|
| `the columns must be exactly id,region,cassette` | remove or rename the extra columns |
| `No feature called ...` / `N features are called ...` | take the name from `inspect`; use `#INDEX` |
| `File not found` | the cassette path is relative to the CSV's folder |
| `orientation must be forward or reverse` | fix the value |
| `no room for the homology next to the region` | the region is too near the end of a linear file, or the homology is longer than the genome |
| `the region to replace does not fit the genome` | positions outside the genome |
| `already exists and is not empty` | choose a new `--out` |

## Example

This is the sample in the repository (`demo/pBAD30_eilR.gb` with `demo/demo_insert.gb` as the cassette); the commands run from the repository root.

```
recombineering check demo/pBAD30_eilR.gb app/agent/examples/recombineering.csv
```
