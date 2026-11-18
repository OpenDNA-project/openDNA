# gibson: Gibson / HiFi assembly of fragments

**Use when:** the user wants to join DNA fragments (a backbone, inserts, parts cut from several files) by Gibson or NEBuilder HiFi assembly: primers that add the overlaps, the checks of every junction, and the assembled plasmid. Same as the Assembly panel's Gibson method. (Golden Gate, Sidewinder, restriction/ligation and the panel's other options are not available on the command line yet.)

## 1. Look at every source first

`inspect FILE.gb` for each file: pick the stretch of each fragment from its features (a gene, a backbone between two sites).

## 2. Build the input: a CSV, one row per fragment, in the order they are joined

```
id,source,region,orientation
backbone,pBAD30_eilR.gb,2000..4000,forward
insert,demo_insert.gb,,forward
```

- `id`: the fragment's name (it names its primers, `backbone-F`, `backbone-R`); different for every row.
- `source`: a GenBank or FASTA file (a relative path is read from the table's folder).
- `region`: `START..END` (1-based, inclusive, START smaller than END) or a feature name or `#INDEX`; empty = the whole file.
- `orientation`: `forward` (default) or `reverse`: the fragment is used as the reverse complement of that stretch.
- Columns are exactly these (`region`, `orientation` optional). **One bad row stops the command** (unlike the batch commands: an assembly needs every fragment); the message names the rows.
- The last fragment joins the first (circular product). `--linear` for a linear one.

## 3. Options

`--min-overlap` (default 20 bp) and `--max-overlap` (40): the overlap each junction gets, the shortest that reaches `--overlap-tm` (52 deg C); `--tm` (60) is the annealing Tm of the primers. Each has the range the window accepts in `settings`; beyond it the value is used as given, with a warning that the result is the user's responsibility.

## 4. Check, show the settings, run

```
gibson check FRAGMENTS.csv
gibson run FRAGMENTS.csv --out NEW_FOLDER [--name my_plasmid]
```

`check` writes nothing. Show `settings` and wait for the user's OK before `run`. JSON: `verdict` (`ok`, `warn`, `bad`: the worst check), `junctions` (each with a level; an end that does not fit is `bad`), `checks` (overlap Tm, GC, whether the overlap is unique in the product, hairpins, product length, and so on), `primers` (name, sequence with its 5' overlap tail, `annealing_length`, `tail_nt` = the bases the primer adds at its 5' end to make the overlap with its neighbour, and Tm), `product` (length, circular, file). `run` writes the assembled plasmid (`assembly.gb`, or `--name`), `primers.csv`, `report.json` and `settings.json`; nothing is overwritten. A junction reads `existing overlap` when the fragment end and the next start are the same bases, which includes the overlap that the primers' tails make (`tail_nt` > 0); if the fragment ends already overlapped, `tail_nt` is 0 and nothing is added. As in the window, the circular product starts a few bases before the first fragment (its overlap with the last): the same circle, turned.

## What to tell the user

The verdict, then every check that is `warn` or `bad` in plain words (an overlap that occurs twice in the product would join the wrong pieces; a low Tm or a hairpin in an overlap lowers efficiency), and the primers to order (`order` turns them into a supplier file). The tool designs and checks the assembly; it does not model the reaction itself.

## Errors you may see

| Message | Fix |
|---|---|
| `the columns must be exactly id,source` | remove or rename the extra columns |
| `the fragments cannot be assembled with a row missing` | fix the rows named, then re-run |
| `The ids must be different` | give each fragment its own id |
| `positions must satisfy` / `No feature called` | check the region against `inspect` |
| `--min-overlap is above --max-overlap` | swap them |

## Example

```
gibson check app/agent/examples/gibson_fragments.csv
```
