# search: find a sequence, a motif or a name

**Use when:** the user asks where something is in a plasmid: a DNA sequence (IUPAC codes allowed), a protein motif, or the name of a feature or enzyme. Same engine as the search bar. Read-only.

## Command

```
search PLASMID.gb QUERY [--fuzzy N] [--single-strand] [--protein] [--details] [--enzymes]
```

- DNA (`GAATTC`, `TTGACA`, `RGATCY`): found as written and as its reverse complement (`--single-strand`: only as written). `--fuzzy N`: up to N mismatches. A query of IUPAC letters shorter than 5 is looked up as a name instead.
- `--protein`: a protein motif (one-letter code, `X` = any) is searched in the translations of both strands.
- A name finds features whose label, gene or product contains it (`--details`: any qualifier).
- `--enzymes`: an enzyme name finds its sites.

JSON: `hits` in order of position: `kind` (`dna`, `protein`, `feature`, `enzyme`), `label`, `start`..`end` (1-based), `strand`, `mismatches`, and `feature_index` for a feature (usable as `#INDEX` in other commands). `note` says what was not searched.

## Example

```
search demo/pBAD30_eilR.gb GAATTC
search demo/pBAD30_eilR.gb araC
search demo/pBAD30_eilR.gb ATGGCT --fuzzy 1
```
