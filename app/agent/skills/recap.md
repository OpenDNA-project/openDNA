# recap: what every track holds at a place

**Use when:** the user has tracks on a sequence (an Open DNA Project with its tracks, or a sequence and track files) and wants to know what they say at a place or at several places: which annotation items are there, the level of each signal (RNA-seq, insertion density...), the read depth. The same numbers as the window's right-click > Recap Tracks Here.

## 1. Recap a place

```
recap SOURCE START..END [--track FILE ...]
```

`SOURCE` is an Open DNA Project (`.odp`: its tracks are loaded, in its order and groups) or a sequence file; `--track` adds a track file (GFF, GTF, BED, bigBed, VCF, bigWig, bedGraph, wig, BAM), as many as needed. The place is `START..END` (1-based, inclusive) or a feature's name.

Instead of a place: `--regions places.bed` recaps each region of a BED file; with neither, the project's highlighted regions.

The answer is JSON, one record per region: `gc` (%), `features` (the sequence's own, with the share of the region they cover), and `tracks` in the view's order, each with its `group` and:

- an annotation track: `items` (name, type, start, end, strand, `covers`: share of the region, `columns`: the BED file's extra columns);
- a signal: `bases_with_value`, `mean`, `min`, `max`, `sum`, read from the values themselves (exact); with `--per-base` and a region of at most 1,000 bp, `values` base by base;
- a BAM file: `mean_depth`, `max_depth`.

## 2. As a table

```
recap SOURCE START..END --tsv
```

One tab-separated row per region and track, and after a group's tracks a row `<group> (mean of N)` with the mean of their numbers: easy to paste into a spreadsheet. The warnings go to stderr.

## What to tell the user

- The place recapped (positions as given), and the tracks with something notable there: items found, the highest and lowest signals. Compare signals by `mean`; `sum` grows with the region's length.
- `bases_with_value` below the region's length: the signal has no value on part of it (not zero: no data).
- Any warning: a track not loaded (its sequence names do not match), a project file not found.

## Example

```
recap demo/pBAD30_eilR.gb 1001..2600 --track demo/hub/pBAD30_eilR/plus.bw --track demo/hub/pBAD30_eilR/sites.bb
recap demo/pBAD30_eilR.gb 1001..1020 --track demo/hub/pBAD30_eilR/plus.bw --per-base --tsv
```

## Errors you may see

| Message | Fix |
|---|---|
| `no region: ...` | give `START..END`, a feature's name, or `--regions FILE.bed` |
| `positions must satisfy 1 <= START < END <= ...` | the place is beyond the sequence: check its length with `inspect` |
| `Track file not found` | check the path of each `--track` |
