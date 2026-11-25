# signal-regions: where a signal is high

**Use when:** the user has a signal track (bigWig, bedGraph, wig: RNA-seq coverage, insertion density, ChIP...) and wants the stretches where it is at least some value: expressed regions, peaks, regions without insertions (with a low threshold on the inverse question, ask the user). The same as the window's right-click on a signal > Regions Above Threshold.

## 1. Find them

```
signal-regions SIGNAL --threshold VALUE [--min-length BP] [--gap BP] [--region START..END] [--seq NAME] [--out NEW.bed]
```

A stretch is a run of bases whose value is at least `--threshold`. Stretches `--gap` bp apart or closer are joined into one; stretches shorter than `--min-length` are left out (after joining). `--region` limits the search (1-based, inclusive); `--seq` names the sequence when the file has several.

The answer lists each stretch with `start`, `end` (1-based), `length`, its `mean` and `max` (read from the values themselves, over the whole stretch: joined bases included, bases with no value not counted). `--out` also writes them as BED (name `<file>_<n>`, the mean as the score), to load as a track or to give to `recap --regions`.

## What to tell the user

- How many stretches, their total length, and the highest ones; the threshold and joining used.
- To see them on the sequence: load the BED file as a track, or use Regions Above Threshold in the window (it can make them highlights or features).

## Example

```
signal-regions demo/hub/pBAD30_eilR/plus.bw --threshold 20 --min-length 100
```

## Errors you may see

| Message | Fix |
|---|---|
| `name the sequence with --seq` | the file has several sequences: pick the user's |
| `give START..END with ...` | the region is beyond the sequence or written another way |
| `cannot read ...` | not a bigWig, bedGraph or wig file |
