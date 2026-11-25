# track-features: a track's items as features of the sequence

**Use when:** the user wants the items of an annotation track (BED, bigBed, GFF, GTF) to become features of their GenBank file: genes of another annotation, transcription start sites, regions found by an experiment. The same as the window's right-click on a track > Copy to Features.

## 1. Copy

```
track-features SEQUENCE.gb TRACK --out NEW.gb [START..END] [--type TYPE] [--name-in QUALIFIER]
```

Every item of `TRACK` on this sequence (or only those overlapping `START..END`, 1-based, or a feature's name) becomes a feature of the copy written to `--out` (a new file: the original is not changed). Each feature gets:

- the item's name, in `/label` (or the qualifier given by `--name-in`: `gene`, `locus_tag`, `product`, `note`);
- a note `from track <file name>`, the item's own note, and each extra column of the BED file as a note of its own (`height=81`, `TU_length=1380`...);
- the type `--type` (default: `CDS` for a track whose name ends in `_CDS`, else `misc_feature`).

An item that is already a feature of the sequence (same type, same place, same strand) is skipped and counted in `skipped`: running it twice adds nothing the second time.

## What to tell the user

- How many features were `added` and `skipped`, their type, and the new file.
- To see them: open the new file; the features table lists them with their notes.

## Example

```
track-features demo/pBAD30_eilR.gb demo/hub/pBAD30_eilR/sites.bb --out example_sites.gb
```

## Errors you may see

| Message | Fix |
|---|---|
| `... is not an annotation track of this sequence` | the file is a signal or BAM file, or its sequence names do not match this sequence (see the warning) |
| `... already exists` | choose a new `--out` |
| `--name-in is one of ...` | use one of the qualifiers listed |
