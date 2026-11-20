# checkup: should each gene of the plasmid be expressed?

**Use when:** the user wants to know whether a construct will work as designed: for each CDS, is there a promoter upstream, a ribosome binding site, a terminator, a start codon, a stop with nothing in between? Same as the Plasmid Checkup window. Read-only. It checks the **design**, not expression levels: real terminators leak and predicted elements are only motifs.

## Command

```
checkup PLASMID.gb [--min-promoter N] [--min-stem N] [--rbs-reach N] [--no-predict] [--no-library]
```

Elements come, in order of trust, from the file's own features (promoter, terminator, RBS, origin), from known parts found in it (the Features Database, as `annotate` does), and from predicted motifs (bacterial sigma-70 promoters, intrinsic terminators, Shine-Dalgarno sites). `--no-predict` / `--no-library` leave a source out. The thresholds are the user's Settings > Plasmid Checkup values: `settings` shows each with where it came from and the range the window accepts; a value given here is used for this run only and never saved. Show the settings and wait for the user's OK before changing one.

JSON: `cds`, one per CDS: `verdict` (`Expressed`, `Likely`, `Unlikely`, `Not expressed`) in its `best_host`, `reasons`, and per host (`Bacteria`; also `Eukaryotic cells` when the plasmid has a eukaryotic promoter) the `promoter`, `terminator` and `rbs` used. `elements` (every promoter, terminator, RBS and origin with its `origin`: annotated, library or predicted), `origins` (each: `label`, `start`, `end`, `strand`, `name`, `host`, `copies_per_cell`, `incompatibility_group`, `note`) and `origin_notes` (for example two origins that may be incompatible).

## What to tell the user

Each CDS's verdict with its reasons in plain words (`relies on predicted elements` means no annotated promoter or terminator backs it: a prediction, to be confirmed). Promoters that need an inducer (araBAD, lac, tet, T7...) are named in the reasons: say which inducer. A CDS marked `Not expressed` may only lack an annotation: ask whether the promoter is missing or just not annotated (`annotate` can find known ones).

## Example

```
checkup demo/pBAD30_eilR.gb
```
