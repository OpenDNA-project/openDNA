# annotate: find known parts in a plasmid

**Use when:** the user wants a plasmid annotated: promoters, terminators, origins, resistance genes, tags and other known parts looked up in openDNA's Features Database (the bundled collections plus the user's own). Same as Detect Common Features.

## Command

```
annotate PLASMID.gb [--kinds promoter,terminator] [--min-dna N] [--min-protein N] [--min-length N] [--add NEW.gb]
```

Both strands are searched, across the origin. A coding part is also searched as protein, to find codon-optimized versions. The thresholds are the user's Settings > Detect Features values (95 percent DNA identity, 90 percent protein, parts of 15 bp or more by default); `settings` shows each with its source and the range the window accepts. A value given here is used for this run only.

JSON: `hits` in order of position: `name`, `kind`, `start`..`end` (1-based), `strand`, `identity` (0..1), `protein_match` (found as protein: the DNA differs), `already_annotated` (the file already has a feature of that kind on exactly those bases, on the same strand for a CDS; a feature that only overlaps does not count, so a part inside or beside another one is listed as new: say so), `source` (where the part comes from), `strength` (a measurement the source reports, with its units, for promoters and terminators).

`--add NEW.gb` writes a copy of the plasmid with the hits that are **not** annotated yet added as features (never overwrites; the original is untouched). Features are added one per part, not grouped into genes.

## What to tell the user

The hits with their identity and source; which ones are new and which were already there. A part found at less than 100 percent identity is a variant of the library part: say so. Parts shorter than 30 bp can match by chance in a long sequence. Ask before using `--add` on anything but a copy.

## Example

```
annotate demo/pBAD30_eilR.gb
```
