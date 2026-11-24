# annotate: find known parts in a plasmid

**Use when:** the user wants a plasmid annotated: promoters, terminators, origins, resistance genes, tags and other known parts looked up in openDNA's Features Database (the bundled collections plus the user's own), and, when it is installed, in pLannotate (the Barrick lab's annotator, a separate tool: McGuffie & Barrick, Nucleic Acids Res 2021). Same as Detect Common Features.

## Command

```
annotate PLASMID.gb [--kinds promoter,terminator] [--min-dna N] [--min-protein N] [--min-length N] [--no-plannotate] [--add NEW.gb]
```

Both strands are searched, across the origin. A coding part is also searched as protein, to find codon-optimized versions. The thresholds are the user's Settings > Detect Features values (95 percent DNA identity, 90 percent protein, parts of 15 bp or more by default); `settings` shows each with its source and the range the window accepts. A value given here is used for this run only.

When pLannotate is installed (`wsltools status`; on Windows it runs inside WSL) it searches the same sequence with its own databases (SnapGene-derived features, FPbase, Swiss-Prot, Rfam), and its features are merged into the list: `found_by` is `openDNA`, `pLannotate` or `openDNA + pLannotate` (both found the same thing in the same place: the name is the Features Database's, `plannotate` gives the other tool's). The top-level `plannotate` says `ran`, `off` (`--no-plannotate` or the setting), `missing` (not installed: a warning says how to install it) or `error`. pLannotate takes seconds, only for sequences up to 50,000 bp, and has no `--kinds` option (its hits are filtered by type afterwards). The thresholds above apply to the Features Database only. When it ran, `plannotate_citation` holds the reference (with its DOI): give it to the user if they will use or publish those results.

JSON: `hits` in order of position: `name`, `kind`, `start`..`end` (1-based), `strand`, `identity` (0..1), `protein_match` (found as protein: the DNA differs), `already_annotated` (the file already has a feature of that kind on exactly those bases, on the same strand for a CDS; a feature that only overlaps does not count, so a part inside or beside another one is listed as new: say so), `annotated_other_way` (the file has that kind on exactly those bases but in the other direction: the file's feature is probably backwards; say so), `found_by`, `source` (where the part comes from; for pLannotate its database after `pLannotate:`), `strength` (a measurement the source reports, with its units, for promoters and terminators).

`--add NEW.gb` writes a copy of the plasmid with the hits that are **not** annotated yet added as features (never overwrites; the original is untouched). Features are added one per part. A CDS, promoter, RBS or primer binding site annotated the other way is turned round in the copy (its name and notes kept; `added.turned_round` counts them), as the window does; other types annotated the other way are left as they are.

## What to tell the user

The hits with their identity, source and who found them (name pLannotate when it did: it is someone else's work); which ones are new and which were already there. A part found at less than 100 percent identity is a variant of the library part: say so. Parts shorter than 30 bp can match by chance in a long sequence. Ask before using `--add` on anything but a copy.

## Example

```
annotate demo/pBAD30_eilR.gb
```
