# openDNA

A desktop sequence editor and viewer for molecular biology: plasmids, genomes, reads and the
analyses run on them.

## Language

### Terminal

**Terminal**:
The tab where command-line programs run against the current document and other files, in one
persistent shell. Everything external openDNA runs goes through it.
_Avoid_: console, shell tab

**Pipeline**:
A ready-to-run chain of command-line steps with a form of options (for example Align reads,
Assemble), listed on the left of the Terminal.
_Avoid_: recipe, tool, workflow, script

**Preset**:
A named set of option values for a pipeline, chosen for a common case so that a user need not
understand each option.
_Avoid_: profile, template, mode

**QC summary**:
The few lines a pipeline run ends with, saying how it went (reads kept, percent mapped, assembly
size and N50) and where the tools' full reports are.
_Avoid_: report, QC report

**Run folder**:
The folder a pipeline run leaves beside its output, holding the inputs, the exact commands and the
log, so the run can be repeated by hand.
_Avoid_: work dir, temp folder

### Recombineering

**Cassette**:
The dsDNA piece that goes into the genome, usually a selection marker. It is a selection or feature
of an open document, amplified with primers whose 5' tails carry the homologies.
_Avoid_: insert, marker

**Homology**:
The stretch at each end of the amplified cassette that is identical to the genome beside the
replaced region: 50 nt by default, one per side. Exact by construction when designing; when
checking, a few mismatches are tolerated so that non-specific sites show up.
_Avoid_: homology arm, flank, overlap (that is Gibson's word)

**Replaced region**:
The stretch of genome between the two homologies, lost when the cassette goes in. Empty for a pure
insertion; with an empty cassette the edit is a deletion.
_Avoid_: target, deleted region

**Genome**:
The open circular document the edit is simulated on. In the assembly panel its box holds what is
kept (from the region's end round to its start); Perform assembly edits the genome itself, so its
origin and features stay where they were.
_Avoid_: chromosome, host

**Design / Check**:
Design builds the cassette's primers and the recombinant genome from a region and a cassette;
check takes a fragment or a primer pair that already exists and reports where each homology lands,
what it would replace and which other places it matches.
