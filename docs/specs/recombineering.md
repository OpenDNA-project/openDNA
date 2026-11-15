# Recombineering simulator: spec

Status: agreed in a design interview, not yet implemented.

## Goal

Simulate λ Red recombineering with a dsDNA cassette. Two ways in, one engine:

- **Design**: give the cassette, the region to replace and the homology length. The tool returns the
  primers and the sequence of the modified genome.
- **Check**: give a finished fragment (or a primer pair and a template). The tool finds where each
  homology lands in the genome, what it would replace, and the resulting genome.

## Decisions

1. **Where it lives.** A new method, "Recombineering", in the assembly panel beside Gibson, Golden
   Gate and Sidewinder. It reuses the panel's flow (fragments, design, simulate, checks) and the
   helpers in `app/assembly.py` (`Design`, `auto_primers`, `simulate`, `placed_features`, `wrapped`,
   `anneal_len`, `oligo_checks`).
   In the panel the two boxes are "Genome" and "Cassette". The genome box holds what is kept (from the
   region's end round to its start), so the checks are the Gibson ones on a two-fragment circle; Perform
   assembly edits the genome itself (`apply`), which keeps its origin and features. The cassette's
   sense/antisense is the existing flip of its box.
2. **Scope of v1.** dsDNA cassette only (a PCR product with homology tails, e.g. KanR from pKD4).
   Replacement, deletion (empty cassette) and insertion (empty replaced region) are the same
   calculation. ssDNA oligo recombineering (MAGE style) is out of scope.
3. **Genome and cassette.** The genome is the open document. The cassette is a selection or feature
   of another open document, as for the other fragments. The primers anchor at the two ends of that
   selection.
4. **Replaced region.** A selection or a whole feature on the genome map, with coordinates editable
   by hand. A "keep start/stop" box serves gene deletions.
5. **Homologies and primers.** 50 nt per side by default, each side adjustable (35 to 100). Each
   primer is a 5' homology tail plus a ~20 nt 3' anchor at Tm about 60 C (`anneal_len`). Tm,
   hairpin and dimer checks come from `oligo_checks`.
6. **Output.** A new document, never an edit in place. Features inside the replaced region are
   removed or clipped, cassette features are transferred (`placed_features`), downstream
   coordinates shift.
7. **Check mode input.** Either a full fragment with its tails, or a primer pair plus a template
   (the cassette is then the amplicon between the primers). Homologies on both sides are searched
   in the genome with `opendna_core.align`; no new search code in the core.
8. **Warnings and mismatch tolerance.** The search tolerates a few mismatches (default up to 5 per
   50 nt arm, adjustable) so that non-specific sites are found, not only the intended one. Each
   hit lists its mismatch count and positions. Two warnings: a homology that is not an exact match
   at the intended site, and a homology that maps to more than one place, near-matches included.
   The intended site is the best hit; the others are listed as non-specific. Nothing else is
   judged (no essential-gene call, no distance warnings). In design mode the arms are cut from the
   genome, so they match exactly by construction.
9. **Circular genomes.** The genome is circular by default, so the replaced region may cross the
   sequence origin. The design interface has a toggle for the cassette orientation (sense or
   antisense; default: same orientation as the replaced gene). A dedicated test covers the origin
   crossing.
10. **One edit at a time.** A design is one cassette and one region. To chain edits, use the output
    document as the next genome. No batch mode.

## Dropped on purpose

- Verification PCR primers and expected amplicon sizes.
- FLP/FRT or Cre/loxP excision and the scar genome.
- Multi-edit batch design.

## Tests and report

Per CLAUDE.md: extend the existing test whose role fits (logic only, no widgets), add entries to
`docs/test-report/tests.json`, and rebuild the report. Cases to cover: replace, delete, insert;
origin crossing; reverse orientation; check mode with a wrong arm, a double-hit arm, and a near-match (a few mismatches) elsewhere
in the genome that must be reported as non-specific.
