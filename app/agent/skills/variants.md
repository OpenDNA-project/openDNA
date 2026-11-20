# variants: where sequencing reads differ from the plasmid

**Use when:** the user has sequencing reads of a plasmid, already mapped to it (a sorted, indexed BAM), and asks what differs from the reference: point mutations, small insertions and deletions, breaks that point to a larger rearrangement, mixed colonies. Same scan as the Variants pane. Mapping reads is a separate step (the Terminal tab's "Align reads" pipelines, or any aligner): this command does not map.

## Command

```
variants PLASMID.gb READS.bam [--seqid NAME] [--min-depth N] [--variant-pct N] [--consensus-pct N] [--min-baseq N] [--min-mapq N] [--apply NEW.gb]
```

- Nothing is left out for being near an end. On a circular plasmid, reads that run over the origin (mapped as two pieces, one clipped at each end of the reference) are told from breaks by their clipped bases, which must continue the circle: such reads count in full (a change at base 1 or at the last base has its real depth) and are not breaks. A clip that does not continue the circle is a break, at the ends too. Reads that are hard-clipped there (their bases are not in the file) are reported as a `clip` whose text says it could not be checked: ask the user whether the reads were mapped with soft clips. A deletion through the origin is one `del` (its text says "through the origin"; `--apply` writes it), and an uncovered start and end are one `gap` through the origin. The BAM's reference must be this plasmid. The sequence is found by name (`--seqid` when the BAM has several); a different length or name is a warning: the reads may not line up, so tell the user.
- Thresholds default to the user's Settings > Alignments (`settings` shows each value, its source and the range the window accepts): at least `--min-depth` reads (5), `--variant-pct` of them differing (20 percent) to list a change, `--consensus-pct` (80 percent) to call it a **consensus** change (the plasmid itself differs from the reference), base and mapping qualities as given. Show the settings and wait for the user's OK to change them.
- `--apply NEW.gb`: also writes a copy of the plasmid with every consensus substitution, insertion and deletion written in (never overwrites). Tell the user which changes were written: the ones with `level` consensus and `can_apply` true; a `variant` (minority) level, a `clip` or a `gap` is never written. After an insertion or deletion the coordinates move: tell the user to map the reads again against the copy to confirm.

JSON: `variants` in order of position: `position` (1-based), `kind` (`snp`, `del`, `ins`; `clip`: reads that stop matching there, the edge of an insertion, inversion or rearrangement; `gap`: no read covers it), `reference`, `alternative`, `reads` of `depth`, `fraction`, `level` (`consensus` or `variant`: a mixed population, for example a mixed colony), `change` in words and `can_apply`.

## What to tell the user

First the consensus point changes and small indels, with positions and read support. Then the `variant` level ones (a minority of reads: a mixed colony, or sequencing error: the fraction says which). Then `clip` and `gap` entries: they mark larger changes (a deletion, insertion, inversion or duplication) that the reads do not spell out; they need the user's attention or another tool (`align` with a long read).

## Errors you may see

| Message | Fix |
|---|---|
| `Cannot read ... as a BAM file` | the file must be a BAM with its index (`.bai`) next to it |
| `holds reads for several sequences` | pass `--seqid` |
| `has no sequence called` | take the name from the list in the message |

## Example

```
variants demo/pBAD30_eilR.gb app/agent/examples/pBAD30_eilR_clone_reads.bam
```
