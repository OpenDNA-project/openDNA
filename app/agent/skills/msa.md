# msa: align several sequences and compare them

**Use when:** the user has two or more sequences (variants of a plasmid, reads, homologous genes) and wants them aligned, with identity to the first one and a consensus. Same built-in aligner as the Multiple alignment tab. Read-only (writes a file only with `--out`).

## Command

```
msa A.gb B.gb [C.gb ...] [--out ALIGNED.fa]
```

One sequence per file (FASTA or GenBank), named by the record's name; duplicate names get a number. The first is the reference of `identity_to_first` (a gap counts as a difference). `--out` writes the alignment as FASTA to a **new** file (an existing file is an error: nothing is overwritten).

JSON: `width`, `consensus`, and per sequence `name`, `file`, `length`, `identity_to_first` (percent) and `aligned` (with `-` for gaps).

## Limits

The built-in aligner is a small progressive one: sequences up to 20,000 bp each, a few tens at most, no large rearrangements (an insertion or deletion longer than its band cannot be found). For bigger jobs the user needs MAFFT or Clustal Omega (openDNA's Terminal tab): say so rather than trusting a poor alignment.

## Example

```
msa demo/pBAD30_eilR.gb demo/pBAD30_RC.gb
```
