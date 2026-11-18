# align: where a query sequence matches the plasmid

**Use when:** the user has a piece of DNA (a sequencing read, an insert, a suspected variant) and wants to know where and how it matches a plasmid, including deletions, insertions, inversions and pieces in another order. Same as the Local Align window. Read-only.

## Command

```
align PLASMID.gb QUERY
```

`QUERY` is the bases themselves (spaces, line breaks, digits and a FASTA header are ignored; U is read as T) or a GenBank/FASTA file (its first record). Both strands are searched; a circular plasmid is searched across its origin.

JSON: `summary` (one line: how many pieces, strands, structural events), and `pieces` in the order they have in the query: `query_start`..`query_end`, `start`..`end` on the plasmid (1-based), `strand` (-1: the query matches the reverse strand), `identity` (0..1), `score`, `mismatches`, `insertions`, `deletions` (counts) and `before`: how this piece follows the previous one (contiguous, deletion from the query, insertion in the query, inversion breakpoint, template jumps back).

## What to tell the user

Read `summary` and `before` first: they say whether the query is one clean match or a rearrangement. Give the identity and the positions. At most 30 pieces are listed (`note` says so). A query over 200,000 bp is refused; a protein is refused with a message.

## Errors you may see

| Message | Fix |
|---|---|
| `Nothing to align` / `No DNA found` | the query is empty |
| `This looks like a protein sequence` | give DNA |
| `Query too long` | cut it into parts |

## Example

```
align demo/pBAD30_eilR.gb ATGGCTGAAGCGCAAAATGATCCCCTGCTGCCGGGATACTCGTTTAATGCCCATCTGGTGGCGGGTTTAACGCCGATTGAGGCCAACGGTTATCTCGATTTTTTTATCGACCGACCGCTGGGAATGAAAGGTTATATTCTCAATCTCACCATTCGCGGTCAGGGGGTGGTGAAAAATCAGGGACGAGAATTTGTCTGCCGACCGGGTGATATTCTGCTGTTCCCGCCAGGAGAGATTCATCACTACGGTCGTCATCCGGAGGCTCGCGAATGGTATCACCAGT
```
