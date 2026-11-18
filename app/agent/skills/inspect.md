# inspect, feature, region: read a plasmid before you act on it

**Use when:** you need to know what is in a GenBank file (names, positions, which feature is a CDS, its protein), or the exact bases of a feature or a stretch. Always do this before building a mutation table.

## Commands

```
inspect FILE.gb          # length, circular, every feature: index, kind, label, start, end, strand, location; CDS: protein
feature FILE.gb NAME     # NAME = feature label (case ignored) or #INDEX: its bases in its own direction, and the protein of a CDS
region  FILE.gb START END   # bases START..END, 1-based inclusive
```

## Reading the output

- `start`/`end`: 1-based, inclusive, on the plasmid, whatever the strand. `strand` is 1 or -1.
- `protein.sequence` starts at the codon given by `/codon_start`; `*` is a stop. Amino-acid number N of a mutation such as `lacZ:D120A` is the Nth letter of this sequence: check that the letter is the one the user expects (here `D`).
- `protein.note` instead of a sequence means openDNA cannot translate that CDS as annotated (a join across the origin, a `/transl_except`): do not guess, tell the user.
- A gene and its CDS often share a label. `feature` picks the CDS; if there are several candidates it fails and lists them, and you pass `#INDEX`.

## Errors you may see

| Message | Fix |
|---|---|
| `File not found` | wrong path |
| `Cannot read ... as GenBank` | not a GenBank file: ask the user |
| `No feature called ...` | take a label from `inspect` |
| `N features are called ...` | use `#INDEX` from the list in the message |
| `Positions must satisfy ...` | stay within the length; ask for a region across the origin in two parts |

## Example

```
inspect demo/pBAD30_eilR.gb
feature demo/pBAD30_eilR.gb araC
region demo/pBAD30_eilR.gb 1 60
```
