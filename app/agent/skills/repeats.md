# repeats: stretches shared inside one sequence or between two

**Use when:** the user asks about repeated sequence (a recombination risk: 20 bases and more can recombine in the cell), or what two plasmids have in common. Same as the Compare view of Local Align. Read-only.

## Command

```
repeats A.gb [B.gb] [--min-len N] [--mismatches N]
```

A file with several records (the Vibrio file, a multi-FASTA) is compared one record at a time: name it `FILE#NAME` (the names are in the error message, or in `split`'s answer). Circular sequences are read as circles: a stretch may run through the origin: its start is inside the sequence and its end may be above the length (start 590, length 30 on a 600 bp circle: end 619, meaning bases 590..600 then 1..19); a stretch is cut where it reaches 1,000 bases past an origin, and a warning says so. With one file the sequence is compared with itself (each pair once, the trivial self-match left out). With two, the stretches they share; two different files with the same bases share their whole length (that stretch is listed first, `identical_sequences` is true). `direct`: same orientation; `inverted`: one is the reverse complement of the other. `--min-len` (default 12; the window accepts 8 to 500), `--mismatches` (default 0; 0 to 3).

JSON: `compared` (the two file names), `repeats`, longest first: `kind`, `start_in_a`..`end_in_a`, `start_in_b`..`end_in_b` (1-based), `length`, `differences`; `note` when the list was capped.

## What to tell the user

The longest repeats and where they are; 20 bases or more is a recombination risk in a cell. Short repeats (12 to 19) are rarely a problem: say so rather than alarming them.

## Example

```
repeats demo/pBAD30_eilR.gb --min-len 20
```
