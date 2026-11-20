# split and merge: files of several records

**Use when:** a file holds several sequences (a multi-FASTA; a GenBank file of several LOCUS records such as a bacterium's two chromosomes, or a genome and its plasmids) and the user wants them apart, or several files and the user wants them in one. Nothing here changes a sequence.

## First: you may not need to split

Every command that takes a file accepts `FILE#NAME`, the record called NAME of a multi-record file (quote it in a shell: `"genome.gb#NC_002506"`). `msa` takes a multi-record file whole. A bare file name is refused when the file holds more than one record, and the error lists the names.

## Commands

```
split FILE [--out NEW_FOLDER]
merge A B [C ...] --out NEW_FILE [--names X,Y,Z]
```

- `split`: one file per record, in a new folder (`FILE_split` next to the file by default; an existing non-empty folder is an error). A GenBank file gives `.gb` files with their features; a FASTA file gives `.fasta` files. A file name is made of the record name, with unsafe characters replaced and a number added when two would be the same (a warning says so). JSON: `out`, `files` (`record`, `file`, `length`, `circular`).
- `merge`: every sequence of every file (a multi-record file counts for all its records, `FILE#NAME` for one) into one new file, in the order given. The extension of `--out` picks the format: `.fasta`/`.fa` writes only the bases (features and circular topology are lost: a warning says so), anything else is GenBank with features and topology kept. `--names A,B,C` sets the record names (one per sequence, in order). Without it a record keeps its own name; when two have the same name the second takes its file's name (`pBAD30_RC` for `pBAD30_RC.gb`), then `_2`, `_3` if that is taken too (a warning says so). The file must not exist. JSON: `file`, `bytes`, `records` (`name`, `source`, `length`, `circular`).

## What to tell the user

Which files or records were written, and every warning. If they only wanted one record analysed, you did not need `split`: use `FILE#NAME`.

## Example

The examples write into the current folder.

```
split demo/msa_demo.fasta --out example_split
merge demo/pBAD30_eilR.gb demo/pBAD30_RC.gb --out example_merged.gb
region demo/msa_demo.fasta#clone_1 1 20
```
