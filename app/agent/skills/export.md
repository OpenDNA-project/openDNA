# export: write the sequence in another file format

**Use when:** the user wants a plasmid as FASTA, GFF3, BED, an NCBI feature table (`tbl`, for table2asn), the proteins of its CDS as FASTA, or as a plain GenBank file (which also converts a Benchling, APE or EMBL file to GenBank).

## Command

```
export FILE FORMAT NEW_FILE
```

`FORMAT`: `genbank`, `fasta` (with `[topology=circular]` in the header for a circular plasmid), `gff3`, `bed` (0-based starts), `tbl`, `proteins` (one FASTA record per CDS, final stop removed, translated with the CDS's own table). `NEW_FILE` must not exist: nothing is overwritten. JSON: `file`, `bytes`, `features`.

## Example

The example writes a file into the current folder.

```
export demo/pBAD30_eilR.gb fasta example_export.fa
```
