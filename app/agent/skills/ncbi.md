# ncbi: search NCBI and fetch a record

**Use when:** the user names a sequence they do not have as a file (a plasmid, a gene, a genome by accession such as `L09137` or `NC_000913.3`) and wants it from NCBI, or asks which records match some words. Same as Import from NCBI. This is the only command that uses the network.

## Privacy: tell the user what is sent

Only the **search words or accession** go to NCBI (eutils.ncbi.nlm.nih.gov), plus the user's saved NCBI API key if they have one (it is never shown). No file, no sequence of theirs is sent. Say so, and show the exact words, before the first call. Do not put anything private (a project name, an unpublished sequence) in the query.

## Commands

```
ncbi search WORDS [--db nuccore|protein] [--limit N]
ncbi get ACCESSION --out NEW_FILE [--db nuccore|protein] [--allow-large]
```

- `search`: records that match, each with `accession` (with version), `title`, `length`, `organism`, `updated`, `topology`, `source` (refseq, genbank, swiss_prot...). Up to 25 by default. Words like `E. coli K-12 complete genome` work; an accession finds that record.
- `get`: writes the record to a **new** file (never overwrites): a GenBank file with its annotation for `nuccore`, FASTA for `protein`. `ACCESSION` is the accession, with or without version (without: the current version). A record already fetched before comes from the local cache (`from_cache`). A record over 20 million bases is refused unless `--allow-large` (a download of minutes): ask the user first. Over 100 million bases it cannot be fetched here.

After `get`, use the file with `inspect`, `kld`, `pcr` and the rest like any other.

## What to tell the user

For a search, the matching records with their accession, length and organism, and which one you would pick and why (a RefSeq record, `NC_`/`NZ_`, is curated). For a fetch, the file written and its size. A record that NCBI could not return, or a connection problem, is an `ERROR:` line: relay it, do not retry in a loop (NCBI limits requests to 3 a second without a key).

## Errors you may see

| Message | Fix |
|---|---|
| `is not the accession of a ... record` | run `ncbi search` and take the accession from it |
| `needs --out` | give the new file |
| `has ... bases: a download of that size` | ask the user, then `--allow-large` |
| `already exists` | choose a new file name |
| `NCBI ...` (connection, too many requests) | wait and try again once; tell the user if it persists |

## Example

```
ncbi search pUC19 --limit 3
ncbi get L09137 --out example_puc19.gb
```
