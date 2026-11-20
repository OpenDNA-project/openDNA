# wsltools: the optional tools that run beside openDNA

**Use when:** the user wants pLannotate (the Barrick lab's plasmid annotator, used by `annotate` when it is installed) or the Terminal's bioinformatics tools (fastp, bwa, minimap2, samtools, unicycler, mafft...) set up, or `annotate` said pLannotate is `missing`. On Windows they run inside WSL; elsewhere in bash. Each is a separate program in its own conda environment, never part of openDNA.

## Command

```
wsltools status
wsltools install [--only plannotate|terminal]
```

`status`: JSON with `wsl` (can it be run), `plannotate` (`state`: `ready`, `no_database`, `missing` or `no_wsl`, and its `version`) and `terminal_tools` (each tool true when installed). It starts WSL: a few seconds.

`install`: downloads and installs what is missing (about 250 MB of pLannotate databases and several hundred MB of programs: BLAST+, DIAMOND, Infernal, the aligners), a few minutes, its text on stderr. It is safe to run again. It never installs WSL itself: that needs an administrator and a restart, so when `wsl` is false tell the user to run `wsl --install` in an administrator terminal, restart, and ask again. **Ask the user before `install`:** it downloads a lot and writes into their WSL home folder.

## What to tell the user

What is installed and what is not. After `install`, that pLannotate now takes part in `annotate`, and that it is the Barrick lab's work (McGuffie & Barrick, Nucleic Acids Res 2021, doi:10.1093/nar/gkab374, https://github.com/barricklab/pLannotate).

## Example

```
wsltools status
```
