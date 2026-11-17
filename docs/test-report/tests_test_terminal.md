# `tests/test_terminal.py`

[← all test files](README.md) · source: [`tests/test_terminal.py`](../../tests/test_terminal.py)

## `test_pipelines_build_valid_bash`

**Result:** ✅ passed in 0.27 s · [line 25](../../tests/test_terminal.py#L25) · **Area:** Terminal and pipelines · **Checked against:** property/invariant · **Strength:** 🟡 adequate

**Why it matters:** A pipeline could type broken or wrong bash into the Terminal (an unfilled option, an unquoted path with a space, a tool it does not check), so a run fails or maps with the wrong settings.

Every preset of every pipeline in tools.toml builds a command with no placeholder left that bash accepts; the Illumina preset checks only the tools it calls, fills nested placeholders, quotes paths and leaves bash's own ${...} alone; the session setup, a wrapped line and Stop's command are valid bash. The flow chart the Terminal draws comes from tools.toml alone (runner.diagram).

| # | Input | Expected output |
|---|---|---|
| 1 | align, preset Illumina paired-end, 6 threads, reads 'my reads_R1.fq' | starts '( need fastp bwa samtools && cd ', has '-w 6 -q 20 -l 50' and 'bwa mem -t 6', the reads path quoted, '${R2:+-I "$R2"}' kept |
| 2 | align preset Nanopore | cleanup filtlong, min_len 50 (default), inputs reads1 + reference |
| 3 | every pipeline x preset | no {name} left; bash -n accepts it |
| 4 | session_setup + wrap("echo 'it''s' \| wc -c") + stop_children(42, 43) + install(['mafft']) | bash -n accepts it |
| 5 | default_output(align, 'pUC 19', 'my reads_R1.fq') | 'pUC_19_my reads.bam' beside the reads |
| 6 | diagram(align, Nanopore with cleanup none); diagram(assemble, Hybrid) | Reads, Reference, skipped 'Read cleanup: none', 'Aligner: minimap2_ont', samtools sort/index/flagstat, 'Open in openDNA'; assemble steps cleanup, unicycler (mode), quast once |

**Review notes:**

- Syntax only: no pipeline is run against real tools (they are not installed in CI).

## `test_open_paths_and_summary`

**Result:** ✅ passed in 0.02 s · [line 61](../../tests/test_terminal.py#L61) · **Area:** Terminal and pipelines · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** A result could open in the wrong place (an alignment as a document, a BAM in the browser), a WSL path could not be turned back into a Windows one, or the QC summary could misreport reads kept or % mapped.

WSL path conversion both ways, open_kind for each kind of file (alignment detected from its content), and the QC summary read from fastp.json, flagstat.txt and quast/report.tsv.

| # | Input | Expected output |
|---|---|---|
| 1 | shell_path('D:\data\x.fq'), windows_path('/mnt/d/data/x.fq') with WSL | /mnt/d/data/x.fq, D:\data\x.fq; a \\wsl$ path unchanged |
| 2 | x.bam, v.vcf.gz, aligned a.fa, plain p.fa, r.gb, fastp.html | track, track, msa, document, document, browser; --msa forces msa |
| 3 | fastp 1000 -> 920 reads, Q30 0.941; flagstat 987 mapped (98.70%); QUAST 3 contigs, 4641652 bp, N50 4600000 | 'fastp: 1,000 reads in, 920 kept (92%), Q30 94%', 'samtools flagstat: 987 reads mapped (98.70%)', 'QUAST: 3 contigs, 4,641,652 bp, N50 4,600,000 bp' |
| 4 | quast/report.html and quast/icarus_viewers/contig_size_viewer.html | only report.html listed as a report |
