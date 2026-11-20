# `tests/test_reference.py`

[← all test files](README.md) · source: [`tests/test_reference.py`](../../tests/test_reference.py)

## `test_off_target_search_covers_every_record`

**Result:** ✅ passed in 0.06 s · [line 79](../../tests/test_reference.py#L79) · **Area:** Off-target search (reference genome) · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** Without it, the primer and guide off-target checks could miss copies on the second or later chromosome, report the open sequence as its own off-target, show green after a search that timed out, or lose/double hits at window edges.

With a fake 3-record genome, checks the primer genome check finds copies on chr2 and chr3, excludes the open record, treats incomplete searches as not conclusive, flags a 3' copy as warning and a clean primer as ok; also checks multi-record GenBank references, the Cas9 guide finder counting every record, and windowed search equals whole search.

| # | Input | Expected output |
|---|---|---|
| 1 | 24 nt primer on a plasmid; genome chr1 (none), chr2 (primer at 501-524), chr3 (reverse complement) | level 'bad'; text mentions chr2: and chr3:, not chr1:, and '501–524'; reference count 3, label starts '3 sequences, ' |
| 2 | Open record is chr2 of the genome file | 'bad' citing chr3 only; with a one-record genome holding just chr2: ('info', 'this sequence is the genome') |
| 3 | search with seconds=-1; search cancelled after 2 progress calls | stopped 'time', searched 0, not complete; stopped 'cancelled', searched < 3, first progress names chr1 of 3 |
| 4 | Cached time-out result used by the check | 'warn' with '0 of 3' and 'not conclusive' |
| 5 | Genome with only a 3' 14 nt copy (after 3 mismatched bases) on chr3; then a genome without any copy | 'warn' citing chr3 only; then 'ok' with '3 sequences' |
| 6 | GenBank file of recA (circular, 400 bp) + recB (primer inside) | count 2, records [recA, recB], recA circular; whole hits only recB; searched 2; complete; sequences() names [recA, recB] |
| 7 | SpCas9 guide with copies (AGG, CGG PAM) on c2 and c3 of a 3-record FASTA | genome count 2, genome_name '3 sequences'; verdict 'bad' |
| 8 | 1000 bp record with primer at 395 (spans 400 window edge) and RC at 700; WINDOW set to 100 | windowed whole and hits equal unwindowed; 2 whole hits at 395 and 700 |

**Review notes:**

- The time-out section injects a result into the private _cache dict, coupling the test to implementation details.
- Random 'dna' filler uses a module-level RNG, so sequence content depends on the order tests run in that module (a chance extra match is unlikely but possible).

## `test_ncbi_genome_search_and_download`

**Result:** ✅ passed in 0.01 s · [line 184](../../tests/test_reference.py#L184) · **Area:** NCBI genome download · **Checked against:** hand-computed · **Strength:** 🟡 adequate

**Why it matters:** Without it, searching NCBI for a reference genome could build wrong URLs, and downloading one could produce a broken FASTA (records joined without newlines), leave partial files, or download the same genome repeatedly.

With the network faked, checks taxon-name and accession searches build the right NCBI Datasets URLs and parse the report, that a blank search is refused, and that a download unpacks chromosome files into one readable FASTA plus GFF, leaves no .part file and is cached on the second request.

| # | Input | Expected output |
|---|---|---|
| 1 | search_genomes('Madeup coli') with a faked report | URL contains 'taxon/Madeup%20coli/dataset_report' and 'reference_only=true'; result accession GCF_000000001.1, length 84, 2 sequences, 'Complete Genome' |
| 2 | search_genomes('gcf_000000001.1') | URL contains 'accession/GCF_000000001.1/dataset_report' (upper-cased) |
| 3 | search_genomes('  ') | NcbiError |
| 4 | download_genome with fake zip (chr1 without final newline, chr2, GFF) | not cached; progress reported; GFF contains 'abcA'; request asks GENOME_FASTA and GENOME_GFF; fasta path <acc>/<acc>.fna; Fasta lists [(chr1,60),(chr2,24)]; no *.part |
| 5 | download_genome again | cached True; still only 1 network call |

**Review notes:**

- Network is faked at ncbi.fetch/open_url, so a change in NCBI's real API or response format would not be caught.
- Unrelated to the multi-record off-target search that names this file; could live with other NCBI tests.
