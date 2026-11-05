# `tests/test_reference.py`

[← all test files](README.md) · source: [`tests/test_reference.py`](../../tests/test_reference.py)

## `test_genome_checks`

**Result:** ✅ passed in 0.41 s · [line 81](../../tests/test_reference.py#L81) · **Area:** Reference genome: primer off-target check · **Checked against:** hand-computed · **Strength:** 🟢 strong

The primer off-target check finds copies on other records, never compares the open record with itself, reports an incomplete/cancelled search as inconclusive (yellow), runs on a worker thread in the dialog, and grades 3' copies yellow and clean primers green.

| # | Input | Expected output |
|---|---|---|
| 1 | 3-record genome (chr2/chr3 hold primer); check primer on a plasmid | level 'bad'; text has chr2:, chr3:, '501–524', not chr1:; ref.count 3, label '3 sequences, ...' |
| 2 | open chr2 as the document | level 'bad', chr3: only, chr2: excluded (self); a single-record file that is the open seq -> ('info','this sequence is the genome') |
| 3 | search with seconds=-1 (no time) and a cancel after 2 progress calls | stopped=='time', searched 0, not complete; cancelled: stopped=='cancelled', searched<3, first progress chr1 with total 3 |
| 4 | cache a time-stopped 0-hit result | level 'warn', '0 of 3', 'not conclusive' (never green no-match) |
| 5 | PrimerDialog._update_report on a worker thread; wait up to 10s | report contains chr2: and chr3: |
| 6 | rewrite genome: 3' 14bp copy on chr3; then clean | warn with chr3: only; then ok with '3 sequences' |

**Review notes:**

- The dialog section (lines 123-130) polls d.report.text() in a wall-clock loop up to 10s for a worker-thread result; a slow machine could time out and the following assertion (line 129) would then fail even though the behaviour is correct. Mutates shared primer_panel settings but via a monkeypatched throw-away QSettings.

## `test_reference_searches_every_record`

**Result:** ✅ passed in 0.01 s · [line 145](../../tests/test_reference.py#L145) · **Area:** Reference genome: multi-record search / windowing · **Checked against:** round-trip/self-consistency · **Strength:** 🟢 strong

A multi-record GenBank reference is read one record at a time (refusing to open as a single Document), the guide finder counts sites across every record, and a long record searched in overlapping windows returns the same sites without loss or duplication.

| # | Input | Expected output |
|---|---|---|
| 1 | concatenated recA(circular)+recB GenBank as one file | Document(both) raises OSError match 'reference genome'; Reference(both).count 2, records ['recA','recB'], recA circular True; search(PRIMER) whole==['recB'], searched 2, complete; sequences() names ['recA','recB'] |
| 2 | 3-record FASTA c1/c2/c3 (c2,c3 hold guide); cas9.scan max_mm 2 with Reference | genome[0]==2, genome_name '3 sequences'; verdict 'bad' |
| 3 | 1000bp seq with PRIMER at 395 (spans window edge 400) and rc(PRIMER) at 700; search with WINDOW=100 vs default | windowed.whole == whole.whole (2 sites); windowed.hits == whole.hits; whole positions == [395,700] |

## `test_ncbi_members_are_joined_with_their_newlines_and_the_reference_is_the_fasta_itself`

**Result:** ✅ passed in 10.11 s · [line 190](../../tests/test_reference.py#L190) · **Area:** GUI: NCBI genome 'use as reference' · **Checked against:** hand-computed · **Strength:** 🟢 strong

Downloading a multi-sequence NCBI assembly joins its member FASTA files (handling a missing final newline) and sets the combined FASTA itself as the reference genome, keeping every record.

| # | Input | Expected output |
|---|---|---|
| 1 | stubbed NCBI report (2-sequence assembly) + zip with chr1 (no trailing newline) and chr2; GenomeDialog search 'Madeup multi'; wait for assemblies | reference_button enabled (2 sequences) |
| 2 | use_as_reference(); wait for 'reference genome' in status | settings primers/genome == the combined .fna path; '2 sequences' in status; Fasta(fasta).sequences() names == ['chr1','chr2'] |

**Review notes:**

- Wall-clock wait loops up to 10s each (lines 212-213, 218-219) for worker-thread search/download; network is stubbed so no real access, but a very slow machine could time out and fail the following assertion.

## `test_other_records`

**Result:** ✅ passed in 0.18 s · [line 229](../../tests/test_reference.py#L229) · **Area:** GUI: multi-contig record handling · **Checked against:** hand-computed · **Strength:** 🟢 strong

Contig names match ignoring 'chr' and case, tracks follow the switched record (a mismatching contig is still loaded after a warning), and saving one record as GenBank notes the others are left behind.

| # | Input | Expected output |
|---|---|---|
| 1 | match_contig cases | ['chr1','chr2'],'chr2'->'chr2'; ['1','2'],'chr2'->'2'; ['Chr1','Chr2'],'chr1'->'Chr1'; ['1','chr1'],'chrX'->None; ['chr1','1'],'chr1'->'chr1'; ['chr7'],'chr1'->None |
| 2 | load 2-contig FASTA; load matching chr1 GFF3; load chr7-only GFF3 | chr1 track no warning; chr7 track loaded (seqid chr7) after a warning naming 'chr7' and 'chr1' |
| 3 | switch_sequence('chr2') | editor.fasta_record 'chr2'; tracks == [('t.gff3','chr2')] (follows); 'only7.gff3' in status bar |
| 4 | 3-contig FASTA; other_records_note; write chr1.gb | note starts ' (only chr1: the other 2 sequences of g.fa'; write ok; status 'the other 2' |

**Review notes:**

- Monkeypatches QMessageBox.warning by direct attribute assignment restored in a finally block for the tracks section (lines 248-258), but an exception in the earlier match_contig asserts (before the try) would not need restoration; fine as written.
