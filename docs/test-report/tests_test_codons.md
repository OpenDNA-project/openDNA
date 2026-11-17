# `tests/test_codons.py`

[← all test files](README.md) · source: [`tests/test_codons.py`](../../tests/test_codons.py)

## `test_recoding_preserves_the_protein_and_obeys_the_constraints`

**Result:** ✅ passed in 0.10 s · [line 79](../../tests/test_codons.py#L79) · **Area:** Codon recoding · **Checked against:** property/invariant · **Strength:** 🟢 strong

**Why it matters:** Without it, recoding a CDS could silently change the protein the user is cloning, or leave a forbidden restriction site, a bad GC window or a hidden RBS in the synthesised gene.

Checks the invariant that a recoded CDS translates to the same protein on either strand and with either method, that the recoded sequence obeys all constraints at once (no forbidden sites including across its two edges, GC window, run length, no RBS), and that an impossible rule is reported with its position instead of being ignored.

| # | Input | Expected output |
|---|---|---|
| 1 | a real CDS of the demo plasmid on each strand, methods MATCHED and HIGHEST, made-up 'uneven' table | sequence changes, same length, same translation under the CDS's genetic code |
| 2 | minus-strand CDS recoded with avoid=8 common enzymes, gc_range 40-60 over 50 bp windows, max_run 5, rbs=True, flanks given | unmet == []; no avoided site anywhere in left+seq+right; no RBS motif; longest homopolymer run <= 5; every 50 bp GC window between 40 and 60 |
| 3 | GCT + ATGTGG + GCT with a forbidden site 'ATGTGG' (Met-Trp: one codon each) | unmet == ['Impossible site at bases 4-9'] and the translation is still AMWA |

## `test_codon_choice_follows_the_table_and_the_edit_is_one_undo_step`

**Result:** ✅ passed in 0.03 s · [line 113](../../tests/test_codons.py#L113) · **Area:** Codon recoding · **Checked against:** in-test reference model · **Strength:** 🟢 strong

**Why it matters:** Without it, recoding could ignore the chosen usage table (defeating its purpose), be irreproducible between runs, or leave the document in a state undo cannot restore.

Checks that 'highest frequency' always picks the most used codon, that 'frequency-matched' reproduces the table's shares statistically, that a zero-frequency codon is avoided unless it is the only option, that a seed makes the result reproducible, that a non-codon selection is refused, and that applying a recoding is a single undo step which a locked sequence refuses.

| # | Input | Expected output |
|---|---|---|
| 1 | 300 random codons recoded with method HIGHEST against the uneven table | every codon equals the most-used synonym for its amino acid, computed in the test from the table |
| 2 | CTG x 3000 with a Leu table CTG 50, others 10 each | CTG share ~0.5 (+-0.04), CTA share ~0.1 (+-0.03) |
| 3 | CTG x 200 with CTG set to 0 / TGG x 10 with TGG set to 0 | CTG never chosen; TGG kept because Trp has no synonym |
| 4 | the same CDS recoded with seed 3 twice, then seed 4 | identical with the same seed, different with another; ValueError 'whole number of codons' for 'ATGAA' |
| 5 | selections of 10 bases, a caret, and 30 bases | 'whole number of codons' message; 'Select a CDS' message; a 30 bp selection read as strand +1 starting at 0 with those bases |
| 6 | applying a recoding to a minus-strand CDS via overwrite_bases, then undo | sequence changes but keeps its length, the CDS reads the recoded bases, all feature kinds/locations/labels unchanged, one undo restores the original bases |
| 7 | the same edit with dna_locked = True | overwrite_bases returns False and the sequence is unchanged |

**Review notes:**

- The frequency-matched tolerance (+-0.04 on 3000 codons) is loose enough that a moderately skewed sampler would pass.

## `test_comparing_a_cds_under_two_genetic_codes`

**Result:** ✅ passed in 0.00 s · [line 176](../../tests/test_codons.py#L176) · **Area:** Genetic code comparison · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** Without it, the comparison that tells a user a gene will be truncated in another host (TGA read as stop instead of Trp) could point at the wrong codon or miscount the final stop.

Checks that comparing a CDS under two genetic codes reports the gained and lost internal stops at the right codon numbers, that the terminal stop is not counted as internal, that identical codes differ nowhere, and that the rendered comparison numbers its rows and marks differences.

| # | Input | Expected output |
|---|---|---|
| 1 | ATG TGA GCT TAA compared Mycoplasma(4) -> standard(1) | proteins 'MWA*' and 'M*A*'; gained == [2], lost == [], differ == [2] |
| 2 | the same compared standard -> Mycoplasma | lost == [2], gained == [] |
| 3 | ATGGCTTGA under codes 1 and 11 | differ, gained and lost all empty |
| 4 | ATGGCTTGA under codes 1 and 4 | differ == [3], gained == [] - the last codon is a terminal stop, not an internal one |
| 5 | comparison_html of ATG + TGA x 12 + TAA, codes 4 and 1 | HTML names both codes, contains the difference marker, a stop symbol and the row number 14 |

**Review notes:**

- The HTML assertion is a weak substring check ('14' could match anything).

## `test_usage_tables_their_files_and_what_they_say_about_a_gene`

**Result:** ✅ passed in 0.15 s · [line 201](../../tests/test_codons.py#L201) · **Area:** Codon usage tables and genetic codes · **Checked against:** round-trip/self-consistency · **Strength:** 🟢 strong

**Why it matters:** Without it, a user's imported or self-made codon table or genetic code could be lost, mis-parsed, or silently saved nowhere, and the adaptation index and checkup notes computed from it could be wrong.

Checks that no table ships and nothing breaks without one; CSV import/export including the genetic code column and the files that must be refused; where user tables live and how failures to write are reported; the CAI and rare-codon maths; the checkup note that uses the active table; and saving, deleting and translating with recoded genetic codes.

| # | Input | Expected output |
|---|---|---|
| 1 | fresh data folder | all_tables() == {}, active() is None, checkup.codon_notes(...) == [] with no error |
| 2 | write_csv of a flat table (AAA=40) then read_csv | name 'My organism', source 'a paper', AAA 40.0, 64 codons |
| 3 | a headerless codon/frequency list using U instead of T | TTT 0 and GGG 63 (accepted) |
| 4 | CSV with 1 row / with codon 'XYZ' / with amino acid 'B' / with no aa column | ValueError '64 codons'; ValueError 'not a codon'; ValueError 'not an amino acid'; aas is None but the table is accepted |
| 5 | write_csv with a saved recoded code (AGA, AGG removed from table 11) | file contains '# code: Export test code', header 'codon,aa,per_thousand,start', 'AGA,X,' and 'ATG,M,'; reading it back restores aas, starts, code_name and the frequencies |
| 6 | duplicate + save_user of a table named 'Recoded host' | Recoded_host.json under OPENDNA_DATA_DIR; listed once with note 'Copy of Original.'; a blank name raises 'needs a name'; a corrupt JSON file in the folder is ignored; after delete_user the active table is cleared |
| 7 | fileio.replace monkeypatched to raise PermissionError | save_user and genetic_code.save_code raise ValueError 'could not be saved...' and nothing is registered |
| 8 | CAI with a table where CTG=50 and CTA=1 | cai('CTGCTGCTG')==1.0, cai('CTACTA')==1/50, cai('CTGCTA')==sqrt(1/50) (geometric mean), cai('ATGTGGTAA') is None; rare_codons finds (3,'CTA'); a never-used codon scores 0.01; adaptiveness gives CTA 0 and omits stops |
| 9 | active table 'Host' (CTG 50, CTA 0.5), checkup.codon_notes | no note for ATG+CTG*20+TAA; one note for CTA*20 containing 'adaptation index 0.01 for Host' and '20 rare codons (CTA)'; a missing active table still returns a list |
| 10 | built-in recoded codes 101/102 and a user code saving/deleting | removed_codons(101)==['TCA','TCG','TAG'], 7 for 102; translate('ATGTCGTAGTCATAA',101)=='MXXX*'; user code id >= FIRST_USER with its JSON file, translate('AGAAGG')=='XR'; duplicate name, built-in id and a 4-letter table all raise; after delete the id and file are gone |
| 11 | genetic code tool setting switched to code 101 and back | protein_problems('ATGTCGAAATCATAA') reports 'uses codons this genetic code does not use: TCA, TCG', and nothing under the default code |

**Review notes:**

- Far too large: CSV parsing, user-file storage, write failures, CAI maths, checkup notes and genetic-code management are six independent features in one test; a single failure stops the rest from running.
- Mutates global state (active table, tool settings) and relies on try/finally to clean up.
