# `tests/test_roundtrip.py`

[← all test files](README.md) · source: [`tests/test_roundtrip.py`](../../tests/test_roundtrip.py)

## `test_real_files_round_trip_and_agree_with_biopython`

**Result:** ✅ passed in 1.90 s · [line 132](../../tests/test_roundtrip.py#L132) · **Area:** GenBank round trip · **Checked against:** independent implementation · **Strength:** 🟢 strong

**Why it matters:** Without it, opening and saving a real plasmid or genome could silently drop or alter features, qualifiers or bases, or write a file other tools read differently.

Loads every real and demo file (the two plasmids, the .dna file, the demo_*.gb files and the 4.6 Mb E. coli genome), saves as GenBank, reloads and compares everything; checks repeated saves are byte-identical and that Biopython reads our output the same. Also checks the .dna file against its GenBank export and that an edit inside a feature moves it, is saved, and undoes.

| # | Input | Expected output |
|---|---|---|
| 1 | Each present file: Document(path) -> save_genbank -> reload | name, definition, topology, bases (case kept) and every feature (kind, location, qualifiers) equal; spans and strands equal; reloaded doc not modified |
| 2 | BW25113 genome file | length > 4,500,000, circular, more than 4000 features |
| 3 | Save, reload, save, reload, save | 1.gb, 2.gb, 3.gb byte-identical |
| 4 | Our saved file read by Biopython (all but the genome) | same upper-case sequence, topology, id/name, feature count, type, location parts and qualifiers in order |
| 5 | pBAD30_eilR.gb and pFG051_Chlor.gb | circular, have features, lengths 5936 and 3780 |
| 6 | pFG051_Chlor.dna vs pFG051_Chlor.gb | both 3780 bp circular, same bases ignoring case; .dna has 15 primers; every CDS has a direction arrow; at least one misc_feature has none |
| 7 | is_directional rules | misc_feature not directional, CDS directional; /directional=yes makes misc_feature directional, /directional=no makes promoter non-directional |
| 8 | Insert GATTACA 10 bases into a >20 bp single-span feature of the .dna, save as GenBank | feature end grows by 7; doc not modified after save; Biopython reads same sequence, feature count and new start/end |
| 9 | undo | feature span back to the original (start, end) |

**Review notes:**

- Skips entirely if the demo files are missing; the .dna sections return early without them.
- The source id 'genome_dna' actually refers to the GenBank file BW25113_2026.gb (GENOME), not a .dna file: misleading name.
- Saving and reloading the 4.6 Mb genome three times makes this one of the slower tests.
- The is_directional checks are unrelated to round trips and could live in a feature-direction test.

**Fixed by the audit of November 2026:**

- 13 Nov 2026: gap closed. The proposed checks were added, pass on the current code and fail under the mutation that slipped through before. Was: 'load, save, reload is identical' and 'Biopython agrees' both compare against the document our own reader built from the ORIGINAL file; the original file is never compared with an independent reading of it. A GenBank reader that loses a qualifier, a feature or bases on open then writes the loss, reloads the loss and Biopython reads the loss: everything agrees. This is exactly the role ('opening and saving a real plasmid or genome could silently drop or alter features, qualifiers'). \| The genome (the only file with /EC_number, most /locus_tag) is excluded from every Biopython comparison.

## `test_unusual_records_round_trip`

**Result:** ✅ passed in 0.08 s · [line 214](../../tests/test_roundtrip.py#L214) · **Area:** GenBank round trip · **Checked against:** independent implementation · **Strength:** 🟢 strong

**Why it matters:** Without it, unusual but legal locations (joins, fuzzy ends, order(), across the origin, between sites), tiny or empty sequences, or lower-case bases could be corrupted on save without anyone noticing.

Saves made-up records that hit the edges of the GenBank format, reloads them with our reader and with Biopython, and checks every location form is written back exactly as given. Also covers topology toggling with undo, tiny/empty sequences, a gap letter and case preservation.

| # | Input | Expected output |
|---|---|---|
| 1 | 400 bp random 'rich' record with join, complement(join), <100..>120, <150..180, 200..>230, order(), complement gene with /pseudo; save, reload, save | snapshot equal; the two saved files byte-identical |
| 2 | Rich record circular and linear, read by Biopython | Biopython agrees on everything; 'fuzzy both ends' has BeforePosition start and AfterPosition end; 'complement of a join' strand -1 |
| 3 | 12 location strings (10..20, complement, joins with mixed strands, <10..20, 10..>20, order, three-part join, join(350..400,1..20) across origin) on a 400 bp circular seq | after reload location, spans and strand equal; location string written exactly as given; Biopython part spans equal ours |
| 4 | Between location 20^21 on an 80 bp linear seq | location '20^21', spans [(20,20)]; Biopython start=end=20 |
| 5 | 105 bp sequence, circular and linear; reload, set_topology(not circular), save | LOCUS line has ' circular '/' linear ' accordingly; reload shows flipped topology; one undo restores it |
| 6 | Sequences 'A', 'g', 'N', 'AC', 'acgtNRYKMSWBDHVU', each linear and circular with a full-length feature | bases, length, topology, name 'tiny' and snapshot preserved; Biopython reads upper-cased bases and 1 feature |
| 7 | Empty sequence (both topologies); 'ACGT' deleted down to nothing | reloads length 0, '' , same topology and name; empty save succeeds and clears modified |
| 8 | 'ACGT-ACGT' | gap letter '-' preserved on reload |
| 9 | Hand-written 60 bp GenBank with mixed case and IUPAC letters | sequence read exactly as typed, slice (8,14) = 'acGTAC'; ORIGIN line written in 10-base blocks with the same case; second save byte-identical; Biopython upper-case equal |
| 10 | insert 'TTT' at 0 | sequence starts 'TTTacg' (case of typed text kept) |

**Review notes:**

- Very long test mixing many unrelated edge cases; a single failure hides the rest. Could be split into location forms, topology, tiny sequences and case.
- The empty-sequence-after-delete case only checks the file exists, not that it reloads.

**Fixed by the audit of November 2026:**

- 13 Nov 2026: gap closed. The proposed checks were added, pass on the current code and fail under the mutation that slipped through before. Was: tests.json case "'ACGT' deleted down to nothing => reloads length 0, '', same topology and name" is not asserted: the shrunk document's file is never reloaded, so a writer that falls back to the original/stale bases (or writes a malformed LOCUS) when length drops to 0 after an edit would pass. Only the new_sequence('empty','') path is reloaded.

## `test_qualifier_text_survives_and_is_written_as_documented`

**Result:** ✅ passed in 0.11 s · [line 340](../../tests/test_roundtrip.py#L340) · **Area:** GenBank qualifiers · **Checked against:** independent implementation · **Strength:** 🟢 strong

**Why it matters:** Without it, notes with quotes, unicode or long text could be mangled, wrapped badly or gain spaces when saved, so other tools (or a reload) would show different annotation text.

Writes a feature with awkward qualifier values and checks the reload equals the documented canonical form (whitespace runs collapsed, translation whitespace removed), lines fit 80 bytes, the file is stable on resave, and Biopython reads the same values. Also checks very long words are never cut and full last lines do not gain a trailing space.

| # | Input | Expected output |
|---|---|---|
| 1 | Qualifiers: doubled quotes, value of only quotes, single '"', newlines/tabs/space runs, ~500 char text, unicode text, quotes near the wrap column, product with =;/()[]{}, translation with whitespace, /pseudo flag | reloaded qualifiers = given with whitespace runs collapsed to one space (translation: all whitespace removed); document holds the same canonical form |
| 2 | The written file | every line at most 80 bytes; contains '日本語' and '/pseudo\n' |
| 3 | Reload and save again | byte-identical |
| 4 | Biopython reads the file | notes equal the canonical values, translation 'MKVLLAGGT*', pseudo [''], product unchanged |
| 5 | Note 'before <300 x> after' and a bare 300-char word | both read back unchanged (word never split) |
| 6 | Note 'sequence: ' + n random bases for n = 40..89 | our reader and Biopython both read the exact note (no trailing space); lines are at most 79 wide unless they are a single unbreakable word |

## `test_primers_survive_save`

**Result:** ✅ passed in 0.03 s · [line 410](../../tests/test_roundtrip.py#L410) · **Area:** Primers in saved files · **Checked against:** round-trip/self-consistency · **Strength:** 🟡 adequate

**Why it matters:** Without it, primers attached to a plasmid (from a .dna file or imported from the database) could vanish or lose their full sequence when the file is saved as GenBank.

Checks primer_bind features passed to save_genbank are written with name and full sequence (including 5' tails) without being added to the document, that a .dna file's primer list survives a GenBank save, and that Editor.save writes only primers imported to the file without adding an undo step.

| # | Input | Expected output |
|---|---|---|
| 1 | save_genbank with fwd 5..24 and rev complement(60..80), notes carrying tailed sequences | reloaded primer_bind features: ('fwd','TTTTTACGTTGCATGCCGATTAG','5..24') and ('rev','cgatccggCATGGCAACTTGGCTAATC','complement(60..80)'); the document itself still has 0 features |
| 2 | pFG051_Chlor.dna saved as GenBank | at least 5 primers; reloaded primer_bind list equals the source's and maps each name to its full sequence |
| 3 | Primer 'myprimer' = bases 5..27 added only to the primer database, Editor.save | no primer_bind in the saved file |
| 4 | import_primers([myprimer]) then Editor.save | import returns 1; saved file has primer_bind features all named myprimer with that sequence, one at '6..27'; doc not modified and the undo history has exactly one step (the import) |

**Review notes:**

- The .dna section is skipped if the demo .dna file is missing.
- Overlaps with test_multi_record's .dna-primers-as-primer_bind check and with primer tests elsewhere.

## `test_fasta_files_are_read_like_biopython`

**Result:** ✅ passed in 0.03 s · [line 453](../../tests/test_roundtrip.py#L453) · **Area:** FASTA reading · **Checked against:** independent implementation · **Strength:** 🟢 strong

**Why it matters:** Without it, opening a multi-sequence FASTA (genomes, Windows line endings, BOM, lower case) could return wrong bases at line breaks or include stray characters, corrupting what the user sees and edits.

Checks the indexed FASTA reader against Biopython for record list, full fetches, clipped fetches and slices crossing line breaks; that the .fai index is written and reused; that line-ending, blank-line, case and long-line variants give the same bases; and that a byte order mark is ignored.

| # | Input | Expected output |
|---|---|---|
| 1 | 4 records (250, 130, 61, 0 bp) wrapped at 60 | sequences() = [(chr1,250),(chr2,130),(chr3,61),(empty_one,0)]; each fetch equals Biopython; end beyond length clipped; Document.from_fasta per record equals Biopython |
| 2 | chr1 slices (0,1),(59,61),(55,130),(249,250),(119,121),(0,250) | equal Biopython slices |
| 3 | Document.from_fasta('chr2') | name chr2, length 130, linear |
| 4 | After opening | multi.fa.fai exists; reopening gives the same list |
| 5 | 8 variants: CRLF, no final newline, trailing blank lines, lower case, CRLF without final newline, CRLF with blank lines, mixed case, one long line | records [(chr1,150),(chr2,77)]; bases equal Biopython and the original ignoring case (no \r, no blank lines); fetch(148,152) clipped to 2 bases |
| 6 | UTF-8 BOM before '>chr1 first' | sequences [(chr1,14),(chr2,12)]; chr1 'ACGTACGTACGTAC'; chr2 starts 'TTTT'; from_fasta works |

## `test_fasta_writing_is_the_string_edit_of_the_file`

**Result:** ✅ passed in 0.02 s · [line 505](../../tests/test_roundtrip.py#L505) · **Area:** FASTA writing · **Checked against:** in-test reference model · **Strength:** 🟢 strong

**Why it matters:** Without it, saving an edited chromosome back into a multi-FASTA could damage the other records, lose CRLF or a BOM, leave a stale index, or write the wrong edited bases.

Applies edits to one record of a FASTA and checks the written file is exactly the original text with only that record replaced (computed independently on a plain string), both to another file and in place; other records stay byte-identical with CRLF or BOM preserved, and exporting does not clear the modified flag.

| # | Input | Expected output |
|---|---|---|
| 1 | Insert GATTACA in chr2, write_fasta to export.fa | Biopython ids [chr1,chr2,chr3,empty_one], chr1 description kept, chr2 = edited sequence, chr1/chr3 unchanged; doc still modified |
| 2 | insert 5 'GGGGGGG', delete 100..110, replace 0..3 by 'tt' on chr2; same edits on a Python string | doc sequence equals the string model; written file bytes equal the expected FASTA text; source file untouched |
| 3 | Build index, then write_fasta in place | file equals expected; new Fasta index returns edited chr2 and unchanged chr1/chr3; no *.tmp left |
| 4 | CRLF file, insert AAAA at 0 of chr2, write in place | bytes before '>chr2' and from '>chr3' identical; chr2 = AAAA+original; Biopython agrees |
| 5 | CRLF file, replace chr2 with 70 G | every line ends CRLF; output starts with chr1 block then '>chr2\r\n' + 60 G + CRLF + 10 G + CRLF |
| 6 | BOM file, replace chr1 with 'GGGG' | exact bytes '\xef\xbb\xbf>chr1 first\nGGGG\n>chr2\nTTTTTTTTTT\nTT\n' |

**Review notes:**

- Partly overlaps test_multi_record::test_records_are_added_reordered_and_removed_byte_for_byte (CRLF/BOM byte preservation), but that one tests reorder/append rather than edit-in-place.
