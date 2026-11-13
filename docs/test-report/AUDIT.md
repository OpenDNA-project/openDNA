# Audit of the unit tests

[← unit-test report](README.md)

Built by `scripts/test_report.py` on 2026-11-13 from `tests.json` and the last run. Every test of the project is listed below, Python (`tests/`) and Rust (`core/src`), with what would break unnoticed without it, what its expected answer is checked against, its strength and the reviewer's notes. The Python suite holds no UI tests: windows, dialogs and drawing are checked by hand before a release.

## Totals

| Suite | Tests | Passed | Failed / not run | Strong | Adequate | Weak | Not reviewed | With review notes |
|---|---|---|---|---|---|---|---|---|
| Python | 94 | 94 | 0 | 68 | 25 | 1 | 0 | 72 |
| Rust | 65 | 65 | 0 | 41 | 22 | 2 | 0 | 44 |
| **All** | 159 | 159 | 0 | 109 | 47 | 3 | 0 | 116 |

## Weak or unreviewed tests

- `tests/test_protein.py::test_mass_pi_charge_and_numbered_layout` (weak): pI and charge are checked only for sign/range on homopolymers: a wrong pKa table or a broken pI search could still pass. An exact pI and mass for a known protein (e.g. against Biopython ProtParam, which is installed in .venv) would make it strong.; Mass is checked for one dipeptide only; a wrong residue mass for any amino acid other than G would go unnoticed.; The index/offset check is a round trip only; offset_of and index_at that were wrong in the same way would pass.
- `core/src/align.rs::tests::max_hits_still_limits_the_number_of_pieces_and_keeps_the_best` (weak): Name promises 'keeps the best' but lines 666-667 only check counts; nothing compares align(..,2) with the top 2 of align(..,10), so returning the 2 worst pieces would pass.
- `core/src/primers.rs::tests::tm_and_tm_range_of_plain_and_degenerate_primers` (weak): (was tm_is_sane) 10 °C window only; exact values are checked in Python by tests/test_core_sequence.py::test_melting_temperature_matches_santalucia_and_biopython. The code comment in the test still points to tests/test_primers.py, which no longer exists.; (was tm_is_sane) Comment says 'more salt or more primer = higher Tm' but primer concentration is never varied; the self-complementary (symmetric) branch is untested.; (was tm_is_sane) Superseded by tests/test_core_sequence.py::test_melting_temperature_matches_santalucia_and_biopython (exact Biopython values, the same input checks, salt/Mg/dNTP/primer direction); it only remains as a fast Rust smoke test.; (was tm_range_of_a_plain_primer_is_its_tm_and_bad_letters_give_none) Duplicates the end of tests/test_core_sequence.py::test_melting_temperature_matches_santalucia_and_biopython (tm(degenerate) is None, tm_range(plain) == (tm, tm), tm_range('ACGTX') is None); the only extra is tm_range('A'), and tm_is_sane already checks tm('ACGTN') is None.; (was tm_range_is_exactly_the_lowest_and_highest_tm_of_the_expanded_oligos) For the two self-complementary primers only bracketing is checked, so an arbitrarily wide range would pass there.; (was tm_range_is_exactly_the_lowest_and_highest_tm_of_the_expanded_oligos) Overlaps tests/test_core_sequence.py::test_melting_temperature_matches_santalucia_and_biopython, which checks bracketing for 'AGCGGRTAACAATYTCACACNGGA' only; this Rust test is the stronger one (exact bounds, 4 primers).

## Every test, by file

### `tests/test_assembly.py` (8 tests)

| Test | Why it matters | Checked against | Strength | Review notes |
|---|---|---|---|---|
| `test_fragment_and_gibson_products` | Without it, Gibson designs could hand the user primers or a predicted product whose junctions do not actually exist in the assembled DNA, and silently accept ends that share no homology. | hand-computed | strong | Long test covering four separate behaviours; a failure in the last block does not say which engine is wrong without reading the loop.<br>Tm threshold 52 is hard-coded here and also in the primer tests. |
| `test_golden_gate_products_and_overhang_checks` | Without it, a Golden Gate design could produce a circle with the wrong overhangs, accept overhangs that cross-ligate or repeat, or fail to warn about leftover enzyme sites, so the cloning would not work at the bench. | hand-computed | strong | Also re-tests the overhangs module directly (OH.assess, pair_risk, junction_checks), which could live in its own test. |
| `test_sidewinder_oligos_toeholds_and_barcodes` | Without it, the ordered oligos could fail to tile the intended product, repeat a toehold so the assembly joins in the wrong place, or reuse barcodes, and the user would only discover it after paying for synthesis. | property/invariant | strong | The opening loop over 60 seeds hides how often a random design is rejected; if the toehold picker regressed to accept almost nothing, the test would still pass on a later seed. |
| `test_products_keep_every_orf_on_its_own_bases` | Without it, an assembled product could annotate a CDS a few bases off, so the user reads a shifted frame and the wrong protein without any error being shown. | independent implementation | strong | Slowest test in the file (many seeds x methods); the Golden Gate loop over 6 seeds adds little once one seed passes. |
| `test_the_product_carries_the_features_and_the_primers_of_its_fragments` | Without it, the simulated product could lose or duplicate the fragments' features and primer_bind annotations, or credit bases to a primer tail that the fragments already shared, misleading the user about what the design actually adds. | hand-computed | strong | Very long; mixes feature transfer, primer annotation, preview colouring and junction wording in one function, so a failure needs the line number to interpret.<br>One assertion checks a full English sentence verbatim, which will break on any rewording. |
| `test_restriction_ligation_products` | Without it, the predicted ligation product could be the wrong arc of the vector or have the insert at the wrong place, so the user would clone a construct different from the one shown. | hand-computed | strong | Depends on two demo GenBank files; if those are edited the expected arcs change silently. |
| `test_restriction_ligation_end_rules_and_warnings` | Without it, the program could pair ends that cannot actually ligate (or refuse ones that can) and stay silent about an insert that can go in backwards or a vector that will re-close on itself. | hand-computed | strong |  |
| `test_assembly_reports_bad_input_instead_of_failing` | Without it, a zero-length or very short fragment could crash the assembly panel with a traceback instead of showing the user a red check saying what is wrong. | property/invariant | adequate | Mostly asserts 'not ok and some bad check', so a wrong but non-crashing message would pass.<br>The repeated-overhang case overlaps with the ('unique','bad') case in test_golden_gate_products_and_overhang_checks. |

### `tests/test_checkup.py` (6 tests)

| Test | Why it matters | Checked against | Strength | Review notes |
|---|---|---|---|---|
| `test_bacterial_cds_verdicts` | Without it, the Checkup could tell a user that a gene will be expressed in E. coli when a terminator blocks it, when it is on the wrong strand or when it has no RBS - or refuse to predict the parts a plain file does not annotate. | hand-computed | strong |  |
| `test_eukaryotic_hosts_and_unit_ends` | Without it, a mammalian expression construct could be judged by bacterial rules - no Kozak check, a polyA or a Pol III promoter ignored - and a shuttle vector would not get a per-host verdict. | hand-computed | strong | The CMV and U6 promoters are random filler named 'CMV promoter' / 'U6 promoter', so this exercises name-based recognition only, not any sequence motif. |
| `test_how_far_a_unit_runs_and_how_it_ends` | Without it, the drawn transcript could stop at the wrong base, miss that it wraps the origin, claim an exact +1 it cannot know, or describe its ending wrongly, misleading the user about what is actually transcribed. | hand-computed | strong |  |
| `test_origins_and_prediction_thresholds` | Without it, the Checkup could name the wrong origin family (or miss that a plasmid has none, or needs a pir host), and the motif predictors could flood a real plasmid with chance hits. | property/invariant | adequate | The <= 8 bounds on a single random sequence are a loose regression threshold, not a real false-positive rate; a predictor that returned zero hits would also pass.<br>Origin recognition here is driven by the feature labels, not the sequence. |
| `test_library_parts_on_a_real_plasmid` | Without it, the built-in parts library could override a file's own annotations, or fail to supply the terminator a bare sequence needs, changing the verdicts on real plasmids. | smoke | adequate | Only three assertions over one file; it shows the library is consulted but not that any particular part is matched correctly.<br>Tied to a demo file whose annotation changes would change the result. |
| `test_unit_features_keep_their_geometry_and_survive_a_file_round_trip` | Without it, a transcription unit saved into a file could be drawn with the wrong +1, the wrong thick/thin split or lost altogether on reload, and old files written under the program's former name would stop working. | round-trip/self-consistency | strong |  |

### `tests/test_codons.py` (4 tests)

| Test | Why it matters | Checked against | Strength | Review notes |
|---|---|---|---|---|
| `test_recoding_preserves_the_protein_and_obeys_the_constraints` | Without it, recoding a CDS could silently change the protein the user is cloning, or leave a forbidden restriction site, a bad GC window or a hidden RBS in the synthesised gene. | property/invariant | strong |  |
| `test_codon_choice_follows_the_table_and_the_edit_is_one_undo_step` | Without it, recoding could ignore the chosen usage table (defeating its purpose), be irreproducible between runs, or leave the document in a state undo cannot restore. | in-test reference model | strong | The frequency-matched tolerance (+-0.04 on 3000 codons) is loose enough that a moderately skewed sampler would pass. |
| `test_comparing_a_cds_under_two_genetic_codes` | Without it, the comparison that tells a user a gene will be truncated in another host (TGA read as stop instead of Trp) could point at the wrong codon or miscount the final stop. | hand-computed | strong | The HTML assertion is a weak substring check ('14' could match anything). |
| `test_usage_tables_their_files_and_what_they_say_about_a_gene` | Without it, a user's imported or self-made codon table or genetic code could be lost, mis-parsed, or silently saved nowhere, and the adaptation index and checkup notes computed from it could be wrong. | round-trip/self-consistency | strong | Far too large: CSV parsing, user-file storage, write failures, CAI maths, checkup notes and genetic-code management are six independent features in one test; a single failure stops the rest from running.<br>Mutates global state (active table, tool settings) and relies on try/finally to clean up. |

### `tests/test_core_alignment.py` (4 tests)

| Test | Why it matters | Checked against | Strength | Review notes |
|---|---|---|---|---|
| `test_align_single_queries_coordinates_indels_and_overhangs` | Without it, Local Align could report wrong coordinates, strands, indel positions or invented alignments when a user maps a read or construct onto a sequence. | hand-computed | strong | Leftover 'for _ in range(100): ... break' loop around pick_unslidable does nothing.<br>The linear over-origin check only asserts the hit is 100 bases and ends <= 1000, not which half. |
| `test_align_finds_every_piece_of_a_structural_variant_query` | Without it, Local Align could merge, drop or misplace pieces of a rearranged construct (insertion, inversion, shuffled or duplicated segments), which is what users rely on to diagnose structural variants. | hand-computed | strong | Overlaps somewhat with the 'ends match, middle does not' section of test_align_single_queries_coordinates_indels_and_overhangs. |
| `test_local_align_query_cleaning_and_junction_texts` | Without it, pasting a FASTA, GenBank ORIGIN or RNA query could be rejected or mangled, and the junction descriptions (insertion, deletion, inversion) in the Local Align report could be wrong. | hand-computed | adequate | Junction texts are checked with startswith only, so the rest of each message (sizes, positions) is not verified.<br>Duplication wording is mentioned in the docstring but no duplication junction is tested. |
| `test_folder_search_finds_a_query_in_its_file_record_place_and_strand` | Without it, Folder Search could fail to find a sequence in the user's files, report the wrong file, record, position or strand, or crash on an unreadable file. | hand-computed | adequate | The cancellation check is weak: the query ACGTACGTACGT would almost certainly not be found in that folder anyway, so a search that ignores cancelled() would still pass.<br>'files >= 5' depends on the contents of demo/. |

### `tests/test_core_formats.py` (4 tests)

| Test | Why it matters | Checked against | Strength | Review notes |
|---|---|---|---|---|
| `test_bam_reads_pileup_coverage_and_variants_match_a_model` | Without it, the read alignment track could show reads at wrong places, miscount depth, ignore duplicate/QC filters or call wrong SNPs and indels from a user's BAM file. | in-test reference model | strong | Very large test covering four functions (reads, pileup, coverage, variants); a failure early hides later sections. Could be split. |
| `test_tracks_gff_bed_parsing_query_density_and_search` | Without it, imported annotation tracks could show features with off-by-one coordinates, wrong names or strands, crash on malformed lines, or miss features when zooming and searching. | in-test reference model | strong |  |
| `test_feature_locations_round_trip` | Without it, opening and saving a GenBank file could silently change or drop complex feature locations (joins, complement, fuzzy ends, between-base sites, remote references). | round-trip/self-consistency | strong |  |
| `test_genome_fasta_fetch_and_search_match_slicing_and_brute_force` | Without it, a reference genome opened from a FASTA file could return wrong bases at line breaks or ends, and motif/sequence search on it could miss hits or report wrong strands. | in-test reference model | strong | Overlaps a little with the FASTA genome part of test_primer_binding.py::test_binding_agrees_with_the_documented_rules_written_in_python and the file-scope part of test_search.py. |

### `tests/test_core_sequence.py` (6 tests)

| Test | Why it matters | Checked against | Strength | Review notes |
|---|---|---|---|---|
| `test_melting_temperature_matches_santalucia_and_biopython` | Without it, the primer melting temperatures shown everywhere in the editor (primer list, PCR design) could drift from the published SantaLucia/Primer3 values and nobody would notice. | independent implementation | strong | The abs=0.5 comparison to the hand formula is made redundant by the abs=1e-6 one on the next line. |
| `test_hairpins_and_dimers` | Without it, the hairpin and primer-dimer warnings in primer checks could report nonsense (positive energies, asymmetric pairs, wrong 3'-end flag) without being caught. | property/invariant | adequate | No absolute dG value is checked against Primer3 or a hand calculation: a scaled or shifted energy model would still pass. |
| `test_gc_bins_match_an_exact_count` | Without it, the GC-content track could show wrong percentages, especially zoomed out where cached 1024-base blocks are used, or misbehave at sequence ends. | in-test reference model | strong | The block-based path is only checked to 0.02 tolerance, so a small rounding bug at block edges would pass. |
| `test_orfs_match_a_brute_force_scan` | Without it, the ORF track could miss open reading frames, report nested duplicates, or place reverse-strand ORFs at wrong coordinates. | in-test reference model | strong | The ORF across the origin is only checked for count and strand, not for its start/end coordinates.<br>Partly overlaps the ORF part of test_genetic_code.py::test_translation_tables_match_biopython (that one covers alternative tables and starts). |
| `test_restriction_cuts_and_digests_match_biopython` | Without it, enzyme cut sites, overhangs and digest fragment sizes shown on the map could be wrong (strand, Type IIS offsets, origin-spanning sites, degenerate sites) and cloning plans built on them would fail. | independent implementation | strong | The real-plasmid part silently does nothing if the sample file is missing (no skip reported).<br>Picks EcoRI/BsaI by index COMMON[14]/COMMON[9]; guarded by a name assertion, but will fail if the enzyme list is reordered. |
| `test_find_parts_identity_coding_and_repeated_parts` | Without it, automatic annotation of known parts could miss mutated or recoded copies, ignore identity thresholds, or report duplicated/overlapping hits. | hand-computed | strong |  |

### `tests/test_documents.py` (8 tests)

| Test | Why it matters | Checked against | Strength | Review notes |
|---|---|---|---|---|
| `test_saving_never_damages_the_original` | Without it, a save that fails (disk path missing, file locked, bad data) could truncate or corrupt the user's only copy of a file or wrongly mark unsaved work as saved. | property/invariant | strong | The rename-wait check only proves the wait on Windows; elsewhere it only checks the result.<br>Timing-based assertions (>= 0.4 s, < 2 s) could be flaky on a loaded machine. |
| `test_corrupt_and_truncated_files_are_a_clear_error` | Without it, opening a truncated or damaged file could crash the app or, worse, open a partial sequence that the user then saves over the original. | property/invariant | strong | 'errors > 100' is a loose count; it does not say which cut points must fail.<br>Parts depend on demo files being present. |
| `test_random_edits_keep_bases_and_features_then_undo_redo_exactly` | Without it, some combination of edits (especially across the origin or with joins) could shift features onto the wrong bases or make undo/redo restore a different document, corrupting annotations silently. | in-test reference model | strong | The models and assertions live in app/selftest.py, not in the test, so the report reader cannot see them here; the test is a one-line wrapper.<br>The app's own self-check runs shorter versions (seeds 1 and 2, 150 edits) of the same machinery: overlapping by design. |
| `test_edited_base_ranges_follow_edits_undo_redo_and_saves` | Without it, the highlight showing which bases changed since the last save could point at the wrong bases after further edits, undo or rotation, misleading the user about what they changed. | hand-computed | strong | Across-origin case only checks the total length of ranges, not their positions. |
| `test_undo_history_steps` | Without it, the History panel could show wrong step names, primer database changes could fall out of sync with base edits on undo/redo, the DNA lock could be bypassed, or undo across records of a multi-record file could go in the wrong order. | hand-computed | strong | Five fairly independent features (annotation step names, primer undo, primer conflicts, DNA lock, FileHistory) in one test; the DNA lock and FileHistory parts could be separate tests. |
| `test_copy_reverse_complement_keeps_annotations_flipped` | Without it, Copy Reverse Complement could paste the right bases but with features on the wrong strand or wrong positions, or break across the plasmid origin. | hand-computed | adequate | Only 'top' gets an exact location check; 'bottom' and the join's parts are checked by strand or prefix only.<br>Across-origin case checks the bases but not that features are carried across the origin.<br>The RC helper is the app's own reverse_complement, so a bug shared by both would not be caught (it is simple enough that this is minor). |
| `test_text_exports_and_the_project_file` | Without it, GFF3/BED/feature-table/protein exports could be malformed (e.g. wrong strand convention for NCBI submission) and Open DNA Project files could lose tracks, display settings or the DNA lock. | hand-computed | adequate | GFF3, BED and proteins_fasta checks on the demo vector are shallow (line counts, prefix, one strand column); a wrong coordinate would pass.<br>Project round trip compares only sequence and feature count, not feature contents.<br>Exports and the project file are two topics merged into one test. |
| `test_gene_groups_and_display_switches` | Without it, gene + CDS pairs might not group by locus_tag, toggling display options could re-read all features from the core (slow on genomes) or wipe Local Align hits, and type counts could go stale after edits. | hand-computed | adequate | Docstring says the switch 'also dropped the Local Align hits' while the test asserts they are kept; it refers to an old bug and reads backwards.<br>Group T_F expects 2 parts although gene F has a gene, a CDS and a misc_feature; the intent (the gene feature is the group head, not a part) is not explained in the test.<br>Does not check collapse_genes actually changes what is shown. |

### `tests/test_enzymes_and_guides.py` (5 tests)

| Test | Why it matters | Checked against | Strength | Review notes |
|---|---|---|---|---|
| `test_cut_finder_the_common_list_and_the_filter` | Without it, the program could report cut positions at the wrong base (or miss sites on the bottom strand, across the origin, or with IUPAC codes), so every digest, map and cloning plan built on them would be wrong. | hand-computed | strong |  |
| `test_what_a_cut_says_about_itself` | Without it, the tooltip a user reads before ordering a digest could state the wrong cut base, the wrong overhang or the wrong strand, which is the only place those numbers are explained. | hand-computed | strong | Asserts on exact English sentences, so harmless rewording of the tooltip breaks the test. |
| `test_guide_candidates_cuts_and_scores` | Without it, the Spacer Finder could list guides with the wrong PAM, the wrong strand or a cut site in the wrong place, so a user would order a guide that cuts somewhere other than shown. | hand-computed | strong |  |
| `test_other_guide_sites_here_and_in_a_genome` | Without it, the off-target counts and the good/warn/bad verdict a user relies on when choosing a guide could be wrong, or the fast Rust search could silently miss sites a careful search would find. | independent implementation | strong |  |
| `test_a_saved_guide_site_is_read_back_and_round_trips_in_the_file` | Without it, a saved guide could come back from a GenBank file with the protospacer, PAM or cut drawn in the wrong place, be unreadable in other software, or break for files written under the program's former name. | round-trip/self-consistency | strong | Six separate concerns (layout, tooltip, parsing, file bytes, legacy key, exports) in one function. |

### `tests/test_feedback.py` (3 tests)

| Test | Why it matters | Checked against | Strength | Review notes |
|---|---|---|---|---|
| `test_scrub_removes_private_text` | Without it a crash report or log could leak the user's name, folder paths, file names or their DNA/protein sequences into a public GitHub issue. | hand-computed | adequate | Mostly 'not in' checks: a scrub that removed too much would pass most rows (only the plain-sentence row, the traceback row and the placeholder counts guard against over-scrubbing).<br>The /home/bobsmith input is scrubbed and checked twice in a row with overlapping assertions (lines 90-91). |
| `test_report_links_are_encoded_capped_and_carry_nothing_private` | Without it the 'Report a bug' / 'Request a feature' links could be malformed, too long for the browser, lose the end of the log, or carry private data to GitHub. | hand-computed | strong | The privacy section partly re-tests scrub (previous test) through the link builder; reasonable as an end-to-end guard. |
| `test_log_file_and_error_hooks` | Without it unhandled errors could go unlogged or unreported, the same error could pop a dialog repeatedly, the log could grow without bound or contain private data, or the error dialog itself could crash the app. | hand-computed | strong | Very long test bundling about ten independent features; splitting it would make a failure point at one feature.<br>The assertion that 'needs+a+report' is not in the URL is an indirect check of the encoding style; its intent is unclear. |

### `tests/test_file_assoc.py` (3 tests)

| Test | Why it matters | Checked against | Strength | Review notes |
|---|---|---|---|---|
| `test_register_writes_only_our_keys_and_status_reads_them_back` | Without it double-clicking a .gb/.fa file might not open openDNA, or registering could overwrite other programs' entries or the user's default-app choice in the registry. | in-test reference model | strong | The fake registry is the test's own model of winreg; the real WinRegistry wrapper is not exercised here, so a bug in that thin layer would go unnoticed. |
| `test_unregister_removes_exactly_what_we_added` | Without it removing the file associations could delete other programs' registry keys or the user's settings, or leave openDNA entries behind. | round-trip/self-consistency | strong |  |
| `test_startup_question_decision_table` | Without it the app could nag users at every start about file associations, ask when running from source or off Windows, ignore 'don't ask again', or fail to repair the association after the exe was moved. | hand-computed | strong |  |

### `tests/test_genes_and_history.py` (2 tests)

| Test | Why it matters | Checked against | Strength | Review notes |
|---|---|---|---|---|
| `test_gene_groups_fold_their_parts_and_survive_edits` | Without it, the feature list and map could fold the wrong features under a gene (or fail to fold them), and after editing the sequence a feature could end up covering bases it never covered. | hand-computed | strong |  |
| `test_history_texts` | Without it, the History panel - the user's only record of what an edit actually did - could describe an edit wrongly or lump several edits into one line, making a mistake impossible to trace back. | regression value | adequate | Entirely string-literal assertions on UI wording: it pins the current text but a wrong-but-plausible message (e.g. an off-by-one position) is only caught where the position is spelled out.<br>One assertion builds its expected text from the history entry it is checking, so that part cannot fail. |

### `tests/test_genetic_code.py` (2 tests)

| Test | Why it matters | Checked against | Strength | Review notes |
|---|---|---|---|---|
| `test_translation_tables_match_biopython` | Without it, translations, protein search and ORFs for non-standard organisms (mitochondria, mycoplasma, bacteria with alternative starts) could use wrong codon assignments. | independent implementation | strong | The line asserting str(Seq(codon).translate()) == aa only checks Biopython against itself.<br>The ORF table/start section partly overlaps test_core_sequence.py::test_orfs_match_a_brute_force_scan (that one covers table 1 only). |
| `test_cds_translation_and_translation_refresh` | Without it, CDS amino acids shown in the views and exported, and the /translation qualifiers rewritten by the refresh command, could ignore /transl_table, /codon_start or /transl_except, or overwrite pseudo/partial genes. | hand-computed | strong | Changes shared tool settings during the test (restored by the fixture). |

### `tests/test_input_parsing.py` (1 test)

| Test | Why it matters | Checked against | Strength | Review notes |
|---|---|---|---|---|
| `test_typed_or_imported_text_is_parsed_or_refused_with_a_reason` | Without it the range box could select the wrong bases, pasted protein text could keep junk characters, and bad names, NCBI keys or primer lists could end up stored in the primer database. | hand-computed | strong | Bundles five unrelated input areas (range box, protein paste, name/key validation, primer lists, primer database) in one test; a failure in the range box hides all the primer checks. Splitting into 2-3 tests would help.<br>Writes good1, good2 and fine_one into the primer database of the module's temp data dir; harmless, but state is shared with anything else in the same process using that dir. |

### `tests/test_multi_record.py` (2 tests)

| Test | Why it matters | Checked against | Strength | Review notes |
|---|---|---|---|---|
| `test_one_record_is_read_and_saved_in_place` | Without it, editing one chromosome or plasmid inside a multi-record GenBank file could corrupt the other records or leave the file half-written when the save fails. | property/invariant | strong | Skips entirely if the demo plasmids are missing.<br>match_contig and other_records_note are main-window helpers unrelated to in-place saving; they would fit better in a tracks or FASTA-save test. |
| `test_records_are_added_reordered_and_removed_byte_for_byte` | Without it, reordering, removing or adding records in a GenBank or FASTA file could alter untouched records, break CRLF/BOM files, reuse a stale index, or create duplicate record names. | property/invariant | strong | The .dna-with-primers part only checks that some primer_bind exists, not names or positions.<br>The test builds the primer_bind list itself, so it does not test how the app's Add Record dialog builds it. |

### `tests/test_online.py` (4 tests)

| Test | Why it matters | Checked against | Strength | Review notes |
|---|---|---|---|---|
| `test_ncbi_client_search_cache_rate_limit_and_errors` | Without it fetching sequences from NCBI could break, send wrong parameters, ignore NCBI's request rate limit (getting the user blocked), cache error pages as sequence files or fail on a temporary 429. | hand-computed | adequate | The rate-limit check uses real wall-clock time (adds about 0.6 s; the '< 0.6 s with a key' branch could be flaky on a very slow machine).<br>The 429 retry checks a single retry only; the retry limit and giving up are not checked. |
| `test_structure_sources_cache_and_auto_order` | Without it the protein 3D view could pick the wrong source, accept invalid input, ignore its cache, hide why every source failed, or break the user's own local folding command. | hand-computed | adequate | Checks on the returned models are mostly presence checks (plddt truthy, 'ATOM' in text, id in detail) rather than content; the ESMFold 0-1 to 0-100 scaling is the one exact value. |
| `test_structure_page_script_safety_and_alphafold_reply_not_followed` | Without it a crafted protein name, PDB file or AlphaFold reply could inject script into the 3D viewer page or make the app fetch a local or foreign file. | hand-computed | strong | The accession_of check is unrelated to page safety; it would fit better with the structure-source test. |
| `test_pdb_entries_and_the_model_read_from_a_pdb_text` | Without it the PDB search list could show duplicates or wrong identity/coverage, linked entries could be badly ordered, and the 3D view could mis-read chains, ligands or secondary structure, or map protein residues to the wrong atoms. | hand-computed | strong | Two unrelated topics (RCSB search/linking and PDB text parsing) in one test; could be split. |

### `tests/test_parts_library.py` (2 tests)

| Test | Why it matters | Checked against | Strength | Review notes |
|---|---|---|---|---|
| `test_collections_sources_and_switches` | Without it, the parts library could show parts from a source the user switched off, lose source/measurement metadata, let users delete shipped parts, or crash on one bad row. | hand-computed | strong | The first two sections rewrite the same two files without clearing; harmless but redundant.<br>Five subtopics in one test (the docstring calls it a merged test). |
| `test_shipped_parts_meet_the_data_quality_rules` | Without it, the shipped parts library could gain junk entries (primer sites, cloning-site slices, unnamed types), parts too short to annotate meaningfully, duplicates, or parts without a citation, cluttering every user's automatic annotation. | property/invariant | adequate | A data-lint check on shipped files rather than a code test; it fails when data changes, not when code breaks. Still useful as a guard.<br>The rules (NEVER, FLOORS) come from scripts/build_parts.py, so a loosened rule there would pass silently.<br>Duplicate check is exact sequence only; reverse-complement duplicates or near-duplicates are not caught. |

### `tests/test_primer_binding.py` (5 tests)

| Test | Why it matters | Checked against | Strength | Review notes |
|---|---|---|---|---|
| `test_binding_agrees_with_the_documented_rules_written_in_python` | Without it, the primer binding sites shown on the map could deviate from the documented seed/extension rules on any of thousands of combinations of mismatches, tails, strands and topologies. | in-test reference model | strong |  |
| `test_binding_mismatches_seed_tails_and_extension` | Without it, specific binding rules (where a mismatch is allowed, how far the annealed part extends, how tails are handled) could change and primers would be shown binding where they do not, or not where they do. | hand-computed | strong | Partly redundant with the random comparison in test_binding_agrees_with_the_documented_rules_written_in_python, but it pins each rule by name, which makes failures easier to read. |
| `test_binding_origin_ends_repeats_and_letters` | Without it, primers on circular plasmids across the origin, repeated or palindromic sites, lower-case/RNA/IUPAC primers or N in the template could produce missing, duplicated or wrong sites. | hand-computed | strong |  |
| `test_binding_matches_the_primer_sites_stored_in_the_real_dna_file` | Without it, the primer sites the app computes could disagree with those SnapGene stored in a real lab file, which is the most direct check that users see the same sites as in their other tool. | independent implementation | strong | Skips if the sample file is missing; depends on one specific private sample file. |
| `test_editor_primer_sites_pcr_product_specificity_and_save` | Without it, the editor could colour the wrong mismatched bases, build a wrong PCR product (tails, origin, mismatches), give wrong non-specific priming warnings, or alter the document and undo history when saving primer sites. | hand-computed | adequate | Bundles four unrelated features (mismatch display, PCR product, non-specific priming check, save) into one test; would read better split.<br>Leftover comment '# the dialog (built, not run) says the same size...' with no code under it.<br>Uses a module-level RNG shared across tests, so the random data depends on which tests ran before. |

### `tests/test_primer_design.py` (4 tests)

| Test | Why it matters | Checked against | Strength | Review notes |
|---|---|---|---|---|
| `test_primer_check_levels_follow_their_limits` | Without it, the green/orange/red verdicts on a primer could fire at the wrong value or ignore the user's own limits, so primers that will not work would look fine (or vice versa) before they are ordered. | hand-computed | strong | Tm is stubbed to a fixed value, so the actual melting-temperature calculation is not exercised here. |
| `test_where_a_primer_comes_from_decides_how_it_is_kept` | Without it, a primer from the shared database could be confused with one stored in the file, be drawn in the wrong colour or hidden by the wrong filter, or an imported primer could fail to be saved or to undo cleanly. | round-trip/self-consistency | strong |  |
| `test_primers_are_designed_for_a_target_at_or_across_the_origin_of_a_circle` | Without it, designing primers for a target at or across a plasmid's origin could return nothing, or primers whose bases are not actually on the template or whose product does not contain the target. | property/invariant | adequate | Only checks that pairs exist and are self-consistent; nothing here checks that the chosen primers are good ones (Tm, specificity), so a designer returning poor pairs would pass. |
| `test_the_files_written_for_ordering` | Without it, the plate file sent to a synthesis company could put oligos in the wrong wells, overflow past well 96, drop the scale or purification, or corrupt names and sequences - a real-money error. | hand-computed | strong |  |

### `tests/test_project.py` (2 tests)

| Test | Why it matters | Checked against | Strength | Review notes |
|---|---|---|---|---|
| `test_release_metadata_and_exe_recipe_are_consistent` | Without it, a release could ship with mismatched version numbers (breaking the updater) or a Windows exe missing data files the app needs (parts library, tool settings, 3D viewer), which only shows up for users of the frozen build. | property/invariant | adequate | Not a test of DNA editing; it is a release-hygiene check. Relevant for packaging but could be marked as such.<br>The bundling check is text matching on the spec ('"name"' appears anywhere), so a name mentioned but not actually bundled would pass.<br>Only catches data files reached through the exact with_name("...") pattern. |
| `test_old_data_folder_is_copied_once` | Without it, users upgrading from the program's old name could lose their primers and settings, or a later start could overwrite their new data with the old copy. | hand-computed | adequate | Lives in test_project.py (release metadata) though it is about primer_db's data migration.<br>The copy-error path (half copy returns False) is not tested. |

### `tests/test_protein.py` (1 test)

| Test | Why it matters | Checked against | Strength | Review notes |
|---|---|---|---|---|
| `test_mass_pi_charge_and_numbered_layout` | Without it the protein panel could show a wrong molecular weight, pI or charge, or clicking in the numbered protein text could select the wrong residue. | hand-computed | weak | pI and charge are checked only for sign/range on homopolymers: a wrong pKa table or a broken pI search could still pass. An exact pI and mass for a known protein (e.g. against Biopython ProtParam, which is installed in .venv) would make it strong.<br>Mass is checked for one dipeptide only; a wrong residue mass for any amino acid other than G would go unnoticed.<br>The index/offset check is a round trip only; offset_of and index_at that were wrong in the same way would pass. |

### `tests/test_reference.py` (2 tests)

| Test | Why it matters | Checked against | Strength | Review notes |
|---|---|---|---|---|
| `test_off_target_search_covers_every_record` | Without it, the primer and guide off-target checks could miss copies on the second or later chromosome, report the open sequence as its own off-target, show green after a search that timed out, or lose/double hits at window edges. | hand-computed | strong | The time-out section injects a result into the private _cache dict, coupling the test to implementation details.<br>Random 'dna' filler uses a module-level RNG, so sequence content depends on the order tests run in that module (a chance extra match is unlikely but possible). |
| `test_ncbi_genome_search_and_download` | Without it, searching NCBI for a reference genome could build wrong URLs, and downloading one could produce a broken FASTA (records joined without newlines), leave partial files, or download the same genome repeatedly. | hand-computed | adequate | Network is faked at ncbi.fetch/open_url, so a change in NCBI's real API or response format would not be caught.<br>Unrelated to the multi-record off-target search that names this file; could live with other NCBI tests. |

### `tests/test_repeats_msa.py` (2 tests)

| Test | Why it matters | Checked against | Strength | Review notes |
|---|---|---|---|---|
| `test_repeats_are_found_where_planted_and_not_in_random_dna` | Without it, the repeat finder could miss planted direct, inverted or tandem repeats, flood random DNA with chance hits, or hang on long low-complexity input. | hand-computed | adequate | Positions are checked within 2 bases and lengths within a range, so small off-by-one errors pass.<br>pair_dots is checked only for being non-empty.<br>In the two-sequence check `near` uses zip without a length check. |
| `test_multiple_alignment_consensus_identity_and_parsing` | Without it, the built-in multiple aligner could lose bases, misplace a gap, or give wrong consensus, identity and column/position mapping in the alignment view. | hand-computed | adequate | identity(1) for sub120 (row 2) is not checked.<br>Parsing accepts rows of different widths (5 and 6) without comment; the test does not say whether that should be refused. |

### `tests/test_roundtrip.py` (6 tests)

| Test | Why it matters | Checked against | Strength | Review notes |
|---|---|---|---|---|
| `test_real_files_round_trip_and_agree_with_biopython` | Without it, opening and saving a real plasmid or genome could silently drop or alter features, qualifiers or bases, or write a file other tools read differently. | independent implementation | strong | Skips entirely if the demo files are missing; the .dna sections return early without them.<br>The source id 'genome_dna' actually refers to the GenBank file BW25113_2026.gb (GENOME), not a .dna file: misleading name.<br>Saving and reloading the 4.6 Mb genome three times makes this one of the slower tests.<br>The is_directional checks are unrelated to round trips and could live in a feature-direction test. |
| `test_unusual_records_round_trip` | Without it, unusual but legal locations (joins, fuzzy ends, order(), across the origin, between sites), tiny or empty sequences, or lower-case bases could be corrupted on save without anyone noticing. | independent implementation | strong | Very long test mixing many unrelated edge cases; a single failure hides the rest. Could be split into location forms, topology, tiny sequences and case.<br>The empty-sequence-after-delete case only checks the file exists, not that it reloads. |
| `test_qualifier_text_survives_and_is_written_as_documented` | Without it, notes with quotes, unicode or long text could be mangled, wrapped badly or gain spaces when saved, so other tools (or a reload) would show different annotation text. | independent implementation | strong |  |
| `test_primers_survive_save` | Without it, primers attached to a plasmid (from a .dna file or imported from the database) could vanish or lose their full sequence when the file is saved as GenBank. | round-trip/self-consistency | adequate | The .dna section is skipped if the demo .dna file is missing.<br>Overlaps with test_multi_record's .dna-primers-as-primer_bind check and with primer tests elsewhere. |
| `test_fasta_files_are_read_like_biopython` | Without it, opening a multi-sequence FASTA (genomes, Windows line endings, BOM, lower case) could return wrong bases at line breaks or include stray characters, corrupting what the user sees and edits. | independent implementation | strong |  |
| `test_fasta_writing_is_the_string_edit_of_the_file` | Without it, saving an edited chromosome back into a multi-FASTA could damage the other records, lose CRLF or a BOM, leave a stale index, or write the wrong edited bases. | in-test reference model | strong | Partly overlaps test_multi_record::test_records_are_added_reordered_and_removed_byte_for_byte (CRLF/BOM byte preservation), but that one tests reorder/append rather than edit-in-place. |

### `tests/test_search.py` (1 test)

| Test | Why it matters | Checked against | Strength | Review notes |
|---|---|---|---|---|
| `test_sequence_search_engine` | Without it, the Find box could miss DNA matches on the other strand or with IUPAC/mismatches, fail across the origin, ignore its options for names, primers, enzymes and translations, or search outside the selection. | hand-computed | adequate | The reverse-strand protein search with 'AXW' computes hits but asserts nothing.<br>IUPAC and fuzzy checks use 'in'/'any', so extra spurious hits would not be noticed.<br>One long test mixing engine, editor bookkeeping and file search; a failure early hides the rest. |

### `tests/test_selftest.py` (1 test)

| Test | Why it matters | Checked against | Strength | Review notes |
|---|---|---|---|---|
| `test_every_selftest_check_passes` | Without it the Help > Troubleshooting > Self-Test checks could start failing in a release without anyone noticing before a user runs them. | property/invariant | adequate | The real checks live in app/selftest.py; this test's strength is theirs and cannot be judged from the test file.<br>Overlap: tests/test_documents.py::test_random_edits_keep_bases_and_features_then_undo_redo_exactly runs the same selftest.run_random_edits machinery with 500 edits; this test re-runs the short in-app versions. Still worth keeping as the guard that the in-app self-test itself passes. |

### `tests/test_updater.py` (5 tests)

| Test | Why it matters | Checked against | Strength | Review notes |
|---|---|---|---|---|
| `test_release_check_offers_only_a_newer_complete_release` | Without it the app could offer users an older, draft, prerelease or incomplete release, or one whose files are hosted off the project's GitHub repository, as an update. | hand-computed | strong | Large multi-section test; a failure in an early section hides the later ones.<br>Error texts are checked by substring ('Refused', 'understood', 'failed'), not exact message. |
| `test_download_comes_only_from_github_and_is_kept_only_when_it_matches_its_checksum` | Without it the updater could download a program from a look-alike or hijacked host, follow a redirect off GitHub, or install a tampered or partial file. | hand-computed | strong | The real network path is not exercised (correct for a unit test); the redirect check calls the handler method directly rather than going through a real redirect. |
| `test_install_swaps_relaunches_and_cleans_up` | Without it a successful update could fail to replace the running exe, fail to restart the program, start the new program with a broken PyInstaller environment, or leave .old/.part files behind (or delete the wrong file). | hand-computed | strong |  |
| `test_install_refuses_or_rolls_back_leaving_the_old_program` | Without it a failed update (no write access, antivirus lock, broken new exe, wrong build type or platform) could leave the user with no working program. | hand-computed | strong |  |
| `test_automatic_check_obeys_the_switches_runs_once_a_day_and_remembers_a_skip` | Without it the app could check for updates when the user turned it off, check on every start, stop checking after the clock goes back, or forget a version the user chose to skip. | hand-computed | adequate | Only checks that a skipped version is stored, not that the automatic check then stays quiet about it (that logic, if in main.py, is untested here). |

### `tests/test_variants.py` (1 test)

| Test | Why it matters | Checked against | Strength | Review notes |
|---|---|---|---|---|
| `test_variants_are_written_into_the_bases_or_annotated_beside_them` | Without it, applying variants found in sequencing reads could corrupt the user's sequence - wrong base, shifted indel, several edits interfering - and the damage would be hard to undo. | in-test reference model | strong | One long function covering apply, ordering, locking and annotation; splitting it would make failures easier to read. |

### `core/src/align.rs` (12 tests)

| Test | Why it matters | Checked against | Strength | Review notes |
|---|---|---|---|---|
| `exact_hit_forward_and_reverse` | Quick Align (local alignment) could report a perfectly matching query at the wrong place, strand or identity and nobody would notice. | hand-computed | strong | Forward case (line 502-505) does not assert hits.len()==1, so spurious extra hits would pass; the reverse case does (line 509). |
| `mismatch_insertion_and_deletion` | A query carrying a SNP and small indels could be split, truncated or reported with the wrong mismatch, so sequencing-read or clone checks against a plasmid would show false differences. | hand-computed | adequate | Positions of the insertion (template 1150) and deletion (1250..1255) are never asserted, only their total sizes (lines 525-526); a misplaced indel would pass.<br>Reverse-strand case (line 530) checks only strand and span, not the mismatch/indel reports. |
| `a_mismatch_near_either_end_of_a_long_query_is_reported_not_trimmed_away` | A mismatch in the last few bases of a query (a typical mutation near a read or primer end) could be silently clipped off the hit, so the user would see a 'perfect' shorter match and miss the difference. | hand-computed | strong | Forward strand only; the end-trimming rule on a reverse-complement query is not checked.<br>Does not assert hits.len()==1 (indexes [0] directly); an extra spurious hit would pass. |
| `a_short_query_with_mismatches_near_both_ends_is_still_one_whole_hit` | A short (45 bp, oligo-sized) query with a few mismatches near its ends could come back as a trimmed or broken hit, hiding the mismatches from the user. | hand-computed | adequate | Mismatch positions and qstart/qend are not asserted, only the count 3; nor is hits.len()==1, despite the name 'one whole hit'.<br>The mismatch at 35 is 10 bases from the end, so only one end (base 9) is really 'near' an end. |
| `no_hit_for_random_query` | Quick Align could start reporting chance matches for an unrelated query, so users would trust false hits. | property/invariant | adequate | Single seed pair only; a threshold regression that lets chance hits through on other seeds could pass. |
| `query_hanging_off_the_template_end_aligns_only_the_overlap` | A query running past the end of a linear sequence could be reported with invented bases or wrong query/template coordinates beyond the end. | hand-computed | strong |  |
| `overhanging_query_wraps_on_a_circle_and_a_tiny_overlap_is_no_hit` | A read spanning the origin of a plasmid could be shown as two pieces (or as one on a linear sequence), and a tiny chance overlap at the end could appear as a hit. | hand-computed | strong | Linear case (line 587) checks only the count 2, not the coordinates of the two halves. |
| `every_piece_of_a_patchwork_query_is_reported_including_the_inverted_one` | Quick Align on a rearranged construct (insertion, inversion, reordering) could drop a piece or report the inverted piece on the wrong strand or place. | hand-computed | strong |  |
| `pieces_on_one_diagonal_with_a_long_unmatched_middle_are_not_merged` | Separate matching pieces with a long unrelated stretch between them could be merged into one long low-identity hit, hiding a large insertion or deletion. | hand-computed | strong |  |
| `max_hits_still_limits_the_number_of_pieces_and_keeps_the_best` | The max-hits setting could stop limiting the number of pieces returned for a patchwork query. | property/invariant | weak | Name promises 'keeps the best' but lines 666-667 only check counts; nothing compares align(..,2) with the top 2 of align(..,10), so returning the 2 worst pieces would pass. |
| `a_piece_repeated_in_the_query_is_found_twice_on_the_same_template_place` | A tandem duplication in a query could be collapsed into one hit, so the user would not see the duplication. | hand-computed | strong |  |
| `an_unmatched_stretch_is_cut_out_but_a_short_one_is_still_bridged` | A short block of differences could split a hit into fragments, or a long unrelated stretch could be swallowed into one misleading hit. | hand-computed | strong | Comment line 689 says the 8 changed bases are 'shown as mismatches' but the mismatch list is not asserted. |

### `core/src/alignments.rs` (1 test)

| Test | Why it matters | Checked against | Strength | Review notes |
|---|---|---|---|---|
| `reads_mismatches_and_coverage` | Imported BAM reads could be drawn with wrong positions, strands, mismatches, indels or soft clips, and the coverage, pileup and variant tracks could show wrong numbers. | hand-computed | strong | Writes fixed file names ge_test_reads.bam/.sam into the shared system temp dir and never deletes them; two concurrent cargo test runs could clobber each other.<br>Overlaps tests/test_core_formats.py::test_bam_reads_pileup_coverage_and_variants_match_a_model, which covers the same BAM decoding from Python with a model check on random reads; this is the faster unit-level copy. |

### `core/src/annotate.rs` (1 test)

| Test | Why it matters | Checked against | Strength | Review notes |
|---|---|---|---|---|
| `finds_parts_on_both_strands_across_the_origin_and_as_protein` | 'Detect common features' could miss known parts on the reverse strand, across the origin, or as codon-changed CDSs, or report the wrong identity. | hand-computed | strong | The promoter checks use contains() (lines 228-229), so extra spurious hits in the circular case would not be detected.<br>Protein-level search is only exercised on the forward strand; reverse-strand protein hits (coordinate mapping at line 164 of annotate.rs) are untested here. |

### `core/src/bio.rs` (3 tests)

| Test | Why it matters | Checked against | Strength | Review notes |
|---|---|---|---|---|
| `finds_orfs_on_both_strands` | The ORF finder could show ORFs at wrong coordinates, ignore the minimum length, or miss bottom-strand ORFs. | hand-computed | strong |  |
| `alternative_starts_and_other_tables` | The ORF finder could use alternative start codons when they are switched off or from the wrong genetic code, or read TGA as a stop in the Mycoplasma code. | hand-computed | strong | Table-4 case (line 102) checks only len()==1, not the coordinates 0..12. |
| `finds_an_orf_that_crosses_the_origin_of_a_ring` | An ORF that crosses the origin of a plasmid could be missed or reported twice. | hand-computed | strong | No case for a ring with a start but no stop in a full lap (the 'full lap' guard, bio.rs line 42). |

### `core/src/cas9.rs` (8 tests)

| Test | Why it matters | Checked against | Strength | Review notes |
|---|---|---|---|---|
| `finds_a_site_with_its_pam_on_the_top_strand` | Guide-site search could miss or misplace a perfect protospacer + NGG site on the top strand. | hand-computed | strong |  |
| `finds_a_site_on_the_bottom_strand` | Guide sites on the bottom strand could be missed or drawn at the wrong place. | hand-computed | strong |  |
| `needs_the_pam` | Sites lacking a PAM could be reported as guide sites (false off-targets), or non-NGG PAMs such as SaCas9's NNGRRT could be ignored. | hand-computed | adequate | The NNGRRT case checks only the count, not position/strand. |
| `counts_mismatches_and_records_where` | Off-target sites could be listed with the wrong number or positions of mismatches, misleading guide specificity checks. | hand-computed | strong |  |
| `a_site_can_run_across_the_origin_of_a_circle` | A guide site crossing the origin of a plasmid could be missed. | hand-computed | adequate | End is only checked as > length (line 235); exact end 48 and strand are not asserted. |
| `several_guides_and_lower_case` | Searching several guides at once, or a lower-case sequence, could return mixed-up or empty results. | hand-computed | adequate | Only relative order is checked (line 246), not exact starts 6 and 35; lower-case guides (as opposed to template) are not tested. |
| `finds_a_site_whose_pam_comes_first` | Cas12a (5' TTTV PAM) sites could be searched in Cas9 orientation and so be missed or wrongly found. | hand-computed | strong |  |
| `mismatch_mask_counts_from_the_guides_5_prime_end_next_to_a_first_pam` | For Cas12a the mismatch positions could be counted from the wrong end, so a PAM-proximal mismatch would be weighted as distal in off-target scoring. | hand-computed | adequate | Does not assert there is exactly one site or its position (line 275 indexes hits[0][0] directly). |

### `core/src/code.rs` (7 tests)

| Test | Why it matters | Checked against | Strength | Review notes |
|---|---|---|---|---|
| `a_registered_code_is_found_by_id_and_a_removed_codon_reads_x` | A recoded organism's or user-made genetic code could fail to register, be looked up wrongly, or translate a reassigned codon as its old amino acid instead of X. | hand-computed | adequate | (was a_registered_code_is_found_by_id_and_a_removed_codon_reads_x) Error cases are is_err() only; the reason is not asserted.<br>(was a_registered_code_is_found_by_id_and_a_removed_codon_reads_x) Replacing an existing custom id (documented 'add (or replace)') and the starts string of a registered code are not tested.<br>(was a_registered_code_is_found_by_id_and_a_removed_codon_reads_x) The Python binding register_genetic_code (lib.rs) is not called by any Python test; this is the only test of the code registry.<br>(was every_table_has_64_codons_and_known_ids) Only lengths are checked; the module doc says 'the tests compare every codon' against Biopython: that is done in Python, not here, so these Rust tests alone would not catch a wrong letter in a table.<br>(was every_table_has_64_codons_and_known_ids) Largely superseded by tests/test_genetic_code.py::test_translation_tables_match_biopython, which compares every codon of every table (a wrong length would fail there too); this only adds a fast length/id smoke check. |
| `the_standard_code_and_the_known_differences` | Translation could get basic codons wrong, ignore table-specific differences (TGA=W in table 4, AGA=stop in table 2), or reject lower-case or RNA input. | hand-computed | adequate |  |
| `a_cds_is_translated_with_its_own_table_and_frame` | A CDS could be translated with the app default instead of its own /transl_table or /codon_start, or a broken CDS (internal stop, incomplete codon) could be translated silently. | hand-computed | adequate | Error cases use is_err() only (lines 366, 375): the error kind/message is not checked, so an error for a different reason would pass. |
| `an_alternative_start_codon_is_read_as_methionine` | A CDS starting with an alternative start such as GTG or TTG could show V or L instead of M as its first residue, or internal GTGs could wrongly become M. | hand-computed | strong |  |
| `bottom_strand_cds_and_transl_except` | A complement CDS could be translated from the wrong strand, or a /transl_except (e.g. selenocysteine) could be ignored. | hand-computed | adequate | Error cases are is_err() only; the reason is not asserted. |
| `joined_cds_reads_its_parts_in_order` | A spliced (join) CDS could be translated from its parts in the wrong order. | hand-computed | adequate | Parts are in ascending order; a complement(join) or out-of-order join is not tested. |
| `pseudo_and_partial_cds_are_skipped` | Pseudogenes and 5'-partial CDSs could get a nonsense /translation instead of being skipped with a reason. | hand-computed | strong |  |

### `core/src/dna_file.rs` (2 tests)

| Test | Why it matters | Checked against | Strength | Review notes |
|---|---|---|---|---|
| `reads_synthetic_file` | Opening a SnapGene .dna file could lose or misplace its features, primers, topology or description. | hand-computed | strong | The primer site (0-based '9-12' -> 10..13) is not where 'GGGG' would actually anneal in this synthetic sequence, so the test fixes the 0-based convention by assertion only; it cannot confirm it against real bases (the code comment says real files were checked). |
| `rejects_other_files` | Opening a non-.dna or a damaged .dna file could give an empty document instead of an error. | hand-computed | adequate | Truncated packets (length past end of data) are not tested. |

### `core/src/document.rs` (12 tests)

| Test | Why it matters | Checked against | Strength | Review notes |
|---|---|---|---|---|
| `features_follow_every_kind_of_edit` | Without it, an insert, delete or replace could move, grow, shrink or drop annotations (plain, complement, joined, fuzzy) so they no longer sit on their bases. | hand-computed | strong | (was edit_moves_features_and_undo_restores_exactly) Partly overlaps modified_flag_follows_edits_undo_and_save (modified/undo) and random_edits_undo_redo_follow_a_snapshot_and_string_model (exact undo/redo); its specific value is the hand-checked '2..13' growth of lacZ.<br>(was deleting_exactly_the_feature_or_more_removes_it_and_undo_restores_it_whole) Cases (30,100) and (50,130) delete the whole feature (0-based 50..100) too, yet line 2083 condition a<=50 && b>=100 sends them to the 'removed' branch; the else branch 'must keep what is left' is never executed by any case, so partial deletion is not tested here.<br>(was complement_feature_keeps_its_bases_and_strand_through_edits) bio == rc(covered) (line 2137) is true by construction of the helper for any strand -1 feature, so it adds nothing.<br>(was complement_join_order_survives_edits_and_reverse_complement) The reverse-complement half overlaps reverse_complement_flips_every_feature_onto_the_same_bases (same join and complement(join) shapes); the outside-edit half is unique.<br>(was fuzzy_locations_survive_edits_and_saves_with_their_markers) Fixed temp file name ge_test_fuzzy.gb (line 2279). |
| `edits_are_validated_and_empty_or_overwrite_edits_behave` | Without it, bad ranges or letters could corrupt the sequence or be silently dropped, and empty or overwrite edits could change bases or features they should not. | hand-computed | adequate | (was rejects_bad_edits) Record content after the refused edits is not compared (only the undo stack).<br>(was rejects_bad_edits) Largely redundant: invalid letters are covered more thoroughly by letters_are_kept_as_typed_and_invalid_ones_are_refused_not_dropped and out-of-range replaces by editing_an_empty_sequence; only the inverted range (5,4) is unique here.<br>(was overwrite_keeps_length_moves_nothing_and_is_one_step_per_call) The overwritten bases themselves are never checked (e.g. seq[10..14]=='TTTT').<br>(was editing_an_empty_sequence) gc_bins on an empty range only checks there are 3 values (expected NaN), not their value. |
| `genbank_text_and_qualifiers_survive_saves_and_edits` | Without it, saving could rewrap, unquote or alter qualifier text and a GenBank file would no longer reload as it was written. | round-trip/self-consistency | strong | (was wrapped_qualifiers_are_stable_across_saves) Fixed file name ge_test_wrapped.gb in the shared temp dir (line 1280), never removed.<br>(was genbank_round_trip) Topology and definition are not compared after the round trip. |
| `feature_edits_and_translation_refresh_are_single_undo_steps` | Without it, a feature change or a /translation refresh could take several undos, or one undo could leave the annotations half restored. | hand-computed | strong |  |
| `copy_paste_and_reverse_complement_keep_features_on_their_bases` | Without it, pasted or reverse-complemented annotations could land on the wrong bases or strand. | hand-computed | strong | (was copy_paste_carries_features) The shift of the existing features by 8 after pasting at 0 is not asserted.<br>(was reverse_complement_twice_is_identity) Overlaps reverse_complement_flips_every_feature_onto_the_same_bases (feature flipping on linear and circular, with undo) and complement_join_order_survives_edits_and_reverse_complement (twice = identity); what is unique here is the hand-computed sequence and coordinates on the shared sample.<br>(was pasted_features_cover_the_pasted_bases_and_existing_ones_shift_correctly) Old h and the complement join k are set up but never checked after the paste, although the name says existing ones shift correctly.<br>(was paste_at_every_boundary_matches_the_string_model_and_is_one_step) For (12,14) only the length of g is checked (line 1897), not its bases. |
| `primer_sites_and_pcr_products` | Without it, a PCR product could have the wrong bases or tails, and saving primer sites could change the document or its undo history. | hand-computed | strong | (was saving_replaces_primer_sites_without_touching_the_document) Comment says 'nothing written' for the bad site but the file content after the failed save is not checked. |
| `gc_bins_exact_and_from_blocks` | The GC-content plot could show wrong values, or stale values after an edit or undo on long sequences that use the block cache. | hand-computed | adequate | The block path (lines 1468-1474) is only fed uniform sequences (all GC or all AT), so a wrong block index, off-by-one block edge or bins swapped would still give 1.0/0.0; a mixed sequence is needed to test it. |
| `circular_edits_across_the_origin` | Without it, edits, copies and features across a plasmid's origin could lose or shift bases. | hand-computed | strong | (was circular_edit_across_the_origin_equals_the_rotated_string_edit) Only 'is some rotation' of the model is checked (line 1832), not the exact rotation, and the two genes' positions after the edit are never asserted. |
| `rotating_and_making_linear` | Without it, rotating a circle or making it linear could move features off their bases or not undo back to the circle. | hand-computed | strong | (was making_linear_is_one_undo_step_that_rotates_and_undo_gives_the_circle_back) Line 2396 says linearize(50) 'cuts the feature cut (bases 31..70)', but nothing asserts what happens to 'cut' (two pieces on the right bases), the most interesting case here. |
| `history_texts_and_the_modified_flag` | Without it, the History panel could describe steps wrongly and the unsaved-changes flag could lie, so changes could be lost on close. | hand-computed | strong |  |
| `random_edits_undo_redo_follow_a_snapshot_and_string_model` | Some combination of edits, pastes, rotations, flips, undo, redo and save could corrupt the bases, desynchronise undo/redo, or get the unsaved-changes flag wrong. | in-test reference model | strong | Feature locations are compared only with snapshots taken from the document itself (self-consistency); no independent model checks where features land after each edit (other drift tests cover that). |
| `edited_ranges_follow_edits_undo_redo_the_origin_and_save` | Without it, the blue marking of edited bases could point at the wrong bases after edits, undo, redo or a save. | in-test reference model | strong | (was edited_ranges_follow_every_edit_undo_and_redo) Line 2477 'if pos > 0 && doc.undo()' silently skips when undo() unexpectedly returns false, so a missing undo step would not fail the test (it would desync later, or not at all).<br>(was edited_ranges_basics_origin_and_save) Docstring (lines 2496-2497) says the origin edit 'marks the new bases on both sides', but lines 2524-2528 only check one range of width 3; its position is never compared with where 'GGG' actually is. |

### `core/src/fasta.rs` (4 tests)

| Test | Why it matters | Checked against | Strength | Review notes |
|---|---|---|---|---|
| `a_fasta_rewritten_in_the_same_clock_tick_is_re_indexed` | After a genome FASTA is rewritten quickly, a stale .fai could be trusted and fetches would return the old bases. | hand-computed | strong | Fixed file names in the shared system temp dir (ge_test_same_tick.fa, line 319); concurrent test runs can collide. |
| `fetch_across_line_breaks_with_crlf` | Fetching a region from a Windows (CRLF) or wrapped genome FASTA could return bases shifted by the line endings. | hand-computed | strong |  |
| `replace_one_sequence_and_search` | Saving one edited record into a multi-record FASTA could damage the other records, or motif search in a FASTA could miss or duplicate hits. | hand-computed | strong | Wrapping at 60 (80 new bases should be 60+20 lines) and the in-place case dst == src, both promised by write_replacing's doc, are not tested. |
| `rejects_uneven_wrapping` | An unevenly wrapped FASTA could be indexed anyway, so every fetch from it would silently return wrong bases. | hand-computed | adequate |  |

### `core/src/locations.rs` (3 tests)

| Test | Why it matters | Checked against | Strength | Review notes |
|---|---|---|---|---|
| `delete_rules` | Deleting bases could leave features shifted, wrongly trimmed, or not removed when fully deleted, especially multi-part complement features. | hand-computed | adequate | No case for a deletion after the range (unchanged) or overlapping only the end. |
| `spans_and_strand` | Feature locations could be flattened into wrong ranges or the wrong strand, which every drawing and export relies on. | hand-computed | adequate | Between (zero-width site) handling is not tested. |
| `clip_window_keeps_parts_in_reading_order` | Copying or rotating across the origin could scramble the order of a feature's parts, lose its strand, or keep parts outside the copied window. | hand-computed | strong | Fuzzy-marker dropping on cut ends (locations.rs line 159) is not tested. |

### `core/src/primers.rs` (5 tests)

| Test | Why it matters | Checked against | Strength | Review notes |
|---|---|---|---|---|
| `tm_and_tm_range_of_plain_and_degenerate_primers` | Without it, the Rust Tm and the Tm range shown for degenerate primers could be implausible or not bracket the oligos they stand for. | property/invariant | weak | (was tm_is_sane) 10 °C window only; exact values are checked in Python by tests/test_core_sequence.py::test_melting_temperature_matches_santalucia_and_biopython. The code comment in the test still points to tests/test_primers.py, which no longer exists.<br>(was tm_is_sane) Comment says 'more salt or more primer = higher Tm' but primer concentration is never varied; the self-complementary (symmetric) branch is untested.<br>(was tm_is_sane) Superseded by tests/test_core_sequence.py::test_melting_temperature_matches_santalucia_and_biopython (exact Biopython values, the same input checks, salt/Mg/dNTP/primer direction); it only remains as a fast Rust smoke test.<br>(was tm_range_of_a_plain_primer_is_its_tm_and_bad_letters_give_none) Duplicates the end of tests/test_core_sequence.py::test_melting_temperature_matches_santalucia_and_biopython (tm(degenerate) is None, tm_range(plain) == (tm, tm), tm_range('ACGTX') is None); the only extra is tm_range('A'), and tm_is_sane already checks tm('ACGTN') is None.<br>(was tm_range_is_exactly_the_lowest_and_highest_tm_of_the_expanded_oligos) For the two self-complementary primers only bracketing is checked, so an arbitrarily wide range would pass there.<br>(was tm_range_is_exactly_the_lowest_and_highest_tm_of_the_expanded_oligos) Overlaps tests/test_core_sequence.py::test_melting_temperature_matches_santalucia_and_biopython, which checks bracketing for 'AGCGGRTAACAATYTCACACNGGA' only; this Rust test is the stronger one (exact bounds, 4 primers). |
| `a_degenerate_primer_binds_every_site_it_matches_on_either_strand` | Without it, a degenerate primer could miss sites it matches, or bind sites it does not, on either strand. | hand-computed | adequate | (was a_degenerate_primer_binds_every_site_it_matches_with_no_mismatch) Near-duplicate of the degenerate-primer block of tests/test_primer_binding.py::test_binding_origin_ends_repeats_and_letters (three variants with R/Y/N, a base outside the code is a mismatch, N in the template matches only N), which also checks the reverse strand; one of the two is redundant.<br>(was a_degenerate_primer_binds_every_site_it_matches_with_no_mismatch) The test binds the same string to two names (let (a, g) = (same, same)): harmless but confusing.<br>(was a_degenerate_primer_binds_every_site_it_matches_with_no_mismatch) Degenerate position 2 lies outside the 15-nt 3' seed, so it only matters through the extension; the negative cases test position 12 (inside the seed) only.<br>(was a_degenerate_reverse_primer_binds_the_opposite_strand_using_complemented_codes) One degenerate letter, one positive case and no negative case (an R facing a base it does not pair with), so code that let every IUPAC letter match anything would pass.<br>(was a_degenerate_reverse_primer_binds_the_opposite_strand_using_complemented_codes) The code comment 'opposite a T on the top strand?' (with a question mark) is wrong: primer base 5 faces top-strand position 34, a C, which pairs with G (allowed by R).<br>(was a_degenerate_reverse_primer_binds_the_opposite_strand_using_complemented_codes) Redundant with tests/test_primer_binding.py::test_binding_origin_ends_repeats_and_letters, which binds a three-letter degenerate primer to the reverse-complemented template (three sites, all reverse strand, 0 mismatches). |
| `annealing_extends_across_a_lone_mismatch` | A primer's displayed binding site could stop at a lone mismatch (showing too much tail) or run on past two mismatches, on either strand. | hand-computed | strong |  |
| `finds_sites_on_both_strands_with_tails_and_mismatches` | Primers with 5' tails, reverse primers, or primers with mismatches could be shown binding wrongly, including a 3'-mismatched primer shown as binding. | hand-computed | strong |  |
| `site_across_the_origin` | A primer binding across the origin of a plasmid could be missed. | hand-computed | strong | No linear counterpart (that the site is not found / is shorter on a linear template). |

### `core/src/search.rs` (6 tests)

| Test | Why it matters | Checked against | Strength | Review notes |
|---|---|---|---|---|
| `finds_an_exact_motif_on_both_strands_in_position_order` | Sequence search could miss bottom-strand hits, report palindromes twice, ignore the top-strand-only option, or be case-sensitive. | hand-computed | strong | Comment line 236 says 'GGATCC too' but GGATCC is never searched. |
| `iupac_letters_in_the_query_match_every_base_they_stand_for` | IUPAC search queries could miss bases they stand for, or an N in the sequence could match everything, giving false or missing hits. | hand-computed | strong |  |
| `allows_a_number_of_mismatches_and_reports_how_many` | Search with mismatches allowed could miss near-matches or report a wrong mismatch count. | hand-computed | strong |  |
| `a_hit_can_run_across_the_origin_of_a_circle` | A motif crossing the origin of a plasmid could be missed by the search. | hand-computed | strong |  |
| `stops_at_the_limit_and_handles_tiny_and_empty_inputs` | A search could run past its hit limit or crash on an empty or very short sequence or query. | hand-computed | adequate | Limit test does not check that the first 3 hits (0, 2, 4) are the ones kept. |
| `finds_a_protein_motif_in_all_six_frames` | Protein motif search could miss hits in some frames, report wrong DNA coordinates, or mishandle X wildcards and stops. | hand-computed | adequate | Name says 'all six frames' but only the three forward frames are checked; the reverse-strand frame number is never asserted (line 309). |

### `core/src/tracks.rs` (1 test)

| Test | Why it matters | Checked against | Strength | Review notes |
|---|---|---|---|---|
| `gff_bed_and_queries` | Loaded GFF/GTF/BED/VCF annotation tracks could show features off by one base, lose escaped names, or miss features in a region. | hand-computed | adequate | Fixed file names ge_test.gff3/.bed/.vcf in the shared temp dir, never removed.<br>Overlaps tests/test_core_formats.py::test_tracks_gff_bed_parsing_query_density_and_search (GFF/GTF/BED with a brute-force query check); only the VCF part is covered here alone. |
