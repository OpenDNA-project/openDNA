# What the tests check, and how

A plain-language catalogue of every automated test in openDNA, written for the bench scientist who wants to know
which claims the program has been held to. `docs/TESTING.md` is the *plan* (what should be covered); this file is the
*inventory* (what each test actually does, and whether its biology is right).

Run everything:

    python -m pytest tests -q          # 271 test functions (many hold several sections), about 4 min
    cd core && cargo test              # 95 tests, about 2 s

The input and expected output of every test, its result in the last run and the audit of November 2026 are in
[docs/test-report/](test-report/README.md).

The tests that need real sequence files use the copies in `demo/` (see `tests/sample_files.py`), so they run on any machine
that has the repository. Three tests skip on a machine without extra files or tools, CI included: the read-mapping part of
`test_gui.py::test_read_colours_local_align_and_map_reads` (needs bwa, samtools and minimap2), and
`test_map_view.py::test_human_genome_chromosomes_listed_drawn_and_zoomed` and `test_multi_genbank.py::test_vibrio_cholerae_two_chromosomes`
(need a human genome and a *V. cholerae* GenBank file that are not in the repository).

**How the tests get their DNA.** Four kinds of input, and the difference matters when you judge a result:

| Input | Used for | Why |
|---|---|---|
| Hand-built sequences (e.g. `CCATGAAAAAATAACC`) | exact, checkable answers | every base is known, so the expected coordinates can be worked out on paper |
| Random sequence with something planted in it | "does it find the needle, and only the needle" | guards against a rule that happens to work on one favourite sequence |
| Real plasmids and a real genome (`demo/`) | opening, saving, primer sites, digests | catches what synthetic data never shows: odd qualifiers, joins, origin-crossing features |
| An independent second calculation (Biopython, or a plain-Python model written inside the test) | the hard arithmetic | the test works the answer out a second way and the two must agree |

That last row is the backbone of the biology testing: melting temperature, restriction digests, translation, GenBank
locations, read pileups and primer binding are all checked against an independent calculation, not against a number
someone once wrote down.

---

## 1. Melting temperature, hairpins and primer-dimers

`tests/test_core_functions.py`, `tests/test_primers.py`, `core/src/primers.rs`

| Test | What it checks | How |
|---|---|---|
| `test_tm_three_known_primers_against_the_santalucia_1998_formula_and_biopython` | Tm follows SantaLucia 1998 nearest-neighbour thermodynamics | The ΔH/ΔS table is written out again inside the test and the formula applied by hand at 1 M Na⁺, 250 nM primer; the program must agree to 1e-6. Biopython's `Tm_NN` must agree within 0.5 °C. Three real primers (M13 −20, an M13/pUC reverse 24-mer, the T7 promoter primer). |
| `test_primers.py::test_tm_matches_biopython` | The same at realistic PCR conditions | 200 random 12–40-mers at three salt / Mg²⁺ / dNTP / primer combinations must match Biopython within 0.01 °C. Palindromes are skipped (self-complementary strands use a different symmetry term). A sequence containing N gives no Tm. |
| `test_tm_rises_with_length_and_gc` | The direction of the physics | Tm rises monotonically over prefixes of length 10→40; at 20 nt, AT-only < mixed < GC-only. Case and U are ignored. |
| `test_tm_salt_magnesium_dntp_and_primer_concentration_directions` | Buffer effects go the right way | More Na⁺, more Mg²⁺ or more primer raises Tm; dNTPs lower it by chelating Mg²⁺. With Mg²⁺ ≤ dNTP the result equals the no-Mg²⁺ case. |
| `test_hairpin_and_self_dimer_of_palindromes_and_none_for_poly_a` | Structures that should form, do | A `GGATCC`-loop-`GGATCC` construct forms a 3′ hairpin; `GAATTCGAATTC` self-dimerises at its 3′ end; poly-A, too-short and empty inputs give 0. |
| `test_dimer_and_hairpin_energy_signs_and_ordering` | Secondary-structure ΔG behaves sensibly | 300 random pairs plus hand-made stems: ΔG ≤ 0 always; the best dimer anywhere is at least as stable as the best 3′ one; dimers are symmetric; longer and GC-richer duplexes are more stable; longer loops are less stable; an internal hairpin is not reported as a 3′ hairpin; more salt is more stable. |
| `primers.rs::tm_is_sane`, `primers.rs::dimers_and_hairpins` | The Rust side agrees | Tm of the 24-mer is 55–65 °C (it is 60.8); a primer against its own full complement gives ΔG < −20 kcal/mol, about −1 per base pair. |

**Checked against the literature and correct:** the 10 unique nearest-neighbour steps and both initiation terms
(SantaLucia 1998); the salt term `0.368·(N−1)·ln[Na⁺]`; the Mg²⁺→Na⁺ conversion `Na_eq = Na + 120·√(Mg − dNTP)`
(von Ahsen 2001); the C/4 rule for non-self-complementary duplexes; the hairpin loop table (SantaLucia & Hicks 2004,
loops 3–10).

## 2. Open reading frames and translation

| Test | What it checks | How |
|---|---|---|
| `test_orfs_minimum_length_stop_coordinates_and_strands` | An ORF runs ATG→stop, the stop belongs to the ORF, and bottom-strand ORFs are reported in top-strand coordinates | `CC ATG AAA AAA TAA CC` gives (2, 14, +1) at a 4-codon minimum and nothing at 5; its reverse complement gives the same span on −1. Each of TAA, TAG, TGA ends an ORF, lower case included. |
| `test_orfs_nested_back_to_back_open_and_missing_start_or_stop` | Only the longest ORF per stop (first ATG); adjacent ORFs both found; nothing without a start or a stop | Tiny hand-built sequences with exact expected tuples. |
| `test_orfs_agree_with_a_brute_force_scan_on_random_sequences` | Six-frame ORF finding is right in general | Random 500 and 3001 nt at minimum lengths 1, 10, 30 must equal a plain six-frame Python scan. In a 6 kb sequence every reported ORF starts with ATG, ends with a stop, and has no internal in-frame stop. |
| `test_orf_across_the_origin_of_a_circular_sequence`, `bio.rs::finds_an_orf_that_crosses_the_origin_of_a_ring` | A gene running through base 1 of a plasmid is still one ORF | `AAATAA` + 30×C + `ATGAAA` as a circle reads `ATG AAA | AAA TAA` and gives exactly one ORF; the Rust test pins the coordinates (36, 48) and the mirror on the reverse strand. A linear copy gives nothing. |
| `test_translation_of_all_64_codons_against_the_standard_table` | The genetic code | All 64 codons against Biopython table 1, forward and reverse strand, upper and lower case; stops become `*`. |
| `test_translation_incomplete_last_codon_and_codons_with_n` | Ragged ends and ambiguity | A trailing partial codon is ignored; frame 2 of `ATGGCCAAAGG` reads W P K; a codon containing N is X. |
| `code.rs::the_standard_code_and_the_known_differences` | The same in Rust | `ATG GCC TAA GNA` → `MA*X`. |

## 3. Restriction enzymes and digests

| Test | What it checks | How |
|---|---|---|
| `test_cuts_and_overhangs_on_a_hand_made_sequence` | Where each enzyme cuts each strand, and the overhang left behind | EcoRI gives top/bottom cuts (11, 15) with a 5′ `AATT` overhang; KpnI (31, 27) with a 3′ `GTAC` overhang; EcoRV (45, 45) blunt. Every common enzyme is then planted in both orientations and must produce its declared offsets, mirrored on the bottom strand. |
| `test_both_strands_and_sites_across_the_origin` | Type IIS enzymes cut *outside* their site, on the correct side | BsaI forward (27, 31), reverse (15, 19), 4-nt 5′ overhang. An EcoRI site split across the origin cuts at (25, 3) on a circle and not at all on a line. A BsaI cut falling past the end of a linear molecule is dropped; on a 7 bp circle it wraps. |
| `test_degenerate_sites_and_both_sided_cutters` | Ambiguous recognition sequences | AccI (`GT^MKAC`) accepts the four legal variants and rejects three illegal ones (M = A/C, K = G/T); BsmBI (17, 21); BglI (`GCCNNNN^NGGC`) gives (9, 6), a 3′ overhang. |
| `test_digest_on_the_real_plasmid_equals_biopython_restriction` | A real digest equals Biopython's | Every common enzyme on the real .dna plasmid, circular and linear; cut positions must equal `Bio.Restriction`, and the fragments must sum to the plasmid length. |
| `test_digest_with_every_enzyme_planted_in_a_random_circle_equals_biopython` | The same under adversarial conditions | Five random 4 kb circles with EcoRI and BsaI sites planted *across the origin*; still equal to Biopython. |
| `test_cut_sites.py::*` | Hover text of cuts, the cut marker, the Cas9 guide display and its storage | Tooltips for EcoRI (5′), KpnI (3′), EcoRV (blunt), BsaI (Type IIS, both strands), a site across a circle's origin, SpCas9 on both strands, Cas12a; the cross at a cut is painted and the origin bar has none; the protospacer / PAM / name colours are painted in both themes; a guide saved by the Spacer Finder reloads byte-identically, Biopython reads it as a `misc_feature`, exports leave the `/OpenDNA_cas9` qualifier out. |
| `test_cut_ui.py::*` | The enzyme pill tooltip, the guide look, the Spacer Finder's automatic check and right-click, the Digest fragment click | Pill pieces for 5′ (EcoRI), 3′ (KpnI), blunt, Type IIS, bottom-strand site, circular origin and a far cut with …; pixmaps drawn in both themes and different per kind; the tooltip carries the pill and Qt renders it; guide arrow is red, hollow-outlined and points up (pixel test), guide text role is black, text fits lane height from font metrics at every size; first 100 guides checked with a faked reference (101st onward not), cancel, no-reference header, sorting maps rows to guides; right-click adds a feature identical to the button's (also for several rows); a Digest fragment click selects its span (linear, across the origin, whole circle), pins its end cuts with the pill tooltip and clears on another selection. |
| `test_view_edges.py::test_enzymes_know_where_each_strand_is_cut` | The cut drawn on screen is on the right strand | EcoRI `G^AATTC` → (3, 7); EcoRV blunt → (5, 5); BsaI forward (9, 13), reverse site (7, 11); EcoRI across a 12 bp origin → (11, 3). |
| `test_bugfixes.py::test_a_cut_that_falls_off_the_end_of_a_linear_sequence_is_no_cut` | A Type IIS enzyme needs room downstream to cut | `GGTCTCA` alone on a linear molecule gives no cut; with a longer tail it cuts at (7, 11); as a circle it cuts. |
| `test_gui.py::test_design_pcr_orfs_digest` (digest part) | The virtual gel | BamHI + PvuI on the real plasmid give 2 fragments whose sizes sum to the plasmid length, and two gel lanes. |

**Checked and correct:** all 47 enzyme entries in `app/enzymes.py` (recognition sequence and both cut offsets), the
AccI degeneracy, and the reverse-strand offsets of the Type IIS enzymes.

## 4. Primer binding — where a primer really anneals

This is the most heavily tested part of the program, because it decides what the Primers panel draws and what a PCR
product looks like. The rule being tested: the 3′ end must match exactly (the last 3 bases), a seed next to it is
checked for mismatches, annealing then extends toward the 5′ end and stops when it has lost too much ground, and
anything left over at the 5′ end is a non-annealing tail (a restriction site, a Gibson overlap).

| Test | What it checks | How |
|---|---|---|
| `test_the_rust_search_agrees_with_the_rules_written_in_python` | The whole rule, in general | 1500 random templates containing repeats, linear and circular, with mutated primers carrying random 5′ tails, seeds 4–20 and mismatch limits 0–3, all compared with an independent Python implementation of the documented rule. |
| `test_every_single_mismatch_position_is_reported_with_count_1_on_both_strands` | What one mismatch does, position by position | For lengths 16, 20, 24, 30 on both strands: a mismatch in the last 3 bases means no binding at all; in the seed it counts as 1; further 5′ it is crossed if enough matching bases lie beyond it; within 3 bases of the 5′ end it simply becomes part of the tail. |
| `test_a_tail_of_any_length_0_to_30_is_not_counted_and_does_not_move_the_site` | A 5′ tail does not shift the binding site | Tails of 0–30 nt made of bases that cannot pair; the reported site and mismatch count must be unchanged. Both strands, both topologies. |
| `test_a_tail_that_matches_the_template_extends_the_site_and_the_count_stays_right` | A 5′ extension that *does* match anneals too | 0–30 extra matching bases; a single mismatch inside the extension is crossed and counted. |
| `test_the_extension_stops_three_below_the_best` | The X-drop extension rule | Ten hand-made match/mismatch patterns × both strands, each worked out by hand. |
| `test_site_across_the_origin_on_both_strands`, `test_a_mismatch_on_each_side_of_the_origin` | A primer that spans base 1 of a plasmid | A 12 + 12 split gives the span (n−12, n+12) on a circle and nothing on a line; rotating the plasmid moves the site accordingly; mismatches on either side of the origin are both counted. |
| `test_tandem_repeat_sites_overlap_but_are_never_duplicated`, `test_palindromic_primer_binds_both_strands_at_the_same_place_once_each` | Repeats and palindromes | `(GATTACA)×8` with a 21-mer gives 6 overlapping sites 7 nt apart; a 22-nt palindrome gives exactly one site per strand, never duplicated. |
| `test_a_site_and_its_mirror_image_are_the_same_site_on_the_other_strand` | Forward and reverse rules are mirror images | 600 random cases: the sites found on the reverse complement must equal the mirrored forward sites. |
| `test_primer_pair_product_linear_with_tails`, `..._across_the_origin_with_a_site_across_it` | The PCR product | Product = forward tail + template between the primers + reverse-complemented reverse tail; a mismatched primer writes *its own* base into the product, as in a real PCR; on a plasmid the product may run through the origin. |
| `test_reverse_primer_upstream_of_the_forward_one_has_no_product_when_linear_but_runs_through_the_origin_on_a_circle` | Orientation decides whether there is a product at all | Two primers pointing away from each other give no product on a linear template and a 228 bp product on a circle; two forward primers never pair. |
| `test_the_pair_with_the_shortest_product_is_picked` | Which pair the program offers | Two forward and two reverse primers; the shortest product (124 bp) wins. |
| `test_non_specific_priming_in_the_sequence_...`, `..._in_the_genome_check` | Mispriming warnings | A 12-nt 3′ match elsewhere is a yellow warning; a second full-length site is red; against a reference genome, unique is green, a 3′ copy yellow, a full reverse-complement copy red, and the sequence is never compared with itself. |
| `test_our_sites_are_the_binding_sites_stored_in_each_real_dna_file` | Agreement with the stored binding sites | Every site found in the real .dna file must be one stored in that file, or sit inside one with the same 3′ base. |
| `test_primer_sites_written_on_save_do_not_change_the_document_or_its_undo_stack` | Saving writes primer sites without touching the document | Two saves are byte-identical; sites shift by +10 after a 10 nt insertion; no undo step is added. |

## 5. Reading and writing files

`tests/test_roundtrip.py`, `tests/test_real_files.py`, `tests/test_bugfixes.py`, `core/src/document.rs`,
`locations.rs`, `dna_file.rs`, `fasta.rs`, `tracks.rs`

The question behind this whole group: **can the program lose or silently alter your annotation?**

| Test | What it checks | How |
|---|---|---|
| `test_load_save_as_genbank_reload_is_identical` | Opening and saving loses nothing | The real .dna file, the real plasmid, the three demo files and the 4.6 Mb genome: name, definition, topology, every base, and every feature (type, location text, qualifiers in order) must come back identical. |
| `test_save_reload_save_is_byte_identical` | No drift between saves | Save, reload, save, reload, save: the files must be byte-for-byte equal, so a file does not creep as you work on it. |
| `test_file_saved_by_us_is_read_the_same_by_biopython`, `test_rich_document_is_read_the_same_by_biopython` | Other programs read our files the same way | Biopython must report the same sequence, topology, feature count, types, coordinates, strands and qualifiers — including `complement(join(...))` and the fuzzy `<100..>120` form. |
| `test_location_forms_survive_a_genbank_round_trip` | Every INSDC location syntax | 12 forms (simple, `complement`, `join`, `complement(join)`, join of complements, mixed-strand join, `<`, `>`, `<>`, `order`, a 3-part join, and `join(350..400,1..20)` across the origin) must survive a save unchanged, and Biopython must agree on the parts. |
| `test_the_gap_in_a_between_location_has_the_same_span_as_biopython` | The "between two bases" form `20^21` | Meaning kept; Biopython reports start = end = 20, the same as us. |
| `test_location_parse_and_format_round_trip` (Rust) | The location grammar | 20 forms including `order`, `bond`, one-of, a remote accession, `gap()` and `200^1` across the origin; 1-based→0-based conversion, strand and the exact text. |
| `test_awkward_qualifier_text_survives_...` and the quote/wrap tests | Qualifier text is never corrupted | Quotes, newlines, Unicode, very long values and valueless flags; a quote is tried at all 120 wrap offsets; a 300-character unbreakable word is never split; a value that ends exactly at the wrap column is read back without a stray space. |
| `test_lower_and_mixed_case_bases_are_kept_as_typed...` | Case is information | Mixed-case sequence (the convention for marking non-coding regions) is preserved through load, slice, edit and save. |
| `dna_file.rs::reads_synthetic_file` | The `.dna` format | Segment `2-10` (1-based) → `2..10`; `directionality=2` with range `18-3` on a 20 bp circle → `complement(join(18..20,1..3))`; binding sites are 0-based, so `9-12` → `complement(10..13)` (the same convention Biopython's .dna reader uses); duplicate "simplified" sites skipped. |
| `fasta.rs`, `test_roundtrip.py` FASTA group | Indexed FASTA and in-place editing | `.fai` lengths and 0-based half-open fetches across 60-column line breaks; CRLF, missing final newline, blank lines, mixed case and one-long-line variants all give the same bases; editing one record rewrites only that record, byte for byte, keeping the file's own line endings and a UTF-8 BOM. |
| `tracks.rs::gff_bed_and_queries`, `test_gff3_...`, `test_bed_with_3_to_12_columns` | The coordinate conventions of each track format | GFF3/GTF are 1-based inclusive (gene `10..100` → internal start 9); BED is 0-based half-open; VCF POS is 1-based (`POS 100 A>G` → (99, 100); `POS 200 ACGT>A` ends at 203, the span of the REF allele). `##FASTA` stops parsing; `%2C` escapes decoded; malformed lines are errors, not crashes. |
| `test_a_truncated_file_is_a_clear_error_never_a_crash_or_a_silent_partial_document` | A damaged file never becomes a silently half-read sequence | The real `.dna` and a demo `.gb` truncated at more than 100 different points: each must either raise a readable error or open complete. |
| `test_a_genbank_file_with_two_records_is_refused...` | Multi-record GenBank is refused rather than silently dropping one | A demo file concatenated with itself raises "2 GenBank records". |
| `test_a_failed_save_leaves_the_original_file_and_the_document_alone` | A failed save never damages what you had | An out-of-range primer, a missing folder and a target that is a directory all raise, leaving the original bytes and leaving no `.tmp` behind. |
| `test_roundtrip.py::test_autosave_and_recovery` | Crash recovery never writes your file | Autosave writes only `<file>.autosave.gb`; the original bytes are untouched; after a simulated crash the app offers recovery once; declining or a stale copy cleans up. A `.dna` file is never written at all. |
| `test_reference.py::test_the_check_sees_copies_on_the_second_and_third_record` | The mis-priming check covers every FASTA record | 3-record FASTA with the primer on chr2 and its reverse complement on chr3: red, hits named `chr2:` and `chr3:`, label "3 sequences". |
| `test_reference.py::test_a_three_prime_copy_on_a_later_record_is_yellow_and_a_clean_primer_green` | Verdict colours per record | 3' copy only on chr3 is yellow; nothing anywhere is green "in 3 sequences". |
| `test_reference.py::test_the_open_record_is_never_compared_with_itself` | No self-match | Open chr2 of the reference: only chr3 reported; a one-record reference that is the open sequence says "this sequence is the genome". |
| `test_reference.py::test_a_long_record_is_searched_in_overlapping_windows_without_losing_or_doubling_a_site` | Windowed search of a long record | A site across a window edge is found once; windowed and whole results are equal. |
| `test_reference.py::test_a_multi_record_genbank_reference_is_read_one_record_at_a_time` | Multi-record GenBank as a reference | Two concatenated records: still refused as a document (message names the reference option), read in full as a reference, topology kept per record. |
| `test_reference.py::test_a_search_that_runs_out_of_time_or_is_cancelled_says_it_is_not_conclusive` | Time bound and Cancel | Zero time stops at once; cancel stops after two records with progress reported; a partial search with no hit is yellow, not green. |
| `test_reference.py::test_the_primer_dialog_searches_the_genome_on_a_worker_thread` | Dialog stays live | The report fills in with hits on chr2 and chr3 after the thread finishes. |
| `test_reference.py::test_the_guide_finder_counts_the_sites_of_every_record` | Spacer Finder over all records | One perfect copy on each of two records adds to 2 perfect sites: red. |
| `test_reference.py::test_ncbi_members_are_joined_with_their_newlines_and_the_reference_is_the_fasta_itself` | NCBI download and "use as reference" | A member without a final newline does not glue headers; a 2-sequence assembly can be the reference; the setting is the FASTA itself ("2 sequences"). |
| `test_reference.py::test_a_contig_name_matches_ignoring_chr_and_case` | Track contig names | chr1 = 1 = Chr1 when unambiguous; no guess when two fit or none. |
| `test_reference.py::test_tracks_follow_the_record_and_a_mismatching_single_contig_is_said` | Tracks on other records | A track whose only contig differs shows a status message; switching record carries the tracks (by contig name) and names those with no data. |
| `test_reference.py::test_saving_one_record_as_genbank_says_the_others_stay_behind` | Save as from a multi-record FASTA | The status bar says "only chr1: the other 2 sequence(s)... not included". |
| `test_add_chromosome.py::test_one_record_genbank_becomes_multi_record` | Import as New Chromosome makes a multi-record GenBank file | A plasmid file gets a second plasmid as a new record: the old bytes are the start of the new file, the added record has the same bases, topology and features, the window shows it, and the first record comes back with its unsaved edit, which Save writes in place. |
| `test_add_chromosome.py::test_names_are_unique_and_can_be_changed` | Record names stay unique | The same plasmid added to itself becomes `<name>_2`; a typed name is used (spaces become `_`); a name the file has is refused and nothing is written. |
| `test_add_chromosome.py::test_multi_fasta_source_and_target` | New sequences in a FASTA file | Ticked records of a multi-FASTA (a clashing name renamed) and a GenBank plasmid (bases only) are appended; the stale `.fai` is rebuilt with every length; a Windows file keeps CRLF only; a last line without a newline is ended, not joined. |
| `test_add_chromosome.py::test_dna_source_keeps_primers_and_dna_target_is_saved_first` | `.dna` files | The primers of a `.dna` source become `primer_bind` features of the new record; an open `.dna` file is saved as GenBank first and never changed. |
| `test_add_chromosome.py::test_cancelled_dialog_writes_nothing` | Cancel | The file's bytes and the window are unchanged. |


## 6. Editing: do annotations stay on the right bases?

`core/src/document.rs`, `core/src/locations.rs`, `tests/test_editing.py`, `app/selftest.py`

| Test | What it checks | How |
|---|---|---|
| `test_editing.py::test_500_random_edits_keep_bases_and_features_then_undo_redo_exactly` | The central guarantee, stress-tested | 10 seeds × linear and circular × 500 random edits. A plain-string model tracks the bases and a per-base identity model tracks which original base each annotated position should still be sitting on; both must agree at every step, and undo/redo must walk back and forth exactly. |
| `document.rs::random_edits_undo_redo_follow_a_snapshot_and_string_model` | The same inside Rust | 600 random operations including paste with features, rotate, reverse-complement and edits across the origin. |
| `insert_at_the_exact_start_or_end_of_a_feature_does_not_grow_it` | The boundary rule | Inserting at a feature's first base shifts the whole feature; inserting just after its last base leaves it alone. |
| `delete_overlapping_the_start_the_end_or_the_middle_keeps_the_remaining_bases`, `deleting_exactly_the_feature_or_more_removes_it...` | Trimming and deletion | Six delete ranges; a feature deleted in full disappears and comes back whole on undo, qualifiers included. |
| `replace_over_part_of_a_feature_follows_the_documented_rules` | Replacement bases join a feature only when strictly inside it | Inside (longer, same, shorter), covering the start, covering the end, covering all. |
| `join_edit_between_parts_grows_the_gap_inside_a_part_grows_that_part` | Multi-exon features | Eight insert positions in `join(21..40,81..100)` with the expected spans; a delete spanning the gap trims both parts. |
| `complement_feature_keeps_its_bases_and_strand_through_edits`, `complement_join_order_survives_edits_and_reverse_complement` | Reverse-strand and multi-part features | Strand stays −1, part order survives, and the feature's own sequence is unchanged by edits elsewhere. |
| `reverse_complement_flips_every_feature_onto_the_same_bases` | Flipping the whole molecule | Four feature shapes, linear and circular: every span mirrors to (n−end, n−start), strand flips, and the feature's own sequence is the reverse complement of what it was. Twice = identity. |
| `feature_crossing_the_origin_keeps_its_bases` | Origin-crossing features survive edits and rotation | `join(191..200,1..10)` on a 200 bp circle: inserts inside either part grow it, and twelve different rotations all keep its sequence and part order. |
| `fuzzy_locations_survive_edits_and_saves_with_their_markers` | Uncertain ends keep their `<` and `>` | `<21..>60` shifts to `<24..>63` and then `<24..>65`; `complement(<101..>150)` likewise. |
| `making_linear_is_one_undo_step_that_rotates_and_undo_gives_the_circle_back` | Cutting a plasmid open | Linearising at base 51 rotates the sequence and sets the topology in a single undo step; features keep their own bases; undo restores the circle. |
| `qualifiers_are_untouched_by_edits` | Notes and translations are never rewritten by an edit | 60 random edits on a join CDS carrying `/note`, `/translation` and `/pseudo`. |

## 7. Reads, coverage and variant calling

`tests/test_core_functions.py` (BAM group), `core/src/alignments.rs`, `tests/test_gui.py`

| Test | What it checks | How |
|---|---|---|
| `test_bam_reads_soft_clips_deletions_insertions_strand_and_ignored_reads` | CIGAR handling | Hand-built reads: unmapped, duplicate and QC-fail reads are hidden; a soft clip does not shift the start; `20M5D25M` spans 200–250 with the deletion at 220–225; `20M3I27M` ends at 347 with the insertion at 320; mate position, mate on another contig and TLEN are read. |
| `test_bam_pileup_hand_count_on_the_designed_reads` | The pileup column counts | A, C, G, T, N, deleted and inserted-before counted by hand at the soft clip, the deletion and the insertion; then the whole contig against a Python model; `min_baseq` drops a Q10 base and `min_mapq` drops a MAPQ-5 read. |
| `test_bam_pileup_and_coverage_equal_the_model_on_random_reads` | The same on messy data | 300 random reads with clips, indels, sequencing errors, flags and qualities, under five filter combinations. |
| `test_variants_planted_snp_at_half_is_called_and_thresholds_exclude_it` | Allele-fraction calling | Ten reads carrying a SNP at 50 %, 10 % and 30 %; the `min_fraction` and `min_depth` boundaries are probed from both sides (depth 10 passes, 11 fails). |
| `test_variants_reference_base_is_compared_and_reported_correctly` | Ref and alt come from the reference you supply | Swapping the reference base at position 120 changes the call; case is ignored; a reference N is never called a SNP. |
| `test_variants_deletions_and_insertions_and_the_random_reads` | Indel calls | A 5 nt deletion is merged into one call carrying its reference bases; a `GGG` insertion is called; soft clips are never variants. |
| `alignments.rs::reads_mismatches_and_coverage` | BAM↔SAM conventions | A read stored 0-based at 30 is written as SAM `POS 31` with CIGAR `5M2I5M`; a G mismatch at 17 on a reverse-strand read is recognised (SAM stores SEQ in forward orientation). |
| `test_gui.py::test_map_reads_and_bam_track` | The whole mapping workflow (needs bwa/samtools) | Simulated paired reads carrying a SNP and a 5 bp deletion are mapped with bwa mem; both are called as consensus changes; applying them turns the reference into the sample; a 4.1 kb insert pair is flagged; mean depth is 12–14× as designed. |
| `test_gui.py::test_read_colour_categories` | IGV-style read colouring | A normal FR pair is grey; insert 4000 is "larger", 40 is "smaller"; LL/RR/RL orientations, mate elsewhere and mate unmapped each get their colour. |

## 8. Local alignment and finding known parts

| Test | What it checks | How |
|---|---|---|
| `test_align_query_on_the_reverse_strand_reports_strand_and_coordinates` | Alignment reports the right span and strand | A 300 bp exact query with a 50 nt non-matching flank, forward and reverse-complemented. |
| `test_align_indels_are_found_at_their_exact_positions` | Gaps are placed exactly | A 5 nt deletion and a `GGG` insertion placed where the gap cannot slide; one mismatch gives identity 399/400. |
| `test_align_query_longer_than_the_template_and_across_the_origin` | Edge cases | A query longer than the template is trimmed; a query spanning the origin aligns on a circle only. |
| `test_align_query_hanging_off_the_end_of_a_linear_template_aligns_only_the_overlap` | A query overhanging a line | The query's first 40 bases are the template's last 40 (and the mirror: its last 40 are the template's first 40), the rest random: one hit stopping at the template end, strand and query part right on both strands; a query made of end + start gives two hits on a line. |
| `test_align_overhanging_query_on_a_circular_template_wraps_across_the_origin` | The same on a circle | End + start of the template aligns completely across the origin (end 1060 > length) on both strands; a random overhang still aligns only its overlap (a few chance bases at most). |
| `test_align_query_whose_ends_match_but_whose_middle_does_not` | Only the extremities match | Same diagonal and a 150 base middle: both ends are returned, the better one first (this used to give only the better one); different diagonals: both ends; a 10 base middle is bridged as one hit with mismatches, a 60 base one is not (two hits). |
| `test_align_end_overlap_below_the_reporting_threshold_gives_no_hit` | Tiny end overlaps | 5, 10 and 15 base overlaps on either end, either strand, line or circle: no hit; exactly 20 bases (score 40) is the first one reported. |
| `align.rs::query_hanging_off_the_template_end_aligns_only_the_overlap`, `overhanging_query_wraps_on_a_circle_and_a_tiny_overlap_is_no_hit` | The same in Rust | Overhang at either end, both strands, circle versus line, and a 15 base overlap giving no hit. |
| `test_every_piece_of_a_structural_variant_query_is_found` (12 cases: 6 layouts x line / circle) | A patchwork query: A + 200 unmatched bp + B + C inverted + D | Random 2000 bp template; pieces far apart, shuffled, close together, A-B-D on one diagonal, touching in the template: exactly 4 hits, each piece's strand, query range and template range right (within 2 bases of a chance extension), also for the query typed as its reverse complement. |
| `test_a_piece_that_runs_over_the_origin_of_a_circle_is_found_with_the_others` | A piece across the origin | Piece D is the template's last 50 + first 50 bases: one hit with `end` past the length on a circle; two halves on a line, the other pieces unaffected. |
| `test_max_hits_keeps_the_best_pieces_and_pieces_are_not_merged_over_a_mismatching_gap` | `max_hits` and no bridging | max_hits 0/1/2 give 0/1/2 hits; with A, B, D on one diagonal no hit is longer than about 100 bases (nothing merged over the 200 base gap). |
| `test_a_duplicated_piece_in_the_query_is_reported_twice_on_the_same_template_place` | Tandem duplication | The same 100 bp twice in the query: two hits with the same template span and different query ranges. |
| `test_a_piece_with_a_few_differences_next_to_an_exact_piece_is_still_two_pieces` | Divergent piece | An exact piece and a piece with a SNP and a 3 base deletion: two hits, the second with its mismatch and deletion listed. |
| `align.rs::every_piece_of_a_patchwork_query_is_reported_including_the_inverted_one`, `pieces_on_one_diagonal_with_a_long_unmatched_middle_are_not_merged`, `max_hits_still_limits_the_number_of_pieces_and_keeps_the_best`, `a_piece_repeated_in_the_query_is_found_twice_on_the_same_template_place`, `an_unmatched_stretch_is_cut_out_but_a_short_one_is_still_bridged` | The same in Rust | Layouts far apart / reordered / close, line and circle; three pieces on one diagonal; `max_hits`; a duplicated piece; an 8 base unmatched stretch is bridged, 40 bases are not. |
| `test_clean_query_accepts_fasta_lowercase_rna_digits_and_gaps` | Query text cleaning | FASTA header, spaces, digits (ORIGIN block), `-`, lowercase, U, IUPAC codes, CRLF, several records (first one used). |
| `test_clean_query_says_what_is_wrong` | Query validation | Empty, header only, digits only, a protein, other letters, 200,001 bp: each raises a clear message; exactly 200,000 bp is accepted. |
| `test_structural_variant_query_gives_one_row_per_piece_in_query_order` | Local Align table for the patchwork | 4 rows ordered by query position, query and template ranges, strand words, identity, the "before this piece" texts (200 bp unmatched, 400 bp of template skipped, inversion breakpoint), the summary line, the purple inverted cell, a row click selects the template span, the track bars say forward / inverted. |
| `test_structural_variant_query_on_a_circular_template_and_typed_as_its_reverse_complement` | Same on a circle, backwards | The reverse-complemented query on a circular template: the same places, order of pieces reversed, strands flipped, summary counts 1 forward / 3 inverted. |
| `test_summary_and_junction_texts_name_insertions_deletions_inversions_and_duplications` | Event wording | Contiguous, insertion, deletion, template jumps back, inversion breakpoint, two reverse pieces, the short way round an origin, and the summary line. |
| `test_paste_from_clipboard_cleans_the_text_and_explains_when_there_is_no_dna` | Paste from clipboard | Empty / whitespace / prose / header-only / protein clipboards leave the query alone with a message; a FASTA with numbered lowercase lines becomes the 60 bases. |
| `test_query_box_input_lowercase_rna_protein_and_huge_text` | The query field | Lowercase RNA is found like DNA; a protein, an empty box and a 200,004 bp query each give a message and no hits (stale hits are cleared). |
| `test_align_clipboard_menu_entry_opens_the_dialog_with_the_clipboard_aligned` | Tools > Align Clipboard | The entry opens the dialog with the query loaded and 4 rows aligned; an empty clipboard opens it with the reason; both Local Align entries are in the menu. |
| `test_find_parts_forward_reverse_and_the_identity_thresholds` | Detecting a known part by DNA identity | A 100 bp part with flanks on both strands; 4 mismatches (0.96) is found at the 0.95 default, 6 (0.94) is not until the threshold drops to 0.90. |
| `test_find_parts_coding_flag_and_min_protein` | Finding a gene that has been codon-optimised | A 60-aa ORF recoded with synonymous codons falls below 0.9 DNA identity but is found by protein identity 1.0; 3 of 60 amino acids changed gives 57/60, which passes at 0.9 and fails at 0.97. |
| `annotate.rs::finds_parts_on_both_strands_across_the_origin_and_as_protein` | The same in Rust, including the origin | A promoter forward, one across the origin, a reverse-strand gene at 98 % identity, and a codon-recoded CDS found by protein. A linear copy gives no origin hit. |

## 9. Assembly planning

`tests/test_assembly.py` (17 tests), `app/assembly.py`, `app/overhangs.py`

| Test | What it checks | How |
|---|---|---|
| `test_fragment_is_a_pcr_product_and_flips_with_its_primers` | A fragment really is what its two primers would amplify | 600 bp template, both primers carrying 6 nt tails: the fragment must equal tail + core + reverse-complemented tail; flipping swaps which primer is forward, and flipping twice restores. |
| `test_gibson_ends_that_do_not_fit_are_reported_then_fixed_by_new_primers` | Gibson needs end homology, and the designer creates it | Three unrelated fragments give "no homology" and a failed circle; after designing primers every overlap is 20–40 bp with Tm ≥ 52 °C, 6 oligos, and the fragments close in order A→B→C→A. |
| `test_gibson_existing_overlap_and_uniqueness` | Homology that already exists is reused; dangerous homology is flagged | A 30 bp shared end needs no new primers; the same 30-mer occurring twice makes the junction ambiguous and fails; a 16 bp overlap still joins but warns (below the 20 bp minimum). |
| `test_every_junction_sits_where_it_says_in_the_product` | The bases a junction claims are the bases actually there | For all three methods, the joining sequence is located in the assembled product (wrapping the circle if need be); Golden Gate junctions show a 4-nt 5′ overhang on both sides. |
| `test_golden_gate_vector_and_inserts_become_one_circle` | A destination vector plus two inserts assembles | The vector is cut, never amplified; three distinct overhangs; no BsaI site left in the product; flipping an insert breaks the design until re-primed; a vector with the wrong number or orientation of sites is refused. |
| `test_golden_gate_reads_the_overhang_each_end_really_carries` | Overhangs are read off the digested ends, not assumed | Fragments built as pad + GGTCTC + spacer + overhang + core…; ends join only when the two overhangs agree, and the mismatching pair is named (TCGC vs TCGA). An internal BsaI site is reported. |
| `test_golden_gate_warns_about_poor_overhang_choices` | The four ways a Golden Gate overhang set goes wrong | Repeated overhang → fails; a palindrome (`GATC`, ligates to itself) → fails; two overhangs one base apart → warning; all-A/T → GC warning; an exact reverse-complement pair → fails. |
| `test_the_papers_overhang_sets` | The Potapov 2018 ligation-fidelity data | 25 unique 4-mers stored correctly; the published low-fidelity set is flagged; `GGAG`/`CTCC` (reverse complements) is "bad" and `GGAG`/`GGAA` (one base apart) is "warn"; the designer prefers published overhangs and never repeats one. |
| `test_sidewinder_and_the_simulated_product` | Oligo-based assembly of short parts | Three 100–120 bp parts give 2 junctions and a 310 bp product carrying each part's features; fragments too long to order as oligos are refused ("at most 165 nt"), as is a single fragment. |
| `test_product_keeps_the_features_of_every_method` | Annotation survives assembly | Gibson and Golden Gate products are circular, the designed length, and carry the source promoter exactly once. |
| `test_gibson_with_new_primers_keeps_every_orf_on_its_bases` | An ORF keeps all of its bases in a Gibson product | Three fragments, each one ORF (ATG, sense codons, TAA) with its CDS over all of it, get new primers for overlaps of 20 to 25 bp (so the part of the next ORF that sits in the neighbour's tail is 10 to 13 bases, every frame offset). In the product each CDS must read exactly the designed bases, keep `codon_start 1` and translate from M to a single final stop. Was a bug: the ORF of fragments 2 and 3 started 10 to 13 bases late. |
| `test_a_circular_gibson_keeps_every_orf_including_the_one_that_closes_the_circle` | The same when the circle is closed | The last fragment's end is cut where it repeats the start of fragment 1; its CDS must still end on its stop codon, wrapping round the origin as a `join`. |
| `test_golden_gate_keeps_every_orf_including_the_bases_in_an_overhang` | An ORF keeps the bases the next fragment's overhang carries | A destination vector and three ORFs, 6 seeds: the overhang is picked where it scores best and takes 0 to 4 bases off the end of an ORF (its stop codon); every CDS must still read its designed bases. |
| `test_gibson_keeps_an_orf_that_starts_inside_an_overlap_the_fragments_already_share` | Overlaps that already exist | Each fragment ends with the first 30 bases of the next ORF, unannotated; every CDS must still start on its ATG in the product. |
| `test_an_orf_on_the_minus_strand_keeps_its_bases_too` | The other side of every junction | An ORF on the minus strand in the middle (its stop at the fragment's left end, its start at the right end), Gibson linear and circular and Golden Gate. |
| `test_sidewinder_keeps_every_orf` | Sidewinder places pieces whole | Three short ORFs; every CDS reads its designed bases. |
| `test_a_feature_inside_an_existing_overlap_is_annotated_once` | Two fragments that both carry an annotation of their overlap | Two stretches of one plasmid overlapping by 30 bp, with a feature inside the overlap: the product has it once, and the features outside the overlap are where they were. |
| `test_designed_primers_are_temporary_until_saved_and_explain_themselves_on_hover` | Designed primers are provisional | Drawn as temporary primers with their tails and a tooltip giving Tm and hairpin; they become permanent only on save, and move to the opposite strands when a fragment is flipped. |
| 18 further tests | The Assembly tab itself | Fragment slots, the junction cards and their chips, the read-only sequence view, defining a fragment from a selection or from two database primers, saving fragments and primers, the product preview, right-click actions, and the tab surviving file loads. |

**Checked and correct:** the four Type IIS specifications (BsaI GGTCTC(1/5), BsmBI CGTCTC(1/5), BbsI GAAGAC(2/6),
SapI GCTCTTC(1/4)) and their overhang lengths; the bottom-strand cut arithmetic as the exact mirror of the top strand;
the 5′-overhang topology (the left fragment's bottom strand carries it, the right fragment's top strand does) and the
rule that each piece contributes only its left overhang, so the junction is scarless; Gibson's split overlap
(20–40 bp, Tm ≥ 52 °C); the overhang rules (a palindrome self-ligates, a reverse-complement pair cross-ligates); and
the DOIs, journals, volumes and pages in `app/citations.py`.

## 10. Enzyme picker, gel, protein and 3D structure

| Test | What it checks | How |
|---|---|---|
| `test_enzymes.py::test_degenerate_sites_and_two_sided_cutters_are_found_where_they_cut` | Wobble letters, enzymes cutting outside and on both sides of their site, bottom-strand sites, origin-crossing sites | AccI on three variants, BglI, BcgI, BsaI both orientations, circular EcoRI; exact top/bottom coordinates. |
| `test_enzymes.py::test_the_common_list_still_cuts_where_it_always_did` | The enzyme table does not drift | EcoRI/KpnI/EcoRV/BsaI/SapI offsets, and the one-line description names the right kind of end (`G^AATTC, 5′ overhang 4 nt`; a Type IIS enzyme as `GGTCTC(1/5)`). |
| `test_enzymes.py::test_the_digest_filters_enzymes_by_name_site_and_number_of_sites` | The Digest picker | Filtering by name or by recognition sequence, "unique cutters only", and a ticked enzyme is remembered when the filter hides it. |
| `test_enzymes.py::test_the_gel_has_five_ladders...`, `test_digest_lanes_are_side_by_side_digests...` | The virtual gel | Five size ladders with their smallest/largest bands and reference band; several digests compared side by side, each with its own enzymes, fragment table and lane; fragment sizes sum to the full length. |
| `test_protein.py::test_mass_pi_charge` | Protein properties | Gly-Gly weighs 132.1191 Da (two glycine residues plus water); poly-lysine has pI > 10 and positive charge at pH 7, poly-aspartate the opposite; translation stops at the first stop codon. |
| `test_protein.py::test_follows_the_selection_and_pins_a_cds`, `test_selecting_residues_selects_codons` | The Protein pane follows the DNA | It shows a selected CDS, a feature inside it, or a bare stretch in a chosen reading frame, and keeps up with edits; dragging across two amino acids selects exactly their six bases. |
| `test_structure.py` (9 tests) | 3D structure lookup | AlphaFold DB by accession, the PDB by sequence search, and ESMFold prediction, each with the network faked; refusals for a too-short, too-long or non-protein sequence; results are cached; "Auto" tries AlphaFold first only when a UniProt accession is known and reports every source that failed. A reply pointing at a local file is refused, and text smuggled into a downloaded PDB cannot run as code in the viewer. |

| `test_genetic_code.py::test_every_codon_of_every_table_matches_biopython_on_both_strands` (one per table, 24 tables) | The embedded NCBI translation tables | All 64 codons of each table against Biopython's `unambiguous_dna_by_id` (used by the test only, never by the app): the Rust `translate`, the Python codon map, and a CDS read backwards on the bottom strand. Tables 27, 28 and 31 (a codon is a stop or an amino acid depending on its place) are not offered. |
| `test_genetic_code.py::test_start_and_stop_sets_of_every_table`, `test_tables_differ_where_they_should` | Start and stop codons per table; the well-known differences | Both sets equal Biopython's; table 4 reads TGA as Trp, table 2 reads AGA/AGG as stops and ATA as Met, table 12 CTG as Ser, table 25 TGA as Gly; table 11 = table 1 amino acids with more start codons; an unknown id falls back to the standard code. |
| `test_genetic_code.py::test_orfs_start_at_alternative_codons_only_when_allowed`, `test_orf_stops_follow_the_table` | ORF finder with GTG / TTG starts and other tables | A GTG ORF is found only with table 11 and alternative starts on (never by default, never in table 1); TTG in table 1; both strands; table 4 does not stop at TGA; table 2 stops at AGA. Rust unit tests do the same in `bio.rs` / `code.rs`. |
| `test_genetic_code.py::test_a_cds_uses_its_own_table_regardless_of_the_setting`, `test_codon_start_and_start_codon_as_methionine`, `test_transl_except_is_applied_and_an_unreproducible_one_is_flagged` | Per-CDS translation | `/transl_table` beats the setting; `/codon_start` 2 shifts the frame; a TTG start is Met unless the 5' end is partial; selenocysteine `/transl_except` is applied, and one that cannot be reproduced is flagged in the tooltip instead of mistranslated silently. |
| `test_genetic_code.py::test_the_setting_persists...`, `test_settings_dialog_shows_stops_and_starts...`, `test_views_and_exports_use_the_chosen_table` | The Genetic code setting | Persists in `settings.ini`; the dialog line lists the table's stop and start codons; the Protein pane title ("Table 4"), the sequence view's CDS lane, the FASTA export and `exports.protein_of` all follow the table. |
| `test_genetic_code.py::test_refresh_updates_stale_translations_and_leaves_the_rest_alone`, `test_refresh_matches_the_translation_drawn_in_the_views`, `test_refresh_uses_the_settings_table...` | `/translation` refresh | Stale translations are recomputed (own table, `/codon_start`, selenocysteine, bottom strand); `/pseudo`, partial-without-frame and unreproducible CDSs and CDSs with no `/translation` are untouched; one undo step; a second run adds none; the result equals what the views draw. |
| `test_genetic_code.py::test_features_menu_toggle_is_persistent...`, `test_save_refreshes_translations_but_autosave_never_does`, `test_save_leaves_translations_alone_when_the_option_is_off`, `test_update_now_acts_on_the_selected_cds_else_all...` | Feature > Keep /translation up to Date | On by default and remembered; save refreshes and reports in the status bar; autosave leaves the document, its version and its history untouched; saving twice adds no step and leaves it unmodified; off means off; "Update now" works on the selected CDS or all, as one undo step. |

**Checked and correct:** all 20 average residue masses and the water mass against ExPASy; the pKa set against EMBOSS
`Epk.dat` (the code credits EMBOSS, not ProtParam, so a pI differing from ProtParam's Bjellqvist values is expected,
not an error); the IUPAC complement table including B↔V and D↔H; the standard genetic code string (ATA = Ile, i.e.
NCBI table 1, not the mitochondrial code, and the other NCBI tables the program offers, codon by codon); AlphaFold's pLDDT bands.

## 11. The application itself

Selected biology-bearing GUI tests (`tests/test_gui.py`, `tests/test_app_logic.py`, `tests/test_import_export.py`).

| Test | What it checks | How |
|---|---|---|
| `test_design_pcr_orfs_digest` | Primer design, PCR and ORFs end to end | Designed primers sit outside the target with Tm near 60 °C; the product equals the template between them; a mismatched primer mutates the product; the plasmid's own CDSs are found as ORFs. |
| `test_clipboard_typing_and_primers` | Copy/paste carries annotation; primers behave | Features travel with copied bases, including between windows; pasted FASTA is cleaned; only ACGTN can be typed; a tailed primer and a 1-mismatch primer are placed correctly and move with edits; only exact sites are written on save. |
| `test_display_filter_translation_menus_zoom_shortcuts` | Translation display and protein motif search | A CDS's codons start at its start codon; a protein motif maps back to the right bases; the 1-based range box ("101..200"), including through the origin. |
| `test_gene_groups`, `test_app_logic.py::test_gene_groups_still_cover_the_same_bases_after_edits` | Grouping a promoter + RBS + CDS as one gene | The gene spans exactly its parts, folds to one arrow on the map, grows when a part moves, and survives inserts and deletes; saved as a shared `/locus_tag`. |
| `test_detect_common_features` | The Features Database | Library parts found on both strands, and a codon-recoded sfGFP found through its protein; the minimum-length setting drops a 30 bp promoter at a threshold of 31. |
| `test_db_tabs.py` (6 tests), `test_parts_library.py::test_local_features_are_a_source_and_can_be_removed` | The Databases menu, filter by values, the reference box | The Features Database is a tab of the window: the Features view (Source, the source's name and the record or file in one field, instead of Location) and, on the right, the feature clicked alone (its bases, the feature over all of them; a hint before any click); filter by Source (by the source's name) shows only Local; the filter and the feature on view survive the database changing; only Local features can be deleted (the one on view goes, the hint comes back) and the file keeps its header; the tab follows the window to another file and reopens as it was; Help shows the database's status. The Primer Database tab lists every primer (one that binds nowhere too) with its note and date, shows the primer clicked alone, sorts by a header click (Length as a number, a second click the other way), and keeps the sort and the primer picked after a delete. Filter by values on a table and a tree: a refill keeps the filter, an unseen value is shown, the header is marked, the search box narrows what All / None tick. The toolbar's reference box lists the last three genomes, newest first, and can switch the check off.. A primer's sites are grouped (on this sequence; in the reference genome, searched once the primer is opened, on a worker thread); the primer's bases that do not pair with the template are marked for drawing in red, on the sequence and in the genome (a planted whole-primer match on the minus strand is found exact); ambiguity codes pair with what they stand for. |
| `test_import_export.py::test_text_exports` | Export formats | GFF3 header and strand column; one BED line per feature part; the NCBI feature table writes minus-strand features as start > stop and drops program-specific qualifiers; a CDS on each strand translates to `MKPG`; FASTA wraps at 70. |
| `test_import_export.py::test_dna_lock` | The DNA lock | Typing, backspace, paste and undo of a base edit are all refused with a message while locked; adding a feature still works. |
| `test_primer_origin.py` (6 tests) | Primers: in this file or only in the database | Roles per origin (database only and designed: temporary style; file only; in both: the file's style wins); pixels of the sequence view in the two colours and the filter switching them; Import to This Sequence is one undo step, redo, tooltip text, the site becomes regular; survives Save and reload as an on-sequence primer, a database-only one is not written; menu entries per primer kind and the bulk helper; a primer made from this sequence belongs to the file; the Display panel's two swatches carry the two colours and tooltips. |
| `test_app_logic.py::test_each_local_align_hit_is_its_own_row_in_the_alignment_view` | One row per alignment | Four hits (two direct copies, the original, one inverted) get four rows with the heights of the track, labels with the inversion mark; clicking a label picks the hit and the table row, picking a table row picks the view's row; condensed and too-short tracks paint. |
| `test_import_export.py::test_hide_a_primer_in_this_file` | Hiding a primer in one file | The site stops being drawn but is still found, still in the database, listed as "(hidden)", and the choice is remembered for that file. |
| `test_import_export.py::test_features_without_direction` | Features with no arrow | `misc_feature` has no direction by default, a CDS does, and either can be changed; only a choice that differs from the default is written to the file. |
| `test_import_export.py::test_tab_image_export` | Export → Image of the Current Tab | PNG, TIFF and SVG of the shown tab: the raster files are twice the tab's size at 2× with 300 DPI written in them and real content; the SVG draws the views as vector paths (the map's cached screen pixmap is not embedded), its text is at screen size inside the scaled group (not resized again for the DPI), it prints at pixels / DPI, and every text is filled `nonzero` (glyphs whose outlines overlap stay solid). |
| `test_map_view.py` (4 tests) | A crowded plasmid map; the position scale | 300 small features on a 6 kb circle and line in a 500 px map: the largest are kept (every one drawn at least as long as every one left out), lanes stay at least 5 px thick, the note counts the rest, a bigger map holds more, a map of 6 features draws them all with the editor's own lanes, and a click over a left-out feature's place finds no feature. The scale, from 100 bp to 10 Mb on a 300 and a 1400 px backbone: steps of 1, 2 or 5 x 10^k, 2 to 40 of them, labels narrower than their spacing, minor ticks dividing the step and 4 px or more apart, a bigger map never coarser; labels 750, 2,500, 15 kb, 2.5 Mb, and 1 for the first base. |
| `test_checkup.py` (23 tests) | Plasmid Checkup: will each CDS be expressed? | Constructs built from a consensus sigma-70 promoter, the rrnB T1 terminator and a strong RBS: an annotated unit expresses its CDS; with nothing annotated the promoter, terminator and RBS are predicted on either strand ("Likely"); a terminator before the CDS blocks it and one inside truncates it; a CDS on the wrong strand or with no RBS; an early stop codon; a lac promoter noted as IPTG-inducible with lacI on the plasmid; CMV (eukaryotic: strong Kozak and polyA expressed, weak Kozak without polyA unlikely) and U6 (Pol III, not translated); a shuttle vector where the bacterial unit reads through a polyA signal and the CMV unit through rrnB T1, each CDS judged per host; a single promoter with no terminator on a circle and on a line; replication origins named (pUC, p15A, ColE1, R6K needing pir, f1 as a phage origin), two ColE1-type origins and a plasmid with no bacterial origin flagged; the panel's own map and sequence showing the units as temporary features only; a unit across the origin added as an mRNA feature; few chance predictions in 50 kb of random sequence; pBAD30 with eilR expressed from araBAD. Units as features: thin over the promoter and thick from +1 on each strand, a unit running all round a circle split at the origin; how units end (intrinsic when a hairpin + T tract is found, rho-dependent? when not, polyA, runs on); the filter's own switch for units (not listed as mRNA, independent of the mRNA type, remembered by projects, the program's qualifier left out of GFF3); both views painting units on a circle and a line; the panel's tree (each unit with the CDSs it carries, then those no unit carries), its results drawn in the document's own views only while it is shown, under the document's filter, selectable and hoverable but not editable, the document's added units not drawn twice; Tools > Plasmid Checkup opening one tab of panes (checkup, map, sequence) and showing it again rather than a second; +1 known to the base only for an annotated TSS and a T7 promoter, estimated for an ordinary promoter; motifs near the start of a linear sequence placed where they are; R6K with pir on the plasmid; the flat one-colour icons. |
| `test_app_logic.py` primer-check group | The traffic-light primer checks | Tm, GC %, 3′ GC clamp, hairpin and dimer verdicts are probed from both sides of every threshold, and the thresholds are then moved to confirm the edges follow. |
| `test_view_help.py::test_every_picker_view_has_a_guide_and_an_options_entry` | The "?" guide texts | Every view of the picker has a guide of 4-8 lines and at most 130 words, and an options entry. |
| `test_view_help.py::test_options_targets_exist` | The cog's targets | Each "settings" target is a real tab of Tools > Settings. |
| `test_view_help.py::test_cog_and_help_buttons_in_every_pane` | The cog | In a real pane, each view's cog is disabled with the tooltip "No options for this view", or opens the right settings tab / the Primers options. |
| `test_view_help.py::test_help_popup_opens` | The "?" button | Clicking it opens the dialog with that view's text; both buttons are reachable by Tab. |
| `test_input_validation.py` (30 cases) | Bad input in every audited field (`docs/INPUT_AUDIT.md`) | A pasted FASTA / numbered block is cleaned in DNA and protein fields and a non-sequence refused; primer name / sequence / note limits and messages, a long selection not pre-filled; Tm can never be -273 °C (Na, primer minimum); a missing genome file refused; feature location blank / nonsense / out of range and bad qualifier keys give messages and add nothing; the locus box survives `,`, `-`, `0`, a 20-digit number and positions past the end; File > New keeps the dialog open on bad bases; search box capped; min / max spin boxes cannot cross; tab names one short line; a giant raster export disabled; save-pair names checked. |
| `test_gui.py::test_genome_is_responsive` | A real genome stays usable | The 4.6 Mb *E. coli* genome scrolls and edits, and primer sites are matched in the background. |

`tests/test_feedback.py`, `tests/test_file_assoc.py` and `tests/test_updater.py` (102 tests) cover no biology: they
test the error log and bug-report link (no sequence, file name or user name ever leaves the machine), Windows file
associations, and the update checker (downloads verified by SHA-256, refused from any host but this repository's,
rolled back if the new program does not start).

---

## Biological accuracy review

Every expected value that encodes biology was re-derived independently. **No test asserts a biologically wrong
result.** The reference values that were re-checked and confirmed: the SantaLucia 1998 thermodynamic table and salt
corrections, the von Ahsen Mg²⁺ conversion, the hairpin loop table, all 47 enzyme recognition sequences and cut
offsets, the standard genetic code and stop set, every other embedded NCBI table (against Biopython's copy), the GenBank/GFF/BED/VCF/SAM coordinate conversions, the Cas9 and
Cas12a cut geometries, the ExPASy residue masses, *E. coli* BW25113 at 4,631,469 bp, and the sizes of the features
used in the plasmid tests.

What the review did turn up, and what was done about it:

**Fixed in this pass (comments that said something false, and one assertion that could never fail):**

1. `tests/test_gui.py` — the test for "a gene whose name is spelled with DNA letters is still found by name" searched
   for `SsrA`, which does not exist in the test plasmid, so the assertion passed no matter what the code did. Replaced
   with a real check: a CDS named `cat` is looked up by name and wins over the DNA that spells it.
2. `tests/test_core_functions.py` — the comment above the origin-crossing ORF test wrote the ORF as five codons; it is
   four (`ATG AAA | AAA TAA`). The assertion was right.
3. `tests/test_core_functions.py` — "all 46 common enzymes"; there are 47.
4. `core/src/search.rs` — a leftover comment described a motif (`ACCG`/`CGGT`) that is not in the sequence being tested.
5. `core/src/search.rs` — a comment said the protein motif starts "at 4 (frame 2)"; it starts at base 3, frame 1, which
   is what the assertion says.
6. `core/src/primers.rs` — a hairpin was described as a `GGGGCCCC` loop; it is a `GCGCGC` stem with a `TTTT` loop.
7. `core/src/document.rs` — a comment said a feature lost its first 5 bases; it loses 1. The assertion was right.

**Two user-facing defects found by the review and fixed:**

8. **A Type IIS enzyme's description hid where it cuts.** `BsaI` was described as `GGTCTC: GGTCTC, 5′ overhang 4 nt` —
   identical in form to a palindromic cutter, with no indication that the cut lies *downstream* of the recognition
   site, which is the entire reason these enzymes are used for Golden Gate. Now written the way suppliers do:
   `GGTCTC: GGTCTC(1/5), 5′ overhang 4 nt`. The stored cut offsets were right all along; only the human-readable line
   (used in the Digest tooltip) was lossy. A test now pins the notation for BsaI, BbsI and SapI.
9. **The Protein pane counted the stop codon as an amino acid.** The header read "5 aa" while the properties line
   directly beneath it read "4 aa" for the same chain (`MKWH*`). The header now agrees: "4 aa + stop".

**Open questions worth your decision (no test or code changed):**

- **`source` does not grow when bases are added at either end.** Confirmed by direct test: on a 200 bp record with
  `source 1..200`, inserting 3 bases at the start gives `source 4..203` and inserting at the end leaves `source 1..200`
  while the sequence is 203 bp. INSDC requires `source` to span the whole record, so a file saved after such an edit is
  technically malformed. This is in the program, not the tests, and happens on linear and circular sequences alike.
- **Inserting exactly at the origin of a plasmid.** On a circle, position 0 and position *n* are the same point. If a
  feature spans the origin (`join(191..200,1..10)`), an insert there is strictly inside it, so it should grow — the way
  an insert one base later does. Instead the new bases fall outside the feature. The current behaviour is pinned by
  tests in `document.rs` and `app/selftest.py`, so changing it means changing those expectations too.
- **Degenerate primers are treated as mismatches.** An IUPAC letter in a primer (say Y where the template has T) counts
  as a mismatch, and is rejected outright in the 3′ triplet. A degenerate primer is synthesised as a pool of oligos, so
  the matching variant anneals perfectly. Primers designed by back-translating a protein will therefore be reported
  with false mismatches, or as not binding at all. This is a deliberate design choice today, pinned by
  `tests/test_primer_binding.py`.
- **GFF3 export writes `.` as the phase of every CDS.** The specification requires 0, 1 or 2, so strict validators will
  reject the file. No test covers phase.
- **The Sidewinder toehold is not single-stranded, so neighbouring fragments have nothing to anneal through.**
  Confirmed by running the designer: each fragment's barcode strand is the reverse complement of its *whole* coding
  strand, so the 10 nt toehold is fully base-paired inside its own fragment, and the next coding strand begins after
  it. Two blunt duplex ends meet. A toehold only works when one duplex presents it single-stranded to a complementary
  single-stranded region on its neighbour. The panel note states the sharing as fact, though `design_sidewinder`'s
  docstring already warns that the layout is simplified and should be checked against the paper before ordering.
- **The Sidewinder barcode arm is attached at the wrong end.** The arm is written as a 3′ extension of the reverse
  complement of the coding strand, which places it at the fragment's *left* end — one junction earlier than the
  junction it is reported against. The result is that one arm dangles off the product's free end and the last junction
  gets no arm at all. Both Sidewinder points need the actual strand layout from the paper to resolve; the method is
  the one part of the program flagged as simplified.
- **The high-fidelity overhang set cannot fail its own orthogonality check.** `pair_risk` deliberately exempts two
  overhangs that are both in the published set from the "one base apart" warning (the paper assembled all 25 together),
  and the set contains five such pairs. The test that asserts the set passes every check therefore proves nothing about
  that rule. The exemption is intentional and documented; the test's value is overstated.

**The main untested areas, in rough order of how much they would matter at the bench:**

1. ~~Alternative start codons (GTG, TTG) and the bacterial genetic code (table 11)~~: now covered. Every embedded NCBI
   table is checked codon by codon against Biopython, and ORFs with alternative starts are tested (`test_genetic_code.py`,
   `bio.rs`, `code.rs`). Still untested: a `/transl_except` stop completed by poly-A (treated as "not reproducible"), and
   tables 27, 28 and 31, which are not offered.
2. Hairpin and dimer ΔG are only checked for sign and ordering, never against a reference value (Primer3's `ntthal`
   or a published ΔG), so a systematic error in the loop table would not be caught.
3. ~~Translation with `/codon_start` other than 1~~: now covered (sequence view codons, FASTA export, `/translation` refresh).
4. Find does not search the DNA of a plasmid across its origin (Local Align does), and no test covers it.
5. Digest fragment *sizes* are checked only by their sum, not individually.
