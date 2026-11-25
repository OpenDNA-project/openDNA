# `tests/test_synthesis_salis.py`

[← all test files](README.md) · source: [`tests/test_synthesis_salis.py`](../../tests/test_synthesis_salis.py)

## `test_the_forest_gives_scikit_learns_answer_for_every_sequence_of_the_paper`

**Result:** ✅ passed in 1.07 s · [line 22](../../tests/test_synthesis_salis.py#L22) · **Area:** Synthesis Success Calculator (Salis) · **Checked against:** independent implementation · **Strength:** 🟢 strong

**Why it matters:** A mistake in reading the exported forest would make every synthesis verdict wrong without any error.

The 1512 exported trees give, for the nine features of each of the 303 paper sequences, the probability scikit-learn gets from the same trees, and the verdict.

| # | Input | Expected output |
|---|---|---|
| 1 | 303 feature rows from the paper's table (+ the repeat metric from the authors' code) | probability within 1e-12 of scikit-learn's, same True/False verdict, for all 303 |
| 2 | the paper's own Synthesized/Canceled labels | the verdict agrees with the label for more than 85 % of rows (it is 269 of 303, 89 %) |

**Review notes:**

- The agreement with the paper's labels is a sanity check only: the model has an F1 of 0.93 on held-out data and some of these fragments trained it.

## `test_melting_temperature_equals_biopythons_nearest_neighbour_value`

**Result:** ✅ passed in 0.05 s · [line 31](../../tests/test_synthesis_salis.py#L31) · **Area:** Synthesis Success Calculator (Salis) · **Checked against:** independent implementation · **Strength:** 🟢 strong

**Why it matters:** A wrong melting temperature changes how many windows count as too cold or too hot, and so the verdict.

The ported nearest-neighbour melting temperature (Tm_NN, dnac1=250, salt correction 7) equals Biopython on 20-mers.

| # | Input | Expected output |
|---|---|---|
| 1 | 2000 random 20-mers | difference below 1e-9 each |
| 2 | poly-A, poly-C, poly-G, poly-T, (AT)x10, (GC)x10 | same |

## `test_the_rust_repeat_finder_agrees_with_the_rules_written_the_plain_way`

**Result:** ✅ passed in 0.04 s · [line 70](../../tests/test_synthesis_salis.py#L70) · **Area:** Synthesis Success Calculator (Salis) · **Checked against:** in-test reference model · **Strength:** 🟢 strong

**Why it matters:** A repeat missed or invented by the fast repeat finder changes two of the nine features.

The Rust port of the authors' maximal-repeat finder returns the same (length, kind, positions) as the same rules written with str.find.

| # | Input | Expected output |
|---|---|---|
| 1 | 60 random sequences of 30 to 160 bases with 0 to 3 planted direct or inverted pieces of 9 to 16 bases | identical sets of (length, inverted, starts) |

**Review notes:**

- The reference model is the same rules transcribed by us, not the authors' program: the check against the authors' own numbers is the feature test below and scripts/salis_verify.py.

## `test_the_rule_based_hairpin_best_gc_is_a_stem_fraction`

**Result:** ✅ passed in 0.00 s · [line 83](../../tests/test_synthesis_salis.py#L83) · **Area:** Synthesis Success Calculator (Salis) · **Checked against:** hand-computed · **Strength:** 🟡 adequate

**Why it matters:** The rule-based hairpin pass feeds the richest_hairpin feature.

A GC-only 11 bp stem scores 1.0, poly-A scores 0, and a sequence shorter than any stem pair scores 0. salis_hairpin_gc also says where that hairpin is (its stem-loop-stem), for the Synthesis view to highlight.

| # | Input | Expected output |
|---|---|---|
| 1 | GCGCGCGCGCGAAACGCGCGCGCGC | (1.0, (0, 25)) |
| 2 | 80 A | (0.0, None) |
| 3 | ACGT | (0.0, None) |

**Review notes:**

- Three cases only; the full rule set is checked on the paper's sequences by the feature test and scripts/salis_verify.py.

## `test_a_window_starting_at_size_minus_20_is_never_counted_as_the_original_does`

**Result:** ✅ passed in 0.00 s · [line 89](../../tests/test_synthesis_salis.py#L89) · **Area:** Synthesis Success Calculator (Salis) · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** The original code never reads the last 20-base window; a port that 'fixed' this would no longer be the published model.

GC_short_l counts windows that start before size - 20 only. Also where the counts were made: the counted windows merged into stretches, and the copies of the longest repeat.

| # | Input | Expected output |
|---|---|---|
| 1 | 40 C then 20 A | 4 (the all-A window at 40 is skipped) |
| 2 | 40 A then 20 C | 25 |
| 3 | 40 C + 20 A, places requested | GC_short_l stretches [(36, 59)]; total_GC [] |
| 4 | GATTACAGGCTA + 20-base spacer + GATTACAGGCTA | longest_repeat 12 at [(0, 12), (32, 44)]; (0, 12) among the repeat_9_metric stretches |

## `test_clean_names_the_first_base_it_cannot_use`

**Result:** ✅ passed in 0.00 s · [line 103](../../tests/test_synthesis_salis.py#L103) · **Area:** Synthesis Success Calculator (Salis) · **Checked against:** hand-computed · **Strength:** 🟡 adequate

**Why it matters:** A sequence with an N or a stray character must be refused with its position, never silently changed.

Whitespace and case are accepted; the first other base is named; an empty text is refused.

| # | Input | Expected output |
|---|---|---|
| 1 | ' ac gt\nAC ' | ACGTAC |
| 2 | ACGNACG | error: Base 4 is 'N' |
| 3 | two spaces | error: no sequence |

## `test_the_features_equal_the_authors_numbers_on_the_shortest_sequences_of_the_paper`

**Result:** ✅ passed in 3.14 s · [line 112](../../tests/test_synthesis_salis.py#L112) · **Area:** Synthesis Success Calculator (Salis) · **Checked against:** independent implementation · **Strength:** 🟢 strong

**Why it matters:** The folding, window and repeat parts together must reproduce the published model on real sequences.

The nine features, the probability and the verdict of the 6 shortest paper sequences equal the authors' own table, the authors' repeat metric and scikit-learn's answer. The places returned with them are consistent: every feature named, inside the sequence, one richest-hairpin stretch exactly when its value is above 0, longest-repeat copies at least that long.

| # | Input | Expected output |
|---|---|---|
| 1 | the 6 shortest of the 303 sequences (under 300 bp) | all nine features equal (richest_hairpin to the table's 9 decimals); probability within 1e-12; places consistent with the values (property) |

**Review notes:**

- Only 6 of the 303 sequences run here (the whole set takes about 80 CPU-minutes): scripts/salis_verify.py runs all of them. Needs ViennaRNA (skipped without it; CI installs it).

## `test_a_fragment_shorter_than_one_folding_window_is_still_scored`

**Result:** ✅ passed in 0.01 s · [line 127](../../tests/test_synthesis_salis.py#L127) · **Area:** Synthesis Success Calculator (Salis) · **Checked against:** smoke · **Strength:** 🔴 weak

**Why it matters:** A short fragment must give an answer, not an error.

A 29-base fragment (no 100-base window to fold) gives size 29 and a probability between 0 and 1.

| # | Input | Expected output |
|---|---|---|
| 1 | 29 bases | size 29, probability in [0, 1] |

**Review notes:**

- Checks that it runs and the numbers are in range, not that they are right.

## `test_a_pasted_text_is_one_fasta_record_or_bare_bases`

**Result:** ✅ passed in 0.00 s · [line 132](../../tests/test_synthesis_salis.py#L132) · **Area:** Synthesis Success Calculator (Salis) · **Checked against:** hand-computed · **Strength:** 🟡 adequate

**Why it matters:** A pasted FASTA header, a second record or a stray digit must never reach the calculator as bases.

The Synthesis tab box accepts bare bases or exactly one FASTA record and names what it refuses.

| # | Input | Expected output |
|---|---|---|
| 1 | >pUC fragment / acgt / ACGT | ACGTACGT |
| 2 | two FASTA records | error: more than one FASTA record |
| 3 | ACGT1 | error: Base 5 is '1' |
| 4 | a header and no bases | error: no sequence |

## `test_the_viennarna_wheel_picked_fits_the_python_and_the_system`

**Result:** ✅ passed in 0.00 s · [line 149](../../tests/test_synthesis_salis.py#L149) · **Area:** Synthesis Success Calculator (Salis) · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** The first-use download of ViennaRNA could fetch a file that cannot load on the user's system, or a musl or source build.

From a PyPI file list, the wheel for this Python and system is chosen: Windows, macOS Intel and Apple Silicon, Linux (manylinux, not musl), and an error names what was looked for.

| # | Input | Expected output |
|---|---|---|
| 1 | cp312 win32 AMD64 | win_amd64 |
| 2 | cp313 darwin arm64 / x86_64 | macosx_11_0_arm64 / macosx_10_15_x86_64 |
| 3 | cp312 linux x86_64 / aarch64 | the manylinux file of that machine, never musllinux |
| 4 | cp39 | error: no ViennaRNA wheel for Python cp39 |
| 5 | win32 arm64 | error naming win32 (arm64) |

**Review notes:**

- Uses the file names of the real 2.8.0 listing, written out here; a new naming scheme on PyPI would not be noticed.

## `test_installing_viennarna_checks_the_checksum_before_unpacking_anything`

**Result:** ✅ passed in 0.01 s · [line 163](../../tests/test_synthesis_salis.py#L163) · **Area:** Synthesis Success Calculator (Salis) · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** A corrupted or swapped download could be unpacked into the data folder and imported.

The installer unpacks a wheel whose SHA-256 matches PyPI's listing and makes it importable, refuses one that does not (leaving no folder), and reports an unreachable PyPI as one message.

| # | Input | Expected output |
|---|---|---|
| 1 | a fake wheel with the right digest | version 9.9 returned; RNA/__init__.py unpacked |
| 2 | digest of 64 zeros | error: does not match the checksum; no folder left |
| 3 | PyPI unreachable | error: Could not download ViennaRNA from PyPI: Could not reach pypi.org |

**Review notes:**

- No real download is made (fetch is replaced), so the real PyPI answer format is only as good as the fake: the key names used are the documented ones.
