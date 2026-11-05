# `tests/test_primers.py`

[← all test files](README.md) · source: [`tests/test_primers.py`](../../tests/test_primers.py)

## `test_tm_matches_biopython`

**Result:** ✅ passed in 0.04 s · [line 21](../../tests/test_primers.py#L21) · **Area:** Primer Tm · **Checked against:** independent implementation · **Strength:** 🟢 strong

Checks that Tm with salt, Mg2+ and dNTP correction matches Biopython's nearest-neighbour model with Owczarzy 2008 correction, and that N gives no Tm.

| # | Input | Expected output |
|---|---|---|
| 1 | 200 random 12-40-mers (seed 1) under (Na 50, Mg 1.5, dNTP 0.6, 50 nM), (50,0,0,250), (100,2,0.8,500); palindromes skipped | equals Biopython Tm_NN(DNA_NN3, saltcorr=5, dnac1=dnac2=C/2) within 0.01 C |
| 2 | tm('ACGTNACGT') | None |

## `test_tm_of_self_complementary_primers_matches_biopython`

**Result:** ✅ passed in 0.00 s · [line 32](../../tests/test_primers.py#L32) · **Area:** Melting temperature · **Checked against:** independent implementation · **Strength:** 🟢 strong

Tm of primers that are their own reverse complement (they pair with themselves: whole concentration instead of C/4, plus the -1.4 e.u. symmetry term), against Biopython's selfcomp mode.

| # | Input | Expected output |
|---|---|---|
| 1 | GAATTCGAATTC, ACGTACGTACGT, GGATCCGGATCC, CGCGAATTCGCG, AAGCTTAAGCTT, GCGGCCGCGCGGCCGC at (Na, Mg, dNTP, primer nM) = (50, 1.5, 0.6, 50), (50, 0, 0, 250), (100, 2, 0.8, 500) | equal to Bio.SeqUtils.MeltingTemp.Tm_NN(DNA_NN3, saltcorr=5, selfcomp=True) within 0.01 °C (e.g. GAATTCGAATTC at the first conditions: 30.95 °C) |
| 2 | the same primers with the non-self-complementary formula | more than 1 °C away, so the branch is shown to matter |

**Fixed by the audit of November 2026:**

- New test added by the audit: test_tm_matches_biopython skips palindromes and the Rust tm test uses none, so the self-complementary branch of tm() was tested nowhere.

## `test_every_primer_the_dna_file_annotates_is_found_on_its_strand_with_the_same_3_prime_end`

**Result:** ✅ passed in 0.00 s · [line 47](../../tests/test_primers.py#L47) · **Area:** Primer binding · **Checked against:** independent implementation · **Strength:** 🟢 strong

Every primer site that SnapGene stored in the real .dna file is found by our binding search, on the same strand, inside the stored site and ending at the same 3' base. Together with test_primer_binding.py::test_our_sites_are_the_binding_sites_stored_in_each_real_dna_file (every site we find is a stored one) this checks both directions against an independent program.

| # | Input | Expected output |
|---|---|---|
| 1 | demo/pFG051_Chlor.dna (copy): its 15 primers and 15 primer_bind features | 15 stored sites and 15 primers (the comparison is never empty) |
| 2 | each stored site, e.g. oAC788 3385..3405 (+), oAC789 3116..3138 (-), oAC734 3764..3780 + 0..2 (-, across the origin) | binding_sites(primer, 0) has a site on that strand inside the stored bases, ending at the same 3' base |
| 3 | oAC813 stored 1233..1296 (-), oAC812 stored 356..410 (+) | found shorter: 1233..1259 and 393..410 (SnapGene also marks a 5' tail that matches by chance); exactly these two are shorter |

**Fixed by the audit of November 2026:**

- Replaced test_primers_bind_where_the_dna_file_says, which compared nothing: it read the GenBank plasmid, whose primer list is empty, so its loop never ran.
