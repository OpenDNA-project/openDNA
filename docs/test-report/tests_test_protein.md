# `tests/test_protein.py`

[← all test files](README.md) · source: [`tests/test_protein.py`](../../tests/test_protein.py)

## `test_mass_pi_charge_and_layout`

**Result:** ✅ passed in 0.00 s · [line 30](../../tests/test_protein.py#L30) · **Area:** Protein properties · **Checked against:** hand-computed · **Strength:** 🟡 adequate

Protein mass, isoelectric point, charge, stop trimming and the numbered text layout give correct values, and residue index <-> text offset round-trips.

| # | Input | Expected output |
|---|---|---|
| 1 | mass('GG') | 132.1191 Da within 1e-3 (2 x 57.0519 + water) |
| 2 | isoelectric_point('KKKKKK'), isoelectric_point('DDDDDD') | > 10 and < 4 |
| 3 | charge('KKKKKK', pH 7), charge('DDDDDD', pH 7) | > 5 and < 0 |
| 4 | until_stop('MK*AA') | 'MK' |
| 5 | layout('A' * 120) | 3 lines; line 1 starts '     1  AAAAAAAAAA AAAAAAAAAA'; line 2 starts with '51' |
| 6 | offset_of(i) then index_at for i in 0, 9, 10, 49, 50, 119 | gives back i |

**Review notes:**

- pI and charge checked only by loose inequalities (> 10, < 4, > 5); a wrong pKa table would still pass (line 34-35)
- offset/index checked only as a round trip, the actual offsets are not asserted (line 41-43)

## `test_follows_the_selection_and_selecting_residues_selects_codons`

**Result:** ✅ passed in 0.01 s · [line 46](../../tests/test_protein.py#L46) · **Area:** GUI: protein view · **Checked against:** hand-computed · **Strength:** 🟡 adequate

The Protein view shows the translation of the selected CDS (or of a DNA selection in a chosen frame), highlights a sub-feature's residues, and selecting residues selects their codons in the editor.

| # | Input | Expected output |
|---|---|---|
| 1 | SEQ = 10 C + 'ATGAAATGGCATTAA' + 10 G; CDS 11..25 'tiny', misc_feature 17..22 'motif'; cursor at 0 | title contains 'Select a CDS' |
| 2 | select the CDS | text contains 'MKWH*', props '4 aa', title 'tiny', frame box disabled, one feature row 'aa 3–4' |
| 3 | select the motif feature | title still 'tiny'; highlighted residues == ['WH'] |
| 4 | select(0, 9) then frame index 3 (-1) | title 'Translation of 1–9', frame enabled, 'PPP'; then 'GGG' |
| 5 | select CDS, insert 'T' at 13, refresh(True) | view.src is None or its version == editor version |
| 6 | select residues K..W in the text (offset_of(1) to offset_of(2)+1) | editor selection (13, 19) (codons AAA TGG); title still 'tiny' |

**Review notes:**

- the edit section (line 70-76) says the changed protein is shown but only checks src is None or version matches; 'is None' makes it pass even if nothing is shown, and the new protein (M*) is never checked
- frame -1 check is only a 'GGG' substring (line 69)

## `test_a_pane_in_the_window`

**Result:** ✅ passed in 0.07 s · [line 92](../../tests/test_protein.py#L92) · **Area:** GUI: protein view · **Checked against:** hand-computed · **Strength:** 🟡 adequate

The Protein view docked in the main window labels a selected ORF with its residue count, not counting the stop, and paints without errors.

| # | Input | Expected output |
|---|---|---|
| 1 | MainWindow with SEQ; select(10, 25) = 'ATGAAATGGCATTAA' | title contains '4 aa + stop' |
| 2 | grab and close the window | no exception reached sys.excepthook |

**Review notes:**

- no docstring
