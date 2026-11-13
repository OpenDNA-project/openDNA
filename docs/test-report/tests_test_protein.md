# `tests/test_protein.py`

[← all test files](README.md) · source: [`tests/test_protein.py`](../../tests/test_protein.py)

## `test_mass_pi_charge_and_numbered_layout`

**Result:** ✅ passed in 0.01 s · [line 16](../../tests/test_protein.py#L16) · **Area:** Protein properties · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** Without it the protein panel could show a wrong molecular weight, pI or charge, or clicking in the numbered protein text could select the wrong residue.

Checks the mass of Gly-Gly against a hand value, pI and charge only in direction for poly-K and poly-D, truncation at the stop, the numbered layout of a 120-residue protein, and the residue index <-> text offset round trip.

| # | Input | Expected output |
|---|---|---|
| 1 | mass('GG') | 132.1191 +/- 0.001 (2 x 57.0519 + water) |
| 2 | isoelectric_point('KKKKKK') / ('DDDDDD') | > 10 / < 4 |
| 3 | charge at pH 7 of KKKKKK / DDDDDD | > 5 / < 0 |
| 4 | until_stop('MK*AA') | 'MK' |
| 5 | layout('A'*120) | 3 lines; first starts '     1  AAAAAAAAAA AAAAAAAAAA'; second starts with 51 |
| 6 | offset_of(i) then index_at for i in 0, 9, 10, 49, 50, 119 | Same index back |

**Review notes:**

- pI and charge are checked only for sign/range on homopolymers: a wrong pKa table or a broken pI search could still pass. An exact pI and mass for a known protein (e.g. against Biopython ProtParam, which is installed in .venv) would make it strong.
- Mass is checked for one dipeptide only; a wrong residue mass for any amino acid other than G would go unnoticed.
- The index/offset check is a round trip only; offset_of and index_at that were wrong in the same way would pass.
- 13 Nov 2026: the weaknesses above were fixed (see Changed by the audit); rated strong since.

**Fixed by the audit of November 2026:**

- 13 Nov 2026: gap closed. The proposed checks were added, pass on the current code and fail under the mutation that slipped through before. Was: Mass pinned only for G: 19 residue masses unchecked. \| pI and charge only by sign/range on homopolymers: any pKa table error passes. \| offset_of/index_at only as a round trip: both wrong the same way passes.
