# `core/src/dna_file.rs`

[← all test files](README.md) · source: [`core/src/dna_file.rs`](../../core/src/dna_file.rs)

## `reads_synthetic_file`

**Result:** ✅ passed · [line 355](../../core/src/dna_file.rs#L355) · **Area:** SnapGene .dna import · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** Opening a SnapGene .dna file could lose or misplace its features, primers, topology or description.

A synthetic .dna file is parsed into sequence, topology, definition, features (with origin-spanning and direction handling) and primers as primer_bind features.

| # | Input | Expected output |
|---|---|---|
| 1 | packets: cookie, DNA '\x01AAAACCCCGGGGTTTTACGT', junk 0x07, features XML (lacZ 2-10 dir 1, ori 18-3 dir 2, box 5-8 no dir, tag misc 11-14 dir 1), primer p1 'ttttGGGG' site 9-12 boundStrand 1 (plus simplified duplicate), notes 'test plasmid' | seq 'AAAACCCCGGGGTTTTACGT', circular, definition 'test plasmid', 5 features |
| 2 | lacZ feature | location 1..10 (0-based), label 'lacZ', note 'a <b> c' (unescaped), codon_start '1', no 'directional' |
| 3 | ori | 'complement(join(18..20,1..3))' |
| 4 | box / tag | directional 'no' / 'yes' |
| 5 | primers | [Primer{p1, 'ttttGGGG'}]; feature 4 kind primer_bind at 'complement(10..13)' (sites 0-based) |

**Review notes:**

- The primer site (0-based '9-12' -> 10..13) is not where 'GGGG' would actually anneal in this synthetic sequence, so the test fixes the 0-based convention by assertion only; it cannot confirm it against real bases (the code comment says real files were checked).

## `says_what_it_could_not_read`

**Result:** ✅ passed · [line 384](../../core/src/dna_file.rs#L384) · **Area:** SnapGene .dna reader · **Checked against:** hand-computed · **Strength:** 🟡 adequate

**Why it matters:** Features or primer sites of a .dna file could be left out of the sequence view without a word.

A clean file has no notes; features without a position (self-closing, or only a gap segment), primer binding sites without a location and sites before any primer are counted in the notes, and the good ones are still read.

| # | Input | Expected output |
|---|---|---|
| 1 | the synthetic sample file | no notes |
| 2 | 3 features (one empty, one with only a gap segment, one good); 1 primer with a site without location and a good site | 2 features kept (the good feature and the good site); notes: '2 feature(s) ... not read.' and '1 primer binding site(s) ... not read.' |
| 3 | a primers packet holding a binding site before any <Primer> | note '1 primer binding site(s)...' |

**Review notes:**

- The GUI dialog that shows the notes is not tested (no widget tests).

## `rejects_other_files`

**Result:** ✅ passed · [line 406](../../core/src/dna_file.rs#L406) · **Area:** SnapGene .dna import · **Checked against:** hand-computed · **Strength:** 🟡 adequate

**Why it matters:** Opening a non-.dna or a damaged .dna file could give an empty document instead of an error.

Non-.dna data and a .dna file without a DNA packet are rejected instead of yielding an empty record.

| # | Input | Expected output |
|---|---|---|
| 1 | b'LOCUS       x' | error |
| 2 | cookie 'SnapGene\0' + notes packet, no DNA packet | error |

**Review notes:**

- Truncated packets (length past end of data) are not tested.

**Fixed by the audit of November 2026:**

- 13 Nov 2026: gap closed. The proposed checks were added, pass on the current code and fail under the mutation that slipped through before. Was: b'LOCUS       x' is rejected as a 'truncated packet' (its length bytes are huge), before the signature check is reached; the second input has a valid cookie and fails for 'no DNA'. So the 'not a .dna file' check is never exercised, and is_err() does not say which check fired.
