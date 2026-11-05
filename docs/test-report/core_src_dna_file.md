# `core/src/dna_file.rs`

[← all test files](README.md) · source: [`core/src/dna_file.rs`](../../core/src/dna_file.rs)

## `reads_synthetic_file`

**Result:** ✅ passed · [line 332](../../core/src/dna_file.rs#L332) · **Area:** SnapGene .dna import · **Checked against:** hand-computed · **Strength:** 🟢 strong

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

## `rejects_other_files`

**Result:** ✅ passed · [line 361](../../core/src/dna_file.rs#L361) · **Area:** SnapGene .dna import · **Checked against:** hand-computed · **Strength:** 🟡 adequate

Non-.dna data and a .dna file without a DNA packet are rejected instead of yielding an empty record.

| # | Input | Expected output |
|---|---|---|
| 1 | b'LOCUS       x' | error |
| 2 | cookie 'SnapGene\0' + notes packet, no DNA packet | error |

**Review notes:**

- Truncated packets (length past end of data) are not tested.
