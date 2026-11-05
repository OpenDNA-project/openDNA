# `core/src/annotate.rs`

[← all test files](README.md) · source: [`core/src/annotate.rs`](../../core/src/annotate.rs)

## `finds_parts_on_both_strands_across_the_origin_and_as_protein`

**Result:** ✅ passed · [line 201](../../core/src/annotate.rs#L201) · **Area:** Feature detection (known parts) · **Checked against:** hand-computed · **Strength:** 🟢 strong

'Detect common features' finds a part forward, a mutated part on the reverse strand with the right identity, a part split across the origin of a circle, and a codon-recoded CDS at the protein level.

| # | Input | Expected output |
|---|---|---|
| 1 | circular plasmid built from random_dna(3000, seed 1): 40-bp promoter at 500 (+25 shift), revcomp of a 300-bp gene with 6 substitutions at 1040, promoter split 15/25 across the origin; parts [promoter, gene], coding [false,true], min_dna 0.95, min_protein 0.9 | contains (part 0, 525, len 40, fwd, DNA) and (0, n-15, 40, fwd, DNA); gene hit start 1065, reverse, DNA level, identity 294/300 (within 1e-6) |
| 2 | same sequence as linear, promoter only | exactly 1 hit |
| 3 | 51-bp CDS 'ATGCTGCTGCTG...TAA' vs synonymous recoding 'ATGTTATTGCTT...TAG' embedded at 200 in 200+51+200 random bp | exactly 1 hit: start 200, len 51, forward, protein=true, identity 1.0 |

**Review notes:**

- The promoter checks use contains() (lines 228-229), so extra spurious hits in the circular case would not be detected.
- Protein-level search is only exercised on the forward strand; reverse-strand protein hits (coordinate mapping at line 164 of annotate.rs) are untested here.
