# `tests/test_enzymes_and_guides.py`

[← all test files](README.md) · source: [`tests/test_enzymes_and_guides.py`](../../tests/test_enzymes_and_guides.py)

## `test_cut_finder_the_common_list_and_the_filter`

**Result:** ✅ passed in 0.00 s · [line 57](../../tests/test_enzymes_and_guides.py#L57) · **Area:** Restriction site finding · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** Without it, the program could report cut positions at the wrong base (or miss sites on the bottom strand, across the origin, or with IUPAC codes), so every digest, map and cloning plan built on them would be wrong.

Checks cut positions for degenerate sites, 3' overhangs, enzymes cutting on both sides, bottom-strand sites and sites spanning a circle's origin; that the shipped common enzyme list has the right cut offsets and supplier-style descriptions; and that the enzyme filter only keeps an enzyme when every selected criterion holds.

| # | Input | Expected output |
|---|---|---|
| 1 | AccI (GTMKAC) on AAA GTAGAC TTT GTCTAC AAA GTGGAC | cuts [(5,7),(14,16)] - GTGGAC does not match |
| 2 | BglI (GCCNNNNNGGC) on TT+site+TT | [(9,6)] - top cut after the bottom cut, a 3' overhang |
| 3 | BcgI (cuts both sides) with 40 bp flanks | [(30,28),(64,62)] |
| 4 | BsaI site on top strand vs GAGACC on the bottom | [(9,13)] and [(15,19)], bottom cuts mirrored, 4 nt overhang |
| 5 | EcoRI on circular ATTCAAAAAAGA | [(11,3)] - site across the origin |
| 6 | the COMMON list | EcoRI cuts ((1,5),), KpnI ((5,1),), EcoRV blunt, BsaI ((7,11),), SapI ((8,11),); find_cuts on TTGAATTCAA gives EcoRI [(3,7)], EcoRV blunt cuts both strands at 5; find_sites gives {EcoRI:[3], BsaI:[5]} |
| 7 | describe() for EcoRI / KpnI / EcoRV / BsaI / BbsI / SapI / BcgI | 'GAATTC: G^AATTC, 5' overhang 4 nt'; KpnI ends '3' overhang 4 nt'; 'blunt' for EcoRV; 'GGTCTC: GGTCTC(1/5), 5' overhang 4 nt'; 'GAAGAC: GAAGAC(2/6),'; 'GCTCTTC: GCTCTTC(1/4), 5' overhang 3 nt'; BcgI mentions 'more place' |
| 8 | EnzymeFilter defaults and single criteria | default == EnzymeFilter(1,1,6) (unique cutters, 6+ bp): AluI rejected, EcoRI kept; palindromic keeps EcoRI not BsaI; type_iis keeps BsaI not EcoRI; ends 5/3/blunt select EcoRI/KpnI/EcoRV; min_size=8 keeps NotI not EcoRI |
| 9 | combined filter (palindromic + 5' + 6 bp) and a count filter (2-3 sites, <=4 bp) | EcoRI kept, KpnI and EcoRV rejected; AluI kept at 3 sites, rejected at 1; EcoRI rejected at 2 sites |

## `test_what_a_cut_says_about_itself`

**Result:** ✅ passed in 0.00 s · [line 119](../../tests/test_enzymes_and_guides.py#L119) · **Area:** Restriction site tooltips · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** Without it, the tooltip a user reads before ordering a digest could state the wrong cut base, the wrong overhang or the wrong strand, which is the only place those numbers are explained.

Checks the generated tooltip text for a 5' overhang, a 3' overhang, a blunt cutter and a Type IIS enzyme, for a site read on the bottom strand and for one crossing a circle's origin, and that every enzyme in the common list produces a tooltip at each of its cuts.

| # | Input | Expected output |
|---|---|---|
| 1 | EcoRI on TTGAATTCAA | names EcoRI, 'G^AATTC', 'bases 3-8', cuts after base 3 (top) / 7 (bottom), "5' overhang of 4 nt: AATT", labelled top strand, no 'Type IIS' |
| 2 | KpnI on TTGGTACCAA | 'GGTAC^C', cuts after base 7 top / 3 bottom, "3' overhang of 4 nt: GTAC" |
| 3 | EcoRV on TTGATATCAA | 'GAT^ATC', cuts after base 5 on both strands, 'Blunt ends', no overhang text |
| 4 | BsaI on TTGGTCTCAAGGAATTTT | 'GGTCTC(1/5)', overhang AGGA, cuts after base 9 / 13, Type IIS sentence naming 1 nt and 5 nt after the site |
| 5 | BsaI site written as GAGACC (bottom strand) | 'bases 13-18, on the bottom strand, where it reads GGTCTC; GAGACC on the top strand', 'Type IIS', 'before the site'; the cut numbers in the text match find_cuts (7,11) |
| 6 | EcoRI across the origin of circular ATTCAAAAAAGA | 'bases 11-4', cuts after base 11 top / 3 bottom, overhang AATT |
| 7 | every enzyme in COMMON, sites concatenated twice, circular | each cut yields a tooltip containing the enzyme name, 'Cuts after base' and 'Recognition site' |

**Review notes:**

- Asserts on exact English sentences, so harmless rewording of the tooltip breaks the test.

## `test_guide_candidates_cuts_and_scores`

**Result:** ✅ passed in 0.00 s · [line 166](../../tests/test_enzymes_and_guides.py#L166) · **Area:** CRISPR guide finding · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** Without it, the Spacer Finder could list guides with the wrong PAM, the wrong strand or a cut site in the wrong place, so a user would order a guide that cuts somewhere other than shown.

Checks SpCas9 guide candidates on both strands with their blunt cut 3 bp from the PAM, reading an offset stretch and a site across a circle's origin, SaCas9's NNGRRT PAM, Cas12a's 5' TTTV PAM with its staggered cut, and that the guide score is a checklist that ranks the best guides first.

| # | Input | Expected output |
|---|---|---|
| 1 | AT filler + GUIDE + TGG + filler, SpCas9 | one candidate: start 40, end 63, strand +1, guide = GUIDE, pam TGG, cut at 57 (3 bp from the PAM) |
| 2 | the reverse-complement of GUIDE+AGG in filler | start 40, end 63, strand -1, pam AGG, cut at 46; the slice equals rc(guide+pam) |
| 3 | a stretch read with offset=50 / a rotated circular sequence with wrap=True | start reported as 60 on the whole sequence; across the origin start = len-12 with end > len; the same sequence as a line yields no candidate |
| 4 | SaCas9 with PAM ATGAAT vs ATGCCT | 21 nt guide with pam ATGAAT; ATGCCT gives no plus-strand candidate (not NNGRRT) |
| 5 | Cas12a, pad + TTTA + 23 nt guide + pad | start 50, end 77, strand +1, pam TTTA, cut 72 on the PAM strand and 77 on the other (18 / 23 bases after the PAM) |
| 6 | Cas12a on the bottom strand (TTTC PAM) | start 50, end 77, strand -1, cuts 50 and 55; slice equals rc(pam+guide) |
| 7 | pad + TTTT + guide, Cas12a | one candidate at 51 - TTTT is no PAM, the TTTG after its first T is |
| 8 | cut_text for Cas12a and SpCas9 | "5' overhang" for Cas12a; '3 bp' for SpCas9 |
| 9 | score_guide on GUIDE / poly-T / poly-G / AT-rich / a balanced guide | (100, no notes); score < 70 with a TTTT note; <= 45 with GC and run notes; <= 60 with a GC note; 100 for the balanced one |
| 10 | a sequence holding a good and a poor guide | candidates returned in descending score, the good GUIDE first |

## `test_other_guide_sites_here_and_in_a_genome`

**Result:** ✅ passed in 0.01 s · [line 219](../../tests/test_enzymes_and_guides.py#L219) · **Area:** CRISPR off-target search · **Checked against:** independent implementation · **Strength:** 🟢 strong

**Why it matters:** Without it, the off-target counts and the good/warn/bad verdict a user relies on when choosing a guide could be wrong, or the fast Rust search could silently miss sites a careful search would find.

Checks off-target counting by mismatch class in the construct and in a reference genome, the verdict each situation gives, that seed mismatches (next to the PAM) are weighted, that a genome identical to the sequence is not scanned twice, and that the Rust search returns exactly what a brute-force Python search returns.

| # | Input | Expected output |
|---|---|---|
| 1 | sequence with the guide, a 1-mismatch site and a 4-mismatch site, scan with max 3 mismatches | here == [0,1,0,0], genome None, verdict 'warn' (the 4-mismatch site is beyond the limit) |
| 2 | reference genome holding a perfect copy and a 2-mismatch copy | genome == [1,0,1,0], genome_name 'ref', verdict 'warn' |
| 3 | reference genome with two perfect copies | genome[0] == 2, verdict 'bad' |
| 4 | a construct holding the guide twice | each spec counts 1 other perfect site here; verdict 'bad' |
| 5 | genome document identical to the sequence being edited | genome is None (not scanned twice); verdict for no spec is 'info' |
| 6 | Cas12a sequence with the guide twice | 2 candidates; scanning the first gives here == [1,0] (PAM counted before the protospacer) |
| 7 | seed_mismatches for bit masks 1, 1<<19, (1<<19)\|(1<<7)\|1, (1<<8)\|(1<<12) over 20 nt | 0, 1, 1, 2 - only the bases next to the PAM count |
| 8 | 4000 random bases, 3 guides, core.find_guide_sites with NGG and 2 mismatches | the set of (position, strand) hits equals a plain Python scan of both strands |

## `test_a_saved_guide_site_is_read_back_and_round_trips_in_the_file`

**Result:** ✅ passed in 0.01 s · [line 284](../../tests/test_enzymes_and_guides.py#L284) · **Area:** Guide site storage and display · **Checked against:** round-trip/self-consistency · **Strength:** 🟢 strong

**Why it matters:** Without it, a saved guide could come back from a GenBank file with the protospacer, PAM or cut drawn in the wrong place, be unreadable in other software, or break for files written under the program's former name.

Checks the layout and tooltip of a saved guide feature on both strands, for Cas12a's two cuts and for a site across the origin; that the custom qualifier parses back (and nonsense is ignored); that saving and reloading is byte-identical and readable by Biopython as a plain misc_feature; and that old OpenGenEdit-keyed files still work while exports omit the private qualifier.

| # | Input | Expected output |
|---|---|---|
| 1 | guide feature on the + and - strand, SpCas9 | tooltip shows 5'-GUIDE-3', 'PAM TGG (NGG)', 'SpCas9', 'blunt cut', '3 bp from the PAM' and the cut base (lead+17 or lead+6) on both strands; layout.proto/pam swap sides with the strand; cut_places agrees |
| 2 | Cas12a guide feature TTTA + 23 nt | tooltip names Cas12a, 'PAM TTTA (TTTV)', "5 nt 5' overhang", cuts after lead+22 top and lead+27 bottom |
| 3 | guide feature written as join(n-9..n,1..13) on a circle | protospacer cells run from n-10 across the origin to 9; the guide appears once in the tooltip |
| 4 | cas9_qualifier / parse_cas9 for every nuclease, and for 'nonsense', 'name=X;pam=NGG', '' | round-trips to the same nuclease; damaged values give None |
| 5 | candidate cut from cas9.candidates vs cs.cas9_layout, both strands | layout.gaps == [candidate.cut] |
| 6 | save -> load -> save again | identical bytes; layout gap at lead+17; protospacer reads GUIDE and PAM reads TGG from the reloaded document |
| 7 | the saved file read with Biopython SeqIO | one misc_feature with label 'gRNA GACTCGTA', the cas9 qualifier, the guide spelled out in /note, start = lead, strand +1; the file contains the new key, no 'OpenGenEdit', and no line over 80 characters |
| 8 | a minus-strand guide saved and reloaded | Biopython reads strand -1 and start = lead; layout PAM at (lead, lead+3) |
| 9 | the same file with /OpenDNA_cas9 renamed to /OpenGenEdit_cas9 | qualifier kept as written, layout still valid with the same cut; neither key appears in the .tbl or .gff3 exports |
| 10 | exports of the new file | CAS9_KEY absent but the plain note (with the guide sequence) kept |

**Review notes:**

- Six separate concerns (layout, tooltip, parsing, file bytes, legacy key, exports) in one function.
