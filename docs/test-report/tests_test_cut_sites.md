# `tests/test_cut_sites.py`

[← all test files](README.md) · source: [`tests/test_cut_sites.py`](../../tests/test_cut_sites.py)

## `test_enzyme_tooltips_by_kind_of_end`

**Result:** ✅ passed in 0.00 s · [line 50](../../tests/test_cut_sites.py#L50) · **Area:** Cut tooltips: end type · **Checked against:** hand-computed · **Strength:** 🟢 strong

The restriction-cut tooltip states the site, the two strand cut positions and the end it leaves for a 5' overhang, a 3' overhang, a blunt enzyme and a Type IIS enzyme (how far outside its site it cuts).

| # | Input | Expected output |
|---|---|---|
| 1 | EcoRI in 'TTGAATTCAA' | 'G^AATTC'; 'Cuts after base 3 on the top strand and after base 7 on the bottom strand.'; '5' overhang of 4 nt: AATT'; 'bases 3-8'; no 'Type IIS' |
| 2 | KpnI in 'TTGGTACCAA' | 'GGTAC^C'; cuts after base 7 top and base 3 bottom; '3' overhang of 4 nt: GTAC' |
| 3 | EcoRV in 'TTGATATCAA' | 'GAT^ATC'; cuts after base 5 on both strands; 'Blunt ends'; no 'overhang of' |
| 4 | BsaI in 'TTGGTCTCAAGGAATTTT' | 'GGTCTC(1/5)'; cuts after base 9 top and base 13 bottom; '5' overhang of 4 nt: AGGA'; 'Type IIS: ... 1 nt after the site on the top strand and 5 nt after ... on the bottom strand.' |

## `test_enzyme_tooltips_by_site_position`

**Result:** ✅ passed in 0.00 s · [line 78](../../tests/test_cut_sites.py#L78) · **Area:** Cut tooltips: site position · **Checked against:** hand-computed · **Strength:** 🟢 strong

The tooltip handles a site on the reverse strand, a site across a circle's origin, a palindromic site listed once, and that every common enzyme has a well-formed tooltip at each of its cuts.

| # | Input | Expected output |
|---|---|---|
| 1 | BsaI in 'TTTTAAAAGGAAGAGACCTT' (GAGACC = reverse-strand GGTCTC, bases 13-18) | 'bases 13-18, on the bottom strand, where it reads GGTCTC; GAGACC on the top strand'; 'Type IIS' and 'before the site'; cuts (7,11), 'after base 7 ... and after base 11 on the bottom' |
| 2 | EcoRI in 'TTGAATTCAA' (palindrome) | 'top strand)' (listed once, top strand) |
| 3 | EcoRI in 'ATTCAAAAAAGA' circular (site spans origin) | 'bases 11-4'; 'Cuts after base 11 on the top strand and after base 3 on the bottom strand.'; '5' overhang of 4 nt: AATT' |
| 4 | every COMMON enzyme in concat of their sites (N->A) x2, circular, each cut | each tooltip contains the name, 'Cuts after base' and 'Recognition site' |

**Review notes:**

- L96-101 the last section loops over every enzyme/cut but only checks substrings exist (name, 'Cuts after base', 'Recognition site'); a wrong coordinate for most enzymes would still pass.

## `test_cas9_and_cas12a_tooltips`

**Result:** ✅ passed in 0.00 s · [line 124](../../tests/test_cut_sites.py#L124) · **Area:** Cas9/Cas12a tooltip and layout · **Checked against:** hand-computed · **Strength:** 🟢 strong

The Cas9 guide tooltip and layout are right on both strands (protospacer, PAM, blunt cut 3 bp from the PAM) and Cas12a gives two staggered cuts with a 5' overhang; cut positions match cut_places.

| # | Input | Expected output |
|---|---|---|
| 1 | guide_doc(strand) SpCas9, strand 1 and -1 | '5'-GACTCGTACGTAGCATCGAT-3''; 'PAM TGG (NGG)'; 'SpCas9'; 'blunt cut', '3 bp from the PAM'; 'after base {cut} on both strands'; cut_places(SP,lead,strand)==(cut,cut); layout.gaps==[cut] |
| 2 | same, proto/pam placement | strand +: proto (lead,lead+20), pam (lead+20,lead+23); strand -: pam (lead,lead+3), proto (lead+3,lead+23) |
| 3 | Cas12a, site 'TTTA'+23 nt | 'Cas12a'; 'PAM TTTA (TTTV)'; 'after base {lead+22} on the top strand and after base {lead+27} on the bottom strand'; '5 nt 5' overhang' |

## `test_guide_site_origin_qualifier_and_spacer_finder`

**Result:** ✅ passed in 0.00 s · [line 154](../../tests/test_cut_sites.py#L154) · **Area:** Cas9 guide: origin, qualifier round-trip · **Checked against:** round-trip/self-consistency · **Strength:** 🟢 strong

A guide site across a circle's origin lays out correctly, the saved /cas9 qualifier round-trips (and a damaged one is ignored), and the Spacer Finder candidate and the view layout agree on the cut position on both strands.

| # | Input | Expected output |
|---|---|---|
| 1 | GUIDE+'TGG' split round the origin of a circular seq; join(start+1..n,1..13) feature | proto_cells[0]==start, proto_cells[-1]==9; tooltip contains GUIDE exactly once |
| 2 | cas9_qualifier(nuc) then parse_cas9 for every NUCLEASE; 'nonsense','name=X;pam=NGG',''  | round-trips to the same Nuclease; the three bad values parse to None |
| 3 | cas9.candidates on FILL*3+GUIDE+TGG+FILL, + and - strand | cas9_layout(...).gaps == [candidate.cut] on both strands |

## `test_hover_reaches_features_cut_marks_and_unique_cutters`

**Result:** ✅ passed in 0.02 s · [line 196](../../tests/test_cut_sites.py#L196) · **Area:** GUI: hover dispatch · **Checked against:** hand-computed · **Strength:** 🟢 strong

hover.describe returns the right text for a feature, a cut mark (callable or plain tip, missing mark safe) and a shown enzyme site; the sequence view locates the cut and the unique-cutter name under the pointer.

| # | Input | Expected output |
|---|---|---|
| 1 | guide_doc feature; set_cut_marks digest (callable tip) and 'x' (plain tip) and missing index | feature tip has 'SpCas9'; callable -> 'tip text'; plain -> 'plain tip'; missing -> '' |
| 2 | shown EcoRI site in 'TTGAATTCAA' | describe('site',name,pos) has 'G^AATTC' |
| 3 | 'TTGAATTCAA'*3+'GGATCC'+T*70, digest mark EcoRI, sequence view | cut_at at top/bottom/joining-line of the cut == ('cut','digest',0); off the cut None; hover_target same |
| 4 | BamHI (only unique cutter) name line | site_label_at above bases == ('site','BamHI',pos); on the base line None |

## `test_the_cut_marker_is_not_the_origin_mark`

**Result:** ✅ passed in 0.04 s · [line 246](../../tests/test_cut_sites.py#L246) · **Area:** GUI: cut cross vs origin bar painting · **Checked against:** GUI state · **Strength:** 🟡 adequate

In both themes a restriction cut is painted as a cross reaching either side of the strand, while a circle's origin mark is a plain vertical bar of a different colour (so the two are not confused).

| # | Input | Expected output |
|---|---|---|
| 1 | 'ACGT'*200 circular, dark False and True; EcoRV cut at bpl+10 on row 1 | cut-colour pixels present in small boxes on both sides above the top strand and below the bottom strand (the cross) |
| 2 | scroll to the last row, origin mark after the last base | no cut-colour-vs-origin crosses at the origin (origin box has no cross); origin-colour bar pixels present; cut colour != origin colour |

**Review notes:**

- L262 the cross is confirmed by counting any 'cut' pixels in 3x7 px boxes; it does not verify a cross shape, only that colour appears on both sides.
- L253 'assert view.rpc >= 3' is a precondition on the test fixture, not a behaviour under test.

## `test_a_saved_guide_site_is_drawn`

**Result:** ✅ passed in 7.22 s · [line 279](../../tests/test_cut_sites.py#L279) · **Area:** GUI: guide site drawing · **Checked against:** GUI state · **Strength:** 🟡 adequate

A saved guide site is drawn with coloured protospacer, PAM, name and cut arrow on both strands in both themes (PAM to the right of the protospacer on the top strand); a guide edited to another size is drawn as a plain feature.

| # | Input | Expected output |
|---|---|---|
| 1 | guide_doc(strand), dark False/True x strand 1/-1, sequence view 1800x300 | cas9_proto, cas9_pam, cas9_label each paint > 50 px |
| 2 | box widths measured column by column over the lane height (max over the ~10 px band, not one pixel row, because the 'PAM' label inside the box hides most of one row with wide fonts) | proto columns > 15*cw and pam columns > 2*cw; pam box is right of proto iff strand > 0 |
| 3 | guide edited to 'lead+1..lead+10' (wrong size) | cs.feature_layout is None; cas9_pam colour painted 0 times (plain feature) |

**Review notes:**

- L302 the width bound uses '> N*view.cw' as a proxy for an N-base box; it tolerates the box being wider or the colour bleeding, so an off-by-a-few-bases box would pass.

## `test_a_saved_guide_in_the_file`

**Result:** ✅ passed in 0.33 s · [line 340](../../tests/test_cut_sites.py#L340) · **Area:** GenBank I/O: guide round-trip · **Checked against:** round-trip/self-consistency · **Strength:** 🟢 strong

A saved guide round-trips byte for byte, Biopython reads it as a plain misc_feature (label, note, program qualifier), a reverse-strand guide keeps strand -1, and exports leave the private /cas9 qualifier out while keeping the readable note.

| # | Input | Expected output |
|---|---|---|
| 1 | saved guide (SpCas9) -> save_genbank a.gb -> reload -> save b.gb | a.gb bytes == b.gb bytes; reloaded layout.gaps==[g.cut]; proto bases == g.guide; pam bases == g.pam |
| 2 | read a.gb with Biopython SeqIO | misc_feature label ['gRNA '+g.guide[:8]]; note has guide and 'PAM '+g.pam; CAS9_KEY qualifier == cas9_qualifier(SP); start/end/strand(1) match; all lines <= 80 chars |
| 3 | guide_doc(-1) saved to r.gb | Biopython strand -1, start == lead; reload layout.pam == [(lead,lead+3)] |
| 4 | exports.tbl_text / gff3_text of the saved editor | CAS9_KEY absent from both; g.guide present in tbl_text (the plain note stays) |
