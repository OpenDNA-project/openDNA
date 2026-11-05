# `tests/test_cas9.py`

[← all test files](README.md) · source: [`tests/test_cas9.py`](../../tests/test_cas9.py)

## `test_candidates_circles_other_nucleases_score_and_seed_mismatches`

**Result:** ✅ passed in 0.00 s · [line 41](../../tests/test_cas9.py#L41) · **Area:** Cas9/Cas12a: candidate guide finding · **Checked against:** hand-computed · **Strength:** 🟢 strong

candidates() finds guides on both strands with correct coordinates, PAM and cut site; handles an offset region, sites across a circular origin, and SaCas9's NNGRRT PAM; score_guide applies its checklist; seed_mismatches counts seed-region mismatches.

| # | Input | Expected output |
|---|---|---|
| 1 | dna(40)+GUIDE+'TGG'+dna(40), SpCas9 | 1 guide; (start,end,strand,guide,pam)==(40,63,1,GUIDE,'TGG'); cut==40+20-3==57 (3bp from PAM) |
| 2 | dna(40)+rc(GUIDE+'AGG')+dna(40) | (40,63,-1,GUIDE,'AGG'); cut==40+3+3==46; bottom[start:end]==rc(guide+pam) |
| 3 | offset=50 on a slice; circular site across origin (wrap); SaCas9 NNGRRT | guide read at its full-seq place (start 60); circular site start==len-12, end>len, no site on the line; SaCas9 guide len 21, pam 'ATGAAT'; 'ATGCCT' rejected |
| 4 | score_guide of GUIDE, TTTT seq, poly-G, poly-AT, clean | (100,[]); <70 with 'TTTT' note; <=45 with GC+run notes; <=60 with GC note; ==100 |
| 5 | candidates on seq with GUIDE + a low-score guide | scores sorted descending; found[0].guide==GUIDE |
| 6 | seed_mismatches(1,20),(1<<19,20),((1<<19)\|(1<<7)\|1,20),((1<<8)\|(1<<12),20) | 0, 1, 1, 2 |

## `test_other_sites_here_and_in_a_genome_and_the_rust_search`

**Result:** ✅ passed in 0.01 s · [line 93](../../tests/test_cas9.py#L93) · **Area:** Cas9: off-target counting and Rust search · **Checked against:** independent implementation · **Strength:** 🟢 strong

scan() counts other matching sites (by mismatch count) in the open sequence and in a reference genome, verdict() grades them, and the Rust find_guide_sites matches an independent Python implementation.

| # | Input | Expected output |
|---|---|---|
| 1 | seq with GUIDE + near(1mm) + far(4mm); scan max_mm=3 | here==[0,1,0,0] (4mm beyond 3 excluded); genome None; verdict 'warn' |
| 2 | genome holding guide perfect + 2mm copy | genome==[1,0,1,0], name 'ref'; verdict 'warn' |
| 3 | genome with two perfect copies; or two perfect sites here | genome[0]==2 verdict 'bad'; both here==1 verdict 'bad' |
| 4 | reference is the open sequence itself; spec None | genome None (not scanned twice); verdict(None,g) 'info' |
| 5 | random 4000bp seq, guides at 100/1000/2500, core.find_guide_sites NGG max 2 | Rust hit set {(start,strand)} == independent Python scan set for each guide |

## `test_the_tool_in_a_pane_with_cas9_cas12a_and_the_reference_genome`

**Result:** ✅ passed in 0.47 s · [line 163](../../tests/test_cas9.py#L163) · **Area:** GUI: Cas9 tool pane · **Checked against:** hand-computed · **Strength:** 🟢 strong

The Spacer Finder tool lists/selects/checks/adds guides, draws Cas9 and Cas12a cuts, narrows to a selection, is a Tools-group pane view, and shares the Primers tool's reference genome (including turning the check off).

| # | Input | Expected output |
|---|---|---|
| 1 | seq with two GUIDE+PAM copies; Cas9Dialog | 2 rows, '2 guide' in status; selecting row selects site; cut mark top==bottom==80+20-3 |
| 2 | start_check with no genome; wait | item(0,6)=='1/0/0/0', (0,7)=='–'; (0,8) has '✗' and 'other perfect site' |
| 3 | add guide as feature | one new misc_feature, label starts 'gRNA ', note has 'SpCas9 guide '+GUIDE and 'cuts 3 bp from the PAM' |
| 4 | select 0..60, reload | rowCount 0; no cut marks |
| 5 | pane.set_view('Spacer Finder') | 'Spacer Finder' in grouped Tools; view is Cas9Dialog; tab labelled |
| 6 | Cas12a mode on TTTA+GUIDE23 seq | 1 row cell 'TTTA '+GUIDE23; cut mark (50+4+18,50+4+23,'Cas12a'); note starts 'Cas12a guide '+GUIDE23 with '5′ overhang' |
| 7 | set primers/genome to ref holding guide once; primer_options_changed; start_check | 'myref' in genome label; cas9.genome() is primer_panel.genome(); (0,7)=='1/0/0/0', (0,8) '✓'; genome_check off -> 'none set' |

**Review notes:**

- start_check uses wait_for() busy loops with 5s timeouts and QTest.qWait (lines 179,242) running real worker-thread searches; slow machines could time out. wait_for asserts its condition so a timeout fails loudly rather than silently.

## `test_cas12a_sites_staggered_cuts_and_other_sites_counted_with_the_pam_first`

**Result:** ✅ passed in 0.00 s · [line 260](../../tests/test_cas9.py#L260) · **Area:** Cas12a: 5' PAM and staggered cuts · **Checked against:** hand-computed · **Strength:** 🟢 strong

Cas12a sites have the TTTV PAM on the 5' side, the two strands cut in different places (5' overhang), the PAM-overlap rule is handled (TTTT not a PAM), and off-target counting puts the PAM before the protospacer.

| # | Input | Expected output |
|---|---|---|
| 1 | FILL+'TTTA'+GUIDE23+FILL, A12 | (start,end,strand,guide,pam)==(50,77,1,GUIDE23,'TTTA'); (cut,bottom)==(50+4+18,50+4+23) (5nt overhang) |
| 2 | FILL+rc('TTTC'+GUIDE23)+FILL | (50,77,-1,GUIDE23,'TTTC'); (cut,bottom)==(50,55); bottom[start:end]==rc(pam+guide) |
| 3 | FILL+'TTTT'+GUIDE23+FILL | site starts at [51] (TTTG after first T is the PAM, TTTT is not); cut_text has '5′ overhang'(A12) and '3 bp'(SpCas9) |
| 4 | two Cas12a sites; scan guides[:1] max_mm 1 | spec.here==[1,0] (other copy found PAM-first, guide itself excluded) |
