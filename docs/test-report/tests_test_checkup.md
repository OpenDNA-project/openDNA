# `tests/test_checkup.py`

[← all test files](README.md) · source: [`tests/test_checkup.py`](../../tests/test_checkup.py)

## `test_bacterial_cds_verdicts`

**Result:** ✅ passed in 0.03 s · [line 74](../../tests/test_checkup.py#L74) · **Area:** Checkup: bacterial CDS verdicts · **Checked against:** in-test reference model · **Strength:** 🟢 strong

In bacteria each CDS gets a verdict: an annotated promoter/RBS/terminator unit expresses it; otherwise elements are predicted from their motifs (either strand); a terminator before or inside the CDS blocks/truncates it; wrong strand or missing RBS lowers it; coding problems and regulation notes are reported.

| # | Input | Expected output |
|---|---|---|
| 1 | filler + Pcons promoter + RBS + gfp CDS + T1 terminator + filler, circular | gfp verdict YES; unit promoter 'Pcons', rbs 'rbs', terminator 'T1', unit.cds == ['gfp'] |
| 2 | only the CDS annotated, unit PROMOTER+RBS+gene+T1 on strand +1 and -1 | promoter, terminator and RBS all 'predicted' on that strand; verdict LIKELY, reason 'relies on predicted elements' |
| 3 | promoter then T1 terminator upstream, then RBS + CDS 'lost' | verdict NO, reason[0] has 'stopped by T1' |
| 4 | terminator T1 inserted inside gene(80), CDS 'cut' spanning it | verdict NO, reason[0] has 'inside the CDS' |
| 5 | CDS 'anti' on reverse strand after a forward promoter; CDS 'bare' with no RBS | anti NO 'no bacterial promoter upstream'; bare UNLIKELY 'no RBS' |
| 6 | broken CDS ATG+(GCT*10)+TAA+(GCT*10)+TAA; lacI gene downstream | broken UNLIKELY 'truncated protein'; lacI reason 'IPTG-inducible (lacI on the plasmid)', verdict YES |

**Review notes:**

- gene()/filler() use fixed random seeds, so a verdict change from a real regression and from a seed/codon change look the same (regression-value fragility).

## `test_eukaryotic_hosts_and_unit_ends`

**Result:** ✅ passed in 0.01 s · [line 140](../../tests/test_checkup.py#L140) · **Area:** Checkup: eukaryotic hosts · **Checked against:** in-test reference model · **Strength:** 🟢 strong

A CMV unit ends at its polyA and needs a Kozak context (weak Kozak / missing polyA are doubts); a U6 (Pol III) CDS is not translated; a eukaryotic promoter adds a eukaryote column; on a shuttle vector each unit stops at its own terminator kind and every CDS is judged per host.

| # | Input | Expected output |
|---|---|---|
| 1 | CMV promoter + Kozak + egfp + polyA, then U6 promoter + odd CDS + polyA | hosts == (BACTERIA, EUKARYOTES); egfp euk YES 'strong Kozak', bacteria NO, host EUKARYOTES; odd euk NO 'Pol III' |
| 2 | CMV promoter + 'TTTTTT' + weak CDS 'ATGA...' + filler, no polyA, linear | euk reasons 'weak Kozak' and 'no polyA'; verdict UNLIKELY |
| 3 | shuttle: Pcons/RBS/kanR/polyA/RBS/after/T1 then CMV/Kozak/egfp/T1b/bGH polyA | Pcons unit terminator 'T1' (polyA ignored), CMV unit terminator 'bGH polyA' (T1 ignored); after bacteria YES; kanR bacteria YES/euk NO; egfp euk YES/bacteria NO |

## `test_unit_extents_and_ends`

**Result:** ✅ passed in 0.05 s · [line 184](../../tests/test_checkup.py#L184) · **Area:** Checkup: unit extent and ending · **Checked against:** in-test reference model · **Strength:** 🟢 strong

A unit with no terminator runs round a circle (far CDS gets a distance note) or to the line end; a circular unit crosses the origin and can be added as an mRNA feature; +1 is exact only for an annotated TSS or a T7 promoter; unit endings are classed intrinsic / rho-dependent? / polyA / runs on.

| # | Input | Expected output |
|---|---|---|
| 1 | promoter/RBS/near CDS + 3500 bp filler + RBS/far CDS, circular and linear | terminator None, extent length = ed.length minus (0 circular else u.a); near and far YES; 'no terminator'; far has 'kb downstream of P' |
| 2 | RBS/CDS/T1 then filler then promoter (circular, unit wraps origin); add unit feature | CDS YES; unit_extent[1] > length; one mRNA feature with 2 spans |
| 3 | T7 promoter, Pcons+annotated TSS, Pplain (all with CDS) | T7 tss_exact, tss 17; Pcons tss_exact at annotated base; Pplain not exact; units.geometry exact {T7:True, Pcons:True, Pplain:False} |
| 4 | Pa+T1, Pb+odd terminator filler(60,32), CMV+bGH polyA; lone promoter+filler | Pa ending 'intrinsic' detail 'hairpin'; Pb 'rho-dependent?'; CMV 'polyA'; labels as expected; lone 'runs on' |

**Review notes:**

- L234 'ending_detail' containing 'hairpin' and the Pb 'rho-dependent?' both depend on the random filler(60,32) seq matching a hairpin motif only by seed; brittle but asserted exactly.

## `test_origins_and_prediction_thresholds`

**Result:** ✅ passed in 0.09 s · [line 243](../../tests/test_checkup.py#L243) · **Area:** Checkup: origins and thresholds · **Checked against:** hand-computed · **Strength:** 🟡 adequate

Replication origins are named with group/host and pir notes; the ColE1-group and no-origin remarks fire; prediction thresholds keep false hits rare in 50 kb of random DNA; SD/motif slices near a linear start are placed from base 0.

| # | Input | Expected output |
|---|---|---|
| 1 | rep_origins pUC/p15A/ColE1/R6K gamma/f1 (filler bodies) | names pUC, p15A, 'pMB1 / ColE1'; f1 host 'Phage'; R6K note has 'pir'; remark '2 origins of the ColE1 group' |
| 2 | R6K gamma ori + pir CDS | origin note 'pir is on the plasmid' |
| 3 | f1 ori + a CDS (no bacterial ori) | note 'no bacterial replication origin' |
| 4 | 50,000 bp random (seed 7) | <= 8 promoters at score 10 and <= 8 terminators at stem 18 |
| 5 | Strand 'AGGAGGTTTTATG'+GCT*20, linear; sd_match(st,10,20) | (6, 0, 6) = AGGAGG matched, slice starts at base 0 |

**Review notes:**

- L264-265 thresholds are '<= 8' on one random sequence: a loose bound that only catches gross over-prediction, not the intended rarity.
- L264-265 the oracle for the count is the code under test itself (predict_promoters/predict_terminators), so a systematic threshold bug would pass as long as it stays under 8.

## `test_library_parts_on_a_real_plasmid`

**Result:** ✅ passed in 0.03 s · [line 272](../../tests/test_checkup.py#L272) · **Area:** Checkup: library parts, real plasmid · **Checked against:** regression value · **Strength:** 🟡 adequate

On the real demo plasmid pBAD30_eilR.gb the checkup uses the built-in part library: the T2 terminator comes from the library and eilR is called expressed behind the araBAD promoter, inducible by arabinose with araC on the plasmid.

| # | Input | Expected output |
|---|---|---|
| 1 | demo/pBAD30_eilR.gb, library enabled (this test only) | some element origin 'library', kind 'terminator' |
| 2 | same | eilR verdict YES, unit promoter 'araBAD promoter'; reason has 'arabinose' and 'araC on the plasmid' |

**Review notes:**

- Depends on the committed demo file and the library contents; no reason given why eilR should be YES beyond those fixtures.

## `test_panel_tree_overlay_and_adding`

**Result:** ✅ passed in 0.03 s · [line 283](../../tests/test_checkup.py#L283) · **Area:** GUI: Checkup panel tree/overlay · **Checked against:** GUI state · **Strength:** 🟡 adequate

The Checkup dialog shows one tree of units and the CDSs they carry (verdict per host) then orphan CDSs; while shown it draws its findings as temporary overlay features under the document's filter; they vanish when hidden; adding the units does not draw them twice; a temp feature can be selected/hovered but not edited/deleted.

| # | Input | Expected output |
|---|---|---|
| 1 | Pcons/RBS/gfp/T1 + filler + lonely CDS; CheckupDialog | top items ['Pcons -> T1 (intrinsic)', NONE]; gfp YES child, lonely NO child; header ['Unit / CDS','Location',BACTERIA,'Notes'] |
| 2 | dialog not shown, then shown | overlay_features empty before; after show exactly one TU overlay in ed.shown, none in ed.features |
| 3 | set_display(units=False) then True | the temp feature leaves ed.shown/visible, then returns |
| 4 | select tree gfp child; select_feature(temp); edit_feature on temp | selection length == len(gene()); current_feature None; describe has 'Pcons'; edit_feature returns False |
| 5 | add_units(); refresh(); hide(); show(); close() | one TU in ed.features; no TU in overlay_features (not drawn twice); overlay empty when hidden |

## `test_checkup_tab_and_pane`

**Result:** ✅ passed in 0.18 s · [line 323](../../tests/test_checkup.py#L323) · **Area:** GUI: Checkup tab/pane wiring · **Checked against:** GUI state · **Strength:** 🟡 adequate

Tools > Plasmid Checkup opens one tab with checkup beside map and sequence (reused, not duplicated); picking checkup in a pane adds the map and sequence its tab lacks, a tab already showing them gets nothing more, and a pane closed afterwards stays closed.

| # | Input | Expected output |
|---|---|---|
| 1 | load demo plasmid; win.open_checkup() twice | tab count +1 then unchanged; pane views [TOOL,'Map','Sequence']; current index back to the checkup tab |
| 2 | new window: Features tab, pane.set_view(TOOL) | ordered panes [TOOL,'Map','Sequence']; after closing Map -> [TOOL,'Sequence'] |
| 3 | standard tab (has Map and Sequence): a Features pane set_view(TOOL) | pane views == sorted ['Sequence','Map',TOOL] (only the pane changes) |

## `test_units_as_features`

**Result:** ✅ passed in 0.12 s · [line 374](../../tests/test_checkup.py#L374) · **Area:** Units as features (units.py) · **Checked against:** hand-computed · **Strength:** 🟢 strong

Transcription-unit features carry geometry on both strands and across the origin (a plain mRNA is not a unit); the display filter has a dedicated units switch (not shown as mRNA), the views draw them without error on circle and line, projects remember the switch, exports omit the private qualifier, and the flat toolbar icons use one colour.

| # | Input | Expected output |
|---|---|---|
| 1 | unit_doc() forward unit Pcons->T1 | tss == len(PROMOTER), stopped, end == len(PROMOTER+RBS+gene+T1); pieces thin over promoter then thick from +1 |
| 2 | reverse unit Prev (no terminator, circular) | not stopped (runs on round the circle); tss at Prev's left end; thin over promoter, thick sums to length-len(PROMOTER); 2 spans |
| 3 | plain gfp mRNA feature | units.geometry None and units.tu_of None |
| 4 | DisplayPanel on circle and line; SequenceView and MapView render | 'mRNA' not in panel.kinds; render raises nothing; 2 TUs in ed.shown |
| 5 | panel.units unchecked; set_display(hidden_kinds={'mRNA'}, units=True) | no TU drawn; project meta display.units False; hiding mRNA type still leaves 2 TUs shown; gff3 lacks 'OpenGenEdit_tu' |
| 6 | icons.filter_icon / history_icon 32x32 | at most 2 opaque colours (one flat colour, antialias aside) |
