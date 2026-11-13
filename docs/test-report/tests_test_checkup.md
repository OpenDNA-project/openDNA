# `tests/test_checkup.py`

[← all test files](README.md) · source: [`tests/test_checkup.py`](../../tests/test_checkup.py)

## `test_bacterial_cds_verdicts`

**Result:** ✅ passed in 0.03 s · [line 77](../../tests/test_checkup.py#L77) · **Area:** Plasmid checkup (bacteria) · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** Without it, the Checkup could tell a user that a gene will be expressed in E. coli when a terminator blocks it, when it is on the wrong strand or when it has no RBS - or refuse to predict the parts a plain file does not annotate.

Checks the verdict given to each CDS in a bacterial context: an annotated promoter/RBS/terminator unit expresses it, predicted motifs on either strand give a 'likely' verdict, terminators before or inside the CDS block or truncate it, a CDS on the wrong strand or without an RBS is downgraded, and coding problems and regulation notes appear.

| # | Input | Expected output |
|---|---|---|
| 1 | filler + consensus sigma-70 promoter + RBS + gene + rrnB T1 + filler, annotated | verdict YES, unit promoter 'Pcons', rbs 'rbs', unit terminator 'T1', unit cds == ['gfp'] |
| 2 | the same parts unannotated except the CDS, on the + and the - strand | promoter, terminator and RBS all predicted on the CDS's strand; verdict LIKELY with reason 'relies on predicted elements' |
| 3 | promoter, T1, filler, RBS, CDS | verdict NO, first reason contains 'stopped by T1' |
| 4 | promoter, RBS, first 120 bp of the gene, T1, rest of the gene, CDS annotated over all of it | verdict NO, first reason contains 'inside the CDS' |
| 5 | promoter + CDS annotated on the minus strand; and a second promoter + CDS with a C-rich spacer instead of an RBS | NO with 'no bacterial promoter upstream'; UNLIKELY with 'no RBS' |
| 6 | a CDS with an internal stop, plus a lacI CDS on the plasmid | UNLIKELY with a 'truncated protein' reason; lacI CDS gets the note 'IPTG-inducible (lacI on the plasmid)' and still verdict YES |

## `test_eukaryotic_hosts_and_unit_ends`

**Result:** ✅ passed in 0.01 s · [line 139](../../tests/test_checkup.py#L139) · **Area:** Plasmid checkup (eukaryotes) · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** Without it, a mammalian expression construct could be judged by bacterial rules - no Kozak check, a polyA or a Pol III promoter ignored - and a shuttle vector would not get a per-host verdict.

Checks that a CMV unit ends at its polyA signal and is judged on Kozak context rather than an RBS, that a weak Kozak and a missing polyA are doubts, that a U6 (Pol III) promoter yields no translation, that a eukaryotic promoter adds a eukaryote column, and that on a shuttle vector each unit ends at its own kind of terminator.

| # | Input | Expected output |
|---|---|---|
| 1 | CMV promoter + GCCACC Kozak + egfp + polyA, then U6 promoter + CDS + polyA | hosts == (BACTERIA, EUKARYOTES); egfp YES in eukaryotes with 'strong Kozak' in the rbs label and NO in bacteria (host EUKARYOTES); the U6 CDS NO with 'Pol III' in the first reason |
| 2 | CMV promoter + TTTTTT + CDS starting ATGA, no polyA, linear | eukaryote reasons include 'weak Kozak' and 'no polyA'; verdict UNLIKELY |
| 3 | shuttle construct: Pcons/RBS/kanR/polyA/RBS/after/T1 then CMV/Kozak/egfp/T1/bGH polyA | Pcons unit ends at T1 (polyA does not stop the bacterial polymerase); CMV unit ends at bGH polyA (T1 does not stop the eukaryotic one) |
| 4 | the same construct's CDS verdicts | 'after' YES in bacteria; kanR YES in bacteria and NO in eukaryotes; egfp YES in eukaryotes and NO in bacteria |

**Review notes:**

- The CMV and U6 promoters are random filler named 'CMV promoter' / 'U6 promoter', so this exercises name-based recognition only, not any sequence motif.

## `test_how_far_a_unit_runs_and_how_it_ends`

**Result:** ✅ passed in 0.04 s · [line 183](../../tests/test_checkup.py#L183) · **Area:** Transcription units · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** Without it, the drawn transcript could stop at the wrong base, miss that it wraps the origin, claim an exact +1 it cannot know, or describe its ending wrongly, misleading the user about what is actually transcribed.

Checks that a promoter without a terminator runs round a circle or to the end of a line (with a distance note on a far CDS), that a unit crossing the origin becomes a two-span mRNA feature, that +1 is only flagged exact where the start base is really known, and that each unit ending is classified and labelled.

| # | Input | Expected output |
|---|---|---|
| 1 | promoter, RBS, near CDS, 3500 bp filler, RBS, far CDS - circular and linear | unit terminator None, extent reaching the end (whole length when circular); both CDS YES; 'no terminator' reason on the near one; 'kb downstream of P' note on the far one |
| 2 | RBS + CDS + T1 + filler + promoter on a circle (promoter last) | CDS YES; unit extent end > sequence length (passes the origin); the added unit feature is one mRNA with 2 spans |
| 3 | T7 promoter / promoter with an annotated TSS feature / plain consensus promoter | T7 tss_exact with tss == 17; Pcons tss_exact at the annotated base; Pplain not exact; after writing unit features, units.geometry reports exact True/True/False |
| 4 | promoter + T1 / promoter + a T-poor 60 bp 'terminator' / CMV promoter + polyA | endings 'intrinsic' (detail mentions hairpin), 'rho-dependent?', 'polyA'; labels 'Pa -> T1 (intrinsic)' and 'CMV promoter -> bGH polyA' |
| 5 | a lone promoter with only filler after it | ending 'runs on', label 'P -> runs on' |

## `test_origins_and_prediction_thresholds`

**Result:** ✅ passed in 0.07 s · [line 241](../../tests/test_checkup.py#L241) · **Area:** Replication origins and prediction noise · **Checked against:** property/invariant · **Strength:** 🟡 adequate

**Why it matters:** Without it, the Checkup could name the wrong origin family (or miss that a plasmid has none, or needs a pir host), and the motif predictors could flood a real plasmid with chance hits.

Checks that annotated rep_origin features are mapped to origin names, hosts and notes (including the R6K pir note and the ColE1-group and no-origin remarks), that the promoter and terminator predictors fire rarely on 50 kb of random DNA, and that a motif near the start of a linear sequence is placed at the right base.

| # | Input | Expected output |
|---|---|---|
| 1 | rep_origin features named pUC ori, p15A ori, ColE1 origin, R6K gamma ori, f1 ori | names 'pUC', 'p15A', 'pMB1 / ColE1'; f1 host 'Phage'; R6K note mentions pir |
| 2 | R6K gamma ori plus a CDS labelled 'pir' | note == 'pir is on the plasmid' |
| 3 | the five-origin construct / a plasmid with only f1 ori and a gene | a note about '2 origins of the ColE1 group'; a note 'no bacterial replication origin' |
| 4 | 50 000 random bases, predict_promoters(threshold 10) and predict_terminators(18) | at most 8 predictions each |
| 5 | sd_match on AGGAGGTTTTATG + GCT x 20 from base 10 to 20, linear | (6, 0, 6) |

**Review notes:**

- The <= 8 bounds on a single random sequence are a loose regression threshold, not a real false-positive rate; a predictor that returned zero hits would also pass.
- Origin recognition here is driven by the feature labels, not the sequence.

## `test_library_parts_on_a_real_plasmid`

**Result:** ✅ passed in 0.04 s · [line 269](../../tests/test_checkup.py#L269) · **Area:** Plasmid checkup (parts library) · **Checked against:** smoke · **Strength:** 🟡 adequate

**Why it matters:** Without it, the built-in parts library could override a file's own annotations, or fail to supply the terminator a bare sequence needs, changing the verdicts on real plasmids.

Runs the Checkup on the annotated demo plasmid pBAD30_eilR and on the same bases with no annotations, checking that the library adds nothing to the annotated one but supplies a terminator to the bare one, and that eilR gets a YES verdict behind the araBAD promoter with an arabinose/araC note.

| # | Input | Expected output |
|---|---|---|
| 1 | demo/pBAD30_eilR.gb opened in an Editor, C.run | no element with origin 'library' |
| 2 | the same bases as a new unannotated circular document | at least one library-sourced terminator element |
| 3 | the eilR CDS verdict in the annotated file | verdict YES, unit promoter 'araBAD promoter'; a reason mentioning arabinose and 'araC on the plasmid' |

**Review notes:**

- Only three assertions over one file; it shows the library is consulted but not that any particular part is matched correctly.
- Tied to a demo file whose annotation changes would change the result.

## `test_unit_features_keep_their_geometry_and_survive_a_file_round_trip`

**Result:** ✅ passed in 0.04 s · [line 284](../../tests/test_checkup.py#L284) · **Area:** Transcription unit features · **Checked against:** round-trip/self-consistency · **Strength:** 🟢 strong

**Why it matters:** Without it, a transcription unit saved into a file could be drawn with the wrong +1, the wrong thick/thin split or lost altogether on reload, and old files written under the program's former name would stop working.

Checks the geometry of unit features on both strands and across the origin (thin over the promoter, thick from +1, 'stopped' or running on), that a plain mRNA is not a unit, that units have their own display switch independent of the mRNA type, and that saving and reloading (including under the legacy /OpenGenEdit_tu key) gives identical geometry while exports omit the private qualifier.

| # | Input | Expected output |
|---|---|---|
| 1 | Pcons/RBS/gfp/T1 + filler + a reverse promoter, unit features added | forward unit: tss == len(PROMOTER), stopped, end == len(PROMOTER+RBS+gene+T1), pieces [(0,tss,thin),(tss,end,thick)] |
| 2 | the reverse-strand unit | not stopped (runs round the circle); tss at Prev's left end; one thin piece over the promoter; thick pieces total length - len(PROMOTER); the feature has 2 spans (crosses the origin) |
| 3 | the plain 'gfp' feature | units.geometry is None and tu_of is None |
| 4 | set_display(units=False) / hidden_kinds={'mRNA'} with units=True | 2 units shown normally; none shown when units off (also in map_features); project_meta display.units is False; hiding the mRNA kind still leaves the 2 units shown |
| 5 | saving to GenBank | the file contains '/OpenDNA_tu=' and no 'OpenGenEdit' |
| 6 | the same file with the key renamed to /OpenGenEdit_tu, reloaded | 2 units whose geometry equals the originals exactly; both shown; the display switch still works |
| 7 | gff3 and tbl exports of both documents | neither TU_KEY nor 'OpenGenEdit_tu' appears |
