# `tests/test_primer_design.py`

[← all test files](README.md) · source: [`tests/test_primer_design.py`](../../tests/test_primer_design.py)

## `test_primer_check_levels_follow_their_limits`

**Result:** ✅ passed in 0.02 s · [line 75](../../tests/test_primer_design.py#L75) · **Area:** Primer quality checks · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** Without it, the green/orange/red verdicts on a primer could fire at the wrong value or ignore the user's own limits, so primers that will not work would look fine (or vice versa) before they are ordered.

Checks that the documented default limits are exactly what the settings declare and that they are stored in the data folder, that Tm, GC, 3' clamp, hairpin and self-dimer levels change exactly at their thresholds and follow the limits when they are moved, and that verdicts roll up with a short primer getting none.

| # | Input | Expected output |
|---|---|---|
| 1 | primer_db.settings() and the CHECKS tool settings | INI format file at data_dir()/primers.ini; PP.settings is primer_db.settings; declared defaults equal DEFAULT_LIMITS exactly |
| 2 | band(v, 52, 65, 7) for 51.9, 52, 65, 65.1, 72, 72.1, 45, 44.9 | warn, ok, ok, warn, warn, bad, warn, bad |
| 3 | Tm values 51.9 / 52 / 58 / 65 / 65.1 / 72 / 72.1 / 45 / 44.9 / 20 | warn, ok, ok, ok, warn, warn, bad, warn, bad, bad |
| 4 | tm_low 55 and tm_margin 2, Tm 53.0 then 52.9 | warn then bad - moving the limit moves the verdict |
| 5 | 20 nt primers with 4,5,6,7,8,12,13,14,15,20 G's | bad, bad, warn, warn, ok, ok, warn, warn, bad, bad; with gc_low 45 / gc_margin 5 the 8/9/10 cases give warn, ok, ok; lower-case bases count the same |
| 6 | primers ending ATATA / ATATG / AGCTA / AGCTG / CGCTG / GCGCG | 0 GC in the last 5 -> warn; 1-3 -> ok; 4-5 -> warn; with clamp_max 4 the 4-GC end becomes ok |
| 7 | a clean primer, a GGATCC...GGATCC loop, a GC-clamped hairpin, and dimer-prone primers | ok/ok; hairpin warn then bad; self-dimer warn then bad |
| 8 | hairpin_bad set 0.01 either side of the primer's own hairpin dG (and the same for dimer3_bad / dimer_bad) | verdict flips warn <-> bad exactly at the crossing; setting only the 'anywhere' limit still gives bad |
| 9 | a 4-base primer / worst() over mixed levels / checks_html | [('Primer','info','Type at least 8 bases.')]; worst gives warn, bad, info and info for an empty list; the HTML carries the ok and bad marks and the value text |

**Review notes:**

- Tm is stubbed to a fixed value, so the actual melting-temperature calculation is not exercised here.

## `test_where_a_primer_comes_from_decides_how_it_is_kept`

**Result:** ✅ passed in 0.02 s · [line 160](../../tests/test_primer_design.py#L160) · **Area:** Primer sources and storage · **Checked against:** round-trip/self-consistency · **Strength:** 🟢 strong

**Why it matters:** Without it, a primer from the shared database could be confused with one stored in the file, be drawn in the wrong colour or hidden by the wrong filter, or an imported primer could fail to be saved or to undo cleanly.

Checks the flags, colour and display-filter behaviour of the four primer origins (database only, file, both, designed here), which primers can be imported into the file, that importing is one undo/redo step with the right history text, and that imported, designed and newly made primers are all written to the file and come back as the file's own.

| # | Input | Expected output |
|---|---|---|
| 1 | editor whose file holds 'mine' and 'both', database holds 'dbonly' and 'both', 'temp' designed in session | four sites; dbonly (db, not file, provisional); mine (file, not provisional); both (db and file, not provisional - the file's style wins); temp temporary and provisional |
| 2 | primer_color of each | dbonly and temp use TEMP_COLOR; mine and both use PRIMER_COLOR |
| 3 | display settings primers/db_primers/file_primers toggled | all four drawn; only dbonly with file_primers off; mine, both and temp with db_primers off; nothing with primers off |
| 4 | feature_table.primers_to_import(ed, ['dbonly','temp','both']) | [('dbonly', DBONLY), ('temp', TEMP)] - 'both' is already in the file |
| 5 | import_primers([('dbonly', DBONLY)]) | exactly one new history step, text 'Imported primer dbonly to this sequence'; the site becomes from_file, not provisional, PRIMER_COLOR, described as 'in this file'; nothing left to import |
| 6 | undo then redo | provisional again with 'not in this file'; redo restores the imported state |
| 7 | importing 'temp', saving, reopening the file | all of mine, both, dbonly, temp written as primer_bind features; after removing dbonly from the database the reloaded file still has all four as regular file primers |
| 8 | a fresh editor on the unsaved document | 'temp' still present and provisional |
| 9 | add_primer(ed, 'madehere', ...) then save | not provisional and written into the file |

## `test_primers_are_designed_for_a_target_at_or_across_the_origin_of_a_circle`

**Result:** ✅ passed in 0.08 s · [line 220](../../tests/test_primer_design.py#L220) · **Area:** Primer design on circular DNA · **Checked against:** property/invariant · **Strength:** 🟡 adequate

**Why it matters:** Without it, designing primers for a target at or across a plasmid's origin could return nothing, or primers whose bases are not actually on the template or whose product does not contain the target.

Checks that on a circular 3000 bp sequence the designer returns pairs for a target at base 1, at the last base and across the origin, that each primer's bases really match the template read round the circle, that the template between primer and target stays within 150 bp, that the reported product size is consistent, and that a linear sequence has nothing before base 1.

| # | Input | Expected output |
|---|---|---|
| 1 | circular 3000 bp, target (0,200) | non-empty pairs |
| 2 | target (2800,3000) and (2900,3100) (across the origin) | non-empty pairs in both cases |
| 3 | each returned pair | fwd bases equal (seq*2)[start:end]; rev bases equal the reverse complement of (seq*2)[start:end] |
| 4 | the gap between each primer and the target | both (start - fwd.end) mod n and (rev.start - end) mod n are <= 150 |
| 5 | the reported size | size == len(fwd) + before + (end-start) + after + len(rev) |
| 6 | the same sequence as a line, targets (0,200) and (500,700) | no pairs for the first (nothing upstream of base 1), pairs for the second |

**Review notes:**

- Only checks that pairs exist and are self-consistent; nothing here checks that the chosen primers are good ones (Tm, specificity), so a designer returning poor pairs would pass.

## `test_the_files_written_for_ordering`

**Result:** ✅ passed in 0.03 s · [line 241](../../tests/test_primer_design.py#L241) · **Area:** Oligo order files · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** Without it, the plate file sent to a synthesis company could put oligos in the wrong wells, overflow past well 96, drop the scale or purification, or corrupt names and sequences - a real-money error.

Checks that the IDT plate export fills wells down the columns, starts a new plate file after 96 oligos, writes scale and purification on every row, handles odd names and degenerate bases in both formats, and refuses an empty order without writing anything.

| # | Input | Expected output |
|---|---|---|
| 1 | 10 oligos exported to o.csv | header ['Well Position','Name','Sequence','Scale','Purification']; wells A1,B1,...,H1,A2,B2 |
| 2 | 100 oligos | two files o_plate1.csv and o_plate2.csv; first has 96 rows ending at H12; second has 4 rows starting at A1 with name P96 |
| 3 | 2 oligos with scale '100nm' and purification 'HPLC' | every row carries ('100nm','HPLC') |
| 4 | name 'fwd, "A" primer' in plate and plain formats | name preserved verbatim in both; plain file is [name, sequence] |
| 5 | name containing a tab and a newline, plain format | written as 'a b c' (whitespace becomes a space) |
| 6 | sequence 'acgtNNRYacgt', plain format | 'ACGTNNRYACGT' - degenerate bases kept, upper-cased |
| 7 | an empty order into an empty directory | ValueError and no file created |
