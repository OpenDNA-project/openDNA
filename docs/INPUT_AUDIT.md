# Input field audit

Every place where the user types or chooses a value was checked by reading the code (`grep` for `QLineEdit`, `QSpinBox`,
`QDoubleSpinBox`, `QPlainTextEdit`, `QComboBox`, `QInputDialog` in `app/`) and by feeding bad input to the fields
programmatically (empty, blank, very long, negative, huge, non-ASCII, wrong alphabet, lowercase, IUPAC, pasted FASTA or
numbered blocks). The tests for the fixes are in `tests/test_input_validation.py`.

Verdicts: **ok** = bad input is refused or answered with a clear message, nothing can break; **weak** = it works but the message
is poor, the value is unbounded or the dialog loses what you typed; **missing** = a traceback, a nonsense result or no check at all.
"Before" is the state when the audit started, "now" is after the fixes of this change.

Fixed values are marked **FIXED** and say what was done. Spin boxes and drop-down lists (non-editable combos) cannot hold a
value outside their range, so they are **ok** unless two of them depend on each other.

## Text fields and dialogs

| # | Where | Expected | Validation before | Bad input before | Enter / OK | Before | Now |
|---|---|---|---|---|---|---|---|
| 1 | Search box (`search_bar.py`) | DNA, protein motif, or a name | none; trimmed | any text is searched; a pasted FASTA loses its line breaks and is searched as a name; 100,000 characters were accepted | Enter = next hit, Shift+Enter = previous | weak | ok: **FIXED** max 5,000 characters |
| 2 | Search options: fuzzy mismatches (spin 1-10) | integer | range | cannot be exceeded; capped at a third of the query | live | ok | ok |
| 3 | Toolbar range box ("101..200", or one base) | bases | 1..length | out of range or not a range: said in the status bar, the text is kept to correct; an end below the start = through the origin on a circle | Enter selects | ok | ok |
| 4 | Edit > Go to position (`QInputDialog.getInt`, 1..length) | base number | range | non-numbers are not accepted by the spin box | OK / Enter | ok | ok |
| 5 | Make linear: break before base (`getInt`, 1..n) | base number | range | same | OK / Enter | ok | ok |
| 6 | Alignment view locus box (`genome_view.py`) | `12000`, `12,000-15,000` or a name | regex for numbers; names looked up | `,` gave a `ValueError` traceback; `99999999999999999999` an `OverflowError` in Qt; `0` and positions past the end were passed on | Enter | missing | ok: **FIXED** a message ("positions run from 1 to N", "is not a position"), nothing is moved |
| 7 | Feature dialog: Name (`feature_dialog.py`) | short text | none | 100,000 characters accepted; trimmed | OK | weak | ok: **FIXED** max 200 (never below an existing label, so editing never cuts one) |
| 8 | Feature dialog: Type (editable combo) | one word | core check "must be one word" | a clear message; dialog stays open | OK | ok | ok |
| 9 | Feature dialog: Location | GenBank syntax inside the sequence | core parser | empty or `banana` gave "bad location ''... Syntax error: Error(Error { input: [], ..." or "outside the 500 bp sequence" | OK | weak | ok: **FIXED** blank and non-location text give "Cannot read the location ... 10..200, complement(...), join(...)"; the field is focused; whitespace is trimmed |
| 10 | Feature dialog: Qualifiers (multi-line `key=value`) | one qualifier per line | none | `=value` and `two words=x` were stored with an empty or invalid key, which writes an unreadable GenBank line | OK | missing | ok: **FIXED** the key must be one word (letters, digits, `_`, `'`, `-`); a message names the line |
| 11 | Feature dialog: Direction, Colour | choice | drop-down / colour dialog | cannot be invalid | OK | ok | ok |
| 12 | Feature table filter (`feature_table.py`) | free text | none needed | case-insensitive substring | live | ok | ok |
| 13 | Group as gene: name (`QInputDialog.getText`) | text | trimmed | an empty name did nothing, silently | OK | weak | ok: **FIXED** "Give the gene a name." |
| 14 | Primer dialog: Name (`primer_panel.py`) | short unique text | none | tabs, line breaks and a 10,000-character name were stored in the database | Save | missing | ok: **FIXED** max 60 (never below an existing name), control characters refused with a message, duplicates already refused |
| 15 | Primer dialog: Sequence | DNA, IUPAC, either case, 8-200 bases | regex validator (IUPAC + spaces) | typing wrong letters did nothing; pasting FASTA, numbers or a header was refused whole; 20,000 bases accepted (the checks stall); a selected genome pre-filled the box | Save | weak | ok: **FIXED** paste is cleaned (header, numbers, line breaks dropped; a tooltip says why when it is not DNA), max 200, a longer selection leaves the box empty, messages say how many bases are missing |
| 16 | Primer dialog: Note | text | none | unbounded | Save | weak | ok: **FIXED** max 200 |
| 17 | Primer dialog: tab order | rows in order | creation order | Sequence > Note > Reverse complement, not the order on screen | | weak | ok: **FIXED** `setTabOrder` |
| 18 | Primer database filter | free text | none needed | substring | live | ok | ok |
| 19 | Primer options: Na⁺/K⁺ (mM) and Primer (nM) | above 0 | minimum 0 | 0 gave a Tm of -273.15 °C | OK | missing | ok: **FIXED** minimum 1, also applied to values read from an old settings file |
| 20 | Primer options: Mg²⁺, dNTPs | 0 or more | range 0..100 | 0 is a real condition | OK | ok | ok |
| 21 | Primer options: mismatches (0-3), 3' bases (6-20), 3' mismatches (0-3) | integers | ranges | cannot be exceeded | OK | ok | ok |
| 22 | Primer options: genome file | path of a FASTA / GenBank file | none | a path that does not exist was saved and the genome check silently found nothing | OK | weak | ok: **FIXED** "There is no file ..." and the dialog stays open |
| 23 | Design Primers: Tm (45-75), window (20-2000) | numbers | ranges | cannot be exceeded | Enter | ok | ok |
| 24 | Design Primers: min / max length | min <= max | two separate ranges | min 30 with max 25 was silently repaired at search time | | weak | ok: **FIXED** moving one pushes the other |
| 25 | Design Primers: Save pair name (`getText`) | text | trimmed, empty cancels | tabs and very long names stored | OK | weak | ok: **FIXED** same rules as the primer dialog |
| 26 | File > New: Name | text | trimmed | 100,000 characters accepted | OK | weak | ok: **FIXED** max 100 |
| 27 | File > New: Sequence (paste box) | DNA, FASTA or numbered blocks | `clean_pasted` after OK | wrong letters: the dialog closed, a warning appeared and everything typed was lost | OK | weak | ok: **FIXED** the message appears under the box and the dialog stays open |
| 28 | Rename tab (`getText`) | short text | trimmed | a 300-character name made the tab bar scroll; line breaks kept | OK | weak | ok: **FIXED** one line, max 60 |
| 29 | Import from NCBI: query | accession or words | empty gives "Type an accession..." | clear message, runs in a thread | Enter / Search | ok | ok: **FIXED** also max 500 characters |
| 30 | Import from NCBI: API key (the e-mail field was removed) | optional | none | a key pasted with a trailing space or line break made every request fail with an unclear NCBI error | | weak | ok: **FIXED** trimmed before use, max 100 / 200 |
| 31 | Reference genome from NCBI: query | organism or GCF_ accession | empty message | clear message | Enter | ok | ok |
| 32 | 3D structure: Accession | UniProt accession | checked when getting the structure | "AlphaFold DB needs a UniProt accession" in the status line | Get structure is the default button | ok | ok |
| 33 | 3D structure: Local command | command with `{fasta}` and `{pdb}` | checked when run | clear message | | ok | ok |
| 34 | Map Reads: Reads, Mate reads, Output | existing FASTQ files, output BAM | file existence only | a missing file just moved the focus: nothing happened and no message; an empty output made a file called `.bam` | Run | weak | ok: **FIXED** a message for each case (no reads, file not found, no output, folder missing) |
| 35 | Map Reads: Threads (1-256) | integer | range | cannot be exceeded | Run | ok | ok |
| 36 | Local Align: query box | DNA / FASTA, up to 200,000 bp | `clean_pasted`; length check | protein or text gets a message instead of a result; empty clipboard has a message | Align | ok | ok |
| 37 | Detect ORFs: minimum length (10-5000) | integer | range | cannot be exceeded | live | ok | ok |
| 38 | Detect ORFs: protein motif | amino acids | regex validator | wrong letters ignored; pasting FASTA was refused whole | live | weak | ok: **FIXED** pasted FASTA / numbers are cleaned |
| 39 | Digest: filter, count menu | text, choice | none needed | substring | live | ok | ok |
| 40 | PCR product: forward / reverse primer lists | choice | non-editable | cannot be invalid | Open product | ok | ok |
| 41 | Spacer Finder: nuclease, region, mismatches (1-4) | choices, integer | non-editable, range | cannot be exceeded | Find guides | ok | ok |
| 42 | Assembly: method, enzyme, count (1-12), Tm spins | choices, numbers | ranges | cannot be exceeded | | ok | ok |
| 43 | Assembly: overlap min / max (15-60 / 15-80) | min <= max | two separate ranges | min 50 with max 20 allowed | | weak | ok: **FIXED** moving one pushes the other |
| 44 | Export image: DPI (72-2400) and size (25-1000 %) | numbers | ranges | 2400 DPI at 1000 % of a large tab asked for a raster of more than 100 million pixels (about 400 MB): out of memory | Export | weak | ok: **FIXED** raster formats over 100 million pixels disable Export with a message (SVG is not limited) |
| 45 | Tools > Settings: every tool parameter (spins / check boxes / lists, `tool_settings.py`) | per-parameter range | each `Param` has a minimum and maximum | cannot be exceeded; Defaults button | OK | ok | ok |
| 46 | Protein view: frame list, Settings of Detect / Checkup | choices | non-editable | cannot be invalid | | ok | ok |
| 47 | Track sequence / Which sequence (`QInputDialog.getItem`, not editable) | choice | non-editable | cannot be invalid | OK | ok | ok |
| 48 | Updater release notes | read-only text | read-only | nothing to enter | | ok | ok |
| 49 | Import → Primers → Primer List (file) | FASTA or name + sequence per line | a line was kept if its 2nd column was IUPAC; FASTA not checked at all | a FASTA record with any letters was stored; lines with no separator vanished silently; names with tabs and 500-letter "primers" were stored | Open | missing | ok: **FIXED** every primer goes through `primer_problem`; each refused line is reported with its number and reason; binary files are refused |
| 50 | Primer database, every way in (dialog, Design, Assembly, list, primers of an opened file) | valid primer | the dialog only | Assembly and opened files could store anything | | weak | ok: **FIXED** `PrimerDB.add` / `add_many` / `update` refuse an invalid primer ("invalid", with the reason) |
| 51 | Import Features to Library (files) | named features of A/C/G/T/N | non-ACGTN sequences dropped silently | a tab or line break in a feature name broke the TSV library | Open | missing | ok: **FIXED** `part_problem`: refused features are listed with the reason |
| 52 | Feature dialog: Name, qualifier values | one line of text | maximum length | a pasted tab or control character was written to the GenBank file | OK | weak | ok: **FIXED** a message names the problem, the dialog stays open |
| 53 | Group as gene name, Rename tab | one line of text | trimmed | a pasted control character was kept | OK | weak | ok: **FIXED** a message |
| 54 | Import from NCBI: API key | letters and digits | trimmed | any character was sent to NCBI | Search / Open | weak | ok: **FIXED** a wrong character is named in the status line as it is typed, and nothing is sent |

## Summary

| | Before | Now |
|---|---|---|
| ok | 26 | 54 |
| weak | 22 | 0 |
| missing | 6 | 0 |
| **Rows** | **54** | **54** |

## What was fixed

Sequence fields take pasted FASTA, numbers and line breaks (`clean_pasted`, new `clean_pasted_protein`); numbers outside a sensible
range are refused with a message (locus box, Primer options); every unbounded text field got a maximum length (search box,
feature name, primer name / sequence / note, new sequence name, NCBI fields, tab name); failing dialogs say what is wrong and
keep what was typed (File > New, Feature dialog, Primer dialog, Map Reads, Primer options); pairs of related spin boxes can no
longer cross (primer length, assembly overlap); the image export refuses a picture too large for memory; the Primer dialog's
tab order follows its rows. OK / Save is the default button of every `QDialogButtonBox` dialog (Qt makes the first accepting button the default, so Enter in
a one-line field presses it; multi-line boxes keep Enter for a new line). Not checked by key events in a real window.

## What was left

* The right-click menu **Paste** in a sequence field goes through the validator only: a pasted FASTA is refused whole there
  (Ctrl+V / Shift+Insert clean it). Qt offers no hook for the menu's paste without replacing the menu.
* The NCBI API key is stored as plain text in `ncbi.ini` in the data folder (it is hidden on screen). It is optional, and
  an encrypted store would be a design change.
* The Rust core says "location ... is outside the N bp sequence" for any location it cannot parse; the dialog now catches the
  common cases (blank, no digits, wrong characters) before it gets there, but `10..20x` still gets the core's wording.
* Free-text filters (feature table, primer database, digest) accept anything by design: they only match substrings.
* File name fields typed by hand (Output BAM, genome path) are checked for existence / folder, not for characters Windows
  forbids; the file dialogs (Browse) do that.
* The visual result of the new tooltips and inline messages was checked with the offscreen platform only, not on screen.
