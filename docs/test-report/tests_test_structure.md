# `tests/test_structure.py`

[← all test files](README.md) · source: [`tests/test_structure.py`](../../tests/test_structure.py)

## `test_sources_cache_auto_order_and_local_command`

**Result:** ✅ passed in 0.07 s · [line 71](../../tests/test_structure.py#L71) · **Area:** 3D structure sources (fake network) · **Checked against:** hand-computed · **Strength:** 🟢 strong

Each structure source (AlphaFold DB, PDB, ESMFold, local command) returns a structure and validates its input; results are cached; Auto tries AlphaFold, then PDB, then ESMFold and reports every failure.

| # | Input | Expected output |
|---|---|---|
| 1 | 46-aa PROTEIN; AlphaFold DB with 'p12345' | has pLDDT, source 'AlphaFold DB', detail contains 'AF-P1-F1', not cached |
| 2 | PDB; ESMFold (fake B-factor 0.5) | no pLDDT, detail '1ABC'; pLDDT, mean B-factor == 50.0 (0-1 scaled to 0-100) |
| 3 | AF without accession; 'MKWV' to PDB; 'A'*401 and 'MKWV1*B' to ESMFold; AF 'BADACC1' | StructureError: 'UniProt accession', '25 amino acids', 'up to 400', 'not amino acids', '404' |
| 4 | Auto with P12345, twice, then use_cache=False | AlphaFold first (URL ends /P12345); second cached with no network; third not cached, network used |
| 5 | Auto without accession, PDB returns no hits | source ESMFold; RCSB search called, AlphaFold not |
| 6 | fetch always raises | StructureError with 'No structure found', 'PDB:' and 'ESMFold:' |
| 7 | local command python fold.py {fasta} {pdb}; 'echo hi'; command writing nothing | source 'Local command', PDB has ATOM; 'needs {fasta}'; 'wrote no PDB' |

## `test_accession_page_and_script_tag_safety`

**Result:** ✅ passed in 0.01 s · [line 126](../../tests/test_structure.py#L126) · **Area:** 3D structure viewer page / security · **Checked against:** hand-computed · **Strength:** 🟢 strong

A UniProt accession is read from db_xref; the 3Dmol viewer page escapes the title and cannot be broken by '</script>' in the data or PDB; an AlphaFold reply pointing to a non-AlphaFold URL is not followed.

| # | Input | Expected output |
|---|---|---|
| 1 | accession_of([('label','x'),('db_xref','UniProtKB/Swiss-Prot:P69905')]); [('db_xref','GeneID:1')] | 'P69905'; '' |
| 2 | viewer_html with detail 'detail </script> here', title 'My <protein>' | contains '$3Dmol', 'ATOM', 'My &lt;protein&gt;'; exactly 2 '</script>' |
| 3 | AlphaFold reply with pdbUrl 'file:///C:/Windows/win.ini' | StructureError 'somewhere else'; only 1 fetch |
| 4 | PDB with 'REMARK </script><script>alert(1)</script>' | no '<script>alert(1)' in page; exactly 2 '</script>' |

## `test_window`

**Result:** ✅ passed in 0.44 s · [line 151](../../tests/test_structure.py#L151) · **Area:** GUI: 3D structure dialog · **Checked against:** GUI state · **Strength:** 🟡 adequate

The structure dialog disables the accession box for ESMFold and shows its privacy notice, fetches in the background and enables Open/Save, writes a viewer page, and shows a source error in its status.

| # | Input | Expected output |
|---|---|---|
| 1 | StructureDialog('tiny', PROTEIN, 'P12345') | Open disabled; NOTICE['ESMFold'] mentions ESM |
| 2 | source ESMFold, start, wait up to 5 s | accession disabled, notice 'sent to api.esmatlas.com'; Open and Save enabled, status mentions ESMFold; view_file() exists and contains '3Dmol' |
| 3 | source AlphaFold DB, accession 'BADACC1', start | status contains '404' |

**Review notes:**

- No docstring.
- Lines 159-160 and 168-169: wait loops end silently after 5 s; a slow worker gives an unexplained assertion failure.
