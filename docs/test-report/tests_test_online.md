# `tests/test_online.py`

[← all test files](README.md) · source: [`tests/test_online.py`](../../tests/test_online.py)

## `test_ncbi_client_search_cache_rate_limit_and_errors`

**Result:** ✅ passed in 6.06 s · [line 153](../../tests/test_online.py#L153) · **Area:** NCBI download · **Checked against:** hand-computed · **Strength:** 🟡 adequate

**Why it matters:** Without it fetching sequences from NCBI could break, send wrong parameters, ignore NCBI's request rate limit (getting the user blocked), cache error pages as sequence files or fail on a temporary 429.

With a fake E-utilities endpoint, checks search results and request parameters, empty and invalid searches, error propagation, the UniProt accession on protein records, the on-disk cache, the size limit, protein FASTA reading, the 3-requests-per-second limit without a key, rejection of non-record answers, and a 429 retry.

| # | Input | Expected output |
|---|---|---|
| 1 | search('nuccore','U1.1', key='K') | One record ('U1.1', 1834, 'circular', 'A made-up plasmid'); request has tool=openDNA, api_key=K, db=nuccore, no email |
| 2 | search 'nothing' / '  ' / 'broken' | [] / NcbiError 'Type an accession' / NcbiError '503' |
| 3 | search('protein','P12345') | uniprot 'P12345'; nucleotide record has uniprot '' |
| 4 | get(record) twice | First: .gb file starting LOCUS, cached False; second: same path, cached True, no request |
| 5 | Record with MAX_BASES+1 bases | NcbiError 'too long' |
| 6 | get a protein record then protein_of(path) | ('P12345.1', 'MKWHAA') |
| 7 | 3 calls without a key / with a key | >= 0.6 s elapsed / < 0.6 s |
| 8 | fetch returning 'Error: nothing here' | NcbiError 'did not return'; nothing cached |
| 9 | fetch raising 429 then returning 'fine' (sleep faked) | call returns 'fine' |

**Review notes:**

- The rate-limit check uses real wall-clock time (adds about 0.6 s; the '< 0.6 s with a key' branch could be flaky on a very slow machine).
- The 429 retry checks a single retry only; the retry limit and giving up are not checked.

## `test_structure_sources_cache_and_auto_order`

**Result:** ✅ passed in 0.10 s · [line 208](../../tests/test_online.py#L208) · **Area:** 3D structure: sources and cache · **Checked against:** hand-computed · **Strength:** 🟡 adequate

**Why it matters:** Without it the protein 3D view could pick the wrong source, accept invalid input, ignore its cache, hide why every source failed, or break the user's own local folding command.

With a fake network, checks each structure source (AlphaFold DB, PDB, ESMFold) and its input checks, the cache, Auto's order, Auto's combined error, and the local command source using a small Python script.

| # | Input | Expected output |
|---|---|---|
| 1 | get(PROTEIN, 'AlphaFold DB', 'p12345') | pLDDT model, source 'AlphaFold DB', detail has AF-P1-F1, not cached |
| 2 | get(PROTEIN, 'PDB') | No pLDDT; detail mentions 1ABC |
| 3 | get(PROTEIN, 'ESMFold') with B-factors 0.5 | pLDDT; mean B-factor 50.0 (0-1 scaled to 0-100) |
| 4 | AlphaFold with no accession; 'MKWV' to PDB; 401 aa to ESMFold; 'MKWV1*B' to ESMFold; accession BADACC1 | StructureError 'UniProt accession' / '25 amino acids' / 'up to 400' / 'not amino acids' / '404' |
| 5 | Auto with accession P12345, clean cache | AlphaFold DB; first request ends /P12345 |
| 6 | Auto again; then use_cache=False | cached, no network / not cached, network used |
| 7 | Auto without accession, PDB search returns nothing | ESMFold; RCSB search was tried; AlphaFold never asked |
| 8 | All fetches raising 'no network' | StructureError with 'No structure found', 'PDB:' and 'ESMFold:' |
| 9 | Local command running a Python script with {fasta} {pdb} | source 'Local command'; PDB text contains ATOM |
| 10 | Command 'echo hi' / a script that writes nothing | StructureError 'needs {fasta}' / 'wrote no PDB' |

**Review notes:**

- Checks on the returned models are mostly presence checks (plddt truthy, 'ATOM' in text, id in detail) rather than content; the ESMFold 0-1 to 0-100 scaling is the one exact value.

## `test_structure_page_script_safety_and_alphafold_reply_not_followed`

**Result:** ✅ passed in 0.00 s · [line 264](../../tests/test_online.py#L264) · **Area:** 3D structure: viewer page safety · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** Without it a crafted protein name, PDB file or AlphaFold reply could inject script into the 3D viewer page or make the app fetch a local or foreign file.

Checks UniProt accession extraction from feature qualifiers, HTML escaping of names, that a '</script>' inside the details or PDB text cannot close the page script, and that an AlphaFold reply naming a non-AlphaFold URL is refused without fetching it.

| # | Input | Expected output |
|---|---|---|
| 1 | accession_of with db_xref 'UniProtKB/Swiss-Prot:P69905' / only 'GeneID:1' | 'P69905' / '' |
| 2 | viewer_html with detail 'detail </script> here' and name 'My <protein>' | Contains $3Dmol and ATOM; name escaped as 'My &lt;protein&gt;'; exactly 2 '</script>' |
| 3 | AlphaFold reply with pdbUrl 'file:///C:/Windows/win.ini' | StructureError 'somewhere else'; only one fetch (the API) |
| 4 | PDB text with 'REMARK </script><script>alert(1)</script>' | No '<script>alert(1)'; exactly 2 '</script>' |

**Review notes:**

- The accession_of check is unrelated to page safety; it would fit better with the structure-source test.

**Fixed by the audit of November 2026:**

- 13 Nov 2026: gap closed. The proposed checks were added, pass on the current code and fail under the mutation that slipped through before. Was: 'an AlphaFold reply naming a non-AlphaFold URL is refused' is tested with file:/// only; a check that just requires https (or any scheme check) passes, so a reply could make the app fetch from any https host.

## `test_pdb_entries_and_the_model_read_from_a_pdb_text`

**Result:** ✅ passed in 0.00 s · [line 293](../../tests/test_online.py#L293) · **Area:** 3D structure: PDB entries and model parsing · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** Without it the PDB search list could show duplicates or wrong identity/coverage, linked entries could be badly ordered, and the 3D view could mis-read chains, ligands or secondary structure, or map protein residues to the wrong atoms.

Checks the RCSB search result list and entry details, PDB links from feature qualifiers, the similarity measure, ordering of linked entries by their best chain, and the PDB text parser: first model only, ligands without water, helix/sheet/loop from alpha-carbon geometry, and mapping protein positions to residues on chains with offset numbering, a missing residue and a substitution.

| # | Input | Expected output |
|---|---|---|
| 1 | USER_AGENT | Starts 'openDNA/' and includes a version |
| 2 | search_pdb(PROTEIN) with hits 1ABC_1, 2XYZ_2, 1ABC_2, '../etc_1' | Entries ['1ABC','2XYZ']; first = (0.95 identity, 0.5 covered, 'Crambin', 'Crambe hispanica', 'A, B', 1.5); 2XYZ identity 0.8 from its score, empty title |
| 3 | pdb_entry(first) | Source PDB; detail has 'PDB 1ABC: Crambin', '95% identity', 'chains A, B'; info dict exact (method 'X-ray diffraction', resolution '1.50 Å', link https://www.rcsb.org/structure/1ABC, ...) |
| 4 | search_pdb('MKWV') | StructureError '25 amino acids' |
| 5 | pdb_links on PDB:1abc, PDB:1ABC, a UniProt xref, PDB:2xyz_A, PDB:bad, note PDB:9ZZZ | ['1ABC', '2XYZ'] |
| 6 | similarity('MKWHLLEGRS','KWHLLEGRS'); with '' | (1.0, 0.9); (0.0, 0.0) |
| 7 | linked_pdb(PROTEIN, ['1ABC'], 'P12345') | ['3DEF','1ABC']; 3DEF identity 1.0, coverage (len-10)/len, chains 'B, D', organism Crambe hispanica; 1ABC identity < 0.5 |
| 8 | linked_pdb(PROTEIN) / accession 'P1' | [] / StructureError '6 to 10' |
| 9 | Synthetic PDB: ideal helix chain A, two adjacent strands B and C, lone strand D, HEM FE and water HETATM, ENDMDL then chain E | Ligands [('FE',(0,0,0),'HEM A900')]; ss A='H'*12, B and C 'EEEEEE', D '------', E absent |
| 10 | residue_map('MKWHLLEGRS') on chain A numbered from 10 missing W with L->A, chain B whole, chain C 'GGGG' | 10 positions; pos0 -> A10,B1; pos2 -> B3 only; pos5 -> A14,B6; nothing on C |

**Review notes:**

- Two unrelated topics (RCSB search/linking and PDB text parsing) in one test; could be split.
