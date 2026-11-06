# `tests/test_ncbi.py`

[← all test files](README.md) · source: [`tests/test_ncbi.py`](../../tests/test_ncbi.py)

## `test_ncbi_calls`

**Result:** ✅ passed in 6.06 s · [line 66](../../tests/test_ncbi.py#L66) · **Area:** NCBI import (fake network) · **Checked against:** hand-computed · **Strength:** 🟢 strong

NCBI search parses summaries and sends the tool name / API key, downloads are cached and size-checked, requests are rate-limited (3/s without key), an invalid answer is not cached, and HTTP 429 is retried.

| # | Input | Expected output |
|---|---|---|
| 1 | search('nuccore', 'U1.1', key='K') with fake esearch/esummary | [('U1.1', 1834, 'circular', 'A made-up plasmid')]; request has tool=openDNA, api_key=K, db=nuccore, no email |
| 2 | terms 'nothing', '  ', 'broken' | []; NcbiError 'Type an accession'; NcbiError '503' |
| 3 | protein search 'P12345' | uniprot 'P12345'; nucleotide record uniprot '' |
| 4 | get(record) twice | .gb file starting 'LOCUS', not cached; second call (path, True) with no request |
| 5 | Record of MAX_BASES+1 bp; protein get | NcbiError 'too long'; protein_of == ('P12345.1', 'MKWHAA') |
| 6 | 3 calls without key; 3 calls with key (real clock) | >= 0.6 s; < 0.6 s |
| 7 | fetch returns 'Error: nothing here' | NcbiError 'did not return', no cache file |
| 8 | fetch raises 429 then returns 'fine' (sleep patched) | call() returns 'fine' |

**Review notes:**

- Lines 93-100: the rate limit is measured with the real clock and real sleeps (0.4 s gap); the `< 0.6` bound with a key could fail on a loaded machine, and the sleeps make this the slowest test in the group (~6 s).

## `test_dialog`

**Result:** ✅ passed in 3.75 s · [line 127](../../tests/test_ncbi.py#L127) · **Area:** GUI: NCBI import dialog · **Checked against:** GUI state · **Strength:** 🟡 adequate

The NCBI dialog shows search results, opens a record in the current empty window and the next one in a new window, sends a protein to the 3D-structure window, and reports 'nothing found' and service errors.

| # | Input | Expected output |
|---|---|---|
| 1 | query 'U1.1', search, wait | 1 row: 'U1.1', '1,834 bp'; Open button enabled |
| 2 | Open | opened in win; untitled, length 1834, has features |
| 3 | second dialog, same record, Open | opened in another window, length 1834 |
| 4 | protein database, 'P12345', Open (structure.get patched) | button '3D structure…'; window protein 'MKWHAA', accession 'P12345', name 'RecName: Full=Tiny' |
| 5 | query 'nothing'; query 'broken' | status 'Nothing found', Open disabled; status contains '503', Find re-enabled |

**Review notes:**

- Lines 120-124: wait() polls up to 5 s and returns silently on timeout; the next assert then fails without saying the worker never finished (flaky under load).
