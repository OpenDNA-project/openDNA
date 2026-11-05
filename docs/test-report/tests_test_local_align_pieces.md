# `tests/test_local_align_pieces.py`

[← all test files](README.md) · source: [`tests/test_local_align_pieces.py`](../../tests/test_local_align_pieces.py)

## `test_query_text_and_junction_texts`

**Result:** ✅ passed in 0.04 s · [line 73](../../tests/test_local_align_pieces.py#L73) · **Area:** Local Align query cleaning + junction/summary text · **Checked against:** in-test reference model · **Strength:** 🟢 strong

clean_query normalises FASTA/lowercase/RNA/digits/gaps and gives clear errors (empty, header-only, no DNA, protein, too-long), and junction_text/summarize describe insertions, deletions, inversions, duplications and the summary line correctly.

| # | Input | Expected output |
|---|---|---|
| 1 | clean_query on FASTA, ORIGIN block, RNA, gap dashes, IUPAC, multi-record, CRLF | '>seq1...\nacgt...'->'ACGTACGTGGCC'; ORIGIN->'ACGTACGTGGCCTT'; 'ACGU UGCA'->'ACGTTGCA'; 'AC-GT--AC'->'ACGTAC'; IUPAC kept; only first record; CRLF handled |
| 2 | clean_query errors: '', whitespace, header-only, '123 456', protein, 'ACGT hello 12', length MAX_QUERY+1 and exactly MAX_QUERY | QueryError matching empty/empty/'only a FASTA header'/'No DNA'/'protein'/'Not DNA'/'too long...200,000'; exactly MAX_QUERY accepted |
| 3 | junction_text(a=hit(0,100,1000,1100), various next hits, n=5000, circular) | contiguous; 'insertion in the query: 200 bp unmatched'; 'deletion from the query: 50 bp of the template skipped'; 'template jumps back'; 'inversion breakpoint'; reverse-reverse contiguous; circular short-way -> 'deletion from the query: 50 bp' |
| 4 | summarize([a],n) and a 3-piece list and [] | '1 piece: 1 forward'; '3 pieces: 2 forward, 1 inverted...; gaps...200 bp unmatched (#1→#2)' with 'insertion in the query (#1→#2)' and 'inversion breakpoint (#2→#3)'; '' for empty |

## `test_structural_variant_query`

**Result:** ✅ passed in 0.41 s · [line 118](../../tests/test_local_align_pieces.py#L118) · **Area:** GUI: Local Align structural variant table · **Checked against:** hand-computed · **Strength:** 🟢 strong

A structural-variant query produces one table row per piece in query order with correct coordinates, strand words, identity, 'before this piece' junction text, summary line, inverted-row colour, row-click selection and alignment bars; and a circular template typed reverse-complemented returns the same four places.

| # | Input | Expected output |
|---|---|---|
| 1 | LocalAlignDialog on T(2000); query A+200+B+rc(C)+D; run | hits (qstart,qend,strand)=[(0,100,1),(300,400,1),(400,500,-1),(500,600,1)]; template [(100,200),(600,700),(1100,1200),(1600,1700)]; table 4 rows; query col '1–100'..'501–600'; template col with thousands commas; strand col FORWARD/FORWARD/INVERTED/FORWARD; all '100.0%'; status '4 hits' |
| 2 | junction cells, summary line, inverted colour, row clicks, alignment bars | row0 'start of the query'; row1 'unmatched query...200 bp...400 bp of the template skipped'; rows2-3 'inversion breakpoint'; summary '4 pieces: 3 forward, 1 inverted...'; inverted row background==INVERTED_COLOR, forward row not; clicking row selects template span (1100,1200)/(1600,1700); genome_view bars set {4 spans}, inverted labelled; clear() empties table/summary |
| 3 | circular T, query rc(QUERY) | hits (start,end,strand)=[(1600,1700,-1),(1100,1200,1),(600,700,-1),(100,200,-1)]; qstart/qend mirrored; summary '4 pieces: 1 forward, 3 inverted...' |

**Review notes:**

- errors fixture (asserting no uncaught Qt-callback exceptions) is good; no issues.

## `test_query_box_and_clipboard`

**Result:** ✅ passed in 0.55 s · [line 176](../../tests/test_local_align_pieces.py#L176) · **Area:** GUI: Local Align clipboard + query box · **Checked against:** GUI state · **Strength:** 🟢 strong

Pasting from the clipboard cleans text and explains when there is no DNA (leaving the query box alone), the query box accepts lowercase/RNA and rejects protein/empty/huge with a message, and the 'Align Clipboard' menu entry opens the dialog with the clipboard aligned.

| # | Input | Expected output |
|---|---|---|
| 1 | paste_clipboard with empty/whitespace/prose/header-only/protein; then a valid ORIGIN-style FASTA of T[:60] | each bad case returns False, message in status ('clipboard is empty'/'no DNA'/'only a FASTA header'/'protein'), query box unchanged ('keep me'); valid -> True, box==T[:60], status 'Query: 60 bp from the clipboard' |
| 2 | run() on T[300:400] lowercase with U-for-T; then protein/empty/too-long | RNA found as hit (300,400,1); each bad input -> message in status, hits [] and table empty |
| 3 | 'Align Clipboard…' QAction trigger with clipboard FASTA then empty | dialog shown (exec monkeypatched), box==QUERY, 4 hits strands [1,1,-1,1], 4 rows; empty clipboard -> dialog still opens, 'clipboard is empty', hits []; Tools menu offers both '&Local Align…' and 'Align &Clipboard…' |
