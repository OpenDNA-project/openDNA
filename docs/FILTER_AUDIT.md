# Filter by values: audit of the tables, and how to validate it

Filter by values (`app/value_filter.py`) lets you right-click a column header and tick which of its values are shown, like a
spreadsheet's filter. It hides rows; it never changes what a table holds, so a table's own code (which reads its rows by number)
is untouched. It is offered only on columns that hold **categories** (a handful of repeated values). Names, sequences, positions and
measured numbers are left to the text filter, the sort and the ranges the tools already have.

## Audit

| Table (file) | Rows | Filter on | Not on, and why |
|---|---|---|---|
| Features view (`feature_table.py`) | features of the file | Type, Strand | Name (unique), Location, Length (numbers) |
| Features Database (`db_tabs.py`, the same view) | every feature of the database | Type, Strand, Source (by the source's name; the field also holds the record) | Name, Length |
| Primers view (`primer_panel.py`) | primers binding here / every primer of the database | Match, Structure, Added | Name, Sequence, Length, Tm (numbers); Note (free text: the text filter) |
| Detect Common Features (`annotation.py`) | parts found | Type, Source | Identity, Strength (numbers), Location |
| Detect ORFs (`tools.py`) | ORFs | Strand, Annotation | Location, Length |
| Restriction Digest fragments (`tools.py`) | fragments | Ends | Size, From, To (numbers) |
| Spacer Finder (`cas9.py`) | guides | Strand, Here, Genome, Verdict | Position, guide, cut, GC, Score |
| Local Align (`local_align.py`) | alignment pieces | Strand | identity, score, positions |
| Plasmid Checkup (`checkup.py`) | elements; origins | Type, From; Replicates in, Copies | names, locations, free-text details |
| NCBI search (`ncbi_dialog.py`) | records found | Organism | Accession, Title, Length, Updated |
| NCBI reference genomes (`ncbi_dialog.py`) | assemblies | Organism, Assembly, Level | Accession, Size, Sequences, Released |
| Primer Design candidates (`primer_design.py`) | primer pairs | (none) | every column is a sequence or a number; sorting does the job |
| Checkup units (tree), Assembly report, Database Status, Shortcuts, Performance, Self-test, History, Welcome lists | | (none) | short fixed lists, or not data to sift |

## Plan for validation

Automated (`tests/test_db_tabs.py`): a table and a tree filtered; a refill keeps the filter and shows a value never seen; the
header mark; the search box with All / None; the Features Database filtered by Source; the filter surviving the database changing.

By hand, once per table of the audit (about 15 minutes in all):

1. Fill the table (open `demo/` files: a plasmid for the Features view, Detect Features, ORFs, Digest, Checkup; the E. coli genome
   for the Spacer Finder with the reference set; an NCBI search for the two NCBI tables).
2. Right-click each column named in the audit: the popup lists its values with their counts, the most useful first; a column not
   named shows no popup.
3. Untick one value, OK: those rows disappear, the header ends with ▼, the other rows keep their order and sort.
4. Click a remaining row: the views select what that row is (the row numbers are unchanged by the filter).
5. Make the table refill (edit the sequence, change a setting of the tool, run the search again): the filter is still on; a value
   that was not there before shows.
6. Filter a second column: rows must pass both. *Clear filter* on one column brings its rows back, the other column still filters.
7. Type in the popup's search box, then *None*: only the listed values are unticked.
8. Light and dark theme: the ▼ and the popup are readable in both.
9. A large file (the 4.6 Mb genome, ~10,000 features): filtering the Features view by Type takes well under a second.

Anything that fails: note the table, the column and the step.
