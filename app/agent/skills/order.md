# order: a primer list as a supplier's order file

**Use when:** the user wants primers ready to order (IDT plate upload, or plain name and sequence). Same as Primers > Export order.

## Input: a CSV with exactly two columns

```
id,sequence
KLD-F1,ATGCGTACGTTAGCAGCT
KLD-R1,TTAGCGGCATCGATCGAA
```

`id` is the primer's name. Take the sequences from a `primers.csv` or `report.jsonl` of an earlier command and write one row per primer. Every primer is checked as the Primers panel checks one (name, A/C/G/T/IUPAC letters, 8 bases at least); one problem and nothing is written.

## Command

```
order PRIMERS.csv NEW_FILE [--format idt|plain] [--scale 25nm] [--purification STD]
```

`idt` (default): `Well Position,Name,Sequence,Scale,Purification`, 96 primers per file (`NEW_FILE_plate2.csv` ... for more). `plain`: name,sequence. Scales: 25nm, 100nm, 250nm, 1um, 2um, 5um, 10um; purifications: STD, PAGE, HPLC, IEHPLC, RNASE, PAGE_HPLC. The user decides the scale and purification: **ask, do not assume.** `NEW_FILE` must not exist.

## Example

The example writes a file into the current folder.

```
order app/agent/examples/order_primers.csv example_order.csv
```
