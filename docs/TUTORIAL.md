# Tutorial: from a sequence to an assembly

About 15 minutes. You will open a sequence, edit it and watch the annotations follow, undo and redo, design
primers, then assemble three fragments three different ways. Everything uses small made-up demo files, so
nothing here touches your own data.

## 0. Setup

Follow *Quick start* in the [README](../README.md), then make the demo fragments (they are also committed in `demo/`):

```powershell
python scripts/make_demo_files.py
python app/main.py demo/demo_insert.gb
```

| File | What it is |
|---|---|
| `demo/demo_vector.gb` | A circular 1,834 bp plasmid: an origin, `AmpR`, and a 120 bp `dropout` between two BsaI sites that point outwards (so BsaI releases the dropout and leaves the other 1,700 bp as the backbone). |
| `demo/demo_insert.gb` | 700 bp: a promoter (`pDemo`) and a `reporter` CDS. |
| `demo/demo_terminator.gb` | 220 bp: a terminator. |

The sequences are random, so the "genes" do nothing biologically; they exist to have something to move around.

## 1. Look around

`demo_insert.gb` opens in the sequence view. The two features are drawn as arrows over the bases; the map and the
feature table show the same thing. Click a feature to select it, or drag across bases to select a range
(the selection has draggable `[` `]` ends). *Help → Keyboard Shortcuts* (F1) lists every shortcut.

## 2. Edit and watch the annotations follow

1. Click between two bases **inside** the `reporter` CDS and type a few bases (or paste some). The arrow grows
   by exactly what you added, and every base of the CDS is still the same base.
2. Click **before** the promoter and type some more. Both features slide to the right, unchanged.
3. Delete a stretch that overlaps the end of `pDemo`: the feature shrinks to what is left of it. Delete all of
   a feature and it disappears (and comes back with undo).
4. Open **History** (toolbar button, Ctrl+H). Each change is described in words ("Inserted 5 bp at 301: …" or "Deleted 3 bp at 40–43: …");
   typing is grouped. Click an earlier row to go back to before it, or click *Current sequence* to return.
5. `Ctrl+Z` / `Ctrl+Y` undo and redo one step at a time.

> The rule for inserts, so nothing surprises you: bases typed **strictly inside** a feature join it; bases typed
> at its very first or last base stay outside. It is tested hundreds of times at random: see
> [TESTING.md](TESTING.md) section 2.

## 3. Add a feature

Select a range, then *Feature → New Feature from Selection…* (Ctrl+T), choose a type, give it a label. It is an
undo step like any edit. *Feature → Edit Feature…* (Ctrl+E) changes it later.

## 4. A plasmid: break the circle

Open `demo/demo_vector.gb` (File → Open). It is circular, so the sequence view shows it twice with a red bar at the origin.

1. *Edit → Circular Sequence* (or the toolbar's **Circular** button) and unpress it.
2. OpenGenEdit asks where to break the circle: that base becomes base 1. Pick one inside a feature to see the
   warning that it will be cut into a two-part feature, and answer *Yes*.
3. The sequence is now linear and rotated. In **History** this is **one** step, and one undo gives you back the
   circle with its original base 1.

## 5. Primers

- *Primers → New Primer…* (Ctrl+Shift+P): type or paste a sequence; binding sites on both strands appear, with
  mismatches, Tm, hairpin and self-dimer checks.
- *Primers → Design Primers for Selection…*: select a stretch first.
- Primers you keep go in the database (*Databases → Primer Database…*, Ctrl+Shift+D, a tab of the window), a `primers.tsv` file you can export.

## 6. The Assembly panel

*Tools → Assembly…* opens the **Assembly tab**, next to your other tabs (Edit, Primers, Tools). It is a tool tab, not a pane: it is
not in the pane drop-down and cannot be split. Close it with its ✕ (its content is kept until you open it again). It needs no file to be open.

The top half is the setup and the map, the bottom half is the same sequence view as in the editor (read-only) on the file behind the selected fragment.

1. **Method:** *Gibson / HiFi*, *Golden Gate* (then pick the enzyme) or *Sidewinder*, and the **number of fragments** (1 to 12).
   Each fragment appears as a gray box, *Fragment 1*, *Fragment 2*…
2. **Click a gray box.** A menu asks where the DNA comes from: *Currently open files*, *Recently opened files* or *Browse…*.
   Choose `demo_vector.gb` for Fragment 1, `demo_insert.gb` for Fragment 2 and `demo_terminator.gb` for Fragment 3.
   A file up to 200 kb is used whole straight away; the box turns into a coloured arrow with a mini map of its features.
3. **Pick what to use.** Click a box to see its file in the bottom half. Drag over the bases (or click a feature) to select a
   stretch and press **Use as fragment**. To amplify it instead, click a primer in the view (the primers of your database and of the file are
   drawn, with their tails) and shift-click the one on the other strand, as in the editor: that selects what they amplify, switches to
   *PCR with two primers* and fills the two lists. Press **Use as fragment**. Any 5′ part of a primer that does not match the template becomes a tail,
   and the fragment is the PCR product, tail included.
4. **Watch the junctions.** The map shows small fragment boxes with a **pill between each pair** (and, for a circular assembly, a last pill
   and an arrow that runs back over the boxes and enters Fragment 1: the circle is closed where you read the fragments). As soon as every box
   is filled each pill says whether the two ends fit: a red pill `✗ ends do not fit` means they do not. For the demo files they do not yet,
   because plain stretches carry no homology. Beside the map, the **product** is drawn as a ring (circle) or a bar (line) in the fragments'
   colours, with a dot on every junction. Click a pill for that junction in full detail, or **Full report…** for every junction, check and oligo.

## 7. Assembly 1: Gibson

With *Gibson / HiFi* and the three demo files in order (vector, insert, terminator), press **Design primers**. OpenGenEdit designs
new primers that carry the neighbour's end as a 5′ tail (any primers you chose earlier are replaced), and re-checks.

What you should see (the files are seeded, so the numbers are exact):

- Status line: **✓ Every junction checks out: 2,754 bp circular product, 6 oligos**. Each box now says `+20` (the tail added).
- Each pill is green and reads, for example, `✓ 20 bp` over `57 °C · 45 %`. Click it and the **junction inspector** opens: the product around
  the joint in the sequence view, on the top strand alone, with the 20 bp joint selected. Each fragment is a bar in its colour, the joint is a
  slate bar and what the primers' tails added is orange (here two bars of 10 nt: each primer brings half of the overlap, which a sentence above
  the sequence says in words); the two primers at the junction are drawn as primers, and every base, primer and feature answers to the pointer
  as in the editor. Under the sequence are the checks of that junction: length, Tm, GC, **unique** (the homology must occur once in the whole
  product; twice is a red error, because the ends could join in the wrong place), hairpin and distinct (no 10-mer shared with another junction).
  A failed check is named in red or orange on its pill; hover a pill for every reason.
- Select a box (click it): the bottom half shows its file with the **designed primers drawn on it** like any primer of the editor (a line along the
  binding site, its tail on a second line with the tail bases in orange), and under the sequence you read each oligo: the tail in lowercase orange,
  the annealing part in capitals, length, Tm and a mark.

**Now rearrange it.** Drag the terminator in front of the insert, or double-click the insert to reverse-complement it: the junctions turn red at
once, because the ends no longer fit. Press **Design primers** again and they are re-made for the new order. Rotating the circle (moving the
last fragment to the front) changes nothing, and the panel agrees. **Perform assembly** (enabled only when every junction checks out) opens the
exact assembled plasmid in a new window with each fragment's features in place. **Save fragments…** writes each PCR product, tails included, to its own
GenBank file, and **Save primers to database** keeps the oligos.

## 8. Assembly 2: Golden Gate

1. Set **Method** to *Golden Gate* and **Enzyme** to *BsaI*.
2. Right-click Fragment 1 (`demo_vector`) and tick **Destination vector**: its box gets a dashed outline and `[vector]`. It will be cut at its two BsaI
   sites and the 1,695 bp outside them kept (the sites must face away from the backbone, as in a real destination vector).
3. Press **Design primers** (for Golden Gate it adds the site, spacer and overhang).

Expect **✓ 2,620 bp circular**, **4 oligos** (the vector is cut, not amplified) and three 4 nt overhangs, for example `AACG`, `GAGG` and `ATCC`.
Each card shows the two ends with the sites and the cut marked, for example `…TTTCCT|GAGG|tgagaccaaa…` over `…tttggtctca|GAGG|GTAAGG…` (lowercase = the
site and spacer the enzyme removes), and these checks: **ends fit** (both sides carry the same overhang), **unique**, **not palindromic**, **orthogonal**
(not the reverse complement of another overhang, not one base apart) and **GC mix**, and **validated** (one of the 25 overhangs
Potapov et al. (2018) built a 24-fragment assembly from). An overhang outside that set is only information, not an error: the paper's full ligation
table is not included, so the check is a rule of thumb (see *Help → Citations*). A fragment with an extra BsaI site gets a red mark on its box
(that is what "domesticate" means). The insert's primers carry `ttt` + the BsaI site + a spacer + the overhang; they are drawn on its sequence like the Gibson ones.

Try breaking it: reverse-complement the insert (its sites now face the wrong way, so the junctions go red), or choose BsmBI (the vector has
no BsmBI sites, and the panel says so).

## 9. Assembly 3: Sidewinder

Sidewinder builds each fragment from two annealed oligos of at most 165 nt, so it wants **short** fragments (about
60-150 bp). Choose *Sidewinder*, set the number of fragments to 3 and give each box a short stretch (select ~100 bp in a file and press
**Use as fragment**). Sidewinder always gives a linear product, and its toeholds come from the product itself, so the junctions always "fit";
each card shows the toehold on both strands, its barcode and the checks (GC, unique, distinct from the other toeholds); the 165 nt limit is in the status line. Select a box to read its two oligos under the sequence. The layout is a simplified reading of the method: check it against the paper
before ordering.

## 10. Check the tool itself

*Help → Troubleshooting → Self-Test… → Run all* runs the built-in checks (edits, annotation drift, primer binding)
and shows ✓/✗ with times. If any fails, please report it with the message shown.

## Where next

- [USER_GUIDE.md](USER_GUIDE.md): everything else (genome view and tracks, read mapping, variants, Local Align, the Python API).
- [TESTING.md](TESTING.md): what is tested and what is not.
