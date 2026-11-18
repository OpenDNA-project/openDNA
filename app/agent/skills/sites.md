# sites: restriction sites and digests

**Use when:** the user asks where an enzyme cuts, which enzymes cut a plasmid once (to linearise or clone), or what fragments a digest gives. Same lists as the Restriction view and the Digest tool. Read-only: nothing is written, so there is no `check`/`run`.

## Command

```
sites PLASMID.gb [--enzymes A,B] [--min-sites N] [--max-sites N] [--digest A,B]
```

- No option: the enzymes that cut **once** and recognise 6 or more bases (the window's default list). `--min-sites`/`--max-sites` widen it (`--max-sites 3`: up to three sites; `--min-sites 2` with no max: every enzyme that cuts twice or more).
- `--enzymes EcoRI,BsaI`: only these, whatever their number of sites (`count` 0 = no site). The names are those openDNA knows (an unknown name lists them).
- `--digest EcoRI,NcoI`: also the fragments when the plasmid is cut by all of them, largest first, with the enzymes at each end. A circular plasmid cut once gives one fragment the length of the plasmid (its `start` is after its `end`: it runs across the origin).

JSON: `enzymes` with, for each, `site`, `notation` (`G^AATTC`; a Type IIS enzyme is `GGTCTC(1/5)`), the `ends` it leaves, `count`, `cuts` (`top_after_base`: the cut falls after that base, 0 = before the first base; `bottom_after_base` differs for sticky ends) and the recognition `sites` (1-based, strand). `digest.fragments` when asked.

## What to tell the user

The enzymes and where they cut, in plain words; for a digest the fragment sizes. An enzyme that cuts more than once is not suitable to linearise. Enzymes outside openDNA's list (about 50 common ones) are not searched: say so if the user names one it does not know.

## Errors you may see

| Message | Fix |
|---|---|
| `No enzyme called ...` | take a name from the list in the message |

## Example

```
sites demo/pBAD30_eilR.gb
sites demo/pBAD30_eilR.gb --enzymes BamHI,XbaI --digest BamHI,XbaI
```
