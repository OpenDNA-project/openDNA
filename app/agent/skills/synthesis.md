# synthesis: will a synthesis company make this fragment quickly?

**Use when:** the user asks whether a piece of DNA can be ordered as a synthetic gene or gBlock, how hard it is to synthesise, or whether a region has sequence that makes synthesis fail (long repeats, extreme GC, strong hairpins). Same as the Local row of the Synthesis tab. Read-only; nothing leaves the computer.

## Command

```
synthesis A.gb [--region START..END]
```

Runs the Synthesis Success Calculator (Halper, Hossain & Salis, ACS Synth Biol 2020, doi:10.1021/acssynbio.9b00460), a random forest trained on 1,076 fragments that real companies tried to make. It needs **ViennaRNA** (`pip install ViennaRNA`; openDNA does not bundle it): without it the command stops with one ERROR line saying so. Only A, C, G and T are accepted (the first other base is named). The model reads a straight piece, never a circle: for a circular file say so (a warning does). It takes about 10 seconds per 1,000 bases: test a region (`--region`) rather than a whole large plasmid.

JSON: `citation`, `bases`, `region` (1-based; an `end` above the length means the region runs through the origin of a circular file, and is read as one straight piece), `predicted_quick_synthesis` (true or false), `probability` (of quick synthesis, 0 to 1) and `features`, the nine numbers the forest reads: `size`, `total_GC` (fraction), `longest_repeat` (0 when no repeat is over 11 bases), `repeat_9_metric`, `richest_hairpin` (GC fraction of the richest hairpin stem), `GC_short_l` (20-base windows with 20% GC or less), `GC_term_l` (low-GC windows in the last and first 50 bases), `Tm_low` and `Tm_high` (20-base windows melting at 40 °C or below, 70 °C or above). `places`: for each of these, the stretches of the file it counted (`start`..`end`, 1-based, the windows merged; empty for `size` and `total_GC`), to say where the problem is.

## What to tell the user

The verdict and the probability, as the calculator's prediction and not a quote: a company's own check decides. Then which of the features is out of the ordinary (a long repeat, many low-GC windows, a very GC-rich hairpin), and where (`places`), since those are what to fix. Always give the citation.

## Example

```
synthesis demo/pBAD30_eilR.gb --region 1..300
```
