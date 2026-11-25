# print: a region of the tracks as an image

**Use when:** the user wants a figure of their tracks over a region (a gene with its RNA-seq, a deletion with its coverage): the Alignment view as an SVG (vector, for a paper) or PNG file. The same image as the window's right-click > Export This Region as Image.

## 1. From a project

```
print PROJECT.odp START..END OUT.svg
```

The project's tracks are drawn as it shows them: the tracks shown and hidden, their order, groups, colours, heights, heatmaps and highlighted regions. `START..END` is 1-based and inclusive, or a feature's name. `OUT` ends in `.svg`, `.png` or `.tif`, and must not exist (nothing is overwritten). `--width` sets the width in pixels (default 1600); the height is what the tracks need.

## 2. From a sequence and track files

```
print SEQUENCE.gb START..END OUT.png --track a.bed --track b.bw
```

Each `--track` file is drawn with the default look (bars for a signal, one colour per annotation track). For another look (heatmaps, groups, colours), set it in the window and save a project, then print the project.

## What to tell the user

- Where the image is and what it shows: the region, `tracks_shown`, and `regions` (highlights drawn).
- Any warning: a track not loaded (its sequence names do not match this sequence), a project file not found.
- A region across the origin of a circular sequence must be printed in two parts.

## Example

```
print demo/pBAD30_eilR.gb 1501..3500 example_region.png --track demo/hub/pBAD30_eilR/plus.bw --track demo/hub/pBAD30_eilR/sites.bb
```

## Errors you may see

| Message | Fix |
|---|---|
| `... already exists` | choose a new file name |
| `the file name must end in .svg, .png or .tif` | change the extension |
| `runs through the origin` | print START..LENGTH and 1..END one after the other |
