# hub: the tracks of a UCSC track hub, listed or copied

**Use when:** the user has a UCSC track hub (the address of its `hub.txt`, a `hub.txt`, or a folder holding one) and wants to know what is in it, or to have its files (bigWig signals, bigBed annotations) on their computer to open in openDNA. The same reading and copying as the window's Import > Features and Tracks > Track Hub, which then loads them.

## 1. List the tracks

```
hub list HUB
```

`HUB` is a URL such as `https://example.org/myhub/hub.txt`, a local `hub.txt`, or the folder holding it. The answer is JSON: `hub` (name, label, `genomes`) and one record per track with data: `name`, `label` (the hub's short label), `group` (the composite or super track it is in, else the hub), `type`, `hidden` (the hub hides it at first), `color`, `url`, and `read` (false for a type openDNA does not read, with `why`).

A hub with several genomes needs `--genome NAME` (one of `hub.genomes`): ask the user which one is their sequence when it is not clear.

## 2. Copy the files

```
hub get HUB --out NEW_FOLDER
```

`--out` must be a folder that does not exist or is empty: nothing is overwritten. Every track openDNA reads is copied (hidden ones too: hidden is only how the hub shows them at first); `files` lists each with its group and whether it starts hidden, `failed` the ones that could not be copied. A remote hub's files can be hundreds of megabytes: say how many tracks there are before `get`, and tell the user it downloads from the hub's server.

## What to tell the user

- How many tracks, in which groups, and how many start hidden. The tracks of a type openDNA does not read are in `warnings`: name them.
- Where the folder is. To see them in openDNA: open their sequence, then Import > Features and Tracks > Track Hub (it reads the hub, keeps its groups, labels and colours), or load the copied files with Signal Track / Annotation Track.
- Any file in `failed`, with its error.

## Example

```
hub list demo/hub
hub get demo/hub --out example_hub
```

## Errors you may see

| Message | Fix |
|---|---|
| `... is not a track hub` | the file's first lines must be `hub NAME` and `genomesFile FILE`: check the address |
| `the hub has several genomes ...` | add `--genome` with one of the names listed |
| `cannot read ...` | the address or a file of the hub cannot be reached: check it, or the network |
| `The output folder ... is not empty` | choose a new folder with `--out` |
