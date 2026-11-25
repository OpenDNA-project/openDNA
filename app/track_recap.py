"""Recap: what every loaded track holds over a region, to compare many tracks at one place. The window (Recap Tracks Here,
Recap Highlighted Regions: recap_dialog) and the `recap` command call `recap` and `rows`, so their numbers agree.

Per region: its GC %, the sequence's own features there, and for each loaded track (in the order of the Alignment view):
- an annotation track: the items that overlap it (name, positions, strand, the share of the region they cover, the BED
  columns past itemRgb), on a circle the items across the origin too (their end past the length, as in the file's view);
- a signal: the bases with a value and the mean, lowest, highest and sum of the values, read from the values themselves (never
  a zoom level); with `per_base`, on a region of at most PER_BASE bases, the value at every base;
- a BAM file: the mean and highest read depth, and the variants found there when the reads were scanned (the window scans
  them; the command line does not).
Positions given out are 1-based and inclusive, as everywhere a user reads them.
"""
import csv
import io
import itertools
import math

from editor import HIDDEN_KINDS

PER_BASE = 1000      # bp: the longest region whose values are given base by base
DEPTH_BINS = 2000    # the highest depth is the highest of this many bins (exact on a region of at most this many bases)
COLUMNS = ("region", "start", "end", "track", "group", "kind", "items", "mean", "min", "max", "sum", "details")


def _overlap(s, e, start, end):
    return max(0, min(e, end) - max(s, start))


def _number(v):
    """A value for JSON: None for NaN (no value), floats rounded to 6 significant figures."""
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return None
    return float(f"{v:.6g}") if isinstance(v, float) else v


def _item(item, start, end):
    s, e, strand, name, kind, _, note, extra = item
    out = {"name": name, "type": kind, "start": s + 1, "end": e, "strand": "+" if strand > 0 else "-" if strand < 0 else ".",
           "covers": round(_overlap(s, e, start, end) / (end - start), 4)}
    if note:
        out["note"] = note
    if extra:
        out["columns"] = extra.split("\t")
    return out


def _track(entry, start, end, per_base, length):
    if entry.kind == "track":
        items, more = entry.track.query(entry.seqid, start, end, 100_000)
        found = [_item(i, start, end) for i in items]
        if getattr(entry, "wraps", False):   # items across the origin end past the length: their part after it is one turn on
            seen = {i[:5] for i in items}
            found += [_item(i, start + length, end + length) for i in entry.track.query(entry.seqid, start + length, end + length, 100_000)[0]
                      if i[:5] not in seen]
        return {"items": found, **({"more": True} if more else {})}
    if entry.kind == "signal":
        covered, total, low, high = entry.signal.stats(entry.seqid, start, end)
        out = {"bases_with_value": int(covered), "mean": _number(total / covered if covered else None), "min": _number(low),
               "max": _number(high), "sum": _number(total)}
        if per_base and end - start <= PER_BASE:
            out["values"] = [_number(mean) for mean, _ in entry.signal.bins(entry.seqid, start, end, end - start)]
        return out
    depths = entry.bam.coverage(entry.seqid, start, end, min(DEPTH_BINS, end - start))
    out = {"mean_depth": _number(sum(depths) / len(depths) if depths else 0.0), "max_depth": _number(max(depths, default=0.0))}
    if entry.variants is not None:
        out["variants"] = [{"position": v.pos + 1, "change": v.text, "level": v.level} for v in entry.variants
                           if v.pos < end and v.pos + max(1, len(v.ref)) > start]
    return out


def recap(editor, regions, per_base=False):
    """[{name, start, end, length, gc, features, tracks}] for `regions` [(name, start, end)] (0-based, end excluded, within the
    sequence). `tracks`: [{track, label, group, kind, ...}] in the view's order (see the module's docstring for what each kind
    gives). Reads files: call it off the window's thread."""
    layout = editor.track_layout
    group_of = {k: g.name for g in layout.groups() for k in g.keys}
    order = {k: i for i, k in enumerate(layout.keys())}
    entries = sorted(editor.loaded(), key=lambda e: order.get(e.name, len(order)))
    out = []
    for name, start, end in regions:
        gc = editor.doc.gc_bins(start, end, 1)[0]
        features = [{"label": f.label, "type": f.kind, "start": f.spans[0][0] + 1, "end": f.spans[-1][1],
                     "strand": "+" if f.strand > 0 else "-",
                     "covers": round(sum(_overlap(s, e, start, end) for s, e in f.spans) / (end - start), 4)}
                    for f in editor.doc.features_in(start, end) if f.kind not in HIDDEN_KINDS]
        tracks = [{"track": e.name, "label": layout.label(e.name), "group": group_of.get(e.name), "kind": e.kind,
                   **_track(e, start, end, per_base, editor.length)} for e in entries]
        out.append({"name": name, "start": start + 1, "end": end, "length": end - start, "gc": _number(gc * 100 if gc == gc else None),
                    "features": features, "tracks": tracks})
    return out


def _share(c):
    return "<1%" if 0 < c < 0.01 else f"{c:.0%}"


def _details(t):
    if t["kind"] == "track":
        return "; ".join(f"{i['name'] or '(no name)'} {i['start']:,}-{i['end']:,} {i['strand']} {_share(i['covers'])}"
                         + (f" [{' '.join(i['columns'])}]" if i.get("columns") else "") for i in t["items"]) + (" …" if t.get("more") else "")
    if t["kind"] == "signal":
        return ", ".join("" if v is None else f"{v:g}" for v in t["values"]) if "values" in t else ""
    return "; ".join(f"{v['position']:,} {v['change']}" for v in t.get("variants", [])) if "variants" in t else "variants not scanned"


def _mean(rows, column):
    values = [r[column] for r in rows if r[column] is not None]
    return _number(sum(values) / len(values)) if values else None


def rows(recaps):
    """The recap as a table: one row per region and track, the sequence's features first, and after each group's tracks a row
    "<group> (mean of N)" with the mean of their numbers (a group of signals: the mean of their means, lowest, highest, sums). A
    group of tracks of several kinds averages its signals only (an item count and a number of bases do not add up). Each row is
    a dict of COLUMNS."""
    out = []
    for r in recaps:
        base = {"region": r["name"], "start": r["start"], "end": r["end"]}
        features = "; ".join(f"{f['label']} ({f['type']})" for f in r["features"])
        out.append({**base, "track": "Features (this sequence)", "group": None, "kind": "features", "items": len(r["features"]),
                    "mean": None, "min": None, "max": None, "sum": None,
                    "details": features if r["gc"] is None else f"GC {r['gc']:.1f} %; {features}"})
        for group, run in itertools.groupby(r["tracks"], key=lambda t: t["group"]):   # (a group's tracks are next to each other)
            members = [{**base, "track": t["label"], "group": t["group"], "kind": t["kind"],
                        "items": len(t["items"]) if t["kind"] == "track" else t.get("bases_with_value"),
                        "mean": t.get("mean", t.get("mean_depth")), "min": t.get("min"), "max": t.get("max", t.get("max_depth")),
                        "sum": t.get("sum"), "details": _details(t)} for t in run]
            out += members
            alike = members if len({m["kind"] for m in members}) == 1 else [m for m in members if m["kind"] == "signal"]
            if group is not None and len(alike) > 1:
                out.append({**base, "track": f"{group} (mean of {len(alike)}" + ("" if alike is members else " signals") + ")",
                            "group": group, "kind": "group", **{c: _mean(alike, c) for c in ("items", "mean", "min", "max", "sum")},
                            "details": ""})
    return out


def tsv(table, sep="\t"):
    """The rows as text with a header line (tab-separated, or `sep`: "," for a CSV file)."""
    buffer = io.StringIO()
    writer = csv.writer(buffer, delimiter=sep, lineterminator="\n")
    writer.writerow(COLUMNS)
    writer.writerows([["" if row[c] is None else row[c] for c in COLUMNS] for row in table])
    return buffer.getvalue()
