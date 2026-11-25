"""Tracks into the genome: a track's items as features of the sequence (Add as Feature, Copy to Features, the
`track-features` command), the stretches where a signal is high (Regions Above Threshold, the `signal-regions` command), a
track's items as BED, and the next or previous item from a position. No window here: the menus and the commands call these.
"""
from editor import site_location, split_at_origin

NAME_FIELDS = ("label", "gene", "locus_tag", "product", "note")   # the qualifier an item's name goes in


def default_type(label):
    """The feature type a track's items become unless chosen: CDS for a track named like "RefSeq_2022_CDS", else misc_feature."""
    stem = label.rsplit(".", 1)[0].lower()
    return "CDS" if stem.endswith(("_cds", "-cds", " cds")) or stem == "cds" else "misc_feature"


def items_in(entry, spans):
    """The items of annotation track `entry` that overlap any of `spans` [(start, end)], each once, in order of start."""
    seen = {}
    for s, e in spans:
        for item in entry.track.query(entry.seqid, s, e, 1 << 40)[0]:
            seen.setdefault(item[:5], item)
    return sorted(seen.values(), key=lambda i: (i[0], i[1]))


def item_feature(item, track_label, kind, name_key, length):
    """(type, location, qualifiers) of one item: its name in `name_key`, a note "from track <label>", its own note, and each
    BED column past itemRgb as a note of its own."""
    s, e, strand, name, _, _, note, extra = item
    qualifiers = ([(name_key, name)] if name else []) + [("note", f"from track {track_label}")]
    qualifiers += [("note", text) for text in ([note] if note else []) + (extra.split("\t") if extra else []) if text]
    return kind, site_location(s, e, strand, length), qualifiers


def feature_keys(features):
    """What tells two features apart for `to_features`: {(type, spans, strand)} of `features`."""
    return {(f.kind, tuple(tuple(s) for s in f.spans), f.strand) for f in features}


def to_features(items, track_label, kind, name_key, have, length):
    """(features to add [(type, location, qualifiers)], how many were skipped): an item matching a feature already there, or
    one added before it (same type, span and strand), is skipped. `have`: feature_keys of the sequence's features, to which the
    new ones are added (one set for several tracks: the same item in two of them is added once)."""
    out, skipped = [], 0
    for item in items:
        key = (kind, tuple(split_at_origin(item[0], item[1], length)), -1 if item[2] < 0 else 1)
        if key in have:
            skipped += 1
            continue
        have.add(key)
        out.append(item_feature(item, track_label, kind, name_key, length))
    return out, skipped


def signal_regions(signal, seqid, start, end, threshold, min_length=1, gap=0):
    """[(start, end, mean, highest)]: the stretches of [start, end) where the signal is at least `threshold`, those closer than
    `gap` bases joined into one, those shorter than `min_length` left out (after joining). The mean and highest are over
    the whole stretch, the joined bases included (bases with no value are not counted). Read from the values themselves."""
    joined = []
    for s, e in signal.runs(seqid, start, end, threshold):
        if joined and s - joined[-1][1] <= gap:
            joined[-1][1] = e
        else:
            joined.append([s, e])
    out = []
    for s, e in joined:
        if e - s >= min_length:
            covered, total, _, high = signal.stats(seqid, s, e)
            out.append((s, e, total / covered if covered else float("nan"), high))
    return out


def regions_bed(found, seqid, name):
    """The regions of `signal_regions` as BED lines (sequence, start, end, "<name>_<n>", the mean as the score column)."""
    return "".join(f"{seqid}\t{s}\t{e}\t{name}_{i}\t{mean:.6g}\t.\n" for i, (s, e, mean, _) in enumerate(found, 1))


def items_bed(items, seqid):
    """Items as BED lines: sequence, start, end, name, score 0, strand, and their colour and own columns when they have them."""
    lines = []
    for s, e, strand, name, _, color, _, extra in items:
        line = [seqid, str(s), str(e), name or ".", "0", "+" if strand > 0 else "-" if strand < 0 else "."]
        if color or extra:
            line += [str(s), str(e), color if color and color.count(",") == 2 else "0"]
        lines.append("\t".join(line + (extra.split("\t") if extra else [])) + "\n")
    return "".join(lines)


def next_item(entry, pos, forward=True):
    """The first item of annotation track `entry` that starts after base `pos` (`forward`), or the last one that starts before
    it; None when there is none. Searches a window that doubles until it finds one (a genome's items stay quick)."""
    width = 10_000
    while True:
        if forward:
            items, _ = entry.track.query(entry.seqid, pos + 1, pos + 1 + width, 1 << 20)
            after = [i for i in items if i[0] > pos]
            if after:
                return min(after, key=lambda i: (i[0], i[1]))
            if not entry.track.query(entry.seqid, pos + 1 + width, 1 << 62, 1)[0]:
                return None
        else:
            items, _ = entry.track.query(entry.seqid, max(0, pos - width), pos, 1 << 20)
            before = [i for i in items if i[0] < pos]
            if before:
                return max(before, key=lambda i: (i[0], i[1]))
            if pos - width <= 0:
                return None
        width *= 2
