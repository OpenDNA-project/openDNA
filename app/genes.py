"""Gene groups ("master annotations"): one `gene` feature per gene, owning its parts
(promoter, RBS, CDS, terminator…).

Stored the standard GenBank way, so other software still reads the file: the gene and its
parts share a `/locus_tag` (NCBI files do this for gene + CDS). Files using only `/gene`
are grouped by that. Nothing else is stored: groups are recomputed from the features.

Every operation here is one undo step (Document.edit_features).
"""
import re
from types import SimpleNamespace

GENE = "gene"
# Order of parts inside a gene, upstream to downstream (others go last).
PART_ORDER = {"promoter": 0, "5'UTR": 1, "RBS": 2, "CDS": 3, "mRNA": 3, "3'UTR": 4, "terminator": 5}


def qualifier(f, key):
    return next((v for k, v in f.qualifiers if k == key), None)


def group_key(f):
    return qualifier(f, "locus_tag") or qualifier(f, "gene")


def compute_groups(features):
    """{key: SimpleNamespace(key, master, parts)} for every `gene` feature with a key, and
    {feature index: key} for the parts. A key shared by several gene features (rare,
    e.g. a duplicated region) groups under the first."""
    groups, part_of = {}, {}
    for f in features:
        key = group_key(f)
        if f.kind == GENE and key and key not in groups:
            groups[key] = SimpleNamespace(key=key, master=f, parts=[])
    for f in features:
        key = group_key(f)
        g = groups.get(key)
        if g and f is not g.master:
            g.parts.append(f)
            part_of[f.index] = key
    for g in groups.values():
        g.parts.sort(key=lambda p: (PART_ORDER.get(p.kind, 9), p.spans[0][0]))
    return groups, part_of


def with_qualifiers(f, **values):
    """f as a (kind, location, qualifiers) tuple with some qualifiers set (None = removed)."""
    quals = [(k, v) for k, v in f.qualifiers if k not in values]
    quals += [(k, v) for k, v in values.items() if v is not None]
    return (f.kind, f.location, quals)


def unique_key(taken, name):
    base = re.sub(r"[^\w.-]+", "_", name.strip()) or "gene"
    key, n = base, 2
    while key in taken:
        key, n = f"{base}_{n}", n + 1
    return key


def span_of(features):
    starts = [s for f in features for s, _ in f.spans]
    ends = [e for f in features for _, e in f.spans]
    return min(starts), max(ends)


def site_location(start, end, strand, length):
    from editor import site_location as location  # here: editor imports this module
    return location(start, end, strand, length)


def make_gene(editor, indices, name):
    """Group the features `indices` under a new gene feature named `name` spanning them.
    Returns the group key. Parts already in another gene move to this one."""
    parts = [editor.by_index[i] for i in indices if i in editor.by_index and editor.by_index[i].kind != GENE]
    if not parts:
        raise ValueError("Select at least one feature (not a gene) to group.")
    key = unique_key(set(editor.groups) | {group_key(f) for f in editor.features if group_key(f)}, name)
    lo, hi = span_of(parts)
    cds = [p for p in parts if p.kind == "CDS"]
    strand = (cds or parts)[0].strand or 1
    master = (GENE, site_location(lo, hi, strand, editor.length), [("gene", name), ("locus_tag", key), ("label", name)])
    changes = [(p.index, with_qualifiers(p, locus_tag=key, gene=name)) for p in parts]
    editor.edit_features(changes + [(None, master)])
    return key


def ungroup(editor, key):
    """Remove the gene feature; its parts stay, without the shared tag."""
    g = editor.groups[key]
    changes = [(p.index, with_qualifiers(p, locus_tag=None)) for p in g.parts]
    editor.edit_features(changes + [(g.master.index, None)])


def delete_gene(editor, key):
    """Delete the gene and all its parts."""
    g = editor.groups[key]
    editor.edit_features([(f.index, None) for f in [g.master] + g.parts])


def fitted_master(editor, key, replaced=None):
    """Change that grows the gene to cover all its parts (None if it already does).
    `replaced`: (index, start, end) of a part about to change, using its new span."""
    g = editor.groups[key]
    spans = [(replaced[1], replaced[2]) if replaced and p.index == replaced[0] else span_of([p]) for p in g.parts]
    if not spans:
        return None
    lo, hi = min(s for s, _ in spans), max(e for _, e in spans)
    m_lo, m_hi = span_of([g.master])
    if m_lo <= lo and hi <= m_hi:
        return None
    lo, hi = min(lo, m_lo), max(hi, m_hi)
    return (g.master.index, (GENE, site_location(lo, hi, g.master.strand or 1, editor.length), g.master.qualifiers))


def retag(features, taken):
    """Pasted features: locus tags already used in the target get fresh ones (so a pasted
    gene is a new group, not merged with the original)."""
    mapping, out = {}, []
    for kind, location, quals in features:
        tag = next((v for k, v in quals if k == "locus_tag"), None)
        if tag and tag in taken:
            if tag not in mapping:
                mapping[tag] = unique_key(taken | set(mapping.values()), tag)
            quals = [(k, mapping[tag] if k == "locus_tag" else v) for k, v in quals]
        out.append((kind, location, quals))
    return out
