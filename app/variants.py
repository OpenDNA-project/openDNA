"""Variants from reads: where the reads of a loaded BAM disagree with the sequence (SNPs,
deletions, insertions), from the background scan in alignments.py. They are drawn as a
track of the Alignment view, under each BAM's coverage; right-click one there (or the
track) to write it into the sequence or to annotate it as a `variation` feature.

Applying writes the changes right to left, so earlier positions stay valid (e.g. to
correct a plasmid map after whole-plasmid sequencing). Each change is its own undo step.
Substitutions are then checked again; after an indel the reads no longer line up with the
new coordinates, so the reads should be mapped again.
"""


def change_text(v):
    """"A > G", "deletion of 5 bp (ACGTA)", "insertion of TT"."""
    if v.kind == "snp":
        return f"{v.ref} > {v.alt}"
    if v.kind == "del":
        return f"deletion of {len(v.ref)} bp" + (f" ({v.ref})" if len(v.ref) <= 10 else "")
    return f"insertion of {v.alt or '?'}"


def can_apply(v):
    """An insertion whose bases the reads did not settle (N) cannot be written in."""
    return v.kind != "ins" or bool(v.alt) and "N" not in v.alt


def apply_variants(editor, entry, chosen):
    """Write `chosen` variants of BAM `entry` into the sequence, last position first.
    Returns how many were applied (0 while the DNA is locked)."""
    chosen = sorted((v for v in chosen if can_apply(v)), key=lambda v: v.pos, reverse=True)
    if not chosen or not editor.dna_editable():
        return 0
    doc = editor.doc
    for v in chosen:
        if v.kind == "snp":
            doc.overwrite(v.pos, v.alt)
        elif v.kind == "del":
            doc.delete(v.pos, v.pos + len(v.ref))
        else:
            doc.insert(v.pos, v.alt)
    shifted = any(v.kind != "snp" for v in chosen)
    # Substitutions keep coordinates: scan the corrected sequence again. Indels move
    # everything after them, so the reads no longer line up: no flags until re-mapped.
    entry.variants, entry.cache = ([] if shifted else None), None
    editor._changed()
    if shifted:
        editor.error.emit(f"Applied {len(chosen)} change(s). Indels moved the coordinates: map the reads again "
                          "(Tools → Map Reads) to check the corrected sequence.")
    return len(chosen)


def variant_feature(v, length, source=""):
    """(kind, location, qualifiers) of a `variation` feature for variant `v` (GenBank's
    /replace holds the sample's bases; an insertion covers the two bases around it)."""
    if v.kind == "snp":
        location = str(v.pos + 1)
    elif v.kind == "del":
        location = f"{v.pos + 1}..{v.pos + len(v.ref)}"
    else:
        location = f"{max(1, v.pos)}..{min(length, v.pos + 1)}"
    note = f"{v.count:,} of {v.depth:,} reads ({100 * v.count / max(1, v.depth):.0f} %)" + (f" in {source}" if source else "")
    return ("variation", location, [("label", change_text(v)), ("replace", "" if v.kind == "del" else v.alt), ("note", note)])


def annotate_variants(editor, entry, chosen):
    """Add a `variation` feature for each of `chosen` (one undo step). Returns how many."""
    chosen = list(chosen)
    if chosen:
        editor.add_features([variant_feature(v, editor.length, entry.name) for v in chosen])
    return len(chosen)
