"""Does an imported track (GFF3 / GTF / BED / VCF) or BAM file belong to the open sequence? `mismatches` lists what does not fit:
another sequence name, items past the end of the sequence, VCF REF bases that are not the bases here, a BAM reference of another
length. The import still goes ahead (the user may know better, e.g. a renamed chromosome): this only warns."""
from pathlib import Path

CHECKED = 5000   # items read for the checks (the first ones of the sequence): enough to tell, fast on any file


def names_match(a, b):
    norm = lambda s: s.lower().removeprefix("chr")
    return bool(a) and bool(b) and (a == b or norm(a) == norm(b))


def match_contig(seqs, *names):
    """The sequence of a track file (`seqs`) that is the open one, or None: the same name as one of `names` (the document's
    name, the FASTA record), else the same name once "chr" and the case are ignored ("chr1" is "1") if only one fits."""
    for n in names:
        if n and n in seqs:
            return n
    plain = lambda x: x.lower().removeprefix("chr")
    wanted = {plain(n) for n in names if n}
    close = [q for q in seqs if plain(q) in wanted]
    return close[0] if len(close) == 1 else None


def mismatches(editor, track, seqid, length=None, path=""):
    """Plain sentences, empty when the track fits. `track`: a opendna_core Track (None for a BAM, then `length` is the length of
    its reference `seqid`)."""
    out = []
    here = editor.fasta_record or editor.doc.name
    name = Path(path).name or "The file"
    if not names_match(seqid, here) and not names_match(seqid, editor.doc.name):
        out.append(f"{name} is for the sequence '{seqid}', and this one is '{here}'.")
    if track is None:
        if length is not None and length != editor.length:
            out.append(f"Its reference '{seqid}' is {length:,} bp; this sequence is {editor.length:,} bp: the reads were mapped to "
                       "another version of it, so they may not line up.")
        return out
    past, more = track.query(seqid, editor.length, 1 << 62, CHECKED)
    if past:
        far = max(e for _, e, *_ in past)
        out.append(f"{len(past):,}{'+' if more else ''} item(s) lie beyond the end of this sequence ({editor.length:,} bp; "
                   f"the track reaches {far:,}).")
    items, _ = track.query(seqid, 0, editor.length, CHECKED)
    wrong, checked, example = 0, 0, None
    for start, end, _, label, kind in items:
        if kind != "variant" or end > editor.length:
            continue
        ref = label.split(" ")[-1].split(">")[0].upper()
        if not ref or "…" in ref or set(ref) - set("ACGTN"):   # shortened or symbolic allele: not comparable
            continue
        checked += 1
        here_bases = editor.doc.sequence(start, end).upper()
        if here_bases != ref:
            wrong += 1
            example = example or f"position {start + 1:,}: the file says {ref}, the sequence has {here_bases}"
    if wrong:
        out.append(f"{wrong:,} of {checked:,} variant(s) checked have REF bases that are not the bases of this sequence "
                   f"({example}): the variants were called against another version of it.")
    return out
