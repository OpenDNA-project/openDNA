"""The commands of the command line for agents that are not table batches (see cli.py, which registers them and holds what they share).
Each one calls the function the matching window calls and answers with one JSON: the settings that shaped it (`settings`) and `warnings` to relay.
Read-only analyses (no `check`/`run` pair, nothing written): `sites`, `guides`, `repeats`, `search`, `align`, `msa` (`--out` writes one FASTA),
`checkup`, `annotate` (`--add` writes a copy), `variants` (`--apply` writes a copy). Files: `export`, `order`. Network: `ncbi`.
Two are batches with `check` then `run`, a table in and files out, like the ones in cli.py: `codon` and `gibson`."""
import json
import re
import sys
from contextlib import contextmanager
from pathlib import Path
from types import SimpleNamespace

C = None   # the cli module (set by register): its helpers (load, given, ...) are shared, not copied
UserError = None


def csv_list(text):
    return [t.strip() for t in (text or "").split(",") if t.strip()]


def answer(settings, **fields):
    return {"settings": settings, "warnings": C.setting_warnings(settings), **fields}


def pick_enzymes(names):
    """The enzymes called `names` (each once, in order); an unknown name is an error that lists the known ones."""
    import enzymes
    by_name = {e.name.lower(): e for e in enzymes.COMMON}
    picked = {}
    for name in names:
        if name.lower() not in by_name:
            raise UserError(f"No enzyme called {name!r}. The enzymes openDNA knows: {', '.join(e.name for e in enzymes.COMMON)}.")
        picked[by_name[name.lower()].name] = by_name[name.lower()]
    return list(picked.values())


def end1(end, n):
    """The 1-based last base of a 0-based end (exclusive) on a circle of n, so a piece that runs over the origin ends where it should."""
    return (end - 1) % n + 1


def region_bounds(doc, text):
    """(start, end) 0-based, end exclusive, of `START..END` (1-based, inclusive) or of the whole sequence when `text` is empty."""
    if not text:
        return 0, doc.length
    first, last = C.target_range(doc, text)
    return first - 1, last


# ---- sites: restriction enzymes and digests ---------------------------------------------------------------------------------------

def cmd_sites(args):
    import enzymes
    import tools
    doc = C.load(args.file)
    chosen = pick_enzymes(csv_list(args.enzymes) + csv_list(args.digest))
    digest = [e.name for e in pick_enzymes(csv_list(args.digest))]
    lo = 1 if args.min_sites is None else args.min_sites
    flt = enzymes.EnzymeFilter(min_sites=lo, max_sites=1 if args.min_sites is None and args.max_sites is None else args.max_sites)   # (the window's default: unique cutters)
    seq = doc.sequence(0, None).upper()
    cuts = enzymes.find_cuts(seq, doc.circular, chosen or None)
    listed = []
    for e in chosen or enzymes.COMMON:
        found = cuts.get(e.name, [])
        if chosen or enzymes.wanted(flt, e, len(found)):
            listed.append({"enzyme": e.name, "site": e.site, "notation": enzymes.notation(e), "ends": enzymes.describe(e).split(", ", 1)[1], "count": len(found),
                           "cuts": [{"top_after_base": t, "bottom_after_base": b} for t, b in found],   # (0 = before the first base)
                           "sites": [{"start": s + 1, "end": end1(en, doc.length), "strand": st} for s, en, st in enzymes.site_spans(seq, doc.circular, e)]})
    out = answer({"filter": {"min_sites": flt.min_sites, "max_sites": flt.max_sites, "min_site_length": flt.min_size,
                             "applies": "no: the enzymes were named" if chosen else "yes"}}, length=doc.length, circular=doc.circular, enzymes=listed)
    if digest:
        positions = {}
        for name in digest:
            for t, _ in cuts.get(name, []):
                positions.setdefault(t, []).append(name)
        out["digest"] = {"enzymes": digest, "fragments": [
            {"size": size, "start": a + 1, "end": b, "left": left, "right": right}
            for size, a, b, left, right in sorted(tools.fragments(doc.length, doc.circular, positions), reverse=True)]}
    return out


# ---- guides: Cas9 / Cas12a ---------------------------------------------------------------------------------------------------------

def cmd_guides(args):
    import cas9
    from cut_sites import NUCLEASES
    doc = C.load(args.file)
    by_name = {n.name.lower(): n for n in NUCLEASES.values()}
    if args.nuclease.lower() not in by_name:
        raise UserError(f"No nuclease called {args.nuclease!r}. Choices: {', '.join(n.name for n in by_name.values())}.")
    nuc = by_name[args.nuclease.lower()]
    a, b = region_bounds(doc, args.region)
    offset, seq, wrap, error = cas9.stretch(doc, (a, b) if args.region else None)
    if error:
        raise UserError(error)
    found = cas9.candidates(seq, doc.circular, nuc, offset, wrap)
    if args.top is not None and args.top < 1:
        raise UserError("--top must be at least 1.")
    shown = found[:args.top or cas9.MAX_CANDIDATES]
    max_mm = C.given(args.max_mm, cas9.DEFAULT_MAX_MM, cas9.MAX_MM_LIMITS, "--max-mm")
    ref, specs = None, []
    if shown:   # (the background genome is only opened when there is a guide to check)
        import primer_panel
        ref = None if args.no_genome else primer_panel.genome()
        specs = cas9.scan(doc.sequence(0, None), doc.circular, doc.name, shown, nuc, max_mm["value"], ref)
    n = doc.length
    guides = []
    for rank, (g, spec) in enumerate(zip(shown, specs), 1):
        guides.append({"rank": rank, **cas9.guide_row(nuc, g, spec, n)})
    genome_name = ref.name if ref is not None else None
    settings = {"cas9": {"nuclease": {"value": nuc.name, "source": "command line"}, "max_mismatches": max_mm,
                         "region": {"value": args.region or "whole sequence", "source": "command line" if args.region else "default"}},
                "background_genome": {"value": genome_name, "source": "none: --no-genome" if args.no_genome else "openDNA settings (Primers > Options)" if ref else "none set: no genome check"}}
    note = f"the best {len(shown)} of {len(found)} are listed (--top)" if len(found) > len(shown) else ""
    return answer(settings, guides_found=len(found), note=note, guides=guides)


# ---- repeats: what two sequences share (the Compare view of Local Align) -----------------------------------------------------------

def cmd_repeats(args):
    import repeats
    a = C.load(args.file)
    b = C.load(args.other) if args.other else a
    lim = repeats.LIMITS
    settings = {"repeats": {"min_length": C.given(args.min_len, repeats.DEFAULT_MIN_LEN, lim["min_len"], "--min-len"),
                            "mismatches": C.given(args.mismatches, repeats.DEFAULT_MISMATCHES, lim["mismatches"], "--mismatches")}}
    sa = a.sequence(0, None)
    res = repeats.pair_repeats(sa, sa if b is a else b.sequence(0, None), settings["repeats"]["min_length"]["value"], settings["repeats"]["mismatches"]["value"])
    return answer(settings, compared={"a": a.name, "b": b.name, "same_file": b is a},
                  repeats=[{"kind": k, "start_in_a": i + 1, "end_in_a": i + n, "start_in_b": j + 1, "end_in_b": j + n, "length": n, "differences": d}
                           for k, i, j, n, d in res.repeats], note=res.capped)


# ---- search: the search bar's engine -----------------------------------------------------------------------------------------------

def cmd_search(args):
    import search
    doc = C.load(args.file)
    opts = search.Options(both=not args.single_strand, fuzzy=args.fuzzy or 0, translations=args.protein, details=args.details, enzymes=args.enzymes, primers=False)
    ed = doc_editor(doc)
    hits, note = search.search_sequence(ed, args.query, opts)
    n = doc.length
    return answer({"search": {"both_strands": {"value": opts.both}, "fuzzy_mismatches": {"value": opts.fuzzy}, "protein_motif_in_translations": {"value": opts.translations},
                              "feature_details": {"value": opts.details}, "enzyme_names": {"value": opts.enzymes}}},
                  query=args.query, note=note, hits=[{"kind": h.kind, "label": h.label, "start": h.start + 1, "end": end1(h.end, n), "strand": h.strand, "mismatches": h.mismatches,
                                                      **({"feature_index": h.feature} if h.feature is not None else {})} for h in hits])


# ---- align: where a query matches the plasmid (the Local Align window) ----------------------------------------------------------------

def query_bases(text):
    """The DNA of a QUERY argument: a file (GenBank or FASTA: its first record) when one is named, else the bases written there."""
    import local_align
    if Path(text).is_file():
        return C.load(text).sequence(0, None).upper()
    try:
        return local_align.clean_query(text)
    except local_align.QueryError as e:
        raise UserError(str(e))


def cmd_align(args):
    import local_align
    doc = C.load(args.file)
    query = query_bases(args.query)
    hits = local_align.find_hits(doc.sequence(0, None), query, doc.circular)
    n = doc.length
    rows = [{"piece": i + 1, "query_start": h.qstart + 1, "query_end": h.qend, "start": h.start + 1, "end": end1(h.end, n), "strand": h.strand, "identity": round(h.identity, 4),
             "score": h.score, "mismatches": len(h.mismatches), "insertions": len(h.insertions), "deletions": len(h.deletions),
             "before": local_align.junction_text(hits[i - 1], h, n, doc.circular) if i else "start of the query"} for i, h in enumerate(hits)]
    note = f"at most {local_align.MAX_HITS} pieces are listed" if len(hits) >= local_align.MAX_HITS else ""
    return answer({"align": {"query_bp": {"value": len(query)}, "max_pieces": {"value": local_align.MAX_HITS}}}, summary=local_align.summarize(hits, n, doc.circular) or "No hit",
                  note=note, pieces=rows)


# ---- msa: several sequences aligned (the Multiple alignment tab) --------------------------------------------------------------------

def new_path(path):
    """`path` as a Path, refused when it exists: nothing is overwritten."""
    if Path(path).exists():
        raise UserError(f"{path} already exists: choose a new file name (nothing is overwritten).")
    return Path(path)


def cmd_msa(args):
    import msa
    docs = [C.load(f) for f in args.files]
    names = msa.unique_names([d.name for d in docs])
    out = new_path(args.out) if args.out else None
    try:
        aln = msa.align(names, [d.sequence(0, None) for d in docs])
    except ValueError as e:
        raise UserError(str(e))
    if out:
        out.write_text(msa.to_fasta(aln.names, aln.rows), encoding="utf-8")
    return answer({"msa": {"aligner": {"value": f"built-in progressive aligner (up to {msa.MAX_BUILTIN} bp each)"}, "identity_against": {"value": names[0]}}},
                  width=aln.width, consensus=aln.consensus(),
                  sequences=[{"name": nm, "length": d.length, "identity_to_first": round(aln.identity(i), 2), "aligned": row}
                             for i, (nm, d, row) in enumerate(zip(aln.names, docs, aln.rows))], **({"file": str(out)} if out else {}))


# ---- codon: recode a CDS for a codon usage table (the Recode window) --------------------------------------------------------------

def gc_range(args):
    import codon_tools as CT
    return None if args.no_gc else (CT.DEFAULT_GC[0] if args.gc_low is None else args.gc_low, CT.DEFAULT_GC[1] if args.gc_high is None else args.gc_high)


def codon_row(doc, row, args, usage, avoid, out_dir):
    import codon_tools as CT
    from assembly import snapshot
    from mutation import find_feature
    rec = {"id": row["id"], "feature": row["feature"], "status": "error"}
    try:
        f = find_feature(doc, row["feature"])
        if f.kind != "CDS":
            raise UserError(f"{row['feature']!r} is a {f.kind}, not a CDS.")
        s = CT.feature_source(doc, f)
        if isinstance(s, str):
            raise UserError(s)
        r = CT.recode_source(s, usage, method=args.method, avoid=avoid, max_run=CT.MAX_RUN if args.max_run is None else args.max_run,
                             rbs=not args.no_rbs, gc_range=gc_range(args), window=CT.GC_WINDOW if args.window is None else args.window, seed=args.seed)
    except (UserError, ValueError, RuntimeError) as e:
        rec["error"] = str(e)
        return rec
    rec.update(status="ok", start=s.start + 1, end=s.end, strand=s.strand, genetic_code=s.table, codons=r.codons, codons_changed=r.changed, gc_before=round(r.gc_before, 1),
               gc_after=round(r.gc_after, 1), warnings=[f"Could not be met: {u}" for u in r.unmet], sequence=r.seq)
    if out_dir is not None:
        made = snapshot(doc, doc.name)
        made.replace(s.start, s.end, CT.recoded_bases(s, r))
        made.write_genbank(str(out_dir / f"{row['file_stem']}.gb"), None)
        rec["file"] = f"{row['file_stem']}.gb"
    return rec


def cmd_codon(args):
    import codon_tools as CT
    import codon_usage
    doc, rows = C.load(args.plasmid), C.read_table(args.targets, "feature")
    C.assign_files(rows)
    tables = codon_usage.all_tables()
    if args.usage and Path(args.usage).is_file():   # a table in a CSV file (the window's Import reads the same file)
        try:
            usage = codon_usage.read_csv(args.usage)
        except ValueError as e:
            raise UserError(f"{args.usage}: {e}")
    else:
        usage = tables.get(args.usage) if args.usage else codon_usage.active()
    if args.usage and usage is None:
        raise UserError(f"No codon usage table called {args.usage!r}. Choices: {', '.join(tables) or 'none (import one in openDNA: Tools > Settings > Codon usage)'}.")
    if usage is None:
        raise UserError("No codon usage table is set: name one with --usage, or choose one in openDNA (Tools > Settings > Codon usage).")
    if args.method not in (CT.MATCHED, CT.HIGHEST):
        raise UserError(f"--method must be {CT.MATCHED!r} or {CT.HIGHEST!r}.")
    avoid = pick_enzymes(csv_list(args.avoid))
    gc = gc_range(args)
    if gc and gc[0] > gc[1]:
        raise UserError("--gc-low is above --gc-high.")
    out_dir = C.make_out_dir(args, args.plasmid, "_codon")
    lim, g = CT.LIMITS, CT.DEFAULT_GC
    settings = C.effective_settings({"codon": {
        "usage_table": {"value": usage.name, "source": "command line" if args.usage else "openDNA settings"}, "method": C.given(args.method if args.method != CT.MATCHED else None, CT.MATCHED),
        "avoid_sites": {"value": [e.name for e in avoid]},
        **({"gc_limits": {"value": "none", "source": "command line: --no-gc"}} if args.no_gc else
           {"gc_low": C.given(args.gc_low, g[0], lim["gc"], "--gc-low"), "gc_high": C.given(args.gc_high, g[1], lim["gc"], "--gc-high"),
            "gc_window": C.given(args.window, CT.GC_WINDOW, lim["window"], "--window")}), "max_run": C.given(args.max_run, CT.MAX_RUN, lim["max_run"], "--max-run"),
        "no_internal_rbs": {"value": not args.no_rbs}, "seed": {"value": args.seed, "source": "command line" if args.seed != 1 else "default (the window's first preview)"}}})
    records = C.run_rows(rows, lambda row: codon_row(doc, row, args, usage, avoid, out_dir), "feature")
    return C.finish(out_dir, settings, records, ["id", "feature", "file", "codons", "codons_changed", "gc_before", "gc_after", "warnings"],
                    lambda r: [[r["id"], r["feature"], r.get("file", ""), r["codons"], r["codons_changed"], r["gc_before"], r["gc_after"], " | ".join(r["warnings"])]],
                    lambda r: bool(r.get("warnings")))


# ---- export: the other file formats (Export menu) ---------------------------------------------------------------------------------

EXPORTS = ("genbank", "fasta", "gff3", "bed", "tbl", "proteins")


def cmd_export(args):
    import exports
    doc = C.load(args.file)
    out = new_path(args.out)
    ed = doc_editor(doc)
    if args.format == "genbank":
        doc.write_genbank(str(out), None)   # (also how a Benchling, APE or EMBL file becomes a plain GenBank file)
    elif args.format == "fasta":
        out.write_text(exports.sequence_fasta(ed), encoding="utf-8")
    else:
        out.write_text({"gff3": exports.gff3_text, "bed": exports.bed_text, "tbl": exports.tbl_text, "proteins": exports.proteins_fasta}[args.format](ed), encoding="utf-8")
    return answer({"export": {"format": {"value": args.format}}}, file=str(out), bytes=out.stat().st_size, features=len(ed.features))


# ---- order: a primer table as a supplier's order file (Primers > Export order) --------------------------------------------------------

def cmd_order(args):
    import primer_db
    rows = C.read_table(args.primers, "sequence")
    fmt = {"idt": primer_db.ORDER_FORMATS[0], "plain": primer_db.ORDER_FORMATS[1]}[args.format]
    problems = [f"line {r['_line']} ({r['id']}): " + why for r in rows if (why := r.get("_error") or primer_db.primer_problem(r["id"], r["sequence"].upper()))]
    if problems:
        raise UserError("nothing was written, the table has problems: " + "; ".join(problems))
    if args.scale not in primer_db.SCALES or args.purification not in primer_db.PURIFICATIONS:
        raise UserError(f"--scale is one of {', '.join(primer_db.SCALES)}; --purification one of {', '.join(primer_db.PURIFICATIONS)}.")
    out = new_path(args.out)
    written = primer_db.export_order([{"name": r["id"], "sequence": r["sequence"]} for r in rows], out, fmt, args.scale, args.purification)
    return answer({"order": {"format": {"value": fmt}, "scale": {"value": args.scale}, "purification": {"value": args.purification}}}, primers=len(rows), files=[str(w) for w in written])


# ---- variants: where sequencing reads (a BAM) differ from the plasmid (the Variants pane) -----------------------------------------

VARIANT_FLAGS = {"min_depth": "--min-depth", "variant_pct": "--variant-pct", "consensus_pct": "--consensus-pct", "min_baseq": "--min-baseq", "min_mapq": "--min-mapq"}


def cmd_variants(args):
    import alignments
    import track_check
    import variants as V
    from opendna_core import Bam
    doc = C.load(args.file)
    if not Path(args.reads).is_file():
        raise UserError(f"File not found: {args.reads}")
    try:
        bam = Bam(str(args.reads))
    except Exception as e:   # (the core raises OSError with its own words)
        raise UserError(f"Cannot read {args.reads} as a BAM file: {e}")
    lengths = dict(bam.sequences())
    if not lengths:
        raise UserError(f"{args.reads} has no sequence in its header.")
    seqid = args.seqid
    if seqid is None:
        seqid = track_check.match_contig(list(lengths), doc.name) or (next(iter(lengths)) if len(lengths) == 1 else None)   # (as the window picks it)
        if seqid is None:
            raise UserError(f"{args.reads} holds reads for several sequences ({', '.join(lengths)}): name the one that is this plasmid with --seqid.")
    elif seqid not in lengths:
        raise UserError(f"{args.reads} has no sequence called {seqid!r} (it has: {', '.join(lengths)}).")
    warnings = track_check.mismatches(SimpleNamespace(fasta_record=None, doc=doc, length=doc.length), None, seqid, lengths[seqid], args.reads)
    with tool_settings_block(alignments.TOOL, VARIANT_FLAGS, args) as settings:
        values = tuple(settings[k]["value"] for k in ("min_depth", "variant_pct", "consensus_pct", "min_baseq", "min_mapq"))
        try:
            found = alignments.scan_variants(bam, seqid, doc.sequence(0, None), values)
        except OSError as e:
            raise UserError(f"Cannot read the reads of {args.reads}: {e}")
    out = answer({"variants": settings}, reads=str(args.reads), sequence=seqid, reference_length=lengths[seqid], variants=[
        {"position": v.pos + 1, "kind": v.kind, "reference": v.ref, "alternative": v.alt, "reads": v.count, "depth": v.depth, "fraction": round(v.count / max(1, v.depth), 3),
         "level": v.level, "change": V.change_text(v), "can_apply": V.can_apply(v)} for v in found])
    out["warnings"] += warnings
    if args.apply:
        chosen = [v for v in found if v.level == "consensus" and V.can_apply(v)]
        target = new_path(args.apply)
        V.apply_to_document(doc, chosen)   # (the document was read for this command alone)
        doc.write_genbank(str(target), None)
        out["applied"] = {"file": str(target), "changes": len(chosen),
                          "note": "only the consensus changes (most reads disagree with the plasmid) that can be written in; after an insertion or deletion the coordinates move: map the reads again to check" if any(v.kind != "snp" for v in chosen) else "substitutions only: the coordinates did not move"}
    return out


# ---- annotate (Detect Features) and checkup (Plasmid Checkup): the tools that read their thresholds from Settings --------------------

def doc_editor(doc):
    """What the modules that take an editor (annotation, checkup, search, exports) read of one, for a document with no window: no selection,
    tracks or primer sites. A module that reads one more attribute needs it here (tests/test_agent_cli.py runs every command)."""
    from editor import HIDDEN_KINDS
    return SimpleNamespace(doc=doc, length=doc.length, circular=doc.circular, features=doc.features(), groups={}, fasta_record=None, tracks=[], primer_sites=[], start=0, end=0,
                           visible_features=lambda start, end, overlays=True: [f for f in doc.features_in(start, end) if f.kind not in HIDDEN_KINDS])


@contextmanager
def tool_settings_block(tool, flags, args):
    """{key: {value, source, min, max}} for the settings of `tool` that have a command-line flag (`flags`: key -> flag; the value is `args.<key>`). A value
    given on the command line is used for the length of the block only (tool_settings.override: it is never saved); otherwise the user's saved value."""
    import tool_settings
    saved = tool_settings.store()
    out = {}
    for p in tool_settings.TOOLS[tool]:
        if p.key not in flags:
            continue
        given = getattr(args, p.key)
        source = "saved in openDNA (Settings)" if saved.contains(f"{tool}/{p.key}") else "default"
        if p.minimum is None:   # a switch (True/False)
            out[p.key] = {"value": tool_settings.value(tool, p.key) if given is None else given, "source": "command line" if given is not None else source}
        else:
            out[p.key] = C.given(given, tool_settings.value(tool, p.key), (p.minimum, p.maximum), flags[p.key])
            if given is None:
                out[p.key]["source"] = source
    with tool_settings.override({(tool, k): v["value"] for k, v in out.items()}):
        yield out


ANNOTATE_FLAGS = {"min_dna": "--min-dna", "min_protein": "--min-protein", "min_length": "--min-length"}


def cmd_annotate(args):
    import annotation
    doc = C.load(args.file)
    ed = doc_editor(doc)
    with tool_settings_block(annotation.TOOL, ANNOTATE_FLAGS, args) as settings:
        hits = annotation.detect(ed, set(csv_list(args.kinds)) or None)
    n = doc.length
    out = answer({"detect_features": settings}, parts_in_library=len(annotation.load_parts()), hits=[
        {"name": h.name, "kind": h.kind, "start": h.start + 1, "end": end1(h.end, n), "strand": h.strand, "identity": round(h.identity, 4), "protein_match": bool(h.protein),
         "already_annotated": h.known, "source": h.source, "strength": h.measured} for h in hits])
    if args.add:
        target = new_path(args.add)
        fresh = [h for h in hits if not h.known]
        if fresh:
            doc.add_features(annotation.features_for(ed, fresh, False))
        doc.write_genbank(str(target), None)
        out["added"] = {"file": str(target), "features": len(fresh), "skipped_already_annotated": len(hits) - len(fresh)}
    return out


CHECKUP_FLAGS = {"min_promoter": "--min-promoter", "min_stem": "--min-stem", "rbs_reach": "--rbs-reach"}


def cmd_checkup(args):
    import checkup
    doc = C.load(args.file)
    args.predict, args.library = (False if args.no_predict else None), (False if args.no_library else None)   # (a switch left alone keeps the saved value)
    with tool_settings_block(checkup.TOOL, {**CHECKUP_FLAGS, "predict": "--no-predict", "library": "--no-library"}, args) as settings:
        res = checkup.run(doc_editor(doc))
    n = doc.length
    span = lambda el: {"start": el.start + 1, "end": end1(el.end, n), "strand": el.strand}
    cds = []
    for c in res.cds:
        f = c.feature
        cds.append({"cds": f.label, "start": f.spans[0][0] + 1, "end": f.spans[-1][1], "strand": f.strand, "verdict": c.verdict, "best_host": c.host, "reasons": c.reasons,
                    "hosts": {h: {"verdict": r.verdict, "promoter": r.unit.promoter.label if r.unit else None, "terminator": r.unit.terminator.label if r.unit and r.unit.terminator else None,
                                  "rbs": r.rbs.element.label if r.rbs else None, "reasons": r.reasons} for h, r in c.hosts.items()}})
    return answer({"checkup": settings}, note=res.note, hosts=list(res.hosts), cds=cds,
                  elements=[{"kind": e.kind, "label": e.label, "origin": e.origin, "detail": e.detail, **span(e)} for e in res.elements],
                  origins=[str(o) for o in res.origins], origin_notes=list(res.origin_notes))


# ---- gibson: Gibson / HiFi assembly of fragments (the Assembly panel's Gibson method) -------------------------------------------------

def cmd_gibson(args):
    import csv
    import assembly as A
    import recombineering as R
    from assembly import Frag, snapshot
    rows = C.read_table(args.fragments, "source", optional=("region", "orientation"))
    bad = [r["_error"] for r in rows if "_error" in r]
    if bad:   # (unlike the batch commands, an assembly needs every fragment: one bad row stops it)
        raise UserError("the fragments cannot be assembled with a row missing: " + " ".join(bad))
    names = [r["id"] for r in rows]
    if len(set(names)) != len(names):
        raise UserError(f"The ids must be different (they name the fragments and their primers): {', '.join(sorted({n for n in names if names.count(n) > 1}))}.")
    if len(rows) < 2:
        raise UserError("An assembly needs at least two fragments.")
    if not re.fullmatch(r"[A-Za-z0-9._-]+", args.name):
        raise UserError(f"--name {args.name!r} is not usable as a file name: letters, digits, . _ - only.")
    folder = Path(args.fragments).resolve().parent
    docs, frags = {}, []
    for r in rows:
        path = Path(r["source"])
        path = str(path if path.is_absolute() else folder / path)
        if path not in docs:
            docs[path] = C.load(path)
        src = docs[path]
        first, last = (1, src.length) if not r.get("region") else C.target_range(src, r["region"])
        orientation = (r.get("orientation") or "forward").lower()
        if orientation not in ("forward", "reverse"):
            raise UserError(f"Row {r['id']!r}: orientation must be forward or reverse, not {orientation!r}.")
        frag = Frag(r["id"], snapshot(src, src.name), first - 1, last, source=Path(path).name)
        if orientation == "reverse":
            frag.flip()
        frags.append(frag)
    lim, dflt = A.GIBSON_LIMITS, A.GIBSON_DEFAULTS
    settings = C.effective_settings({"gibson": {
        "circular": {"value": not args.linear, "source": "command line" if args.linear else "default"},
        "min_overlap": C.given(args.min_overlap, dflt["min_overlap"], lim["min_overlap"], "--min-overlap"), "max_overlap": C.given(args.max_overlap, dflt["max_overlap"], lim["max_overlap"], "--max-overlap"),
        "overlap_tm": C.given(args.overlap_tm, dflt["overlap_tm"], lim["overlap_tm"], "--overlap-tm"), "annealing_tm": C.given(args.tm, R.DEFAULT_TM, R.LIMITS["annealing_tm"], "--tm")}})
    g = {k: v["value"] for k, v in settings["gibson"].items()}
    if g["min_overlap"] > g["max_overlap"]:
        raise UserError("--min-overlap is above --max-overlap.")
    out_dir = C.make_out_dir(args, args.fragments, "_gibson")
    circular = g["circular"]
    A.auto_primers(frags, 0, circular, min_ov=g["min_overlap"], max_ov=g["max_overlap"], min_tm=g["overlap_tm"], tm=g["annealing_tm"])
    try:
        d = A.analyse(frags, 0, circular, min_ov=g["min_overlap"], min_tm=g["overlap_tm"], tm=g["annealing_tm"], max_ov=100)
    except ValueError as e:
        raise UserError(str(e))
    product = A.simulate(frags, d)
    import primer_panel
    primers = [{"fragment": f.name, "name": p[0], "sequence": p[1], "annealing_length": p[2], "annealing_tm": primer_panel.melting_temp(p[1][-p[2]:]), "role": role}
               for f in frags for role, p in (("forward", f.fwd), ("reverse", f.rev)) if p]
    levels = [lv for _, lv, _ in d.checks]
    verdict = "bad" if not d.ok else "warn" if "warn" in levels else "ok"   # (d.ok: every junction fits and no check is red, as the panel decides)
    result = {"settings": settings, "warnings": C.setting_warnings(settings), "verdict": verdict,
              "fragments": [{"id": f.name, "source": f.source, "start": f.start + 1, "end": f.end, "orientation": "reverse" if f.flipped else "forward", "length": f.core.length} for f in frags],
              "junctions": [{"text": t, "level": lv} for t, lv in zip(d.junctions, d.jlevels)], "checks": [{"check": c, "level": lv, "text": t} for c, lv, t in d.checks],
              "primers": primers, "product": {"length": product.length, "circular": product.circular}}
    if out_dir is not None:
        A.snapshot(product, args.name).write_genbank(str(out_dir / f"{args.name}.gb"), None)
        result["product"]["file"] = f"{args.name}.gb"
        with open(out_dir / "primers.csv", "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["fragment", "role", "name", "sequence", "annealing_length", "annealing_tm"])
            w.writerows([p["fragment"], p["role"], p["name"], p["sequence"], p["annealing_length"], C.tm_of(p)] for p in primers)
        (out_dir / "report.json").write_text(json.dumps(result, indent=1), encoding="utf-8")
        result["out"] = str(out_dir)
    for w in result["warnings"]:
        print(f"WARNING: {w}", file=sys.stderr)
    return result


# ---- ncbi: search and fetch records (Import from NCBI) -------------------------------------------------------------------------------

def cmd_ncbi(args):
    import ncbi
    import shutil
    db = args.db
    if args.action == "get" and not args.out:
        raise UserError("`ncbi get` needs --out NEW_FILE.")
    if args.limit is not None and args.limit < 1:
        raise UserError("--limit must be at least 1.")
    out = new_path(args.out) if args.action == "get" else None   # (before the network is used)
    key = ncbi.saved_key().strip()
    settings = {"ncbi": {"database": {"value": args.db, "source": "command line" if args.db != "nuccore" else "default"},
                         "api_key": {"value": "set" if key else "none", "source": "saved in openDNA (Import from NCBI)" if key else "default: 3 requests a second"},
                         "sends": {"value": f"only {args.query!r} to eutils.ncbi.nlm.nih.gov"}}}
    try:
        found = ncbi.search(db, args.query, key, ncbi.MAX_RESULTS if args.limit is None else args.limit)
        records = [{"accession": r.accession, "title": r.title, "length": r.length, "organism": r.organism, "updated": r.updated, "topology": r.topology, "source": r.source}
                   for r in found]
        if args.action == "search":
            return answer(settings, records=records)
        q = args.query.strip().lower()
        wanted = [r for r in found if r.accession.lower() == q or ("." not in q and r.accession.lower().split(".")[0] == q)]
        if not wanted:
            raise UserError(f"{args.query!r} is not the accession of a {args.db} record (`ncbi search` lists what matches).")
        record = wanted[0]
        if record.db == "nuccore" and record.length > ncbi.WARN_BASES and not args.allow_large:
            raise UserError(f"{record.accession} has {record.length:,} bases: a download of that size takes minutes; add --allow-large if the user wants it.")
        path, cached = ncbi.get(record, key)
        shutil.copyfile(path, out)
    except ncbi.NcbiError as e:
        raise UserError(str(e))
    return answer(settings, record=records[found.index(record)], file=str(out), bytes=out.stat().st_size, from_cache=cached)


def register(sub, cli):
    """Add the commands to the argument parser `sub`; `cli` is the cli module."""
    global C, UserError
    C, UserError = cli, cli.UserError
    s = sub.add_parser("sites", help="restriction sites (unique 6+ cutters by default) and the fragments of a digest")
    s.add_argument("file")
    s.add_argument("--enzymes", help="comma-separated names: only these (counted whatever their number of sites)")
    s.add_argument("--min-sites", type=int, help="list enzymes that cut at least this many times (default 1)")
    s.add_argument("--max-sites", type=int, help="list enzymes that cut at most this many times (default 1: unique cutters)")
    s.add_argument("--digest", help="comma-separated enzymes cut together: the fragments")
    s.set_defaults(run=cmd_sites)
    s = sub.add_parser("guides", help="Cas9/Cas12a guides ranked, with their off-target counts and a mark (Spacer Finder)")
    s.add_argument("file")
    s.add_argument("--nuclease", default="SpCas9", help="SpCas9 (default), SaCas9 or Cas12a")
    s.add_argument("--region", help="START..END (1-based) to look in (default: the whole sequence)")
    s.add_argument("--max-mm", type=int, help="count other sites up to this many mismatches (default 3)")
    s.add_argument("--top", type=int, help="how many of the best guides to report and check (default: 300)")
    s.add_argument("--no-genome", action="store_true", help="do not count sites in the background genome")
    s.set_defaults(run=cmd_guides)
    s = sub.add_parser("repeats", help="stretches shared by two files (or by one with itself): direct and inverted repeats")
    s.add_argument("file")
    s.add_argument("other", nargs="?", help="a second file (default: the first one again)")
    s.add_argument("--min-len", type=int, help="shortest shared stretch in bp (default 12)")
    s.add_argument("--mismatches", type=int, help="differences allowed inside a stretch (default 0)")
    s.set_defaults(run=cmd_repeats)
    s = sub.add_parser("search", help="a DNA sequence (with IUPAC codes), a protein motif or a name, in the sequence and its features")
    s.add_argument("file")
    s.add_argument("query")
    s.add_argument("--fuzzy", type=int, help="mismatches allowed (default 0)")
    s.add_argument("--single-strand", action="store_true", help="do not look for the reverse complement")
    s.add_argument("--protein", action="store_true", help="also look for a protein motif in the six-frame translations")
    s.add_argument("--details", action="store_true", help="names match any qualifier of a feature, not just its name")
    s.add_argument("--enzymes", action="store_true", help="an enzyme name finds its sites")
    s.set_defaults(run=cmd_search)
    s = sub.add_parser("align", help="where a query (bases, or a file) matches the plasmid: pieces, strands, mismatches, structural changes (Local Align)")
    s.add_argument("file")
    s.add_argument("query", help="DNA text, or a GenBank/FASTA file")
    s.set_defaults(run=cmd_align)
    s = sub.add_parser("msa", help="multiple alignment of two or more sequences (one per file), with identity and consensus")
    s.add_argument("files", nargs="+")
    s.add_argument("--out", help="also write the alignment as FASTA to this new file")
    s.set_defaults(run=cmd_msa)
    s = sub.add_parser("codon", help="recode CDS features for a codon usage table keeping the protein (check: report only; run: also write plasmids)")
    s.add_argument("action", choices=["check", "run"])
    s.add_argument("plasmid")
    s.add_argument("targets", help="CSV with the columns id,feature")
    s.add_argument("--out", help="run: a new output folder (default: PLASMID_codon next to the plasmid)")
    s.add_argument("--usage", help="a codon usage table: its name in openDNA, or a CSV file (default: the one chosen in openDNA)")
    s.add_argument("--method", default="Frequency matched", help='"Frequency matched" (default) or "Highest frequency"')
    s.add_argument("--avoid", help="comma-separated enzymes whose sites must not appear")
    s.add_argument("--gc-low", type=int, help="lowest GC percent in every window (default 35)")
    s.add_argument("--gc-high", type=int, help="highest GC percent (default 65)")
    s.add_argument("--no-gc", action="store_true", help="no GC limits")
    s.add_argument("--window", type=int, help="bases per GC window (default 50)")
    s.add_argument("--max-run", type=int, help="longest run of one base (default 5)")
    s.add_argument("--no-rbs", action="store_true", help="allow an internal ribosome binding site")
    s.add_argument("--seed", type=int, default=1, help="the random draw of a frequency-matched recoding (default 1; the same seed gives the same sequence)")
    s.set_defaults(run=cmd_codon, returns_code=True)
    s = sub.add_parser("export", help="write the sequence as genbank, fasta, gff3, bed, tbl (NCBI feature table) or proteins (FASTA of the CDS)")
    s.add_argument("file")
    s.add_argument("format", choices=EXPORTS)
    s.add_argument("out", help="a new file (never overwritten)")
    s.set_defaults(run=cmd_export)
    s = sub.add_parser("order", help="a primer table (columns id,sequence) as an IDT plate or plain name,sequence order file")
    s.add_argument("primers")
    s.add_argument("out", help="a new file (never overwritten)")
    s.add_argument("--format", choices=["idt", "plain"], default="idt")
    s.add_argument("--scale", default="25nm")
    s.add_argument("--purification", default="STD")
    s.set_defaults(run=cmd_order)
    s = sub.add_parser("variants", help="where the reads of a BAM file differ from the plasmid (SNPs, insertions, deletions, breaks); --apply writes a corrected copy")
    s.add_argument("file")
    s.add_argument("reads", help="a sorted, indexed BAM file mapped to this plasmid")
    s.add_argument("--seqid", help="the BAM's name for this plasmid (default: the one that matches its name, or the only one)")
    s.add_argument("--min-depth", dest="min_depth", type=int, help="fewest reads to count a position (default: the Alignments setting, 5)")
    s.add_argument("--variant-pct", dest="variant_pct", type=int, help="percent of reads that must differ (setting default 20)")
    s.add_argument("--consensus-pct", dest="consensus_pct", type=int, help="percent of reads that make a change a consensus change (setting default 80)")
    s.add_argument("--min-baseq", dest="min_baseq", type=int, help="lowest base quality counted (setting default 13)")
    s.add_argument("--min-mapq", dest="min_mapq", type=int, help="lowest mapping quality counted (setting default 0)")
    s.add_argument("--apply", help="write a copy of the plasmid with the consensus changes applied to this new GenBank file")
    s.set_defaults(run=cmd_variants)
    s = sub.add_parser("annotate", help="find known parts (promoters, terminators, origins, resistance genes...) in the plasmid (Detect Features); --add writes a copy with them annotated")
    s.add_argument("file")
    s.add_argument("--kinds", help="comma-separated feature types to look for (promoter, terminator, CDS, rep_origin...)")
    s.add_argument("--min-dna", dest="min_dna", type=int, help="minimum DNA identity in percent (the Detect Features setting, 95)")
    s.add_argument("--min-protein", dest="min_protein", type=int, help="minimum protein identity of a CDS in percent (90)")
    s.add_argument("--min-length", dest="min_length", type=int, help="ignore parts shorter than this many bp (15)")
    s.add_argument("--add", help="write a copy of the plasmid with the parts not annotated yet added as features, to this new GenBank file")
    s.set_defaults(run=cmd_annotate)
    s = sub.add_parser("checkup", help="for each CDS: should it be expressed? promoters, RBS, terminators, transcription units (Plasmid Checkup)")
    s.add_argument("file")
    s.add_argument("--min-promoter", dest="min_promoter", type=int, help="predicted promoter: minimum score (the Plasmid Checkup setting, 10)")
    s.add_argument("--min-stem", dest="min_stem", type=int, help="predicted terminator: minimum stem score (18)")
    s.add_argument("--rbs-reach", dest="rbs_reach", type=int, help="an RBS counts within this many bp of the start codon (20)")
    s.add_argument("--no-predict", action="store_true", help="use only annotated features and library parts, no motif prediction")
    s.add_argument("--no-library", action="store_true", help="do not use library parts")
    s.set_defaults(run=cmd_checkup)
    s = sub.add_parser("gibson", help="Gibson / HiFi assembly of fragments in a table (source file, region, orientation): primers with overlaps, checks and the assembled plasmid")
    s.add_argument("action", choices=["check", "run"])
    s.add_argument("fragments", help="CSV with the columns id,source (and optionally region,orientation), one row per fragment, in order")
    s.add_argument("--out", help="run: a new output folder (default: next to the table, FRAGMENTS_gibson)")
    s.add_argument("--name", default="assembly", help="run: name of the assembled plasmid's file (default assembly)")
    s.add_argument("--linear", action="store_true", help="the product is linear (default: circular, the last fragment joins the first)")
    s.add_argument("--min-overlap", dest="min_overlap", type=int, help="shortest overlap in bp (default 20)")
    s.add_argument("--max-overlap", dest="max_overlap", type=int, help="longest overlap in bp (default 40)")
    s.add_argument("--overlap-tm", dest="overlap_tm", type=float, help="the overlap's lowest Tm in deg C (default 52)")
    s.add_argument("--tm", type=float, help="annealing Tm of the primers in deg C (default 60)")
    s.set_defaults(run=cmd_gibson)
    s = sub.add_parser("ncbi", help="search NCBI Nucleotide or Protein, or fetch one record as a GenBank (nucleotide) or FASTA (protein) file; sends only the query to NCBI")
    s.add_argument("action", choices=["search", "get"])
    s.add_argument("query", help="search words, or (get) an accession such as NC_000913.3 or J01749")
    s.add_argument("--db", choices=["nuccore", "protein"], default="nuccore")
    s.add_argument("--limit", type=int, help="search: how many records (default 25)")
    s.add_argument("--out", help="get: the new file to write")
    s.add_argument("--allow-large", action="store_true", help="get: a record of more than 20 million bases")
    s.set_defaults(run=cmd_ncbi)
