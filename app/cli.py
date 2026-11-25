"""Command line for agents (and scripts): the same functions as the windows, with no window.

    oDNA-cli.exe <command> ...                      (installed: next to openDNA.exe)
    .venv\\Scripts\\python.exe -m app.cli <command>  (from the repository, run from its root)

Commands: `skills list|show NAME` (the recipes in app/agent/), `inspect FILE`, `feature FILE NAME`, `region FILE START END`;
batches with a table, `check` then `run`: `kld`, `pcr`, `recombineering`, `codon`, `gibson`; read-only analyses (app/cli_tools.py): `sites`, `guides`,
`repeats`, `search`, `align`, `msa`, `checkup`, `annotate`, `variants`, `synthesis`; files: `export`, `order`; network: `ncbi`, `hub` (a UCSC track hub); tracks: `recap`, `print`, `track-features`, `signal-regions`; optional tools beside openDNA (pLannotate, the Terminal's): `wsltools`.
Everything prints JSON (skills print Markdown). An error is one line `ERROR: ...` on stderr and exit code 1, never a traceback.
Positions are 1-based and inclusive, as in the GenBank file and HGVS (`region f.gb 10 20` is bases 10 to 20).
"""
import argparse
import json
import os
import sys
from pathlib import Path

# Some modules import Qt; nothing is ever shown. On Windows the "minimal" platform with the system's fonts: "offscreen" finds no
# font there, and the text of an image drawn by `print` came out as boxes.
os.environ.setdefault("QT_QPA_PLATFORM", "minimal:enable_fonts" if os.name == "nt" else "offscreen")
if os.name == "nt" and os.environ["QT_QPA_PLATFORM"].startswith("offscreen"):   # (offscreen asked for, as the tests do: its fonts)
    os.environ.setdefault("QT_QPA_FONTDIR", os.path.join(os.environ.get("WINDIR", r"C:\Windows"), "Fonts"))
sys.path.insert(0, str(Path(__file__).resolve().parent))   # the modules import each other by flat name

import cli_tools
from user_error import UserError   # a mistake the agent can fix: its message is printed as it is

AGENT_DIR = Path(__file__).with_name("agent")   # the recipes; the PyInstaller spec bundles the folder next to the modules


def qt_setup():
    """What the modules that read settings need: a Qt application with the program's name (so the data folder is the app's)."""
    from PySide6.QtWidgets import QApplication
    import project
    app = QApplication.instance() or QApplication(sys.argv[:1])
    app.setApplicationName(project.NAME)
    return app


def file_records(path):
    """[name] of the records of a FASTA file or of a GenBank file with several LOCUS lines; None for an ordinary one-record file."""
    from opendna_core import Fasta
    import multigb
    from reference import is_fasta_file
    if is_fasta_file(path):
        return [n for n, _ in Fasta(str(path)).sequences()]
    return [n for n, *_ in multigb.records(path)] if multigb.is_multi_genbank(path) else None


def said(doc, path):
    """`doc`, after a WARNING line on stderr for each thing the file held that was not read (doc.notes)."""
    for note in doc.notes:
        print(f"WARNING: {path}: {note}", file=sys.stderr)
    return doc


def load(path):
    """The Document of a file. A file of several records (multi-FASTA, multi-record GenBank) is read one record at a time: `FILE#NAME`
    (the file need not be split); a bare file name is accepted when it holds exactly one record."""
    from opendna_core import Document
    import multigb
    from reference import is_fasta_file
    path = str(path)
    record = None
    if not Path(path).is_file() and "#" in path:
        path, record = path.rsplit("#", 1)
    if not Path(path).is_file():
        raise UserError(f"File not found: {path}")
    try:
        names = file_records(path)
        if names is None:
            if record is not None:
                raise UserError(f"{path} holds one record: write the file name without #{record}.")
            return said(Document(path), path)
        listed = ", ".join(names[:8]) + (", ..." if len(names) > 8 else "")
        if record is None:
            if len(names) != 1:
                raise UserError(f"{path} holds {len(names)} records ({listed}): name one as {path}#NAME, or write one file per record with `split`.")
            record = names[0]
        if record not in names:
            raise UserError(f"{path} has no record called {record!r}. Its records: {listed}.")
        return Document.from_fasta(path, record) if is_fasta_file(path) else multigb.read_record(path, record)
    except UserError:
        raise
    except Exception as e:   # the core raises ValueError / OSError with its own words
        raise UserError(f"Cannot read {path}: {e}")


def feature_dna(doc, f):
    """The bases of a feature read in its own direction (a join: the parts in order), or None for a part that wraps the origin."""
    from mutation import reverse_complement
    if any(s < 0 or e > doc.length for s, e in f.spans):
        return None
    seq = "".join(doc.sequence(s, e) for s, e in f.spans)
    return reverse_complement(seq) if f.strand < 0 else seq


def protein(doc, f):
    """{frame, table, length, sequence} of a CDS, or {"note": why none}."""
    import genetic_code
    import mutation
    plan = genetic_code.cds_plan(doc, f)
    dna = feature_dna(doc, f)
    if dna is None or plan.why:
        return {"note": plan.why or "crosses the origin: not translated here"}
    aa = mutation.translate(dna[plan.frame:].upper(), genetic_code.codon_map(plan.table))
    return {"frame": plan.frame + 1, "table": plan.table, "length": len(aa.rstrip("*")), "sequence": aa}


def cmd_inspect(args):
    from mutation import arc
    doc = load(args.file)
    out = {"file": str(args.file), "name": doc.name, "length": doc.length, "circular": doc.circular, "features": []}
    for f in doc.features():
        item = {"index": f.index, "kind": f.kind, "label": f.label, **span_fields(doc, *arc(f, doc.length, doc.circular)), "strand": f.strand, "location": f.location}
        if f.kind == "CDS":
            item["protein"] = protein(doc, f)
        out["features"].append(item)
    return out


def cmd_feature(args):
    from mutation import arc, find_feature
    doc = load(args.file)
    f = find_feature(doc, args.name)
    dna = feature_dna(doc, f)
    if dna is None:
        raise UserError(f"{args.name!r} has a part outside the sequence.")
    out = {"index": f.index, "kind": f.kind, "label": f.label, **span_fields(doc, *arc(f, doc.length, doc.circular)), "strand": f.strand, "length": len(dna), "sequence": dna}
    if f.kind == "CDS":
        out["protein"] = protein(doc, f)
    return out


def cmd_region(args):
    doc = load(args.file)
    n = doc.length
    if doc.circular and 1 <= args.end < args.start <= n:   # across the origin: START..N then 1..END
        return {"start": args.start, "end": args.end, "length": n - args.start + 1 + args.end, "crosses_origin": True,
                "sequence": doc.sequence(args.start - 1, n) + doc.sequence(0, args.end)}
    if not 1 <= args.start <= args.end <= n:
        raise UserError(f"Positions must satisfy 1 <= START <= END <= {n}" + (" (or START > END across the origin of a circular sequence)." if doc.circular else " (this sequence is linear)."))
    return {"start": args.start, "end": args.end, "length": args.end - args.start + 1, "sequence": doc.sequence(args.start - 1, args.end)}


def given(value, default, limits=None, name="this value"):
    """A setting for the agent to show: the value in use, where it came from, and the range the windows accept (`limits`). A value
    given beyond that range is used as given, with a warning that the result is the sender's responsibility."""
    out = {"value": default if value is None else value, "source": "default" if value is None else "override (command line)"}
    if limits:
        out["min"], out["max"] = limits
        if value is not None and not limits[0] <= value <= limits[1]:
            out["warning"] = (f"{name} = {value} is outside what openDNA's windows accept ({limits[0]} to {limits[1]}): it is used as given, "
                              "and the result is your responsibility.")
    return out


def setting_warnings(settings):
    """The warnings found anywhere in a settings block."""
    if isinstance(settings, dict):
        return ([settings["warning"]] if "warning" in settings else []) + [w for v in settings.values() for w in setting_warnings(v)]
    return []


def effective_settings(extra):
    """Every setting that changes the answer, with its value and where it came from: for the agent to show before `run`. `extra`: the
    settings of the command itself."""
    import codon_usage
    import genetic_code
    import primer_panel
    import tool_settings
    from primer_db import BINDING, PCR
    saved = tool_settings.store()
    sourced = lambda tool, keys: {k: {"value": tool_settings.value(tool, k), "source": "saved in openDNA" if saved.contains(f"{tool}/{k}") else "default"}
                                  for k in keys}
    usage = codon_usage.active()
    seed, mismatches = primer_panel.homology_settings()
    genome = primer_panel.genome_path()
    return {
        **extra,
        "pcr_conditions": sourced(PCR, primer_panel.conditions()),
        "genetic_code_table": {"value": genetic_code.table_id(), "source": "openDNA settings"},
        "codon_usage_table": {"value": usage.name if usage else None,
                              "source": "openDNA settings" if usage else "none: ties between codons are broken alphabetically"},
        "primer_binding": sourced(BINDING, ("seed", "exact_3p", "mismatches", "min_tm")),
        "primer_checks": {"homology_seed": {"value": seed, "source": "openDNA settings"}, "homology_mismatches": {"value": mismatches, "source": "openDNA settings"},
                          "background_genome": {"value": genome, "source": "openDNA settings (Primers > Options)" if genome else "none set: no genome check"}},
    }


RESERVED = {"con", "prn", "aux", "nul", *(f"com{i}" for i in range(1, 10)), *(f"lpt{i}" for i in range(1, 10))}   # names Windows will not make a file of


def read_table(path, *columns, optional=()):
    """Rows of a CSV with the columns id and `columns` (and the `optional` ones; nothing else). A wrong header stops the run (the table
    cannot be read at all). A malformed row does not: it comes back with an `_error` saying what is wrong, and the other rows go on."""
    import csv
    import io
    if not Path(path).is_file():
        raise UserError(f"File not found: {path}")
    raw = Path(path).read_bytes()
    try:
        text = raw.decode("utf-8-sig")
    except UnicodeDecodeError:   # a table saved by Excel
        text = raw.decode("cp1252", errors="replace")
    reader = csv.reader(io.StringIO(text))
    header, need = [h.strip() for h in next(reader, [])], ["id", *columns]
    if len(set(header)) != len(header) or not set(need) <= set(header) or not set(header) <= set(need) | set(optional):
        raise UserError(f"{path}: the columns must be exactly {','.join(need)}" + (f" (and optionally {','.join(optional)})" if optional else "")
                        + f"; found: {', '.join(header) or 'none'}.")
    rows = []
    for line, cells in enumerate(reader, 2):   # (line 1 is the header)
        if not any(c.strip() for c in cells):
            continue
        row = {"_line": line, **dict.fromkeys(header, ""), **dict(zip(header, (c.strip() for c in cells)))}
        row["id"] = row["id"] or f"line{line}"
        if len(cells) != len(header):
            row["_error"] = (f"line {line} of {path} has {len(cells)} cells but the header has {len(header)} ({','.join(header)}): the cells found are "
                             f"{[c.strip() for c in cells]!r}. A cell that contains a comma must be in double quotes. This row was skipped; the others were processed.")
        elif empty := [c for c in need if not row[c]]:
            row["_error"] = f"line {line} of {path}: {', '.join(empty)} is empty. This row was skipped; the others were processed."
        rows.append(row)
    if not rows:
        raise UserError(f"{path}: no rows.")
    return rows


def assign_files(rows):
    """`file_stem` of every row: its id made safe for a file name (letters, digits, . _ -), and never the name of another row, in any
    letter case (a Windows or macOS disk would make them one file): `m1`, `M1` -> m1, M1_2. A changed name is said in the row's warnings."""
    import re
    taken = set()
    for row in rows:
        stem = re.sub(r"[^A-Za-z0-9._-]", "_", row["id"]).strip(".") or f"line{row['_line']}"
        stem = f"_{stem}" if stem.split(".")[0].lower() in RESERVED else stem
        name, n = stem, 1
        while name.lower() in taken:
            n += 1
            name = f"{stem}_{n}"
        taken.add(name.lower())
        row["file_stem"] = name
        if name != row["id"]:
            row["file_note"] = f"The id {row['id']!r} is not usable as a file name here (unsafe characters, or the same name as another row): the file is {name}.gb."


def run_rows(rows, row_fn, key):
    """A record for every row, whatever happens: a malformed row, or one that fails in an unforeseen way (a bug: it is reported with its
    traceback), is an error record and the batch goes on."""
    import traceback
    records = []
    for row in rows:
        if "_error" in row:
            rec = {"id": row["id"], key: row.get(key, ""), "status": "error", "error": row["_error"]}
        else:
            try:
                rec = row_fn(row)
            except Exception as e:   # one row's unforeseen failure must not end the batch
                rec = {"id": row["id"], key: row.get(key, ""), "status": "error", "detail": traceback.format_exc(),
                       "error": f"Unexpected {type(e).__name__}: {e} (a bug in openDNA: please report it with this row)"}
        if rec["status"] == "ok" and row.get("file_note"):
            rec.setdefault("warnings", []).append(row["file_note"])
        records.append(rec)
    return records


def make_out_dir(args, source, suffix):
    """None for `check`; for `run` a new, empty folder (--out, or next to `source`): nothing is ever overwritten."""
    if args.action != "run":
        return None
    out_dir = Path(args.out) if args.out else Path(source).with_name(Path(source).stem + suffix)
    if out_dir.exists() and any(out_dir.iterdir()):
        raise UserError(f"The output folder {out_dir} already exists and is not empty: choose a new one with --out (nothing is overwritten).")
    out_dir.mkdir(parents=True, exist_ok=True)
    return out_dir


def finish(out_dir, settings, records, header, csv_rows, warned):
    """The answer of a batch: settings, summary and every record; with `out_dir` also primers.csv (`header`, and the lines `csv_rows`
    gives for each record), report.jsonl and settings.json, after checking that every file a record names is there. (result, exit code):
    1 when a row failed, the others are still reported. Every error and settings warning is also printed on stderr."""
    import csv
    if out_dir is not None:
        for r in records:
            if r["status"] == "ok" and r.get("file") and not (out_dir / r["file"]).is_file():
                r.update(status="error", error=f"The file {r['file']} was expected in {out_dir} and is not there.")
    failed = [r for r in records if r["status"] == "error"]
    result = {"settings": settings, "warnings": setting_warnings(settings),
              "summary": {"rows": len(records), "ok": len(records) - len(failed), "errors": len(failed), "with_warnings": sum(map(warned, records)),
                          **({"files_written": sum(1 for r in records if r.get("file"))} if out_dir is not None else {})}, "rows": records}
    if out_dir is not None:
        with open(out_dir / "primers.csv", "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(header)
            for r in records:
                if r["status"] == "ok":
                    w.writerows(csv_rows(r))
        (out_dir / "report.jsonl").write_text("".join(json.dumps(r) + "\n" for r in records), encoding="utf-8")
        (out_dir / "settings.json").write_text(json.dumps(settings, indent=1), encoding="utf-8")
        result["out"] = str(out_dir)
    for w in result["warnings"]:
        print(f"WARNING: {w}", file=sys.stderr)
    for r in failed:
        print(f"ERROR [row {r['id']}]: {r['error']}" + (f"\n{r['detail']}" if r.get("detail") else ""), file=sys.stderr)
    return result, 1 if failed else 0


def primer_records(rows):
    """The primers of a design as records: rows are (name, sequence, annealing length, ...)."""
    import primer_panel
    return [{"name": n, "sequence": s, "annealing_length": a, "annealing_tm": primer_panel.melting_temp(s[-a:])} for n, s, a, *_ in rows]


def tm_of(primer):
    return round(primer["annealing_tm"], 1) if primer["annealing_tm"] else ""


def primer_pair_line(r, *lead):
    """primers.csv line of a record with one forward and one reverse primer (a deletion, or a mutation that changes nothing, has none)."""
    f, rv = (r["primers"] + [{"sequence": "", "annealing_tm": None}] * 2)[:2]
    return [*lead, f["sequence"], tm_of(f), rv["sequence"], tm_of(rv)]


def kld_row(source, row, args, usage, out_dir):
    """One mutation: the result record, and the edited plasmid with its primers written to `out_dir` (when there is one)."""
    import kld
    import mutation
    from assembly import snapshot
    rec = {"id": row["id"], "mutation": row["mutation"], "status": "error"}
    try:
        r = mutation.resolve(source, row["mutation"], usage)
        seq = source.sequence(0, None).upper()
        d = None if r.effect == "no change" else kld.design_selected(
            seq[r.shift:] + seq[:r.shift], [(r.start, r.end, r.text)], source.circular,   # (a change through the origin is designed on the plasmid turned by r.shift)
            kld.DEFAULT_TM if args.tm is None else args.tm, kld.DEFAULT_GROUP if args.group is None else args.group,
            kld.DEFAULT_WARN if args.warn is None else args.warn)
    except (UserError, ValueError) as e:   # (design raises ValueError with a sentence: too short, too many new bases, ...)
        rec["error"] = str(e)
        return rec
    notes = d.notes if d else []
    rec.update(status="ok", feature=r.feature, hgvs_dna=r.hgvs_dna, hgvs_protein=r.hgvs_protein, effect=r.effect,
               warnings=r.warnings + [t for level, t in notes if level == "warn"], notes=[t for level, t in notes if level != "warn"],
               primers=primer_records(d.rows) if d else [])
    if out_dir is not None:   # (a mutation that changes nothing still gets its file: the plasmid as it is)
        doc = snapshot(source, source.name)
        if d:
            if r.shift:
                doc.rotate(r.shift)
            doc.replace(r.start, r.end, r.text)
            doc = kld.with_primers(doc, d)
            if r.shift:   # back to the file's origin: its first base stays first, unless the change replaced it (then the first base after the change)
                j = (-r.shift) % source.length
                first = (r.start + len(r.text) + max(0, j - r.end) if j >= r.start else j) % doc.length
                if first:
                    doc.rotate(first)
        doc.write_genbank(str(out_dir / f"{row['file_stem']}.gb"), None)
        rec["file"] = f"{row['file_stem']}.gb"
    return rec


def cmd_kld(args):
    import codon_usage
    import kld
    source, rows = load(args.plasmid), read_table(args.mutations, "mutation")
    assign_files(rows)
    out_dir = make_out_dir(args, args.plasmid, "_kld")
    lim = kld.LIMITS
    settings = effective_settings({"kld": {"tm": given(args.tm, kld.DEFAULT_TM, lim["tm"], "--tm"), "group": given(args.group, kld.DEFAULT_GROUP, lim["group"], "--group"),
                                           "warn_insert": given(args.warn, kld.DEFAULT_WARN, lim["warn_insert"], "--warn")}})
    usage = codon_usage.active()
    records = run_rows(rows, lambda row: kld_row(source, row, args, usage, out_dir), "mutation")
    return finish(out_dir, settings, records,
                  ["id", "mutation", "file", "hgvs_dna", "hgvs_protein", "forward_primer", "forward_tm", "reverse_primer", "reverse_tm", "warnings"],
                  lambda r: [primer_pair_line(r, r["id"], r["mutation"], r.get("file", ""), r["hgvs_dna"], r["hgvs_protein"]) + [" | ".join(r["warnings"])]],
                  lambda r: bool(r.get("warnings")))


def target_range(doc, text):
    """(first, last), 1-based inclusive, of a feature (name or #INDEX) or of START..END. On a circular sequence a target through the origin has
    `last` above the length (START..END with START > END, or a feature that crosses it): the bases first..length then 1..last - length."""
    import re
    from mutation import arc, find_feature
    m = re.fullmatch(r"(\d+)\.\.(\d+)", text)
    if not m:
        return arc(find_feature(doc, text), doc.length, doc.circular)
    first, last = int(m.group(1)), int(m.group(2))
    if doc.circular and 1 <= last < first <= doc.length:
        return first, last + doc.length
    if not 1 <= first < last <= doc.length:
        raise UserError(f"{text!r}: positions must satisfy 1 <= START < END <= {doc.length}" + (" (or START > END across the origin)." if doc.circular else " (this sequence is linear)."))
    return first, last


def span_fields(doc, first, last):
    """The `start`, `end` (1-based, `end` below `start` when the target runs through the origin) of a target from `target_range`."""
    n = doc.length
    return {"start": first, "end": last - n if last > n else last, **({"crosses_origin": True} if last > n else {})}


def pcr_row(doc, row, args):
    """The best primer pairs for one target (a feature name or START..END, 1-based inclusive): the record of the row."""
    import primer_design
    import primer_panel
    rec = {"id": row["id"], "target": row["target"], "status": "error"}
    try:
        first, last = target_range(doc, row["target"])
        opts = {k: v for k, v in dict(tm_target=args.tm, min_len=args.min_len, max_len=args.max_len, window=args.window).items() if v is not None}
        pairs = primer_design.design_region(doc, first - 1, last, **opts)
    except UserError as e:
        rec["error"] = str(e)
        return rec
    if not pairs:
        rec["error"] = "No pair found: widen the search window (--window), the length range (--min-len, --max-len) or the Tm (--tm)."
        return rec
    ed = primer_design.editor_like(doc)
    out = []
    for rank, p in enumerate(pairs, 1):
        checks = {k: primer_panel.primer_checks(ed, p[k]["seq"], "wait") for k in ("fwd", "rev")}
        out.append({"rank": rank, "product_size": p["size"], "penalty": round(p["penalty"], 2), "notes": p["notes"],
                    "mark": primer_panel.worst(checks["fwd"] + checks["rev"]),
                    **{k: {"sequence": p[k]["seq"], "tm": round(p[k]["tm"], 1), "gc": round(p[k]["gc"]), "start": p[k]["start"] + 1, "end": p[k]["end"],
                           "checks": [{"check": n, "level": lv, "text": t} for n, lv, t in checks[k]]} for k in ("fwd", "rev")}})
    rec.update(status="ok", **span_fields(doc, first, last), pairs=out)
    if pairs.more:
        rec["warnings"] = [f"Only the best {len(pairs)} pairs are listed: other pairs were found."]
    if pairs.window is not None:
        rec.setdefault("warnings", []).append(f"The search window was narrowed to {pairs.window} bases on each side: on this circular plasmid a wider one would make the two sides meet.")
    if all(p["mark"] == "bad" for p in out):   # (the window is the user's choice: only offer a wider one)
        window = opts.get("window", primer_design.DEFAULT_WINDOW)
        rec.setdefault("warnings", []).append(f"Every pair is marked bad with --window {window}. Tell the user, and offer to try a wider search window"
                           f"{' (for example --window 30)' if window < 30 else ''}: the primers may then start or end outside the target, so the product is longer than the target. "
                           f"Do not widen it without their OK.")
    return rec


def cmd_pcr(args):
    import primer_design
    doc, rows = load(args.plasmid), read_table(args.targets, "target")
    out_dir = make_out_dir(args, args.plasmid, "_pcr")
    lim = primer_design.LIMITS
    settings = effective_settings({"pcr_design": {
        "tm": given(args.tm, primer_design.DEFAULT_TM, lim["tm"], "--tm"), "min_len": given(args.min_len, primer_design.DEFAULT_MIN_LEN, lim["min_len"], "--min-len"),
        "max_len": given(args.max_len, primer_design.DEFAULT_MAX_LEN, lim["max_len"], "--max-len"),
        "window": given(args.window, primer_design.DEFAULT_WINDOW, lim["window"], "--window")}})
    records = run_rows(rows, lambda row: pcr_row(doc, row, args), "target")
    result, code = finish(out_dir, settings, records,
                          ["id", "target", "rank", "forward_primer", "forward_tm", "reverse_primer", "reverse_tm", "product_size", "penalty", "mark", "notes"],
                          lambda r: [[r["id"], r["target"], p["rank"], p["fwd"]["sequence"], p["fwd"]["tm"], p["rev"]["sequence"], p["rev"]["tm"], p["product_size"],
                                      p["penalty"], p["mark"], " | ".join(p["notes"])] for p in r["pairs"]],
                          lambda r: any(p["mark"] in ("warn", "bad") for p in r.get("pairs", [])))
    if not args.full:   # (report.jsonl, written above, keeps every check)
        for r in records:
            for p in r.get("pairs", []):
                for k in ("fwd", "rev"):
                    checks = p[k].pop("checks")
                    p[k]["checks_passed"] = sum(c["level"] == "ok" for c in checks)
                    p[k]["issues"] = [c for c in checks if c["level"] != "ok"]   # (only what is not fine: with --full, all the checks)
    return result, code


def rec_row(genome, row, args, folder, cassettes, out_dir):
    """One recombineering design: the cassette (a GenBank file, or `none` for a deletion) replaces a region of the genome."""
    import recombineering as R
    from assembly import transfer_primers
    from opendna_core import Document
    from types import SimpleNamespace
    rec = {"id": row["id"], "region": row["region"], "cassette": row["cassette"], "status": "error"}
    writing = out_dir is not None and not args.primers_only   # (only a written genome needs the recombinant: a bacterial genome is costly to copy)
    try:
        first, last = target_range(genome, row["region"])
        orientation = row.get("orientation", "").lower() or "forward"
        if orientation not in ("forward", "reverse"):
            raise UserError(f"orientation must be forward or reverse, not {orientation!r}.")
        if row["cassette"].lower() in ("none", "-"):   # a deletion
            cassette = Document.new_sequence("none", "", False)
        else:
            path = Path(row["cassette"])
            path = str(path if path.is_absolute() else folder / path)
            if path not in cassettes:
                cassettes[path] = load(path)
            cassette = cassettes[path]
        left = R.HOMOLOGY if args.homology is None else args.homology
        d = (R.design if writing else R.design_primers)(genome, cassette, first - 1, last, left, left, orientation == "reverse",
                                                         R.DEFAULT_TM if args.tm is None else args.tm)
    except (UserError, ValueError) as e:
        rec["error"] = str(e)
        return rec
    warnings = []
    covered = [(first - 1, min(last, genome.length))] + ([(0, last - genome.length)] if last > genome.length else [])   # (a region through the origin is two pieces)
    for f in genome.features():   # (by the parts of a feature: one that joins the two ends of the plasmid is not "all of it")
        if f.kind == "source" or not any(s < b and a < e for s, e in f.spans for a, b in covered):
            continue
        if all(any(a <= s and e <= b for a, b in covered) for s, e in f.spans):
            warnings.append(f"The region removes the {f.kind} {f.label!r} ({f.location}).")
        else:
            warnings.append(f"The region cuts the {f.kind} {f.label!r} ({f.location}): only part of it is replaced.")
    rec.update(status="ok", **span_fields(genome, first, last), replaced_bp=last - first + 1, cassette_bp=cassette.length, orientation=orientation,
               left_homology=d.left, right_homology=d.right, warnings=warnings, primers=primer_records(d.rows))
    if writing:
        if d.rows:
            transfer_primers(d.product, [SimpleNamespace(fwd=d.rows[0][:2], rev=d.rows[1][:2])])
        d.product.write_genbank(str(out_dir / f"{row['file_stem']}.gb"), None)
        rec["file"] = f"{row['file_stem']}.gb"
    return rec


def cmd_recombineering(args):
    import recombineering as R
    genome = load(args.genome)
    rows = read_table(args.table, "region", "cassette", optional=("orientation",))
    assign_files(rows)
    out_dir = make_out_dir(args, args.genome, "_rec")
    lim = R.LIMITS
    settings = effective_settings({"recombineering": {
        "homology": given(args.homology, R.HOMOLOGY, lim["homology"], "--homology"), "annealing_tm": given(args.tm, R.DEFAULT_TM, lim["annealing_tm"], "--tm"),
        "write_genomes": given(False if args.primers_only else None, True)}})
    cassettes = {}
    folder = Path(args.table).resolve().parent
    records = run_rows(rows, lambda row: rec_row(genome, row, args, folder, cassettes, out_dir), "region")
    return finish(out_dir, settings, records,
                  ["id", "region", "cassette", "orientation", "file", "forward_primer", "forward_tm", "reverse_primer", "reverse_tm", "warnings"],
                  lambda r: [primer_pair_line(r, r["id"], r["region"], r["cassette"], r["orientation"], r.get("file", "")) + [" | ".join(r["warnings"])]],
                  lambda r: bool(r.get("warnings")))


def recipes():
    """{name: path} of the recipes (the index is `index`)."""
    found = {p.stem: p for p in sorted((AGENT_DIR / "skills").glob("*.md"))}
    return {"index": AGENT_DIR / "index.md", **found}


def cmd_skills(args):
    names = recipes()
    if args.action == "list":
        return "\n".join(f"{n}: " + next((ln.lstrip("# ").strip() for ln in p.read_text(encoding="utf-8").splitlines() if ln.strip()), "") for n, p in names.items())
    if args.name not in names:
        raise UserError(f"No recipe called {args.name!r}. Choices: {', '.join(names)}.")
    return names[args.name].read_text(encoding="utf-8")


def parser():
    p = argparse.ArgumentParser(prog="oDNA-cli", description=__doc__.split("\n\n")[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="command", required=True)
    s = sub.add_parser("skills", help="the recipes for agents: `skills list`, `skills show index` first")
    s.add_argument("action", choices=["list", "show"])
    s.add_argument("name", nargs="?", default="index")
    s.set_defaults(run=cmd_skills, text=True)
    s = sub.add_parser("inspect", help="features of a GenBank file with positions, strand and the protein of each CDS")
    s.add_argument("file")
    s.set_defaults(run=cmd_inspect)
    s = sub.add_parser("feature", help="the sequence (and protein) of one feature, by label or #index")
    s.add_argument("file")
    s.add_argument("name")
    s.set_defaults(run=cmd_feature)
    s = sub.add_parser("kld", help="KLD primers and edited plasmids for a table of HGVS mutations (check: report only; run: also write files)")
    s.add_argument("action", choices=["check", "run"])
    s.add_argument("plasmid")
    s.add_argument("mutations", help="CSV with the columns id,mutation")
    s.add_argument("--out", help="run: a new output folder (default: PLASMID_kld next to the plasmid)")
    s.add_argument("--tm", type=float, help="annealing Tm in deg C (default 60)")
    s.add_argument("--group", type=int, help="changes closer than this many bp share one primer pair (default 30)")
    s.add_argument("--warn", type=int, help="warn for insertions above this many bp (default 120)")
    s.set_defaults(run=cmd_kld, returns_code=True)
    s = sub.add_parser("pcr", help="PCR primer pairs for a table of targets (a feature name or START..END), with the Primers panel's checks")
    s.add_argument("action", choices=["check", "run"])
    s.add_argument("plasmid")
    s.add_argument("targets", help="CSV with the columns id,target")
    s.add_argument("--out", help="run: a new output folder (default: PLASMID_pcr next to the plasmid)")
    s.add_argument("--tm", type=float, help="target Tm in deg C (default 60)")
    s.add_argument("--min-len", type=int, help="shortest primer (default 18)")
    s.add_argument("--max-len", type=int, help="longest primer (default 25)")
    s.add_argument("--window", type=int, help="how far from the target ends a primer may start, in bp (default 0: on the ends)")
    s.add_argument("--full", action="store_true", help="print every check of every primer (default: a summary: sequence, Tm, GC, position and the checks that are not fine)")
    s.set_defaults(run=cmd_pcr, returns_code=True)
    s = sub.add_parser("recombineering", help="lambda Red: primers (homology tails) and the recombinant genome for cassettes that replace regions")
    s.add_argument("action", choices=["check", "run"])
    s.add_argument("genome")
    s.add_argument("table", help="CSV with the columns id,region,cassette (and optionally orientation)")
    s.add_argument("--out", help="run: a new output folder (default: GENOME_rec next to the genome)")
    s.add_argument("--homology", type=int, help="bp of homology on each side of the cassette (default 50)")
    s.add_argument("--tm", type=float, help="annealing Tm in deg C (default 60)")
    s.add_argument("--primers-only", action="store_true", help="run: do not write the recombinant genomes (a bacterial genome is megabytes each)")
    s.set_defaults(run=cmd_recombineering, returns_code=True)
    cli_tools.register(sub, sys.modules[__name__])
    s = sub.add_parser("region", help="the bases START..END (1-based, inclusive)")
    s.add_argument("file")
    s.add_argument("start", type=int)
    s.add_argument("end", type=int)
    s.set_defaults(run=cmd_region)
    return p


def main(argv=None):
    args = parser().parse_args(argv)
    for stream in (sys.stdout, sys.stderr):   # a Windows console is not UTF-8: primer notes and recipes have symbols
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(encoding="utf-8")
    try:
        if args.command != "skills":
            qt_setup()   # (genetic_code reads the settings file: the same table choice as the windows)
        result = args.run(args)
    except UserError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1
    code = 0
    if getattr(args, "returns_code", False):   # a batch: failed rows make the exit code 1, the rest is still reported
        result, code = result
    print(result if getattr(args, "text", False) else json.dumps(result, indent=1))
    return code


if __name__ == "__main__":
    sys.exit(main())
