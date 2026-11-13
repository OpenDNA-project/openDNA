"""Build the unit-test report in docs/test-report/ (see its README.md).

For every test (Python: each `def test_*` in tests/test_*.py; Rust: each `#[test]` fn in core/src/*.rs) the report gives what goes in,
what is expected to come out, the kind of reference the answer is checked against, a strength verdict, review notes and the
result of the last run. The written descriptions live in docs/test-report/tests.json; a test that has none (a new test, say) is
still listed, with its docstring and its assert lines taken from the code, and marked as such.

    python -m pytest tests -q --junitxml=junit.xml                     # QT_QPA_PLATFORM=offscreen on a machine without a screen
    cd core && cargo test > ../cargo.txt 2>&1; cd ..
    python scripts/test_report.py --junit junit.xml --cargo cargo.txt  # writes docs/test-report/*.md and report.html

Without --junit / --cargo the result column says "not run".
"""
import argparse
import ast
import datetime
import json
import re
import subprocess
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "docs" / "test-report"
ICON = {"passed": "✅ passed", "failed": "❌ failed", "skipped": "⏭️ skipped", "not run": "· not run"}
STRENGTH = {"strong": "🟢 strong", "adequate": "🟡 adequate", "weak": "🔴 weak"}


# ---- the tests in the code ---------------------------------------------------------------------------------------------------
def python_tests():
    out = []
    for path in sorted((ROOT / "tests").glob("test_*.py")):
        source = path.read_text(encoding="utf-8")
        for node in ast.parse(source).body:
            if isinstance(node, ast.FunctionDef) and node.name.startswith("test"):
                checks = [ast.unparse(n).split("\n")[0] for n in ast.walk(node) if isinstance(n, ast.Assert)]
                checks += [ast.unparse(n.items[0].context_expr) for n in ast.walk(node)
                           if isinstance(n, ast.With) and "pytest.raises" in ast.unparse(n.items[0].context_expr)]
                out.append({"id": f"tests/{path.name}::{node.name}", "file": f"tests/{path.name}", "name": node.name,
                            "line": node.lineno, "doc": ast.get_docstring(node) or "", "checks": checks})
    return out


def rust_tests():
    out = []
    for path in sorted((ROOT / "core" / "src").glob("*.rs")):
        lines = path.read_text(encoding="utf-8").splitlines()
        for i, line in enumerate(lines):
            if line.strip() != "#[test]":
                continue
            j = i + 1
            while not re.match(r"\s*fn ", lines[j]):
                j += 1
            name = re.match(r"\s*fn (\w+)", lines[j]).group(1)
            doc, k = [], i - 1
            while k >= 0 and lines[k].strip().startswith("///"):
                doc.insert(0, lines[k].strip()[3:].strip())
                k -= 1
            body, depth, m = [], 0, j
            while m < len(lines):
                depth += lines[m].count("{") - lines[m].count("}")
                body.append(lines[m])
                if depth <= 0 and "{" in "".join(body):
                    break
                m += 1
            checks = [b.strip() for b in body if re.search(r"\bassert(_eq|_ne)?!", b)]
            out.append({"id": f"core/src/{path.name}::tests::{name}", "file": f"core/src/{path.name}", "name": name,
                        "line": j + 1, "doc": " ".join(doc), "checks": checks})
    return out


# ---- results of the last run -------------------------------------------------------------------------------------------------
def junit_results(path):
    results = {}
    for case in ET.parse(path).getroot().iter("testcase"):
        file = "tests/" + case.get("classname", "").split(".")[-1] + ".py"
        status = "passed"
        message = ""
        for child in case:
            if child.tag in ("failure", "error"):
                status, message = "failed", (child.get("message") or "")[:300]
            elif child.tag == "skipped":
                status, message = "skipped", (child.get("message") or "")[:300]
        results[f"{file}::{case.get('name')}"] = (status, float(case.get("time", 0)), message)
    return results


def cargo_results(path):
    results = {}
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        m = re.match(r"test (\w+)::tests::(\w+) \.\.\. (ok|FAILED|ignored)", line)
        if m:
            status = {"ok": "passed", "FAILED": "failed", "ignored": "skipped"}[m.group(3)]
            results[f"core/src/{m.group(1)}.rs::tests::{m.group(2)}"] = (status, None, "")
    return results


# ---- writing -----------------------------------------------------------------------------------------------------------------
def cell(text):
    return str(text).replace("|", "\\|").replace("\n", "<br>")


def page_name(file):
    return file.replace("/", "_").replace(".py", "").replace(".rs", "") + ".md"


def write(tests, described, results):
    OUT.mkdir(parents=True, exist_ok=True)
    by_file = defaultdict(list)
    for t in tests:
        by_file[t["file"]].append(t)
    index = []
    totals = defaultdict(int)
    for file, group in by_file.items():
        lines = [f"# `{file}`", "", f"[← all test files](README.md) · source: [`{file}`](../../{file})", ""]
        counts = defaultdict(int)
        for t in group:
            d = described.get(t["id"])
            status, seconds, message = results.get(t["id"], ("not run", None, ""))
            counts[status] += 1
            totals[status] += 1
            strength = d.get("strength", "") if d else ""
            totals["strength:" + (strength or "not reviewed")] += 1
            counts["strength:" + (strength or "not reviewed")] += 1
            took = f" in {seconds:.2f} s" if seconds else ""
            lines += [f"## `{t['name']}`", "",
                      f"**Result:** {ICON[status]}{took} · [line {t['line']}](../../{file}#L{t['line']})"
                      + (f" · **Area:** {d['area']}" if d and d.get("area") else "")
                      + (f" · **Checked against:** {d['oracle']}" if d and d.get("oracle") else "")
                      + (f" · **Strength:** {STRENGTH.get(strength, strength)}" if strength else ""), ""]
            if message:
                lines += [f"> {cell(message)}", ""]
            if d:
                if d.get("role"):
                    lines += [f"**Why it matters:** {d['role']}", ""]
                lines += [d.get("purpose", ""), "", "| # | Input | Expected output |", "|---|---|---|"]
                lines += [f"| {i} | {cell(c.get('input', ''))} | {cell(c.get('expected', ''))} |" for i, c in enumerate(d.get("cases", []), 1)]
                lines.append("")
                if d.get("verdict") == "gap":
                    lines += ["**Does not fully prove its claim:**", ""] + [f"- {cell(x)}" for x in d.get("gaps", [])] + [""]
                    if d.get("proposal"):
                        lines += [f"**Proposed bulletproof test:** {d['proposal']}", ""]
                if d.get("issues"):
                    lines += ["**Review notes:**", ""] + [f"- {cell(x)}" for x in d["issues"]] + [""]
                if d.get("fixed"):
                    lines += ["**Fixed by the audit of November 2026:**", ""] + [f"- {cell(x)}" for x in d["fixed"]] + [""]
            else:
                lines += ["*No written description yet: the docstring and the checks below are taken from the code as they are.*", ""]
                if t["doc"]:
                    lines += [t["doc"], ""]
                if t["checks"]:
                    lines += ["Checks (input → expected output, as written in the test):", "", "```" + ("python" if file.endswith(".py") else "rust")]
                    lines += t["checks"] + ["```", ""]
        (OUT / page_name(file)).write_text("\n".join(lines), encoding="utf-8")
        verdicts = " · ".join(f"{STRENGTH.get(k[9:], k[9:])} {v}" for k, v in sorted(counts.items()) if k.startswith("strength:"))
        index.append(f"| [`{file}`]({page_name(file)}) | {len(group)} | {counts['passed']} | {counts['failed']} | "
                     f"{counts['skipped']} | {counts['not run']} | {verdicts} |")
    n_py = sum(t["file"].endswith(".py") for t in tests)
    strengths = " · ".join(f"{STRENGTH.get(k[9:], k[9:])}: {v}" for k, v in sorted(totals.items()) if k.startswith("strength:"))
    readme = [
        "# Unit-test report", "",
        "What every automated test of openDNA feeds in, what it expects to come out, what the answer is checked against, and",
        "how much a pass is worth. Written so that a user can check, test by test, that the tests check relevant things and check",
        "them well. **[AUDIT.md](AUDIT.md) lists every test on one page with why it matters, its strength and the review notes.**", "",
        "## How to read a test entry", "",
        "- **Result**: the outcome of the last run this report was built from (✅ passed, ❌ failed, ⏭️ skipped, · not run).",
        "- **Checked against**: where the expected answer comes from. *hand-computed*: worked out on paper and written in the test;"
        " *independent implementation*: another program (mostly Biopython) computes it; *in-test reference model*: a plain model"
        " written inside the test (e.g. the sequence as a Python string); *property/invariant*: a rule that must always hold"
        " (e.g. undo restores the original); *round-trip/self-consistency*: save then reload gives the same thing; *GUI state*:"
        " what the window shows; *smoke*: only that nothing crashes.",
        "- **Strength**: 🟢 strong = exact values or an independent check, edge cases included; 🟡 adequate = the main behaviour,"
        " some edges or exact values missing; 🔴 weak = a clearly wrong program could still pass.",
        "- **Why it matters**: what would break unnoticed for a user if the test did not exist. A test whose reason you cannot"
        " agree with is a candidate for removal.",
        "- **Input / Expected output**: the concrete cases, in the order the test runs them. Many tests are *topical*: one test"
        " function holds several sections, each a row.",
        "- **Review notes**: what the reviewer found could be stronger. **Fixed by the audit**: what was changed.", "",
        f"## Totals ({len(tests)} tests: {n_py} Python, {len(tests) - n_py} Rust)", "",
        f"Last run: ✅ {totals['passed']} passed · ❌ {totals['failed']} failed · ⏭️ {totals['skipped']} skipped · · {totals['not run']} not run",
        "", f"Review: {strengths}", "",
        "## Test files", "",
        "| Test file | Tests | Passed | Failed | Skipped | Not run | Strength (reviewed) |", "|---|---|---|---|---|---|---|"] + index + [
        "", "## Rebuilding this report", "",
        "The descriptions are in [tests.json](tests.json); the pages are written by `scripts/test_report.py`, which also lists any",
        "test that has no description yet (with its docstring and asserts taken from the code), so a new test never goes missing:", "",
        "```bash", "QT_QPA_PLATFORM=offscreen python -m pytest tests -q --junitxml=junit.xml",
        "(cd core && cargo test > ../cargo.txt 2>&1)", "python scripts/test_report.py --junit junit.xml --cargo cargo.txt", "```", ""]
    (OUT / "README.md").write_text("\n".join(readme), encoding="utf-8")
    return totals


def write_audit(tests, described, results):
    """AUDIT.md: every test on one page with its role, reference, strength and review notes, built from tests.json and the run."""
    def suite(t):
        return "Python" if t["file"].endswith(".py") else "Rust"

    def strength(t):
        return (described.get(t["id"]) or {}).get("strength") or "not reviewed"

    kinds = ["strong", "adequate", "weak", "not reviewed"]
    lines = ["# Audit of the unit tests", "", "[← unit-test report](README.md)", "",
             f"Built by `scripts/test_report.py` on {datetime.date.today().isoformat()} from `tests.json` and the last run."
             " Every test of the project is listed below, Python (`tests/`) and Rust (`core/src`), with what would break"
             " unnoticed without it, what its expected answer is checked against, its strength and the reviewer's notes."
             " The Python suite holds no UI tests: windows, dialogs and drawing are checked by hand before a release.", "",
             "## Totals", "", "| Suite | Tests | Passed | Failed / not run | " + " | ".join(k.capitalize() for k in kinds) + " | With review notes |",
             "|---|---|---|---|" + "---|" * len(kinds) + "---|"]
    for name in ("Python", "Rust", None):
        group = [t for t in tests if name is None or suite(t) == name]
        res = [results.get(t["id"], ("not run",))[0] for t in group]
        noted = sum(bool((described.get(t["id"]) or {}).get("issues")) for t in group)
        lines.append(f"| {name or '**All**'} | {len(group)} | {res.count('passed')} | {len(res) - res.count('passed')} | "
                     + " | ".join(str(sum(strength(t) == k for t in group)) for k in kinds) + f" | {noted} |")
    lines += ["", "## Weak or unreviewed tests", ""]
    flagged = [t for t in tests if strength(t) in ("weak", "not reviewed")]
    lines += [f"- `{t['id']}` ({strength(t)}): " + cell("; ".join((described.get(t["id"]) or {}).get("issues") or ["no written description"]))
              for t in flagged] or ["None."]
    gaps = [t for t in tests if (described.get(t["id"]) or {}).get("verdict") == "gap"]
    lines += ["", f"## Gaps and proposed tests ({len(gaps)} of {len(tests)} tests)", "",
              "Tests that do not fully prove what they claim: what a bug could slip past, and a test that would catch it.", ""]
    for t in gaps:
        d = described[t["id"]]
        lines += [f"### `{t['id']}`", ""] + [f"- {cell(g)}" for g in d.get("gaps", [])] + [""]
        if d.get("proposal"):
            lines += [f"**Proposed test:** {cell(d['proposal'])}", ""]
    lines += ["", "## Every test, by file", ""]
    by_file = defaultdict(list)
    for t in tests:
        by_file[t["file"]].append(t)
    for file, group in by_file.items():
        lines += [f"### `{file}` ({len(group)} test{'s' if len(group) > 1 else ''})", "",
                  "| Test | Why it matters | Checked against | Strength | Proves its claim | Review notes |", "|---|---|---|---|---|---|"]
        for t in group:
            d = described.get(t["id"]) or {}
            status = results.get(t["id"], ("not run",))[0]
            name = f"`{t['name']}`" + ("" if status == "passed" else f" ({status})")
            lines.append(f"| {name} | {cell(d.get('role') or d.get('purpose') or t['doc'] or '')} | {cell(d.get('oracle', ''))} | "
                         f"{strength(t)} | {d.get('verdict') or 'not checked'} | {cell('<br>'.join(d.get('issues') or []))} |")
        lines.append("")
    (OUT / "AUDIT.md").write_text("\n".join(lines), encoding="utf-8")


def markdown_html(text):
    """The Markdown that AUDIT.md uses (headings, paragraphs, pipe tables, - and 1. lists, **bold**, `code`, [links](...)) as HTML."""
    def inline(t):
        t = t.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
        t = re.sub(r"`([^`]+)`", r"<code>\1</code>", t)
        t = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", t)
        t = re.sub(r"(?<![\w*])\*([^*\s][^*]*?)\*(?![\w*])", r"<i>\1</i>", t)
        t = re.sub(r"\[([^\]]+)\]\(([^)\s]+)\)", r'<a href="\2">\1</a>', t)
        return t.replace("\\|", "|").replace("&lt;br&gt;", "<br>")

    def cells(line):
        return [c.strip() for c in re.split(r"(?<!\\)\|", line.strip().strip("|"))]

    out, lines, i = [], text.splitlines(), 0
    while i < len(lines):
        line = lines[i]
        if not line.strip():
            i += 1
        elif m := re.match(r"(#{1,4}) (.*)", line):
            out.append(f"<h{len(m.group(1))}>{inline(m.group(2))}</h{len(m.group(1))}>")
            i += 1
        elif line.startswith("|") and i + 1 < len(lines) and re.match(r"\|[-|: ]+\|?$", lines[i + 1].strip()):
            head, i, body = cells(line), i + 2, []
            while i < len(lines) and lines[i].startswith("|"):
                body.append(cells(lines[i]))
                i += 1
            out.append("<table><thead><tr>" + "".join(f"<th>{inline(c)}</th>" for c in head) + "</tr></thead><tbody>"
                       + "".join("<tr>" + "".join(f"<td>{inline(c)}</td>" for c in r) + "</tr>" for r in body) + "</tbody></table>")
        elif re.match(r"\s*(- |\d+\. )", line):
            tag = "ol" if re.match(r"\s*\d+\. ", line) else "ul"
            items = []
            while i < len(lines) and (re.match(r"\s*(- |\d+\. )", lines[i]) or (lines[i].startswith("  ") and lines[i].strip())):
                m = re.match(r"(\s*)(?:- |\d+\. )(.*)", lines[i])
                if m and m.group(1):          # an indented item: a list inside the last item
                    items[-1][1].append(m.group(2))
                elif m:
                    items.append([m.group(2), []])
                else:                         # a continuation line of the last item
                    items[-1][0] += " " + lines[i].strip()
                i += 1
            out.append(f"<{tag}>" + "".join(f"<li>{inline(a)}" + ("<ul>" + "".join(f"<li>{inline(x)}</li>" for x in sub) + "</ul>" if sub else "")
                                            + "</li>" for a, sub in items) + f"</{tag}>")
        else:
            para = []
            while i < len(lines) and lines[i].strip() and not re.match(r"(#{1,4} |\||\s*- |\s*\d+\. )", lines[i]):
                para.append(lines[i].strip())
                i += 1
            out.append(f"<p>{inline(' '.join(para))}</p>")
    return "\n".join(out)


def write_html(tests, described, results, path, bare=False, repo="", ref="main"):
    """One self-contained page (search, filters, every test opening in place), from the same data as the Markdown pages."""
    rows = []
    for t in tests:
        d = described.get(t["id"]) or {}
        status, seconds, message = results.get(t["id"], ("not run", None, ""))
        rows.append({"id": t["id"], "file": t["file"], "name": t["name"], "line": t["line"], "suite": "py" if t["file"].endswith(".py") else "rs",
                     "result": status, "seconds": seconds, "message": message, "described": bool(d), "area": d.get("area", ""), "role": d.get("role", ""), "verdict": d.get("verdict", ""), "gaps": d.get("gaps", []), "proposal": d.get("proposal", ""), "evidence": d.get("evidence", ""),
                     "purpose": d.get("purpose", ""), "cases": d.get("cases", []), "oracle": d.get("oracle", ""), "strength": d.get("strength", ""),
                     "issues": d.get("issues", []), "fixed": d.get("fixed", []), "doc": t["doc"], "checks": t["checks"]})
    audit = OUT / "AUDIT.md"
    data = {"generated": datetime.date.today().isoformat(), "repo": repo, "ref": ref, "tests": rows,
            "audit": markdown_html(audit.read_text(encoding="utf-8")) if audit.exists() else ""}
    page = (Path(__file__).with_name("test_report_template.html").read_text(encoding="utf-8")
            .replace("/*DATA*/", json.dumps(data, ensure_ascii=False).replace("</", "<\\/")))
    if not bare:   # a page to open from disk; the published artifact gets the page alone (its host adds the document around it)
        page = '<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n<meta name="viewport" content="width=device-width, initial-scale=1">\n' + page.replace("<div class=\"wrap\">", "</head>\n<body>\n<div class=\"wrap\">", 1) + "\n</body>\n</html>\n"
    Path(path).write_text(page, encoding="utf-8")


def git_ref():
    try:
        ref = subprocess.run(["git", "rev-parse", "--abbrev-ref", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    except OSError:
        ref = ""
    return ref if ref and ref != "HEAD" else "main"


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--junit", help="pytest --junitxml file")
    parser.add_argument("--cargo", help="saved output of `cargo test`")
    parser.add_argument("--html", default=str(OUT / "report.html"), help="where to write the HTML page (default: docs/test-report/report.html)")
    parser.add_argument("--bare", action="store_true", help="HTML without <!doctype>/<head>/<body>, for hosts that add them")
    parser.add_argument("--repo", default="https://github.com/OpenDNA-project/openDNA", help="repository URL for the source links")
    args = parser.parse_args()
    described = {t["id"]: t for t in json.loads((OUT / "tests.json").read_text(encoding="utf-8"))["tests"]}
    results = {}
    if args.junit:
        results.update(junit_results(args.junit))
    if args.cargo:
        results.update(cargo_results(args.cargo))
    tests = python_tests() + rust_tests()
    known = {t["id"] for t in tests}
    stale = sorted(set(described) - known)
    totals = write(tests, described, results)
    write_audit(tests, described, results)
    write_html(tests, described, results, args.html, args.bare, args.repo, git_ref())
    print(f"{len(tests)} tests ({sum(t['file'].endswith('.py') for t in tests)} Python, {sum(t['file'].endswith('.rs') for t in tests)} Rust); "
          f"{len(known & set(described))} with a written description")
    print("results:", {k: v for k, v in totals.items() if not k.startswith("strength:")})
    print("strength:", {k[9:]: v for k, v in totals.items() if k.startswith("strength:")})
    if stale:
        print("descriptions of tests that no longer exist (rename or remove them in tests.json):", *stale, sep="\n  ")


if __name__ == "__main__":
    main()
