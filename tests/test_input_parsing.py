"""What the user types or imports (no window involved): the range box, pasted protein text, names and the NCBI key, primers and
primer lists. Each bad input is refused with a reason and never reaches the database."""
import os
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ["OPENDNA_DATA_DIR"] = tempfile.mkdtemp(prefix="opendna_test_")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "app"))

import pytest

import validate
from editor import clean_pasted_protein
from main import parse_range
from primer_db import db, primer_problem
from primer_panel import import_primers, parse_primer_list


def test_typed_or_imported_text_is_parsed_or_refused_with_a_reason(tmp_path):
    """The range box reads the usual ways of writing a stretch (and refuses the rest); pasted protein text is cleaned or refused;
    names and the NCBI key are checked; a primer and a primer list are checked line by line, every refused line named with its
    number and reason, and nothing invalid reaches the primer database whichever way it comes (a list, add, add_many, update)."""
    # --- the range box
    assert parse_range("101..200", 3780, True) == (100, 200)
    assert parse_range(" 101-200 ", 3780, True) == parse_range("101 to 200", 3780, True) == parse_range("200..101", 3780, False)
    assert parse_range("1,000..2,000", 3780, True) == (999, 2000)
    assert parse_range("3700..50", 3780, True) == (3699, 3830), "through the origin of a circular sequence"
    assert parse_range("3700..50", 3780, False) == (49, 3700), "a linear one: the two ends in order"
    assert parse_range("42", 3780, True) == (42, 42), "one number: the cursor after that base"
    for bad in ("", "abc", "1..2..3", "0..10", "100..4000"):
        with pytest.raises(ValueError):
            parse_range(bad, 3780, True)

    # --- pasted protein text: header, numbers and spaces dropped, anything else refused
    assert clean_pasted_protein(">p\nMKT AYI\n10 AKQ") == "MKTAYIAKQ"
    with pytest.raises(ValueError):
        clean_pasted_protein("MKT123?!")

    # --- names and the NCBI key
    assert validate.text_problem("fine name", "name", 20) == ""
    assert "control character" in validate.text_problem("a\tb", "name", 20) and "too long" in validate.text_problem("x" * 30, "name", 20)
    assert validate.api_key_problem("0123abcdef") == "" and "letters and digits" in validate.api_key_problem("abc-def!")

    # --- one primer: each problem named
    for name, seq, expect in [
        ("", "ACGTACGTAC", "name"),
        ("a\tb", "ACGTACGTAC", "control character"),
        ("ok_short", "ACGTAC", "at least 8"),
        ("ok_long", "A" * 201, "at most"),
        ("ok_bad", "ACGTACGTXX", "not DNA"),
    ]:
        assert expect in primer_problem(name, seq), (name, expect)

    # --- a bogus primer list: only real primers get in, every refused line named
    bogus = tmp_path / "bogus.tsv"
    bogus.write_text("name\tsequence\n"
                     "good1\tACGTACGTACGTACGTACGT\n"
                     "words\thello world, this is not dna\n"
                     "short\tACGT\n"
                     "just one column with no separator\n"
                     "\tACGTACGTACGTACGT\n"
                     "uracil\tACGUACGUACGUACGU\n"
                     "\n# a comment\n"
                     "good2,acgtacgtacgtACGTAC\n", encoding="utf-8")
    rejected = []
    report = import_primers(str(bogus), None, rejected)
    assert db().find("good1") and db().find("good2")
    assert [n for n, _, _ in rejected] == [3, 4, 5, 6, 7], rejected
    reasons = dict((n, why) for n, _, why in rejected)
    assert "not DNA" in reasons[3] and "at least 8" in reasons[4] and "No tab, comma or semicolon" in reasons[5]
    assert "name" in reasons[6].lower() and "U" in reasons[7]
    assert "5 lines not imported" in report and "line 3" in report
    for name in ("words", "short", "uracil"):
        assert db().find(name) is None

    # --- a FASTA list, and a binary file
    rejected = []
    out = parse_primer_list(">fa1 desc\nACGTACGT\nACGTAC\n>fa2\nACGT-XYZ-ACGT\n>\nACGTACGTACGT\n", rejected)
    assert out == [("fa1", "ACGTACGTACGTAC")]
    assert [n for n, _, _ in rejected] == [4, 6] and "not DNA" in rejected[0][2]
    binary = tmp_path / "x.tsv"
    binary.write_bytes(b"name\tsequence\nab\x00cd\tACGTACGTACGT\n")
    with pytest.raises(ValueError, match="not a text file"):
        import_primers(str(binary))

    # --- the database itself refuses what the list would refuse
    assert db().add("tab\there", "ACGTACGTACGT") == "invalid" and db().add("ok_name", "ACGTXX") == "invalid"
    report = db().add_many([("fine_one", "ACGTACGTACGT"), ("bad\tone", "ACGTACGTACGT"), ("toolong", "A" * 500)])
    assert report["added"] == 1 and [n for n, _ in report["invalid"]] == ["bad\tone", "toolong"]
    with pytest.raises(ValueError):
        db().update("fine_one", "fine_one", "NOT DNA AT ALL", "")
