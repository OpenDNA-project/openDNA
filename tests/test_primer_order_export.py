"""Export for Ordering: the IDT plate file and the plain name/sequence file."""
import csv
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "app"))

import pytest

import primer_db


def rows(n, name="P{}"):
    return [{"name": name.format(i), "sequence": "ACGTACGTAC", "note": "", "added": ""} for i in range(n)]


def read(path):
    with open(path, encoding="utf-8", newline="") as fh:
        return list(csv.reader(fh))


def test_plate_wells_fill_down_the_columns_before_moving_right(tmp_path):
    [out] = primer_db.export_order(rows(10), tmp_path / "o.csv")
    table = read(out)
    assert table[0] == ["Well Position", "Name", "Sequence", "Scale", "Purification"]
    assert [r[0] for r in table[1:]] == ["A1", "B1", "C1", "D1", "E1", "F1", "G1", "H1", "A2", "B2"]


def test_more_than_96_oligos_continue_on_a_second_plate_file_starting_again_at_a1(tmp_path):
    written = primer_db.export_order(rows(100), tmp_path / "o.csv")
    assert [w.name for w in written] == ["o_plate1.csv", "o_plate2.csv"]
    first, second = read(written[0]), read(written[1])
    assert len(first) - 1 == 96 and first[-1][0] == "H12"
    assert len(second) - 1 == 4 and second[1][0] == "A1" and second[1][1] == "P96"


def test_scale_and_purification_are_written_on_every_row(tmp_path):
    [out] = primer_db.export_order(rows(2), tmp_path / "o.csv", scale="100nm", purification="HPLC")
    assert {tuple(r[3:]) for r in read(out)[1:]} == {("100nm", "HPLC")}


def test_a_name_with_a_comma_and_a_quote_survives_the_round_trip(tmp_path):
    odd = [{"name": 'fwd, "A" primer', "sequence": "ACGTACGT", "note": "", "added": ""}]
    [plate] = primer_db.export_order(odd, tmp_path / "p.csv")
    [plain] = primer_db.export_order(odd, tmp_path / "q.csv", primer_db.ORDER_FORMATS[1])
    assert read(plate)[1][1] == 'fwd, "A" primer'
    assert read(plain) == [['fwd, "A" primer', "ACGTACGT"]]


def test_a_tab_or_line_break_in_a_name_becomes_a_space(tmp_path):
    odd = [{"name": "a\tb\nc", "sequence": "ACGTACGT", "note": "", "added": ""}]
    [out] = primer_db.export_order(odd, tmp_path / "p.csv", primer_db.ORDER_FORMATS[1])
    assert read(out) == [["a b c", "ACGTACGT"]]


def test_degenerate_bases_are_exported_as_they_are(tmp_path):
    odd = [{"name": "deg", "sequence": "acgtNNRYacgt", "note": "", "added": ""}]
    [out] = primer_db.export_order(odd, tmp_path / "p.csv", primer_db.ORDER_FORMATS[1])
    assert read(out) == [["deg", "ACGTNNRYACGT"]]


def test_exporting_no_primers_is_refused_and_writes_nothing(tmp_path):
    with pytest.raises(ValueError):
        primer_db.export_order([], tmp_path / "o.csv")
    assert not os.listdir(tmp_path)
