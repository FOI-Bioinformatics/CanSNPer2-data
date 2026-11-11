"""Tests for validate.py."""

import pytest

import validate
from conftest import SNP_HEADER, SNP_ROWS


def problems_of(scheme):
    problems, _warnings = validate.validate(scheme)
    return problems


def test_minimal_scheme_is_valid(make_scheme):
    assert validate.validate(make_scheme()) == ([], [])


def test_repository_schemes_are_valid():
    for scheme in validate.discover():
        problems, _warnings = validate.validate(scheme)
        assert problems == [], (scheme.name, problems)


def test_non_ascii_digit_position_is_reported_not_raised(make_scheme):
    # "1²" passes str.isdigit() but int() rejects it; this used to end the run
    rows = list(SNP_ROWS)
    rows[1] = "A\ts\tTest2026\tG\t1²\tT\tC"
    problems = problems_of(make_scheme(snp_rows=rows))
    assert any("position '1²' is not a positive integer" in p for p in problems)


def test_positions_with_leading_zeros_count_as_the_same_position(make_scheme):
    rows = list(SNP_ROWS)
    rows[1] = "A\ts\tTest2026\tG\t0100\tT\tC"
    problems = problems_of(make_scheme(snp_rows=rows))
    assert any("position 100 on 'G' is annotated 2 times" in p for p in problems)


def test_file_that_is_not_utf8_is_reported_not_raised(make_scheme):
    raw = (SNP_HEADER + "\n".join(SNP_ROWS) + "\n").encode("utf-8")
    raw = raw.replace(b"\ts\t", b"\ts\xe9\t", 1)
    problems = problems_of(make_scheme(raw_snp=raw))
    assert any("not valid UTF-8" in p for p in problems)


def test_byte_order_mark_is_reported(make_scheme):
    raw = b"\xef\xbb\xbf" + (SNP_HEADER + "\n".join(SNP_ROWS) + "\n").encode()
    problems = problems_of(make_scheme(raw_snp=raw))
    assert any("byte order mark" in p for p in problems)


def test_node_with_two_parents_is_rejected(make_scheme):
    tree = "R\nR\tA\nR\tC\nR\tA\tB\nR\tC\tB\n"
    problems = problems_of(make_scheme(tree=tree))
    assert any("node 'B' is given more than one parent: A, C" in p for p in problems)


def test_known_issue_is_reported_as_a_warning(make_scheme):
    rows = list(SNP_ROWS)
    rows[1] = "A\ts\tTest2026\tG\t200\t/\tC"
    problems, warnings = validate.validate(make_scheme(snp_rows=rows, known_issues=["A"]))
    assert problems == []
    assert any("derived_base '/'" in w for w in warnings)


def test_known_issue_naming_no_snp_is_a_problem(make_scheme):
    problems = problems_of(make_scheme(known_issues=["NOT.A.SNP"]))
    assert any("'NOT.A.SNP', which is not in" in p for p in problems)


def test_known_issue_for_a_snp_with_valid_bases_is_a_problem(make_scheme):
    problems = problems_of(make_scheme(known_issues=["A"]))
    assert any("'A', whose bases are valid" in p for p in problems)


def test_unknown_config_key_is_reported():
    problems = validate.check_config({"x": {"known_issue": ["A"]}, "y": {"wip": True}})
    assert len(problems) == 1
    assert "[x] has an unknown key 'known_issue'" in problems[0]


def test_config_entry_that_is_not_a_table_is_reported():
    assert "must be a table" in validate.check_config({"x": 1})[0]


def test_discover_stops_with_a_message_on_an_unknown_key(tmp_path, monkeypatch, capsys):
    source = tmp_path / "source"
    source.mkdir()
    config = tmp_path / "schemes.toml"
    config.write_text('[x]\nkind = "snp"\nknwon_issues = []\n')
    monkeypatch.setattr(validate, "SOURCE_DIR", str(source))
    monkeypatch.setattr(validate, "CONFIG", str(config))
    with pytest.raises(SystemExit) as stop:
        validate.discover()
    assert stop.value.code == 2
    assert "unknown key 'knwon_issues'" in capsys.readouterr().err
