"""Tests for validate.py."""

import os

import pytest

import validate
from conftest import SNP_HEADER, SNP_ROWS, readme_text


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
    source = tmp_path / "schemes"
    source.mkdir()
    config = tmp_path / "schemes.toml"
    config.write_text('[x]\nkind = "snp"\nknwon_issues = []\n')
    monkeypatch.setattr(validate, "SCHEMES_DIR", str(source))
    monkeypatch.setattr(validate, "CONFIG", str(config))
    with pytest.raises(SystemExit) as stop:
        validate.discover()
    assert stop.value.code == 2
    assert "unknown key 'knwon_issues'" in capsys.readouterr().err


# -- the genomes an SNP scheme was made from ---------------------------------

GENOMES_HEADER = "genome\tnode\trole\n"


def test_a_scheme_without_a_genomes_file_is_valid(make_scheme):
    assert problems_of(make_scheme()) == []


def test_a_genomes_file_is_checked_when_present(make_scheme):
    genomes = (GENOMES_HEADER + "GCA_000000002.1\tB\trepresentative\n"
               "GCF_000000003.2\tB\tmember\nGCA_000000004.1\t\toutgroup\n")
    assert problems_of(make_scheme(genomes=genomes)) == []


@pytest.mark.parametrize("row,message", [
    ("strain-7\tB\trepresentative", "not a versioned assembly accession"),
    ("GCA_000000002\tB\trepresentative", "not a versioned assembly accession"),
    ("GCA_000000002.1\tZ\trepresentative", "node 'Z' does not appear"),
    ("GCA_000000002.1\tB\ttype strain", "role 'type strain' is not one of"),
    ("GCA_000000002.1\tB\toutgroup", "outgroup sits at no node"),
    ("GCA_000000002.1\tB", "expected 3 columns"),
])
def test_a_malformed_genomes_row_is_a_problem(make_scheme, row, message):
    problems = problems_of(make_scheme(genomes=GENOMES_HEADER + row + "\n"))
    assert any(message in problem for problem in problems), problems


def test_one_assembly_under_two_accessions_is_a_problem(make_scheme):
    genomes = GENOMES_HEADER + "GCA_000000002.1\tB\tmember\nGCF_000000002.1\tA\tmember\n"
    problems = problems_of(make_scheme(genomes=genomes))
    assert any("one assembly is listed 2 times" in problem for problem in problems)


def test_a_genomes_file_without_its_header_is_a_problem(make_scheme):
    problems = problems_of(make_scheme(genomes="GCA_000000002.1\tB\tmember\n"))
    assert any("unexpected header" in problem for problem in problems)


# -- the optional ranked canSNP set ------------------------------------------

SNP_SET_HEADER = "node\trank\tgenome\tposition\tderived_base\tancestral_base\tscore\n"
SNP_SET = SNP_SET_HEADER + (
    "R\t1\tG\t100\tA\tG\t520\n"
    "R\t2\tG\t1500\tC\tT\t510\n"
    "A\t1\tG\t200\tT\tC\t500\n"
    "B\t1\tG\t300\tC\tT\t495\n"
    "B\t2\tG\t2500\tG\tA\t480.5\n"
    "C\t1\tG\t400\tG\tA\t470\n"
)


def test_a_scheme_without_a_snp_set_is_valid(make_scheme):
    assert validate.validate(make_scheme()) == ([], [])


def test_a_snp_set_is_checked_when_present(make_scheme):
    scheme = make_scheme(snp_set=SNP_SET, readme=readme_text(per_node="up to 2"))
    assert validate.validate(scheme) == ([], [])


@pytest.mark.parametrize("old, new, message", [
    ("R\t1\tG\t100\tA\tG", "R\t1\tG\t101\tA\tG", "rank 1 of node 'R' is not the SNP"),
    ("B\t2\t", "B\t3\t", "ranks of node 'B' are 1, 3"),
    ("B\t2\tG\t2500", "B\t2\tG\t1500", "position 1500 on 'G' is listed 2 times"),
    ("B\t2\tG\t2500\tG\tA", "B\t2\tG\t2500\tA\tA", "derived and ancestral base"),
    ("B\t2\tG\t2500\tG\tA", "B\t2\tG\t2500\tN\tA", "not both one of A, C, G, T"),
    ("R\t2\tG", "X\t2\tG", "node 'X' does not appear"),
    ("R\t2\tG", "R\t2\tH", "genome 'H' is not listed"),
    ("\t480.5", "\thigh", "score 'high' is not a number"),
    ("C\t1\tG\t400\tG\tA\t470\n", "", "node 'C' has no canSNP of rank 1"),
])
def test_a_snp_set_that_breaks_the_rules_is_a_problem(make_scheme, old, new, message):
    assert old in SNP_SET
    problems, _warnings = validate.validate(make_scheme(snp_set=SNP_SET.replace(old, new)))
    assert any(message in problem for problem in problems), problems


def test_a_snp_set_without_its_header_is_a_problem(make_scheme):
    problems, _warnings = validate.validate(make_scheme(snp_set=SNP_SET[len(SNP_SET_HEADER):]))
    assert problems


# -- one folder per scheme under schemes/ ------------------------------------

def _scheme_folder(root, folder, name):
    d = root / folder
    d.mkdir(parents=True)
    for kind in ("tree", "snp", "references"):
        (d / f"{name}_{kind}.txt").write_text("")


def test_discover_finds_a_scheme_by_its_folder(tmp_path, monkeypatch):
    schemes = tmp_path / "schemes"
    _scheme_folder(schemes, "abc", "abc")
    monkeypatch.setattr(validate, "REPO", str(tmp_path))
    monkeypatch.setattr(validate, "SCHEMES_DIR", str(schemes))
    monkeypatch.setattr(validate, "CONFIG", str(tmp_path / "none.toml"))
    found = validate.discover()
    assert [s.name for s in found] == ["abc"]
    assert found[0].snp == os.path.join("schemes", "abc", "abc_snp.txt")


def test_discover_reports_a_snp_file_in_the_wrong_folder(tmp_path, monkeypatch, capsys):
    schemes = tmp_path / "schemes"
    _scheme_folder(schemes, "abc", "xyz")
    (schemes / "loose_snp.txt").write_text("")
    monkeypatch.setattr(validate, "REPO", str(tmp_path))
    monkeypatch.setattr(validate, "SCHEMES_DIR", str(schemes))
    monkeypatch.setattr(validate, "CONFIG", str(tmp_path / "none.toml"))
    with pytest.raises(SystemExit) as stop:
        validate.discover()
    assert stop.value.code == 2
    err = capsys.readouterr().err
    assert "schemes/abc/xyz_snp.txt" in err
    assert "schemes/loose_snp.txt" in err


# -- the README of each scheme ------------------------------------------------


def test_a_conforming_readme_is_valid(make_scheme):
    assert problems_of(make_scheme()) == []


def test_a_missing_readme_is_reported(make_scheme):
    assert any("has no README.md" in p for p in problems_of(make_scheme(readme=False)))


def test_facts_must_match_the_scheme_files(make_scheme):
    readme = readme_text(name="other", nodes="5", reference="strain X (GCA_000000009.1)",
                         per_node="up to 3", changed="19 Nov 2026")
    problems = "\n".join(problems_of(make_scheme(readme=readme)))
    assert "Scheme is `other`" in problems
    assert "Nodes is 5" in problems
    assert "Reference genome names no accession" in problems
    assert "canSNPs per node is 'up to 3'" in problems
    assert "Last changed" in problems


def test_status_must_be_draft_or_verified(make_scheme):
    problems = problems_of(make_scheme(readme=readme_text(status="final")))
    assert any("Status" in p and "draft or verified" in p for p in problems)


def test_sections_must_all_be_present_in_order(make_scheme):
    sections = ["Description", "Validation", "Origin", "Known issues", "References"]
    problems = "\n".join(problems_of(make_scheme(readme=readme_text(sections=sections))))
    assert "sections must be" in problems


def test_an_empty_section_is_reported(make_scheme):
    readme = readme_text().replace("Text about known issues.", "")
    assert any("'Known issues' is empty" in p for p in problems_of(make_scheme(readme=readme)))


def test_canSNPs_per_node_follows_the_ranked_file(make_scheme):
    snp_set = ("node\trank\tgenome\tposition\tderived_base\tancestral_base\tscore\n"
               "R\t1\tG\t100\tA\tG\t1\nA\t1\tG\t200\tT\tC\t1\nA\t2\tG\t250\tT\tC\t1\n"
               "B\t1\tG\t300\tC\tT\t1\nC\t1\tG\t400\tG\tA\t1\n")
    assert problems_of(make_scheme(snp_set=snp_set, readme=readme_text(per_node="up to 2"))) == []
    problems = problems_of(make_scheme(snp_set=snp_set))
    assert any("canSNPs per node is '1'" in p for p in problems)


def test_a_draft_is_rejected_only_when_verified_schemes_are_required(make_scheme):
    scheme = make_scheme(readme=readme_text(status="draft"))
    assert validate.validate(scheme)[0] == []
    problems, _ = validate.validate(scheme, require_verified=True)
    assert any("is a draft" in p for p in problems)
