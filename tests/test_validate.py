"""Tests for validate.py."""

import os

import pytest

import validate
from conftest import CANSNP_HEADER, CANSNP_ROWS, readme_text


def row(node, position, derived, ancestral, sequence="seq1", genome="G", rank=1,
        reference="Test2026", score="500", typed="yes"):
    """One canSNP row, with the columns of validate.CANSNP_HEADER."""
    return "\t".join([node, str(rank), reference, genome, sequence, str(position),
                      derived, ancestral, score, typed])



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
    rows = list(CANSNP_ROWS)
    rows[1] = row("A", "1²", "T", "C")
    problems = problems_of(make_scheme(cansnp_rows=rows))
    assert any("sequence_position '1²' is not a positive integer" in p for p in problems)


def test_positions_with_leading_zeros_count_as_the_same_position(make_scheme):
    rows = list(CANSNP_ROWS)
    rows[1] = row("A", "0100", "T", "C")
    problems = problems_of(make_scheme(cansnp_rows=rows))
    assert any("position 100 on 'G' sequence 'seq1' is annotated 2 times" in p for p in problems)


def test_old_header_names_the_conversion_command(make_scheme):
    old = ("snp_id\tstrain\treference\tgenome\tposition\tderived_base\tancestral_base\n")
    raw = (old + "R\ts\tTest2026\tG\t100\tA\tG\n").encode()
    problems = problems_of(make_scheme(raw_cansnps=raw))
    old_header = [p for p in problems if "CanSNPer2-database --convert" in p]
    assert len(old_header) == 1 and "_canSNPs.txt" in old_header[0]


def test_uniqueness_on_genome_sequence_position(make_scheme):
    # the same position on another sequence is another SNP
    rows = list(CANSNP_ROWS)
    rows[1] = row("A", 100, "T", "C", sequence="seq2")
    assert problems_of(make_scheme(cansnp_rows=rows)) == []
    rows[1] = row("A", 100, "T", "C")
    problems = problems_of(make_scheme(cansnp_rows=rows))
    assert any("position 100 on 'G' sequence 'seq1' is annotated 2 times" in p for p in problems)


def test_an_empty_sequence_name_is_a_problem(make_scheme):
    rows = list(CANSNP_ROWS)
    rows[1] = row("A", 200, "T", "C", sequence="")
    problems = problems_of(make_scheme(cansnp_rows=rows))
    assert any("sequence is empty" in p for p in problems)


def test_file_that_is_not_utf8_is_reported_not_raised(make_scheme):
    raw = (CANSNP_HEADER + "\n".join(CANSNP_ROWS) + "\n").encode("utf-8")
    raw = raw.replace(b"\tTest2026\t", b"\tTest\xe92026\t", 1)
    problems = problems_of(make_scheme(raw_cansnps=raw))
    assert any("not valid UTF-8" in p for p in problems)


def test_byte_order_mark_is_reported(make_scheme):
    raw = b"\xef\xbb\xbf" + (CANSNP_HEADER + "\n".join(CANSNP_ROWS) + "\n").encode()
    problems = problems_of(make_scheme(raw_cansnps=raw))
    assert any("byte order mark" in p for p in problems)


def test_node_with_two_parents_is_rejected(make_scheme):
    tree = "R\nR\tA\nR\tC\nR\tA\tB\nR\tC\tB\n"
    problems = problems_of(make_scheme(tree=tree))
    assert any("node 'B' is given more than one parent: A, C" in p for p in problems)


def test_known_issue_is_reported_as_a_warning(make_scheme):
    rows = list(CANSNP_ROWS)
    rows[1] = row("A", 200, "/", "C")
    problems, warnings = validate.validate(make_scheme(cansnp_rows=rows, known_issues=["A"]))
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


def test_a_genome_at_an_inner_node_is_a_problem(make_scheme):
    genomes = GENOMES_HEADER + "GCA_000000002.1\tA\trepresentative\n"
    problems = problems_of(make_scheme(genomes=genomes))
    assert any("node 'A' is not a terminal node" in problem for problem in problems), problems


def test_an_outgroup_is_not_held_to_the_terminal_node_rule(make_scheme):
    genomes = GENOMES_HEADER + "GCA_000000004.1\t\toutgroup\n"
    assert problems_of(make_scheme(genomes=genomes)) == []


def test_a_terminal_node_without_a_genome_is_valid(make_scheme):
    genomes = GENOMES_HEADER + "GCA_000000002.1\tB\trepresentative\n"
    assert problems_of(make_scheme(genomes=genomes)) == []


def test_one_assembly_under_two_accessions_is_a_problem(make_scheme):
    genomes = GENOMES_HEADER + "GCA_000000002.1\tB\tmember\nGCF_000000002.1\tC\tmember\n"
    problems = problems_of(make_scheme(genomes=genomes))
    assert any("one assembly is listed 2 times" in problem for problem in problems)


def test_a_genomes_file_without_its_header_is_a_problem(make_scheme):
    problems = problems_of(make_scheme(genomes="GCA_000000002.1\tB\tmember\n"))
    assert any("unexpected header" in problem for problem in problems)


SOURCE_HEADER = "genome\tnode\trole\tbasis\tsource\n"


def test_genomes_with_basis_and_source_are_valid(make_scheme):
    genomes = (SOURCE_HEADER + "GCA_000000002.1\tB\tmember\tpublished\tFoster2008\n"
               "GCF_000000003.2\tC\trepresentative\ttyped\tGTDB r232, CanSNPer2 2.3.0\n"
               "GCA_000000004.1\t\toutgroup\tpublished\tFoster2008\n")
    assert problems_of(make_scheme(genomes=genomes)) == []


@pytest.mark.parametrize("row,message", [
    ("GCA_000000002.1\tB\tmember\tverified\tFoster2008", "basis 'verified' is not one of"),
    ("GCA_000000002.1\tB\tmember\tpublished\t", "source is empty"),
    ("GCA_000000002.1\tB\tmember\tpublished", "expected 5 columns"),
])
def test_a_malformed_basis_or_source_is_a_problem(make_scheme, row, message):
    problems = problems_of(make_scheme(genomes=SOURCE_HEADER + row + "\n"))
    assert any(message in problem for problem in problems), problems


# -- the ranks, scores and typed marks of the canSNP file --------------------

SET_ROWS = [
    row("R", 100, "A", "G", rank=1, score="520", typed="yes"),
    row("R", 1500, "C", "T", rank=2, score="510", typed="no"),
    row("A", 200, "T", "C", rank=1, score="500", typed="yes"),
    row("B", 50, "C", "T", sequence="seq2", rank=1, score="495", typed="yes"),
    row("B", 2500, "G", "A", sequence="seq2", rank=2, score="480.5", typed="no"),
    row("C", 150, "G", "A", sequence="seq2", rank=1, score="470", typed="yes"),
]
SET_TEXT = "\n".join(SET_ROWS)
UP_TO_2 = readme_text(per_node="up to 2")


def test_a_scheme_with_one_canSNP_per_node_is_valid(make_scheme):
    assert validate.validate(make_scheme()) == ([], [])


def test_several_canSNPs_per_node_are_valid(make_scheme):
    scheme = make_scheme(cansnp_rows=SET_ROWS, readme=UP_TO_2)
    assert validate.validate(scheme) == ([], [])


def test_a_publication_canSNP_may_have_no_score_and_a_reference_of_NA(make_scheme):
    rows = list(SET_ROWS)
    rows[1] = row("R", 1500, "C", "T", rank=2, reference="NA", score="", typed="no")
    assert validate.validate(make_scheme(cansnp_rows=rows, readme=UP_TO_2)) == ([], [])


@pytest.mark.parametrize("old, new, message", [
    ("B\t2\t", "B\t3\t", "ranks of node 'B' are 1, 3"),
    ("B\t2\t", "B\t0\t", "rank '0' is not a positive integer"),
    ("B\t2\t", "B\t02\t", "rank '02' is not a positive integer"),
    ("B\t2\t", "B\tsecond\t", "rank 'second' is not a positive integer"),
    ("B\t2\t", "B\t1\t", "ranks of node 'B' are 1, 1"),
    ("B\t2\tTest2026\tG\tseq2\t2500", "B\t2\tTest2026\tG\tseq2\t50",
     "position 50 on 'G' sequence 'seq2' is annotated 2 times"),
    ("B\t2\tTest2026\tG\tseq2\t2500\tG\tA", "B\t2\tTest2026\tG\tseq2\t2500\tA\tA",
     "derived and ancestral base"),
    ("B\t2\tTest2026\tG\tseq2\t2500\tG\tA", "B\t2\tTest2026\tG\tseq2\t2500\tN\tA",
     "derived_base 'N' is not one of A, C, G, T"),
    ("R\t2\tTest2026\tG", "X\t2\tTest2026\tG", "node 'X' does not appear"),
    ("R\t2\tTest2026\tG", "R\t2\tTest2026\tH", "genome 'H' is not listed"),
    ("R\t2\tTest2026", "R\t2\t", "reference is empty"),
    ("\t480.5\t", "\thigh\t", "score 'high' is not a number"),
    ("\t480.5\tno", "\t480.5\tmaybe", "typed 'maybe' is not yes or no"),
    ("\t510\tno", "\t510\tyes", "node 'R' has 2 canSNPs marked typed yes"),
    ("\t520\tyes", "\t520\tno", "node 'R' has no canSNP marked typed yes"),
])
def test_a_canSNP_file_that_breaks_the_rules_is_a_problem(make_scheme, old, new, message):
    assert old in SET_TEXT, old
    rows = SET_TEXT.replace(old, new).split("\n")
    problems, _warnings = validate.validate(make_scheme(cansnp_rows=rows, readme=UP_TO_2))
    assert any(message in problem for problem in problems), problems


def test_ranks_out_of_order_are_a_problem(make_scheme):
    rows = list(SET_ROWS)
    rows[0], rows[1] = rows[1], rows[0]
    problems, _warnings = validate.validate(make_scheme(cansnp_rows=rows, readme=UP_TO_2))
    assert any("ranks of node 'R' are 2, 1, not 1, 2, 3 ... in order" in p for p in problems)


def test_a_node_without_a_canSNP_is_a_problem(make_scheme):
    rows = [r for r in CANSNP_ROWS if not r.startswith("C\t")]
    problems = problems_of(make_scheme(cansnp_rows=rows))
    assert any("node 'C' has no canSNP" in problem for problem in problems), problems


def test_a_canSNP_file_without_its_header_is_a_problem(make_scheme):
    raw = ("\n".join(SET_ROWS) + "\n").encode()
    assert problems_of(make_scheme(raw_cansnps=raw, readme=UP_TO_2))


def test_a_canSNP_file_with_a_missing_column_is_a_problem(make_scheme):
    rows = list(CANSNP_ROWS)
    rows[0] = "R\t1\tTest2026\tG\tseq1\t100\tA\tG\t520"
    problems = problems_of(make_scheme(cansnp_rows=rows))
    assert any("expected 10 columns, found 9" in problem for problem in problems), problems


# -- one folder per scheme under schemes/ ------------------------------------

def _scheme_folder(root, folder, name):
    d = root / folder
    d.mkdir(parents=True)
    for kind in ("tree", "canSNPs", "references"):
        (d / f"{name}_{kind}.txt").write_text("")


def test_discover_finds_a_scheme_by_its_folder(tmp_path, monkeypatch):
    schemes = tmp_path / "schemes"
    _scheme_folder(schemes, "abc", "abc")
    monkeypatch.setattr(validate, "REPO", str(tmp_path))
    monkeypatch.setattr(validate, "SCHEMES_DIR", str(schemes))
    monkeypatch.setattr(validate, "CONFIG", str(tmp_path / "none.toml"))
    found = validate.discover()
    assert [s.name for s in found] == ["abc"]
    assert found[0].cansnps == os.path.join("schemes", "abc", "abc_canSNPs.txt")


def test_discover_reports_a_canSNP_file_in_the_wrong_folder(tmp_path, monkeypatch, capsys):
    schemes = tmp_path / "schemes"
    _scheme_folder(schemes, "abc", "xyz")
    (schemes / "loose_canSNPs.txt").write_text("")
    monkeypatch.setattr(validate, "REPO", str(tmp_path))
    monkeypatch.setattr(validate, "SCHEMES_DIR", str(schemes))
    monkeypatch.setattr(validate, "CONFIG", str(tmp_path / "none.toml"))
    with pytest.raises(SystemExit) as stop:
        validate.discover()
    assert stop.value.code == 2
    err = capsys.readouterr().err
    assert "schemes/abc/xyz_canSNPs.txt" in err
    assert "schemes/loose_canSNPs.txt" in err


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


def test_canSNPs_per_node_follows_the_largest_rank(make_scheme):
    rows = list(CANSNP_ROWS)
    rows.insert(2, row("A", 250, "T", "C", rank=2, typed="no"))
    assert problems_of(make_scheme(cansnp_rows=rows, readme=readme_text(per_node="up to 2"))) == []
    problems = problems_of(make_scheme(cansnp_rows=rows))
    assert any("canSNPs per node is '1'" in p for p in problems)


def test_a_draft_is_rejected_only_when_verified_schemes_are_required(make_scheme):
    scheme = make_scheme(readme=readme_text(status="draft"))
    assert validate.validate(scheme)[0] == []
    problems, _ = validate.validate(scheme, require_verified=True)
    assert any("is a draft" in p for p in problems)


# -- the README of a taxonomy scheme ------------------------------------------

def _taxonomy_scheme(tmp_path, readme):
    (tmp_path / "tax_tree.txt").write_text("R\nR\tA\nR\tB\n")
    (tmp_path / "tax_genomes.txt").write_text("g1.fasta\tA\ng2.fasta\tB\n")
    if readme is not False:
        (tmp_path / "README.md").write_text(readme)
    return validate.Scheme("tax", kind="taxonomy", tree=str(tmp_path / "tax_tree.txt"),
                           genomes=str(tmp_path / "tax_genomes.txt"),
                           readme=str(tmp_path / "README.md"))


def test_a_taxonomy_scheme_needs_a_readme(tmp_path):
    problems, _ = validate.validate(_taxonomy_scheme(tmp_path, False))
    assert any("has no README.md" in p for p in problems)


def test_a_taxonomy_readme_names_no_reference_genome(tmp_path):
    readme = readme_text(name="tax", nodes="3", reference="none", per_node="none")
    assert validate.validate(_taxonomy_scheme(tmp_path, readme))[0] == []
    readme = readme_text(name="tax", nodes="3", reference="none", per_node="1")
    problems, _ = validate.validate(_taxonomy_scheme(tmp_path, readme))
    assert any("canSNPs per node is '1'" in p for p in problems)


def test_a_draft_taxonomy_scheme_is_rejected_when_verified_is_required(tmp_path):
    readme = readme_text(name="tax", status="draft", nodes="3", reference="none",
                         per_node="none")
    scheme = _taxonomy_scheme(tmp_path, readme)
    assert validate.validate(scheme)[0] == []
    problems, _ = validate.validate(scheme, require_verified=True)
    assert any("is a draft" in p for p in problems)


def test_a_taxonomy_scheme_defaults_to_its_folder():
    scheme = validate.Scheme("tax", kind="taxonomy")
    assert scheme.tree == os.path.join("schemes", "tax", "tax_tree.txt")
    assert scheme.genomes == os.path.join("schemes", "tax", "tax_genomes.txt")
    assert scheme.readme == os.path.join("schemes", "tax", "README.md")


# -- accessions of the reference genomes --------------------------------------

def test_a_genbank_id_repeated_on_two_rows_is_a_problem(make_scheme):
    references = ("genome\tstrain\tgenbank_id\trefseq_id\tassembly_name\n"
                  "G\tstrain G\tGCA_000000001.1\tGCF_000000001.1\tASM1v1\n"
                  "H\tstrain H\tGCA_000000001.1\tGCF_000000002.1\tASM2v1\n")
    scheme = make_scheme(references=references)
    problems = problems_of(scheme)
    relative = scheme.references
    assert (f"{relative}:3: genbank_id 'GCA_000000001.1' of genome 'H' repeats "
            f"the one of genome 'G' at {relative}:2") in problems


@pytest.mark.parametrize("accession", ["GCA_000000001.1", "GCF_000000001.1",
                                       "GCA_000000001.123"])
def test_a_genbank_id_follows_the_pattern_of_cansnper2(make_scheme, accession):
    references = ("genome\tstrain\tgenbank_id\trefseq_id\tassembly_name\n"
                  f"G\tstrain G\t{accession}\tGCF_000000001.1\tASM1v1\n")
    problems = problems_of(make_scheme(references=references))
    assert not any(f"genbank_id '{accession}'" in p for p in problems), problems


@pytest.mark.parametrize("accession", ["", "GCA_00000001.1", "GCA_000000001",
                                       "GCX_000000001.1", "GCA_000000001.\u0661"])
def test_a_malformed_genbank_id_is_a_problem(make_scheme, accession):
    references = ("genome\tstrain\tgenbank_id\trefseq_id\tassembly_name\n"
                  f"G\tstrain G\t{accession}\tGCF_000000001.1\tASM1v1\n")
    problems = problems_of(make_scheme(references=references))
    assert any(f"genbank_id '{accession}' is not of the form" in p
               for p in problems), problems


def test_the_genbank_id_pattern_is_the_one_cansnper2_uses():
    source_check = pytest.importorskip("CanSNPer2.modules.SourceCheck")
    assert validate.GENBANK_RE.pattern == source_check.ACCESSION_RE.pattern
