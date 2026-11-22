"""Tests for check_bases.py: SNP bases against the reference genomes.

The genome files are placed in the cache folder, so no test downloads anything."""

import gzip

import pytest

import check_bases
from conftest import CANSNP_ROWS

# Two sequences of 250 bases. The SNPs of the test scheme are at seq1:100, seq1:200,
# seq2:50 and seq2:150, counted along each sequence. The keys of BASES are the
# positions along the two sequences joined in file order (seq2 starts at 251).
BASES = {100: "G", 200: "C", 300: "T", 400: "A"}   ## the ancestral bases of CANSNP_ROWS


def genome(swap=False, change=None):
    joined = ["N"] * 500   ## a filler that matches no SNP base where the sequences are swapped
    for pos, base in {**BASES, **(change or {})}.items():
        joined[pos - 1] = base
    seq1, seq2 = "".join(joined[:250]), "".join(joined[250:])
    records = [("seq2", seq2), ("seq1", seq1)] if swap else [("seq1", seq1), ("seq2", seq2)]
    return "".join(f">{name} test\n{seq}\n" for name, seq in records)


@pytest.fixture
def cache(tmp_path):
    folder = tmp_path / "cache"
    folder.mkdir()

    def put(text):
        path = folder / "GCA_000000001.1_ASM1v1_genomic.fna.gz"
        path.write_bytes(gzip.compress(text.encode()))
        return str(folder)
    return put


def never_download(url, target):
    raise AssertionError(f"no download expected, asked for {url}")


def test_a_scheme_whose_reference_has_the_snp_bases_passes(make_scheme, cache):
    problems, warnings, report = check_bases.check_scheme(make_scheme(), cache(genome()), retrieve=never_download)
    assert problems == [] and warnings == []
    assert "G: seq1 250 bp (2 SNPs), seq2 250 bp (2 SNPs)" in report


def test_a_reference_with_its_sequences_in_another_order_passes(make_scheme, cache):
    problems, _warnings, report = check_bases.check_scheme(make_scheme(), cache(genome(swap=True)), retrieve=never_download)
    assert problems == []
    assert "G: seq2 250 bp (2 SNPs), seq1 250 bp (2 SNPs)" in report


def test_a_missing_sequence_is_a_problem(make_scheme, cache):
    rows = list(CANSNP_ROWS)
    rows[1] = "A\t1\tTest2026\tG\tplasmid\t200\tT\tC\t500\tyes"
    problems, _warnings, _report = check_bases.check_scheme(make_scheme(cansnp_rows=rows), cache(genome()), retrieve=never_download)
    assert len(problems) == 1
    assert "sequence 'plasmid'" in problems[0] and "seq1, seq2" in problems[0]


def test_a_single_wrong_base_fails(make_scheme, cache):
    problems, _warnings, _report = check_bases.check_scheme(make_scheme(), cache(genome(change={200: "A"})), retrieve=never_download)
    assert len(problems) == 1
    assert "_canSNPs.txt:3: canSNP A at seq1 position 200 of G" in problems[0]
    assert "has A" in problems[0]


def test_a_known_issue_is_a_warning(make_scheme, cache):
    scheme = make_scheme(known_issues=["A"])
    problems, warnings, _report = check_bases.check_scheme(scheme, cache(genome(change={200: "A"})), retrieve=never_download)
    assert problems == []
    assert len(warnings) == 1 and "canSNP A" in warnings[0]


def test_a_position_beyond_the_reference_fails(make_scheme, cache):
    rows = list(CANSNP_ROWS)
    rows[3] = "C\t1\tTest2026\tG\tseq2\t900\tG\tA\t470\tyes"
    problems, _warnings, _report = check_bases.check_scheme(make_scheme(cansnp_rows=rows), cache(genome()), retrieve=never_download)
    assert any("position 900" in p and "seq2" in p and "250 bases" in p for p in problems)


def test_the_ranked_canSNPs_are_checked_as_well(make_scheme, cache):
    rows = list(CANSNP_ROWS)
    rows.insert(2, "A\t2\tTest2026\tG\tseq1\t220\tC\tG\t490\tno")
    problems, _warnings, _report = check_bases.check_scheme(make_scheme(cansnp_rows=rows), cache(genome()), retrieve=never_download)
    assert len(problems) == 1
    assert "_canSNPs.txt:4: canSNP A (rank 2) at seq1 position 220 of G" in problems[0]


def test_a_missing_genome_is_downloaded_once(make_scheme, tmp_path):
    folder = tmp_path / "empty"
    folder.mkdir()
    asked = []

    def retrieve(url, target):
        asked.append(url)
        with open(target, "wb") as out:
            out.write(gzip.compress(genome().encode()))
    problems, _warnings, _report = check_bases.check_scheme(make_scheme(), str(folder), retrieve=retrieve)
    assert problems == []
    assert asked == ["https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/000/001/"
                     "GCA_000000001.1_ASM1v1/GCA_000000001.1_ASM1v1_genomic.fna.gz"]
    check_bases.check_scheme(make_scheme(), str(folder), retrieve=never_download)


def test_a_failed_download_is_a_problem_not_a_crash(make_scheme, tmp_path):
    folder = tmp_path / "empty"
    folder.mkdir()

    def retrieve(url, target):
        raise OSError("no network")
    problems, _warnings, _report = check_bases.check_scheme(make_scheme(), str(folder), retrieve=retrieve)
    assert len(problems) == 1 and "could not get GCA_000000001.1" in problems[0] and "no network" in problems[0]
    assert list(folder.iterdir()) == []      ## nothing half written is left to be used later
