"""Tests for build.py. They build real databases, so they need CanSNPer2 2.3.0 or later."""

import importlib.metadata
import json
import os
import shutil
import sqlite3

import pytest

import build
from conftest import CANSNP_HEADER, CANSNP_ROWS, readme_text

try:
    build.require_cansnper2()
except RuntimeError as error:
    pytest.skip(f"needs CanSNPer2: {error}", allow_module_level=True)

from CanSNPer2.modules.DatabaseConnection import SourceFileError  # noqa: E402


@pytest.fixture
def built(make_scheme, tmp_path):
    """Build the conftest scheme and return (scheme, path to the database)."""
    def make(**kwargs):
        scheme = make_scheme(**kwargs)
        target = str(tmp_path / f"built_{len(os.listdir(tmp_path))}.db")
        warnings = build.build(scheme, target)
        return scheme, target, warnings
    return make


# -- building and structure -------------------------------------------------

def test_minimal_scheme_builds_cleanly(built):
    scheme, target, warnings = built()
    assert warnings == []
    assert build.check_structure(target) == []
    assert build.check_counts(scheme, target) == []


def test_a_snp_on_an_unknown_genome_is_not_built_silently(built):
    # CanSNPer2 with source file checks (CanSNPer2#66) refuses the files; an
    # earlier CanSNPer2 skips the row with a logged warning, which build()
    # returns and check_counts sees as a gap
    rows = list(CANSNP_ROWS)
    rows[3] = "C\t1\tTest2026\tOTHER\tseq2\t150\tG\tA\t470\tyes"
    try:
        scheme, target, warnings = built(cansnp_rows=rows)
    except SourceFileError as error:
        assert "genome 'OTHER' is not listed" in str(error)
        return
    assert warnings, "the skipped row must be reported"
    assert build.check_counts(scheme, target) == [
        "snp_annotation has 3 rows, the sources define 4"]


def test_check_counts_sees_a_missing_row(built):
    scheme, target, _warnings = built()
    assert build.check_counts(scheme, target) == []
    connection = sqlite3.connect(target)
    connection.execute("DELETE FROM snp_annotation WHERE snp_id = 'C'")
    connection.commit()
    connection.close()
    assert build.check_counts(scheme, target) == [
        "snp_annotation has 3 rows, the sources define 4"]


def test_check_counts_counts_every_canSNP_of_a_node(built):
    rows = list(CANSNP_ROWS)
    rows.insert(2, "A\t2\tTest2026\tG\tseq1\t220\tC\tG\t490\tno")
    scheme, target, warnings = built(cansnp_rows=rows, readme=readme_text(per_node="up to 2"))
    assert warnings == []
    assert build.check_counts(scheme, target) == []
    connection = sqlite3.connect(target)
    assert connection.execute("SELECT COUNT(*) FROM snp_annotation").fetchone()[0] == 5
    connection.close()


def test_check_counts_counts_every_canSNP_of_a_node(built):
    rows = list(CANSNP_ROWS)
    rows.insert(2, "A\t2\tTest2026\tG\tseq1\t220\tC\tG\t490\tno")
    scheme, target, warnings = built(cansnp_rows=rows, readme=readme_text(per_node="up to 2"))
    assert warnings == []
    assert build.check_counts(scheme, target) == []
    connection = sqlite3.connect(target)
    assert connection.execute("SELECT COUNT(*) FROM snp_annotation").fetchone()[0] == 5
    connection.close()


def test_node_with_two_parents_is_not_built_silently(built):
    # validate.py rejects this source; CanSNPer2 with source file checks
    # refuses it too, and with an earlier CanSNPer2 the database check is the
    # second line
    tree = "R\nR\tA\nR\tC\nR\tA\tB\nR\tC\tB\n"
    try:
        _scheme, target, _warnings = built(tree=tree)
    except SourceFileError as error:
        assert "node 'B' is given more than one parent" in str(error)
        return
    assert "node name 'B' is used by more than one node" in build.check_structure(target)


def test_more_than_one_root_is_rejected_for_a_snp_scheme(built):
    _scheme, target, _warnings = built()
    connection = sqlite3.connect(target)
    connection.execute("INSERT INTO nodes (id, name) VALUES (99, 'root')")
    connection.execute("INSERT INTO tree (parent, child, rank_i) VALUES (99, 99, 1)")
    connection.commit()
    connection.close()
    problems = build.check_structure(target)
    assert any(p.startswith("2 roots where one is expected") for p in problems)
    # the flextaxd taxonomy layout has an orphan root beside the scheme root
    assert build.check_structure(target, single_root=False) == []


# -- comparison --------------------------------------------------------------

def test_comparison_sees_a_snp_on_the_wrong_node(built, tmp_path):
    _scheme, target, _warnings = built()
    swapped = str(tmp_path / "swapped.db")
    shutil.copy(target, swapped)
    connection = sqlite3.connect(swapped)
    ids = dict(connection.execute("SELECT name, id FROM nodes"))
    connection.executescript(f"""
        UPDATE snp_annotation SET node_id = -1 WHERE node_id = {ids['A']};
        UPDATE snp_annotation SET node_id = {ids['A']} WHERE node_id = {ids['B']};
        UPDATE snp_annotation SET node_id = {ids['B']} WHERE node_id = -1;""")
    connection.commit()
    connection.close()
    differences = build.compare(build.content(target), build.content(swapped))
    assert [table for table, _, _ in differences] == ["snp_annotation"]


def test_comparison_sees_a_repeated_node_name(built, tmp_path):
    _scheme, target, _warnings = built()
    repeated = str(tmp_path / "repeated.db")
    shutil.copy(target, repeated)
    connection = sqlite3.connect(repeated)
    connection.execute("INSERT INTO nodes (name) VALUES ('B')")
    connection.commit()
    connection.close()
    differences = dict((table, right) for table, _, right in
                       build.compare(build.content(target), build.content(repeated)))
    assert differences == {"duplicated node names": ["B"]}


def test_comparison_reads_a_database_with_joined_positions(built, tmp_path):
    # a committed database built before CanSNPer2 2.3.0 has a single position
    # column; it must be compared, and differ, rather than stop the build
    _scheme, target, _warnings = built()
    joined = str(tmp_path / "joined.db")
    shutil.copy(target, joined)
    connection = sqlite3.connect(joined)
    # the table is recreated, since a UNIQUE constraint names the sequence column
    connection.executescript("""
        CREATE TABLE snp_annotation_joined (
            node_id INTEGER, snp_id TEXT, position INTEGER, ancestral_base TEXT,
            derived_base TEXT, reference TEXT, date TEXT, genome_i INTEGER);
        INSERT INTO snp_annotation_joined
            SELECT node_id, snp_id, sequence_position, ancestral_base, derived_base,
                   reference, date, genome_i FROM snp_annotation;
        DROP TABLE snp_annotation;
        ALTER TABLE snp_annotation_joined RENAME TO snp_annotation;""")
    connection.commit()
    connection.close()
    differences = build.compare(build.content(joined), build.content(target))
    assert [table for table, _, _ in differences] == ["snp_annotation"]
    _table, only_joined, only_built = differences[0]
    assert len(only_joined) == len(only_built) == 4
    assert all("None" in row for row in only_joined)


# -- scheme_info -------------------------------------------------------------

def scheme_info_of(path):
    from CanSNPer2.modules.SchemeInfo import read_scheme_info
    return read_scheme_info(path)


@pytest.fixture
def clean_git(monkeypatch):
    """Give every scheme a committed source state without a git checkout."""
    monkeypatch.setattr(build, "source_commit", lambda scheme: "a" * 40)
    monkeypatch.setattr(build, "dirty_sources", lambda scheme: [])


def test_the_database_records_the_scheme(built, clean_git):
    _scheme, target, _warnings = built()
    info = scheme_info_of(target)
    assert info["scheme"] == "test"
    assert info["status"] == "verified"
    assert info["source"] == "FOI-Bioinformatics/CanSNPer2-data"
    # the branch of a checkout goes stale once it is merged; commit names the sources
    assert "ref" not in info
    assert info["commit"] == "a" * 40
    assert info["format"] == "2.4"
    assert "known_issues" not in info or info["known_issues"] == ""


def test_the_status_is_read_from_the_readme(built, clean_git):
    _scheme, target, _warnings = built(readme=readme_text(status="draft"))
    assert scheme_info_of(target)["status"] == "draft"


def test_known_issues_are_read_from_schemes_toml(make_scheme, clean_git, tmp_path, monkeypatch):
    from CanSNPer2.modules.FetchScheme import LocalSchemeSource
    config = tmp_path / "toml"
    config.mkdir()
    (config / "schemes.toml").write_text('[test]\nknown_issues = ["A", "C"]\n')
    monkeypatch.setattr(build, "scheme_source", lambda: LocalSchemeSource(str(config)))
    assert build.scheme_info(make_scheme())["known_issues"] == "A,C"


def test_uncommitted_sources_are_marked_in_the_commit(make_scheme, clean_git, monkeypatch):
    monkeypatch.setattr(build, "dirty_sources", lambda scheme: ["schemes/test/test_canSNPs.txt"])
    assert build.scheme_info(make_scheme())["commit"] == "a" * 40 + "+uncommitted"


def test_the_comparison_ignores_the_build_date_but_not_the_status(built, clean_git, tmp_path):
    _scheme, target, _warnings = built()
    for key, value, differs in [("built", "2000-01-01", False),
                                ("built_with", "2.3.1", False),
                                ("status", "draft", True),
                                ("commit", "b" * 40, True),
                                ("known_issues", "A", True),
                                ("format", "2.5", True)]:
        changed = str(tmp_path / f"changed_{key}.db")
        shutil.copy(target, changed)
        connection = sqlite3.connect(changed)
        connection.execute("INSERT OR REPLACE INTO scheme_info (key, value) VALUES (?, ?)",
                           (key, value))
        connection.commit()
        connection.close()
        differences = build.compare(build.content(target), build.content(changed))
        assert bool(differences) == differs, (key, differences)
        if differs:
            assert [table for table, _, _ in differences] == ["scheme_info"]


@pytest.mark.parametrize("column,value", [("rank", 2), ("score", "480"), ("typed", 0)])
def test_a_database_that_differs_only_in_a_canSNP_column_is_not_equal(built, clean_git,
                                                                     tmp_path, column, value):
    """The comparison covers the rank, the score and the typed flag.

    Without them a database whose canSNPs are ranked or marked differently, which is
    what the scheme calls, would compare as equal to its rebuild.
    """
    _scheme, target, _warnings = built()
    changed = str(tmp_path / f"changed_{column}.db")
    shutil.copy(target, changed)
    connection = sqlite3.connect(changed)
    connection.execute(f"UPDATE snp_annotation SET {column} = ? WHERE snp_id = 'A'", (value,))
    connection.commit()
    connection.close()
    differences = build.compare(build.content(target), build.content(changed))
    assert [table for table, _, _ in differences] == ["snp_annotation"], differences


def test_a_database_without_scheme_info_differs_from_a_rebuild(built, clean_git, tmp_path):
    _scheme, target, _warnings = built()
    bare = str(tmp_path / "bare.db")
    shutil.copy(target, bare)
    connection = sqlite3.connect(bare)
    connection.execute("DROP TABLE scheme_info")
    connection.commit()
    connection.close()
    differences = build.compare(build.content(bare), build.content(target))
    assert [table for table, _, _ in differences] == ["scheme_info"]


# -- git state ---------------------------------------------------------------

def git(repo, *arguments):
    import subprocess
    subprocess.run(["git", "-C", str(repo), *arguments], check=True, capture_output=True,
                   env={**os.environ, "GIT_AUTHOR_NAME": "t", "GIT_AUTHOR_EMAIL": "t@t",
                        "GIT_COMMITTER_NAME": "t", "GIT_COMMITTER_EMAIL": "t@t"})


def head(repo):
    import subprocess
    return subprocess.run(["git", "-C", str(repo), "rev-parse", "HEAD"],
                          capture_output=True, text=True, check=True).stdout.strip()


@pytest.fixture
def checkout(make_scheme, tmp_path):
    """The conftest scheme committed to a git repository in tmp_path."""
    scheme = make_scheme()
    git(tmp_path, "init", "-q", "-b", "main")
    git(tmp_path, "add", "-A")
    git(tmp_path, "commit", "-q", "-m", "sources")
    return scheme, tmp_path


def test_the_source_commit_is_the_last_commit_of_the_scheme_files(checkout):
    scheme, repo = checkout
    sources = head(repo)
    (repo / "unrelated.txt").write_text("x\n")
    git(repo, "add", "unrelated.txt")
    git(repo, "commit", "-q", "-m", "unrelated")
    assert head(repo) != sources
    assert build.source_commit(scheme, repo=str(repo)) == sources


def test_a_changed_cansnp_file_is_uncommitted(checkout):
    scheme, repo = checkout
    assert build.dirty_sources(scheme, repo=str(repo)) == []
    with open(scheme.path(scheme.cansnps), "a") as handle:
        handle.write("D\t1\tTest2026\tG\tseq2\t250\tA\tC\t1\tyes\n")
    assert build.dirty_sources(scheme, repo=str(repo)) == ["test_canSNPs.txt"]


def test_a_changed_readme_is_uncommitted(checkout):
    scheme, repo = checkout
    with open(scheme.path(scheme.readme), "a") as handle:
        handle.write("More text.\n")
    assert build.dirty_sources(scheme, repo=str(repo)) == ["README.md"]


# -- manifests ---------------------------------------------------------------

def read_manifest(scheme):
    with open(build.manifest_path(scheme)) as handle:
        return json.load(handle)


def write_json(path, data):
    with open(path, "w") as handle:
        json.dump(data, handle)


def test_the_manifest_records_the_database_and_its_format(make_scheme, clean_git):
    scheme = make_scheme()
    assert build.process(scheme, write=True)
    manifest = read_manifest(scheme)
    assert manifest["database_sha256"] == build.file_digest(build.database_path(scheme))
    assert manifest["format"] == "2.4"
    assert manifest["source_commit"] == "a" * 40
    assert "uncommitted_sources" not in manifest


def test_the_manifest_hashes_the_genomes_file_when_there_is_one(make_scheme, clean_git):
    scheme = make_scheme()
    assert build.process(scheme, write=True)
    assert sorted(read_manifest(scheme)["sources"]) == sorted(
        [scheme.tree, scheme.cansnps, scheme.references])
    genomes = ("genome\tnode\trole\nGCA_000000002.1\tB\trepresentative\n")
    scheme = make_scheme(genomes=genomes)
    assert build.process(scheme, write=True)
    sources = read_manifest(scheme)["sources"]
    assert len(sources) == 4 and scheme.genomes in sources


def test_the_manifest_hashes_the_genomes_file_when_there_is_one(make_scheme, clean_git):
    scheme = make_scheme()
    assert build.process(scheme, write=True)
    assert sorted(read_manifest(scheme)["sources"]) == sorted(
        [scheme.tree, scheme.cansnps, scheme.references])
    genomes = "genome\tnode\trole\nGCA_000000002.1\tB\trepresentative\n"
    scheme = make_scheme(genomes=genomes)
    assert build.process(scheme, write=True)
    sources = read_manifest(scheme)["sources"]
    assert len(sources) == 4 and scheme.genomes in sources


def test_manifest_lifecycle(make_scheme, clean_git, monkeypatch):
    scheme = make_scheme()
    assert build.process(scheme, write=True)
    first = read_manifest(scheme)
    assert build.process(scheme, write=False)

    # an unchanged rebuild leaves the database and the manifest as they were,
    # not stamped with this run's tool versions
    monkeypatch.setattr(build, "tool_versions",
                        lambda scheme: {"cansnper2": "9.9.9", "python": "0"})
    assert build.process(scheme, write=True)
    assert read_manifest(scheme) == first

    # a manifest written before database_sha256 was recorded is completed;
    # the build description is kept
    stale = dict(first)
    del stale["database_sha256"]
    write_json(build.manifest_path(scheme), stale)
    assert not build.process(scheme, write=False)
    assert build.process(scheme, write=True)
    assert read_manifest(scheme) == first

    # a committed source edit, here to the reference column, which the database
    # does not store, changes the source commit and so the database
    rows = list(CANSNP_ROWS)
    rows[0] = "R\t1\tAnother2026\tG\tseq1\t100\tA\tG\t520\tyes"
    with open(scheme.path(scheme.cansnps), "w") as handle:
        handle.write(CANSNP_HEADER + "".join(row + "\n" for row in rows))
    monkeypatch.setattr(build, "source_commit", lambda scheme: "b" * 40)
    assert not build.process(scheme, write=False)
    assert build.process(scheme, write=True)
    second = read_manifest(scheme)
    assert second["sources"] != first["sources"]
    assert second["source_commit"] == "b" * 40
    assert second["database_sha256"] == build.file_digest(build.database_path(scheme))
    assert scheme_info_of(build.database_path(scheme))["commit"] == "b" * 40
    assert build.process(scheme, write=False)


def test_missing_manifest_fails_the_comparison(make_scheme, clean_git):
    scheme = make_scheme()
    assert build.process(scheme, write=True)
    os.remove(build.manifest_path(scheme))
    assert not build.process(scheme, write=False)


def test_a_database_that_does_not_match_its_manifest_fails(make_scheme, clean_git, capsys):
    scheme = make_scheme()
    assert build.process(scheme, write=True)
    manifest = read_manifest(scheme)
    manifest["database_sha256"] = "0" * 64
    write_json(build.manifest_path(scheme), manifest)
    capsys.readouterr()
    assert not build.process(scheme, write=False)
    assert "does not match the SHA-256 of test.db" in capsys.readouterr().out


def test_write_refuses_uncommitted_sources(make_scheme, clean_git, monkeypatch, capsys):
    scheme = make_scheme()
    monkeypatch.setattr(build, "dirty_sources", lambda scheme: ["schemes/test/test_canSNPs.txt"])
    monkeypatch.setattr(build, "discover", lambda only: [scheme])
    monkeypatch.setattr("sys.argv", ["build.py", "--write"])
    assert build.main() == 1
    out = capsys.readouterr().out
    assert ("FAIL schemes/test/test_canSNPs.txt has uncommitted changes; commit the "
            "sources first, then run build.py --write") in out
    assert not os.path.exists(build.database_path(scheme))
    assert not os.path.exists(build.manifest_path(scheme))


def test_write_refuses_sources_without_a_commit(make_scheme, clean_git, monkeypatch, capsys):
    # outside a git checkout, or for files never committed, no commit names the sources
    scheme = make_scheme()
    monkeypatch.setattr(build, "source_commit", lambda scheme: None)
    monkeypatch.setattr(build, "discover", lambda only: [scheme])
    monkeypatch.setattr("sys.argv", ["build.py", "--write"])
    assert build.main() == 1
    out = capsys.readouterr().out
    assert ("FAIL no commit of the source files of test was found; run build.py "
            "--write in a git checkout in which they are committed") in out
    assert not os.path.exists(build.database_path(scheme))
    assert not os.path.exists(build.manifest_path(scheme))
# -- versions ----------------------------------------------------------------

def test_version_tuple():
    assert build.version_tuple("2.2.0") == (2, 2, 0)
    assert build.version_tuple("2.10.1.dev3") == (2, 10, 1)
    assert build.version_tuple("unknown") == ()


def test_older_cansnper2_is_refused(monkeypatch):
    # 2.2 counts positions along the joined reference sequences
    monkeypatch.setattr(importlib.metadata, "version", lambda name: "2.2.2")
    with pytest.raises(RuntimeError, match="2.3.0 or later is required, found 2.2.2"):
        build.require_cansnper2()


def test_cansnper2_without_scheme_info_is_refused(monkeypatch):
    # a 2.3.0 development commit from before scheme_info has the same version
    real = build.importlib.import_module

    def without_scheme_info(name, *args, **kwargs):
        if name == "CanSNPer2.modules.SchemeInfo":
            raise ImportError(name)
        return real(name, *args, **kwargs)
    monkeypatch.setattr(build.importlib, "import_module", without_scheme_info)
    with pytest.raises(RuntimeError, match="does not write scheme_info"):
        build.require_cansnper2()
