"""Tests for build.py. They build real databases, so they need CanSNPer2 2.2.0 or later."""

import importlib.metadata
import json
import os
import shutil
import sqlite3

import pytest

import build
from conftest import SNP_HEADER, SNP_ROWS

try:
    build.require_cansnper2()
except RuntimeError as error:
    pytest.skip(f"needs CanSNPer2: {error}", allow_module_level=True)


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


def test_loader_warnings_are_collected_and_rows_counted(built):
    # a SNP on a genome without a reference is skipped by the loader with a
    # logged warning only; build() returns it and check_counts sees the gap
    rows = list(SNP_ROWS)
    rows[3] = "C\ts\tTest2026\tOTHER\t400\tG\tA"
    scheme, target, warnings = built(snp_rows=rows)
    assert warnings, "the skipped row must be reported"
    assert build.check_counts(scheme, target) == [
        "snp_annotation has 3 rows, the sources define 4"]


def test_node_with_two_parents_gives_a_repeated_name(built):
    # validate.py rejects this source; the database check is the second line
    tree = "R\nR\tA\nR\tC\nR\tA\tB\nR\tC\tB\n"
    _scheme, target, _warnings = built(tree=tree)
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


# -- manifests ---------------------------------------------------------------

def read_manifest(scheme):
    with open(build.manifest_path(scheme)) as handle:
        return json.load(handle)


def test_manifest_lifecycle(make_scheme, monkeypatch):
    scheme = make_scheme()
    monkeypatch.setattr(build, "source_commit", lambda: "commit-1")
    assert build.process(scheme, write=True)
    first = read_manifest(scheme)
    assert first["source_commit"] == "commit-1"
    assert os.path.exists(build.database_path(scheme))
    assert build.process(scheme, write=False)

    # an unchanged rebuild leaves the manifest as it was, not stamped with
    # this run's commit and tool versions
    monkeypatch.setattr(build, "source_commit", lambda: "commit-2")
    monkeypatch.setattr(build, "tool_versions",
                        lambda scheme: {"cansnper2": "9.9.9", "python": "0"})
    assert build.process(scheme, write=True)
    assert read_manifest(scheme) == first

    # the strain column is not stored, so the database still matches; the
    # comparison mode must still notice that the manifest is out of date
    rows = list(SNP_ROWS)
    rows[0] = "R\tanother strain\tTest2026\tG\t100\tA\tG"
    with open(scheme.path(scheme.snp), "w") as handle:
        handle.write(SNP_HEADER + "".join(row + "\n" for row in rows))
    assert not build.process(scheme, write=False)

    # --write updates only the source record; the build description is kept
    assert build.process(scheme, write=True)
    second = read_manifest(scheme)
    assert second["sources"] != first["sources"]
    assert second["source_commit"] == "commit-2"
    for key in ("built", "versions", "counts"):
        assert second[key] == first[key]
    assert build.process(scheme, write=False)


def test_missing_manifest_fails_the_comparison(make_scheme):
    scheme = make_scheme()
    assert build.process(scheme, write=True)
    os.remove(build.manifest_path(scheme))
    assert not build.process(scheme, write=False)


# -- versions ----------------------------------------------------------------

def test_version_tuple():
    assert build.version_tuple("2.2.0") == (2, 2, 0)
    assert build.version_tuple("2.10.1.dev3") == (2, 10, 1)
    assert build.version_tuple("unknown") == ()


def test_older_cansnper2_is_refused(monkeypatch):
    monkeypatch.setattr(importlib.metadata, "version", lambda name: "2.1.0")
    with pytest.raises(RuntimeError, match="2.2.0 or later is required, found 2.1.0"):
        build.require_cansnper2()
