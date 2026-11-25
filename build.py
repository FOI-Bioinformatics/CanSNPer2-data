#!/usr/bin/env python3
"""Build the CanSNPer2 databases in this repository from the scheme sources.

By default every scheme is built into a temporary directory and compared, row by
row, against the database committed under database/. Nothing in the working tree
is touched, which is what CI runs on a pull request:

    python build.py

Pass --write to install the result, archiving the database it replaces and
writing a provenance manifest beside it:

    python build.py --write
    python build.py --write --scheme bacillus_anthracis

--write refuses a scheme whose source files, README or schemes.toml have
uncommitted changes: the database records the commit of its sources in its
scheme_info table, so the sources are committed first and built after.

A scheme marked wip in schemes.toml, one whose build needs tools that continuous
integration does not install, is skipped unless it is named with --scheme or
--include-wip is given.

Requires CanSNPer2 2.3.0 or later, whose own loader builds a SNP scheme,
stores each SNP as a reference sequence and a position on that sequence, and
records the scheme in a scheme_info table. The flextaxd command is only needed
to build a taxonomy scheme.
"""

import argparse
import datetime
import hashlib
import importlib
import importlib.metadata
import json
import logging
import os
import re
import shutil
import sqlite3
import subprocess
import sys
import tempfile

from validate import CONFIG, REPO, discover, read_rows, tree_nodes_and_parents, validate

# Versions before 2.2.0 add an orphan node named 'root' beside the scheme root,
# so their databases neither match the committed ones nor pass check_structure.
# Versions before 2.3.0 count SNP positions along the reference sequences
# joined in file order and do not read the per-sequence SNP files. Development
# commits of 2.3.0 from before the scheme_info table carry the same version, so
# require_cansnper2 also checks that the module writing that table is present.
#: The schemes keep their canSNPs in <name>_canSNPs.txt, which CanSNPer2 reads from 2.4.0.
MINIMUM_CANSNPER2 = (2, 4, 0)

# The repository recorded as the source of every database built here.
DATA_REPOSITORY = "FOI-Bioinformatics/CanSNPer2-data"

# scheme_info keys that describe the build rather than the scheme: the date and
# the CanSNPer2 version, which the manifest records as well. They are left out
# of the comparison with the committed database.
BUILD_ONLY_KEYS = ("built", "built_with")


def database_path(scheme):
    return os.path.join(REPO, scheme.database_dir, f"{scheme.name}.db")


def manifest_path(scheme):
    return os.path.join(REPO, scheme.database_dir, f"{scheme.name}.manifest.json")


def archive_dir(scheme):
    return os.path.join(REPO, scheme.database_dir, "archive")


# --------------------------------------------------------------------------
# Building
# --------------------------------------------------------------------------

def version_tuple(version):
    """The leading numeric part of a version string, e.g. '2.2.0.dev1' -> (2, 2, 0)."""
    match = re.match(r"(\d+(?:\.\d+)*)", version)
    return tuple(int(part) for part in match.group(1).split(".")) if match else ()


def cansnper2_version():
    """The version of the CanSNPer2 that will be imported.

    The version of the package itself, not of the installed distribution: the two
    disagree for a checkout on PYTHONPATH or an editable install, and it is the
    code that will run that has to be able to read the schemes.
    """
    try:
        return importlib.import_module("CanSNPer2").__version__
    except (ImportError, AttributeError):
        pass
    return importlib.metadata.version("CanSNPer2")


def require_cansnper2():
    try:
        installed = cansnper2_version()
    except (importlib.metadata.PackageNotFoundError, ImportError):
        raise RuntimeError("CanSNPer2 is not installed")
    if version_tuple(installed) < MINIMUM_CANSNPER2:
        wanted = ".".join(str(part) for part in MINIMUM_CANSNPER2)
        raise RuntimeError(f"CanSNPer2 {wanted} or later is required, "
                           f"found {installed}")
    try:
        importlib.import_module("CanSNPer2.modules.SchemeInfo")
    except ImportError:
        raise RuntimeError(f"CanSNPer2 {installed} does not write scheme_info; "
                           f"install CanSNPer2 2.4.0 or later")


def scheme_source():
    """This repository as CanSNPer2 reads a local folder of schemes.

    CanSNPer2-download takes the status of a scheme from its README and the
    known issues from schemes.toml through this class, so a database built
    here records them in the same way.
    """
    from CanSNPer2.modules.FetchScheme import LocalSchemeSource
    return LocalSchemeSource(REPO)


def scheme_info(scheme):
    """The scheme_info keys of a database built from this scheme.

    ref is not written: the branch a database is built on goes stale once it is
    merged, and commit already identifies the sources.
    """
    source = scheme_source()
    commit = source_commit(scheme) or "unknown"
    if commit != "unknown" and dirty_sources(scheme):
        # the convention of CanSNPer2 for a checkout with uncommitted changes
        commit += "+uncommitted"
    return {
        "scheme": scheme.name,
        "status": source.status(scheme.readme) or "unknown",
        "source": DATA_REPOSITORY,
        "commit": commit,
        "known_issues": source.known_issues(scheme.name),
    }


def cansnper2_loads_genomes():
    """Whether the CanSNPer2 that will be imported loads a genomes file into the
    database (the genomes option of CanSNPer2-database, with read typing)."""
    try:
        importlib.import_module("CanSNPer2.modules.SchemeGenomes")
    except ImportError:
        return False
    return True


def build_snp_scheme(scheme, target):
    """Build a canSNP database with CanSNPer2's own loader.

    CanSNPer2.SNPDatabase parses command line arguments at import time, so the
    loader class is imported directly instead of going through that module.
    The genomes file of the scheme, when it has one, is loaded too: read typing
    with local references takes its representative genomes from the database.
    A CanSNPer2 that cannot load it builds the database without it, as before.
    """
    require_cansnper2()
    from CanSNPer2.modules.InitiateDatabase import CanSNPDatabase

    options = {}
    genomes = scheme.path(scheme.genomes)
    if os.path.exists(genomes):
        if cansnper2_loads_genomes():
            options["genomes"] = genomes
        else:
            # not a warning: build() takes the warnings of a build as problems,
            # and a database without genomes is what this CanSNPer2 builds
            logging.getLogger(__name__).info(
                f"CanSNPer2 {cansnper2_version()} cannot load the genomes file "
                f"{scheme.genomes}, the database is built without it")

    CanSNPDatabase(
        database=target,
        create=True,
        tree=scheme.path(scheme.tree),
        annotation=scheme.path(scheme.cansnps),
        references=scheme.path(scheme.references),
        source_type=scheme.source_type,
        verbose=False,
        scheme_info=scheme_info(scheme),
        **options,
    )


def flextaxd_command():
    """Locate the flextaxd console script.

    Prefer the one installed alongside the interpreter running this script, so
    that a virtual environment is honoured even when its bin directory is not
    on PATH.
    """
    local = os.path.join(os.path.dirname(sys.executable), "flextaxd")
    if os.path.exists(local):
        return local
    found = shutil.which("flextaxd")
    if found:
        return found
    raise RuntimeError("flextaxd was not found; install flextaxd to build "
                       "taxonomy schemes")


def build_taxonomy_scheme(scheme, target):
    """Build a taxonomy database with the flextaxd command.

    CanSNPer2-database has no option for a genome-to-node file, so a scheme of
    this kind goes through flextaxd directly. --force suppresses the interactive
    confirmation that flextaxd asks for before creating a database.
    """
    command = [
        flextaxd_command(), "--force",
        "--database", target,
        "--taxonomy_file", scheme.path(scheme.tree),
        "--taxonomy_type", scheme.source_type,
        "--genomeid2taxid", scheme.path(scheme.genomes),
    ]
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"flextaxd failed for {scheme.name}:\n"
                           f"{result.stdout}\n{result.stderr}")


class WarningCollector(logging.Handler):
    """Keep the warnings a loader logs while it builds one database."""

    def __init__(self):
        super().__init__(level=logging.WARNING)
        self.messages = []

    def emit(self, record):
        self.messages.append(record.getMessage())


def build(scheme, target):
    """Build one scheme into target and return the warnings its loader logged.

    The CanSNPer2 loader skips a row it cannot use (a short row, a SNP that is
    not in the tree, a genome without a reference) with a logged warning only,
    so the warnings are collected here whatever the console log level is.
    """
    root = logging.getLogger()
    collector = WarningCollector()
    previous_level = root.level
    root.setLevel(min(previous_level, logging.WARNING))
    root.addHandler(collector)
    try:
        if scheme.kind == "taxonomy":
            build_taxonomy_scheme(scheme, target)
        else:
            build_snp_scheme(scheme, target)
    finally:
        root.removeHandler(collector)
        root.setLevel(previous_level)
    return collector.messages


def check_structure(target, single_root=True):
    """Check the tree in a freshly built database for orphans and cycles.

    flextaxd's own validate_tree is not used: it requires every node to descend
    from node id 1, which does not hold for databases built by CanSNPer2 before
    2.2.0 or by flextaxd, where the scheme root is a self-linked sibling of an
    orphan node 1 named 'root'. CanSNPer2 2.2.0 makes the scheme root node 1.
    The checks below hold for both layouts: every self-linked node is treated
    as a root. With single_root, which applies to SNP schemes built by
    CanSNPer2 2.2.0, exactly one root is allowed.

    Returns a list of problems, empty when the tree is sound.
    """
    problems = []
    connection = sqlite3.connect(target)
    try:
        names = dict(connection.execute("SELECT id, name FROM nodes"))
        edges = [(parent, child) for parent, child in connection.execute(
            "SELECT DISTINCT parent, child FROM tree") if parent != child]
        self_linked = {row[0] for row in connection.execute(
            "SELECT parent FROM tree WHERE parent = child")}
        duplicated = [row[0] for row in connection.execute(
            "SELECT name FROM nodes GROUP BY name HAVING count(*) > 1 ORDER BY name")]
    finally:
        connection.close()

    # The tree reader creates a second node of the same name rather than
    # rejecting a node that a source file places under two parents.
    for name in duplicated:
        problems.append(f"node name '{name}' is used by more than one node")

    if not self_linked:
        problems.append("no root: no node is linked to itself")
    elif single_root and len(self_linked) > 1:
        listed = ", ".join(sorted(names.get(n, str(n)) for n in self_linked))
        problems.append(f"{len(self_linked)} roots where one is expected: {listed}")

    parents = {}
    for parent, child in edges:
        parents.setdefault(child, set()).add(parent)

    for child, child_parents in sorted(parents.items()):
        if len(child_parents) > 1:
            listed = ", ".join(sorted(names.get(p, str(p)) for p in child_parents))
            problems.append(f"node '{names.get(child, child)}' has "
                            f"{len(child_parents)} parents: {listed}")

    # Everything must be reachable from a self-linked node; those are the roots.
    children = {}
    for parent, child in edges:
        children.setdefault(parent, set()).add(child)
    reached, frontier = set(self_linked), list(self_linked)
    while frontier:
        for child in children.get(frontier.pop(), ()):
            if child not in reached:
                reached.add(child)
                frontier.append(child)

    for orphan in sorted(set(names) - reached):
        problems.append(f"node '{names.get(orphan, orphan)}' is not reachable "
                        f"from any root")

    # A cycle shows up as a node reachable from itself.
    for start in sorted(names):
        seen, stack = set(), [start]
        while stack:
            for child in children.get(stack.pop(), ()):
                if child == start:
                    problems.append(f"node '{names.get(start, start)}' is part "
                                    f"of a cycle")
                    stack = []
                    break
                if child not in seen:
                    seen.add(child)
                    stack.append(child)

    return problems


def check_counts(scheme, target):
    """Check that every row of a SNP scheme's sources reached the database.

    A node may hold several canSNPs, so snp_annotation holds one row per canSNP
    rather than one per node.

    The loader skips rows it cannot use, so a database that is missing rows
    can still be built without an error; for a new scheme there is no
    committed database for the comparison to catch it against.
    """
    if scheme.kind != "snp":
        return []
    _, _, tree_rows = read_rows(scheme.path(scheme.tree))
    tree_nodes, _ = tree_nodes_and_parents(tree_rows)
    expected = {
        "nodes": len(tree_nodes),
        "snp_annotation": len(read_rows(scheme.path(scheme.cansnps))[2]) - 1,
        "snp_references": len(read_rows(scheme.path(scheme.references))[2]) - 1,
    }
    found = row_counts(target)
    return [f"{table} has {found.get(table, 0)} rows, the sources define {count}"
            for table, count in expected.items() if found.get(table, 0) != count]


# --------------------------------------------------------------------------
# Comparison
# --------------------------------------------------------------------------

def content(path):
    """Read a database into comparable sets, keyed on names rather than row ids.

    snp_annotation.date holds the build date and the surrogate keys depend on
    insertion order, so neither is part of the comparison, and neither are the
    scheme_info keys in BUILD_ONLY_KEYS. Each SNP is compared
    together with the name of the node it is attached to, and node names that
    occur more than once are listed, since a set of names would merge them.
    """
    connection = sqlite3.connect(path)
    query = connection.execute
    try:
        tables = {row[0] for row in query(
            "SELECT name FROM sqlite_master WHERE type = 'table'")}
        result = {
            "nodes": {row[0] for row in query("SELECT name FROM nodes")},
            "duplicated node names": {row[0] for row in query(
                "SELECT name FROM nodes GROUP BY name HAVING count(*) > 1")},
            "edges": set(query("""SELECT p.name, c.name, r.rank FROM tree t
                                  JOIN nodes p ON p.id = t.parent
                                  JOIN nodes c ON c.id = t.child
                                  LEFT JOIN rank r ON r.rank_i = t.rank_i""")),
            "rank": {row[0] for row in query("SELECT rank FROM rank")},
        }
        # The genomes table gained a "reference" column in later flextaxd
        # versions, so only the columns common to every vintage are compared,
        # and the role that CanSNPer2 stores for the genomes of a canSNP scheme.
        genome_columns = {row[1] for row in query("PRAGMA table_info(genomes)")}
        role = "g.role" if "role" in genome_columns else "NULL"
        result["genomes"] = set(query(f"""SELECT n.name, g.genome, {role} FROM genomes g
                                          LEFT JOIN nodes n ON n.id = g.id"""))
        if "snp_annotation" in tables:
            # A database built before CanSNPer2 2.3.0 has a single position
            # column, counted along the joined reference sequences. It is read
            # with no sequence, so that it differs from a rebuild instead of
            # stopping the comparison.
            columns = {row[1] for row in query("PRAGMA table_info(snp_annotation)")}
            if "sequence_position" in columns:
                location = "a.sequence, a.sequence_position"
            else:
                location = "NULL, a.position"
            # A database built before the canSNP file holds one canSNP per node and
            # none of these three columns; it is read without them, so that it differs
            # from a rebuild instead of stopping the comparison.
            if "typed" in columns:
                ranking = "a.rank, a.score, a.typed"
            else:
                ranking = "NULL, NULL, NULL"
            result["snp_annotation"] = set(query(
                f"""SELECT n.name, a.snp_id, {location}, a.derived_base,
                          a.ancestral_base, a.reference, g.genome, {ranking}
                   FROM snp_annotation a
                   LEFT JOIN nodes n ON n.id = a.node_id
                   LEFT JOIN snp_references g ON g.id = a.genome_i"""))
            result["snp_references"] = set(query(
                """SELECT genome, strain, genbank_id, refseq_id, assembly_name
                   FROM snp_references"""))
        if "snp_annotation" in tables or "scheme_info" in tables:
            # A database without the table, built before CanSNPer2 2.3.0,
            # gives an empty set and so differs from a rebuild.
            result["scheme_info"] = set()
            if "scheme_info" in tables:
                result["scheme_info"] = {
                    (key, value) for key, value in query("SELECT key, value FROM scheme_info")
                    if key not in BUILD_ONLY_KEYS}
        return result
    finally:
        connection.close()


def compare(committed, built):
    """Return a list of (table, only_in_committed, only_in_built) differences."""
    differences = []
    for table in sorted(set(committed) | set(built)):
        left = committed.get(table, set())
        right = built.get(table, set())
        only_left = sorted(str(v) for v in left - right)
        only_right = sorted(str(v) for v in right - left)
        if only_left or only_right:
            differences.append((table, only_left, only_right))
    return differences


def report_differences(differences, limit=10):
    for table, only_committed, only_built in differences:
        print(f"    {table}: {len(only_committed)} only in the committed "
              f"database, {len(only_built)} only in the rebuild")
        for value in only_committed[:limit]:
            print(f"      committed only: {value}")
        for value in only_built[:limit]:
            print(f"      rebuild only:   {value}")
        hidden = (len(only_committed) + len(only_built)) - min(
            limit, len(only_committed)) - min(limit, len(only_built))
        if hidden > 0:
            print(f"      ... and {hidden} more")


# --------------------------------------------------------------------------
# Provenance
# --------------------------------------------------------------------------

def file_digest(path):
    with open(path, "rb") as handle:
        return hashlib.sha256(handle.read()).hexdigest()


def git(repo, *arguments):
    """The output of a git command in repo, or None when it fails."""
    try:
        result = subprocess.run(["git", "-C", repo] + list(arguments),
                                capture_output=True, text=True)
    except OSError:
        return None
    return result.stdout if result.returncode == 0 else None


def inside(repo, paths):
    """The paths that lie in repo; git refuses a path outside its checkout."""
    top = os.path.realpath(repo)
    return [path for path in paths
            if os.path.commonpath([top, os.path.realpath(path)]) == top]


def source_commit(scheme, repo=REPO):
    """The last commit that changed the source files of the scheme.

    This, rather than HEAD, is recorded, so that committing the rebuilt
    database, or a change to another scheme, does not make the database
    describe a different commit. None when it cannot be found.
    """
    paths = inside(repo, [scheme.path(relative) for relative in scheme.files])
    if not paths:
        return None
    commit = (git(repo, "log", "-1", "--format=%H", "--", *paths) or "").strip()
    return commit or None


def recorded_inputs(scheme):
    """Every file whose content the database records: the source files, and
    the README and schemes.toml that give its status and known issues."""
    paths = [scheme.path(relative) for relative in scheme.files]
    if scheme.readme:
        paths.append(scheme.path(scheme.readme))
    paths.append(CONFIG)
    return paths


def dirty_sources(scheme, repo=REPO):
    """The files of recorded_inputs that have uncommitted changes or are not
    committed at all, relative to the top of the checkout."""
    paths = inside(repo, recorded_inputs(scheme))
    if not paths:
        return []
    status = git(repo, "status", "--porcelain", "--untracked-files=all", "--", *paths)
    if status is None:
        return []
    return sorted(line[3:].strip() for line in status.splitlines() if line.strip())


def tool_versions(scheme):
    """Record the versions that produced a build.

    CanSNPer2's version is the one the gate reads, the version of the package
    that is imported: a checkout on PYTHONPATH or an editable install builds
    with a version the installed distribution metadata does not know, and a
    manifest naming the other one would credit a version that did not build the
    database. flextaxd is read from the distribution metadata, as its package
    does not expose a __version__ attribute, and it is only recorded for a
    scheme that was built with it.
    """
    versions = {"python": sys.version.split()[0]}
    try:
        versions["cansnper2"] = cansnper2_version()
    except importlib.metadata.PackageNotFoundError:
        versions["cansnper2"] = "unknown"
    if scheme.kind == "taxonomy":
        try:
            versions["flextaxd"] = importlib.metadata.version("flextaxd")
        except importlib.metadata.PackageNotFoundError:
            versions["flextaxd"] = "unknown"
    return versions


def row_counts(path):
    connection = sqlite3.connect(path)
    try:
        tables = sorted(row[0] for row in connection.execute(
            "SELECT name FROM sqlite_master WHERE type = 'table'"))
        return {table: connection.execute(
            f"SELECT count(*) FROM {table}").fetchone()[0] for table in tables}
    finally:
        connection.close()


def source_digests(scheme):
    return {relative: file_digest(scheme.path(relative)) for relative in scheme.files}


def read_manifest(scheme):
    path = manifest_path(scheme)
    if not os.path.exists(path):
        return None
    with open(path) as handle:
        return json.load(handle)


def save_manifest(scheme, manifest):
    with open(manifest_path(scheme), "w") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")


def database_record(scheme, target):
    """The manifest keys that describe the database file itself: its SHA-256
    and, for a database that records one, its format."""
    record = {"database_sha256": file_digest(target)}
    if scheme.kind == "snp":
        from CanSNPer2.modules.SchemeInfo import read_scheme_info
        info_format = read_scheme_info(target).get("format")
        if info_format:
            record["format"] = info_format
    return record


def write_manifest(scheme, target):
    """Record what produced a newly installed database, as text beside the binary."""
    manifest = {
        "scheme": scheme.name,
        "kind": scheme.kind,
        "built": datetime.date.today().isoformat(),
        "source_commit": source_commit(scheme),
        "sources": source_digests(scheme),
        "versions": tool_versions(scheme),
        "counts": row_counts(target),
    }
    manifest.update(database_record(scheme, target))
    save_manifest(scheme, manifest)
    return manifest


def manifest_is_current(scheme, manifest, target):
    """True if the manifest records the current sources and the database file."""
    if manifest.get("sources", {}) != source_digests(scheme):
        return False
    if manifest.get("source_commit") != source_commit(scheme):
        return False
    record = database_record(scheme, target)
    return all(manifest.get(key) == value for key, value in record.items())


def refresh_manifest(scheme, manifest, target):
    """Update the record of the sources and of the database file for a
    database that a rebuild left unchanged, for example a manifest written
    before database_sha256 was recorded. built, versions and counts still
    describe the earlier build, which produced the database, and are kept.
    """
    manifest["sources"] = source_digests(scheme)
    manifest["source_commit"] = source_commit(scheme)
    manifest.pop("uncommitted_sources", None)
    manifest.update(database_record(scheme, target))
    save_manifest(scheme, manifest)
    return manifest


def sources_changed(scheme):
    """True if the sources no longer hash to what the manifest recorded."""
    manifest = read_manifest(scheme)
    if manifest is None:
        return True
    return manifest.get("sources", {}) != source_digests(scheme)


def install(scheme, target):
    """Archive the database being replaced, then move the new one into place."""
    destination = database_path(scheme)
    if os.path.exists(destination):
        archive = archive_dir(scheme)
        os.makedirs(archive, exist_ok=True)
        stamp = datetime.date.today().strftime("%y%m%d")
        archived = os.path.join(archive, f"{scheme.name}_{stamp}.db")
        suffix = 1
        while os.path.exists(archived):
            suffix += 1
            archived = os.path.join(archive, f"{scheme.name}_{stamp}-{suffix}.db")
        shutil.move(destination, archived)
        print(f"    archived the previous database to "
              f"{os.path.relpath(archived, REPO)}")
    shutil.move(target, destination)
    print(f"    wrote {os.path.relpath(destination, REPO)}")


# --------------------------------------------------------------------------

def process(scheme, write):
    """Build one scheme. Returns True if it is in the state we want."""
    print(f"{scheme.name}{' (built on request)' if scheme.wip else ''}:")

    problems, warnings = validate(scheme)
    for warning in warnings:
        print(f"    known issue: {warning}")
    if problems:
        print(f"    FAIL validation ({len(problems)} problem"
              f"{'s' if len(problems) != 1 else ''}); run validate.py for detail")
        for problem in problems[:5]:
            print(f"      - {problem}")
        return False

    if write:
        # The database records the source commit, which does not describe
        # files with uncommitted changes.
        dirty = dirty_sources(scheme)
        if dirty:
            listed = ", ".join(dirty)
            print(f"    FAIL {listed} {'has' if len(dirty) == 1 else 'have'} "
                  f"uncommitted changes; commit the sources first, then run "
                  f"build.py --write")
            return False
        if source_commit(scheme) is None:
            print(f"    FAIL no commit of the source files of {scheme.name} was "
                  f"found; run build.py --write in a git checkout in which they "
                  f"are committed")
            return False

    with tempfile.TemporaryDirectory(prefix=f"cansnper2-{scheme.name}-") as workdir:
        target = os.path.join(workdir, f"{scheme.name}.db")
        try:
            loader_warnings = build(scheme, target)
        except Exception as error:
            print(f"    FAIL build: {error}")
            return False
        if loader_warnings:
            print(f"    FAIL the loader reported {len(loader_warnings)} warning"
                  f"{'s' if len(loader_warnings) != 1 else ''}")
            for message in loader_warnings[:10]:
                print(f"      - {message}")
            return False

        structure = (check_structure(target, single_root=scheme.kind == "snp")
                     + check_counts(scheme, target))
        if structure:
            print("    FAIL structure of the built database")
            for problem in structure:
                print(f"      - {problem}")
            return False
        print("    built, tree structure sound")

        committed = database_path(scheme)
        matched = False
        if not os.path.exists(committed):
            print("    no committed database to compare against")
        else:
            # The comparison always runs, because it decides whether a write
            # needs to replace the database at all. scheme.verify only decides
            # whether a difference is a failure.
            differences = compare(content(committed), content(target))
            if not differences:
                matched = True
                print("    matches the committed database")
            elif scheme.verify:
                print("    FAIL the rebuild does not match the committed database")
                report_differences(differences)
                if not write:
                    return False
            else:
                print("    differs from the committed database; not treated as a "
                      "failure for this scheme (verify = false)")
                report_differences(differences, limit=3)

        manifest = read_manifest(scheme)
        manifest_name = os.path.basename(manifest_path(scheme))
        if not write:
            # The comparison above cannot see a source edit that does not reach
            # the database, so the recorded hashes are checked separately.
            if os.path.exists(committed):
                if manifest is None:
                    print(f"    FAIL {manifest_name} is missing; run build.py --write")
                    return False
                if manifest.get("sources", {}) != source_digests(scheme):
                    print(f"    FAIL {manifest_name} does not record the current "
                          f"sources; run build.py --write")
                    return False
                if manifest.get("database_sha256") != file_digest(committed):
                    print(f"    FAIL {manifest_name} does not match the SHA-256 of "
                          f"{os.path.basename(committed)}; run build.py --write")
                    return False
            return True

        if matched:
            # The committed database already holds this content. Replacing
            # it would change the SQLite bytes, and so produce a diff, for
            # no change in meaning.
            print("    database already up to date, left in place")
            if manifest is not None:
                if manifest_is_current(scheme, manifest, committed):
                    # Rewriting it would put this run's date and tool versions
                    # on a database that an earlier build produced.
                    print(f"    {manifest_name} already up to date, left in place")
                    return True
                refresh_manifest(scheme, manifest, committed)
                print(f"    updated the record of the sources and the database "
                      f"in {manifest_name}")
                return True
        else:
            install(scheme, target)
        manifest = write_manifest(scheme, committed)
        print(f"    wrote {manifest_name} "
              f"({', '.join(f'{k}={v}' for k, v in sorted(manifest['counts'].items()))})")
    return True


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--scheme", action="append", metavar="NAME",
                        help="build only this scheme; repeatable")
    parser.add_argument("--write", action="store_true",
                        help="install the built databases and write manifests")
    parser.add_argument("--include-wip", action="store_true",
                        help="also build schemes marked wip in schemes.toml")
    parser.add_argument("--changed-only", action="store_true",
                        help="skip schemes whose sources still hash to their manifest")
    parser.add_argument("--verbose", action="store_true",
                        help="show the logging from the build tools")
    args = parser.parse_args()

    level = logging.INFO if args.verbose else logging.ERROR
    logging.basicConfig(level=level, format="%(levelname)s %(name)s: %(message)s")
    # build() lowers the root logger to WARNING to collect loader warnings; the
    # console keeps the level asked for, and the warnings are reported by process().
    for handler in logging.getLogger().handlers:
        handler.setLevel(level)

    schemes = discover(args.scheme)
    if not args.scheme and not args.include_wip:
        for scheme in schemes:
            if scheme.wip:
                print(f"{scheme.name}: built on request only, skipping "
                      f"(use --scheme {scheme.name} or --include-wip)")
        schemes = [s for s in schemes if not s.wip]
    if args.changed_only:
        skipped = [s for s in schemes if not sources_changed(s)]
        for scheme in skipped:
            print(f"{scheme.name}: sources unchanged since the manifest, skipping")
        schemes = [s for s in schemes if sources_changed(s)]
        if not schemes:
            print(f"\nNothing to do; {len(skipped)} scheme(s) up to date.")
            return 0

    failed = [s.name for s in schemes if not process(s, args.write)]
    if failed:
        print(f"\n{len(failed)} of {len(schemes)} schemes failed: "
              f"{', '.join(failed)}")
        return 1
    print(f"\nAll {len(schemes)} schemes "
          f"{'are installed and up to date' if args.write else 'match their sources'}.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
