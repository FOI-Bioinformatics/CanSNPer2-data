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

A scheme marked as work in progress in schemes.toml is skipped unless it is named
with --scheme or --include-wip is given.

Requires CanSNPer2 2.2.0 or later (the master branch), whose own loader builds
a SNP scheme. The flextaxd command is only needed to build a taxonomy scheme, which at
present means the work-in-progress francisellaceae scheme.
"""

import argparse
import datetime
import hashlib
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

from validate import REPO, discover, read_rows, tree_nodes_and_parents, validate

# Earlier versions add an orphan node named 'root' beside the scheme root, so
# their databases neither match the committed ones nor pass check_structure.
MINIMUM_CANSNPER2 = (2, 2, 0)


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


def require_cansnper2():
    try:
        installed = importlib.metadata.version("CanSNPer2")
    except importlib.metadata.PackageNotFoundError:
        raise RuntimeError("CanSNPer2 is not installed")
    if version_tuple(installed) < MINIMUM_CANSNPER2:
        wanted = ".".join(str(part) for part in MINIMUM_CANSNPER2)
        raise RuntimeError(f"CanSNPer2 {wanted} or later is required, "
                           f"found {installed}")


def build_snp_scheme(scheme, target):
    """Build a canSNP database with CanSNPer2's own loader.

    CanSNPer2.SNPDatabase parses command line arguments at import time, so the
    loader class is imported directly instead of going through that module.
    """
    require_cansnper2()
    from CanSNPer2.modules.InitiateDatabase import CanSNPDatabase

    CanSNPDatabase(
        database=target,
        create=True,
        tree=scheme.path(scheme.tree),
        annotation=scheme.path(scheme.snp),
        references=scheme.path(scheme.references),
        source_type=scheme.source_type,
        verbose=False,
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
        "snp_annotation": len(read_rows(scheme.path(scheme.snp))[2]) - 1,
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
    insertion order, so neither is part of the comparison. Each SNP is compared
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
            # The genomes table gained a "reference" column in later flextaxd
            # versions, so only the columns common to every vintage are compared.
            "genomes": set(query("""SELECT n.name, g.genome FROM genomes g
                                    LEFT JOIN nodes n ON n.id = g.id""")),
        }
        if "snp_annotation" in tables:
            result["snp_annotation"] = set(query(
                """SELECT n.name, a.snp_id, a.position, a.derived_base,
                          a.ancestral_base, a.reference, g.genome
                   FROM snp_annotation a
                   LEFT JOIN nodes n ON n.id = a.node_id
                   LEFT JOIN snp_references g ON g.id = a.genome_i"""))
            result["snp_references"] = set(query(
                """SELECT genome, strain, genbank_id, refseq_id, assembly_name
                   FROM snp_references"""))
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


def source_commit():
    try:
        result = subprocess.run(["git", "-C", REPO, "rev-parse", "HEAD"],
                                capture_output=True, text=True)
        if result.returncode == 0:
            return result.stdout.strip()
    except OSError:
        pass
    return None


def dirty_sources(scheme):
    """Return this scheme's source files that have uncommitted changes."""
    paths = [scheme.path(relative) for relative in scheme.files]
    result = subprocess.run(["git", "-C", REPO, "status", "--porcelain", "--"] + paths,
                            capture_output=True, text=True)
    if result.returncode != 0:
        return []
    return [line[3:].strip() for line in result.stdout.splitlines() if line.strip()]


def tool_versions(scheme):
    """Record the installed versions that produced a build.

    The distribution metadata is the reliable source: flextaxd's package does
    not expose a __version__ attribute. flextaxd is only recorded for a scheme
    that was built with it.
    """
    versions = {"python": sys.version.split()[0]}
    tools = [("cansnper2", "CanSNPer2")]
    if scheme.kind == "taxonomy":
        tools.append(("flextaxd", "flextaxd"))
    for name, distribution in tools:
        try:
            versions[name] = importlib.metadata.version(distribution)
        except importlib.metadata.PackageNotFoundError:
            versions[name] = "unknown"
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


def write_manifest(scheme, target):
    """Record what produced a newly installed database, as text beside the binary."""
    manifest = {
        "scheme": scheme.name,
        "kind": scheme.kind,
        "built": datetime.date.today().isoformat(),
        "source_commit": source_commit(),
        "uncommitted_sources": dirty_sources(scheme),
        "sources": source_digests(scheme),
        "versions": tool_versions(scheme),
        "counts": row_counts(target),
    }
    save_manifest(scheme, manifest)
    return manifest


def refresh_manifest_sources(scheme, manifest):
    """Record new source hashes for a database that a rebuild left unchanged.

    This happens when a source edit does not reach the database, for example
    the strain column of a SNP file, which CanSNPer2 does not store. The
    database, and so built, versions and counts, still describe the earlier
    build; only the record of which sources it corresponds to is updated.
    """
    manifest["sources"] = source_digests(scheme)
    manifest["source_commit"] = source_commit()
    manifest["uncommitted_sources"] = dirty_sources(scheme)
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
    print(f"{scheme.name}{' (work in progress)' if scheme.wip else ''}:")

    problems, warnings = validate(scheme)
    for warning in warnings:
        print(f"    known issue: {warning}")
    if problems:
        print(f"    FAIL validation ({len(problems)} problem"
              f"{'s' if len(problems) != 1 else ''}); run validate.py for detail")
        for problem in problems[:5]:
            print(f"      - {problem}")
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
            return True

        if matched:
            # The committed database already holds this content. Replacing
            # it would change the SQLite bytes, and so produce a diff, for
            # no change in meaning.
            print("    database already up to date, left in place")
            if manifest is not None:
                if manifest.get("sources", {}) == source_digests(scheme):
                    # Rewriting it would put this run's commit, date and tool
                    # versions on a database that an earlier build produced.
                    print(f"    {manifest_name} already up to date, left in place")
                    return True
                manifest = refresh_manifest_sources(scheme, manifest)
                print(f"    updated the source hashes in {manifest_name}")
                note_uncommitted(manifest)
                return True
        else:
            install(scheme, target)
        manifest = write_manifest(scheme, committed)
        print(f"    wrote {manifest_name} "
              f"({', '.join(f'{k}={v}' for k, v in sorted(manifest['counts'].items()))})")
        note_uncommitted(manifest)
    return True


def note_uncommitted(manifest):
    if manifest["uncommitted_sources"]:
        print("    note: sources have uncommitted changes, so "
              "source_commit does not fully describe this build")


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--scheme", action="append", metavar="NAME",
                        help="build only this scheme; repeatable")
    parser.add_argument("--write", action="store_true",
                        help="install the built databases and write manifests")
    parser.add_argument("--include-wip", action="store_true",
                        help="also build schemes marked as work in progress")
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
                print(f"{scheme.name}: work in progress, skipping "
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
