#!/usr/bin/env python3
"""Validate the canSNP scheme sources in this repository.

Checks the text files under source/ for format problems and for the cross-file
invariants that the database loaders assume but do not enforce. Run with no
arguments to check every scheme:

    python validate.py
    python validate.py --scheme francisella_tularensis

Exits non-zero if any scheme has a problem.
"""

import argparse
import codecs
import inspect
import os
import re
import sys
import tomllib
from collections import Counter

REPO = os.path.dirname(os.path.abspath(__file__))
SOURCE_DIR = os.path.join(REPO, "source")
CONFIG = os.path.join(REPO, "schemes.toml")

SNP_HEADER = ["snp_id", "strain", "reference", "genome", "position",
              "derived_base", "ancestral_base"]
REFERENCE_HEADER = ["genome", "strain", "genbank_id", "refseq_id", "assembly_name"]

BASES = set("ACGT")
GENBANK_RE = re.compile(r"^GCA_\d{9}\.\d{1,2}$")
REFSEQ_RE = re.compile(r"^GCF_\d{9}\.\d{1,2}$")
# ASCII digits only: str.isdigit() also accepts characters such as superscripts,
# which int() then rejects.
POSITION_RE = re.compile(r"^[0-9]+$")


class Scheme:
    """One canSNP or taxonomy scheme and the source files that define it."""

    def __init__(self, name, kind="snp", tree=None, snp=None, genomes=None,
                 references=None, verify=True, source_type="CanSNPer",
                 known_issues=None, wip=False, database_dir="database"):
        self.name = name
        self.kind = kind
        self.verify = verify
        self.source_type = source_type
        self.known_issues = set(known_issues or ())
        self.wip = wip
        self.database_dir = database_dir
        self.tree = tree or os.path.join("source", f"{name}_tree.txt")
        self.snp = snp or os.path.join("source", f"{name}_snp.txt")
        self.references = references or os.path.join("source", f"{name}_references.txt")
        self.genomes = genomes

    @property
    def files(self):
        if self.kind == "taxonomy":
            return [self.tree, self.genomes]
        return [self.tree, self.snp, self.references]

    def path(self, relative):
        return os.path.join(REPO, relative)

    def __repr__(self):
        return f"<Scheme {self.name} ({self.kind})>"


def load_config():
    if not os.path.exists(CONFIG):
        return {}
    with open(CONFIG, "rb") as handle:
        return tomllib.load(handle)


def config_keys():
    """The keys a schemes.toml entry may use: the Scheme constructor's arguments."""
    return [name for name in inspect.signature(Scheme.__init__).parameters
            if name not in ("self", "name")]


def check_config(config):
    """Return a problem for every entry or key in schemes.toml that is not understood.

    An unknown key used to reach Scheme() and stop the script with a TypeError,
    and a misspelt optional key would otherwise be ignored without notice.
    """
    problems = []
    allowed = config_keys()
    for name, entry in sorted(config.items()):
        if not isinstance(entry, dict):
            problems.append(f"schemes.toml: '{name}' must be a table ([{name}])")
            continue
        for key in sorted(set(entry) - set(allowed)):
            problems.append(f"schemes.toml: [{name}] has an unknown key '{key}'; "
                            f"recognised keys: {', '.join(allowed)}")
    return problems


def discover(only=None):
    """Return the declared schemes, discovering SNP schemes by convention."""
    config = load_config()
    config_problems = check_config(config)
    if config_problems:
        for problem in config_problems:
            print(f"error: {problem}", file=sys.stderr)
        raise SystemExit(2)
    names = set()
    for entry in sorted(os.listdir(SOURCE_DIR)):
        if entry.endswith("_snp.txt"):
            names.add(entry[: -len("_snp.txt")])
    names.update(config)

    schemes = []
    for name in sorted(names):
        schemes.append(Scheme(name, **config.get(name, {})))

    missing = [s for s in schemes if any(not os.path.exists(s.path(f)) for f in s.files)]
    if missing:
        for scheme in missing:
            absent = [f for f in scheme.files if not os.path.exists(scheme.path(f))]
            print(f"error: scheme '{scheme.name}' declares missing files: "
                  f"{', '.join(absent)}", file=sys.stderr)
        raise SystemExit(2)

    if only:
        known = {s.name for s in schemes}
        unknown = set(only) - known
        if unknown:
            print(f"error: no such scheme: {', '.join(sorted(unknown))}. "
                  f"Known schemes: {', '.join(sorted(known))}", file=sys.stderr)
            raise SystemExit(2)
        schemes = [s for s in schemes if s.name in only]
    return schemes


def read_rows(path, problems=None, relative=None):
    """Return the file's lines split on tabs, with the raw bytes for format checks.

    A file that is not valid UTF-8 is reported to problems, when given, and read
    with the undecodable bytes replaced, so the other checks still run.
    """
    with open(path, "rb") as handle:
        raw = handle.read()
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError as error:
        if problems is None:
            raise
        problems.append(f"{relative or path}: not valid UTF-8 ({error.reason} at "
                        f"byte {error.start})")
        text = raw.decode("utf-8", errors="replace")
    lines = text.split("\n")
    if lines and lines[-1] == "":
        lines.pop()
    return raw, lines, [line.split("\t") for line in lines]


def check_format(problems, relative, raw, lines):
    """Checks that apply to every scheme source file regardless of its columns."""
    if raw.startswith(codecs.BOM_UTF8):
        problems.append(f"{relative}: starts with a UTF-8 byte order mark, which "
                        f"becomes part of the first field")
    if b"\r" in raw:
        count = raw.count(b"\r")
        problems.append(f"{relative}: contains {count} carriage return(s); "
                        f"files must use LF line endings")
    if raw and not raw.endswith(b"\n"):
        problems.append(f"{relative}: missing a final newline")
    for number, line in enumerate(lines, start=1):
        if line != line.rstrip("\t "):
            problems.append(f"{relative}:{number}: trailing whitespace")
        if not line.strip():
            problems.append(f"{relative}:{number}: blank line")
        if ";" in line:
            problems.append(f"{relative}:{number}: contains a semicolon, which the "
                            f"tree reader treats as a column separator")


def check_header(problems, relative, rows, expected):
    if not rows:
        problems.append(f"{relative}: file is empty")
        return False
    if rows[0] != expected:
        problems.append(f"{relative}:1: unexpected header\n"
                        f"      expected: {chr(9).join(expected)}\n"
                        f"      found:    {chr(9).join(rows[0])}")
        return False
    return True


def tree_nodes_and_parents(rows):
    """Derive the node set and each node's parents from the cumulative paths."""
    nodes = set()
    parents = {}
    for row in rows:
        fields = [f for f in row if f.strip()]
        nodes.update(fields)
        for index in range(1, len(fields)):
            parents.setdefault(fields[index], set()).add(fields[index - 1])
    return nodes, parents


def check_tree(problems, scheme):
    relative = scheme.tree
    raw, lines, rows = read_rows(scheme.path(relative), problems, relative)
    check_format(problems, relative, raw, lines)
    if not rows:
        problems.append(f"{relative}: file is empty")
        return set()

    for number, row in enumerate(rows, start=1):
        for index, field in enumerate(row):
            if not field.strip() and index < len(row) - 1:
                problems.append(f"{relative}:{number}: empty column {index + 1} "
                                f"between populated columns")

    roots = {row[0] for row in rows if row and row[0].strip()}
    if len(roots) > 1:
        problems.append(f"{relative}: every line must start at the same root node, "
                        f"found {len(roots)}: {', '.join(sorted(roots))}")

    nodes, parents = tree_nodes_and_parents(rows)
    for node, node_parents in sorted(parents.items()):
        if len(node_parents) > 1:
            problems.append(f"{relative}: node '{node}' is given more than one "
                            f"parent: {', '.join(sorted(node_parents))}")
    return nodes


def check_references(problems, scheme):
    relative = scheme.references
    raw, lines, rows = read_rows(scheme.path(relative), problems, relative)
    check_format(problems, relative, raw, lines)
    if not check_header(problems, relative, rows, REFERENCE_HEADER):
        return set()

    genomes = set()
    for number, row in enumerate(rows[1:], start=2):
        if len(row) != len(REFERENCE_HEADER):
            problems.append(f"{relative}:{number}: expected "
                            f"{len(REFERENCE_HEADER)} columns, found {len(row)}")
            continue
        genome, _strain, genbank, refseq, _assembly = row
        genomes.add(genome)
        # CanSNPer2 recovers the reference name from "<genome>_<query>.xmfa" by
        # splitting on the first underscore.
        if "_" in genome:
            problems.append(f"{relative}:{number}: genome '{genome}' contains an "
                            f"underscore, which CanSNPer2 cannot use in a "
                            f"reference name")
        if not GENBANK_RE.match(genbank):
            problems.append(f"{relative}:{number}: genbank_id '{genbank}' is not of "
                            f"the form GCA_000000000.0")
        if not REFSEQ_RE.match(refseq):
            problems.append(f"{relative}:{number}: refseq_id '{refseq}' is not of "
                            f"the form GCF_000000000.0")

    duplicates = [g for g, n in Counter(r[0] for r in rows[1:]).items() if n > 1]
    for genome in sorted(duplicates):
        problems.append(f"{relative}: genome '{genome}' is listed more than once")
    return genomes


def check_snp(problems, warnings, scheme, tree_nodes, reference_genomes):
    relative = scheme.snp
    raw, lines, rows = read_rows(scheme.path(relative), problems, relative)
    check_format(problems, relative, raw, lines)
    if not check_header(problems, relative, rows, SNP_HEADER):
        return

    snp_ids = []
    positions = Counter()
    base_issues = set()
    for number, row in enumerate(rows[1:], start=2):
        if len(row) != len(SNP_HEADER):
            problems.append(f"{relative}:{number}: expected {len(SNP_HEADER)} "
                            f"columns, found {len(row)}")
            continue
        snp_id, _strain, _reference, genome, position, derived, ancestral = row
        snp_ids.append(snp_id)

        if not POSITION_RE.match(position) or int(position) < 1:
            problems.append(f"{relative}:{number}: position '{position}' is not a "
                            f"positive integer")
        else:
            # Keyed on the number, so that 0100 and 100 count as the same position.
            positions[(genome, int(position))] += 1
        # A SNP listed in known_issues is already published with a malformed
        # base column. Report it every run, but do not fail the build on it.
        base_report = warnings if snp_id in scheme.known_issues else problems
        base_problems = []
        if derived not in BASES:
            base_problems.append(f"{relative}:{number}: derived_base '{derived}' is "
                                 f"not one of A, C, G, T")
        if ancestral not in BASES:
            base_problems.append(f"{relative}:{number}: ancestral_base '{ancestral}' "
                                 f"is not one of A, C, G, T")
        if derived == ancestral and derived in BASES:
            base_problems.append(f"{relative}:{number}: derived and ancestral base "
                                 f"are both '{derived}'")
        if base_problems:
            base_issues.add(snp_id)
            base_report.extend(base_problems)
        if genome not in reference_genomes:
            problems.append(f"{relative}:{number}: genome '{genome}' is not listed "
                            f"in {scheme.references}")

    # A known_issues entry that names no SNP, or a SNP whose bases are now
    # valid, is stale: it would exempt a later error in that SNP from failing.
    for snp_id in sorted(scheme.known_issues - set(snp_ids)):
        problems.append(f"schemes.toml: known_issues for '{scheme.name}' lists "
                        f"'{snp_id}', which is not in {relative}")
    for snp_id in sorted((scheme.known_issues & set(snp_ids)) - base_issues):
        problems.append(f"schemes.toml: known_issues for '{scheme.name}' lists "
                        f"'{snp_id}', whose bases are valid; remove the entry")

    for snp_id, count in sorted(Counter(snp_ids).items()):
        if count > 1:
            problems.append(f"{relative}: snp_id '{snp_id}' appears {count} times")
    for (genome, position), count in sorted(positions.items()):
        if count > 1:
            problems.append(f"{relative}: position {position} on '{genome}' is "
                            f"annotated {count} times")

    annotated = set(snp_ids)
    for snp_id in sorted(annotated - tree_nodes):
        problems.append(f"{relative}: snp_id '{snp_id}' does not appear in "
                        f"{scheme.tree}")
    for node in sorted(tree_nodes - annotated):
        problems.append(f"{scheme.tree}: node '{node}' has no SNP in {scheme.snp}")


def check_genomes(problems, scheme, tree_nodes):
    relative = scheme.genomes
    raw, lines, rows = read_rows(scheme.path(relative), problems, relative)
    check_format(problems, relative, raw, lines)

    seen = Counter()
    for number, row in enumerate(rows, start=1):
        if len(row) != 2:
            problems.append(f"{relative}:{number}: expected 2 columns "
                            f"(genome, node), found {len(row)}")
            continue
        genome, node = row
        seen[genome] += 1
        if node not in tree_nodes:
            problems.append(f"{relative}:{number}: node '{node}' does not appear "
                            f"in {scheme.tree}")
    for genome, count in sorted(seen.items()):
        if count > 1:
            problems.append(f"{relative}: genome '{genome}' is listed {count} times")


def validate(scheme):
    """Return the problems and the non-fatal warnings found in one scheme."""
    problems, warnings = [], []
    tree_nodes = check_tree(problems, scheme)
    if scheme.kind == "taxonomy":
        check_genomes(problems, scheme, tree_nodes)
    else:
        reference_genomes = check_references(problems, scheme)
        check_snp(problems, warnings, scheme, tree_nodes, reference_genomes)
    return problems, warnings


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--scheme", action="append", metavar="NAME",
                        help="validate only this scheme; repeatable")
    args = parser.parse_args()

    schemes = discover(args.scheme)
    failed = 0
    warned = 0
    for scheme in schemes:
        problems, warnings = validate(scheme)
        label = f"{scheme.name} (work in progress)" if scheme.wip else scheme.name
        if problems:
            failed += 1
            print(f"FAIL {label} ({len(problems)} problem"
                  f"{'s' if len(problems) != 1 else ''})")
        else:
            print(f"ok   {label}")
        for problem in problems:
            print(f"  - {problem}")
        if warnings:
            warned += len(warnings)
            for warning in warnings:
                print(f"  ! known issue: {warning}")

    if warned:
        print(f"\n{warned} known issue(s) reported; see known_issues in schemes.toml.")
    if failed:
        print(f"{failed} of {len(schemes)} schemes failed validation.")
        return 1
    print(f"All {len(schemes)} schemes valid.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
