#!/usr/bin/env python3
"""Validate the canSNP scheme sources in this repository.

Checks the text files under schemes/ for format problems and for the cross-file
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
# Each SNP scheme has a folder schemes/<name>/ holding <name>_tree.txt and the
# other files of the scheme.
SCHEMES_DIR = os.path.join(REPO, "schemes")
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
                 references=None, snp_set=None, verify=True, source_type="CanSNPer",
                 known_issues=None, wip=False, database_dir="database", readme=None):
        self.name = name
        self.kind = kind
        self.verify = verify
        self.source_type = source_type
        self.known_issues = set(known_issues or ())
        self.wip = wip
        self.database_dir = database_dir
        folder = os.path.join("schemes", name)
        self.tree = tree or os.path.join(folder, f"{name}_tree.txt")
        self.snp = snp or os.path.join(folder, f"{name}_snp.txt")
        self.references = references or os.path.join(folder, f"{name}_references.txt")
        # A taxonomy scheme needs its genome-to-node file. An SNP scheme may
        # list the genomes it was made from, which is optional and checked when
        # present.
        if genomes is None and kind == "snp":
            genomes = os.path.join(folder, f"{name}_genomes.txt")
        self.genomes = genomes
        # An SNP scheme may also list more than one canSNP per node, ranked,
        # which is optional and checked when present.
        if snp_set is None and kind == "snp":
            snp_set = os.path.join(folder, f"{name}_canSNPs.txt")
        self.snp_set = snp_set
        # An SNP scheme is described by a README.md in its folder, in the
        # layout given in schemes/README.md.
        if readme is None and kind == "snp":
            readme = os.path.join(folder, "README.md")
        self.readme = readme

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
    names, misplaced = set(), []
    for root, _dirs, files in os.walk(SCHEMES_DIR):
        for entry in sorted(files):
            if not entry.endswith("_snp.txt"):
                continue
            name = entry[: -len("_snp.txt")]
            if os.path.normpath(root) == os.path.join(os.path.normpath(SCHEMES_DIR), name):
                names.add(name)
            else:
                misplaced.append(os.path.relpath(os.path.join(root, entry), REPO))
    if misplaced:
        # A file outside its own folder would otherwise not be checked or built.
        for path in sorted(misplaced):
            print(f"error: {path} is not in the folder of its scheme; an SNP "
                  f"scheme <name> belongs in schemes/<name>/", file=sys.stderr)
        raise SystemExit(2)
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


SNP_GENOME_HEADER = ["genome", "node", "role"]
SNP_GENOME_ROLES = {
    "representative": "its sequence was in the alignment the canSNPs were found in",
    "member": "placed at the node through a representative, as dereplication folds it in",
    "outgroup": "aligned, but placed at no node",
}


def check_snp_genomes(problems, scheme, tree_nodes):
    """The genomes an SNP scheme was made from, when the scheme lists them.

    Each row names an assembly by its versioned accession, so that the genomes
    can be fetched again, the node it sits at, and its role. An outgroup sits
    at no node, so its node column is empty.
    """
    relative = scheme.genomes
    if not relative or not os.path.exists(scheme.path(relative)):
        return
    raw, lines, rows = read_rows(scheme.path(relative), problems, relative)
    check_format(problems, relative, raw, lines)
    if not check_header(problems, relative, rows, SNP_GENOME_HEADER):
        return

    assemblies = {}
    for number, row in enumerate(rows[1:], start=2):
        if len(row) != 3:
            problems.append(f"{relative}:{number}: expected 3 columns "
                            f"(genome, node, role), found {len(row)}")
            continue
        genome, node, role = row
        if not (GENBANK_RE.match(genome) or REFSEQ_RE.match(genome)):
            problems.append(f"{relative}:{number}: '{genome}' is not a versioned assembly "
                            f"accession such as GCA_000007125.1")
            continue
        if role not in SNP_GENOME_ROLES:
            problems.append(f"{relative}:{number}: role '{role}' is not one of "
                            f"{', '.join(SNP_GENOME_ROLES)}")
        elif role == "outgroup" and node:
            problems.append(f"{relative}:{number}: an outgroup sits at no node, so its "
                            f"node column must be empty")
        elif role != "outgroup" and node not in tree_nodes:
            problems.append(f"{relative}:{number}: node '{node}' does not appear in "
                            f"{scheme.tree}")
        # GCA_ and GCF_ accessions with the same digits are one assembly.
        assembly = genome[4:].split(".")[0]
        assemblies.setdefault(assembly, []).append((number, genome))
    for listed in assemblies.values():
        if len(listed) > 1:
            places = ", ".join(f"{genome} on line {number}" for number, genome in listed)
            problems.append(f"{relative}: one assembly is listed {len(listed)} times: {places}")


SNP_SET_HEADER = ["node", "rank", "genome", "position", "derived_base", "ancestral_base",
                  "score"]
SCORE_RE = re.compile(r"^-?[0-9]+(\.[0-9]+)?$")


def check_snp_set(problems, scheme, tree_nodes, reference_genomes):
    """The canSNPs of each node, ranked, when the scheme lists more than one.

    CanSNPer2 calls one SNP per node, from the SNP file, so rank 1 of each node
    must be the SNP the SNP file gives it. The other ranks are further canSNPs
    of the same node, for assays that type several, and every position may
    mark only one node.
    """
    relative = scheme.snp_set
    if not relative or not os.path.exists(scheme.path(relative)):
        return
    raw, lines, rows = read_rows(scheme.path(relative), problems, relative)
    check_format(problems, relative, raw, lines)
    if not check_header(problems, relative, rows, SNP_SET_HEADER):
        return

    ranks = {}
    first = {}
    positions = Counter()
    for number, row in enumerate(rows[1:], start=2):
        if len(row) != len(SNP_SET_HEADER):
            problems.append(f"{relative}:{number}: expected {len(SNP_SET_HEADER)} "
                            f"columns, found {len(row)}")
            continue
        node, rank, genome, position, derived, ancestral, score = row
        if node not in tree_nodes:
            problems.append(f"{relative}:{number}: node '{node}' does not appear in "
                            f"{scheme.tree}")
        if not POSITION_RE.match(rank) or int(rank) < 1:
            problems.append(f"{relative}:{number}: rank '{rank}' is not a positive integer")
        else:
            ranks.setdefault(node, []).append(int(rank))
        if genome not in reference_genomes:
            problems.append(f"{relative}:{number}: genome '{genome}' is not listed in "
                            f"{scheme.references}")
        if not POSITION_RE.match(position) or int(position) < 1:
            problems.append(f"{relative}:{number}: position '{position}' is not a "
                            f"positive integer")
            continue
        positions[(genome, int(position))] += 1
        if derived not in BASES or ancestral not in BASES:
            problems.append(f"{relative}:{number}: the bases '{derived}' and '{ancestral}' "
                            f"are not both one of A, C, G, T")
        elif derived == ancestral:
            problems.append(f"{relative}:{number}: derived and ancestral base are both "
                            f"'{derived}'")
        if not SCORE_RE.match(score):
            problems.append(f"{relative}:{number}: score '{score}' is not a number")
        if rank == "1":
            first[node] = (genome, int(position), derived, ancestral)

    for node, listed in sorted(ranks.items()):
        if listed != list(range(1, len(listed) + 1)):
            problems.append(f"{relative}: the ranks of node '{node}' are "
                            f"{', '.join(map(str, listed))}, not 1, 2, 3 ... in order")
    for (genome, position), count in sorted(positions.items()):
        if count > 1:
            problems.append(f"{relative}: position {position} on '{genome}' is listed "
                            f"{count} times")

    # Rank 1 must agree with the SNP file, which is what CanSNPer2 calls.
    # Undecodable bytes in the SNP file are reported by check_snp already.
    _raw, _lines, snp_rows = read_rows(scheme.path(scheme.snp), [], scheme.snp)
    called = {}
    for row in snp_rows[1:]:
        if len(row) == len(SNP_HEADER) and POSITION_RE.match(row[4]):
            called[row[0]] = (row[3], int(row[4]), row[5], row[6])
    for node in sorted(set(called) - set(first)):
        problems.append(f"{relative}: node '{node}' has no canSNP of rank 1")
    for node in sorted(set(called) & set(first)):
        if first[node] != called[node]:
            problems.append(f"{relative}: rank 1 of node '{node}' is not the SNP "
                            f"{scheme.snp} gives it")


# The README of an SNP scheme: a title, a table of these facts in this order,
# and these sections in this order. schemes/README.md describes them.
README_FACTS = ["Scheme", "Status", "Organisms", "Reference genome", "Nodes",
                "canSNPs per node", "Origin", "Last changed"]
README_SECTIONS = ["Description", "Origin", "Validation", "Known issues",
                   "References", "Changes"]
STATUSES = ("draft", "verified")
ACCESSION_RE = re.compile(r"GC[AF]_\d{9}\.\d{1,2}")
DATE_RE = re.compile(r"^\d{4}-\d{2}-\d{2}$")


def parse_readme(text):
    """Split a scheme README into its title, its facts table and its sections.

    Returns (title, facts, sections): facts is a list of (key, value) pairs in
    the order of the table, and sections a list of (heading, body) pairs for
    the second-level headings.
    """
    title, facts, sections = None, [], []
    lines = text.splitlines()
    index = 0
    while index < len(lines) and not lines[index].strip():
        index += 1
    if index < len(lines) and lines[index].startswith("# "):
        title = lines[index][2:].strip()
        index += 1
    # The facts table is the first table after the title.
    while index < len(lines) and not lines[index].startswith("|") \
            and not lines[index].startswith("## "):
        index += 1
    rows = []
    while index < len(lines) and lines[index].startswith("|"):
        rows.append([cell.strip() for cell in lines[index].strip().strip("|").split("|")])
        index += 1
    # Skip the empty header row and the separator row.
    for cells in rows[2:]:
        if len(cells) == 2:
            facts.append((cells[0], cells[1]))
    current = None
    for line in lines[index:]:
        if line.startswith("## "):
            current = [line[3:].strip(), []]
            sections.append(current)
        elif current is not None:
            current[1].append(line)
    return title, facts, [(heading, "\n".join(body).strip()) for heading, body in sections]


def readme_status(scheme):
    """The Status given in a scheme's README, or None."""
    if not scheme.readme or not os.path.exists(scheme.path(scheme.readme)):
        return None
    with open(scheme.path(scheme.readme), encoding="utf-8") as handle:
        return dict(parse_readme(handle.read())[1]).get("Status")


def check_readme(problems, scheme, tree_nodes, require_verified=False):
    """Check the README of an SNP scheme against the layout and the scheme files."""
    relative = scheme.readme
    path = scheme.path(relative)
    if not os.path.exists(path):
        problems.append(f"{scheme.name} has no README.md; expected {relative}, "
                        f"in the layout described in schemes/README.md")
        return
    with open(path, encoding="utf-8") as handle:
        title, facts, sections = parse_readme(handle.read())
    if not title:
        problems.append(f"{relative}: must start with a '# ' title")
    keys = [key for key, _value in facts]
    if keys != README_FACTS:
        problems.append(f"{relative}: the facts table must have the rows "
                        f"{', '.join(README_FACTS)}, in that order; found "
                        f"{', '.join(keys) or 'none'}")
    values = dict(facts)

    scheme_value = values.get("Scheme")
    if scheme_value is not None and scheme_value != f"`{scheme.name}`":
        problems.append(f"{relative}: Scheme is {scheme_value}, expected `{scheme.name}`")
    status = values.get("Status")
    if status is not None and status not in STATUSES:
        problems.append(f"{relative}: Status is '{status}'; it must be draft or verified")
    if require_verified and status == "draft":
        problems.append(f"{relative}: {scheme.name} is a draft; only verified "
                        f"schemes belong on master")
    nodes = values.get("Nodes")
    if nodes is not None and nodes != str(len(tree_nodes)):
        problems.append(f"{relative}: Nodes is {nodes}, but the tree file has "
                        f"{len(tree_nodes)} nodes")
    reference = values.get("Reference genome")
    if reference is not None:
        accessions = set()
        for row in read_rows(scheme.path(scheme.references))[2][1:]:
            accessions.update(row[2:4])
        named = set(ACCESSION_RE.findall(reference))
        if not named & accessions:
            problems.append(f"{relative}: Reference genome names no accession from "
                            f"{scheme.references}")
    per_node = values.get("canSNPs per node")
    if per_node is not None:
        expected = "1"
        if scheme.snp_set and os.path.exists(scheme.path(scheme.snp_set)):
            ranks = [int(row[1]) for row in read_rows(scheme.path(scheme.snp_set))[2][1:]
                     if len(row) > 1 and row[1].isdigit()]
            if ranks and max(ranks) > 1:
                expected = f"up to {max(ranks)}"
        if per_node != expected:
            problems.append(f"{relative}: canSNPs per node is '{per_node}', "
                            f"expected '{expected}' from the scheme files")
    changed = values.get("Last changed")
    if changed is not None and not DATE_RE.match(changed):
        problems.append(f"{relative}: Last changed is '{changed}'; use YYYY-MM-DD")

    headings = [heading for heading, _body in sections]
    if headings != README_SECTIONS:
        problems.append(f"{relative}: the sections must be "
                        f"{', '.join(README_SECTIONS)}, in that order; found "
                        f"{', '.join(headings) or 'none'}")
    for heading, body in sections:
        if heading in README_SECTIONS and not body:
            problems.append(f"{relative}: section '{heading}' is empty; write "
                            f"'None known.' or similar where there is nothing to report")


def validate(scheme, require_verified=False):
    """Return the problems and the non-fatal warnings found in one scheme.

    With require_verified, an SNP scheme whose README gives its Status as draft
    is a problem, as on master, which holds verified schemes only.
    """
    problems, warnings = [], []
    tree_nodes = check_tree(problems, scheme)
    if scheme.kind == "taxonomy":
        check_genomes(problems, scheme, tree_nodes)
    else:
        reference_genomes = check_references(problems, scheme)
        check_snp(problems, warnings, scheme, tree_nodes, reference_genomes)
        check_snp_genomes(problems, scheme, tree_nodes)
        check_snp_set(problems, scheme, tree_nodes, reference_genomes)
        if not scheme.wip:
            check_readme(problems, scheme, tree_nodes, require_verified)
    return problems, warnings


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--scheme", action="append", metavar="NAME",
                        help="validate only this scheme; repeatable")
    parser.add_argument("--require-verified", action="store_true",
                        help="fail on a scheme whose README gives its Status as "
                             "draft; used for master, which holds verified schemes")
    args = parser.parse_args()

    schemes = discover(args.scheme)
    failed = 0
    warned = 0
    for scheme in schemes:
        problems, warnings = validate(scheme, args.require_verified)
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
