import os
import sys

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)

from validate import Scheme  # noqa: E402

TREE = "R\nR\tA\nR\tA\tB\nR\tC\n"
REFERENCES = ("genome\tstrain\tgenbank_id\trefseq_id\tassembly_name\n"
              "G\tstrain G\tGCA_000000001.1\tGCF_000000001.1\tASM1v1\n")
CANSNP_HEADER = ("node\trank\treference\tgenome\tsequence\tsequence_position\t"
                 "derived_base\tancestral_base\tscore\ttyped\n")
CANSNP_ROWS = [
    "R\t1\tTest2026\tG\tseq1\t100\tA\tG\t520\tyes",
    "A\t1\tTest2026\tG\tseq1\t200\tT\tC\t500\tyes",
    "B\t1\tTest2026\tG\tseq2\t50\tC\tT\t495\tyes",
    "C\t1\tTest2026\tG\tseq2\t150\tG\tA\t470\tyes",
]


def readme_text(name="test", status="verified", nodes="4", reference="strain G (GCA_000000001.1)",
                per_node="1", changed="2026-11-19", sections=None):
    """A scheme README in the required layout, for the small scheme below."""
    if sections is None:
        sections = ["Description", "Origin", "Validation", "Known issues",
                    "References", "Changes"]
    facts = [("Scheme", f"`{name}`"), ("Status", status), ("Organisms", "*Testus testi*"),
             ("Reference genome", reference), ("Nodes", nodes),
             ("canSNPs per node", per_node), ("Origin", "test data"),
             ("Last changed", changed)]
    lines = ["# Test canSNP scheme", "", "| | |", "|---|---|"]
    lines += [f"| {key} | {value} |" for key, value in facts]
    for section in sections:
        lines += ["", f"## {section}", "", f"Text about {section.lower()}."]
    return "\n".join(lines) + "\n"


@pytest.fixture
def make_scheme(tmp_path):
    """Write a small SNP scheme (R -> A -> B, R -> C) and return a Scheme for it.

    Paths are absolute, so the scheme lives outside the repository; the
    database and manifest go to tmp_path/database. Several canSNPs of a node
    are given as rows with ranks 2 and up in cansnp_rows.
    """
    def make(tree=TREE, cansnp_rows=None, references=REFERENCES, name="test",
             raw_cansnps=None, genomes=None, readme=None, **options):
        files = {
            "tree": tmp_path / f"{name}_tree.txt",
            "cansnps": tmp_path / f"{name}_canSNPs.txt",
            "references": tmp_path / f"{name}_references.txt",
        }
        files["tree"].write_text(tree)
        files["references"].write_text(references)
        if raw_cansnps is not None:
            files["cansnps"].write_bytes(raw_cansnps)
        else:
            rows = CANSNP_ROWS if cansnp_rows is None else cansnp_rows
            files["cansnps"].write_text(CANSNP_HEADER + "".join(row + "\n" for row in rows))
        genomes_path = tmp_path / f"{name}_genomes.txt"
        if genomes is not None:
            genomes_path.write_text(genomes)
        # A conforming README unless the test gives its own text, or False for none.
        readme_path = tmp_path / "README.md"
        if readme is None:
            readme = readme_text(name=name)
        if readme is not False:
            readme_path.write_text(readme)
        database_dir = tmp_path / "database"
        database_dir.mkdir(exist_ok=True)
        return Scheme(name, tree=str(files["tree"]), cansnps=str(files["cansnps"]),
                      references=str(files["references"]), genomes=str(genomes_path),
                      readme=str(readme_path), database_dir=str(database_dir), **options)
    return make
