import os
import sys

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)

from validate import Scheme  # noqa: E402

TREE = "R\nR\tA\nR\tA\tB\nR\tC\n"
REFERENCES = ("genome\tstrain\tgenbank_id\trefseq_id\tassembly_name\n"
              "G\tstrain G\tGCA_000000001.1\tGCF_000000001.1\tASM1v1\n")
SNP_HEADER = "snp_id\tstrain\treference\tgenome\tposition\tderived_base\tancestral_base\n"
SNP_ROWS = [
    "R\ts\tTest2026\tG\t100\tA\tG",
    "A\ts\tTest2026\tG\t200\tT\tC",
    "B\ts\tTest2026\tG\t300\tC\tT",
    "C\ts\tTest2026\tG\t400\tG\tA",
]


@pytest.fixture
def make_scheme(tmp_path):
    """Write a small SNP scheme (R -> A -> B, R -> C) and return a Scheme for it.

    Paths are absolute, so the scheme lives outside the repository; the
    database and manifest go to tmp_path/database.
    """
    def make(tree=TREE, snp_rows=None, references=REFERENCES, name="test",
             raw_snp=None, genomes=None, **options):
        files = {
            "tree": tmp_path / f"{name}_tree.txt",
            "snp": tmp_path / f"{name}_snp.txt",
            "references": tmp_path / f"{name}_references.txt",
        }
        files["tree"].write_text(tree)
        files["references"].write_text(references)
        if raw_snp is not None:
            files["snp"].write_bytes(raw_snp)
        else:
            rows = SNP_ROWS if snp_rows is None else snp_rows
            files["snp"].write_text(SNP_HEADER + "".join(row + "\n" for row in rows))
        genomes_path = tmp_path / f"{name}_genomes.txt"
        if genomes is not None:
            genomes_path.write_text(genomes)
        database_dir = tmp_path / "database"
        database_dir.mkdir(exist_ok=True)
        return Scheme(name, tree=str(files["tree"]), snp=str(files["snp"]),
                      references=str(files["references"]), genomes=str(genomes_path),
                      database_dir=str(database_dir), **options)
    return make
