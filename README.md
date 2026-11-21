# CanSNPer2-data

The canSNP schemes and databases used by CanSNPer2
(https://github.com/FOI-Bioinformatics/CanSNPer2) to type bacterial genomes.
Each scheme is kept as tab-separated text under `schemes/`, with a README that
describes where it comes from, how it was checked and its known issues. The
SQLite databases under `database/` are built from that text.

## Branches

| Branch | Contents |
| ------ | -------- |
| `master` | Verified schemes. |
| `dev` | The verified schemes and draft schemes that have not yet been verified. New schemes are added here by pull request and moved to `master` once verified. |

## Repository layout

| Path | Contents |
| ---- | -------- |
| `schemes/<name>/` | One folder per canSNP scheme, holding its files as tab-separated text and a `README.md` describing it. This is what you edit. See `schemes/README.md` for the formats. |
| `database/` | The SQLite databases generated from `schemes/`, with a provenance manifest beside each one. |
| `database/archive/` | The database each rebuild replaced. |
| `schemes.toml` | Declares schemes that are not plain SNP schemes, and records known data issues. |

## Changing a scheme

Edit the text files under `schemes/<name>/`, then check and rebuild:

```sh
python validate.py          # check the scheme files
python check_bases.py       # check the SNP bases against the reference genomes
python build.py             # rebuild and compare, changing nothing
python build.py --write     # install the rebuilt databases and manifests
```

Commit the changed scheme files together with the database and manifest that
`--write` produced. Continuous integration runs the first three commands on every
pull request, so a scheme file that breaks an invariant, a SNP whose reference
genome has neither its derived nor its ancestral base, or a database that was not
regenerated after its sources changed, fails the build.

`check_bases.py` downloads the GenBank assembly of every reference genome from
NCBI into `.references/` (once; `--cache` names another folder) and checks the
base at every position of `<name>_snp.txt` and `<name>_canSNPs.txt`. Each SNP
names a sequence of the assembly (a chromosome or plasmid, by its GenBank
accession) and a position on it, so the order of the sequences in the assembly
file does not matter. A sequence that the assembly lacks is reported as a
problem. It prints the sequences of each reference with the number of SNPs on
each.

Each scheme folder has a `README.md` in a fixed layout, with a Status of `draft`
or `verified`; `schemes/TEMPLATE.md` is a copy to start from. When a scheme
changes, update its README as well, including Last changed and Changes. On pull
requests into `master`, continuous integration also runs
`python validate.py --require-verified`, which fails on a draft scheme.

`validate.py` needs only Python 3.11 or later. `build.py` also needs CanSNPer2
2.3.0 or later; until 2.3.0 is released it is installed from the `dev_2.3.0`
branch. Building a taxonomy scheme also needs flextaxd; SNP schemes do not:

```sh
pip install git+https://github.com/FOI-Bioinformatics/CanSNPer2.git@dev_2.3.0
```

Since CanSNPer2 2.3.0 the scheme files and the databases give each SNP as a
reference sequence (the GenBank accession of a chromosome or plasmid) and a
1-based position on that sequence. Earlier versions counted positions along the
reference sequences joined in the order of the assembly file. CanSNPer2 2.3.0
refuses scheme files and databases in that format; `CanSNPer2-database
--convert` converts them (see `schemes/README.md`). The databases built in the
earlier format, which CanSNPer2 2.2 reads, are kept under `database/archive/`.

No aligner is needed to build a database. Typing a sample needs one of the
aligners CanSNPer2 supports: progressiveMauve, minimap2 or MUMmer4.

The scripts have their own tests, which CI also runs:

```sh
python -m pytest tests
```

## Reserving node names

Node names for a publication can be reserved before the canSNPs and the tree
are made public: contact the maintainers. The README of a scheme gives the next
free node name where the scheme numbers its nodes.
