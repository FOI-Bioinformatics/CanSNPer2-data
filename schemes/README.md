# Schemes

Each scheme in this repository has its own folder here, holding the scheme as
tab-separated text and a `README.md` describing it. The SQLite databases under
`database/` are built from these files and are not edited directly.

There are two kinds of scheme. An SNP scheme is a tree of canSNPs that
CanSNPer2 types genomes with. A taxonomy scheme is a tree with genomes assigned
to its nodes and no canSNPs; CanSNPer2 cannot type with it.

All files use LF line endings, carry no trailing whitespace and end with a
newline. Carriage returns matter: a stray `\r` can end up inside a node name.
`.gitattributes` and `validate.py` both enforce this.

A semicolon may not appear anywhere in these files. The CanSNPer tree reader
converts tabs to semicolons before splitting, so a semicolon inside a name is
read as a column break.

## SNP schemes

A scheme is named in lower case with underscores, after the organisms it
covers (`yersinia_pestis`, `brucella_classical`). Its folder,
`schemes/<name>/`, holds four files that are required:

- `<name>_tree.txt`, the tree;
- `<name>_snp.txt`, one canSNP per node;
- `<name>_references.txt`, the reference genomes the positions are given on;
- `README.md`, describing the scheme;

and two that are optional:

- `<name>_genomes.txt`, the genomes the scheme was made from;
- `<name>_canSNPs.txt`, more than one canSNP per node.

Adding the folder with the required files is all that is needed to declare an
SNP scheme; `validate.py` and `build.py` discover it by name. A `_snp.txt` file
outside the folder of its name is reported as an error.

Node names may contain letters, digits, `.`, `/`, `-`, `_` and spaces, but no
semicolon and no tab. A node name is unique within a scheme.

## `<name>_tree.txt`

One line per path from the root of the scheme to a node, with the nodes of the
path in tab-separated columns:

```
A/B.Br.001
A/B.Br.001	B.Br.003
A/B.Br.001	B.Br.003	B.Br.004
```

Every line starts at the same root node. A node may be named on many lines but
must always appear under the same parent. The tree reader does not reject a
node placed under two parents: it creates a second node with the same name, and
the SNP is attached to only one of them. `validate.py` rejects such a file, and
`build.py` rejects a database with a repeated node name.

## `<name>_snp.txt`

One line per canSNP, with this header:

```
snp_id	strain	reference	genome	sequence	sequence_position	derived_base	ancestral_base
```

`snp_id` must name a node in the tree file, and every node in the tree file must
have exactly one SNP here. `genome` must name a row in the references file.
`sequence` is the GenBank accession of the reference sequence (chromosome,
plasmid or contig) that holds the SNP, as it is named in the assembly file of
that genome, for example `AE008917.1` for chromosome I of Brucella melitensis
16M. `sequence_position` is a positive integer, one-based, on that sequence.
Each sequence of a reference genome thus has its own coordinates, and the order
of the sequences in the assembly file does not matter.
`derived_base` and `ancestral_base` are single characters from A, C, G and T and
must differ from each other. No two SNPs may annotate the same position on the
same sequence of a reference genome.

Before CanSNPer2 2.3.0 the file had a single `position` column, counted along
the sequences of the reference genome joined in the order of its assembly file.
CanSNPer2 2.3.0 and `validate.py` refuse that format. A file in it can be
converted with

```
CanSNPer2-database --convert --annotation <name>_snp.txt \
    --references <name>_references.txt --refdir <folder> -o <outdir>
```

which checks each SNP against its reference genome (`<folder>/<genome>.fna`,
downloaded from NCBI when missing) and writes the converted file to `<outdir>`.
The same command converts a `<name>_canSNPs.txt` file.

A SNP whose bases are published in another form can be listed under
`known_issues` in `schemes.toml`, which reports it as a warning instead of an
error. The README of the scheme describes the issue.

## `<name>_references.txt`

One line per reference genome, with this header:

```
genome	strain	genbank_id	refseq_id	assembly_name
```

`genome` is the short name used in the SNP file and may not contain an
underscore, because CanSNPer2 recovers the reference name from alignment file
names by splitting on the first one. `genbank_id` and `refseq_id`
must be well formed accessions: `genbank_id` of the form `GCA_000000000.0`
(`GCF_000000000.0` is accepted as well, as CanSNPer2 does) and `refseq_id` of
the form `GCF_000000000.0`. `CanSNPer2-download` builds NCBI paths from the
digit triplets, so a malformed accession is only noticed when a download fails.
Two rows may not give the same `genbank_id`: CanSNPer2 downloads one reference
per row and refuses a references file that repeats one.

## `<name>_genomes.txt`

Optional for an SNP scheme. The genomes a scheme was made from, so that it can
be checked again against them, and so that a typing result can be related to
known genomes. A taxonomy scheme uses a file of the same name in another
format; see Taxonomy schemes.
One line per genome, with this header:

```
genome	node	role
```

`genome` is a versioned assembly accession, `GCA_000000000.0` or
`GCF_000000000.0`, so that the genome can be fetched again; one assembly is
listed once, under either accession. `node` names a node in the tree file.
`role` is one of:

- `representative`: its sequence was in the alignment the canSNPs were found
  in, as a representative genome of a RepGenR2 run is;
- `member`: placed at the node through a representative, as dereplication folds
  it in, without being aligned itself;
- `outgroup`: aligned, but placed at no node, so that it gives the root its
  polarity. Its `node` is empty.

Two columns may follow, both or neither, to record where each genome comes
from:

```
genome	node	role	basis	source
```

`basis` is `published`, when a publication places the genome at the node, or
`typed`, when CanSNPer2 placed it there itself. `source` names the publication
by the key of the scheme README, or, for `typed`, the genome collection and
CanSNPer2 version, such as `GTDB r232, CanSNPer2 2.3.0`. A `typed` genome
confirms that the scheme calls the node; it is not independent evidence for the
node. Where a node has several typed genomes, one is listed: the type material
if there is one, otherwise the GTDB representative, otherwise the genome with
the fewest contigs and the best CheckM2 completeness and contamination. Its role
is `representative`, in the sense that it stands for the node.

A node may list no genome. CanSNPer2 does not read this file, and the database
does not depend on it. `CanSNPCreator export --genome-map ... --alignment ...`
writes it, and `CanSNPCreator audit` reads it to check every canSNP in those
genomes.

## `<name>_canSNPs.txt`

Optional. More than one canSNP per node, ranked within the node, for assays
that type several canSNPs of a branch. One line per canSNP, with this header:

```
node	rank	genome	sequence	sequence_position	derived_base	ancestral_base	score
```

`node` names a node in the tree file, and `rank` counts from 1 within the node,
in order of the file. `genome`, `sequence`, `sequence_position`, `derived_base`
and `ancestral_base` follow the SNP file: a position is 1-based on the named
reference sequence, and a position marks one node only. `score` is the total the SNP
scored when the scheme was made, which ranks it and is not comparable between
schemes scored with different settings.

CanSNPer2 calls one SNP per node and reads only the SNP file, so rank 1 of each
node must be the SNP that file gives the node, and every node of the SNP file
must have a rank 1. `CanSNPCreator export --snps-per-node N` (or `all`) writes
this file, and `CanSNPCreator audit` checks every canSNP in it.

## Taxonomy schemes

A taxonomy scheme is declared in `schemes.toml` with `kind = "taxonomy"`. Its
folder holds `<name>_tree.txt`, in the format above, `<name>_genomes.txt` and
`README.md`, and no SNP or references file. Its `<name>_genomes.txt` has two
columns and no header: a genome file name and the node the genome belongs to.
Each genome is listed once, and each node must appear in the tree.

`build.py` builds a taxonomy scheme with the flextaxd command, which continuous
integration does not install. Such a scheme is marked `wip = true` in
`schemes.toml`, so that `build.py` builds it only when it is named with
`--scheme`.

## Status

The README of every scheme gives its Status, `draft` or `verified`:

- A new scheme is added to `dev` by pull request, with Status `draft`.
- A scheme is `verified` when two things hold, and its README describes both,
  under Origin and under Validation:
  - every canSNP is published, or its origin is recorded in the README, as the
    canSNPs that are not published must be;
  - the scheme has been typed with CanSNPer2 against genomes whose lineage is
    known independently, and the calls agree with it.
- A verified scheme is moved to `master`, which holds verified schemes only.

A taxonomy scheme has no canSNPs; it is verified when its genome assignments
have been checked against an independent classification.

## `README.md`

Every scheme has a `README.md` in its folder, in this layout, so that each
scheme is described in the same way. `TEMPLATE.md` in this folder is a copy to
start from.

```markdown
# Brucella suis canSNP scheme

| | |
|---|---|
| Scheme | `brucella_suis` |
| Status | draft |
| Organisms | *Brucella suis* (NCBI taxon 29461), *Brucella canis* (taxon 36855) |
| Reference genome | *B. suis* 1330 (GCF_000007505.1) |
| Nodes | 31 |
| canSNPs per node | up to 5 |
| Origin | made with CanSNPCreator 0.5.1 |
| Last changed | 2026-11-19 |

## Description
## Origin
## Validation
## Known issues
## References
## Changes
```

The facts table has these rows, in this order:

- **Scheme**: the folder name, in backticks.
- **Status**: `draft` or `verified`; see Status above.
- **Organisms**: the species the scheme covers, with NCBI taxon IDs.
- **Reference genome**: strain and assembly accession of each reference genome
  the canSNP positions are given on; `none` for a taxonomy scheme.
- **Nodes**: the number of nodes in the tree file.
- **canSNPs per node**: `1`, or `up to N` when `<name>_canSNPs.txt` ranks
  more than one; `none` for a taxonomy scheme.
- **Origin**: in a few words, the publications the canSNPs come from, or the
  tool that made the scheme.
- **Last changed**: the date the scheme files last changed, `YYYY-MM-DD`: the
  tree, SNP, references, genomes or canSNPs file. A change to the README alone
  does not count.

The sections follow in this order, and none may be empty:

- **Description**: what the scheme resolves, and how its nodes are named.
- **Origin**: for a published scheme, which publications the canSNPs come
  from; for a scheme made with CanSNPCreator, the genomes and how they were
  chosen, the outgroup, the representatives, how the tree was built, and the
  settings.
- **Validation**: what has been checked, with which CanSNPer2 version and
  settings, and the result; also what has not been checked.
- **Known issues**: problems known in the data, or `None known.`
- **References**: full citations, with DOIs where available.
- **Changes**: one dated line per change to the scheme or its README, with the
  pull request where there is one.

A section may have `###` subsections. `validate.py` checks the layout, and that
Scheme, Nodes, Reference genome and canSNPs per node agree with the scheme
files. With `--require-verified`, as continuous integration runs it on
`master`, it also rejects a scheme whose Status is `draft`.

## Checking a change

The top-level `README.md` describes how to check and rebuild a scheme after a
change: `python validate.py`, `python build.py` and `python -m pytest tests`.
