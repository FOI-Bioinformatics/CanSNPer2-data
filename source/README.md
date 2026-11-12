# Scheme sources

These tab-separated text files are the authored form of every canSNP scheme in
this repository. The SQLite databases under `database/` are generated from them
and are not edited directly.

All files use LF line endings, carry no trailing whitespace and end with a
newline. Carriage returns matter: a stray `\r` can end up inside a node name. `.gitattributes` and `validate.py` both enforce this.

A semicolon may not appear anywhere in these files. The CanSNPer tree reader
converts tabs to semicolons before splitting, so a semicolon inside a name is
read as a column break.

## Schemes

A SNP scheme is named after its organism and consists of three files:
`<organism>_tree.txt`, `<organism>_snp.txt` and `<organism>_references.txt`,
and optionally `<organism>_genomes.txt`, listing the genomes it was made from,
and `<organism>_canSNPs.txt`, listing more than one canSNP per node.
Adding those three files is all that is needed to declare a new scheme; both
scripts discover it by name.

## `<organism>_tree.txt`

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

## `<organism>_snp.txt`

One line per canSNP, with this header:

```
snp_id	strain	reference	genome	position	derived_base	ancestral_base
```

`snp_id` must name a node in the tree file, and every node in the tree file must
have exactly one SNP here. `genome` must name a row in the references file.
`position` is a positive integer, one-based, on that reference genome.
`derived_base` and `ancestral_base` are single characters from A, C, G and T and
must differ from each other. No two SNPs may annotate the same position on the
same reference genome.

Three Yersinia canSNPs (`1.ANT2.b`, `1.ORI2.c`, `1.ORI3.b`) are published with
`/` in the `derived_base` column and are listed under `known_issues` in
`schemes.toml`. CanSNPer2 compares an aligned base against both columns, so those
three can never be called derived. In tests with simulated sequence, samples on
these branches were placed at a node above them, and a sample at `1.ORI2.d`, which
descends from `1.ORI2.c`, was placed at `1.ORI2.b`. No warning was given. The
correct bases have not been established; CO92 carries the ancestral base at all
three positions, so it cannot supply them.

## `<organism>_references.txt`

One line per reference genome, with this header:

```
genome	strain	genbank_id	refseq_id	assembly_name
```

`genome` is the short name used in the SNP file and may not contain an
underscore, because CanSNPer2 recovers the reference name from alignment file
names by splitting on the first one. `genbank_id` and `refseq_id`
must be well formed accessions, `GCA_000000000.0` and `GCF_000000000.0`:
`CanSNPer2-download` builds NCBI paths from the digit triplets, so a malformed
accession is only noticed when a download fails.

## `<organism>_genomes.txt`

Optional. The genomes a scheme was made from, so that it can be checked again
against them, and so that a typing result can be related to known genomes.
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

A node may list no genome. CanSNPer2 does not read this file, and the database
does not depend on it. `CanSNPCreator export --genome-map ... --alignment ...`
writes it, and `CanSNPCreator audit` reads it to check every canSNP in those
genomes.

## `<organism>_canSNPs.txt`

Optional. More than one canSNP per node, ranked within the node, for assays
that type several canSNPs of a branch. One line per canSNP, with this header:

```
node	rank	genome	position	derived_base	ancestral_base	score
```

`node` names a node in the tree file, and `rank` counts from 1 within the node,
in order of the file. `genome`, `position`, `derived_base` and `ancestral_base`
follow the SNP file: positions are 1-based over the concatenated reference
genome, and a position marks one node only. `score` is the total the SNP
scored when the scheme was made, which ranks it and is not comparable between
schemes scored with different settings.

CanSNPer2 calls one SNP per node and reads only the SNP file, so rank 1 of each
node must be the SNP that file gives the node, and every node of the SNP file
must have a rank 1. `CanSNPCreator export --snps-per-node N` (or `all`) writes
this file, and `CanSNPCreator audit` checks every canSNP in it.

## Checking a change

```sh
python validate.py
python build.py
```

`validate.py` checks the files described above. `build.py` rebuilds every
database and compares it against the committed one, so it also reports a scheme
whose database has not been regenerated.
