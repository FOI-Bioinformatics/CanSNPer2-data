# Francisellaceae taxonomy scheme

| | |
|---|---|
| Scheme | `francisellaceae` |
| Status | draft |
| Organisms | *Francisellaceae* (NCBI taxon 34064) |
| Reference genome | none |
| Nodes | 67 |
| canSNPs per node | none |
| Origin | compiled at FOI |
| Last changed | 2026-11-19 |

## Description

A taxonomy of the family *Francisellaceae*, from the family to the genera
*Francisella*, *Allofrancisella*, *Parafrancisella*, *Pseudofrancisella*,
*Caedibacter*, *Fangia* and others, and within *F. tularensis* to the
subspecies and the groups A.I, A.II and B2. It is a taxonomy scheme: it holds
a tree and the genomes assigned to its nodes, and no canSNPs. Internal nodes
without a taxon name carry codes such as `F2` and `C1001`.

`francisellaceae_genomes.txt` has two columns and no header: a genome file
name and the node the genome belongs to. It lists 505 genomes, 245 of them
at `B2`.

## Origin

The tree and the genome assignments were added in 2021 (commit 7412278) and
were compiled at FOI. How the tree was inferred is not recorded in this
repository.

## Validation

`validate.py` checks that the tree has a single root, that every node has one
parent, and that every genome is assigned to a node of the tree. The genome
assignments have not been checked otherwise.

## Known issues

- CanSNPer2 cannot type against this database: it holds a taxonomy and genome
  assignments but no SNP annotations, and CanSNPer2 rejects it for lacking the
  `snp_annotation` and `snp_references` tables.
- CanSNPer2-database has no option for a genome-to-node file, so the scheme is
  built with the flextaxd command, which continuous integration does not
  install; `build.py` builds it only on request
  (`python build.py --scheme francisellaceae`). flextaxd master stops with
  `NameError: name 'cnp' is not defined` when reading a CanSNPer tree, and the
  committed database was built with a local, unpublished correction to
  `flextaxd/modules/ReadTaxonomyCanSNPer.py`. The flextaxd version in the
  manifest does not reflect that correction.
- No database has been published for this scheme, so a rebuild has no
  independent baseline; `build.py` reports a difference rather than failing on
  one.

## References

None; the taxonomy is not taken from a publication.

## Changes

- 2021-11-02: tree and genome assignments added.
- 2026-11-10: database built with flextaxd and CanSNPer2 2.1.0, under
  `work_in_progress/`.
- 2026-11-19: moved to `schemes/francisellaceae/`, with its database in
  `database/`; README in the standard layout.
