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

A canSNP scheme for *B. suis* and *B. canis*. *B. canis* is included because
it lies within *B. suis*. Nodes are named after the biovar or species that the
genomes below them carry in their NCBI BioSample records (`Bsu.bv1` to
`Bsu.bv5`, `Bsu.canis`). Nodes below a named node keep the rest of their code,
and nodes without a clear label keep their code (`Bsu.1.1`). The names rest on
the 54 genomes with a recorded biovar and should be read as indicative.

## Origin

Made with CanSNPCreator 0.5.1 from the current RefSeq assemblies of *B. suis*
(120 assemblies) and *B. canis* (75 assemblies), with *B. melitensis* 16M
(GCF_000007125.1) as the outgroup and *B. suis* 1330 as the reference. Forty
representatives were chosen from a tree over all 195 genomes (20,000 SNP sites,
FastTree), cut into 40 clusters. A tree over the representatives was built
with RepGenR2 and FastTree, and nodes with less than 90% support were
collapsed. The scheme keeps the 31 nodes that have a canSNP
(`--snps-per-node 5 --prune`), with up to five ranked canSNPs per node in
`brucella_suis_canSNPs.txt`. The `reference` column of the SNP file is `NA`, as
the canSNPs are not taken from a publication.

## Validation

With CanSNPer2 2.2.2 (`--aligner minimap2`, default settings), the 155 genomes
not used as representatives were each called at the node the CanSNPCreator
typing report gives them. Against a tree over all genomes, that node agreed
with the genome's place in the tree for 65.8% of them, was deeper for 21.9% and
higher for 11.6%. Of the 40 representatives, 36 were called at their node and
three, at a node without canSNPs that was removed, at its parent. The scheme
has not been checked against published *B. suis* typing.

## Known issues

- One divergent representative at the root (GCF_000292125.1) is not called.
- Two genomes carry a biovar label that does not fit their node (one bv2 in
  `Bsu.bv1.1`, one bv1 in `Bsu.bv2.1.1.1.1`).

## References

None; the canSNPs are not taken from a publication. CanSNPCreator:
https://github.com/FOI-Bioinformatics/CanSNPCreator

## Changes

- 2026-11-18: added (#14).
- 2026-11-19: files moved to `schemes/brucella_suis/`; README in the standard
  layout.
