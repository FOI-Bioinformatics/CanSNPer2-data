# Brucella anthropi canSNP scheme

| | |
|---|---|
| Scheme | `brucella_anthropi` |
| Status | draft |
| Organisms | *Brucella anthropi* (NCBI taxon 529), *Brucella cytisi* (taxon 407152) |
| Reference genome | *B. anthropi* ATCC 49188 (GCF_000017405.1) |
| Nodes | 42 |
| canSNPs per node | up to 5 |
| Origin | made with CanSNPCreator 0.5.1 |
| Last changed | 2026-11-21 |

## Description

A canSNP scheme for *B. anthropi* (formerly *Ochrobactrum anthropi*) and
*B. cytisi*. In the tree the scheme was made from, the two species are sister
groups, and *B. anthropi* divides into two lineages, each with 100% bootstrap
support. The nodes are named after them: `Bant.cytisi`; `Bant.anthropi`;
`Bant.anthropi.T`, the lineage of the type strain (26 genomes, 12 of the 17
with a recorded source clinical); `Bant.anthropi.T.E`, a group of five genomes
within it, all from Europe; and `Bant.anthropi.env` (43 genomes, 22 of the 31
with a recorded source environmental). Nodes below a named node keep the rest
of their code. The names describe the genomes at hand and are not established
designations.

## Origin

Made with CanSNPCreator 0.5.1 from the RefSeq assemblies of *B. anthropi* and
*B. cytisi*. Three assemblies filed as *B. anthropi* were left out because
their average nucleotide identity (ANI) to the type strain placed them in
another species, which leaves 76 genomes. One genome filed as *B. cytisi*
(WTB2) falls within *B. anthropi*, in agreement with its ANI. The reference is
the *B. anthropi* type strain ATCC 49188. The outgroup is *B. tritici*
(GCF_008932295.1), the closest species outside the two by ANI (91.5-91.8% to
both type strains). Fifty-seven representatives were chosen by RepGenR2
dereplication, and the tree over them was built with IQ-TREE (GTR+F+I+R4, 1,000
ultrafast bootstrap replicates). Nodes with less than 90% bootstrap support
were collapsed. The scheme keeps the 42 nodes that have a canSNP
(`--snps-per-node 5 --prune`), with up to five ranked canSNPs per node in
`brucella_anthropi_canSNPs.txt`. The `reference` column of the SNP file is
`NA`, as the canSNPs are not taken from a publication.

The positions in `brucella_anthropi_snp.txt` and
`brucella_anthropi_canSNPs.txt` are given on each reference sequence of
ATCC 49188: 35 of the 42 scheme SNPs lie on chromosome 1 (CP000758.1,
2,887,297 bp) and 7 on chromosome 2 (CP000759.1, 1,895,911 bp); none lie on the
four plasmids. Until 2026-11-21 the positions were counted along the
concatenated genome, as CanSNPer2 before 2.3.0 required.

## Validation

With CanSNPer2 2.2.2 (`--aligner minimap2`, default settings), 75 of the 76
genomes were called at the node the CanSNPCreator typing report gives them, or
at its nearest ancestor in the scheme where that node has no canSNP; this
includes all 19 genomes not used as representatives. The scheme has not been
checked against published typing of these species. The results above
were obtained with CanSNPer2 2.2.2, before the positions were converted to
positions per reference sequence. Typing with CanSNPer2 2.3.0 after the
conversion gave the same calls with all three aligners.

## Known issues

- One *B. cytisi* representative (GCF_058245555.1) sits directly at
  `Bant.cytisi` and carries only two derived canSNPs on its path, `Bant` and
  `Bant.cytisi`, so it is not called with the default
  `--min_required_hits 3`.

## References

None; the canSNPs are not taken from a publication. CanSNPCreator:
https://github.com/FOI-Bioinformatics/CanSNPCreator

## Changes

- 2026-11-18: added (#15).
- 2026-11-19: files moved to `schemes/brucella_anthropi/`; README in the
  standard layout.
- 2026-11-21: SNP positions given per reference sequence (CP000758.1 and
  CP000759.1) with `CanSNPer2-database --convert`; database rebuilt with
  CanSNPer2 2.3.0.
