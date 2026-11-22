# Bacillus anthracis canSNP scheme

| | |
|---|---|
| Scheme | `bacillus_anthracis` |
| Status | verified |
| Organisms | *Bacillus anthracis* (NCBI taxon 1392) |
| Reference genome | Ames (GCF_000007845.1) |
| Nodes | 13 |
| canSNPs per node | 1 |
| Origin | published (Van Ert et al. 2007) |
| Last changed | 2026-11-22 |

## Description

The canSNP scheme of Van Ert et al. (2007), which divides *B. anthracis* into
the A and B branches and their sub-lineages. Nodes carry the published names,
from the root `A/B.Br.001` to the `A.Br` and `B.Br` branch points.

### Node names

The next free node names in the published *B. anthracis* nomenclature are
`A.Br.177`, `B.Br.012` and `C.Br.003`. This scheme holds 13 of the named
nodes.

## Origin

All 13 canSNPs and their positions on the Ames genome are taken from Van Ert et
al. (2007).

Van Ert et al. (2007) sequenced seven genomes to find the canSNPs: the
assemblies that the Institute for Genomic Research (TIGR) deposited and the Ames
Ancestor sequence (NC_007530). Each defines a lineage of that paper: Ames
Ancestor at `A.Br.001`, Australia 94 at `A.Br.003`, Vollum at `A.Br.007`,
Western North America USA6153 at `A.Br.009`, Kruger B at `B.Br.001` and
CNEVA-9066 at `B.Br.004`, with A1055, the genome of the C lineage, as the
outgroup, since the scheme covers the A and B branches only. The paper names its
seventh genome A0155 in the supplementary text and A1055 in Figure 1 and Table
2; A1055 is taken here.

The genomes file `bacillus_anthracis_genomes.txt` lists the five of those that
sit at a terminal node, as `representative`, and A1055 as the `outgroup`.
Australia 94 (`GCF_000167335.1`) is not listed: `A.Br.003` carries further
nodes below it, and the genomes file holds the genomes of the terminal nodes.
Every one of the five terminal nodes has a genome.

## Validation

`validate.py` checks the files, and `build.py` rebuilds the database and
compares it with the committed one (CanSNPer2 2.2.0). When the build tooling
was added (#13), the database rebuilt from these files was identical in nodes,
edges, SNP annotations and references to the one committed before.

Typed with CanSNPer2 2.2.2 (`--aligner minimap2`, default settings) on
2026-11-19:

- **One representative per species of the genus.** Of the 70 *Bacillus_A*
  species representatives of GTDB r232 that have a sequence, only
  *B. anthracis* Vollum was called, at `A.Br.007`; the path
  `A/B.Br.001`, `A.Br.006`, `A.Br.007` is the one that defines the
  A.Br.Vollum group in Van Ert et al. (2007). The other 69 species, among
  them *B. cereus* and *B. paranthracis*, received no call. The
  representative of *B. wiedmannii* could not be typed, as its assembly is
  suppressed at NCBI.
- **The species cluster.** All 947 genomes of the GTDB r232 species
  `s__Bacillus_A anthracis` were typed. Of the 497 that NCBI records as
  *B. anthracis*, 413 (83%) were called. None of the 410 genomes that NCBI
  records as another species, 363 of them *B. cereus* and 44
  *B. thuringiensis*, received a call.

Of the 84 *B. anthracis* genomes without a call, 83 carry fewer than three
derived canSNPs, which the default `--min_required_hits 3` requires; see Known
issues for the rest. The calls have not been compared with an independent
typing of these genomes.

## Known issues

- Of the 947 genomes in the GTDB species cluster, 59 that NCBI records as
  *B. anthracis* are ancestral at every canSNP of the scheme, including the
  root `A/B.Br.001`, so the scheme places them outside the A/B branch. Their
  identity has not been examined further.
- One genome (GCA_001278435.1) is derived at `A.Br.004`, `A.Br.003`,
  `A.Br.002` and `A.Br.001` but ancestral at `A.Br.006`, which the tree places
  above them. CanSNPer2 needs an unbroken path from the root, so the genome
  receives no call. Either the position of `A.Br.006` in the tree or its
  alleles should be checked.

## References

- Van Ert MN, Easterday WR, Huynh LY, et al. Global genetic population
  structure of *Bacillus anthracis*. PLoS ONE. 2007;2(5):e461.
  https://doi.org/10.1371/journal.pone.0000461

## Changes

- 2023-03-13: added from the earlier CanSNPer data; tree columns separated by
  tabs instead of semicolons.
- 2026-11-11: database rebuilt with CanSNPer2 2.2.0 (#13).
- 2026-11-19: files moved to `schemes/bacillus_anthracis/` (#16); README added.
- 2026-11-19: next free node names moved here from the repository README.
- 2026-11-19: typed against the GTDB r232 species representatives of the genus
  and the whole species cluster; the result is under Validation.
- 2026-11-21: SNP positions given per reference sequence (AE016879.1)
  with `CanSNPer2-database --convert`; database rebuilt with CanSNPer2
  2.3.0.
- 2026-11-21: added `bacillus_anthracis_genomes.txt`, the seven genomes of Van Ert
  et al. (2007).
- 2026-11-22: `bacillus_anthracis_genomes.txt` gained the columns `basis` and `source` (all
  rows `published`).
- 2026-11-22: `bacillus_anthracis_canSNPs.txt` replaces the SNP file, and
  `bacillus_anthracis_genomes.txt` keeps only the genomes of the terminal nodes
  (Australia 94, at `A.Br.003`, was dropped).
