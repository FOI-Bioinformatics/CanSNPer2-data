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
| Last changed | 2023-03-13 |

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

## Validation

`validate.py` checks the files, and `build.py` rebuilds the database and
compares it with the committed one (CanSNPer2 2.2.0). When the build tooling
was added (#13), the database rebuilt from these files was identical in nodes,
edges, SNP annotations and references to the one committed before. The scheme
has not been tested here against a set of typed *B. anthracis* genomes.

## Known issues

None known.

## References

- Van Ert MN, Easterday WR, Huynh LY, et al. Global genetic population
  structure of *Bacillus anthracis*. PLoS ONE. 2007;2(5):e461.
  doi:10.1371/journal.pone.0000461

## Changes

- 2023-03-13: added from the earlier CanSNPer data; tree columns separated by
  tabs instead of semicolons.
- 2026-11-11: database rebuilt with CanSNPer2 2.2.0 (#13).
- 2026-11-19: files moved to `schemes/bacillus_anthracis/` (#16); README added.
- 2026-11-19: next free node names moved here from the repository README.
