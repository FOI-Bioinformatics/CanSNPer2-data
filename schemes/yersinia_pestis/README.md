# Yersinia pestis canSNP scheme

| | |
|---|---|
| Scheme | `yersinia_pestis` |
| Status | verified |
| Organisms | *Yersinia pestis* (NCBI taxon 632) |
| Reference genome | CO92 (GCF_000009065.1) |
| Nodes | 100 |
| canSNPs per node | 1 |
| Origin | published (Morelli et al. 2010) |
| Last changed | 2026-11-04 |

## Description

The canSNP scheme of Morelli et al. (2010) for the global population of *Y.
pestis*. Nodes carry the published branch names (for example `0.PE4.a`,
`1.ORI2.c`, `2.MED1.a`); the Roman numerals (`III` to `XVII`) are the internal
nodes of that phylogeny, and the root is `root`.

## Origin

All 100 canSNPs and their positions on the CO92 genome are taken from Morelli
et al. (2010).

## Validation

`validate.py` checks the files, and `build.py` rebuilds the database and
compares it with the committed one (CanSNPer2 2.2.0). When the build tooling
was added (#13), the database rebuilt from these files was identical in nodes,
edges, SNP annotations and references to the one committed before. Samples
simulated on the branches named under Known issues were typed with CanSNPer2;
the result is described there. The scheme has not been tested here against a
set of typed *Y. pestis* genomes.

## Known issues

Three canSNPs (`1.ANT2.b`, `1.ORI2.c`, `1.ORI3.b`) are published with `/` in
the `derived_base` column and are listed under `known_issues` in
`schemes.toml`. CanSNPer2 compares an aligned base against both columns, so
those three can never be called derived. In tests with simulated sequence,
samples on these branches were placed at a node above them, and a sample at
`1.ORI2.d`, which descends from `1.ORI2.c`, was placed at `1.ORI2.b`. No
warning was given. The correct bases have not been established; CO92 carries
the ancestral base at all three positions, so it cannot supply them. An attempt
to read the derived bases from RefSeq assemblies did not recover them, so the
three rows are left as published until the bases are taken from the original
description of the SNPs.

## References

- Morelli G, Song Y, Mazzoni CJ, et al. *Yersinia pestis* genome sequencing
  identifies patterns of global phylogenetic diversity. Nat Genet.
  2010;42(12):1140-1143. doi:10.1038/ng.705

## Changes

- 2023-03-13: added from the earlier CanSNPer data.
- 2026-11-04: files normalised to LF line endings without trailing whitespace.
- 2026-11-11: database rebuilt with CanSNPer2 2.2.0 (#13).
- 2026-11-19: files moved to `schemes/yersinia_pestis/` (#16); README added.
