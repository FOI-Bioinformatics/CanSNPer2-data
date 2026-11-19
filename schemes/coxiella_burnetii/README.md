# Coxiella burnetii canSNP scheme

| | |
|---|---|
| Scheme | `coxiella_burnetii` |
| Status | draft |
| Organisms | *Coxiella burnetii* (NCBI taxon 777) |
| Reference genome | RSA 493 (GCF_000007765.2) |
| Nodes | 11 |
| canSNPs per node | 1 |
| Origin | published (Karlsson et al. 2014), ported from CanSNPer 1.0.10 |
| Last changed | 2026-11-19 |

## Description

The canSNP scheme for *C. burnetii* that CanSNPer 1.0.10 shipped, with the
ten canSNPs of Karlsson et al. (2014) for the major groups `C.1` to `C.11`.
The tree, node names and alleles are those of CanSNPer 1.0.10.

## Origin

Converted from `coxiella_burnetii_tree.txt` and `coxiella_burnetii_snp.txt` of
CanSNPer 1.0.10 (https://github.com/adrlar/CanSNPer). CanSNPer 1 gave the
positions on AE016828.2, an earlier version of the RSA 493 chromosome. They
were moved to the current chromosome, AE016828.3 (NC_002971.4) in
GCF_000007765.2, by a unique match of 40 bp on each side; the shifts are
between 0 and +206 bp and no base changed. Karlsson et al. (2014, Table 3)
give the same positions and alleles on AE016828.2 as CanSNPer 1.

## Validation

`validate.py` checks the files, and every derived and ancestral pair was
checked against the base of RSA 493 at its position. RSA 493 is derived
exactly along its own path, `Root`, `C.7`, `C.9`, `C.10`.

With CanSNPer2 2.2.2 (`--aligner minimap2`), 156 *C. burnetii* RefSeq
assemblies were typed. With the default `--min_required_hits 3`, 151 were
called (C.11 60, C.10 40, C.4 20, C.3 18, C.8 9, C.9 3, C.5 1) and 5 were not.
With `--min_required_hits 2`, four of those five were called (`C.6` three
times, `C.1` once). For 124 of the 125 genomes that a CanSNPCreator typing
report also placed, the call is the same; the exception, GCF_000820465.1, is
placed at `C.1` by the report and at `C.3` by CanSNPer2, which finds the
derived base of `C.3`. Karlsson et al. (2014, Table 4) assign strains to
groups, and an earlier typing run agreed with those assignments.

## Known issues

- `Root` is not a published canSNP. CanSNPer 1 notes that it was found in an
  alignment with *Rickettsiella grylli*, used to root the tree.
- Karlsson et al. (2014) could not determine the polarity of `C.1` and `C.2`
  for lack of an outgroup, and made one assay for the two. That assay is
  `C.1` here, so `C.1` may sit on either branch, and the scheme has no `C.2`.
- `C.6` and `C.1` carry two derived canSNPs on their path including `Root`,
  so CanSNPer2 does not call them with its default `--min_required_hits 3`.
  Use `--min_required_hits 2` to type to these nodes.
- Karlsson et al. (2014) name the branches `C.12` to `C.15` but give no assay
  for them, so genomes on those branches stop at a node above.
- GCF_000826165.1 has N at five of the ten positions and is not called.

## References

- Karlsson2014: Karlsson E, Macellaro A, Byström M, et al. Eight new genomes
  and synthetic controls increase the accessibility of rapid melt-MAMA SNP
  typing of *Coxiella burnetii*. PLoS ONE. 2014;9(1):e85417.
  doi:10.1371/journal.pone.0085417
- CanSNPer 1: Lärkeryd A, Myrtennäs K, Karlsson E, et al. CanSNPer: a
  hierarchical genotype classifier of clonal pathogens. Bioinformatics.
  2014;30(12):1762-1764. doi:10.1093/bioinformatics/btu113

## Changes

- 2026-11-19: ported from CanSNPer 1.0.10, with positions moved from
  AE016828.2 to AE016828.3.
