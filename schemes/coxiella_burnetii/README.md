# Coxiella burnetii canSNP scheme

| | |
|---|---|
| Scheme | `coxiella_burnetii` |
| Status | verified |
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

Typed with CanSNPer2 2.2.2 (`--aligner minimap2`) on 2026-11-19 against all
300 genomes of the genus *Coxiella* in GTDB r232: the 156 of the species
cluster `s__Coxiella burnetii`, and 144 of twelve other species, which are
tick endosymbionts, among them *Candidatus* Coxiella mudrowiae.

| | Default (`--min_required_hits 3`) | `--min_required_hits 2` |
|---|---|---|
| *C. burnetii*, 156 genomes | 151 called | 155 called |
| other *Coxiella* species, 144 genomes | none called | none called |

The calls of the species cluster at the default threshold are `C.11` 60 times,
`C.10` 40, `C.4` 20, `C.3` 18, `C.8` 9, `C.9` 3 and `C.5` once. The four
genomes that only `--min_required_hits 2` calls are placed at `C.6` three
times and `C.1` once, the two nodes directly below the root; see Known issues.

Genomes of the same strain received the same call, and the strains of the
known genomic groups fell on separate nodes:

| Strain | Genomes | Node |
|---|---|---|
| Nine Mile clone, including RSA 493, RSA 439 phase II with its passage series, RSA 514 and RSA 363 | 17 | `C.10` |
| Dugway 5J108-111, 7D77-80 and 7E65-68 | 3 | `C.6` |
| CbuG Q212 and Scurry Q217 | 2 | `C.8` |
| CbuK Q154 | 1 | `C.4` |
| Henzerling, Z3055 | 2 | `C.11` |
| MSU Goat Q177 | 1 | `C.4` |
| Namibia | 1 | `C.5` |

Dugway, the most divergent of these strains, is placed at `C.6`, which
descends directly from the root.

Of the 125 genomes that a CanSNPCreator typing report also placed, 122
received the same call. The three others do not contradict it: two are the
Dugway genomes that the default threshold leaves uncalled, and
GCF_000820465.1 is called at `C.3` where the report gives its parent `C.1`.
Karlsson et al. (2014, Table 4) assign strains to groups, and an earlier
typing run agreed with those assignments.

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
- GCF_000826165.1 has N at five of the ten positions and is not called. It is
  the only genome of the species cluster that `--min_required_hits 2` does not
  call either.
- Three genomes carry the organism name of RSA 493 but are not placed at
  `C.10`: two are called at `C.8` and one is GCF_000826165.1 above. They are
  probably deposits that reuse the name of the reference strain; this has not
  been examined further.

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
  AE016828.2 to AE016828.3 (#20).
- 2026-11-19: typed against all 300 genomes of the genus in GTDB r232; Status
  set to verified (#26).
- 2026-11-19: moved to `master` as a verified scheme.
