# Classical Brucella canSNP scheme

| | |
|---|---|
| Scheme | `brucella_classical` |
| Status | verified |
| Organisms | *Brucella*, classical species (NCBI taxon 234) |
| Reference genome | *B. melitensis* 16M (GCF_000007125.1) |
| Nodes | 24 |
| canSNPs per node | 1 |
| Origin | published (Foster et al. 2008, 2012), ported from CanSNPer 1.0.10 |
| Last changed | 2026-11-22 |

## Description

The canSNP scheme for the classical *Brucella* species that CanSNPer 1.0.10
shipped. The classical species are *B. melitensis*, *B. abortus*, *B. suis*,
*B. canis*, *B. ovis*, *B. neotomae* and the species from marine mammals,
*B. ceti* and *B. pinnipedialis*; the scheme does not cover the atypical
species or the former genus *Ochrobactrum*, now part of *Brucella*. Seven canSNPs from Foster et al. (2008) separate the species-level
lineages (named after the gene and position, for example `abc_205` for
*B. abortus*, `cysW_288` for *B. melitensis*, `abc_246` for *B. suis* and
*B. canis*). Sixteen canSNPs from Foster et al. (2012) resolve clades within
them (named after the assay and clade letter, for example `2366A`, `8872G`,
`10621J`). The tree, node names and alleles are those of CanSNPer 1.0.10,
except for the two corrections listed under Known issues.

## Origin

Converted from `brucella_tree.txt` and `brucella_snp.txt` of CanSNPer 1.0.10
(https://github.com/adrlar/CanSNPer). CanSNPer 1 gave positions on the two
16M chromosomes, and they are given here in the same way: 22 canSNPs on
chromosome I (AE008917.1, NC_003317.1, 2,117,144 bp) and 2 on chromosome II
(AE008918.1, NC_003318.1). From 2026-11-19 to 2026-11-21 the positions were
counted along the concatenated genome, chromosome I followed by chromosome II,
as CanSNPer2 before 2.3.0 required. The CanSNPer 1 sequence of 16M is identical
to GCF_000007125.1.

The positions were checked against the publications by typing nine *Brucella*
genomes in windows of 4 bp around each published position. Each Foster et al.
(2012) position has exactly one variable site nearby, always at offset +1 or
-3, and that site is the CanSNPer 1 position, so the CanSNPer 1 positions are
used. The two branch-J assays of that paper point at the same site,
`10621J`. The Foster et al. (2008) positions match CanSNPer 1 exactly.

The genomes file `brucella_classical_genomes.txt` lists, as `representative`,
the six genomes the canSNPs were discovered from, each at its deepest node.
Foster et al. (2008) compared gene fragments in silico against six whole
genomes: *B. melitensis* 16M (`8872G`), *B. abortus* 2308 (`2366A`),
*B. abortus* 9-941 (`3740D`), *B. suis* 1330 (`8306I`), *B. ovis* 63/290
(`cysW_234`) and *B. canis* RM6/66 (`omp25_256`). Foster et al. (2012)
discovered its 1,000 microarray SNPs from four of these (2308, 9-941, 16M and
1330). All six have a current assembly and are listed. The scheme lists no
outgroup; the root canSNP is from CanSNPer 1 without a published polarity.

## Validation

`validate.py` checks the files, and every derived and ancestral pair was
checked against the base of 16M at its position. 16M is derived exactly along
its own path, `Root`, `cysW_288`, `1688H`, `6214F`, `2995G`, `8872G`.

### The genus of GTDB r232

Typed with CanSNPer2 2.2.2 (`--aligner minimap2`) on 2026-11-19 against all
1,785 genomes of the genus *Brucella* in GTDB r232: the 1,337 of the species
cluster `s__Brucella melitensis`, which holds the classical species, and 448
of 43 other species, which are the former genus *Ochrobactrum*.

| | Default (`--min_required_hits 3`) | `--min_required_hits 2` |
|---|---|---|
| the species cluster, 1,337 genomes | 304 called | 1,237 called |
| other *Brucella* species, 448 genomes | none called | none called |

Not one genome of the 448 outside the cluster received a call at either
threshold. The default threshold calls only 304 of the cluster because the
species-level nodes lie directly below the root; see Known issues.

At `--min_required_hits 2`, the species-level node each genome reached is this,
against the species in its NCBI record:

| NCBI species | Genomes | Node reached |
|---|---|---|
| *B. melitensis* | 701 | 697 `cysW_288`, 2 `abc_246`, 1 `abc_205`, 1 no call |
| *B. abortus* | 378 | 362 `abc_205`, 11 `cysW_288`, 1 `abc_246`, 4 no call |
| *B. suis* | 97 | 45 `abc_246`, 1 `cysW_288`, 51 no call |
| *B. canis* | 32 | 31 `abc_246`, 1 no call |
| *B. ceti* | 24 | 24 `trpE_290` |
| *B. ovis* | 18 | 18 `cysW_234` |
| *B. pinnipedialis* | 8 | 8 `trpE_290` |
| *B. neotomae* | 6 | 5 `rpoB_2673`, 1 no call |
| *B. microti*, *B. inopinata*, *B. vulpis*, *B. amazoniensis* | 13 | 1 `abc_246`, 12 no call |
| recorded only as *Brucella* sp. | 60 | 18 `cysW_288`, 10 `trpE_290`, 1 `abc_246`, 1 `abc_205`, 30 no call |

Fifteen genomes were called on the branch of another species than their NCBI
record names. Each was compared by average nucleotide identity (ANI) with
*B. melitensis* 16M, *B. abortus* 2308, *B. suis* 1330 and *B. canis* ATCC
23365, and in every case the ANI agrees with the call rather than with the
name: the eleven deposited as *B. abortus* and called at `cysW_288` are
closest to 16M (99.79 to 99.91%), the two deposited as *B. melitensis* and
called below `abc_246` are closest to *B. suis* or *B. canis*, the one
deposited as *B. melitensis* and called at `5754E` is 99.96% to 2308, and
strain 104M, deposited as *B. abortus*, is 99.99% to *B. suis* 1330. No call
contradicts the identity of its genome.

Of the 100 genomes of the cluster that `--min_required_hits 2` does not call,
94 carry fewer than two derived canSNPs: 51 are *B. suis*, whose biovars 2 and
5 no canSNP of the scheme marks, and 12 are the atypical species *B. microti*,
*B. inopinata*, *B. vulpis* and *B. amazoniensis*, which the scheme does not
cover.

### Reference and type strains

Earlier, 24 *Brucella* reference and type strains were typed, and every call
agreed with the species of the strain:

| Strain | Species and biovar | Call, default | Call, `--min_required_hits 2` |
|---|---|---|---|
| 544 | *B. abortus* bv1 | `2366A` | `2366A` |
| 86/8/59 | *B. abortus* bv2 | `5754E` | `5754E` |
| Tulya | *B. abortus* bv3 | NA | `abc_205` |
| 292 | *B. abortus* bv4 | `5754E` | `5754E` |
| B3196 | *B. abortus* bv5 | NA | `abc_205` |
| 870 | *B. abortus* bv6 | NA | `abc_205` |
| 63/75 | *B. abortus* bv7 | NA | `abc_205` |
| C68 | *B. abortus* bv9 | NA | `abc_205` |
| 16M | *B. melitensis* bv1 | `8872G` | `8872G` |
| 63/9 (ATCC 23457) | *B. melitensis* bv2 | NA | `cysW_288` |
| Ether | *B. melitensis* bv3 | NA | `cysW_288` |
| 1330 | *B. suis* bv1 | `8306I` | `8306I` |
| Thomsen (ATCC 23445) | *B. suis* bv2 | NA | NA |
| 686 | *B. suis* bv3 | NA | `abc_246` |
| 40 | *B. suis* bv4 | `10621J` | `10621J` |
| 513 | *B. suis* bv5 | NA | NA |
| RM6/66 (ATCC 23365) | *B. canis* | `omp25_256` | `omp25_256` |
| 63/290 (ATCC 25840) | *B. ovis* | NA | `cysW_234` |
| 5K33 (NCTC 10084) | *B. neotomae* | NA | `rpoB_2673` |
| B1/94 | *B. ceti* | NA | `trpE_290` |
| B2/94 | *B. pinnipedialis* | NA | `trpE_290` |
| CCM 4915 | *B. microti* | NA | NA |
| F60 | *B. vulpis* | NA | NA |
| BO1 | *B. inopinata* | NA | NA |

The scheme has not been checked against published *Brucella* canSNP typing of
these genomes below species level.

## Known issues

- Nodes at depth one or two, `abc_205`, `cysW_288`, `abc_246`, `cysW_234`,
  `trpE_290` and `rpoB_2673`, carry two derived canSNPs on their path including
  `Root`, so CanSNPer2 does not call them with its default
  `--min_required_hits 3`. Use `--min_required_hits 2` to type to these nodes.
  Over the species cluster of GTDB r232 this is the difference between 304 and
  1,237 of 1,337 genomes called.
- No canSNP of the scheme marks the lineages of *B. suis* biovars 2 and 5, so
  genomes of those biovars stop above `abc_246`, nor the atypical species
  *B. microti*, *B. inopinata*, *B. vulpis* and *B. amazoniensis*, which stop
  at `Root`.
- Fifteen genomes of the species cluster are deposited under a species name
  that neither their call nor their ANI supports; the scheme places them with
  their ANI. See Validation.
- `Root` is not a published canSNP: CanSNPer 1 gives no reference for it.
- `2922C` and `991C` give the ancestral base on the other strand in CanSNPer
  1.0.10 (C and G). They are corrected here to G and C; the corrected alleles
  agree with Foster et al. (2012) and with seven genomes.
- For `rpoB_2673` and `cysW_288`, the allele letters of Foster et al. (2008)
  disagree with the genomes. Typing *B. neotomae* 5K33 and 16M supports the
  polarity of CanSNPer 1, which is kept.
- `10621J` is derived in *B. abortus* 9-941 as well as in *B. canis*: a
  homoplasy or an incongruent marker. It does not change the call of 9-941
  (`3740D`).
- In an earlier comparison over 1,337 genomes, nine *B. abortus* genomes were
  derived at `2922C` and `991C` but not at `1562B` and `4748B`, which the tree
  places above them. The nesting of the C clade under the B clade has not been
  checked further.

## References

- Foster2008: Foster JT, Okinaka RT, Svensson R, et al. Real-time PCR assays
  of single-nucleotide polymorphisms defining the major *Brucella* clades. J
  Clin Microbiol. 2008;46(1):296-301. https://doi.org/10.1128/JCM.01496-07
- Foster2012: Foster JT, Price LB, Beckstrom-Sternberg SM, et al. Genotyping
  of *Brucella* species using clade specific SNPs. BMC Microbiol.
  2012;12:110. https://doi.org/10.1186/1471-2180-12-110
- CanSNPer 1: Lärkeryd A, Myrtennäs K, Karlsson E, et al. CanSNPer: a
  hierarchical genotype classifier of clonal pathogens. Bioinformatics.
  2014;30(12):1762-1764. https://doi.org/10.1093/bioinformatics/btu113

## Changes

- 2026-11-19: ported from CanSNPer 1.0.10, with the alleles of `2922C` and
  `991C` corrected and positions given on the concatenated 16M genome (#20).
- 2026-11-19: renamed from `brucella` to `brucella_classical`, so that
  `brucella` can name a scheme for the whole genus (#23).
- 2026-11-19: typed against all 1,785 genomes of the genus in GTDB r232;
  Status set to verified (#28).
- 2026-11-19: moved to `master` as a verified scheme.
- 2026-11-21: SNP positions given per reference sequence (AE008917.1 and
  AE008918.1) with `CanSNPer2-database --convert`; database rebuilt with
  CanSNPer2 2.3.0.
- 2026-11-21: added `brucella_classical_genomes.txt`, the six genomes the
  canSNPs were discovered from (Foster et al. 2008, 2012).
- 2026-11-22: `brucella_classical_genomes.txt` gained the columns `basis` and
  `source` (all rows `published`). The leaves `rpoB_2673` and `trpE_290` have no
  genome; no GTDB genome is called there.
- 2026-11-22: `brucella_classical_canSNPs.txt` replaces the SNP file; every genome of the
  genomes file already sat at a terminal node.
