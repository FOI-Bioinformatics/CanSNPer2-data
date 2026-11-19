# Brucella canSNP scheme

| | |
|---|---|
| Scheme | `brucella` |
| Status | draft |
| Organisms | *Brucella*, classical species (NCBI taxon 234) |
| Reference genome | *B. melitensis* 16M (GCF_000007125.1) |
| Nodes | 24 |
| canSNPs per node | 1 |
| Origin | published (Foster et al. 2008, 2012), ported from CanSNPer 1.0.10 |
| Last changed | 2026-11-19 |

## Description

The canSNP scheme for the classical *Brucella* species that CanSNPer 1.0.10
shipped. Seven canSNPs from Foster et al. (2008) separate the species-level
lineages (named after the gene and position, for example `abc_205` for
*B. abortus*, `cysW_288` for *B. melitensis*, `abc_246` for *B. suis* and
*B. canis*). Sixteen canSNPs from Foster et al. (2012) resolve clades within
them (named after the assay and clade letter, for example `2366A`, `8872G`,
`10621J`). The tree, node names and alleles are those of CanSNPer 1.0.10,
except for the two corrections listed under Known issues.

## Origin

Converted from `brucella_tree.txt` and `brucella_snp.txt` of CanSNPer 1.0.10
(https://github.com/adrlar/CanSNPer). CanSNPer 1 gave positions on the two
16M chromosomes; they are given here on the concatenated genome, chromosome I
(NC_003317.1, 2,117,144 bp) followed by chromosome II (NC_003318.1). The
CanSNPer 1 sequence of 16M is identical to GCF_000007125.1.

The positions were checked against the publications by typing nine *Brucella*
genomes in windows of 4 bp around each published position. Each Foster et al.
(2012) position has exactly one variable site nearby, always at offset +1 or
-3, and that site is the CanSNPer 1 position, so the CanSNPer 1 positions are
used. The two branch-J assays of that paper point at the same site,
`10621J`. The Foster et al. (2008) positions match CanSNPer 1 exactly.

## Validation

`validate.py` checks the files, and every derived and ancestral pair was
checked against the base of 16M at its position. 16M is derived exactly along
its own path, `Root`, `cysW_288`, `1688H`, `6214F`, `2995G`, `8872G`.

With CanSNPer2 2.2.2 (`--aligner minimap2`), 24 *Brucella* reference and type
strains were typed:

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

Every call agrees with the species of the strain. With the default
`--min_required_hits 3`, 17 strains are not called; see Known issues. *B. suis*
bv2 and bv5 and the atypical species are derived at `Root` only, as no canSNP
in the scheme marks their lineages. The 76 *B. anthropi* and *B. cytisi*
genomes of the `brucella_anthropi` scheme received no call at either setting.

## Known issues

- Nodes at depth one or two, `abc_205`, `cysW_288`, `abc_246`, `cysW_234`,
  `trpE_290` and `rpoB_2673`, carry two derived canSNPs on their path including
  `Root`, so CanSNPer2 does not call them with its default
  `--min_required_hits 3`. Use `--min_required_hits 2` to type to these nodes.
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
  Clin Microbiol. 2008;46(1):296-301. doi:10.1128/JCM.01496-07
- Foster2012: Foster JT, Price LB, Beckstrom-Sternberg SM, et al. Genotyping
  of *Brucella* species using clade specific SNPs. BMC Microbiol.
  2012;12:110. doi:10.1186/1471-2180-12-110
- CanSNPer 1: Lärkeryd A, Myrtennäs K, Karlsson E, et al. CanSNPer: a
  hierarchical genotype classifier of clonal pathogens. Bioinformatics.
  2014;30(12):1762-1764. doi:10.1093/bioinformatics/btu113

## Changes

- 2026-11-19: ported from CanSNPer 1.0.10, with the alleles of `2922C` and
  `991C` corrected and positions given on the concatenated 16M genome.
