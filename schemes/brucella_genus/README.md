# Brucella genus canSNP scheme

| | |
|---|---|
| Scheme | `brucella_genus` |
| Status | draft |
| Organisms | *Brucella*, the genus as GTDB r232 delimits it: the classical species and the former genus *Ochrobactrum* |
| Reference genome | *B. melitensis* 16M (GCF_000007125.1) |
| Nodes | 283 |
| canSNPs per node | up to 3 |
| Origin | made with CanSNPCreator 0.9.0 and RepGenR2 2.3.0 |
| Last changed | 2026-11-24 |

## Description

A canSNP scheme for the whole genus *Brucella* as GTDB r232 delimits it, from
the species of the former genus *Ochrobactrum* down to the species and biovars
of the classical brucellae. The scheme was built in two tiers and joined at the
node `BR.classical`, which holds the 1,337 genomes of the GTDB species cluster
`s__Brucella melitensis`, the classical species.

Above `BR.classical` (112 nodes) a node is named after a GTDB species, as
`BR.anthropi` or `BR.intermedia`, where every genome below it belongs to that
species. Below it (171 nodes) the nodes are named after the species and biovars
of 23 reference strains, as `BR.classical.abortus`, `BR.classical.abortus_bv1`
or `BR.classical.ceti`. Nodes below a named node keep the rest of their code,
as `BR.classical.melitensis.1.2`. The names describe the genomes at hand and
are not established designations.

*B. suis* is not one clade in this scheme, as *B. canis* falls within it; its
four parts are named after the biovar of the reference strain each holds
(`suis_bv1` to `suis_bv4`). See Known issues for the labels that name no node.

## Origin

Made from all 1,785 genomes of the genus *Brucella* in GTDB r232, in two tiers
built with RepGenR2 2.3.0 and CanSNPCreator 0.9.0. Positions were taken
against *B. melitensis* 16M (GCF_000007125.1) and are given on its GenBank
sequences, chromosome I (AE008917.1) and chromosome II (AE008918.1).

**Genus tier.** The 1,785 genomes were dereplicated with skder 1.3.6 at 99%
average nucleotide identity, keeping as representatives the 66 genomes of
species and reference strains listed for the build, which left 212
representatives. The outgroup is *Paenochrobactrum pullorum*
(GCF_041929845.1). Across the genus most genomes share well under half of the
16M sequence, so the tree of this tier was built from core genes rather than
from SNPs against 16M: the 213 genomes were annotated with bakta 1.12.2 (light
database 6.0), the pangenome was built with PPanGGOLiN 2.3.2, and the 300 most
informative of 1,392 core genes (614,278 alignment columns) were aligned. A
tree of each gene was built with IQ-TREE 3.1.2 (GTR+G, `--fast`, aBayes
support), and the species tree was summarised from the 300 gene trees with
weighted ASTRAL (ASTER 1.25). Nodes with a local posterior probability below
0.9 were collapsed (64 of 210). The canSNPs of this tier were taken from an
alignment of the 213 genomes against 16M (minimap2 2.31, bcftools 1.23.1),
keeping columns called in at least half of the genomes.

**Classical tier.** The 1,337 genomes of `s__Brucella melitensis` were
dereplicated at 99.99% identity, keeping the 23 biovar and species reference
strains as representatives, which left 197 representatives. The outgroup is
*B. teleogrylli* (GCF_006376685.1). They were aligned against 16M as above,
keeping columns called in at least 95% of the genomes (419,842 variable
columns), and the tree was built with IQ-TREE 3.1.2 (GTR+G, 1,000 ultrafast
bootstrap replicates). Nodes with less than 90% bootstrap support were
collapsed (15 of 197). The classical genomes form one clade of the genus tree,
which `repgenr graft` checked before the two tiers were joined.

**canSNPs.** In each tier, SNPs were extracted from the alignment
(`--min-called 0.9`), scored by CanSNPCreator against the 16M annotation and
the tandem repeats found by Tandem Repeats Finder 4.10.0-rc.2, and up to three
canSNPs per node were exported, preferring canSNPs called in more genomes,
then the higher score, at least 1,000 bases apart (`--snps-per-node 3
--prune`). Nodes without a canSNP were left out, and their children hang from
the nearest node that has one (11 in the genus tier, 10 in the classical
tier). `CanSNPCreator merge` joined the tiers at `BR.classical`: the 23 genus
tier nodes below it were replaced by the 171 of the classical tier, and the
canSNPs of `BR.classical` itself are those of the genus tier. Of the 809
canSNPs, 649 lie on chromosome I and 160 on chromosome II; 257 nodes have three
canSNPs, 12 have two and 14 have one. The first canSNP of each node is marked
`typed`. The `reference` column of the canSNP file is `FOIunpublished`, as the
canSNPs are not taken from a publication.

### Genomes

`brucella_genus_genomes.txt` lists the genomes of the 102 terminal nodes: 221
representatives, 975 genomes folded into them by dereplication (`member`), and
the two outgroups. Genomes that sit at an inner node are not listed; they were
used to find the canSNPs.

## Validation

`validate.py` and `check_bases.py` pass: the derived or ancestral base of all
809 canSNPs is the base of 16M at its position.

### Typing the genomes the scheme was made from

`CanSNPCreator audit` typed the representatives of each tier from its own
alignment. All were called at their node when the canSNPs of their path are
taken together. By the rule of CanSNPer2, one canSNP per node, 178 of the 182
representatives of the genus tier and 182 of the 184 of the classical tier
were called at their node; the others were called one node higher or not
called.

### The genus of GTDB r232

All 1,785 genomes were typed with CanSNPer2 2.4.0 (`--aligner minimap2`) on
2026-11-24, with this scheme and with `brucella_classical` for comparison.

| | `brucella_genus`, default (`--min-required-hits 3`) | `brucella_genus`, `--min-required-hits 2` | `brucella_classical`, default | `brucella_classical`, `--min-required-hits 2` |
|---|---|---|---|---|
| the species cluster, 1,337 genomes | 1,337 called | 1,337 called | 304 called | 1,237 called |
| other *Brucella* species, 448 genomes | 440 called | 443 called | none called | none called |

No genome of the species cluster was called outside `BR.classical`, and no
other genome was called inside it.

Each genome was compared with the node the scheme places it at, which for a
genome folded into a representative is the representative's node. At the
default threshold:

| | Genomes | At its node | Above it | Below it | On another branch | No call |
|---|---|---|---|---|---|---|
| representatives, species cluster | 197 | 195 | 1 | 1 | 0 | 0 |
| representatives, other species | 172 | 167 | 3 | 0 | 0 | 2 |
| members, species cluster | 1,140 | 547 | 319 | 35 | 239 | 0 |
| members, other species | 276 | 232 | 23 | 1 | 14 | 6 |

Most calls of members on another branch are within *B. melitensis*: for 158 of
the 239 the call and the representative's node part below
`BR.classical.melitensis_bv2.1.1.1.1`. A member was folded into its
representative at 99.99% identity, which still allows some hundred
differences, so its representative's node is an approximation of where it
belongs. These calls have not been checked against the members' own sequence.

The named node reached by the genomes of the species cluster, against the
species in their NCBI record, at the default threshold (the nearest named node
on the path of the call):

| NCBI species | Genomes | Named node of the call |
|---|---|---|
| *B. melitensis* | 701 | 462 `melitensis_bv2`, 177 `melitensis_bv3`, 58 `melitensis_bv1`, 1 `abortus`, 1 `suis_bv1`, 1 `suis_bv3`, 1 `BR.classical` |
| *B. abortus* | 378 | 303 `abortus`, 39 `abortus_bv1`, 17 `abortus_bv3`, 4 `abortus_bv9`, 2 `abortus_bv7`, 10 `melitensis_bv2`, 1 `melitensis`, 2 `suis_bv1` |
| *B. suis* | 97 | 46 `suis_bv2`, 35 `suis_bv1`, 6 `suis_bv4`, 3 `suis_bv3`, 1 `inopinata`, 1 `melitensis_bv1`, 5 `BR.classical` |
| *B. canis* | 32 | 32 `canis` |
| *B. ceti* | 24 | 15 `ceti`, 9 `BR.classical` |
| *B. ovis* | 18 | 18 `ovis` |
| *B. pinnipedialis* | 8 | 8 `pinnipedialis` |
| *B. neotomae* | 6 | 6 `neotomae` |
| *B. microti* | 5 | 5 `microti` |
| *B. inopinata* | 4 | 3 `inopinata`, 1 `suis_bv4` |
| *B. vulpis* | 2 | 2 `BR.classical` |
| *B. amazoniensis* | 2 | 2 `neotomae` |
| recorded only as *Brucella* sp. | 60 | 24 `inopinata`, 18 `melitensis_bv2`, 4 `ceti`, 4 `pinnipedialis`, 2 `neotomae`, 1 `ovis`, 1 `suis_bv3`, 1 `abortus_bv1`, 5 `BR.classical` |

The ten genomes filed as *B. abortus* that are called in `melitensis_bv2` are
also called on the *B. melitensis* branch (`cysW_288`) by `brucella_classical`
at `--min-required-hits 2`; its README discusses the genomes called on the
branch of another species than their NCBI record names. Of the 23 reference
strains, 18 were called at the node named after their label. *B. abortus*
biovars 2, 4, 5 and 6 were called at `BR.classical.abortus` and *B. suis*
biovar 5 at `BR.classical`; see Known issues.

The scheme has not been checked against published typing of these genomes, or
against genomes outside GTDB r232.

## Known issues

- 26 nodes carry fewer than three canSNPs (14 one, 12 two), most of them deep
  in *B. melitensis* biovar 2 and in *B. pseudogrignonensis*. A genome is
  called at such a node with the default `--min-required-hits 3` only when
  enough canSNPs of the path above it are derived.
- Labels that name no node, because their reference strains share every node
  with a strain of another label: *B. abortus* biovars 4, 5 and 6 and *B. suis*
  biovar 5 in the classical tier, and 17 GTDB species in the genus tier,
  among them *B. anthropi*-A, *B. daejeonensis*, *B. grignonensis* and
  *B. rhizosphaerae*. The node of *B. abortus* biovar 2 had no canSNP and was
  left out, so its reference strain is called at `BR.classical.abortus`.
- *B. suis* is paraphyletic with respect to *B. canis*, so no node is named
  `suis`.
- Eight of the 448 genomes outside the species cluster are not called at the
  default threshold, three of them below `BR.ovis-A` and three below
  `BR.anthropi.1.1.1.1`; five remain uncalled at `--min-required-hits 2`.
- The genus tier was built from 300 core genes and its canSNPs from an
  alignment against 16M, which covers less than half of most genomes of the
  former *Ochrobactrum*. Eleven nodes of the species tree had no canSNP in
  that alignment and were left out.

## References

None; the canSNPs are not taken from a publication. CanSNPCreator:
https://github.com/FOI-Bioinformatics/CanSNPCreator; RepGenR2:
https://github.com/FOI-Bioinformatics/repgenr2

## Changes

- 2026-11-24: first version.
