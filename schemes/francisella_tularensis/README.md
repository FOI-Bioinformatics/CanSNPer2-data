# Francisella tularensis canSNP scheme

| | |
|---|---|
| Scheme | `francisella_tularensis` |
| Status | verified |
| Organisms | *Francisella tularensis* (NCBI taxon 263) |
| Reference genome | OSU18 (GCF_000014605.1), FSC200 (GCF_000168775.2), FTNF002-00 (GCF_000017785.1), LVS (GCF_000009245.1), SCHU S4 (GCF_000008985.1) |
| Nodes | 328 |
| canSNPs per node | 1 |
| Origin | published canSNPs and canSNPs from FOI, compiled at FOI |
| Last changed | 2026-11-22 |

## Description

A hierarchy of canSNPs for *F. tularensis*, from the root `T/N.1` to the
subspecies and their lineages: `A.I` and `A.II` (subspecies *tularensis*),
`M.1` (subspecies *mediasiatica*) and the `B.` nodes (subspecies *holarctica*),
which make up 301 of the 328 nodes. Nodes carry the names used in the
publications below. The canSNP positions are given on five reference genomes,
each chosen to be close to the lineages whose canSNPs it carries.

### Node names

The next free node name in subspecies *holarctica* is `B.317`; `B.316` is
reserved.

## Origin

The scheme was compiled at FOI from published canSNPs and canSNPs that FOI has
not published. The `reference` column of the canSNP file gives the source of each
canSNP:

| Source | canSNPs |
|---|---|
| FOI, unpublished | 136 |
| Kevin2020 | 83 |
| Vogler2009 | 22 |
| Dwibedi2016 | 21 |
| Birdsell2014 | 15 |
| Wittwer2018 | 11 |
| Svensson2009 | 8 |
| Chanturia2011 | 7 |
| Gyuranecz2012 | 6 |
| Schulze2016 | 5 |
| Myrtennas2016 | 4 |
| Dwibedi2017 | 3 |
| Karadenizli2015 | 2 |
| Karlsson2013 | 2 |
| Birdsell2019 | 1 |
| Busch2016 | 1 |
| Sissonen2015 | 1 |

### Published canSNPs not included

These published canSNPs are not in the scheme. Most identify the same branch
as a canSNP that is included; for the others, the reason was not recorded.

| canSNP | Source | Identifies |
|---|---|---|
| `B.14` | Vogler2009 | the same branch as `B.24` |
| `B.15` | Svensson2009 | the same branch as `B.1` and Ftind47 |
| `B.17` | Svensson2009 | the same branch as `B.4` and Ftind38 |
| `B.18` | Svensson2009 | the same branch as `B.11` |
| `B.19` | Svensson2009 | the same branch as `B.12` and Ftind33 |
| `B.20` | Svensson2009 | the same branch as `B.42` |
| `B.40` | Karlsson2013 | the same branch as `B.7` |
| `B.41` | Karlsson2013 | the same branch as `B.10` |
| `B.43` | Karlsson2013 | the same branch as `B.27` |
| `B.69` | Sissonen2015 | the same branch as `B.11` |
| `B.310` | Kevin2020 | the same branch as `B.215` |
| `B.311` | Kevin2020 | the same branch as `B.168` |
| `B.312` | Kevin2020 | the same branch as `B.66` |
| `B.313` | Kevin2020 | the same branch as `B.190` |
| `A.2` | Svensson2009 | not recorded |
| `A.3` | Svensson2009 | not recorded |
| `A.4` | Svensson2009 | not recorded |
| `A.5` | Svensson2009 | not recorded |
| `A.6` | Svensson2009 | not recorded |
| `A.II.5` | Vogler2009 | not recorded |
| `A.II.7` | Vogler2009 | not recorded |

### Genomes

The genomes file `francisella_tularensis_genomes.txt` lists 219 assemblies at 109 of the 161 leaves. 180 are `published`, placed at the node by a publication (Kevin2020 144, Wittwer2018 13, Dwibedi2016 11, Karadenizli2015 6, scheme reference 3, Schulze2016 3), and 39 are `typed`, the genome of the GTDB r232 species cluster that CanSNPer2 calls at a leaf no publication covers (chosen as described in `schemes/README.md`).

- Each published genome was typed with CanSNPer2 2.3.0 (minimap2) and is listed only where the call is its node or a node below it. Kevin2020 genomes left out: five called at a shallower node (B.232 twice, B.284, B.292, B.300), one called on another lineage (B.237 called B.264), one called at an inner node (B.88 called B.89), and five at nodes B.310, B.312 and B.313, which are not in the scheme tree.
- The labels B.88 to B.91 of Wittwer2018 (Table 1) are one step lower than the scheme's: its genomes of B.88, B.89, B.90 and B.91 carry the scheme's SNPs B.89, B.90, B.91 and B.92. The text of that paper and the scheme list the same five subclades (B.87/B.88, B.89, B.90, B.91, B.92). The genomes are listed at the scheme's nodes.
- 52 leaves have no genome: A.I.6, A.I.16, B.32, B.38, B.96, B.101, B.105, B.123, B.132, B.134, B.141, B.142, B.143, B.144, B.145, B.148, B.149, B.152, B.154, B.157, B.158, B.161, B.167, B.169, B.172, B.173, B.174, B.175, B.176, B.178, B.181, B.182, B.183, B.184, B.185, B.186, B.187, B.189, B.194, B.198, B.200, B.204, B.206, B.208, B.209, B.211, B.214, B.216, B.221, B.225, B.227, B.232.

## Validation

`validate.py` checks the files, and `build.py` rebuilds the database and
compares it with the committed one (CanSNPer2 2.2.0). When the build tooling
was added (#13), the database rebuilt from these files was identical in nodes,
edges, SNP annotations and references to the one committed before.

Typed with CanSNPer2 2.2.2 (`--aligner minimap2`, default settings) on
2026-11-19:

- **One representative per species of the genus.** Of the 26 *Francisella*
  species representatives of GTDB r232, only *F. tularensis* NIH B-38, a
  subsp. *tularensis* strain, was called, at `A.II.4`. The other 25 species,
  among them *F. hispaniensis*, *F. philomiragia* and *F. persica*, received
  no call.
- **The species cluster.** All 937 genomes of the GTDB r232 species
  `s__Francisella tularensis` were typed, and the call of each was compared
  with the subspecies in its NCBI record:

| Subspecies in the NCBI record | Genomes | Called on the `A` branch | on `B` | at `M.1` or below | No call |
|---|---|---|---|---|---|
| subsp. *tularensis* | 48 | 48 | 0 | 0 | 0 |
| subsp. *holarctica* | 309 | 0 | 306 | 0 | 3 |
| subsp. *mediasiatica* | 16 | 0 | 0 | 16 | 0 |
| subsp. *novicida* | 42 | 0 | 0 | 0 | 42 |
| no subspecies recorded | 521 | 38 | 477 | 1 | 5 |
| recorded only as *Francisella* | 1 | 0 | 0 | 0 | 1 |

  No genome was called on a branch that contradicts its recorded subspecies.
  The subspecies *novicida* has no branch in the scheme, so no call is the
  expected result for those 42 genomes.

In November 2026, the same 937 genomes were also typed with the Mauve,
minimap2 and MUMmer4 aligners, and each genome received the same final call
with all three, so the calls do not depend on the aligner.

The calls below subspecies level have not been compared with an independent
typing of these genomes.

## Known issues

- The full reference for one published source, Dwibedi2017 (`B.65`, `B.66` and
  `B.81`), is not yet recorded below. These canSNPs are not in Dwibedi2016
  (`B.44` to `B.64`), and no further paper of that author on *Francisella*
  was found. Busch2016, Karadenizli2015 and Kevin2020 were added on
  2026-11-21 and should be confirmed. Myrtennas2016 was added on the same
  date; its Table 3 gives `B.76` to `B.79` at the SCHUS4.2 positions of this
  scheme.

## References

References marked * were matched to the source key by first author, year and
subject and should be confirmed.

- Birdsell2014: Birdsell DN, Yaglom H, Rodriguez E, et al. *Francisella
  tularensis* subsp. *tularensis* group A.I, United States. Emerg Infect Dis.
  2014;20(5):861-865. https://doi.org/10.3201/eid2005.131559
- Birdsell2019*: Phylogenetic analysis of *Francisella tularensis* group A.II
  isolates from 5 patients with tularemia, Arizona, USA, 2015-2017. Emerg
  Infect Dis. 2019;25(5). https://doi.org/10.3201/eid2505.180363
- Busch2016*: Busch A, Thomas P, Myrtennäs K, et al. High-quality draft genome
  sequence of *Francisella tularensis* subsp. *holarctica* strain 08T0073
  isolated from a wild European hare. Genome Announc. 2017;5(12):e01577-16.
  https://doi.org/10.1128/genomea.01577-16
- Chanturia2011: Chanturia G, Birdsell DN, Kekelidze M, et al.
  Phylogeography of *Francisella tularensis* subspecies *holarctica* from the
  country of Georgia. BMC Microbiol. 2011;11:139. https://doi.org/10.1186/1471-2180-11-139
- Dwibedi2016: Dwibedi C, Birdsell D, Lärkeryd A, et al. Long-range dispersal
  moved *Francisella tularensis* into Western Europe from the East. Microb
  Genom. 2016;2(12):e000100. https://doi.org/10.1099/mgen.0.000100
- Gyuranecz2012: Gyuranecz M, Birdsell DN, Splettstoesser W, et al.
  Phylogeography of *Francisella tularensis* subsp. *holarctica*, Europe.
  Emerg Infect Dis. 2012;18(2):290-293. https://doi.org/10.3201/eid1802.111305
- Karadenizli2015*: Karadenizli A, Forsman M, Şimşek H, et al. Genomic analyses
  of *Francisella tularensis* strains confirm disease transmission from
  drinking water sources, Turkey, 2008, 2009 and 2012. Euro Surveill.
  2015;20(21):21136. https://doi.org/10.2807/1560-7917.es2015.20.21.21136
- Karlsson2013: Karlsson E, Svensson K, Lindgren P, et al. The phylogeographic
  pattern of *Francisella tularensis* in Sweden indicates a Scandinavian origin
  of Eurosiberian tularaemia. Environ Microbiol. 2013;15(2):634-645.
  https://doi.org/10.1111/1462-2920.12052
- Kevin2020*: Kevin M, Girault G, Caspar Y, et al. Phylogeography and genetic
  diversity of *Francisella tularensis* subsp. *holarctica* in France
  (1947-2018). Front Microbiol. 2020;11:287. https://doi.org/10.3389/fmicb.2020.00287
- Myrtennas2016: Myrtennäs K, Marinov K, Johansson A, et al. Introduction and
  persistence of tularemia in Bulgaria. Infect Ecol Epidemiol. 2016;6:32838.
  https://doi.org/10.3402/iee.v6.32838
- Schulze2016: Schulze C, Heuner K, Myrtennäs K, et al. High and novel genetic
  diversity of *Francisella tularensis* in Germany and indication of
  environmental persistence. Epidemiol Infect. 2016;144(14):3025-3036.
  https://doi.org/10.1017/S0950268816001175
- Sissonen2015*: Phylogeography of *Francisella tularensis* subspecies
  *holarctica* in Finland, 1993-2011. Infect Dis (Lond). 2015;47(10):701-706.
  https://doi.org/10.3109/23744235.2015.1049657
- Svensson2009: Svensson K, Granberg M, Karlsson L, et al. A real-time PCR
  array for hierarchical identification of *Francisella* isolates. PLoS ONE.
  2009;4(12):e8360. https://doi.org/10.1371/journal.pone.0008360
- Vogler2009: Vogler AJ, Birdsell D, Price LB, et al. Phylogeography of
  *Francisella tularensis*: global expansion of a highly fit clone. J
  Bacteriol. 2009;191(8):2474-2484. https://doi.org/10.1128/JB.01786-08
- Wittwer2018: Wittwer M, Altpeter E, Pilo P, et al. Population genomics of
  *Francisella tularensis* subsp. *holarctica* and its implication on the
  eco-epidemiology of tularemia in Switzerland. Front Cell Infect Microbiol.
  2018;8:89. https://doi.org/10.3389/fcimb.2018.00089

## Changes

- 2020-04-29: added (#1, #2).
- 2020-05-07: base errors corrected; B.314 added.
- 2020-06-25: A.II.3 and A.II.4 corrected; A.II.8 added (#5).
- 2020-08-27: B.315 added (#8).
- 2026-11-04: files normalised to LF line endings without trailing whitespace.
- 2026-11-11: database rebuilt with CanSNPer2 2.2.0 (#13).
- 2026-11-19: files moved to `schemes/francisella_tularensis/` (#16); README
  added.
- 2026-11-19: next free node name and the published canSNPs not included
  moved here from the repository README.
- 2026-11-19: typed against the GTDB r232 species representatives of the genus
  and the whole species cluster; the result is under Validation.
- 2026-11-21: SNP positions given per reference sequence (the chromosomes of the five reference genomes)
  with `CanSNPer2-database --convert`; database rebuilt with CanSNPer2
  2.3.0.
- 2026-11-21: added the references of Busch2016, Karadenizli2015, Kevin2020 and
  Myrtennas2016.
- 2026-11-22: added `francisella_tularensis_genomes.txt`: 219 genomes, 180 published and 39 typed, with the columns `basis` and `source`.
- 2026-11-22: `francisella_tularensis_canSNPs.txt` replaces the SNP file: the same 328
  canSNPs, each rank 1 of its node and typed, with the `reference` column kept and
  the `strain` column dropped.
