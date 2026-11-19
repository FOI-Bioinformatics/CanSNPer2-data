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
| Last changed | 2026-11-04 |

## Description

A hierarchy of canSNPs for *F. tularensis*, from the root `T/N.1` to the
subspecies and their lineages: `A.I` and `A.II` (subspecies *tularensis*),
`M.1` (subspecies *mediasiatica*) and the `B.` nodes (subspecies *holarctica*),
which make up 301 of the 328 nodes. Nodes carry the names used in the
publications below. The canSNP positions are given on five reference genomes,
each chosen to be close to the lineages whose canSNPs it carries.

## Origin

The scheme was compiled at FOI from published canSNPs and canSNPs that FOI has
not published. The `reference` column of the SNP file gives the source of each
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

## Validation

`validate.py` checks the files, and `build.py` rebuilds the database and
compares it with the committed one (CanSNPer2 2.2.0). When the build tooling
was added (#13), the database rebuilt from these files was identical in nodes,
edges, SNP annotations and references to the one committed before.

In November 2026, 937 *F. tularensis* genomes were typed with CanSNPer2 using
the Mauve, minimap2 and MUMmer4 aligners, and each genome received the same
final call with all three. This shows that the calls do not depend on the
aligner; it does not show that they are correct. The calls have not been
compared here with an independent typing of those genomes.

## Known issues

- The full references for six sources (Busch2016, Dwibedi2017,
  Karadenizli2015, Kevin2020, Myrtennas2016 and the FOI canSNPs) are not yet
  recorded below.

## References

References marked * were matched to the source key by first author, year and
subject and should be confirmed.

- Birdsell2014: Birdsell DN, Yaglom H, Rodriguez E, et al. *Francisella
  tularensis* subsp. *tularensis* group A.I, United States. Emerg Infect Dis.
  2014;20(5):861-865. doi:10.3201/eid2005.131559
- Birdsell2019*: Phylogenetic analysis of *Francisella tularensis* group A.II
  isolates from 5 patients with tularemia, Arizona, USA, 2015-2017. Emerg
  Infect Dis. 2019;25(5). doi:10.3201/eid2505.180363
- Chanturia2011: Chanturia G, Birdsell DN, Kekelidze M, et al.
  Phylogeography of *Francisella tularensis* subspecies *holarctica* from the
  country of Georgia. BMC Microbiol. 2011;11:139. doi:10.1186/1471-2180-11-139
- Dwibedi2016: Dwibedi C, Birdsell D, Lärkeryd A, et al. Long-range dispersal
  moved *Francisella tularensis* into Western Europe from the East. Microb
  Genom. 2016;2(12):e000100. doi:10.1099/mgen.0.000100
- Gyuranecz2012: Gyuranecz M, Birdsell DN, Splettstoesser W, et al.
  Phylogeography of *Francisella tularensis* subsp. *holarctica*, Europe.
  Emerg Infect Dis. 2012;18(2):290-293. doi:10.3201/eid1802.111305
- Karlsson2013: Karlsson E, Svensson K, Lindgren P, et al. The phylogeographic
  pattern of *Francisella tularensis* in Sweden indicates a Scandinavian origin
  of Eurosiberian tularaemia. Environ Microbiol. 2013;15(2):634-645.
  doi:10.1111/1462-2920.12052
- Schulze2016: Schulze C, Heuner K, Myrtennäs K, et al. High and novel genetic
  diversity of *Francisella tularensis* in Germany and indication of
  environmental persistence. Epidemiol Infect. 2016;144(14):3025-3036.
  doi:10.1017/S0950268816001175
- Sissonen2015*: Phylogeography of *Francisella tularensis* subspecies
  *holarctica* in Finland, 1993-2011. Infect Dis (Lond). 2015;47(10):701-706.
  doi:10.3109/23744235.2015.1049657
- Svensson2009: Svensson K, Granberg M, Karlsson L, et al. A real-time PCR
  array for hierarchical identification of *Francisella* isolates. PLoS ONE.
  2009;4(12):e8360. doi:10.1371/journal.pone.0008360
- Vogler2009: Vogler AJ, Birdsell D, Price LB, et al. Phylogeography of
  *Francisella tularensis*: global expansion of a highly fit clone. J
  Bacteriol. 2009;191(8):2474-2484. doi:10.1128/JB.01786-08
- Wittwer2018: Wittwer M, Altpeter E, Pilo P, et al. Population genomics of
  *Francisella tularensis* subsp. *holarctica* and its implication on the
  eco-epidemiology of tularemia in Switzerland. Front Cell Infect Microbiol.
  2018;8:89. doi:10.3389/fcimb.2018.00089

## Changes

- 2020-04-29: added (#1, #2).
- 2020-05-07: base errors corrected; B.314 added.
- 2020-06-25: A.II.3 and A.II.4 corrected; A.II.8 added (#5).
- 2020-08-27: B.315 added (#8).
- 2026-11-04: files normalised to LF line endings without trailing whitespace.
- 2026-11-11: database rebuilt with CanSNPer2 2.2.0 (#13).
- 2026-11-19: files moved to `schemes/francisella_tularensis/` (#16); README
  added.
