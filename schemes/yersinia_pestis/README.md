# Yersinia pestis canSNP scheme

| | |
|---|---|
| Scheme | `yersinia_pestis` |
| Status | verified |
| Organisms | *Yersinia pestis* (NCBI taxon 632) |
| Reference genome | CO92 (GCF_000009065.1) |
| Nodes | 97 |
| canSNPs per node | 1 |
| Origin | published (Morelli et al. 2010) |
| Last changed | 2026-11-21 |

## Description

The canSNP scheme of Morelli et al. (2010) for the global population of *Y.
pestis*. Nodes carry the published branch names (for example `0.PE4.a`,
`1.ORI2.b`, `2.MED1.a`); the Roman numerals (`III` to `XVII`) are the internal
nodes of that phylogeny, and the root is `root`.

## Origin

The canSNPs and their positions on the CO92 genome are taken from Morelli et
al. (2010). Of the 100 published canSNPs, 97 are kept; the three without a
derived base were removed on 2026-11-20 (see Known issues).

The genomes file `yersinia_pestis_genomes.txt` lists the genomes of Morelli et
al. (2010, Supplementary Table 1), each by its current assembly accession and
at the node that table gives it. Fifteen *Y. pestis* genomes are listed as
`representative` and *Y. pseudotuberculosis* IP32953 as the `outgroup`. Two of
the table's genomes are not listed: IP674 has no assembly, and FV-1 was left
out of the genome comparison the canSNPs came from, for suspected sequencing
errors. The node of `91001` is `0.PE4.c` as the table
gives it, although CanSNPer2 does not call it there because `91001` is
ancestral at the root canSNP (see Known issues).

## Validation

`validate.py` checks the files, and `build.py` rebuilds the database and
compares it with the committed one (CanSNPer2 2.2.0). When the build tooling
was added (#13), the database rebuilt from these files was identical in nodes,
edges, SNP annotations and references to the one committed before. Samples
simulated on the branches that were later removed were typed with CanSNPer2;
the result is described under Known issues.

Of the 17 genomes Morelli et al. list with their node (Supplementary Table 1),
16 are in GTDB r232 (IP674 is not; KIM is the KIM10+ assembly GCF_000006645.1,
chromosome NC_004088 as in the paper). Reading the bases at the canSNP positions
of these 16 places each at the node the paper gives it, with every
derived canSNP on that node's path; Pestoides F, node `0.PE2 a` in the paper, is
`0.PE2a.a` here (checked 2026-11-20).

Typed with CanSNPer2 2.2.2 (`--aligner minimap2`, default settings) on
2026-11-19:

- **One representative per species of the genus.** Of the 28 *Yersinia*
  species representatives of GTDB r232, only *Y. pestis* NCTC 5923 was
  called, at `1.ORI1.a`. The other 27 species received no call.
- **The species cluster.** All 1,039 genomes of the GTDB r232 species
  `s__Yersinia pestis` were typed. That cluster also holds
  *Y. pseudotuberculosis*, which the scheme separates from *Y. pestis*: of
  the 921 genomes that NCBI records as *Y. pestis*, 855 (93%) were called,
  and of the 113 recorded as *Y. pseudotuberculosis*, 111 received no call.
  The two that were called are 99.8% identical to CO92 by average nucleotide
  identity, so they are *Y. pestis* deposited under the wrong name. The four
  *Y. wautersii* genomes and one further genome received no call.

Of the 66 *Y. pestis* genomes without a call, 56 carry fewer than three
derived canSNPs, which the default `--min_required_hits 3` requires; see Known
issues for the rest. The calls have not been compared with an independent
typing of these genomes.

## Known issues

Three canSNPs, `1.ANT2.b` (CO92 position 4,191,752), `1.ORI2.c` (2,030,090)
and `1.ORI3.b` (434,961), were published with `/` instead of a derived base, and
were removed with their nodes on 2026-11-20; `1.ORI2.d`, the child of
`1.ORI2.c`, now descends directly from `1.ORI2.b`. CanSNPer2 compares an aligned
base against the ancestral and derived base, so the three could never be called
derived: in tests with simulated sequence, samples on these branches were
placed at a node above them, and a sample at `1.ORI2.d` was placed at
`1.ORI2.b`, without a warning. CanSNPer2 with source file checks
(FOI-Bioinformatics/CanSNPer2#66) refuses a scheme with a derived base other
than A, C, G or T.

The three positions lie at the start of runs of A in CO92 (`TCAGG[A]AAACT`,
`TTTTT[A]AAAAA`, `GCCCT[A]AAAAA`), so `/` most likely marks a one-base deletion
in the run rather than a substitution: the 933 markers Morelli et al. genotyped
for their minimum spanning tree include eight frameshifts (Supplementary Fig. 1
and Table 6 of the paper), and CanSNPer2 types substitutions only. The derived state could not be read from
genomes: of the 1,039 genomes of the GTDB r232 species cluster, none is derived
at `1.ORI2.d`, one is derived at `1.ORI2.b` and carries the CO92 base at
`1.ORI2.c`, and those derived at `1.ANT2.a` or `1.ORI3.d` all carry the CO92
base at the other two positions. Seven genomes lack one A at 434,961, one of them
also at 4,191,752, but these belong to other lineages (`2.MED`, `III`, `IX`), so
the deletions are not specific to the removed branches. No genome of the set is
derived at a node below a removed one, so its calls under Validation are the
same with the 97-node scheme. The nodes can be restored once their derived
states are taken from the original description of the SNPs.

Ten *Y. pestis* genomes of the species cluster are ancestral at the `root`
canSNP (CO92 position 257,340) although they are derived at `0.PE3.a` or
below. CanSNPer2 needs an unbroken path from the root, so these genomes
receive no call. The canSNP of the root node is therefore not derived in every
*Y. pestis*.

## References

- Morelli G, Song Y, Mazzoni CJ, et al. *Yersinia pestis* genome sequencing
  identifies patterns of global phylogenetic diversity. Nat Genet.
  2010;42(12):1140-1143. https://doi.org/10.1038/ng.705

## Changes

- 2023-03-13: added from the earlier CanSNPer data.
- 2026-11-04: files normalised to LF line endings without trailing whitespace.
- 2026-11-11: database rebuilt with CanSNPer2 2.2.0 (#13).
- 2026-11-19: files moved to `schemes/yersinia_pestis/` (#16); README added.
- 2026-11-19: typed against the GTDB r232 species representatives of the genus
  and the whole species cluster; the result is under Validation.
- 2026-11-20: removed `1.ANT2.b`, `1.ORI2.c` and `1.ORI3.b`, published without
  a derived base; `1.ORI2.d` now descends from `1.ORI2.b`.
- 2026-11-21: SNP positions given per reference sequence (AL590842.1)
  with `CanSNPer2-database --convert`; database rebuilt with CanSNPer2
  2.3.0.
- 2026-11-21: added `yersinia_pestis_genomes.txt`, the genomes of Morelli et
  al. (2010, Supplementary Table 1), with IP32953 as the outgroup.
