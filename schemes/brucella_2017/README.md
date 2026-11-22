# Brucella genus canSNP scheme (2017, revised 2023)

| | |
|---|---|
| Scheme | `brucella_2017` |
| Status | draft |
| Organisms | *Brucella* (NCBI taxon 234), the genus including the atypical species |
| Reference genome | 2308 (GCF_000054005.1) |
| Nodes | 43 |
| canSNPs per node | up to 3 |
| Origin | tree and genome assignment from FOI, unpublished; canSNPs made with CanSNPCreator |
| Last changed | 2026-11-22 |

## Description
A canSNP hierarchy for the genus *Brucella*, from the root `B.1`, which every
*Brucella* carries, down to the biovars of the classical species. `B.3` holds
the early branching lineages (`B.4`, `I.1`), and `B.2` the rest, which divides
into the atypical and classical groups: `A.` for *B. abortus*, `Me.` for *B.
melitensis*, `S.` and `Ca.` for *B. suis* and *B. canis*, `Ce.` and `Pi.` for
*B. ceti* and *B. pinnipedialis*, `O.` and `Pa.` for *B. ovis* and *B.
papionis*, `N.` for *B. neotomae*, `Mi.` for *B. microti* and `V.` for *B.
vulpis*. The node names are those of the 2017 work and are kept unchanged.

### Node names
The tree is the 2023 revision of the 2017 hierarchy. It adds `B.2a`, `B.2b` and
`B.4` and drops `Bsp.1`, `IBsp.1`, `S.3` and `SCa.5` of the 2017 version, whose
`V.1` and `IBsp.1` each had two parents, which no scheme can have.

## Origin
The tree and the placement of the genomes are unpublished FOI work: a hierarchy
drawn up in 2017 and revised in September 2023, with 57 genomes assigned to its
nodes.

The canSNPs are not those of the 2017 work. That selection, a table of 42
positions, is given in coordinates that match no assembly of *Brucella*: of
1,134 genomes tried, the best agreed with the alleles at 55% of the sites,
which is chance. The positions are columns of an alignment that was not kept,
so they cannot be placed on a reference sequence, which is what a scheme needs.

The canSNPs here were made again from the same genomes with CanSNPCreator
0.7.0: the 56 *Brucella* genomes of the 2023 assignment were aligned to 2308
with RepGenR2's `simple` typer (minimap2, projected onto the reference), and
*Ochrobactrum teleogrylli* (GCF_006376685.1) was aligned as the outgroup, which
is what gives the root canSNP its polarity. A site is used where at least 90%
of the genomes have a base (`extract --min-called 0.9`); the outgroup covers
only 42% of the reference, so a strict core would have left `B.2a`, `B.2b` and
`B.3` without a canSNP. Of 387,212 core SNPs, 195,490 are canonical for a node,
and the three best of each node are listed, scored against the reference's
genes (the NCBI annotation of GCA_000054005.1) and tandem repeats (TRF).

### Genomes
`brucella_2017_genomes.txt` lists a genome of each of the 24 terminal nodes, 52
in all, and *O. teleogrylli* as the `outgroup`. Four genomes of the 2023
assignment sit at an inner node and are not listed: `GCF_000177155.1`,
`GCF_002994775.1` and `GCF_024970695.1` at `B.2`, and `GCF_000292145.1` at
`Me.5`. They were used to find the canSNPs.

## Validation
Partly validated. The scheme passes `validate.py` and `check_bases.py` of this
repository. The 56 *Brucella* genomes of the 2023 assignment were typed against
it with CanSNPer2 2.4.0 and minimap2: 53 are called at the node they were
assigned to, and the outgroup is `NA`, as it should be. The three that are not
are the `B.2` genomes, which are called `NA` at the default
`--min_required_hits 3` because `B.2` lies two nodes from the root; with
`--min_required_hits 2` all three are called `B.2`.

Before it is moved to master it needs what the other schemes have: a run
against the GTDB species cluster of the genus, and genomes typed that were not
used to find the canSNPs.

## Known issues
- The canSNPs do not come from the 2017 table, whose coordinates could not be
  placed on a reference (see Origin). A canSNP of this scheme and a canSNP of
  the 2017 work at the same node are not the same position.
- `B.2a` and `B.2b` have a canSNP only because the core is relaxed to 90% of
  the genomes; with a strict core neither can be typed.
- The 2023 assignment places four genomes at an inner node, so those nodes have
  no genome of their own in the genomes file.
- The root `B.1` is the node every *Brucella* carries. Its canSNP rests on one
  outgroup genome, *O. teleogrylli*, which covers 42% of the reference.
- `B.2` and `B.3` lie fewer than three nodes from the root, so CanSNPer2's
  default `--min_required_hits 3` returns `NA` for a genome that belongs to
  either. `--min_required_hits 2` is needed to call them. The same applies to
  the other schemes whose shallowest nodes sit close to the root.

## References
- None. The tree, the genome assignment and the 2017 canSNP table are
  unpublished FOI work, recorded as `FOIunpublished` in the canSNP file.

## Changes
- 2026-11-22: added, as a draft: the 2023 tree and genome assignment with
  canSNPs made again with CanSNPCreator 0.7.0 against 2308. The genomes of the
  scheme were typed against it; see Validation.
