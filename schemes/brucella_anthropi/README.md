# Brucella anthropi scheme (with B. cytisi)

Made with CanSNPCreator 0.5.1 from the RefSeq assemblies of *Brucella
anthropi* (formerly *Ochrobactrum anthropi*, NCBI taxon 529) and *Brucella
cytisi* (taxon 407152). Three assemblies filed as *B. anthropi* were left out
because their average nucleotide identity (ANI) to the type strain placed them
in another species, which leaves 76 genomes. The reference is the *B. anthropi*
type strain ATCC 49188 (GCF_000017405.1). The outgroup is *B. tritici*
(GCF_008932295.1), the closest species outside the two by ANI (91.5-91.8% to
both type strains). Fifty-seven representatives were chosen by RepGenR2
dereplication, and the tree over them was built with IQ-TREE (GTR+F+I+R4,
1,000 ultrafast bootstrap replicates). Nodes with less than 90% bootstrap
support were collapsed. The scheme keeps the 42 nodes that have a canSNP, with
up to five ranked canSNPs per node in `brucella_anthropi_canSNPs.txt`. The
`reference` column of the SNP file is `NA`, as the canSNPs are not taken from a
publication.

In this tree, *B. cytisi* and *B. anthropi* are sister groups, and *B.
anthropi* divides into two lineages, each with 100% bootstrap support. The
nodes are named after them: `Bant.cytisi`; `Bant.anthropi`;
`Bant.anthropi.T`, the lineage of the type strain (26 genomes, 12 of the 17
with a recorded source clinical); `Bant.anthropi.T.E`, a group of five genomes
within it, all from Europe; and `Bant.anthropi.env` (43 genomes, 22 of the 31
with a recorded source environmental). Nodes below a named node keep the rest
of their code. The names describe the genomes at hand and are not established
designations. One genome filed as *B. cytisi* (WTB2) falls within
*B. anthropi*, in agreement with its ANI.

With CanSNPer2 2.2.2 (`--aligner minimap2`), 75 of the 76 genomes were called
at the node the CanSNPCreator typing report gives them, or at its nearest
ancestor in the scheme where that node has no canSNP; this includes all 19
genomes not used as representatives. One *B. cytisi* representative
(GCF_058245555.1) sits directly at `Bant.cytisi` and carries only two derived
canSNPs on its path, `Bant` and `Bant.cytisi`, so it is not called with the
default `--min_required_hits 3`. The scheme has not been checked against
published typing of these species.
