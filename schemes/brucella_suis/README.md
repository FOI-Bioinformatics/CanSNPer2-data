# Brucella suis scheme (with B. canis)

Made with CanSNPCreator 0.5.1 from the current RefSeq assemblies of
*Brucella suis* (NCBI taxon 29461, 120 assemblies) and *Brucella canis*
(taxon 36855, 75 assemblies), with *B. melitensis* 16M (GCF_000007125.1) as the
outgroup and *B. suis* 1330 (GCF_000007505.1) as the reference. *B. canis* is
included because it lies within *B. suis*. Forty representatives were chosen
from a tree over all 195 genomes, and the scheme keeps the 31 nodes that have a
canSNP, with up to five ranked canSNPs per node in `brucella_suis_canSNPs.txt`.
The `reference` column of the SNP file is `NA`, as the canSNPs are not taken
from a publication.

Nodes are named after the biovar or species that the genomes below them carry
in their NCBI BioSample records (`Bsu.bv1` to `Bsu.bv5`, `Bsu.canis`); nodes
below a named node keep the rest of their code, and nodes without a clear label
keep their code (`Bsu.1.1`). The names rest on the 54 genomes with a recorded
biovar and should be read as indicative.

With CanSNPer2 2.2.2 (`--aligner minimap2`), the 155 genomes not used as
representatives were each called at the node the CanSNPCreator typing report
gives them. Against a tree over all genomes, that node agreed with the
genome's place in the tree for 65.8% of them, was deeper for 21.9% and higher
for 11.6%. The scheme has not been checked against published *B. suis* typing.
