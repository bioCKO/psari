## Psari: Pseudogene discovery

Psari is (for the moment) a single perl script that processes tblastn output to produce a GFF file describing pseudogene candidates in a genome.  Query sequences used for the tblastn are amino acid sequences from high-confidence genes *within the target genome*, while the subject sequence is the genome itself.  Psari uses BioPerl to parse tblastn output, and does its own tiling of blast hits to discover candidates.

Psari is currently very much tied to the Norway spruce genome project, but I will be generalising the interface in the near future.

