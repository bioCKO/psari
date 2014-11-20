# Psari: Pseudogene discovery

Psari is a pipeline that processes tblastn output to produce a GFF file describing pseudogene candidates in a genome.  Query sequences used for the tblastn are amino acid sequences from high-confidence genes *within the target genome*, while the subject sequence is the genome itself.  The basic steps of the pipeline are:

1. Discover candidate pseudogenes using tblastn of the high-confidence coding regions to the whole-genome sequence (script `psari_tblastn.sh`)
2. Tile the tblastn hits and create a GFF describing candidate pseudogenic regions; candidates are annotated in the GFF with similarity to the original gene, etc. (script `psari.pl`)
3. Screen out pseudogenes that overlap with high-confidence genes (script `psari_screen-HC-overlaps.sh`, with alternate script `psari_screen-HC-overlaps_alternate.sh` that may be a bit more complete but requires more memory)
4. Annotate the GFF entries for each screened pseudogene candidate with the occurrence of SNV and indel variants (script `psari_annotateBedWithVariants.sh`)

The script `run_psari.sh` gives an overview of running the pipeline, once the tblastn results have been gathered.

Psari requires BioPerl, BEDtools and awk.

Psari is currently very much tied to the Norway spruce genome project, but I hope to generalise the interface in the near future.
