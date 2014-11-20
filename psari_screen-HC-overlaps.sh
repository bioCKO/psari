#!/bin/bash

#SBATCH -A b2010042
#SBATCH -o HC-overlap.out 
#SBATCH -e HC-overlap.err
#SBATCH -J HC-overlap.job
#SBATCH -p core
#SBATCH -t 3:00:00
#SBATCH --mail-user douglas.scofield@ebc.uu.se
#SBATCH --mail-type=ALL

GFFs=("$@")

#  NOTE:  I worry that a pseudogenic_block may not overlap with
#         a high-confidence region while its parent pseudogenic_region
#         does.  See region-HC-overlap.sh for a version of this that
#         looks for overlaps in regions and removed regions and blocks
#         belonging to overlapping regions.  The grep within that takes
#         a lot of memory to run so it may not fit on a node.  This
#         approach is much much faster and probably good enough for
#         downstream work that primarily focuses on regions.

TmpID=${SLURM_JOB_ID:-$$}

HCGFF=sprucev01_HQ.gff

Stats=$HCGFF.overlap-stats.txt
cat /dev/null > $Stats

for gff in ${GFFs[@]} ; do
    TotILines=$(grep -v '^#' $gff | wc -l | cut -f1 -d" ")
    CandILines=$(grep -v '^#' $gff | grep 'pseudogenic_region' | wc -l | cut -f1 -d" ")
    echo -n "working $gff, which has $TotILines input lines..."
    HCOverlap=${gff%.gff}.HC-overlap.gff
    bedtools intersect -u -a $gff -b $HCGFF > $HCOverlap
    chmod 440 $HCOverlap
    TotOLines=$(wc -l $HCOverlap | cut -f1 -d" ")
    CandOLines=$(grep 'pseudogenic_region' $HCOverlap | wc -l | cut -f1 -d" ")
    Screened=${gff%.gff}.HC-screened.gff
    bedtools intersect -v -a $gff -b $HCGFF > $Screened
    chmod 440 $Screened
    TotSLines=$(wc -l $Screened | cut -f1 -d" ")
    CandSLines=$(grep 'pseudogenic_region' $Screened | wc -l | cut -f1 -d" ")
    echo " total overlap $TotOLines, screened $TotSLines; candidate overlap $CandOLines, screened $CandSLines"
    echo "$gff totalinput $TotILines" >> $Stats
    echo "$gff totaloverlap $TotOLines" >> $Stats
    echo "$gff totalscreened $TotSLines" >> $Stats
    echo "$gff candidateinput $CandILines" >> $Stats
    echo "$gff candidateoverlap $CandOLines" >> $Stats
    echo "$gff candidatescreened $CandSLines" >> $Stats
done
