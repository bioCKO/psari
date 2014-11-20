#!/bin/bash -x

#SBATCH -A b2010042
#SBATCH -o HC-overlap.out 
#SBATCH -e HC-overlap.err
#SBATCH -J HC-overlap.job
#SBATCH -p core
#SBATCH -t 3:00:00
#SBATCH --mail-user douglas.scofield@ebc.uu.se
#SBATCH --mail-type=ALL

GFFs=("$@")


TmpID=${SLURM_JOB_ID:-$$}

HCGFF=sprucev01_HQ.gff

for gff in ${GFFs[@]} ; do
    Regions=${gff%.gff}.pseudogenic_region.gff
    grep pseudogenic_region $gff > $Regions
    HCRegions=${gff%.gff}.HC-overlaps.txt
    bedtools intersect -u -a $Regions -b $HCGFF | sed -e 's/^.*ID=\([^;]\+\);.*$/\1[:;]/g' > $HCRegions
    chmod 440 $HCRegions
    rm -f $Regions  # no longer need these, since we have the overlapping IDs
    # now use these IDs to remove matching pseudogenic_region and pseudogenic_block 
    Screened=${gff%.gff}.HC-screened.gff
    grep -v -f $HCRegions $gff > $Screened
    chmod 440 $Screened
done
