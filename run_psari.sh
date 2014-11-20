#! /bin/bash -l

echo -n "Time started: " 
date

set -e
set -x

#SBATCH -A b2010042
#SBATCH -o pseudo.out
#SBATCH -e pseudo.err
#SBATCH -J pseudo.job
#SBATCH -p core
#SBATCH -t 1-00:00:00
#SBATCH --mail-user douglas.scofield@ebc.uu.se
#SBATCH --mail-type=ALL

module load bioinfo-tools
module load BioPerl/1.6.922

TmpID=${SLURM_JOB_ID:-$$}

BlastInput=$1

Prefix=${BlastInput%.txt.gz}.psari

echo -n "Time started psari.pl: " 
date

DiscoveryOutput=$Prefix.prelim.gff3

./psari.pl $BlastInput > $DiscoveryOutput

echo -n "Time started psari_remove-HC-overlaps.sh: " 
date

ScreeningOutput=$Prefix.HC-screened.gff
ScreeningOverlapOutput=$Prefix.HC-overlap.gff  # those that do overlap HC genes

./psari_screen-HC-overlaps.sh $DiscoveryOutput

# The alternate screening script is perhaps more complete but may not run because
# it includes a grep that may gobble up too much memory
# ./psari_screen-HC-overlaps_alternate.sh $DiscoveryOutput  

echo -n "Time started psari_annotateBedWithVariants: " 
date

PsariOutput=$Prefix.full.gff3

./annotateBedWithVariants.sh $ScreeningOutput | gzip -c > $PsariOutput.gz

# convenience scripts produced by the annotate script
rm -f psari_join-annotation.pl psari_produce-annotation.awk  


echo -n "Time ended: " 
date
echo
echo "Output in $PsariOutput.gz" 

