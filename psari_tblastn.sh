#! /bin/bash -l

echo -n "Time started: " 
date

set -e
set -x

#SBATCH -A b2010042
#SBATCH -o tblastn.out
#SBATCH -e tblastn.err
#SBATCH -J tblastn.JOB
#SBATCH -p node -n 8
#SBATCH -t 3-00:00:00
#SBATCH --mail-user douglas.scofield@plantphys.umu.se
#SBATCH --mail-type=ALL

TmpID=${SLURM_JOB_ID:-$$}

module load bioinfo-tools
module load blast/2.2.24+

Query=$1
QueryName=$( basename $Query )
if [ "$Query" = "" ] ; then
	echo "Usage: $0 file.fa"
	exit 1
fi

# Target=/proj/b2010042/fosmidpools/data/MAIN_FP_FOLDERS/SCAFFOLDS/k51_650_2500_renamed/FosmidPools.all450.fa
# Target=/proj/b2010042/assembly/ASSEMBLYLOCK_2012_JULY/indices/blast/fosmidpools/picea_abies.fp.july2012.fa
Target=/proj/b2010042/nobackup/douglas/pseudogenes/picea_abies.rna-scaffolded.nov2012.fa
TargetName=$( basename $Target )

Output=pseudo.aa.${QueryName}-to-${TargetName}.$TmpID.tblastn.verbose.txt

tblastn -query $Query -db $Target -evalue 0.00001 -max_intron_length 100000 -window_size 1000 -num_threads 8 | gzip -c > $Output.gz

echo
echo "Output in $Output.gz" 
echo -n "Time finished: " 
date

