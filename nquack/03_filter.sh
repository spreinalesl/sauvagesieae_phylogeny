#!/bin/bash

#SBATCH --job-name=filter
#SBATCH --partition=medium
#SBATCH --cpus-per-task=1
#SBATCH --mem=100MB
#SBATCH --mail-user=spreinalesl@alumni.usp.br
#SBATCH --mail-type=END,FAIL

###########
# NOTES
###########
## Filter low quality
## To remove regions of low quality based on the MAPQ score. MAPQ is equal to
## -10*log-base-10(Pr(mapping position is wrong)), rounded to the nearest integer
## In this case, use a very stringent approach and remove any sites with that have
##  a 10% chance or more of being mapped to the wrong location

###########
# COMMANDS
###########

## Paths
bampath=/home/sreinale/scratch/private/nquack/example/sauv/bamfiles

# outpath
mkdir -p /home/sreinale/scratch/private/nquack/example/sauv/filtered
outpath=/home/sreinale/scratch/private/nquack/example/sauv/filtered

## Define variable to list .bam files and iterate over each element using the SLURM_ARRAY_TASK_ID variable
#bam_list=($(ls * | grep bam$))
#SAMP_ID=(${bam_list[$SLURM_ARRAY_TASK_ID]})

## Activate conda environment
source activate nquack

## Remove low quality positions
cd $bampath
SAMP_LIST=$(ls * | grep bam)
for SAMP_ID in $SAMP_LIST; do
 samtools view -b -q 15 $SAMP_ID > $outpath/$SAMP_ID
 echo "END samtools for" $SAMP_ID
done
