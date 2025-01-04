#!/bin/bash

#SBATCH --job-name=fasTOnex
#SBATCH --partition=medium
#SBATCH --cpus-per-task=1
#SBATCH --mem=5MB
#SBATCH --mail-user=s.reinales@kew.org
#SBATCH --mail-type=END,FAIL

###########
# NOTES
###########
## Run pxs2nex function of phyx package to convert alignments from fasta to nexus

###########
# COMMANDS
###########

## Input fasta
fpath=~/scratch/private/phylonetworks/alignments/trimmed_alignm

## Iterate fasta files
cd $fpath

for i in `ls *.fas`; do
 /home/sreinale/scratch/apps/conda/bin/pxs2nex -s $i -o $i.nexus
 echo "Fasta $i converted to nexus"
done

echo DONE!
