#!/bin/bash

#SBATCH --job-name=burnin
#SBATCH --partition=medium
#SBATCH --cpus-per-task=1
#SBATCH --mem=50MB
#SBATCH --mail-user=s.reinales@kew.org
#SBATCH --mail-type=END,FAIL

###########
# NOTES
###########
## This script read .t files form MrBayes runs (4) containing 1000 trees each. The, It excludes 500 samples at the beginning of the file
## and subsample every 25 trees for a total of 20 samples per loci, resulting in 20*522=10440 samples. Selected trees are then converted
## in a newick format

###########
# COMMANDS
###########

## Input oath
fpath=~/scratch/private/phylonetworks/alignments/trimmed_alignm

## Out path
outpath=~/scratch/private/phylonetworks/mbresults

## List log files
cd $fpath
logfiles=$(ls | grep nexus$)

## Iterate nexus files and sample trees
for j in $logfiles; do
 echo $j
 /home/sreinale/scratch/apps/conda/bin/pxlog -t $j.run*.t -b 500 -n 25 | /home/sreinale/scratch/apps/conda/bin/pxt2new >> $outpath/lesstrees_allloci.newick
done

echo DONE!
