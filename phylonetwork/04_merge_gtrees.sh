#!/bin/bash

#SBATCH --job-name=merge
#SBATCH --partition=medium
#SBATCH --cpus-per-task=1
#SBATCH --mem=5MB
#SBATCH --mail-user=s.reinales@kew.org
#SBATCH --mail-type=END,FAIL

###########
# NOTES
###########


###########
# COMMAND #
###########

## Input path
cd ~/scratch/private/phylonetworks/iqtree

## Merge all ML trees into one file
cat *.treefile >> lesstrees_allloci.newick

echo"Number of gene trees:"
ls * | grep treefile | wc -l
