#!/bin/bash

#SBATCH --job-name=qcont_sum
#SBATCH --partition=medium
#SBATCH --cpus-per-task=1
#SBATCH --mem=100MB
#SBATCH --mail-user=s.reinales@kew.org
#SBATCH --mail-type=END,FAIL

###########
# NOTES
###########
## Run R script 02_tree_qualitycontrol_sum.R to sumarize amd plot sdsf values across trees. Use the table generated
## by 02_tree_qualitycontrol_gen.sh as input.

###########
# COMMANDS
###########

## Script path
SPATH=~/scratch/private/phylonetworks/scripts

## Activate conda environment
source activate renv

## Run r script
cd $SPATH

Rscript --vanilla 02_tree_qualitycontrol_sum.R

echo END run!
