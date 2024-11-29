#!/bin/bash

#SBATCH --job-name=cf
#SBATCH --partition=medium
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --mail-user=s.reinales@kew.org
#SBATCH --mail-type=END,FAIL

###########
# NOTES
###########
# This script calculates concordance factors for a combined sample of trees. We have 10440 trees.
# A sample of 100.000 trees for this project requires allocation of roughly 17Gb. Then, we used 3G

###########
# COMMAND #
###########

## Out path
outpath=~/scratch/private/phylonetworks/mbresults

## Activate conda environment
source activate juliaFORnetw

## Run Julia script
echo BEGIN!

/home/sreinale/scratch/apps/conda/envs/juliaFORnetw/julia/bin/./julia 04_CF_calculation.jl > $outpath/CF_calculation.outerr 2>&1

echo END!
