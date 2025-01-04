#!/bin/bash

#SBATCH --job-name=CFs
#SBATCH --partition=medium
#SBATCH --cpus-per-task=1
#SBATCH --mem=700MB
#SBATCH --mail-user=s.reinales@kew.org
#SBATCH --mail-type=END,FAIL

###########
# NOTES
###########
# This script calculates concordance factors for 592 gene trees estimated with iqtree.
# 700MB were required for 592 gene trees. Time 3min

###########
# COMMAND #
###########

## Out path
mkdir -p ~/scratch/private/phylonetworks/snaq
outpath=~/scratch/private/phylonetworks/snaq

## Run Julia script
echo BEGIN!

/home/sreinale/apps/julia/bin/./julia 04_CF_calculation.jl > $outpath/CF_calculation.outerr 2>&1

echo END!
