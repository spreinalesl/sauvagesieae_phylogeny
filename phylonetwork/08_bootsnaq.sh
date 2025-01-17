#!/bin/bash

#SBATCH --job-name=bootsnaq
#SBATCH --partition=long
#SBATCH --array=1-100%5
#SBATCH --cpus-per-task=6
#SBATCH --mem=41G
#SBATCH --mail-user=s.reinales@kew.org
#SBATCH --mail-type=END,FAIL

###########
# NOTES
###########
## Run 100 bootstrap iterations for the best network, based on the collection of trees resulted of the Ultrafast BT from Iqtree,
## for each locus for a total of 592000 gene trees (100 per locus)
## Each iteration has 5 independent runs

###########
# COMMANDS
###########

## Out path
outpath=~/scratch/private/phylonetworks/bootsnaq

## Define variables
NRUNS=5
#SEED=($(echo {620..719})[$SLURM_ARRAY_TASK_ID]])
SEED=($SLURM_ARRAY_TASK_ID)

## Use array variable to run snaq in parallel
echo BEGIN SNaQ iteration $SEED

/home/sreinale/apps/julia/bin/./julia bootsnaq.jl $NRUNS $SEED > $outpath/bootsnaq$SLURM_ARRAY_TASK_ID.outerr 2>&1

echo SNaQ FINISHED!
