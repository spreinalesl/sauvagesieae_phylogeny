#!/bin/bash

#SBATCH --job-name=nQuack
#SBATCH --partition=medium
#SBATCH --cpus-per-task=1
#SBATCH --mem=100MB
#SBATCH --mail-user=spreinalesl@alumni.usp.br
#SBATCH --mail-type=END,FAIL

###########
# NOTES
###########
## Select the best model for each sample using BIC

###########
# COMMANDS
###########
cd /home/sreinale/scratch/private/nquack/example/sauv

## Activate conda environment
source activate nquack

## Run r script
Rscript --vanilla 09_summary_boot.R
