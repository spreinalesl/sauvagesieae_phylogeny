#!/bin/bash

#SBATCH --job-name=nQuack09
#SBATCH --partition=short
#SBATCH --cpus-per-task=1
#SBATCH --mem=100MB
#SBATCH --mail-user=xxx@alumni.usp.br
#SBATCH --mail-type=END,FAIL

###########
# NOTES
###########
## Select the best model for each sample using BIC

###########
# COMMANDS
###########
cd ~/projects/nquack

## Activate conda environment
source activate nquack

## Run r script
Rscript --vanilla 09_summary_boot.R
