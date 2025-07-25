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

## Create out directory
mkdir -p /home/sreinale/scratch/private/nquack/example/sauv/models
mkdir -p /home/sreinale/scratch/private/nquack/example/sauv/models/interpret

## Activate conda environment
source activate nquack

## Run r script
Rscript --vanilla 07_model_select.R
