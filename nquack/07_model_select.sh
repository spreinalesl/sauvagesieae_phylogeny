#!/bin/bash

#SBATCH --job-name=nQuack07
#SBATCH --partition=short
#SBATCH --cpus-per-task=1
#SBATCH --mem=100MB
#SBATCH --mail-user=xxx@alumni.usp.br
#SBATCH --mail-type=END,FAIL

###########
# NOTES
###########
## Select the best model for each sample using BIC. Run this script with min 100MB
## to load the environment

###########
# COMMANDS
###########
cd ~/projects/nquack

## Create out directory
mkdir -p ~/projects/nquack/models
mkdir -p ~/projects/nquack/models/interpret

## Activate conda environment
source activate nquack

## Run r script
Rscript --vanilla 07_model_select.R
