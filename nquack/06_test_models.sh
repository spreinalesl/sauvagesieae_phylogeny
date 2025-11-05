#!/bin/bash

#SBATCH --job-name=nQuack06
#SBATCH --partition=long
#SBATCH --array=0-79%5
#SBATCH --cpus-per-task=10
#SBATCH --mem=2G
#SBATCH --mail-user=xxx@alumni.usp.br
#SBATCH --mail-type=END,FAIL

###########
# NOTES
###########
## Test 18 model types for all three distributions (Normal, Beta and Beta-binomial),
## totalising 54 models
## First, set analysis options in the R script

###########
# COMMANDS
###########

## Create out directory
mkdir -p ~/projects/nquack/models

## Activate conda environment
source activate nquack

## Create task ID
cd ~/projects/nquack/processed
file_list=($(ls * | grep -e '_dd.csv'))
SAMP_ID=(${file_list[$SLURM_ARRAY_TASK_ID]})

echo $SAMP_ID

## Run r script
Rscript --vanilla ../06_test_models.R $SAMP_ID
