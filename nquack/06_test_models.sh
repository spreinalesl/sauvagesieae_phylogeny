#!/bin/bash

#SBATCH --job-name=nQuack
#SBATCH --partition=medium
#SBATCH --array=0-59%10
#SBATCH --cpus-per-task=10
#SBATCH --mem=200MB
#SBATCH --mail-user=spreinalesl@alumni.usp.br
#SBATCH --mail-type=END,FAIL

###########
# NOTES
###########
## Test 18 model types for all three distributions (Normal, Beta and Beta-binomial),
## totalising 54 models

###########
# COMMANDS
###########

## Create out directory
mkdir -p /home/sreinale/scratch/private/nquack/example/sauv/models

## Activate conda environment
source activate nquack

## Create task ID
cd /home/sreinale/scratch/private/nquack/example/sauv/processed
file_list=($(ls * | grep -e '_dd.csv'))
SAMP_ID=(${file_list[$SLURM_ARRAY_TASK_ID]})

echo $SAMP_ID

## Run r script
Rscript --vanilla ../06_test_models.R $SAMP_ID
