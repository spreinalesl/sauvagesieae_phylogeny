#!/bin/bash

#SBATCH --job-name=nQuack08
#SBATCH --partition=long
#SBATCH --array=0-59%10
#SBATCH --cpus-per-task=1
#SBATCH --mem=750MB
#SBATCH --mail-user=spreinalesl@alumni.usp.br
#SBATCH --mail-type=END,FAIL

###########
# NOTES
###########
## Bootstrap under the best model
## first, set the analysis options in the R script

###########
# COMMANDS
###########

## Create out directory
mkdir -p /home/sreinale/scratch/private/nquack/example/sauv/bootstrap

## Activate conda environment
source activate nquack

## Create task ID
cd /home/sreinale/scratch/private/nquack/example/sauv/processed
file_list=($(ls * | grep -e '_dd.csv'))
SAMP_ID=(${file_list[$SLURM_ARRAY_TASK_ID]})

## Run r script
Rscript --vanilla ../08_bootstrap.R $SAMP_ID
