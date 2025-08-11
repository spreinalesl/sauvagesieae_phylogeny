#!/bin/bash

#SBATCH --job-name=nQuack05
#SBATCH --partition=medium
#SBATCH --array=0-59%5
#SBATCH --cpus-per-task=1
#SBATCH --mem=100MB
#SBATCH --mail-user=spreinalesl@alumni.usp.br
#SBATCH --mail-type=END,FAIL

###########
# NOTES
###########
## Filter data by total and allele coverage, and remove noise using the
## R script 05_process_data.R

###########
# COMMANDS
###########

## Paths
inpath=/home/sreinale/scratch/private/nquack/example/sauv/processed

## Activate conda environment
source activate nquack

## Define variable to list all files and iterate over each element using the SLURM_ARRAY_TASK_ID variable
cd $inpath
file_list=($(ls | grep txt))
SAMP_ID=(${file_list[$SLURM_ARRAY_TASK_ID]})

## Run r script
Rscript --vanilla ../05_process_data.R $SAMP_ID
