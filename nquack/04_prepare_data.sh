#!/bin/bash

#SBATCH --job-name=nQuack
#SBATCH --partition=medium
#SBATCH --cpus-per-task=1
#SBATCH --mem=250MB
#SBATCH --mail-user=xxx@alumni.usp.br
#SBATCH --mail-type=END,FAIL

###########
# NOTES
###########
## Prepare data using samtools. The .bam file will convert to a tab-separated text file
## using the R script 04_prepare_data.R. Seems that this function should not be run into
## an array job because it generates a temp directory that causes conflict when different
## samples are runing in the same node.

###########
# COMMANDS
###########

## Paths
cd ~/projects/nquack

## Create a directory for the output
mkdir -p ~/projects/nquack/processed

## Activate conda environment
source activate nquack

## Run r script
Rscript --vanilla 04_prepare_data.R
