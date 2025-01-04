#!/bin/bash

#SBATCH --job-name=mb
#SBATCH --partition=long
#SBATCH --array=0
#SBATCH --cpus-per-task=32
#SBATCH --mem=2G
#SBATCH --mail-user=s.reinales@kew.org
#SBATCH --mail-type=END,FAIL

###########
# NOTES
###########
## Run MrBayes for each locus
## Memory requirements info from https://wiki.eduuni.fi/display/cscfgi/MrBayes+runtime+environment
## The 39% of loci failed with 250MB (mainly ochn loci). Then, we re-run all loci with 750MB.

###########
# COMMAND #
###########

## Input fasta
fpath=~/scratch/private/phylonetworks/alignments/trimmed_alignm

## Go to input-file folder
cd $fpath

## Define variable to list all alignment and iterate over each element using the SLURM_ARRAY_TASK_ID variable
gene_list=($(ls * | grep mb$))
gene=(${gene_list[$SLURM_ARRAY_TASK_ID]})
#gnm=$(echo $gene | sed 's/_align\.fas\.nexus//')

## Run mb templates
echo Run mb for $gene

mpirun -n ${SLURM_CPUS_PER_TASK} -oversubscribe /home/sreinale/scratch/apps/addprograms/MrBayes/src/mb $gene > $gene.outerr 2>&1

echo mb END!
