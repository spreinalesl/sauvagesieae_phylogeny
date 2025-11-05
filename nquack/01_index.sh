#!/bin/bash

#SBATCH --job-name=bwa
#SBATCH --partition=medium
#SBATCH --cpus-per-task=1
#SBATCH --mem=5MB
#SBATCH --mail-user=xxx@alumni.usp.br
#SBATCH --mail-type=END,FAIL

###########
# NOTES
###########
## Indexing the reference sequence

###########
# COMMANDS
###########
cd ~/projects/nquack/reference_seqs

## Activate conda environment
source activate nquack

bwa-mem2 index Cespedesia_spathulata_reference.fasta

echo END!
