#!/bin/bash

#SBATCH --job-name=CFs
#SBATCH --partition=long
#SBATCH --cpus-per-task=4
#SBATCH --mem=30G
#SBATCH --mail-user=xxx@kew.org
#SBATCH --mail-type=END,FAIL

###########
# NOTES
###########

## This script calculates the gene concordance factor (gCF) for each branch based on a cnocatenated matrix,
## as the fraction of decisive gene trees concordant with this branch, and the site concordance factor (sCF) as the fraction
## of decisive alignment sites supporting that branch
## concatenate tree required 12 threads and 26G to run
## gene trees required 4 threads and 16G
## sCF required 20G

###########
# COMMANDS
###########

# set dir
DATASET=rmout_informprlgs

## Create out folder
mkdir -p ~/projects/concordance_factors
mkdir -p ~/projects/concordance_factors/$DATASET 

ALN_DIR=~/projects/concordance_factors/$DATASET

cd $ALN_DIR

# infer a concatenation-based species tree with 1000 ultrafast bootstrap and
# an edge-linked partition model in which all partitions share the same set of branch lengths but
# allowing each partition to have its own evolution rate in order to see how each partition inform the topology
/home/sreinale/apps/iqtree3/bin/iqtree3 -p $ALN_DIR --prefix concat -B 1000 -T AUTO --threads-max ${SLURM_CPUS_PER_TASK}

echo END species tree!

# infer locus trees
/home/sreinale/apps/iqtree3/bin/iqtree3 -S $ALN_DIR --prefix loci -T AUTO --threads-max ${SLURM_CPUS_PER_TASK}

echo  END locus trees!

# compute gene concordance factors
/home/sreinale/apps/iqtree3/bin/iqtree3 -t concat.treefile --gcf loci.treefile --prefix gCF

echo END gCF!

# compute site concordance factor using likelihood with v2.2.2
/home/sreinale/apps/iqtree3/bin/iqtree3 -t concat.treefile -p $ALN_DIR --scfl 100 --prefix sCF

echo END sCF!
