#!/bin/bash

#SBATCH --job-name=sCF
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=3G
#SBATCH --mail-user=s.reinales@kew.org
#SBATCH --mail-type=END,FAIL

###########
# NOTES
###########
## This script calculates the gene concordance factor (gCF) for each branch of the species tree as the fraction 
## of decisive gene trees concordant with this branch and the site concordance factor (sCF) as the fraction 
## of decisive alignment sites supporting that branch
## For interpretation and nice plot of results see http://www.robertlanfear.com/blog/files/concordance_factors.html

###########
# COMMANDS
###########
## Go to path
cd ~/scratch/private/str_mssdata_prlgs

## Calculate gcf base on the species tree from wAstral and the ggtrees infered with iqtree (all in the same file)
#/home/sreinale/scratch/apps/addprograms/iqtree2/bin/iqtree2 -nt ${SLURM_CPUS_PER_TASK} -t sauv_smp30_prlgbothgns_outliers_rm_exons_BS33_sptree.tre --gcf sauv_smp30_prlgbothgns_outliers_rm_exons_gtrees_nm_ed.tre --prefix gCF

#echo END gcf!

## Calculated sCF base on the species tree from wAstral and the alignments of each locus. Scores will be appended to the existing label separated by a /
## in the output file from gcf
/home/sreinale/scratch/apps/addprograms/iqtree2/bin/iqtree2 -nt ${SLURM_CPUS_PER_TASK} -t sauv_smp30_prlgbothgns_outliers_rm_exons_BS33_sptree.tre -s ./names_ed -m GTR+I+G --scfl 100 --prefix sCF

echo END scfl!
