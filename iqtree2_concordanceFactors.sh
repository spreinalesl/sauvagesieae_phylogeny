#!/bin/bash

## This script calculates the gene concordance factor (gCF) for each branch of the species tree as the fraction 
## of decisive gene trees concordant with this branch and the site concordance factor (sCF) as the fraction 
## of decisive alignment sites supporting that branch 

## Set names
DATASET=inform_prlgs
GTREES=sauv_putative_orth_exons_gtrees_nm_ed.tre
SPTREE=sauv_putative_orth_exons_BS33_sptree_annot.tre
ALPATH="./names_ed"

## Go to path
cd /home/spreinalesl/Dropbox/Sandra/Manuscritos/NGS_phyloOchn/analyses/results/trees/concordance_factors/$DATASET

## Calculate gcf base on the species tree from wAstral and the ggtrees infered with iqtree (all in the same file)
iqtree2 -nt 8 -t $SPTREE --gcf $GTREES --prefix gCF

echo END gcf!

## Calculated sCF base on the species tree from wAstral and the alignments of each locus. Scores will be appended to the existing label separated by a /
## in the output file from gcf
#iqtree2 -nt 8 -t $SPTREE -s $ALPATH -m GTR+I+G --scfl 100 --prefix sCF

#echo END scfl!
