#!/bin/bash

#SBATCH --job-name=prepfiles
#SBATCH --partition=short
#SBATCH --cpus-per-task=1
#SBATCH --mem=3MB
#SBATCH --mail-user=xxx@kew.org
#SBATCH --mail-type=END,FAIL

###########
# NOTES
###########
## Merge fboot trees for all gene trees into a single file

###########
# COMMANDS
###########

## Input path
cd ~/projects/phylonetworks/iqtree

## Outpath
mkdir -p ~/projects/phylonetworks/bootsnaq
outpath=~/projects/phylonetworks/bootsnaq

## Untar tree files
tar -xzvf bootrees.tar.gz

echo "Number of tree files:"
ls * | grep ufboot | wc -l

## Iterate files and subselect boot trees
logfiles=$(ls | grep ufboot)

for j in $logfiles; do
 echo $j
 /home/sreinale/apps/conda/bin/pxlog -t $j -b 0 -n 10 >> $outpath/fboottrees_allloci.newick
done

echo DONE!

## Merge all bootstrap trees into one file
#cat *.ufboot >> $outpath/fboottrees_allloci.newick

## Remove bootrees
rm -f -v *.ufboot

## Check previous steps
echo TEST!
head $outpath/fboottrees_allloci.newick
echo DONE!
