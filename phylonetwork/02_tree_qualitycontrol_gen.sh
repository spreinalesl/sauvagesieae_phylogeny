#!/bin/bash

#SBATCH --job-name=qcont
#SBATCH --partition=short
#SBATCH --cpus-per-task=1
#SBATCH --mem=5MB
#SBATCH --mail-user=s.reinales@kew.org
#SBATCH --mail-type=END,FAIL

###########
# NOTES
###########
## This script read the log files form MrBayes runs to check convergence

###########
# COMMANDS
###########

## Input oath
fpath=~/scratch/private/phylonetworks/alignments/trimmed_alignm

## Out path
outpath=~/scratch/private/phylonetworks/mbresults

## List log files
cd $fpath
logfiles=$(ls | grep outerr)

## Iterate log files
for i in $logfiles; do
 echo "Reading $i"
 #gnm=$(echo $i | sed 's/_align\.fas\.nexus\.log//')
 sdsf=$(grep "Average standard deviation of split frequencies =" $i | cut -d "=" -f 2)
 echo "$i $sdsf" >> $outpath/topostability.txt
done

echo DONE!
