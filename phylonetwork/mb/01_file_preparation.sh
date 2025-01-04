#!/bin/bash

#SBATCH --job-name=makeTemp
#SBATCH --partition=medium
#SBATCH --cpus-per-task=1
#SBATCH --mem=5MB
#SBATCH --mail-user=s.reinales@kew.org
#SBATCH --mail-type=END,FAIL

###########
# NOTES
###########
## Create the template for MrBayes

###########
# COMMANDS
###########

## Input fasta
fpath=~/scratch/private/phylonetworks/alignments/trimmed_alignm

## Path to template
tempath=~/scratch/private/phylonetworks/scripts

## process each nexus file to create the script and move them to the appropriate dir
cd $fpath

for i in `ls *.nexus`; do
    # create the mb script from the template file and just append the extension .mb
    touch $i.mb
    cat $tempath/01_template.mb > $i.mb
    # replace the placeholder ALIGNMENT with the alignment filename
    sed -i "s/ALIGNMENT/$i/g" $i.mb
    echo "creating mb template file for $i"
done

echo RUN END!
