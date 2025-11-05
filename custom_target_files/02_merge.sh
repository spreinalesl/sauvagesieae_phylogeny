#!/bin/bash

###########
# NOTES
###########
## Merge files containing assemblies in order to increase the number of Angiosperms353 genes to be blasted
# against CKDK scaffolds 

###########
# COMMANDS
###########

## Define paths
dpath=~/projects/custom_target/edited_fasta

## List of names
cd $dpath
ls *-SOAPdenovo-Trans-translated_ed.fas > names.txt
sed -i 's/-SOAPdenovo-Trans-translated_ed.fas//' names.txt
paste -sd "-" names.txt > names_ed.txt
names=$(cat names_ed.txt)

## Merge files
cat *-translated_ed.fas > ./merged/$names.fas

## Remove intermediate file
rm names*
