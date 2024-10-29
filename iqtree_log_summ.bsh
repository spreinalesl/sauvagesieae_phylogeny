#!/bin/bash

## Read log files
cd /home/spreinalesl/Documents/iqtree_summ

## tar files (one per dataset)
tarfiles=$(ls | grep tar)
echo $tarfiles

for j in $tarfiles; do
 
 ## Dataset name
 DATASET=$(echo $j | sed 's/_iqtree.tar.gz//')
 
 ## Extract log files
 tar -xzf $j
 
 ## Log files (one per gene)
 logfiles=$(ls | grep fas)

 ## Iterate logfiles and extract information

 for i in $logfiles; do
  ## Extract info
  FNM=$(grep "Input file name" $i | cut -d " " -f 4) ## File name
  GNM=$(echo $FNM | cut -d "_" -f 2) ## Gene name
  NSAMP=$(grep "Input data" $i | cut -d " " -f 3) ## Number of samples in this tree
  NSITES=$(grep "Input data" $i | cut -d " " -f 6) ## Number of sites
  PINFS=$(grep "Number of parsimony informative sites" $i | cut -d " " -f 6) ## Number of prs info sites
  IBLENG=$(grep "Sum of internal branch lengths" $i | cut -d " " -f 6) ## Summ of internal branch lengths
  TBLENG=$(grep "Total tree length" $i | cut -d " " -f 8) ## Summ of all branch lengths
 
  ## Write info in a file
  echo "$DATASET;$FNM;$GNM;$NSAMP;$NSITES;$PINFS;$IBLENG;$TBLENG" >> ./iqtree_summ.txt
 done
 
 rm *.fas.iqtree
 echo $j DONE!
 
done
