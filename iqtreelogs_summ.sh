#!/bin/bash

## Read log files
cd ~/projects/iqtree_summ

## tar files (one per dataset)
tarfiles=$(ls | grep tar)
echo $tarfiles

for j in $tarfiles; do
 
 ## Dataset name
 DATASET=$(echo $j | sed 's/_logfiles.tar.gz//')

 ## Extract log files
 tar -xzf $j
 
 ## Log files (one per gene)
 logfiles=$(ls | grep fas.log)

 ## Iterate logfiles and extract information

 for i in $logfiles; do
  ## Extract info
  FNM=$(grep -m 1 "Reading alignment file" $i | cut -d " " -f 4) ## File name
  GNM=$(echo $FNM | cut -d "_" -f 2) ## Gene name
  NSAMP=$(grep -m 1 "Alignment has" $i | cut -d " " -f 3) ## Number of samples in this tree
  NSITES=$(grep -m 1 "Alignment has" $i | cut -d " " -f 6) ## Number of sites
  PINFS=$(grep -m 1 "parsimony-informative" $i | cut -d " " -f 1) ## Number of prs info sites
   
  ## Write info in a file
  echo "$DATASET;$FNM;$GNM;$NSAMP;$NSITES;$PINFS" >> ./sum_iqtreelogs.txt
 done
 
 rm *.fas.log
 echo $j DONE!
 
done
