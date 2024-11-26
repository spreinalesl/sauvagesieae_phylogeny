#!/bin/bash

## Read log files
cd /home/spreinalesl/Documents/subsampling

## iterate fasta files
# Replace ID for Adenarake to merge a353 and ochn data
for file in `ls *.FNA`; do
 if grep -q -e ERR5034038 $file; then
  sed -i 's/ERR5034038/SRR10951608/' $file
  echo $file edited
 fi
done

echo First loop END!


# Replace ID for Ochna_membranacea to merge a353 and ochn data
for file in `ls *.FNA`; do
 if grep -q -e ERR5034076 $file; then
  sed -i 's/ERR5034076/SRR14767553/' $file
  echo $file edited
 fi
done

echo Second loop END!
