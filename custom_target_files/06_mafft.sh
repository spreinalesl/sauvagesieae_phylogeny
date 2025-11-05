#!/bin/sh

###########
# COMMAND #
###########
path=~/projects/custom_target/alignments

cd $path

for i in `ls *.fas`; do
  mafft --localpair --adjustdirectionaccurately --maxiterate 1000 $i > align_$i;
done
