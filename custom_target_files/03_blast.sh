#!/bin/bash

###########
# NOTES
###########
# Commands based on this documentation https://manpages.ubuntu.com/manpages/xenial/en/man1/blastall_old.1.html
# bl2seq  performs  a  comparison  between  two  sequences using either the blastn or blastp algorithm.
# Both sequences must be either nucleotides or proteins.
# "blastn" compares a nucleotide query sequence against a nucleotide sequence database.
# -i query sequence(s) > sequence to compare
# -j subject sequence(s) > sequence to use as database or ref seq
# -T Produce HTML output
# -D 1 Output format: tabular
# -outfmt 10
# by default bl2seq output contains just hits with e-value < 10. For blastn -evalue 0.99 saves just hits with e-value < 0.99

###########
# COMMANDS
###########

## Define paths
ipath=~/projects/custom_target/edited_fasta/merged/*.fas
jpath=~/projects/custom_target/edited_fasta/CKDK-SOAPdenovo-Trans-assembly_ed.fas
rpath=~/projects/custom_target/blast/blastn_CKDK-ALL.csv

## Run blast using bl2seq
#bl2seq -i $ipath -j $jpath -p "blastn" -D 1 -o blastn_CKDK-All.txt

## Using blastn
blastn -query $ipath -subject $jpath -evalue 0.99 -outfmt 10 -out $rpath
