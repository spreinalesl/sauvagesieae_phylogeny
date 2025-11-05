#!/bin/bash

#SBATCH --job-name=bwa
#SBATCH --partition=medium
#SBATCH --array=0-79%10
#SBATCH --cpus-per-task=10
#SBATCH --mem=12G
#SBATCH --mail-user=xxx@alumni.usp.br
#SBATCH --mail-type=END,FAIL

###########
# NOTES
###########
## Mapping reads to reference

###########
# COMMANDS
###########
inpath=~/projects/nquack
targetpath=~/projects/nquack/reference_seqs
readspath=~/projects/nquack/reads
outpath=~/projects/nquack/bamfiles

## Define variable to list all samples and iterate over each element using the SLURM_ARRAY_TASK_ID variable
sampling=($(cat $inpath/sampling.txt))
samp=(${sampling[$SLURM_ARRAY_TASK_ID]})

echo $samp

## Untar reads files in the temp local disk
cp -v $readspath/${samp}_1.trim.fastq.gz $TMPDIR
cp -v $readspath/${samp}_2.trim.fastq.gz $TMPDIR

## Go to temp disk
cd $TMPDIR
gunzip -v ${samp}_1.trim.fastq.gz
gunzip -v ${samp}_2.trim.fastq.gz

## Activate conda environment
source activate nquack

## Mapping
bwa-mem2 mem -t 10 -M $targetpath/*.fasta ${samp}_1.trim.fastq ${samp}_2.trim.fastq > ${samp}.sam

echo "END mapping for" $samp

## SAM to BAM
samtools view -S -b ${samp}.sam > ${samp}.bam

echo "END bam for" $samp

## SORT
samtools sort ${samp}.bam -o $outpath/${samp}.bam

echo "END sort bam for" $samp
