#!/bin/bash

#SBATCH --job-name=trimm
#SBATCH --partition=medium
#SBATCH --array=0-19
#SBATCH --cpus-per-task=8
#SBATCH --mem=4G
#SBATCH --mail-user=spreinalesl@alumni.usp.br
#SBATCH --mail-type=END,FAIL

###########
# NOTES
###########
# Sliding Window Size: length of window in bp. SLIDINGWINDOW:4:20 to use a sliding window of size 4 that will remove bases if their phred score is below 20.
# By considering multiple bases, a single poor quality base will not cause the removal of high quality data later in the read.
# LEADING:20 would delete all bases below a quality threshold of 20 or that are N, beginning at the first base and continuing until the first base that is at least a quality of 20
# and is not N. Q20 means you will accept minimum 53 of quality score per base that corresponds to 0.01 of probability of error.
# MINLEN:30 to trim reads shorter than 30.
# According to google 1 node, 8 tasks or processors, 8gb for memory and 1 hour for walltime should be enough to run Trimmomatic (but haw many reads?).
# This analysis took 5h32min using 16 theads and just 2.98G to process 192 files with ca. 3M reads each.
# BUT using an array is extremely efficient. It is not necesary to use more than 4 or 8 threads.
# This code will perform the following:
# Remove adapters using the palindrome mode (designed for reading through a short fragment into adapter seq)
# Remove leading low quality or N bases (below quality 20)
# Remove trailing low quality or N bases (below quality 20)
# Scan the read with a 4-base wide sliding window, cutting when the average quality per base drops below 20
# Drop reads below the 30 bases long

###########
# COMMANDS
###########

## Define names, files and folders
dpath=/home/sreinale/projects/rbgk/users_area/sreinales/rawdata/addseqs_schnkit
rpath=/home/sreinale/scratch/private/nquack/example/sauv/reads
adapters=/home/sreinale/apps/conda/share/trimmomatic/adapters

## Define variable to read the sample list and iterate over each element using the SLURM_ARRAY_TASK_ID variable
sample_list=($(cat $rpath/sampling.txt))
sample=(${sample_list[$SLURM_ARRAY_TASK_ID]})

## Untar rawdata
cd $TMPDIR
tar -xzvf $dpath/allfiles.tar.gz ${sample}_1.fastq.gz
tar -xzvf $dpath/allfiles.tar.gz ${sample}_2.fastq.gz

## Run trimmomatic
java -jar /home/sreinale/apps/conda/share/trimmomatic/trimmomatic.jar PE \
     -threads ${SLURM_CPUS_PER_TASK} \
     -phred33 \
     ${sample}_1.fastq.gz \
     ${sample}_2.fastq.gz \
     ${sample}_1.trim.fastq.gz \
     ${sample}_1.unp.fastq.gz \
     ${sample}_2.trim.fastq.gz \
     ${sample}_2.unp.fastq.gz \
     ILLUMINACLIP:$adapters/TruSeq3-PE-2.fa:2:30:10:2:True \
     SLIDINGWINDOW:4:20 \
     LEADING:20 \
     TRAILING:20 \
     MINLEN:30 > trimmomatic.outerr 2>&1

echo END trimmomatic!

mv -v ${sample}_1.trim.fastq.gz $rpath
mv -v ${sample}_2.trim.fastq.gz $rpath
