#!/bin/bash

#SBATCH --job-name=snaq
#SBATCH --partition=long
#SBATCH --cpus-per-task=21
#SBATCH --mem=57G
#SBATCH --mail-user=s.reinales@kew.org
#SBATCH --mail-type=END,FAIL

###########
# NOTES
###########
# This script runs SNaQ for mb trees with h0 (no hibridisation), h1 to h5, each one after the previous.
# H0 ran with 33G for 25h
# H1 ran with 41G for 15h
# H2 ran with 42G for 20h
# H3 ran with 41G for 32h
# H4 ran with 45G for
# H5 did not rin with 57G :(

###########
# COMMAND #
###########

## Out path
outpath=~/scratch/private/phylonetworks/snaq

## Set script and file names
outh=h5
jscript=network_estimation_$outh.jl

## Run Julia script
echo BEGIN!

/home/sreinale/apps/julia/bin/./julia $jscript > $outpath/network_estimation_$outh.outerr 2>&1

echo END!
