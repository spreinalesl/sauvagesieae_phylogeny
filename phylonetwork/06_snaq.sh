#!/bin/bash

#SBATCH --job-name=snaq
#SBATCH --partition=long
#SBATCH --cpus-per-task=21
#SBATCH --mem=37G
#SBATCH --mail-user=s.reinales@kew.org
#SBATCH --mail-type=END,FAIL

###########
# NOTES
###########
# This script runs SNaQ for mb trees with h0 (no hibridisation), h1, h2 and h3, each one after the previous.
# H0 ran with 33G for 24h
# H1 ran with 35G for 14h
# H2 ran with 37G for 53h
# H3 ran with 40G for 28h

###########
# COMMAND #
###########

## Out path
outpath=~/scratch/private/phylonetworks/snaq

## Set script and file names
outh=h1
jscript=network_estimation_$outh.jl

## Run Julia script
echo BEGIN!

/home/sreinale/apps/julia/bin/./julia $jscript > $outpath/network_estimation_$outh.outerr 2>&1

echo END!
