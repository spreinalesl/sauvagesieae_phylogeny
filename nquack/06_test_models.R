## Load packages
library(nQuack)

## Paths
inpath <- "~/projects/nquack/processed/"
outpath <- "~/projects/nquack/models/"

## Number of cores
ncores <- 10

## Text pattern of .csv files
tpattern <- "_dd.csv"

## Iterate files using an array job
sample <- commandArgs(trailingOnly = T)
cat("Begin models for", sample, "\n")

## Test 18 model types for all three distributions (Normal, Beta and Beta-binomial),
## totalising 54 models
temp <- as.matrix(read.csv(paste0(inpath, sample)))
  
## Test models
sample <- gsub(pattern = tpattern, replacement = "", sample)

out1 <- quackNormal(xm = temp, samplename = sample, cores = ncores, parallel = T)
cat("End normal", "\n")
  
out2 <- quackBeta(xm = temp, samplename = sample, cores = ncores, parallel = T)
cat("End beta", "\n")

out3 <- quackBetaBinom(xm = temp, samplename = sample, cores = ncores, parallel = T)
cat("End betabinom", "\n")

## Write output
allout <- rbind(out1, out2, out3)
write.csv(allout, file = paste0(outpath, sample, ".csv"),
          row.names = FALSE)
cat("Writing table", sample, "\n")