## Load packages
library(nQuack)

## Paths
inpath <- "~/projects/nquack/filtered/"
outpath <- "~/projects/nquack/processed/"

## Iterate files using an array job
samp_list <- list.files(pattern = "bam", inpath)
samp <- gsub(".bam", "", samp_list)
cat("Test", samp,  "\n")

## Iterate .bam files to convert them tables
for(i in 1:length(samp)){
 prepare_data(samp[i], inpath, outpath)
 cat("Processing", samp[i], "\n")
}