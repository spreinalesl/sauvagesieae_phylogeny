###########################################
## This script read trees and change 
## sample names
###########################################

## Packages
library(ape)

## Input paths
namesfile <- "sppidTOsppname.csv"
fpath <- paste("~/projects/rbgk/users_area/sreinales/results/sauvagesieae/hybpiper/fasta_files",
                "/", probename, "/", type, "/", sep = "")

## Read data and iterate files
ffiles  <- list.files(fpath, pattern = "\\.FNA")

# Iterate files
for(i in ffiles){

  ## Read tree
 

  ## Change sample names from a ID to spp names file
  smp_names <- read.csv(paste0(namespath, namesfile, collapse = ""), as.is = T)

  fnames <- unlist(lapply(X = names(fasta),
                                   FUN = function(x) strsplit(x, split = " ")[[1]][1]))
  for(n in 1:length(fnames)){
    names(fasta)[n] <- smp_names$abbrev[which(smp_names$seqID == fnames[n])]
  }
}
