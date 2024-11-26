#################################
## Read fasta files and remove
## no selected samples for phylonetworks
#################################

library(ape)

## input path
ffpath <- "/home/spreinalesl/Documents/subsampling"

## List of selected taxa
# IDs for Adenarake muriculata and Ochna membranaceae were previously changed in order to merge a353 and ochn data
selspp <- read.csv("~/Dropbox/Sandra/Manuscritos/NGS_phyloOchn/analyses/phylonetworks/sampling.csv",
                   sep = ";")

## Iterate fasta files
warning_list1 <- c()
fflist <- list.files(pattern = "FNA", ffpath, full.names = T)

for(i in fflist){
  ff <- read.FASTA(i)
  fftips <- gsub(labels(ff), pattern = " .+", replacement = "")
  
  ## Extract positions
  pst <- which(fftips %in% selspp$ID)
  
  if(length(pst) > 0){
    ## Write new fasta
    ff <- ff[pst]
    
    write.dna(ff, format = "fasta", nbcol = 1, colsep = "", colw = 60, indent = 0, blocksep = 0,
              file = paste(ffpath, "/ffed", gsub(i, pattern = ffpath, replacement = ""), sep = ""))
    
    cat("Editing", i, "\n")
    
  } else{
    warning_list1 <- rbind(warning_list1, gsub(i, pattern = ffpath, replacement = ""))
  }
}

cat("No samples selected for gene", warning_list1, sep = "\n")