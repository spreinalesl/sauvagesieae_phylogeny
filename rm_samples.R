#################################################
## This script removes samples with less than 
## XX% of bp and genes recovered 
##################################################

## Packages
library(ape)

## Paths
fpath <- "~/projects/hybpiper/fasta_files/schn/exons"
tablespath <- "~/projects/name_list"
smpTOrm <- "samplesTOremove33_all.csv"
outfastadir <- paste(fpath, "smprm", sep = "/")
outnm <- "rm30_smpANDgns_schn_summ.csv"

## Set min number of sampling to remove a gene
# for a353 35 samples is the min samples number
# for schn 46 samples is the min sample number
minsmp <- 45
  
## Read tables
smpnm <- read.csv(paste(tablespath, smpTOrm, sep = "/"), sep = "\t")
smpnm <- unique(smpnm$sample)

## List of fasta files to edit
fastafl <- list.files(pattern = "FNA", fpath)

## Save summary
fsumm <- data.frame(gene = rep(NA, length(fastafl)), nsample_org = NA, nsample_ed = NA, smprm = NA, gnrm = NA)

## Iterate fasta files
for(i in 1:length(fastafl)){
  
  ## Gene name
  gnm <- unlist(strsplit(fastafl[i], split = "\\."))[1]
  
  ## Read tree
  fastaf <- read.FASTA(paste(fpath, fastafl[i], sep = "/"))
  
  ## Extract info
  fsumm$gene[i] <- gnm
  fsumm$nsample_org[i] <- length(fastaf)
  
  ## If there is no sequences in this file, jump to the next iteration
  if(length(fastaf) < 1) {
    fsumm$gnrm[i] <- "Yes"
    next
  }
  
  ## If original fasta file has less than XX samples, jump to the next iteration
  if(length(fastaf) < minsmp) {
    fsumm$gnrm[i] <- "Yes"
    next
  }
  
  ## Extract tip labels
  tips <- unlist(lapply(strsplit(labels(fastaf), split = " "), function(x) x[1]))
  
  ## Positions to remove
  rmpst <- which(tips %in% smpnm)
  
  ## Samples to remove
  fsumm$smprm[i] <- length(rmpst)
  
  ## If there is no sequences to remove write original 
  ## fasta file and jump to next iteration
  if(length(rmpst) == 0) {
    fsumm$nsample_ed[i] <- length(fastaf)
    write.FASTA(fastaf, file = paste(outfastadir, "/", gnm, ".FNA", sep = ""))
    next
  }
  
  ## Edit fasta file
  fastaf <- fastaf[-rmpst]
  
  ## Extract info
  fsumm$nsample_ed[i] <- length(fastaf)
  
  ## Write edited fasta file just for genes with more than XX% of samples
  if(length(fastaf) > minsmp){
    write.FASTA(fastaf, file = paste(outfastadir, "/", gnm, ".FNA", sep = ""))
  } else{ fsumm$gnrm[i] <- "Yes" }
  
  cat("Editing", gnm, "gene", "\n")
}

## Write summ table
write.csv(fsumm, paste(tablespath, "/", outnm, sep = ""))