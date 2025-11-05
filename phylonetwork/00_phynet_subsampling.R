#################################
## Read fasta files and remove
## no selected samples for phylonetworks
#################################

library(ape)

setwd("")

## input path
ffpath <- "subsampling"

## List of selected taxa
selspp <- read.csv("supplementary_tables/TableS1_sampling.csv",
                   sep = ",")

unique(selspp$spp_name)

## Iterate fasta files
warning_list1 <- c()
fflist <- list.files(pattern = "FNA", ffpath, full.names = T)
length(fflist)

for(i in fflist){
  ff <- read.FASTA(i)
  fftips <- gsub(labels(ff), pattern = " .+", replacement = "")
  
  ## Extract positions
  pst <- which(fftips %in% selspp$ID)
  
  if(length(pst) > 0){
    ## Write new fasta
    ff <- ff[pst]
    
    ## Change tip names
    for(j in fftips){
      names(ff)[which(grepl(pattern = j, labels(ff)))] <- selspp$spp_name[which(selspp$ID == j)]  
    }
    
    ## Write new fasta file
    write.dna(ff, format = "fasta", nbcol = 1, colsep = "", colw = 60, indent = 0, blocksep = 0,
              file = paste(ffpath, "/ffed", gsub(i, pattern = ffpath, replacement = ""), sep = ""))
    
    cat("Editing", i, "\n")
    
  } else{
    warning_list1 <- rbind(warning_list1, gsub(i, pattern = ffpath, replacement = ""))
  }
}

cat("No samples selected for gene", warning_list1, sep = "\n")

########################
## Confirmation test
########################
fflist <- list.files(pattern = "FNA", "subsampling/ffed", full.names = T)

tipnames <- vector()
fftable <- data.frame(gene = NA, tipnum = NA)

for(i in fflist){
  ff <- read.FASTA(i)
  tipnames <- c(tipnames, names(ff))
  
  fftable <- rbind(fftable, data.frame(gene = i, tipnum = length(names(ff))))
}

tipnames <- unique(tipnames)
length(tipnames)

table(sort(tipnames) == sort(unique(selspp$spp_name)))

write.csv(fftable, "subsampling/ffed/fftable.csv", row.names = F)