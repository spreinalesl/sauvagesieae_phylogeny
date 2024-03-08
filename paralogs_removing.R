###########################################
## This script check gene trees for genes
## containing paralogous sequences as tagged
## by Hybpiper, in case that paralogs are
## in different clades, they are removed from
## the original fasta matrix
###########################################

## Packages
library(ape)

## Input paths
apgtpth <- "/home/spreinalesl/Dropbox/Sandra/Manuscritos/NGS_phyloOchn/analyses/results/hybpiper/paralogs/paralogs_all/a353/trees/"
apgtl <- "/home/spreinalesl/Dropbox/Sandra/Manuscritos/NGS_phyloOchn/analyses/results/hybpiper/paralogs/paralogs_all/a353/paralogs_above_threshold_report.txt"
apnwf <- "/home/spreinalesl/Dropbox/Sandra/Manuscritos/NGS_phyloOchn/analyses/results/hybpiper/paralogs/paralogs_all/a353/paralogs_all/"

# output files
outmpath <- "/home/spreinalesl/Dropbox/Sandra/Manuscritos/NGS_phyloOchn/analyses/results/hybpiper/paralogs/paralogs_all/a353/paralogs_all/paralogs_rm/"

## Read paralogs_retrieve output and extract a list of loci with
## paralogous sequences
plist <- readLines(apgtl)
p <- grep(pattern = "genes with paralogs", x = plist)
loci_w_prlgs <- plist[(p+1):length(plist)] 

## Iterate tree files containing all paralogous sequences
ffiles  <- list.files(apgtpth, pattern = "\\.tre")
warning_list <- c()
mono_prlgs <- c()
non_mono_prlgs <- c()

for(i in ffiles){

  ## Extract locus name
  lnm <- gsub(pattern = "_.+", replacement = "", x = i)
  
  if(!(lnm %in% loci_w_prlgs)){
    warning_list <- rbind(warning_list, lnm)
  } else{
    
    ## Read tree
    t <- read.tree(paste(apgtpth, i, sep = ""))
    
    ## Extract list of paralogous sequences
    ## But first, test that just paralog names have point (.)
    pnm <- grep(pattern = "\\.", x = t$tip.label, value = T) 
    tst <- unique(gsub(pattern = "^.+\\.", replacement = "", x = pnm))
    allwdv <- c("main", 0:9)

    if(!unique(tst %in% allwdv)){
      stop("Tip labels should not contain points (.)")
    } else{
      locin <- unique(gsub(pattern = "\\..+", replacement = "", x = pnm))
      
      ## Iterate each sample with paralogous sequences, check if all sequences form a clade
      for(j in locin){
        ptips <- grep(pattern = j, x = pnm, value = T) 
        mono_test <- is.monophyletic(t, tips = ptips)
        
        if(!mono_test){
          ## if paralogs do not form a clade, remove them for the original fasta file
          
          #################################################################
          
          ## Removed sequences
          non_mono_prlgs <- c()
            
        } else{
          ## Write in a log file
          mono_prlgs <- rbind(mono_prlgs, ptips)
          
        }
        
        ## Test if main sequences and the other sequences recnstruct the same relations
        
      }
    }
  }
}

## Print warnings
cat("Those loci were not tagged as containing paralogs by Hybpiper" , warning_list,
    "", non_mono_prlgs,
    "Those paralogs form a clade, then, they are keept in the original fasta file", mono_prlgs, sep = "\n")