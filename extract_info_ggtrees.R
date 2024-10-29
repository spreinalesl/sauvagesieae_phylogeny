########################
## Extract support info
## from gene trees
########################

library(ape)

## Set path
tpath <- "~/Documents/iqtree_summ/gtrees"
pattern1 <- "sauv_"
pattern2 <- "_gtrees.tre"

## iterate datasets (each file contains several gene trees)
tfiles <- list.files(pattern = "tre", tpath)

gtrees_summ <- list()
for(i in 1:length(tfiles)){
  
  ## Dataset name
  dataset <- gsub(pattern = pattern1, tfiles[i], replacement = "")
  dataset <- gsub(pattern = pattern2, dataset, replacement = "")
  
  ## Read trees
  trs <- read.tree(paste(tpath, tfiles[i], sep = "/"))
  
  ## Number of genes in this dataset
  ngs <- length(trs)
  
  ## Iterate gene trees 
  dfsum <- list()
  
  for(j in 1:ngs){
    
    ## Save data
    dfsum[[j]] <- data.frame(ntips = length(trs[[j]]$tip.label),
                           tbrachn_sum = sum(trs[[j]]$edge.length),
                           tbrachn_mean = mean(trs[[j]]$edge.length),
                           tbrachn_median = median(trs[[j]]$edge.length),
                           tbrachn_var = var(trs[[j]]$edge.length, na.rm = T),
                           btsupp = mean(as.numeric(trs[[j]]$node.label), na.rm = T))
    
    cat("Read tree", tfiles[i], j, "-", ngs, "\n")
  }
  
  ## Collapse info for all trees into one data frame
  dfsum <- do.call(rbind, dfsum)
  
  ## Add columns
  dfsum$dataset <- dataset
  dfsum$ngenes <- ngs
  
  ## Save data
  gtrees_summ[[i]] <- dfsum
}

## Merge data for all datasets
gtrees_summ <- do.call(rbind, gtrees_summ)

## Write table
write.csv(gtrees_summ, "~/Documents/iqtree_summ/gtrees_summ.csv", row.names = F)
