########################
## Extract support info
## from gene trees
########################

library(ape)

setwd("")

## Set path
tpath <- "iqtree_summ"
pattern1 <- "sauv_"
pattern2 <- "_gtrees.tre"

## iterate datasets (each file contains several gene trees)
tfiles <- list.files(pattern = "tre$", tpath)

gtrees_summ <- list()
allnames <- list()
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
  
  tipnames <- list()
  for(j in 1:ngs){
    
    ## Extract tip names
    tipnames[[j]] <- trs[[j]]$tip.label
    
    ## Save data
    dfsum[[j]] <- data.frame(ntips = length(trs[[j]]$tip.label),
                             tbrachn_sum = sum(trs[[j]]$edge.length, na.rm = T),
                             tbrachn_mean = mean(trs[[j]]$edge.length, na.rm = T),
                             tbrachn_median = median(trs[[j]]$edge.length, na.rm = T),
                             tbrachn_var = var(trs[[j]]$edge.length, na.rm = T),
                             btsupp = mean(as.numeric(trs[[j]]$node.label), na.rm = T))
    
    cat("Read tree", tfiles[i], j, "-", ngs, "\n")
  }
  
  ## Merge tip names
  allnames[[i]] <- unique(do.call(c, tipnames))
  
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

## Extract tip names
allnames <- unique(do.call(c, allnames))
  
## Write tables
write.csv(gtrees_summ, paste(tpath, "gtrees_summ.csv", sep = "/"), row.names = F)
writeLines(allnames, paste(tpath, "tipnames.csv", sep = "/"), sep = "\n")

system("cat iqtree_summ/*.csv > merge_iqtree_gtrees_sum.csv")

################################
## Get summary for table1
################################
dat <- read.csv(paste(tpath, "merge_iqtree_gtrees_sum.csv", sep = "/"))

## Prepare data
dat$probe[which(grepl("g[0-9]", dat$flname))] <- "schn"
dat$probe[which(is.na(dat$probe))] <- "a353"

dat$dataset <- as.factor(dat$dataset)
dat$probe <- as.factor(dat$probe)

## Extract info
## nloci
nloci <- aggregate(dat$ngenes ~ dat$dataset + dat$probe, data = dat, FUN = length)

## ntips
ntips <- aggregate(dat$ntips ~ dat$dataset, data = dat, FUN = median)
ntips$`dat$ntips` <- round(ntips$`dat$ntips`, digits = 1)

ntipsp <- aggregate(dat$ntips ~ dat$dataset + dat$probe, data = dat, FUN = median)

## nsites
nsitest <- aggregate(dat$nsites ~ dat$dataset, data = dat, FUN = sum)
nsitestp <- aggregate(dat$nsites ~ dat$dataset + dat$probe, data = dat, FUN = sum)

nsitespl <- aggregate(dat$nsites ~ dat$dataset, data = dat, FUN = median)
nsitesp <- aggregate(dat$nsites ~ dat$dataset + dat$probe, data = dat, FUN = median)

## Pinformative sites
npfsitest <- aggregate(dat$npars_inf_sites ~ dat$dataset, data = dat, FUN = sum)
npfsitestp <- aggregate(dat$npars_inf_sites ~ dat$dataset + dat$probe, data = dat, FUN = sum)