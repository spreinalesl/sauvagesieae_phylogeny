library(phytools)

## List of trees to compare
tpath <- "/home/spreinalesl/Documents"
tfiles <- list.files(pattern = "\\.tre", tpath, full.names = T)
tpattern <- "/rooted_sauv_"
tpattern1 <- "_BS33_sptree_annot.tre"

## Iterate tree files
alltrees <- list()
tnames <- vector()

for(i in 1:length(tfiles)){
  
  ## Tree name
  tname <- gsub(pattern = tpath, replacement = "", tfiles[i])
  tname <- gsub(pattern = tpattern, replacement = "", tname)
  tname <- gsub(pattern = tpattern1, replacement = "", tname)
  tnames[i] <- tname
  
  ## Read tree
  alltrees[[i]] <- read.tree(tfiles[i])
  cat(is.rooted(alltrees[[i]]), "\n")
}

## Plot trees using phytools
# Define trees
tr1 <- 2
tr2 <- 7

## Support pattern
tsupp <- "q1="

# Generate cophylo
compt <- cophylo(alltrees[[tr1]], alltrees[[tr2]])

## Edit node labels for both trees in the cophylo (Label edition is only needed for Astral trees which include complex annotations)
# pp1 values are in the first position of the vector containing all node label information
# qq1 are in the 7º position

for(i in 1:2){
  
  ## Extract pp1 or qq1 values for each node
  tlabels <- strsplit(compt$trees[[i]]$node.label, split = "-")
  
  ## Position for target support
  pst <- grep(pattern = tsupp, tlabels[[2]])
  
  ## Extract pst 'position'
  tlabels <- unlist(lapply(tlabels[c(2:length(tlabels))], FUN = function(x) x[pst]))
  tlabels <- gsub(tlabels, pattern = tsupp, replacement = "")
  
  ## Convert support value in numeric and replace node labels
  compt$trees[[i]]$node.label <- round(as.numeric(tlabels), digits = 2)
}

## Plot cophylo
pdf(paste("~/Documents/comparing_", tnames[tr1], "_", tnames[tr2], ".pdf", sep = ""), width = 15, height = 24)

plot(compt, fsize = 0.7)

nodelabels.cophylo(compt$trees[[1]]$node.label[2:Nnode(compt$trees[[1]])],
                   2:Nnode(compt$trees[[1]]) + Ntip(compt$trees[[1]]), frame = "none",
                   cex = 0.7, which = "left")

nodelabels.cophylo(compt$trees[[2]]$node.label[2:Nnode(compt$trees[[2]])],
                   2:Nnode(compt$trees[[2]]) + Ntip(compt$trees[[2]]), frame = "none",
                   cex = 0.7, which = "right")

mtext(text = tnames[tr1], side = 3, adj = 0, line = -3 , cex = 2)
mtext(text = tnames[tr2], side = 3, adj = 1, line = -3 , cex = 2)
mtext(text = paste("Node labels =", gsub(pattern = "=", "", tsupp)), side = 1, adj = 0, line = -3 , cex = 2)

dev.off()
