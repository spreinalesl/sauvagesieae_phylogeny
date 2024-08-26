###################################################
## This script read the ML tree and extract
## info about sequence identifications
## of each paralog (e.g. .main, .0, .1) to write
## a map file for Astral-Pro analyses
###################################################

## Required Packages
library(ape)

## Define paths
fpath <- "~/projects/rbgk/users_area/sreinales/results/sauvagesieae/trees/all_paralogs/iqtree/BScollapsed_trees/"
rpath <- "~/projects/rbgk/users_area/sreinales/results/sauvagesieae/trees/paralogs_map/"

## Tree pattern
tpattern1 <- "BS33_trimal_"
tpattern2 <- "_paralogs_all_align.fas.treefile"

## List tree files
tree_files <- list.files(fpath, pattern = tpattern1)
head(tree_files)

## Iterate gene trees
seq_name <- vector()
species_name <- vector()

for(i in tree_files){
  tree <- read.tree(paste(fpath, i, sep = ""))
  
  ## Test for paralog names
  gnames <- grep(pattern = "\\.[0-9]|main", x = tree$tip.label, value = T)
  
  if(length(gnames) > 0){
    seq_name <- sort(gnames)
    species_name <- gsub(pattern = "\\.[0-9]|\\.main", x = seq_name, replacement = "") 
  
  }else{
    seq_name <- sort(tree$tip.label)
    species_name <- seq_name 
  }
  
  ## Extract gene name
  geneID <- gsub(pattern = tpattern1, i, replacement = "")
  geneID <- gsub(pattern = tpattern2, geneID, replacement = "")
  
  ## Write a table
  writeLines(text = paste(seq_name, species_name, sep = "\t"), sep = "\n",
             con = paste(rpath, geneID, "_multiseqID.map", sep = ""))
  
  cat("Processing gen", geneID, "\n")
}
