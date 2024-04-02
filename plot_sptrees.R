#######################################
## Plot gene trees and species trees
#######################################

## Tree paths
tpath <- "~/Documents"
tfiles <- list.files(pattern = ".tre$", tpath)
tpattern <- "_sauv_all_exons"
tpattern1 <- "_sptree_annot.tre"
tpattern2 <- "rooted_"
tpattern3 <- "plus_models_"

## Data path
dpath <- "~/Documents/hybpiper_all_summary_ed.csv" 

#####################################
## Plot species tree
#####################################
## Required packages
library(treeio) ## Installed using BiocManager::install("treeio")
library(ggtree) ## Installed using remotes::install_github('YuLab-SMU/ggtree')
library(ggimage) ## To plot pies on the tree
library(gridExtra) ## To arrange all plots into one pdf
library(ggpubr) ## To plot one tree per page into one pdf
library(ggdist) ## To plot densities in the boxplots
library(ggnewscale) ## To include two heatmaps in the same plot, adding a new scale

## Plot names and title
ptitle <- "wASTRAL_all_exons"

## Read table with information about gene recovery and sequence length
sdata <- read.table(dpath, header = T, sep = ",")

# Extract info about number of gene recovered
df1 <- tapply(X = sdata$GenesWithSeqs, INDEX = sdata$abbrev, FUN = sum)
df1 <- as.data.frame(df1, row.names = names(df1))
colnames(df1) <- "GenesWithSeqs"

# Extract info about seq length
df2 <- tapply(X = sdata$RecoveryLenght, INDEX = sdata$abbrev, FUN = sum)
df2 <- as.data.frame(df2, row.names = names(df2))
colnames(df2) <- "RecoveryLenght"

## Function to format quartet information from astral tree annotation
# The normalized quartet score QS represents the proportion of input gene tree quartet trees satisfied
# by the species tree). The higher this number, the less discordant your gene trees are.
# modified from: https://raw.githubusercontent.com/sidonieB/scripts/master/plot_Astral_trees_v2.R
QS_calc <- function(AS_A) {
  Q1 <- as.numeric(AS_A@data$q1) * 100
  Q <- as.data.frame(Q1)
  Q$Q2 <- as.numeric(AS_A@data$q2) * 100
  Q$Q3 <- as.numeric(AS_A@data$q3) * 100
  Q$node <- AS_A@data$node
  return(Q)
}

# LPP (local posterior probability): the probability this branch is the true branch given the set
# of gene trees (computed based on the quartet score and assuming ILS).
LP_calc <- function(LP_A) {
  P1 <- as.numeric(LP_A@data$pp1)
  P <- as.data.frame(P1)
  P$P2 <- as.numeric(LP_A@data$pp2)
  P$P3 <- as.numeric(LP_A@data$pp3)
  P$node <- LP_A@data$node
  return(P)
}

## Iterate tree files, read trees from Astral and export plots to pdf
#pdf(paste(tpath, "/", "ASTRAL_all_exons.pdf", sep = ""), width = 7, height = 9) # All plots for each tree in the same page
tree_summ_q <- data.frame()
tree_summ_p <- data.frame()

for(i in 1:length(tfiles)){
  
  ## Tree name
  tname <- gsub(pattern = tpattern, replacement = "", tfiles[i])
  tname <- gsub(pattern = tpattern1, replacement = "", tname)
  tname <- gsub(pattern = tpattern2, replacement = "", tname)
  tname <- gsub(pattern = tpattern3, replacement = "", tname)
  
  ## Read tree
  mtree <- read.astral(paste(tpath, "/", tfiles[i], sep = ""))
  is.rooted(mtree)
  
  ## Extract quartet and local probability values using QS_calc and  LP_calc functions
  Q_t_all <- QS_calc(AS_A = mtree)
  P_t_all <- LP_calc(LP_A = mtree)
  
  ## Make a label with the three QS separated by underscores
  for (i in 1:length(Q_t_all$node)){
    Q_t_all$QS[i] <- paste(round(Q_t_all$Q1[i]), "_", round(Q_t_all$Q2[i]), "_", round(Q_t_all$Q3[i]), sep = "")
  }
  
  for (j in 1:length(P_t_all$node)){
    P_t_all$LP[j] <- paste(round(P_t_all$P1[j], digits = 1), "_", round(P_t_all$P2[j], digits = 1), "_", round(P_t_all$P3[j], digits = 1), sep = "")
  }
  
  ## Save summary
  Q_t_all$tree <- tname
  P_t_all$tree <- tname
  tree_summ_q <- rbind(tree_summ_q, Q_t_all)
  tree_summ_p <- rbind(tree_summ_p, P_t_all)
  
  ## Plot tree with node labels
  t_nl <- ggtree(mtree@phylo, ladderize = T, branch.length = "none") +
           geom_text(aes(label = node), size = 2, hjust=-.3) +
           geom_tiplab(size = 2, hjust=-.2) +
           xlim_tree(70) +
           ggtitle(paste("Node labels:", ptitle, tname))
  
  ## Plot a basic tree
  t_all <- ggtree(mtree@phylo, ladderize = T, branch.length = "none") +
            geom_tiplab(size = 2, hjust=-.2) +
            xlim_tree(70) +
            ggtitle(paste("Basic tree:", ptitle, tname))
  
  ## Make pie charts (takes some time if tree is big)
  pies_t_all <- nodepie(Q_t_all, cols = 1:3, color = c(Q1='darkblue', Q2='green', Q3='grey'), alpha = 0.8) # change style as needed, alpha is for transparency

  ## Insert pies to the basic tree
  t_pies <- inset(t_all, pies_t_all, width = 0.04, height = 0.04)
  
  ## create plot including text annotations
  # Q scores
  Q_t_all <- as.data.frame(Q_t_all)
  
  t_textq <- ggtree(mtree@phylo, ladderize = T, branch.length = "none") %<+% Q_t_all + 
              geom_nodelab(aes(x = branch, label = QS), vjust = -0.5, size = 2) +
              geom_tiplab(size = 2.5, hjust = -0.05) + xlim_tree(70) + 
              ggtitle(paste("Q scores tree:", ptitle, tname))
  
  # Local probability LPP
  P_t_all <- as.data.frame(P_t_all)
  
  t_textp <- ggtree(mtree@phylo, ladderize = T, branch.length = "none") %<+% P_t_all + 
              geom_nodelab(aes(x = branch, label = LP), vjust = -0.5, size = 2) +
              geom_tiplab(size = 2.5, hjust = -0.05) + xlim_tree(70) + 
              ggtitle(paste("Local probability tree:", ptitle, tname))
  
  ## Tree with Q scores as pies and LPP as annotations on the branches
  t_complete <- inset(t_textp, pies_t_all, width = 0.04, height = 0.04)
  
  ## tcomplete plus recovery info as heatmaps
  t_hmap_l <- gheatmap(t_complete, df2, offset = 14, width = 0.1, colnames = F) +
                       scale_fill_viridis_c(option = "C", name = "TotalSeqLenght", direction = -1)
  
  t_hmap_g <- t_hmap_l + new_scale_fill()
  
  t_hmap_g <- gheatmap(t_hmap_g, df1, offset = 10, width = 0.1, colnames = F) +
                       scale_fill_viridis_c(option = "D", name = "NumGenesWithSeqs", direction = -1)
  
  ## arrange all the plots in the same page
  #grid.arrange(t_nl, t_textq, t_pies, t_complete) # All plots for each tree in the same page
  #multi.page <- ggarrange(t_nl, t_textq, t_textp, t_pies, t_complete, nrow = 1, ncol = 1) # for one plot per page
  #ggexport(multi.page, filename = paste(tpath, "/", ptitle, "_", tname, ".pdf", sep = ""), width = 480, height = 520)
  
  plots <- list(t_nl, t_textq, t_textp, t_pies, t_complete, t_hmap_g)
  ggsave(filename = paste(tpath, "/", ptitle, "_", tname, ".pdf", sep = ""), 
         plot = marrangeGrob(plots, nrow = 1, ncol = 1),
         width = 15, height = 20)
  
  cat("Plotting tree", tname, "\n")
}
#dev.off()

## Export summary table
write.csv(tree_summ_q, file = paste(tpath, "/", "qscores.csv", sep = ""), row.names = F)
write.csv(tree_summ_p, file = paste(tpath, "/", "local_prob.csv", sep = ""), row.names = F)
