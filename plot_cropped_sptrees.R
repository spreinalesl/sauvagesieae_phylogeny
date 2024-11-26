#####################################
## Plot trees
#####################################
## Required packages
library(treeio) ## Installed using BiocManager::install("treeio")
library(ggplot2)
library(ggtree) ## Installed using BiocManager::install("ggtree")
library(ggimage) ## To plot pies on the tree
library(gridExtra) ## To arrange all plots into one pdf
library(ggpubr) ## To plot one tree per page into one pdf
library(ggdist) ## To plot densities in the boxplots
library(ggnewscale) ## To include two heatmaps in the same plot, adding a new scale

## Tree paths
tpath <- "~/Documents/plot_trees"
tfiles <- list.files(pattern = ".tre$", tpath)
tpattern <- "rooted_sauv_"
tpattern1 <- "_BS33_sptree_annot.tre"

## Data path
dpath <- "~/Documents/plot_trees/hybpiper_all_summary_ed.csv" 

## Read table with information about gene recovery and sequence length
sdata <- read.csv(dpath, header = T)

# Extract info about number of gene recovered
df1 <- tapply(X = sdata$genesWseqsFasta, INDEX = sdata$sppName, FUN = sum)
df1 <- as.data.frame(df1, row.names = names(df1))
colnames(df1) <- "GenesWithSeqs"

# Extract info about seq length
df2 <- tapply(X = sdata$total_captured_bp, INDEX = sdata$sppName, FUN = sum)
df2 <- as.data.frame(df2, row.names = names(df2))
colnames(df2) <- "RecoveryLenght"

# Extract info about paralogs
df3 <- tapply(X = sdata$TotalParalogWarnings, INDEX = sdata$sppName, FUN = sum)
df3 <- as.data.frame(round((df3/df1)*100, digits = 1), row.names = names(df3))
colnames(df3) <- "PctParalogs"

# Extract info about AD
df4 <- tapply(X = sdata$AD, INDEX = sdata$sppName, FUN = mean)
df4 <- as.data.frame(df4, row.names = names(df4))
colnames(df4) <- "AD"

# Extract info about LH
df5 <- tapply(X = sdata$LH, INDEX = sdata$sppName, FUN = mean)
df5 <- as.data.frame(df5, row.names = names(df5))
colnames(df5) <- "LH"

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
tree_summ_q <- data.frame()
tree_summ_p <- data.frame()


pdf(paste(tpath, "/", "croppedtrees.pdf", sep = ""))#, width = 12, height = 20) # All plots for each tree in the same page
for(i in 1:length(tfiles)){
  
  ## Tree name
  tname <- gsub(pattern = tpattern, replacement = "", tfiles[i])
  tname <- gsub(pattern = tpattern1, replacement = "", tname)
  
  ## Plot names and title
  ptitle <- paste("MSCtree:", gsub(pattern = "[[:digit:]]", tname, replacement = ""))
  
  ## Read tree
  mtree <- read.astral(paste(tpath, "/", tfiles[i], sep = ""))
  is.rooted(mtree)
  
  ## Crop tips
  tipsTOkeep <- c("S_sprengelii_Oliveira085", "S_linearifolia_Reinales208", "S_ramosa_Cardenas44000", "Ak_muriculata_Maguire60447", "S_vellozii_Zavatin158", "S_semicylindrifolia_Cardoso3844",
                  "S_insignis_Cardoso2656", "Ad_bicarpellata_Pipoly7299", "T_silvana_Liesner22006", "W_insignis_Arbeláez1109", "Poe_pumila_Hoffmann1636", "Indov_calophylla_Beaman11839",
                  "Indos_involucrata_Palane6484", "N_serrata_Khairuddin31754", "Euth_minor_Lanas789", "Schll_angustifolia_Bramley100471", "R_regalis_Idrobo871", "G_antioquiensis_Reinales058",
                  "K_disticha_Neill15849", "C_spathulata_Cardoso3396", "F_felicis_Haba18", "B_grandiflorus_Ducke23783", "Ochna_membranacea_Konomou261", "L_mogolensis_Freire-FierrosCFCR12462",
                  "P_auriculata_Pipoly6867", "Bonn_stricta_Benko-Iseppon1726")
  
  tipsTOrm <- mtree@phylo$tip.label[!(mtree@phylo$tip.label %in% tipsTOkeep)]
  
  #tipsTOrm <- c("S_amoena_CampbellP22564", "T_spectabilis_Nee31125", "T_apiculata_Liesner24854", "S_africana_Nawa119",
  #              "T_pendula_Maguire42286", "T_spathulata_Huber4305", "T_floribunda_Steyermark124565", "T_tremuloidea_Croat59545",
  #              "T_spathulata_Liesner17587", "T_breweriana_Steyermark124326", "T_linearis_Tate555", "Ochna_serrulata_OneKPJLM2013-53",
  #              "Camp_reticulatum_Cable3433", "O_castaneifolia_Sasaki1809", "Elv_calophyllea_Kubitzki84-349", "L_octandra_Forzza3713",
  #              "L_mogolensis_Freire-FierrosCFCR12462", "P_auriculata_Chase3088", "Hype_perforatum_OneKP132735")
  
  #tipsTOkeep <- mtree@phylo$tip.label[!(mtree@phylo$tip.label %in% tipsTOrm)]
  
  if(length(unique(tipsTOkeep %in% mtree@phylo$tip.label)) > 1) { stop("There is a tip in the list wich is not in the tree", tfiles[i], "\n") }
  
  mtree <- tidytree::keep.tip(mtree, tipsTOkeep)
  
  ## Extract quartet and local probability values using QS_calc and  LP_calc functions
  Q_t_all <- QS_calc(AS_A = mtree)
  P_t_all <- LP_calc(LP_A = mtree)
  
  ## Make a label with the three QS separated by underscores
  for (k in 1:length(Q_t_all$node)){
    Q_t_all$QS[k] <- paste(round(Q_t_all$Q1[k]), "_", round(Q_t_all$Q2[k]), "_", round(Q_t_all$Q3[k]), sep = "")
  }
  
  for (j in 1:length(P_t_all$node)){
    P_t_all$LP[j] <- paste(round(P_t_all$P1[j], digits = 1), "_", round(P_t_all$P2[j], digits = 1), "_", round(P_t_all$P3[j], digits = 1), sep = "")
  }
  
  ## Save summary
  Q_t_all$tree <- tname
  P_t_all$tree <- tname
  tree_summ_q <- rbind(tree_summ_q, Q_t_all)
  tree_summ_p <- rbind(tree_summ_p, P_t_all)
  
  ## Plot a basic tree
  t_all <- ggtree(mtree@phylo, ladderize = T, branch.length = "none") +
            geom_tiplab(size = 6, hjust=-.2) +
            xlim_tree(70)
  
  ## Make pie charts (takes some time if tree is big)
  pies_t_all <- nodepie(Q_t_all, cols = 1:3, color = c(Q1='darkblue', Q2='green', Q3='grey'), alpha = 0.8) # change style as needed, alpha is for transparency

  ## Insert pies to the basic tree
  t_pies <- inset(t_all, pies_t_all, width = 0.3, height = 0.3)
  
  ## Q scores
  Q_t_all <- as.data.frame(Q_t_all)
  
  ## create plot including text annotations
  t_textq <- ggtree(mtree@phylo, ladderize = T, branch.length = "none") %<+% Q_t_all + 
              geom_nodelab(aes(x = branch, label = QS), vjust = -0.5, size = 2.5) +
              geom_tiplab(size = 6, hjust = -0.05) + xlim_tree(70)
  
  # Local probability LPP
  P_t_all <- as.data.frame(P_t_all)
  #P_t_all <- data.frame(P1 = round(P_t_all$P1, digits = 1), node = P_t_all$node, tree = P_t_all$tree)
  #P_t_all$P1[which(P_t_all$P1 == 1.0)] <- NA 
  
  t_textp <- ggtree(mtree@phylo, ladderize = T, branch.length = "none") %<+% P_t_all + 
              geom_nodelab(aes(x = branch, label = LP), vjust = -0.5, hjust = 1.2, size = 3) +
              geom_tiplab(size = 2.5, hjust = -0.05) + xlim_tree(70) +
              ggtitle(ptitle)
  
  ## Tree with Q scores as pies and LPP as annotations on the branches
  t_complete <- inset(t_textp, pies_t_all, width = 0.3, height = 0.3)
  
  ## tcomplete plus recovery info as heatmaps
  #t_hmap_g <- gheatmap(t_complete, df1, offset = 10, width = 0.15, colnames = F) +
  #              scale_fill_viridis_c(option = "D", name = "NumGenesWithSeqs", direction = -1)
  
  #t_hmap_l <- t_hmap_g + new_scale_fill()
  
  #t_hmap_l <- gheatmap(t_hmap_l, df2, offset = 14, width = 0.15, colnames = F) +
  #                    scale_fill_viridis_c(option = "D", name = "TotalSeqLenght", direction = -1)
  
  #t_hmap_p <- t_hmap_l + new_scale_fill()
  
  #t_hmap_p <- gheatmap(t_hmap_p, df3, offset = 18, width = 0.15, colnames = F) +
  #                  scale_fill_viridis_c(option = "D", name = "PctParalogs", direction = -1)
  
  #t_hmap_ad <- t_hmap_p + new_scale_fill()
  
  #t_hmap_ad <- gheatmap(t_hmap_ad, df4, offset = 22, width = 0.15, colnames = F) +
  #                scale_fill_viridis_c(option = "D", name = "AD", direction = -1)
  
  #t_hmap_lh <- t_hmap_ad + new_scale_fill()
  
  #t_hmap_lh <- gheatmap(t_hmap_lh, df5, offset = 26, width = 0.15, colnames = F) +
  #                 scale_fill_viridis_c(option = "D", name = "LH", direction = -1)
  
  ## Plot
  #plot(t_hmap_lh)
  plot(t_complete)
  
  cat("Plotting tree", tname, "\n")
}
dev.off()
