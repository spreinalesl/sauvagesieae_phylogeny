#######################################
## Plot gene trees and species trees
#######################################

## Tree paths
tpath <- "~/Documents/plot_trees"
tfiles <- list.files(pattern = ".tre$", tpath)
tpattern <- "rooted_sauv_"
tpattern1 <- "_BS33_sptree_annot.tre"

## Data path
dpath <- "~/Dropbox/Sandra/Manuscritos/NGS_phyloOchn/analyses/results/hybpiper/assembly_stats/hybpiper_all_summary_ed.csv" 

#####################################
## Plot species tree
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

## Read table with information about gene recovery and sequence length
sdata <- read.table(dpath, header = T, sep = ";", dec = ",")

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
#pdf(paste(tpath, "/", "wASTRAL_all_trees.pdf", sep = ""), width = 15, height = 30) # All plots for each tree in the same page
tree_summ_q <- data.frame()
tree_summ_p <- data.frame()

for(i in 1:length(tfiles)){
  
  ## Tree name
  tname <- gsub(pattern = tpattern, replacement = "", tfiles[i])
  tname <- gsub(pattern = tpattern1, replacement = "", tname)
  
  ## Plot names and title
  ptitle <- paste("wASTRAL", tname, sep = "_")
  
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
           ggtitle(paste("Node labels:", ptitle))
  
  ## Plot a basic tree
  t_all <- ggtree(mtree@phylo, ladderize = T, branch.length = "none") +
            geom_tiplab(size = 2, hjust=-.2) +
            xlim_tree(70) +
            ggtitle(paste("Basic tree:", ptitle))
  
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
              ggtitle(paste("Q scores tree:", ptitle))
  
  # Local probability LPP
  P_t_all <- as.data.frame(P_t_all)
  
  t_textp <- ggtree(mtree@phylo, ladderize = T, branch.length = "none") %<+% P_t_all + 
              geom_nodelab(aes(x = branch, label = LP), vjust = -0.5, size = 2) +
              geom_tiplab(size = 2.5, hjust = -0.05) + xlim_tree(70) + 
              ggtitle(paste("Local probability tree:", ptitle))
  
  ## Tree with Q scores as pies and LPP as annotations on the branches
  t_complete <- inset(t_textp, pies_t_all, width = 0.04, height = 0.04)
  
  ## tcomplete plus recovery info as heatmaps
  t_hmap_l <- gheatmap(t_complete, df2, offset = 14, width = 0.1, colnames = F) +
                      scale_fill_viridis_c(option = "C", name = "TotalSeqLenght", direction = -1)
  
  t_hmap_g <- t_hmap_l + new_scale_fill()
  
  t_hmap_g <- gheatmap(t_hmap_g, df1, offset = 10, width = 0.1, colnames = F) +
                     scale_fill_viridis_c(option = "D", name = "NumGenesWithSeqs", direction = -1)
  
  t_hmap_p <- t_hmap_g + new_scale_fill()
  
  t_hmap_p <- gheatmap(t_hmap_p, df3, offset = 18, width = 0.1, colnames = F) +
                    scale_fill_viridis_c(option = "A", name = "PctParalogs", direction = -1)
  
  ## Plot one tree per page
  #plots <- list(t_nl, t_textq, t_textp, t_pies, t_complete, t_hmap_p)
  #ggsave(filename = paste(tpath, "/", ptitle, ".pdf", sep = ""), 
  #       plot = marrangeGrob(plots, nrow = 1, ncol = 1),
  #       width = 15, height = 20)
  
  ## Alternative arrange all plots per dataset in the same page
  #grid.arrange(t_textq, t_complete) # All plots for each tree in the same page
  
  cat("Plotting tree", tname, "\n")
}
#dev.off()

## Export summary table
write.csv(tree_summ_q, file = paste(tpath, "/", "qscores.csv", sep = ""), row.names = F)
write.csv(tree_summ_p, file = paste(tpath, "/", "local_prob.csv", sep = ""), row.names = F)

##########################################################
## Plot summary
## Comparison among methods
##########################################################
tree_summ_q$prop_q2_q3 <- tree_summ_q$Q2/tree_summ_q$Q3

## Plots
## quartet support
q1 <- ggplot(tree_summ_q, aes(tree, Q1)) + 
  geom_boxplot(width = .2, outlier.shape = NA) + 
  geom_jitter(aes(fill = "grey"), cex = 0.5, width = .05, alpha = .1) +
  scale_color_manual(values = "grey") +
  scale_fill_manual(values = "black") +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5, hjust = 1), legend.position = "none")

q2 <- ggplot(tree_summ_q, aes(tree, Q2)) + 
  geom_boxplot(width = .2, outlier.shape = NA) + 
  geom_jitter(aes(fill = "grey"), cex = 0.5, width = .05, alpha = .1) +
  scale_color_manual(values = "grey") +
  scale_fill_manual(values = "black") +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5, hjust = 1), legend.position = "none")

q3 <- ggplot(tree_summ_q, aes(tree, Q3)) + 
  geom_boxplot(width = .2, outlier.shape = NA) + 
  geom_jitter(aes(fill = "grey"), cex = 0.5, width = .05, alpha = .1) +
  scale_color_manual(values = "grey") +
  scale_fill_manual(values = "black") +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5, hjust = 1), legend.position = "none")

prop <- ggplot(tree_summ_q, aes(tree, log(prop_q2_q3, base = 10))) + 
  geom_boxplot(width = .2, outlier.shape = NA) + 
  geom_jitter(aes(fill = "grey"), cex = 0.5, width = .05, alpha = .1) +
  scale_color_manual(values = "grey") +
  scale_fill_manual(values = "black") +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5, hjust = 1), legend.position = "none")

# alternative
#ggplot(tree_summ_q, aes(tree, Q1)) +
# geom_boxplot(width = .2) +
# geom_point(shape = 21, colour = "black", fill = "grey", size = 1, stroke = 0.2, position = "jitter") +
# theme_bw() +
# theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5, hjust = 1), legend.position = "none")

## Export pdf
pdf(paste(tpath, "/", "q_plots.pdf", sep = ""), width = 9, height = 9)
ggarrange(q1, q2, q3, prop,
          #labels = c("A", "B"),
          font.label = list(size = 12, color = "black", face = "bold", family = "serif"),
          ncol = 2, nrow = 2, legend = "none")
dev.off()

## local probability
p1 <- ggplot(tree_summ_p, aes(tree, P1)) + 
  geom_boxplot(width = .2, outlier.shape = NA) + 
  geom_jitter(aes(fill = "grey"), cex = 0.5, width = .05, alpha = .1) +
  scale_color_manual(values = "grey") +
  scale_fill_manual(values = "black") +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5, hjust = 1), legend.position = "none")

p2 <- ggplot(tree_summ_p, aes(tree, P2)) + 
  geom_boxplot(width = .2, outlier.shape = NA) + 
  geom_jitter(aes(fill = "grey"), cex = 0.5, width = .05, alpha = .1) +
  scale_color_manual(values = "grey") +
  scale_fill_manual(values = "black") +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5, hjust = 1), legend.position = "none")

p3 <- ggplot(tree_summ_p, aes(tree, P3)) + 
  geom_boxplot(width = .2, outlier.shape = NA) + 
  geom_jitter(aes(fill = "grey"), cex = 0.5, width = .05, alpha = .1) +
  scale_color_manual(values = "grey") +
  scale_fill_manual(values = "black") +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5, hjust = 1), legend.position = "none")

## Export pdf
pdf(paste(tpath, "/", "p_plots.pdf", sep = ""), width = 9, height = 9)
ggarrange(p1, p2, p3,
          #labels = c("A", "B"),
          font.label = list(size = 12, color = "black", face = "bold", family = "serif"),
          ncol = 2, nrow = 2, legend = "none")
dev.off()

######################################
## Test clades
######################################
## Take care, the function MRCA does not accept (-) in the tip labels
clades <- list(blastemanthus = c("B_gemmiflorus_Prance15495", "B_sp_Amaralsn", "B_gemmiflorus_Castaño3445", "B_grandiflorus_Ducke23783"),
               cespclade = c("F_felicis_Haba445", "F_felicis_Haba18", "R_regalis_Reinales159", "R_splendida_Reinales128"),
               poePLUSwala = c("Poe_pumila_Maas5732", "W_insignis_Kawasaki243", "W_multiflora_Berry5926"),
               asianclade = c("Indov_calophylla_Brooke9671", "Schll_angustifolia_Bramley100471"),
               cortyleria = c("T_silvana_Liesner22006", "Ad_bicarpellata_Huber10214", "Ad_bicarpellata_Pipoly7299"),
               sauvagesiaA = c("S_elegantissima_Lovo428", "S_glandulosa_rb_Cardoso3942", "S_insignis_Cardoso2656"),
               sauvagesiaB = c("S_sprengelii_Oliveira085", "S_aliciae_Sastre2339_sprl", "S_longifolia_Harley56621", "S_linearifolia_Cardenas21665", "S_ramosa_Cardenas44000", "S_guianensis_Maguire32554"))

cladesumm <- data.frame()

for(i in 1:length(tfiles)){
  
  ## Tree name
  tname <- gsub(pattern = tpattern, replacement = "", tfiles[i])
  tname <- gsub(pattern = tpattern1, replacement = "", tname)

  ## Read tree
  mtree <- read.astral(paste(tpath, "/", tfiles[i], sep = ""))
  
  for(j in 1:length(clades)){
   
    ## Extract MCRC info for some important clades
    cnode <- MRCA(as_tibble(mtree), clades[j][[1]])$node
    coffspring <- offspring(as_tibble(mtree), cnode, tiponly = F, self_include = T)
    ctips <- coffspring$label[-which(is.na(coffspring$label))]
    dnodes <- coffspring$node[which(is.na(coffspring$label))]
    
    ## Test
    cat(tname, "- Clade:", names(clades[j]), "- Test nodes:", (cnode %in% dnodes), "\n")
    
    ## Extract qscores
    q1 <- mtree@data[which(mtree@data$node %in% dnodes), c("node", "q1")]
    qscores <- merge(q1, mtree@data[which(mtree@data$node %in% dnodes), c("node", "q2")], by = "node")
    qscores <- merge(qscores, mtree@data[which(mtree@data$node %in% dnodes), c("node", "q3")], by = "node")
    
    ## Extract local probabilities
    p1 <- mtree@data[which(mtree@data$node %in% dnodes), c("node", "pp1")]
    localp <- merge(p1, mtree@data[which(mtree@data$node %in% dnodes), c("node", "pp2")], by = "node")
    localp <- merge(localp, mtree@data[which(mtree@data$node %in% dnodes), c("node", "pp3")], by = "node")
    
    ## Save summary for clades
    ct <- merge(qscores, localp, by = "node")
    ct$tree <- tname
    ct$clade <- names(clades[j])
    ct$tips <- paste0(sort(ctips), collapse = ";")
    cladesumm <- rbind(cladesumm, ct)
    
    cat("Processing clade", names(clades)[j], "\n")
  }
  
  cat("Tree", tfiles[i], "done!", "\n")
}

## Plots
cladesumm_org <- cladesumm
cladesumm$q1 <- as.numeric(cladesumm$q1)
cladesumm$q2 <- as.numeric(cladesumm$q2)
cladesumm$q3 <- as.numeric(cladesumm$q3)

## Export report
write.csv(cladesumm, file = paste(tpath, "/", "clades_sum.csv", sep = ""), row.names = F)

## colors for clades
cols <- c("#fde725", "#90d743", "#0d0887", "#f37819", "#35b779", "#9c179e", "#000004")

## Colors for dataset
#cols2 <- c("#0d0887", "#fde725", "#90d743", "#9c179e")

# by tree by clade
q1_tree <- ggplot(cladesumm, aes(x = tree, y = q1, fill = clade)) + 
             geom_boxplot() +
             labs(fill = "clade") +
             geom_point(position = position_jitterdodge(jitter.width = 0.05), alpha = .4) +
             theme_bw() +
             theme(axis.title.x = element_blank(),
                   axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust = 1),
                   legend.position = "none")

p1_tree <- ggplot(cladesumm, aes(x = tree, y = pp1, fill = clade)) + 
             geom_boxplot() +
             labs(fill = "clade") +
             geom_point(position = position_jitterdodge(jitter.width = 0.05), alpha = .4) +
             theme_bw() +
             theme(axis.title.x = element_blank(),
                   axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust = 1),
                   legend.position = "none")

## Export pdf
pdf(paste(tpath, "/", "trees-clades.pdf", sep = ""), width = 12, height = 10)
ggarrange(q1_tree, p1_tree,
          font.label = list(size = 12, color = "black", face = "bold", family = "serif"),
          ncol = 1, nrow = 2,
          common.legend = TRUE, legend = "bottom")
dev.off()

# by clade by tree
q1_clade <- ggplot(cladesumm, aes(x = clade, y = q1, fill = tree)) + 
             geom_boxplot() +
             labs(fill = "tree") +
             geom_point(position = position_jitterdodge(jitter.width = 0.05), alpha = .4) +
             theme_bw() +
             theme(axis.title.x = element_blank(),
             axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust = 1),
             legend.position = "none")

p1_clade <- ggplot(cladesumm, aes(x = clade, y = pp1, fill = tree)) + 
             geom_boxplot() +
             labs(fill = "tree") +
             geom_point(position = position_jitterdodge(jitter.width = 0.05), alpha = .4) +
             theme_bw() +
             theme(axis.title.x = element_blank(),
             axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust = 1),
             legend.position = "none")

## Export pdf
pdf(paste(tpath, "/", "clades-trees.pdf", sep = ""), width = 12, height = 10)
ggarrange(q1_clade, p1_clade,
          font.label = list(size = 12, color = "black", face = "bold", family = "serif"),
          ncol = 1, nrow = 2,
          common.legend = TRUE, legend = "bottom")
dev.off()