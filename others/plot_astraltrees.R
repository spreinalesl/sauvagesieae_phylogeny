#######################################
## Plot gene trees and species trees
#######################################

## Tree paths
tpath <- "/home/spreinalesl/Documents/wastral/rooted"
tfiles <- list.files(pattern = ".tre$", tpath)
tpattern <- "_sauv_all_exons_sptree_annot.tre"
tpattern2 <- "rooted_"
tpattern3 <- "plus_models_"

#####################################
## Plot species tree
#####################################
## Required packages
library(ggplot2)
library(treeio) ## Installed using BiocManager::install("treeio")
library(ggtree) ## remotes::install_github('YuLab-SMU/ggtree')
library(ggimage) ## To plot pies on the tree
library(gridExtra) ## To arrange plots into one pdf
library(ggpubr) ## Plot one plot per page into one pdf

## Plot names and title
ptitle <- "wASTRAL_all_exon"

## Function to format quartet information from astral tree annotation
# modified from: https://raw.githubusercontent.com/sidonieB/scripts/master/plot_Astral_trees_v2.R
# The normalized quartet score QS represents the proportion of input gene tree quartet trees satisfied
# by the species tree). The higher this number, the less discordant your gene trees are.
# qs are values show quartet support for the main topology (LR|SO): q1, first alternative (RS|LO):q2, and second alternative (RO|LS):q3.
QS_calc <- function(AS_A) {
  Q1 <- as.numeric(AS_A@data$q1) * 100
  Q <- as.data.frame(Q1)
  Q$Q2 <- as.numeric(AS_A@data$q2) * 100
  Q$Q3 <- as.numeric(AS_A@data$q3) * 100
  Q$node <- AS_A@data$node
  return(Q)
}

# LPP (local posterior probability): the probability this branch is the true branch given the set of gene trees (computed based on the quartet score and assuming ILS).
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
  
  ## Plot basic tree
  t_all <- ggtree(mtree@phylo, ladderize = T, branch.length = "none") +
            geom_tiplab(size = 2, hjust=-.2) +
            xlim_tree(70) +
            ggtitle(paste("Basic tree:", ptitle, tname))
  
  ## Make pie charts (takes some time if tree is big)
  pies_t_all <- nodepie(Q_t_all, cols = 1:3, color = c(Q1='darkblue', Q2='green', Q3='grey'), alpha = 0.8) # change style as needed, alpha is for transparency

  ## Insert pies to the basic tree
  t_pies <- inset(t_all, pies_t_all, width = 0.04, height = 0.04)
  
  ## Make barplots instead of pies (It doesn't work)
  bars <- nodebar(Q_t_all, cols = 1:3, position = 'dodge', color = c(Q1 = 'darkblue', Q2 = 'green', Q3 = 'gray'), alpha = 0.8)
  t_bars <- inset(t_all, bars, x = 'node', width = 0.04, height = 0.04)
  
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
  
  ## Q scores and LPP
  t_complete <- inset(t_textp, pies_t_all, width = 0.04, height = 0.04)
  
  ## arrange all the plots
  #grid.arrange(t_nl, t_textq, t_pies, t_complete) # All plots for each tree in the same page
  #multi.page <- ggarrange(t_nl, t_textq, t_textp, t_pies, t_complete, nrow = 1, ncol = 1) # for one plot per page
  #ggexport(multi.page, filename = paste(tpath, "/", ptitle, "_", tname, ".pdf", sep = ""), width = 480, height = 520)
  
  plots <- list(t_nl, t_textq, t_textp, t_pies, t_complete)
  ggsave(filename = paste(tpath, "/", ptitle, "_", tname, ".pdf", sep = ""), 
         plot = marrangeGrob(plots, nrow = 1, ncol = 1),
         width = 15, height = 20)
  
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
tree_summ_q$prop <- tree_summ_q$Q1/(tree_summ_q$Q2 + tree_summ_q$Q3) 
  
## Plots
p1 <- ggplot(tree_summ_p, aes(tree, P1, colour = "grey")) + 
       geom_boxplot(width = .2, outlier.shape = NA) + 
       geom_jitter(aes(fill = "grey"), width = .05, alpha = .3) +
       stat_halfeye(aes(fill = "grey"), adjust = .9, width = .7, .width = 0, justification = -.2, point_colour = NA, alpha = .7) +
       scale_color_manual(values = "grey") +
       scale_fill_manual(values = "grey") +
       theme_bw() +
       theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5, hjust = 1))

p2 <- ggplot(tree_summ_p, aes(tree, P2, colour = "grey")) + 
       geom_boxplot(width = .2, outlier.shape = NA) + 
       geom_jitter(aes(fill = "grey"), width = .05, alpha = .3) +
       stat_halfeye(aes(fill = "grey"), adjust = .9, width = .7, .width = 0, justification = -.2, point_colour = NA, alpha = .7) +
       scale_color_manual(values = "grey") +
       scale_fill_manual(values = "grey") +
       theme_bw() +
       theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5, hjust = 1))

p3 <- ggplot(tree_summ_p, aes(tree, P3, colour = "grey")) + 
       geom_boxplot(width = .2, outlier.shape = NA) + 
       geom_jitter(aes(fill = "grey"), width = .05, alpha = .3) +
       stat_halfeye(aes(fill = "grey"), adjust = .9, width = .7, .width = 0, justification = -.2, point_colour = NA, alpha = .7) +
       scale_color_manual(values = "grey") +
       scale_fill_manual(values = "grey") +
       theme_bw() +
       theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5, hjust = 1))

prop <- ggplot(tree_summ_q, aes(tree, prop, colour = "grey")) + 
         geom_boxplot(width = .2, outlier.shape = NA) + 
         geom_jitter(aes(fill = "grey"), width = .05, alpha = .3) +
         stat_halfeye(aes(fill = "grey"), adjust = .9, width = .7, .width = 0, justification = -.2, point_colour = NA, alpha = .7) +
         scale_color_manual(values = "grey") +
         scale_fill_manual(values = "grey") +
         theme_bw() +
         theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5, hjust = 1))


## Pdf
pdf(paste(tpath, "/", "q_plot.pdf", sep = ""), width = 9, height = 9)
 ggarrange(q1, q2, q3, prop,
           #labels = c("A", "B"),
           font.label = list(size = 12, color = "black", face = "bold", family = "serif"),
           ncol = 2, nrow = 2, legend = "none")
dev.off()

pdf(paste(tpath, "/", "p_plot.pdf", sep = ""), width = 9, height = 9)
 ggarrange(p1, p2, p3,
           #labels = c("A", "B"),
           font.label = list(size = 12, color = "black", face = "bold", family = "serif"),
           ncol = 2, nrow = 2, legend = "none")
dev.off()

######################################
## Test clades
######################################
## Take care, the function MRCA does not accept (-) in the tip labels
clades <- list(blastemanthus = c("B_sp_Amaralsn", "B_grandiflorus_Ducke23783"),
               cespclade = c("F_felicis_Haba445", "F_felicis_Haba18", "R_regalis_Reinales159"),
               poePLUSwala = c("Poe_pumila_Maas5732", "W_insignis_Kawasaki243"),
               asianclade = c("Indov_calophylla_Brooke9671", "Sch_henningsii_Hoogland8954", "Schll_angustifolia_Bramley100471"),
               cortyleria = c("T_silvana_Liesner22006", "Ad_bicarpellata_Huber10214"),
               sauvagesiaA = c("S_glandulosa_rb_Cardoso3942", "S_alpestris_Mota993", "S_insignis_Cardoso2656"),
               sauvagesiaB = c("S_sprengelii_Oliveira085", "S_longifolia_Harley56621", "S_ramosa_Cardenas44000"))

cladesumm <- data.frame()

for(i in 1:length(tfiles)){
  
  ## Tree name
  tname <- gsub(pattern = tpattern, replacement = "", tfiles[i])
  tname <- gsub(pattern = tpattern2, replacement = "", tname)
  tname <- gsub(pattern = tpattern3, replacement = "", tname)
  
  ## Read tree
  mtree <- read.astral(paste(tpath, "/", tfiles[i], sep = ""))
  
  for(j in 1:length(clades)){
    ## Extract MCRC info for some important clades
    cnode <- MRCA(as_tibble(mtree), clades[j][[1]])$node
    ctips <- offspring(as_tibble(mtree), cnode, tiponly = F, self_include = T)$label
    ctips <- ctips[-which(is.na(ctips))]
    
    ## Extract qscores
    q1 <- mtree@data$q1[which(mtree@data$node == cnode)] 
    qscores <- as.data.frame(q1)
    qscores$q2 <- mtree@data$q2[which(mtree@data$node == cnode)]
    qscores$q3 <- mtree@data$q3[which(mtree@data$node == cnode)]
    
    ## Extract local probabilities
    p1 <- mtree@data$pp1[which(mtree@data$node == cnode)]
    localp <- as.data.frame(p1)
    localp$p2 <- mtree@data$pp2[which(mtree@data$node == cnode)]
    localp$p3 <- mtree@data$pp3[which(mtree@data$node == cnode)]
    
    ## Save summary for clades
    ct <- cbind(qscores, localp)
    ct$node <- cnode
    ct$clade <- names(clades)[j]
    ct$tree <- tname
    ct$tips <- paste0(sort(ctips), collapse = ";")
    cladesumm <- rbind(cladesumm, ct)
    
    cat("Processing clade", names(clades)[j], "\n")
  }

  cat("Tree", tfiles[i], "done!", "\n")
}

## Export report
write.csv(cladesumm, file = paste(tpath, "/", "clades_sum.csv", sep = ""), row.names = F)

## Plot summary
cols <- c("#fde725", "#90d743", "#0d0887", "#f37819", "#35b779", "#9c179e", "#000004")

cladesumm$q1 <- round(as.numeric(cladesumm$q1), digits = 2)
cladesumm$p1 <- round(as.numeric(cladesumm$p1), digits = 2)

q1 <- ggplot(cladesumm, aes(tree, q1, colour = clade)) + 
        geom_boxplot(aes(colour = clade)) +
        geom_point(aes(colour = clade), alpha = .4,
             position = position_jitterdodge(jitter.width = 0.08)) +
        scale_color_manual(values = cols) +
        scale_fill_manual(values = cols) +
        theme_bw() +
        theme(axis.title.x = element_blank(),
              axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust = 1),
              legend.position = "none")

p1 <- ggplot(cladesumm, aes(tree, p1, colour = clade)) + 
       geom_boxplot(aes(colour = clade)) +
       geom_point(aes(colour = clade), alpha = .4,
             position = position_jitterdodge(jitter.width = 0.08)) +
       scale_color_manual(values = cols) +
       scale_fill_manual(values = cols) +
       theme_bw() +
       theme(axis.title.x = element_blank(),
             axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust = 1), 
             legend.position = "none")

## Plot for tree inference
cladesumm$treeinf[grep(pattern = "raxml", cladesumm$tree)] <- "raxml"
cladesumm$treeinf[grep(pattern = "iqtree", cladesumm$tree)] <- "iqtree"

ti <- ggplot(cladesumm, aes(treeinf, q1, colour = clade)) + 
       geom_boxplot(aes(colour = clade)) +
       geom_point(aes(colour = clade), alpha = .4,
             position = position_jitterdodge(jitter.width = 0.08)) +
       scale_color_manual(values = cols) +
       scale_fill_manual(values = cols) +
       theme_bw() +
       theme(axis.title.x = element_blank(),
             axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust = 1),
             legend.position = "none")

## Plot for alignment trimming method
cladesumm$trim[grep(pattern = "trimal", cladesumm$tree)] <- "trimal"
cladesumm$trim[grep(pattern = "taper", cladesumm$tree)] <- "taper"

trim <- ggplot(cladesumm, aes(trim, q1, colour = clade)) + 
         geom_boxplot(aes(colour = clade)) +
         geom_point(aes(colour = clade), alpha = .4,
              position = position_jitterdodge(jitter.width = 0.08)) +
         scale_color_manual(values = cols) +
         scale_fill_manual(values = cols) +
         theme_bw() +
         theme(axis.title.x = element_blank(),
               axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust = 1),
               legend.position = "none")

## Plot for BS support
cladesumm$bs[which(!grepl(pattern = "BS033", cladesumm$tree))] <- "BS10"
cladesumm$bs[grep(pattern = "BS033", cladesumm$tree)] <- "BS33"

bs <- ggplot(cladesumm, aes(bs, q1, colour = clade)) + 
       geom_boxplot(aes(colour = clade)) +
       geom_point(aes(colour = clade), alpha = .4,
             position = position_jitterdodge(jitter.width = 0.08)) +
       scale_color_manual(values = cols) +
       scale_fill_manual(values = cols) +
       theme_bw() +
       theme(axis.title.x = element_blank(),
             axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust = 1),
             legend.position = "none")

## General for tree inference
treeinfg <- ggplot(cladesumm, aes(treeinf, q1, colour = treeinf)) + 
             geom_boxplot(width = .2, outlier.shape = NA) + 
             geom_point(aes(colour = treeinf), alpha = .4,
                  position = position_jitterdodge(jitter.width = 0.08)) +           
             theme_bw() +
             theme(axis.title.x = element_blank(),
                   axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust = 1),
                   legend.position = "none")

## Export pdf
pdf(paste(tpath, "/", "clade_summ.pdf", sep = ""), width = 9, height = 7)
  ggarrange(q1, p1,
           #labels = c("A", "B", "C"),
            font.label = list(size = 12, color = "black", face = "bold", family = "serif"),
            ncol = 2, nrow = 1,
            common.legend = TRUE, legend = "bottom")
dev.off()

pdf(paste(tpath, "/", "clade_summ2.pdf", sep = ""), width = 10, height = 7)
ggarrange(ti, trim, bs, 
          #labels = c("A", "B", "C"),
          font.label = list(size = 12, color = "black", face = "bold", family = "serif"),
          ncol = 3, nrow = 1,
          common.legend = TRUE, legend = "bottom")
dev.off()

pdf(paste(tpath, "/", "clade_summ3.pdf", sep = ""), width = 7, height = 7)
ggarrange(treeinfg,
          #labels = c("A", "B", "C"),
          font.label = list(size = 12, color = "black", face = "bold", family = "serif"),
          ncol = 1, nrow = 1,
          common.legend = F)
dev.off()
