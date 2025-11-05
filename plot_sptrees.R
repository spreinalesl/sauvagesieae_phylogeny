#######################################################
## PLOT SPECIES TREES
## Due to the continuous changes in ggplot2,
## this code does not work with ggplot2 v.4.0.
## So, please use an environment with
## the following versions of the required packages
## to run this code
## treeio v.1.10.0
## ggplot2 v.3.3.3 
## ggtree v.2.0.0
## ggimage v.0.2.8
## ggnewscale v.0.4.4
## gridExtra v.2.3
#######################################################

## Required packages
library(treeio) ## reads ASTRAL trees 
library(ggplot2) ## Manages trees
library(ggtree) ## plots ASTRAL trees
library(ggimage) ## plots pies on the tree
library(ggnewscale) ## includes two heatmaps in the same plot, adding a new scale
library(gridExtra) ## arranges all plots into one pdf

setwd("")

## Tree paths
tpath <- "plot_trees"
tfiles <- list.files(pattern = ".tre$", tpath)
tpattern <- "rooted_sauv_"
tpattern1 <- "_BS33_sptree_annot.tre"
nmpath <- "supplementary_tables/TableS3_hybpiper_summary.csv" ## This table contains sample IDs and species names

## Read table containing species names
sdata <- read.csv(nmpath)

########################################################
## Function to format quartet information from ASTRAL
## tree annotation the normalized quartet score QS
## represents the proportion of input gene tree quartet
## trees satisfied by the species tree. The higher this number,
## the less discordant your gene trees are. Modified from: 
## https://raw.githubusercontent.com/sidonieB/scripts/master/plot_Astral_trees_v2.R
########################################################
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

########################################################
## Extract info about clades
## Take care, the function MRCA does not accept (-)
## in the tip labels
########################################################
clades <- list(blastemanthus = c("B_sp_Amaralsn", "B_grandiflorus_Ducke23783"),
               cespclade = c("F_felicis_Haba445", "F_felicis_Haba18", 
                             "R_regalis_Reinales159"),
               poePLUSwala = c("Poe_pumila_Maas5732", "W_insignis_Kawasaki243"),
               asianclade = c("Indov_calophylla_Brooke9671", 
                              "Sch_henningsii_Morawetz13-2287",
                              "Schll_angustifolia_Bramley100471"),
               cortyleria = c("T_silvana_Liesner22006", "Ad_bicarpellata_Huber10214"),
               sauvagesiaA = c("S_glandulosa_rb_Cardoso3942", "S_alpestris_Mota993",
                               "S_insignis_Cardoso2656"),
               sauvagesiaB = c("S_sprengelii_Oliveira085", "S_longifolia_Harley56621",
                               "S_ramosa_Cardenas44000"))

########################################################
## Iterate tree files, read trees from ASTRAL
## and export tables and plots to pdf
########################################################
## Save Q and LPP values for each tree
tree_summ_q <- data.frame()
tree_summ_p <- data.frame()

## Save info about clades
cladesumm <- data.frame()
plotclades <- TRUE
changeNames <- FALSE

## Iterate species trees
for(i in 1:length(tfiles)){
  
  ## Tree name
  tname <- gsub(pattern = tpattern, replacement = "", tfiles[i])
  tname <- gsub(pattern = tpattern1, replacement = "", tname)
  tname <- gsub(pattern = "[0-9]", replacement = "", tname)
  
  ## Plot names and title
  ptitle <- paste(letters[i], ") ", tname, sep = "")
  
  ## Read tree
  mtree <- read.astral(paste(tpath, "/", tfiles[i], sep = ""))
  is.rooted(mtree)
  
  ## Extract quartet and local probability values using QS_calc and  LP_calc functions
  Q_t_all <- QS_calc(AS_A = mtree)
  P_t_all <- LP_calc(LP_A = mtree)
  
  ## Make a label with the three QS separated by underscores
  for (qi in 1:length(Q_t_all$node)){
    Q_t_all$QS[qi] <- paste(round(Q_t_all$Q1[qi]), "_", round(Q_t_all$Q2[qi]), "_",
                           round(Q_t_all$Q3[qi]), sep = "")
  }
  
  for (pi in 1:length(P_t_all$node)){
    P_t_all$LP[pi] <- paste(round(P_t_all$P1[pi], digits = 1), "_", round(P_t_all$P2[pi], digits = 1), "_",
                           round(P_t_all$P3[pi], digits = 1), sep = "")
  }
  
  ## Save summary
  Q_t_all$tree <- tname
  P_t_all$tree <- tname
  tree_summ_q <- rbind(tree_summ_q, Q_t_all)
  tree_summ_p <- rbind(tree_summ_p, P_t_all)
  
  ## Extract Q and LPP values by clades
  if(plotclades){
    ## Iterate clades and extract Q and LPP information for each clade
    for(j in 1:length(clades)){
      ## Clade members
      clnms <- clades[j][[1]]
      
      ## First check tip labels
      if(grepl("-", paste0(clnms, collapse = ", "))){
        clnms <- gsub("-", "", clnms)
        mtree@phylo$tip.label <- gsub("-", "", mtree@phylo$tip.label)
      }
      
      ## Extract MCRC info for some important clades
      cnode <- MRCA(as_tibble(mtree), clnms)$node
      coffspring <- offspring(as_tibble(mtree), cnode, tiponly = F, self_include = T)
      ctips <- coffspring$label[-which(is.na(coffspring$label))]
      dnodes <- coffspring$node[which(is.na(coffspring$label))]
      
      ## Test
      if(!cnode %in% dnodes){
        cat("Test nodes failed for clade", names(clades[j]),
            "in tree", tname, "\n")
      }
      
      ## Extract qscores
      q1 <- mtree@data[which(mtree@data$node %in% dnodes), c("node", "q1")]
      qscores <- merge(q1, mtree@data[which(mtree@data$node %in% dnodes),
                                      c("node", "q2")], by = "node")
      qscores <- merge(qscores, mtree@data[which(mtree@data$node %in% dnodes),
                                           c("node", "q3")], by = "node")
      
      ## Extract local probabilities
      p1 <- mtree@data[which(mtree@data$node %in% dnodes), c("node", "pp1")]
      localp <- merge(p1, mtree@data[which(mtree@data$node %in% dnodes),
                                     c("node", "pp2")], by = "node")
      localp <- merge(localp, mtree@data[which(mtree@data$node %in% dnodes),
                                         c("node", "pp3")], by = "node")
      
      ## Save summary for clades
      ct <- merge(qscores, localp, by = "node")
      ct$tree <- tname
      ct$clade <- names(clades[j])
      ct$tips <- paste0(sort(ctips), collapse = ";")
      cladesumm <- rbind(cladesumm, ct)
      
      cat("Processed: clade", names(clades)[j], "for tree", tname, "\n")
    }  
  }
  
  ## Change tip names
  if(changeNames){
    for(k in  1:length(mtree@phylo$tip.label)){
      
      if(grepl("-", mtree@phylo$tip.label[k])){ 
        mtree@phylo$tip.label[k] <- gsub("-","",mtree@phylo$tip.label[k]) 
      }
      
      pst <- which(sdata$abbrev == mtree@phylo$tip.label[k])
      
      if(length(pst) >= 1){
        test <- length(unique(sdata$spp_name_voucher[pst]))
        if(length(test) == 1){
          newnm <- unique(sdata$spp_name_voucher[pst])
          if(grepl("-", newnm)) { newnm <- gsub("-", "\\\u00ad", newnm) }
          if(grepl("_sprl", newnm)){ newnm <- gsub("_sprl", "*", newnm) }
          mtree@phylo$tip.label[k] <- newnm
        }
      }else{
        cat("Name", mtree@phylo$tip.label[k], "did not find in tree", tname, "\n")
      }
    }  
  }
  
  ## Plot a basic tree
  t_all <- ggtree(mtree@phylo, ladderize = T, branch.length = "none") +
            geom_tiplab(size = 2, hjust=-.2) +
            xlim_tree(70) +
            ggtitle(paste("Basic tree:", ptitle))
  
  ## Plot tree with node labels
  t_nl <- ggtree(mtree@phylo, ladderize = T, branch.length = "none") +
           geom_text2(aes(label = node), size = 2, hjust=-.3) +
           geom_tiplab(size = 2, hjust=-.2) +
           xlim_tree(70) +
           ggtitle(paste("Node labels:", ptitle))
  
  ## Make pie charts (takes some time if tree is big)
  pies_t_all <- nodepie(Q_t_all, cols = 1:3, 
                        color = c(Q1='darkblue', Q2='green3', Q3='grey'),
                        alpha = 0.8)

  ## Insert pies to the basic tree
  t_pies <- inset(t_all, pies_t_all, width = 0.08, height = 0.08) # 0.42
  
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
  t_complete <- inset(t_textp, pies_t_all, width = 0.08, height = 0.08) +
                 ggtitle(paste0("          \n          \n          \n          ",
                                ptitle))
    
  ## Plot one tree per page
  #plots <- list(t_nl, t_textq, t_textp, t_complete)
  #ggsave(filename = paste(tpath, "/", tname, ".pdf", sep = ""), 
  #       plot = marrangeGrob(plots, nrow = 1, ncol = 1),
  #       width = 15, height = 20)
  
  ## Alternative, arrange all plots per dataset in the same page
  ## open an image window (e.g. pdf) right before the loop
  #grid.arrange(t_textq, t_complete) # All plots for each tree in the same page
  
  ## Plot just the last tree
  ggsave(filename = paste(tpath, "/", "FigureS", i, "_speciesTrees_", tname, ".pdf", sep = ""),
         plot = t_complete,
         width = 12, height = 20)
  
  cat("####### Tree:", tname, " plotted #######", "\n")
}

## Export summary tables Q and LPP
write.csv(tree_summ_q, 
          file = paste(tpath, "/", "qscores.csv", sep = ""),
          row.names = F)
write.csv(tree_summ_p,
          file = paste(tpath, "/", "local_prob.csv", sep = ""),
          row.names = F)

## Merge tables
## Merge Q and LPP tables
table(tree_summ_p$node == tree_summ_q$node) ## Validation test
cladesumm_all <- merge(tree_summ_q, tree_summ_p, by = c("tree","node"))
cladesumm_all <- merge(cladesumm_all, cladesumm[,c("node", "tree", "clade", "tips")],
                       by = c("tree","node"), all = TRUE)

# Adjust tree names
cladesumm_all$tree[which(cladesumm_all$tree == "rmout_infprlgs")] <- "0rmout_infprlgs"
cladesumm_all$tree[which(cladesumm_all$tree == "all_exons")] <- "1all_exons"
cladesumm_all$tree[which(cladesumm_all$tree == "rm_out")] <- "2rm_out"
cladesumm_all$tree[which(cladesumm_all$tree == "str_mssdata")] <- "3str_mssdata"
cladesumm_all$tree[which(cladesumm_all$tree == "str_mssdata_prlgs")] <- "4str_mssdata_prlgs"
cladesumm_all$tree[which(cladesumm_all$tree == "putative_orth")] <- "5putative_orthlgs"
cladesumm_all$tree[which(cladesumm_all$tree == "inform_paralogs")] <- "6inform_prlgs"
cladesumm_all$tree[which(cladesumm_all$tree == "all_paralogs")] <- "7all_prlgs"
cladesumm_all$tree[which(cladesumm_all$tree == "prlgsTen")] <- "8prlgs10"

## Adjust colnames
colnames(cladesumm_all) <- c("tree","node","q1","q2","q3","QS",
                             "pp1","pp2","pp3","LP","clade","tips")

## Remove dataset paralogs10
cladesumm_all <- subset(cladesumm_all, subset = !(cladesumm_all$tree == "8prlgs10"))

## Write summary table
write.csv(cladesumm_all,
          file = paste(tpath, "/", "clades_sum.csv", sep = ""), 
          row.names = F)

###################################################
## Plot recovery data on a species tree
## in this case, species tree from dataset
## rmout_infprlgs
###################################################
## Data path
dpath <- "supplementary_tables/TableS3_hybpiper_summary.csv" 
adlh <- "supplementary_tables/TableS4_hybphaser_summ_table_species.csv"
tname <- "rooted_sauv_0rmout_infprlgs_BS33_sptree_annot.tre"

###################################
## Read data
###################################
## Gene recovery and sequence length
sdata <- read.csv(dpath, header = T, as.is = T)

## AD and LH 
adlhdata <- read.csv(adlh, header = T, as.is = T)

## Merge tables
for(i in unique(sdata$sampleID)){
  pst <- grep(pattern = i, strsplit(x = adlhdata$sampleID, split = "-"))
  sdata[which(sdata$sampleID == i), "AD"] <- adlhdata$allele_divergence_bothprobes_rm_mssdata[pst]
  sdata[which(sdata$sampleID == i), "LH"] <- adlhdata$locus_heterozygosity_bothprobes_rm_mssdata[pst]
}

# Fix hyphen in species names
for(i in 1:length(sdata$sppNameVoucher)){
  if(grepl("-", sdata$sppNameVoucher[i])){
    sdata$sppNameVoucher[i] <- gsub("-", "\\\u00ad", sdata$sppNameVoucher[i])
  }
  
  if(grepl("_sprl", sdata$sppNameVoucher[i])){
    sdata$sppNameVoucher[i] <- gsub("_sprl", "*", sdata$sppNameVoucher[i])
  }
}

## Extract info
# gene recovered
df1 <- data.frame(GenesWithSeqs = tapply(X = sdata$genesWseqs,
                                         INDEX = sdata$sppNameVoucher, FUN = sum))

# total sequence length
df2 <- data.frame(RecoveryLenght = tapply(X = sdata$totalCapturedbp, 
                                          INDEX = sdata$sppNameVoucher, FUN = sum))

# percentage of paralogs
sdata$TotalParalogWarnings <- sdata$ParalogWarningsLong + sdata$ParalogWarningsDepth
df3 <- tapply(X = sdata$TotalParalogWarnings, INDEX = sdata$sppNameVoucher, FUN = sum)
df3 <- data.frame(PctParalogs = round((df3/df1$GenesWithSeqs)*100, digits = 1))

# Allele divergence
df4 <- data.frame(AD = tapply(X = sdata$AD, 
                              INDEX = sdata$sppNameVoucher, FUN = mean))

# Locus heterozygosity
df5 <- data.frame(LH = tapply(X = sdata$LH, 
                              INDEX = sdata$sppNameVoucher, FUN = mean))

# probe
df6 <- data.frame(probe = tapply(X = sdata$dataset, INDEX = sdata$sppNameVoucher, 
                                 FUN = function(x) paste0(x, collapse = "-")))
df6 <- data.frame(species = row.names(df6), probe = df6$probe)
df6$probe <- gsub(pattern = "new_", "", df6$probe)
df6$probe <- gsub(pattern = "old_", "", df6$probe)
df6$probe[which(df6$probe == "a353")] <- "A353"
df6$probe[which(df6$probe == "ochn")] <- "OCHN"
df6$probe[which(df6$probe == "a353-ochn")] <- "BOTH"

###################################
## Read tree
###################################
mtree <- read.astral(paste(tpath, tname, sep = "/"))
is.rooted(mtree)

## Change tip names
for(j in  1:length(mtree@phylo$tip.label)){
  pst <- which(sdata$sppNameAbbrev == mtree@phylo$tip.label[j])
  
  if(length(pst) >= 1){
    test <- length(unique(sdata$sppNameVoucher[pst]))
    if(length(test) == 1){
      newnm <- unique(sdata$sppNameVoucher[pst])
      if(grepl("-", newnm)) { newnm <- gsub("-", "\\\u00ad", newnm) } # change en-dash by hyphen
      if(grepl("_sprl", newnm)){ newnm <- gsub("_sprl", "*", newnm) }
      mtree@phylo$tip.label[j] <- newnm 
    }
  }else{
    cat("Name", mtree@phylo$tip.label[j], "did not find", "\n")
  }
}

## Remove duplicated samples (just those that were retrieved as a clade)
tipsTOrm <- c("Sauvagesia_aliciae_Sastre2339", "Sauvagesia_rubra_Cardoso3942",
              "Sauvagesia_glandulosa_Cardoso2741", "Sauvagesia_paganuccii_Cardoso3578",
              "Godoya_obovata_Neill14978", "Krukoviella_disticha_Neill15849",
              "Fleurydora_felicis_Haba18")

## Remove rouge taxa
tipsTOrm <- c(tipsTOrm, "Sauvagesia_nudicaulis_Maguire36283",
              "Sauvagesia_roraimensis_Wehlburg5647", "Sauvagesia_amoena_CampbellP22564",
              "Tyleria_aristata_Liesner16010", "Sauvagesia_guianensis_araracuarensis_Cardiel4660",
              "Tyleria_aristata_Maguire42313", "Sauvagesia_paniculata_Cardoso2673",
              "Tyleria_tremuloidea_Croat59545", "Tyleria_grandiflora_Pedro.J10394",
              "Tyleria_pendula_Maguire42286", "Tyleria_spectabilis_Nee31125", 
              "Tyleria_apiculata_Liesner24854", "Tyleria_floribunda_Steyermark124565",
              "Philacra_auriculata_Pipoly6867")

mtree <- drop.tip(mtree, tipsTOrm)

###################################
## Extract quartet and local probability
## values using QS_calc and  LP_calc functions
###################################
Q_t_all <- QS_calc(AS_A = mtree)
P_t_all <- LP_calc(LP_A = mtree)

## Make a label with the three QS separated by underscores
for (i in 1:length(Q_t_all$node)){
  Q_t_all$QS[i] <- paste(round(Q_t_all$Q1[i]), "_", round(Q_t_all$Q2[i]), "_",
                         round(Q_t_all$Q3[i]), sep = "")
}

for (j in 1:length(P_t_all$node)){
  P_t_all$LP[j] <- paste(round(P_t_all$P1[j], digits = 1), "_", round(P_t_all$P2[j], digits = 1), "_",
                         round(P_t_all$P3[j], digits = 1), sep = "")
}

###################################
## Plot trees 
###################################
## A basic tree
t_all <- ggtree(mtree@phylo, ladderize = T, branch.length = "none") +
          geom_tiplab(size = 2.8, family = "serif", hjust = -0.04) +
          xlim_tree(70)

## Make pie charts (takes some time if tree is big)
pies_t_all <- nodepie(Q_t_all, cols = 1:3, 
                      color = c(Q1='darkblue', Q2='green3', Q3='grey'),
                      alpha = 0.8)

## Insert pies to the basic tree
t_pies <- inset(t_all, pies_t_all, width = 0.096, height = 0.096)

## Probe annotation on a basic tree
t_complete2 <- t_pies %<+% df6 +
  geom_tippoint(aes(shape = probe)) + 
  theme(legend.position = "right",
        legend.text = element_text(family = "serif"))

## Recovery, AD and LH annotation on t_complete using heatmaps
# Number of genes
colnames(df1) <- "NG"
t_hmap_g <- gheatmap(t_complete2, df1, offset = 15, width = 0.12, 
                     colnames = T, colnames_position = "top", 
                     font.size = 3.5, colnames_offset_y = 0.2) +
              scale_fill_viridis_c(option = "D", name = "Num. Genes\nwith seqs",
                                   direction = -1)

t_hmap_l <- t_hmap_g + new_scale_fill()

# Total sequence length
colnames(df2) <- "TSL"
t_hmap_l <- gheatmap(t_hmap_l, df2, offset = 18.5, width = 0.12,
                     colnames = T, colnames_position = "top",
                     font.size = 3.5, colnames_offset_y = 0.2) +
             scale_fill_viridis_c(option = "D", name = "Total seq\nlenght (bp)",
                                  direction = -1)

t_hmap_p <- t_hmap_l + new_scale_fill()

# Percentage of paralogs
colnames(df3) <- "%PR"
t_hmap_p <- gheatmap(t_hmap_p, df3, offset = 22, width = 0.12,
                     colnames = T, colnames_position = "top", 
                     font.size = 3.5, colnames_offset_y = 0.2) +
             scale_fill_viridis_c(option = "D", name = "%Paralogs",
                                  direction = -1)

t_hmap_ad <- t_hmap_p + new_scale_fill()

# Allele divergence AD
t_hmap_ad <- gheatmap(t_hmap_ad, df4, offset = 25.5, width = 0.12,
                      colnames = T, colnames_position = "top",
                      font.size = 3.5, colnames_offset_y = 0.2) +
               scale_fill_viridis_c(option = "D", name = "AD", 
                                    direction = -1)

t_hmap_lh <- t_hmap_ad + new_scale_fill()

# Locus heterozygosity LH
t_hmap_lh <- gheatmap(t_hmap_lh, df5, offset = 29, width = 0.12,
                      colnames = T, colnames_position = "top",
                      font.size = 3.5, colnames_offset_y = 0.2) +
              scale_fill_viridis_c(option = "D", name = "LH",
                                   direction = -1)
## Clade annotation
t_clades <- t_hmap_lh +
  geom_strip("Sauvagesia_sprengelii_Santarem_Cardoso4081",
             "Sauvagesia_guianensis_Liesner24969",
             barsize = 2, color = "black", fontsize = 5, 
             label = "Sauvagesia\ncladeB", family = "serif",
             offset = 35.5, offset.text = 0.6) +
  geom_strip("Adenarake_muriculata_Maguire60447",
             "Sauvagesia_africana_Nawa119",
             barsize = 2, color = "black", fontsize = 5, 
             label = "Adenarake", family = "serif",
             offset = 35.5, offset.text = 0.6) +
  geom_strip("Sauvagesia_rubra_Baptista570",
             "Sauvagesia_longipes_Pipoly11078",
             barsize = 2, color = "black", fontsize = 5, 
             label = "Sauvagesia\ncladeA", family = "serif",
             offset = 35.5, offset.text = 0.6) +
  geom_strip("Indovethia_calophylla_Brooke9671",
             "Schuurmansia_henningsii_Hoogland8954",
             barsize = 2, color = "black", fontsize = 5, 
             label = "SE Asian\nclade", family = "serif",
             offset = 35.5, offset.text = 0.6) +
  geom_strip("Tyleria_bicarpellata_Pipoly7299",
             "Tyleria_linearis_Tate555",
             barsize = 2, color = "black", fontsize = 5, 
             label = "Tyleria", family = "serif",
             offset = 35.5, offset.text = 0.6) +
  geom_strip("Wallacea_insignis_Kawasaki243",
             "Poecilandra_pumila_Maas5732",
             barsize = 2, color = "black", fontsize = 5, 
             label = "Wallacea +\nPoecilandra", family = "serif",
             offset = 35.5, offset.text = 0.6) +
  geom_strip("Rhytidanthera_regalis_Reinales161",
             "Fleurydora_felicis_Haba445",
             barsize = 2, color = "black", fontsize = 5, 
             label = "Cespedesia\nclade +\nFleurydora", family = "serif",
             offset = 35.5, offset.text = 0.6) +
  geom_strip("Blastemanthus_gemmiflorus_Prance15495",
             "Blastemanthus_grandiflorus_Ducke23783",
             barsize = 2, color = "black", fontsize = 5, 
             label = "Blastemanthus", family = "serif",
             offset = 35.5, offset.text = 0.6) +
  geom_strip("Ochna_membranacea_Konomou261",
             "Hypericum_perforatum_OneKP132735",
             barsize = 2, color = "black", fontsize = 5, 
             label = "Outgroup", family = "serif",
             offset = 35.5, offset.text = 0.6)

###################################
## Export plot
###################################
tname <- gsub(pattern = tpattern, replacement = "sptree_", tname)
tname <- gsub(pattern = tpattern1, replacement = ".pdf", tname)
pdf(paste(tpath, tname, sep = "/"), width = 14, height = 20)
 t_clades
dev.off()