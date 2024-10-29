##################################################
## This script plot summary information from gene
## trees for each dataset using the output table
## of the bash script 'iqtree_log_summ'
## and q1 values from astral trees
##################################################

## Packages
library(ggplot2)
library(gridExtra) ## for multiple plots
library(viridis) ## heatmap colors
library(ggpubr) ## for multiple plots aligned

## Set path
outpath <- "~/Documents/iqtree_summ"

## Read gtrees summary table
data <- read.table("~/Documents/iqtree_summ/merge_iqtree_gtrees_sum.csv",
                   header = T, sep = ",")

## Read q1 scores from astral by clades
cladesumm <- read.csv("~/Dropbox/Sandra/Manuscritos/NGS_phyloOchn/analyses/results/trees/plots/clades_sum.csv",
                      sep = ",", as.is = T)
## Overall q1 per dataset
totalq1 <- aggregate(cladesumm$q1 ~ cladesumm$tree, FUN = function(x) round(median(x, na.rm = T), digits = 4))
colnames(totalq1) <- c("dataset", "median_q1")

## Merge tables
data <- merge(data, totalq1, by.y = "dataset")

######################
## Plot iqtree summ
######################

bntips <- ggplot(data, aes(x = dataset, y = ntips, fill = dataset, alpha = 1.7)) + 
           geom_boxplot(width = 0.3, outlier.size = 1) +
           theme_bw() +
           theme(axis.title.x = element_blank(),
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 legend.position = "none")

bnsites <- ggplot(data, aes(x = dataset, y = log(nsites), fill = dataset, alpha = 1.7)) + 
            geom_boxplot(width = 0.3, outlier.size = 0.5) +
            theme_bw() +
            theme(axis.title.x = element_blank(),
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  legend.position = "none")

prsminf <- ggplot(data, aes(x = dataset, y = log(npars_inf_sites), fill = dataset, alpha = 1.7)) + 
            geom_boxplot(width = 0.3, outlier.size = 0.5) +
            theme_bw() +
            theme(axis.title.x = element_blank(),
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  legend.position = "none")

totbrl_sum <- ggplot(data, aes(x = dataset, y = log(tbrachn_sum), fill = dataset, alpha = 1.7)) + 
               geom_boxplot(width = 0.3, outlier.size = 0.5) +
               theme_bw() +
               labs(x = "", y = "log(Sum branch lenghts)") +
               theme(axis.text.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     legend.position = "none")

totbrl_median <- ggplot(data, aes(x = dataset, y = tbrachn_median, fill = dataset, alpha = 1.7)) + 
                  geom_boxplot(width = 0.3, outlier.size = 0.5) +
                  theme_bw() +
                  labs(x = "", y = "Median branch lenghts") +
                  theme(axis.text.x = element_blank(),
                        axis.ticks.x = element_blank(),
                        legend.position = "none")

totbrl_mean <- ggplot(data, aes(x = dataset, y = tbrachn_mean, fill = dataset, alpha = 1.7)) + 
                geom_boxplot(width = 0.3, outlier.size = 0.5) +
                theme_bw() +
                theme(axis.title.x = element_blank(),
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                legend.position = "none")

totbrl_var <- ggplot(data, aes(x = dataset, y = log(tbrachn_var), fill = dataset, alpha = 1.7)) + 
               geom_boxplot(width = 0.3, outlier.size = 0.5) +
               theme_bw() +
               labs(x = "", y = "log(Branch lenghts variance)") +
               theme(axis.text.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     legend.position = "none")

bootsupport <- ggplot(data, aes(x = dataset, y = btsupp, fill = dataset, alpha = 1.7)) + 
                geom_boxplot(width = 0.3, outlier.size = 0.5) +
                theme_bw() +
                labs(x = "", y = "Median Bootstap support") +
                theme(axis.text.x = element_blank(),
                      axis.ticks.x = element_blank(),
                      legend.position = "none")

## Export boxplots
pdf(paste(outpath, "/", "gtrees_sum_boxplots.pdf", sep = ""), width = 9, height = 7)
ggarrange(bntips, bnsites, prsminf, totbrl_median, totbrl_var, bootsupport,
          ncol = 3, nrow = 2, align = "hv", common.legend = T, legend = "bottom")
dev.off()

##########################
## Heatmap
##########################
ntips <- ggplot(data, aes(x = dataset, y = gene, fill = ntips)) +
          geom_tile(aes(fill = ntips)) +
          labs(x = "", y = "") +
          theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                axis.text.x = element_text(size = 6),
                legend.key.size = unit(0.4, 'cm'),
                legend.text = element_text(size = 4), legend.title = element_text(size = 6)) +
           scale_x_discrete(labels = c("all_exons", "all_paralogs", "inform\nprlgs",
                              "putative\northlgs", "rm_out", "str\nmssdata",  "str_mssdata\nprlgs")) +
           scale_fill_viridis(direction = -1, option = "A", alpha = 0.9, name = "nsamples")

nsites <- ggplot(data, aes(x = dataset, y = gene, fill = nsites)) +
           geom_tile(aes(fill = nsites)) +
           labs(x = "", y = "") +
           theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                 axis.text.x = element_text(size = 6),
                 legend.key.size = unit(0.4, 'cm'),
                 legend.text = element_text(size = 4), legend.title = element_text(size = 6)) +
           scale_x_discrete(labels = c("all_exons", "all_paralogs", "inform\nprlgs",
                                       "putative\northlgs", "rm_out", "str\nmssdata",  "str_mssdata\nprlgs")) +
           scale_fill_viridis(direction = -1, option = "A", alpha = 0.9, transform = "log", name = "log(nsites)")

## Export occupancy
pdf(paste(outpath, "/", "dataset_occupancy_heatmaps.pdf", sep = ""), width = 9)
ggarrange(ntips, nsites, ncol = 2, nrow = 1, align = "hv")
dev.off()

######################
## Plot astral q1 values
######################
## Plot overall Q1
q1 <- ggplot(cladesumm, aes(x = tree, y = q1, fill = tree, alpha = 1.7)) + 
       geom_boxplot(width = 0.3, outlier.size = 0.5) +
       theme_bw() +
       labs(x = "", y = "Q1") +
       theme(axis.text.x = element_blank(),
             axis.ticks.x = element_blank(),
             legend.position = "none")

lpp1 <- ggplot(cladesumm, aes(x = tree, y = pp1, fill = tree, alpha = 1.7)) + 
         geom_boxplot(width = 0.3, outlier.size = 0.5) +
         theme_bw() +
         labs(x = "", y = "LPP1") +
         theme(axis.text.x = element_blank(),
               axis.ticks.x = element_blank(),
               legend.position = "none")

## prop q2/q3 
cladesumm$prop_q2_q3 <- cladesumm$q2/cladesumm$q3

# Remove extreme value for each dataset
cladesumm$prop_q2_q3[which(cladesumm$prop_q2_q3 == max(cladesumm$prop_q2_q3[which(cladesumm$tree == "all_exons")], na.rm = T))] <- NA
cladesumm$prop_q2_q3[which(cladesumm$prop_q2_q3 == max(cladesumm$prop_q2_q3[which(cladesumm$tree =="all_paralogs")], na.rm = T))] <- NA
cladesumm$prop_q2_q3[which(cladesumm$prop_q2_q3 == max(cladesumm$prop_q2_q3[which(cladesumm$tree == "inform_prlgs")], na.rm = T))] <- NA
cladesumm$prop_q2_q3[which(cladesumm$prop_q2_q3 == max(cladesumm$prop_q2_q3[which(cladesumm$tree == "putative_orthlgs")], na.rm = T))] <- NA
cladesumm$prop_q2_q3[which(cladesumm$prop_q2_q3 == max(cladesumm$prop_q2_q3[which(cladesumm$tree == "rm_out")], na.rm = T))] <- NA
cladesumm$prop_q2_q3[which(cladesumm$prop_q2_q3 == max(cladesumm$prop_q2_q3[which(cladesumm$tree == "str_mssdata")], na.rm = T))] <- NA
cladesumm$prop_q2_q3[which(cladesumm$prop_q2_q3 == max(cladesumm$prop_q2_q3[which(cladesumm$tree == "str_mssdata_prlgs")], na.rm = T))] <- NA

# plot
qprp <- ggplot(cladesumm, aes(x = tree, y = prop_q2_q3, fill = tree, alpha = 1.7)) + 
         geom_boxplot(width = 0.3, outlier.size = 0.5) +
         theme_bw() +
         labs(x = "", y = "Q2/Q3") +
         geom_hline(yintercept = 1, linetype = "dotted", col = "red") +
         theme(axis.text.x = element_blank(),
               axis.ticks.x = element_blank(),
               legend.position = "none")

## Plots by clades
q1_clade <- ggplot(cladesumm, aes(x = clade, y = q1, fill = tree, alpha = 1.7)) + 
             geom_boxplot(outlier.size = 0.5) +
             theme_bw() +
             labs(x = "", y = "Q1") +
             theme(axis.text.x = element_text(size = 5),
             legend.position = "none")

lpp1_clade <- ggplot(cladesumm, aes(x = clade, y = pp1, fill = tree, alpha = 1.7)) + 
               geom_boxplot(outlier.size = 0.5) +
               theme_bw() +
               labs(x = "", y = "LPP1") +
               theme(axis.text.x = element_text(size = 5),
                     legend.position = "none")

qprp_clade <- ggplot(cladesumm, aes(x = clade, y = prop_q2_q3, fill = tree, alpha = 1.7)) + 
               geom_boxplot(outlier.size = 0.5) +
               theme_bw() +
               labs(x = "", y = "Q2/Q3") +
               geom_hline(yintercept = 1, linetype = "dotted", col = "red") +
               theme(axis.text.x = element_text(size = 5),
                     legend.position = "none")

##########################
## Extract legend
##########################
provplot <- ggplot(data, aes(x = dataset, y = ntips, fill = dataset)) + 
  geom_boxplot(width = 0.3) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "bottom")

## Function from: http://www.sthda.com/english/wiki/wiki.php?id_contents=7930
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

extlegend <- get_legend(provplot)

##########################
## Plot arrange
##########################
pdf("~/Documents/iqtree_summ/astraltrees_q1_lpp_summ.pdf")
grid.arrange(grobs = list(q1, q1_clade, lpp1, lpp1_clade, qprp, qprp_clade, extlegend),
             heights = c(3,3,3,1), widths = c(3,5),
             layout_matrix = rbind(c(1,2), c(3,4), c(5,6), c(7,7)))
dev.off()

##########################
## Table number of genes per dataset
##########################
ngenes <- tapply(data$gene, as.factor(data$dataset), FUN = length)  