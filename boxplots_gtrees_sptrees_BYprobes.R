##################################################
## This script plot summary information from gene
## trees for the dataset rmout_infprlgs analysed
## by probe set
##################################################

## Packages
library(ggplot2)
library(gridExtra) ## for multiple plots
library(viridis) ## heatmap colors
library(ggpubr) ## for multiple plots aligned
library(gtable) ## for complex layout
library(grid) ## for complex layout
library(cowplot) ## scatterplot with densities

## Set working directory
setwd("")

## Set path
outpath <- "plot_trees"

## Read summary tables for species trees
## Use the output of the script plot_sptrees.R
tree_summ_q <- read.csv("plot_trees/qscores.csv")
tree_summ_p <- read.csv("plot_trees/local_prob.csv")

## Merge tables
table(tree_summ_p$node == tree_summ_q$node) ## Validation test
tree_summ <- merge(tree_summ_q, tree_summ_p, by = c("tree","node"))

## Subset table to compare just rmout_infprlgs complete + rmout_infprlgsA353 and rmout_infprlgsOCHN
tree_summ <- tree_summ[which(tree_summ$tree == "rmout_infprlgs" | 
                                     tree_summ$tree == "a" |
                                     tree_summ$tree == "schn"),]

tree_summ$tree[which(tree_summ$tree == "rmout_infprlgs")] <- "0rmout_infprlgs"
tree_summ$tree[which(tree_summ$tree == "a")] <- "A353"
tree_summ$tree[which(tree_summ$tree == "schn")] <- "OCHN"

## Concordance factors
cfs <- read.csv("plot_trees/summ_CFs_datasets.csv")
cfs <- cfs[which(cfs$tree == "0rmout_infprlgs" | 
                 cfs$tree == "A353" |
                 cfs$tree == "OCHN"),]

## Read gtrees summary table
data <- read.table("plot_trees/merge_iqtree_gtrees_sum.csv",
                   header = T, sep = ",")
data <- subset(data, subset = (!(data$dataset == "0rmout_infprlgs_rmlociGTE")))
data$probe[grepl(pattern = "g", data$gene)] <- "OCHN"
data$probe[which(is.na(data$probe))] <- "A353"
data2 <- data[which(data$dataset == "0rmout_infprlgs"),]
data3 <- data2[which(data2$probe == "A353"),]
data3$dataset <- "A353"
data4 <- data2[which(data2$probe == "OCHN"),]
data4$dataset <- "OCHN"
data5 <- rbind(data2,data3,data4)

###################################################
## Plots
###################################################
## nsites vs bootstrap support (all datasets by probes)
ccols <- c("#E69F00", "#404788FF")

pdf(paste0(outpath, "glengthBS_alldatasets.pdf"))
nsbs <- ggplot(data, aes(x = log(nsites), y = btsupp, color = probe)) +
         geom_point(size = 0.5) +
         geom_smooth(method = lm, se = FALSE, fullrange = TRUE) +
         scale_color_manual(values = ccols) +
         scale_fill_manual(values = ccols) +
         theme_bw() +
         labs(x = "gene length - log(nsites)", y = "Mean bootstrap support")
nsbs + facet_wrap( ~ dataset)
dev.off()

## nsites vs bootstrap support (just rmout_infprlgs)
data2$nsites <- log(data2$nsites)
nsbs2 <- ggplot(data2, aes(x = nsites, y = btsupp, color = probe)) +
          geom_point(size = 0.5) +
          geom_smooth(method = lm, se = FALSE, fullrange = TRUE) +
          scale_color_manual(values = ccols) +
          scale_fill_manual(values = ccols) +
          theme_bw() +
          labs(x = "gene length - log(nsites)",
               y = "Mean bootstrap support") +
          theme(legend.position = "inside",
                legend.position.inside = c(0.8, 0.15),
                legend.title = element_blank())

# densities
xplot <- ggdensity(data2, "nsites", fill = "probe") +
          scale_color_manual(values = ccols) +
          scale_fill_manual(values = ccols) +
          theme_transparent() +
          theme(legend.position = "none",
                axis.title = element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank())

yplot <- ggdensity(data2, "btsupp", fill = "probe") +
          scale_color_manual(values = ccols) +
          scale_fill_manual(values = ccols) +
          theme_transparent() +
          theme(legend.position = "none",
                axis.title = element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank()) +
          rotate()

# Arranging the plot using cowplot
pdf(paste0(outpath, "glengthBS_rmout_infprlgs.pdf"))
plot_grid(xplot, NULL, nsbs2, yplot, ncol = 2, align = "hv", 
          rel_widths = c(1.5, 0.8), rel_heights = c(0.8, 1.5))
dev.off()

## Gene trees
bntips <- ggplot(data5, aes(x = dataset, y = ntips, fill = dataset, alpha = 1.7)) + 
           geom_boxplot(width = 0.3, outlier.size = 1) +
           theme_bw() +
           theme(axis.title.x = element_blank(),
                 axis.text.x = element_blank(),
                 legend.position = "none")

###########################
## Species trees
###########################
## Compare Q1 and LLP1
q1 <- ggplot(tree_summ, aes(x = tree, y = Q1, fill = tree, alpha = 1.7)) + 
       geom_boxplot(width = 0.3, outlier.size = 0.3) +
       theme_bw() +
       labs(x = "", y = "Q1") +
       theme(axis.title.x = element_text(size = 6),
             legend.position = "none")

lpp1 <- ggplot(tree_summ, aes(x = tree, y = P1, fill = tree, alpha = 1.7)) + 
         geom_boxplot(width = 0.3, outlier.size = 0.3) +
         theme_bw() +
         labs(x = "", y = "LPP1") +
         theme(axis.title.x = element_blank(),
               axis.text.x = element_blank(),
               legend.position = "none")

## gCF
gCF <- ggplot(cfs, aes(x = tree, y = gCF, fill = tree, alpha = 1.7)) + 
        geom_boxplot(width = 0.3, outlier.size = 0.3) +
        theme_bw() +
        labs(x = "", y = "gCF") +
        theme(axis.title.x = element_text(size = 6),
              legend.position = "none")

## Arrage plots
pdf(paste0(outpath, "boxplots_rmout_infprlgsA353-OCHN.pdf"))
ggarrange(bntips, lpp1, q1, gCF,
          ncol = 2, nrow = 2, align = "hv",
          legend = "none")
dev.off()
