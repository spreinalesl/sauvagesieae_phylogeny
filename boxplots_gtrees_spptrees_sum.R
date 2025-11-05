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
library(gtable) ## for complex layout
library(grid) ## for complex layout
library(conover.test) ## For a post-hoc Conover-Iman test

## Set working directory
setwd("")

## Set path
outpath <- "plot_trees"

## Read gtrees summary table
data <- read.table(paste(outpath, "merge_iqtree_gtrees_sum.csv", sep = "/"),
                   header = T, sep = ",")
data <- subset(data, subset = (!(data$dataset == "0rmout_infprlgs_rmlociGTE")))

## Read q1 scores from astral by clades
cladesumm <- read.csv(paste(outpath, "clades_sum.csv", sep = "/"),
                      sep = ",")

## Concordance factors
cfs <- read.csv(paste(outpath, "summ_CFs_datasets.csv", sep = "/"))

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
pdf(paste(outpath, "gtrees_sum_boxplots.pdf", sep = "/"), width = 9, height = 7)
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
           scale_x_discrete(labels = c("rmout\ninfprlgs", "all_exons", "rm_out", "str\nmssdata",
                                       "str_mssdata\nprlgs", "putative\northlgs", "inform\nprlgs",
                                       "all_paralogs")) +
           scale_fill_viridis(direction = -1, option = "A", alpha = 0.9, name = "nsamples")

nsites <- ggplot(data, aes(x = dataset, y = gene, fill = nsites)) +
           geom_tile(aes(fill = nsites)) +
           labs(x = "", y = "") +
           theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                 axis.text.x = element_text(size = 6),
                 legend.key.size = unit(0.4, 'cm'),
                 legend.text = element_text(size = 4), legend.title = element_text(size = 6)) +
           scale_x_discrete(labels = c("rmout\ninfprlgs", "all_exons", "rm_out", "str\nmssdata",
                                       "str_mssdata\nprlgs", "putative\northlgs", "inform\nprlgs",
                                       "all_paralogs")) +
           scale_fill_viridis(direction = -1, option = "A", alpha = 0.9, transform = "log", name = "log(nsites)")

## Export occupancy
pdf(paste(outpath, "dataset_occupancy_heatmap.pdf", sep = "/"), width = 9)
ggarrange(ntips, nsites, ncol = 2, nrow = 1, align = "hv")
dev.off()

##########################
## Conover test gene trees
##########################
## Iterate variables
data <- data[,c("dataset", "nsites", "tbrachn_median", "btsupp")]

for(i in colnames(data)[-1]){
  
  ## Remove outliers before, because they affect KW test
  # use boxplot output which detects outliers using quantiles (Q1-1.5*IQR) 
  datasub <- data.frame()
  for(j in unique(data$dataset)){
    bxpl <- boxplot(data[which(data$dataset == j), c(i)], plot = F)
    datasub <- rbind(datasub,
                     subset(data[which(data$dataset == j), c("dataset", i)],
                            subset = !(data[which(data$dataset == j), c(i)] %in% bxpl$out)))
  }
  
  ## Plot test
  #boxplot(datasub[,c(i)] ~ datasub$dataset, main = i)
  
  ## Kuskal Wallis test
  kw <- kruskal.test(datasub[,c(i)] ~ datasub$dataset)
  cat("Kw for", i, kw$p.value, "\n")
  
  ## Export results
  capture.output(file = paste(outpath, "conoverTest.txt", sep = "/"), append = T,
                 cat("\n", "\n", "#################### ", i, " ####################", "\n", "\n"),
                 paste0("n=", length(datasub[,c(i)])), kw)
  
  ## Post-hoc Conover-Iman test for significant different variables
  if(kw$p.value < 0.05){
    capture.output(file = paste(outpath, "conoverTest.txt", sep = "/"), append = T,
                   cat("\n", "\n"),
                   conover.test(datasub[,c(i)], datasub$dataset, method = "bonferroni"),
                   cat("\n", "\n"))
  }
}

######################
## Export gtrees figure
######################
## Align boxplots
g1 <- ggplotGrob(bnsites)
g2 <- ggplotGrob(totbrl_median)
g3 <- ggplotGrob(bootsupport)
g <- rbind(g1, g2, g3, size = "first")
g$widths <- unit.pmax(g1$widths, g2$widths, g3$widths)

## Complex layout
pdf(paste(outpath, "gtrees_summ_stats.pdf", sep = "/"), width = 9, height = 7)
ggarrange(ntips, g, ncol = 2, widths = c(1.2, 0.8))
dev.off()

######################
## Plot astral q1 values
######################
## Plot overall Q1
q1 <- ggplot(cladesumm, aes(x = tree, y = q1, fill = tree, alpha = 1.7)) + 
       geom_boxplot(width = 0.3, outlier.size = 0.3) +
       theme_bw() +
       labs(x = "", y = "Q1") +
       theme(axis.title.x = element_text(size = 6),
             axis.text.x = element_blank(),
             axis.ticks.x = element_blank(),
             legend.position = "none")

lpp1 <- ggplot(cladesumm, aes(x = tree, y = pp1, fill = tree, alpha = 1.7)) + 
         geom_boxplot(width = 0.3, outlier.size = 0.3) +
         theme_bw() +
         labs(x = "", y = "LPP1") +
         theme(axis.title.x = element_text(size = 6),
               axis.text.x = element_blank(),
               axis.ticks.x = element_blank(),
               legend.position = "none")

## prop q2/q3 
cladesumm$prop_q2_q3 <- cladesumm$q2/cladesumm$q3

# plot
qprp <- ggplot(cladesumm, aes(x = tree, y = prop_q2_q3, fill = tree, alpha = 1.7)) + 
         geom_boxplot(width = 0.3, outlier.size = 0.3) +
         theme_bw() +
         labs(x = "", y = "Q2/Q3") +
         geom_hline(yintercept = 1, linetype = "dotted", col = "red") +
         theme(axis.title.x = element_text(size = 6),
               axis.text.x = element_blank(),
               axis.ticks.x = element_blank(),
               legend.position = "none")

## Plots by clades
q1_clade <- ggplot(cladesumm, aes(x = clade, y = q1, fill = tree, alpha = 1.7)) + 
             geom_boxplot(outlier.size = 0.3) +
             theme_bw() +
             labs(x = "", y = "") +
             theme(axis.text.x = element_text(size = 5),
                   legend.position = "none")

lpp1_clade <- ggplot(cladesumm, aes(x = clade, y = pp1, fill = tree, alpha = 1.7)) + 
               geom_boxplot(outlier.size = 0.3) +
               theme_bw() +
               labs(x = "", y = "") +
               theme(axis.text.x = element_text(size = 5),
                     legend.position = "none")

qprp_clade <- ggplot(cladesumm, aes(x = clade, y = prop_q2_q3, fill = tree, alpha = 1.7)) + 
               geom_boxplot(outlier.size = 0.3) +
               theme_bw() +
               labs(x = "", y = "") +
               geom_hline(yintercept = 1, linetype = "dotted", col = "red") +
               theme(axis.text.x = element_text(size = 5),
                     legend.position = "none")

######################
## Plot astral gCF values
######################
gCF <- ggplot(cfs, aes(x = tree, y = gCF, fill = tree, alpha = 1.7)) + 
        geom_boxplot(width = 0.3, outlier.size = 0.3) +
        theme_bw() +
        labs(x = "", y = "gCF") +
        theme(axis.title.x = element_text(size = 6),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              legend.position = "none")

# by clades
gCFs <- subset(cfs, subset = !(is.na(cfs$clade)))
gCF_clade <- ggplot(gCFs, aes(x = clade, y = gCF, fill = tree, alpha = 1.7)) + 
              geom_boxplot(outlier.size = 0.3) +
              theme_bw() +
              labs(x = "", y = "") +
              theme(axis.text.x = element_text(size = 5),
                    legend.position = "none")

#######################
## Export plot
#######################
pdf(paste(outpath, "astraltrees_q1_lpp_gCF_summ.pdf", sep = "/"), height = 10, width = 12)
ggarrange(q1, q1_clade, lpp1, lpp1_clade, gCF, gCF_clade, widths = c(0.8, 1.5),
          ncol = 2, nrow = 3, align = "hv", common.legend = T, legend = "bottom")
dev.off()

##########################
## Table number of genes per dataset
##########################
ngenes <- tapply(data$gene, as.factor(data$dataset), FUN = length)

write.table(ngenes, paste(outpath, "ngenesPERdataset.txt", sep = "/"))

##########################
## Conover test for 
## species trees
##########################
## Iterate variables
#cladesumm <- cfs
colnames(cladesumm)[9] <- "dataset"
variables <- c("q1","pp1")#,"gCF")
data <- cladesumm

for(i in variables){
 
  ## Remove outliers before, because they affect KW test
  # use boxplot output which detects outliers using quantiles (Q1-1.5*IQR) 
  datasub <- data.frame()
  for(j in unique(data$dataset)){
    bxpl <- boxplot(data[which(data$dataset == j), c(i)], plot = F)
    datasub <- rbind(datasub,
                     subset(data[which(data$dataset == j), c("dataset", i)],
                            subset = !(data[which(data$dataset == j), c(i)] %in% bxpl$out)))
  }
  
  ## Plot test
  #boxplot(datasub[,c(i)] ~ datasub$dataset, main = i)
  
  ## Kuskal Wallis test between datasets
  kw <- kruskal.test(datasub[,c(i)] ~ datasub$dataset)
  cat("Kw for", i, kw$p.value, "\n")
  
  ## Export results
  capture.output(file = paste(outpath, "conoverTest_sptrees.txt", sep = "/"), append = T,
                 cat("\n", "\n", "#################### ", i, " ####################", "\n", "\n"),
                 paste0("n=", length(datasub[,c(i)])), kw)
  
  ## Post-hoc Conover-Iman test for significant different variables
  if(kw$p.value < 0.05){
    capture.output(file = paste(outpath, "conoverTest_sptrees.txt", sep = "/"), append = T,
                   cat("\n", "\n"),
                   conover.test(datasub[,c(i)], datasub$dataset, method = "bonferroni"),
                   cat("\n", "\n"))
  }
  
  ###########################
  ## Kruskal Wallis by clade
  ## and dataset
  ###########################
  for(k in unique(data$clade)){
    datasub2 <- subset(data, subset = data$clade == k)
    
    ## rm outliers
    datasub3 <- data.frame()
    for(j in unique(datasub2$dataset)){
      bxpl <- boxplot(datasub2[which(datasub2$dataset == j), c(i)], plot = F)
      datasub3 <- rbind(datasub3,
                       subset(datasub2[which(datasub2$dataset == j), c("dataset", i)],
                              subset = !(datasub2[which(datasub2$dataset == j), c(i)] %in% bxpl$out)))
    }
    
    ## Kruskal Wallis within each clade
    kw <- kruskal.test(datasub3[,c(i)] ~ datasub3$dataset)
    cat("Kw for", i, "in clade", k, kw$p.value, "\n")
    
    ## Plot test
    #boxplot(datasub2[,c(i)] ~ datasub2$dataset, main = i)
    
    ## Export results
    capture.output(file = paste(outpath, "conoverTest_sptreesByClade.txt", sep = "/"), append = T,
                   cat("\n", "\n", "#################### ", i, "clade", k, " ####################", "\n", "\n"),
                   paste0("n=", length(datasub3[,c(i)])), kw)
    
    ## Post-hoc Conover-Iman test for significant different variables
    if(kw$p.value < 0.05){
      capture.output(file = paste(outpath, "conoverTest_sptreesByClade.txt", sep = "/"), append = T,
                     cat("\n", "\n"),
                     conover.test(datasub3[,c(i)], datasub3$dataset, method = "bonferroni"),
                     cat("\n", "\n"))
    }
  }
}

#####################################
## Kuskal Wallis test between clades
#####################################
for(i in variables){
  for(j in unique(data$dataset)){
    ## Kruskal Wallis between clades
    datasub3 <- subset(data, subset = data$dataset == j)
    
    datasub <- data.frame()
    for(k in unique(datasub3$clade)){
      bxpl <- boxplot(datasub3[which(datasub3$clade == k), c(i)], plot = F)
      datasub <- rbind(datasub,
                       subset(datasub3[which(datasub3$clade == k), c("clade", i)],
                              subset = !(datasub3[which(datasub3$clade == k), c(i)] %in% bxpl$out)))
    }
    
    ## Kruskal Wallis
    kw_btw <- kruskal.test(datasub[,c(i)] ~ datasub$clade)
    boxplot(datasub[,c(i)] ~ datasub$clade, main = paste(i,j, sep = " - "))
    cat("Kw for", i, j, kw_btw$p.value, "\n")
    
    ## Export results
    capture.output(file = paste(outpath, "conoverTest_sptrees_betweenClades.txt", sep = "/"), append = T,
                   cat("\n", "\n", "#################### ", i, j, " ####################", "\n", "\n"),
                   paste0("n=", length(datasub[,c(i)])), kw_btw)
    
    ## Post-hoc Conover-Iman test for significant different variables
    if(kw_btw$p.value < 0.05){
      capture.output(file = paste(outpath, "conoverTest_sptrees_betweenClades.txt", sep = "/"), append = T,
                     cat("\n", "\n"),
                     conover.test(datasub[,c(i)], datasub$clade, method = "bonferroni"),
                     cat("\n", "\n"))
    }
  }
}