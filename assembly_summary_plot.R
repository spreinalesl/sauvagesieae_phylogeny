##################################################
## This script summarises assembly results
## using HybPiper in a plot
##################################################

## Packages
library(ggplot2)
library(ggdist)
library(ggpubr)
#library(viridis)

## Set working directory
setwd("")

##########################################
## Capture and assembly
##########################################

## Read table
dpath <- "./supplementary_tables"
filename <- "TableS3_hybpiper_summary.csv"
data_all <- read.csv(paste(dpath, filename, sep = "/"), as.is = T)

## Figures path
fpath <- "./figures"

## Create new columns for probe set
# probeset
data_all$probes <- gsub("new_|old_", "", data_all$dataset)

# PctParalogs
data_all$PctParalogs <- ((data_all$ParalogWarningsDepth + data_all$ParalogWarningsLong)/
                          data_all$genesWseqs)*100

# Proportion of genes with seqs
data_all$pctgenesWseqsFasta <- (data_all$genesWseqs/data_all$GenesMapped)*100

##############################
## boxplot + jitter just for
## new data (Sauvagesieae)
##############################
## Subset for new data
data <- subset(data_all, subset = grepl("new", data_all$dataset))

#ccols <- c("#404788FF", "#29AF7FFF")
ccols <- c("#E69F00", "#404788FF")

## Percentage of reads on target
pcota <- ceiling(mean(data$PctOnTarget[which(data$probes == "a353")]))
pcoto <- ceiling(mean(data$PctOnTarget[which(data$probes == "ochn")]))

pc <- ggplot(data, aes(probes, PctOnTarget, colour = probes)) + 
        geom_boxplot(width = .2, outlier.shape = NA) + 
        geom_jitter(aes(fill = probes), width = .05, alpha = .3) +
        stat_halfeye(aes(fill = probes), adjust = .9, width = .7, .width = 0, justification = -.2,
                     point_colour = NA, alpha = .7) +
        scale_color_manual(values = ccols) +
        scale_fill_manual(values = ccols) +
        theme_bw() +
        theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
              axis.ticks.x = element_blank(), legend.position="none") +
        labs(y = "% Reads on target") #+
        #geom_text(label = "Mean % reads \n on target", colour = "darkgrey", x = 0.9, y = 50, size = 3) +
        #geom_text(label = paste0(pcota, "%", collapse = ""), colour = ccols[1], x = 0.9, y = 45, size = 3) +
        #geom_text(label = paste0(pcoto, "%", collapse = ""), colour = ccols[2], x = 0.9, y = 42, size = 3)

## Genes recovered
# Mapped genes
mpa <- ceiling(mean(data$GenesMapped[which(data$probes == "a353")]))
mpo <- ceiling(mean(data$GenesMapped[which(data$probes == "ochn")]))

# genes with sequences
gwsa <- ceiling(mean(data$genesWseqs[which(data$probes == "a353")]))
gwso <- ceiling(mean(data$genesWseqs[which(data$probes == "ochn")]))

# Table
stable <- matrix(data = c(mpa, gwsa, mpo, gwso), nrow = 2, ncol = 2,
                 dimnames = list(c("GenesMapped", "GenesWSeqs")))
stable <- ggtexttable(stable, cols = NULL,
                      theme = ttheme(rownames.style = rownames_style(color = "darkgrey", size = 7.5, face =  "plain"),
                                     tbody.style = tbody_style(color = sort(c(ccols, ccols), decreasing = T),
                                                               fill = "white", size = 7.5)))

stable <- tab_add_border(stable, from.row = 1, to.row = tab_nrow(stable), from.column = 1,
                           to.column = tab_ncol(stable), linetype = 1, linewidth = 1.5, linecolor = "black")

# Plot
pctgws <- ggplot(data, aes(probes, pctgenesWseqsFasta, colour = probes)) + 
           geom_boxplot(width = .2, outlier.shape = NA) + 
           geom_jitter(aes(fill = probes), width = .05, alpha = .3) +
           stat_halfeye(aes(fill = probes), adjust = .7, width = .7, .width = 0, justification = -.2,
                            point_colour = NA, alpha = .7) +
           scale_color_manual(values = ccols) +
           scale_fill_manual(values = ccols) +
           theme_bw() +
           theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(), legend.position="none")+
           annotation_custom(ggplotGrob(stable), xmin = 2, xmax = 2.1, ymax = 25) +
           labs(y = "% Genes with seqs")

## Number of paralogs
plg <- ggplot(data, aes(probes, PctParalogs, colour = probes)) + 
         geom_boxplot(width = .2, outlier.shape = NA) + 
         geom_jitter(aes(fill = probes), width = .05, alpha = .3) +
         stat_halfeye(aes(fill = probes), adjust = .9, width = .7, .width = 0, justification = -.2,
                      point_colour = NA, alpha = .7) +
         scale_color_manual(values = ccols) +
         scale_fill_manual(values = ccols) +
         theme_bw() +
         theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
               axis.ticks.x = element_blank(), legend.position="none") +
         labs(y = "% Genes with paralog warnings")

## total length of sequences recovered by hybpiper
## Total length for a353 820790
## Total length for ochn 657352
bpa353 <- round(mean(data$totalCapturedbp[which(data$probes == "a353")]), digits = 1)
bpochn <- round(mean(data$totalCapturedbp[which(data$probes == "ochn")]), digits = 1)

stable2 <- matrix(data = c("Recovery\nlength (bp)", bpa353, bpochn), nrow = 3, ncol = 1)

stable2 <- ggtexttable(stable2, cols = NULL,
                      theme = ttheme(tbody.style = tbody_style(color = c("darkgrey", ccols),
                                                               fill = "white", size = 7.5)))

stable2 <- tab_add_border(stable2, from.row = 1, to.row = tab_nrow(stable2), from.column = 1,
                         to.column = tab_ncol(stable2), linetype = 1, linewidth = 1.5, linecolor = "black")

rl <- ggplot(data, aes(probes, percCapturedbp, colour = probes)) + 
        geom_boxplot(width = .2, outlier.shape = NA) + 
        geom_jitter(aes(fill = probes), width = .05, alpha = .3) +
        stat_halfeye(aes(fill = probes), adjust = .9, width = .7, .width = 0, justification = -.2,
                     point_colour = NA, alpha = .7) +
        scale_color_manual(values = ccols) +
        scale_fill_manual(values = ccols) +
        theme_bw() +
        theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
              axis.ticks.x = element_blank(), legend.position="none") +
       labs(y = "% Base pairs recovered") +
       annotation_custom(ggplotGrob(stable2), xmin = 0.7, xmax = 1.1, ymin = 63)
            
##############################
## Open pdf and arrange multiple
## ggplots on a single page
##############################
pdf(paste(fpath, "assemply_summ2.pdf", sep = "/"), width = 8, height = 7.5)
#png(paste0(fpath, "assemply_summ.png", sep = ""))
ggarrange(pc, rl, pctgws, plg,
          labels = c("A", "B", "C", "D"),
          font.label = list(size = 12, color = "black", face = "bold", family = "serif"),
          ncol = 2, nrow = 2, align = "hv",
          common.legend = TRUE, legend = "bottom")

#annotate_figure(allplots, top = text_grob("Capture and assembly summary", color = "black",
#                                          face = "bold", size = 12, family = "serif"))
dev.off()

##################################################################
#### Another plots, New vs old data
##################################################################
## Assemply_summ_oldvsnew3 doesn't include data for non Sauvagesieae
# taxa as Ouratea, Ochna, etc

## Remove some data
# Ochna, Ouratea, Bonnetia, Campylospermum, Elvasia
data <- data_all
data <- data[-grep(pattern = "O_|Ochna|Bonn|Camp|Elv", x = data$sppNameVoucher),]
data$source[grepl("new",data$dataset)] <- "newdata"
data$source[grepl("old_ochn",data$dataset)] <- "Schn2020"
data$source[grepl("old_a353",data$dataset)] <- "Shah2021"

## Number of genes with seqs according the seq_lengths.tsv file (A353 Vs Ochn)
## and (New vs old data)
pc2 <- ggplot(data, aes(source, PctOnTarget, colour = probes)) + 
         geom_boxplot(aes(colour = probes)) +
         geom_point(aes(colour = probes), alpha = .3,
                    position = position_jitterdodge(jitter.width = 0.08)) +
         scale_color_manual(values = ccols) +
         scale_fill_manual(values = ccols) +
         theme_bw() +
         theme(axis.title.x = element_blank(),
               axis.text.x = element_blank(),
               legend.position="none") +
         labs(y = "% Reads on target")

## Number of genes with seqs according the seq_lengths.tsv file (A353 Vs Ochn)
## and (New vs old data)
ngenes <- ggplot(data, aes(source, pctgenesWseqsFasta, colour = probes)) + 
            geom_boxplot(aes(colour = probes)) +
            geom_point(aes(colour = probes), alpha = .3,
                       position = position_jitterdodge(jitter.width = 0.08)) +
            scale_color_manual(values = ccols) +
            scale_fill_manual(values = ccols) +
            theme_bw() +
            theme(axis.title.x = element_blank(),
                  legend.position ="none") +
            labs(y = "% Genes with seqs")

## Total length of sequences recovered by hybpiper (A353 Vs Ochn) and (New vs old data)
rl2 <- ggplot(data, aes(source, percCapturedbp, colour = probes)) + 
         geom_boxplot(aes(colour = probes)) +
         geom_point(aes(colour = probes), alpha = .3,
                    position = position_jitterdodge(jitter.width = 0.08)) +
         scale_color_manual(values = ccols) +
         scale_fill_manual(values = ccols) +
         theme_bw() +
         theme(axis.title.x = element_blank(), 
               axis.text.x = element_blank(),
               legend.position ="none") +
         labs(y = "% Base pairs recovered")

## Paralogs (A353 Vs Ochn) and (New vs old data)
prlg2 <- ggplot(data, aes(source, PctParalogs, colour = probes)) + 
          geom_boxplot(aes(colour = probes)) +
          geom_point(aes(colour = probes), alpha = .3,
                     position = position_jitterdodge(jitter.width = 0.08)) +
          scale_color_manual(values = ccols) +
          scale_fill_manual(values = ccols) +
          theme_bw() +
          theme(axis.title.x = element_blank(),
                legend.position="none") +
          labs(y = "% Genes with paralog warnings")

## Pdf
pdf(paste(fpath, "assemply_summ_oldvsnew_without_Outgroup_ed.pdf", sep = "/"),
    width = 9, height = 6)
#png(paste0(fpath, "assemply_summ_oldvsnew.png", sep = ""))
ggarrange(pc2, rl2, ngenes, prlg2,
         #labels = c("A", "B", "C"),
          font.label = list(size = 12, color = "black", face = "bold", family = "serif"),
          ncol = 2, nrow = 2, align = "hv",
          common.legend = TRUE, legend = "right")

#annotate_figure(anotherplots, top = text_grob("Capture and assembly results: old vs new data", color = "black",
#                                               face = "bold", size = 12, family = "serif"))
dev.off()