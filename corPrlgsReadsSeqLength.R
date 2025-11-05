## Packages
library(ggplot2)
library(ggdist)
library(ggpubr)

##########################################
## Read capture and assembly data
##########################################
dpath <- "supplementary_tables"
filename <- "TableS3_hybpiper_summary.csv"
data <- read.csv(paste(dpath, filename, sep = "/"), as.is = T)

## Probeset
data$probe <- gsub("new_|old_", "", data$dataset)

## PctParalogs
data$Totalparalogs <- data$ParalogWarningsDepth + data$ParalogWarningsLong 
data$PctParalogs <- (data$Totalparalogs/data$genesWseqs)*100

##########################################
## Read gene trees summary table
##########################################
gtpath <- paste(dpath, "TableS5_hybphaserANDhybpiper_summ_table_genes.csv", sep = "/")
gtdata <- read.csv(gtpath, header = T, as.is = T)

##########################################
## ## Plots
##########################################
cols <- c("#E69F00", "#404788FF")

## Total sequence length recovered
seql <- ggplot(data, aes(x = totalCapturedbp, y = PctParalogs, col = probe)) + 
         geom_point(aes(fill = probe), shape = 21, size = 1.5) +
         scale_color_manual(values = cols) +
         scale_fill_manual(values = cols) +
         labs(y = "% Genes with paralog warnings",
              x = "Total sequence length recovered (bp)") +
         theme_bw()

## Number of reads
nreads <- ggplot(data, aes(x = NumReads, y = PctParalogs, col = probe)) + 
           geom_point(aes(fill = probe), shape = 21, size = 1.5) +
           geom_smooth(method = "lm", se = T, fullrange = T, alpha = 0.2) +
           scale_color_manual(values = cols) +
           scale_fill_manual(values = cols) +
           labs(x = "Total number of reads", y = "") +
           theme_bw()

## Gene length
gtdata <- gtdata[-which(gtdata$geneMeanLength == max(gtdata$geneMeanLength)),]
gtdata <- gtdata[-which(gtdata$percentageParalogs == max(gtdata$percentageParalogs, na.rm = T)),]
glength <- ggplot(gtdata, aes(x = geneMeanLength, y = percentageParalogs, col = probe)) +
            geom_point(aes(fill = probe), shape = 21, size = 1.5) +
            geom_smooth(method = "lm", se = T, fullrange = T, alpha = 0.2) +
            scale_color_manual(values = cols) +
            scale_fill_manual(values = cols) +
            theme(axis.title.y = element_blank(), 
                  axis.text.y = element_blank()) +
            labs(x = "Mean sequence length per gene (bp)") +
            theme_bw()
 
## Export plots
pdf(paste(dpath, "corr_prlgs_seqLength_reads.pdf", sep = "/"),
    height = 3.5)
ggarrange(seql, nreads,
          labels = c("a)", "b)"),
          font.label = list(size = 12, color = "black",
                            family = "serif"),
          ncol = 2, nrow = 1, align = "hv",
          common.legend = TRUE, legend = "right")
dev.off()

#########################################################
##################### Alternative #######################
## Use the length of the gene alignment
iqtreed <- read.table(paste(dpath, "merge_iqtree_gtrees_sum.csv", sep = "/"),
                      header = T, sep = ",")
iqtreed <- subset(iqtreed, subset = iqtreed$dataset == "0rmout_infprlgs")

new_gtdata <- data.frame(gene = gtdata$gene,
                         probe = gtdata$probe,
                         geneLength = NA,
                         pctPrlgs = gtdata$percentageParalogs)

for(i in 1:nrow(new_gtdata)){
  pst <- which(iqtreed$gene == new_gtdata$gene[i])
  if(length(pst) > 0){
    new_gtdata$geneLength[i] <- iqtreed$nsites[pst]
  }
}

glength <- ggplot(new_gtdata, aes(x = log(geneLength), y = pctPrlgs, col = probe)) +
  geom_point(aes(fill = probe), shape = 21, size = 2) +
  geom_smooth(method = "lm", se = F, fullrange = T) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_blank()) +
  labs(x = "Gene length (log)") +
  theme_bw()