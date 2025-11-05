#####################################################################################################################
#### Script to plot phylogenetic trees with all paralogous sequences retrieved from HybPiper 2.0 
#### The script will plot the trees for ALL loci (not only loci with paralogs) in a single PDF called paralog_trees.pdf
#### samples with no paralog warning are plotted normally
#### samples with paralog warnings are flagged as follow:
#### The "main" sequence is flagged in orange and the other alternative sequences (.0, .1, .2, etc.) in blue
#### This script assumes that:
######### - the trees files are in the Newick format and have a "paralogs.tre" suffix
######### - the trees were reconstructed with FastTree (see https://github.com/mossmatters/HybPiper/wiki/Paralogs for more details);
###  other program would work too, but you would need to change the 'pattern' arguments in list.files() and str_remove()
#### Modified from: @Leo-Paul Dagallier - 1st Nov. 2022 available on
### https://github.com/LPDagallier/Phylogenomics_Workshop/blob/master/PHYLOGENY_RECONSTRUCTION/R/plot_hybpiper_paralog_trees.R
#####################################################################################################################

library(ape)
library(stringr)
library(ggplot2)
library(ggtree)
library(gridExtra)

setwd("")

##  Paths and names
probe <- "a353"
fpath <- paste("~/projects/hybpiper/paralogs/", probe, "/trees", sep = "")

## Read all tree files
files <- list.files(path = fpath, pattern = "\\.tre$", full.names = T, recursive=FALSE)
filenames <- list.files(path = fpath, pattern = "\\.tre$", full.names = F, recursive=FALSE)

## Open pdf screen
pdf(paste(fpath, "/paralog_trees.pdf", sep = ""), onefile = TRUE, height = 7, width = 14)
warning_list1 <- c()
warning_list2 <- c()
pb <- txtProgressBar(min = 0, max = length(files), style = 3, width = 100, char = "=") 

## Iterate tree files
for(i in 1:length(files)){
  # print(i)
  t<-read.tree(files[i])
  t_trimmed <- drop.tip(t, tip = t$tip.label[str_which(t$tip.label, pattern = "\\.", negate = T)])
  nsamp <- length(t$tip.label)
  nsamp_trimmed <- length(t_trimmed$tip.label)
  palette <- c("#DD8D29", "#46ACC8", "")
  
  if (nsamp_trimmed > 2) {
    # whole trees:
    lab_size = max(min(c(5,70/nsamp)), 0.2)
    line_size = min(c(0.5, 30/nsamp))
    title = str_remove(string = filenames[i], pattern = "_paralogs.*")
    color = t$tip.label
    color[str_which(string = t$tip.label,pattern = ".main")] <- palette[1]
    color[str_which(string = t$tip.label,pattern = "\\.[:digit:]")] <- palette[2]
    color[which(color != "#46ACC8" & color != "#DD8D29")] <- NA
    p1 <- ggtree(t, layout="equal_angle", size = line_size) + geom_tippoint(color = color) + labs(title = title, size = 0.5)
    p2 <- ggtree(t, size = line_size) + geom_tippoint(color = color, alpha = c(0.7), size = min(5, 100/nsamp)) + geom_tiplab(size = lab_size, color = "black") + geom_treescale(color = "black") + ggexpand(5/nsamp, side = "h")
    
    # trimmed trees:
    lab_size = max(min(c(5,70/nsamp_trimmed)), 0.2)
    line_size = min(c(0.5, 30/nsamp_trimmed))
    title = paste0(str_remove(string = filenames[i], pattern = "_paralogs.*"), "_trimmed")
    color = t_trimmed$tip.label
    color[str_which(string = t_trimmed$tip.label,pattern = ".main")] <- palette[1]
    color[str_which(string = t_trimmed$tip.label,pattern = "\\.[:digit:]")] <- palette[2]
    color[which(color != "#46ACC8" & color != "#DD8D29")] <- NA
    p3 <- ggtree(t_trimmed, layout="equal_angle", size = line_size) + geom_tippoint(color = color) + labs(title = title, size = 0.5)
    p4 <- ggtree(t_trimmed, size = line_size) + geom_tippoint(color = color, alpha = c(0.7), size = min(5, 100/nsamp_trimmed)) + geom_tiplab(size = lab_size, color = "black") + geom_treescale(color = "black") + ggexpand(5/nsamp_trimmed, side = "h")
    
  } else {
    if(nsamp_trimmed > 0){
      
      warning_list1 <- rbind(warning_list1, filenames[i])
      
      # whole trees:
      lab_size = max(min(c(5,70/nsamp)), 0.2)
      line_size = min(c(0.5, 30/nsamp))
      title = str_remove(string = filenames[i], pattern = "_paralogs.*")
      color = t$tip.label
      color[str_which(string = t$tip.label,pattern = ".main")] <- palette[1]
      color[str_which(string = t$tip.label,pattern = "\\.[:digit:]")] <- palette[2]
      color[which(color != "#46ACC8" & color != "#DD8D29")] <- NA
      p1 <- ggtree(t, layout="equal_angle", size = line_size) + geom_tippoint(color = color) + labs(title = title, size = 0.5)
      p2 <- ggtree(t, size = line_size) + geom_tippoint(color = color, alpha = c(0.7), size = min(5, 100/nsamp)) + geom_tiplab(size = lab_size, color = "black") + geom_treescale(color = "black") + ggexpand(5/nsamp, side = "h")
      p3 <- ggplot() + theme_void() 
      p4 <- ggplot() + theme_void() 
    } else{
      warning_list2 <- rbind(warning_list2, filenames[i])
    }
  }
  
  # arrange all the plots:
  grid.arrange(p1, p3, p2, p4)
  
  # set progress
  setTxtProgressBar(pb, i)
}
dev.off()
close(pb)
cat("Less than 3 paralogous sequences were found in the following trees (just complete trees were plotted):" , warning_list1,
    "There are no paralogous sequences in the following trees:", warning_list2, sep = "\n")