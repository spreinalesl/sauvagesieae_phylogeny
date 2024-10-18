###########################################
## This script estimates the
## generalized Robinson–Foulds distances
## between all combinations of two trees
###########################################

## Package
library(ape) # to manipulate trees
library(TreeDist) # to generalized Robinson–Foulds distances
library(reshape2) # to plot matrix
library(ggplot2) # to plot matrix

## List of trees to compare
tpath <- "/home/spreinalesl/Documents/plot_trees"
tfiles <- list.files(pattern = "\\.tre", tpath, full.names = T)
tpattern <- "/sauv_"
tpattern1 <- "_BS33_sptree_annot.tre"

## Iterate tree files
alltrees <- list()
tnames <- vector()

for(i in 1:length(tfiles)){
  
  ## Tree name
  tname <- gsub(pattern = tpath, replacement = "", tfiles[i])
  tname <- gsub(pattern = tpattern, replacement = "", tname)
  tname <- gsub(pattern = tpattern1, replacement = "", tname)
  tnames[i] <- tname
  
  ## Read tree
  ttree <- read.tree(tfiles[i])
  
  ## Unrooting tree (Astral produces species trees rooted in the first sample present in the gene trees)
  ttree <- unroot(ttree)
  
  ## Save tree in a list
  alltrees[[i]] <- ttree
  names(alltrees)[i] <- tname
  
  cat(tname, "is.rooted", is.rooted(ttree), "ntips", length(ttree$tip.label), "\n")
}

## Create all combinations among input trees
comb <- expand.grid(1:length(tfiles), 1:length(tfiles))

## Calculate TreeDistance for each pair of trees
rfmatrix <- matrix(nrow = 7, ncol = 7, dimnames = list(tnames, tnames))

for(i in 1:nrow(comb)){
  ## Select trees
  tre1 <- comb$Var1[i]
  tre2 <- comb$Var2[i]
  
  ## Test
  cat(dimnames(rfmatrix)[[1]][tre1] == names(alltrees)[tre1], 
      dimnames(rfmatrix)[[2]][tre2] == names(alltrees)[tre2], "\n")
  
  ## Save generalized RF distance in the matrix
  rfmatrix[tre1,tre2] <- TreeDistance(alltrees[[tre1]], alltrees[[tre2]])
}

## Some adjustments before plot
rfmatrix[which(upper.tri(rfmatrix))] <- NA
df <- data.frame(dataset = dimnames(rfmatrix)[[1]], rfmatrix)
df <- melt(df, id = "dataset")
df$value <- round(df$value, digits = 3)

## Plot matrix for generalized RF distances
pdf(paste(tpath, "heatmap_gRFdistances.pdf", sep = "/"))
ggplot(df, aes(x = dataset, y = variable, fill = value)) +
  geom_tile()+
  scale_fill_gradient(low ="#FFFF88", high = "#FF0000", na.value = "white",
                      name="Generalized\nRF distances", limit = c(0, 0.2)) +
  coord_fixed() +
  geom_text(aes(dataset, variable, label = value), color = "black", size = 4) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks = element_blank(),
        legend.justification = c(1, 0),
        legend.position = c(0.4, 0.7),
        legend.direction = "horizontal") +
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1, title.position = "top", title.hjust = 0.5))
dev.off()