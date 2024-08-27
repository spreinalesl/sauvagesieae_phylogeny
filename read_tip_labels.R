library(ape)
lfiles <- list.files(pattern = "tre$", "/home/spreinalesl/Documents", full.names = T)

tipslab <- list()
for(i in 1:length(lfiles)){
  name <- gsub(x = lfiles[i], pattern = "/home/spreinalesl/Documents/rooted_", replacement = "")
  name <- gsub(x = name, pattern = "_BS33_sptree_annot.tre", replacement = "")
  
  tipslab[[i]] <- as.data.frame(x = sort(read.tree(lfiles[i])$tip.label), nm = name)
#  write.table(as.data.frame(x = sort(read.tree(lfiles[i])$tip.label), nm = name), row.names = F,
#              paste("/home/spreinalesl/Documents/trees/tables/", name, ".csv", sep = ""))
}

tipslab <- do.call(cbind, tipslab[-2])
write.table(tipslab, row.names = F, "/home/spreinalesl/Documents/trees/tables/all.csv")
