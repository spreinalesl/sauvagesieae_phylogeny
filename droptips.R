library(ape)

tfiles <- list.files(pattern = "\\.tre", "~/Dropbox/Sandra/Manuscritos/NGS_phyloOchn/analyses/results/phylonetworks/starting_tree", full.names = T)

subsetspp <- read.table("~/Dropbox/Sandra/Manuscritos/NGS_phyloOchn/analyses/results/phylonetworks/sampling/sampling_sppnames.txt",
                      header = T)

smpdata <- read.csv("~/Dropbox/Sandra/Manuscritos/NGS_phyloOchn/analyses/sampling/sppidTOsppname.csv")

for(i in tfiles){
  
  ## Dataset name
  dnm <- gsub("/home/spreinalesl/Dropbox/Sandra/Manuscritos/NGS_phyloOchn/analyses/results/phylonetworks/starting_tree/sauv_",
              i, replacement = "")
  dnm <- gsub("_BS33_sptree.tre", dnm, replacement = "")
  
  ## Read tree
  tre <- read.tree(i)
  is.rooted(tre)
  
  ## Change tip labels
  for(j in 1:length(tre$tip.label)){
    newID <- smpdata$seqID[which(smpdata$abbrev == tre$tip.label[j])]
    if(length(newID) > 1){
      cat(newID, "\n")
      newID <- grep("SRR", newID, value = T)
    }
    tre$tip.label[j] <- newID
  }
  
  ## Drop tips
  tipTOdel <- which(!(tre$tip.label %in% subsetspp$ID))
  tre <- drop.tip(tre, tip = tipTOdel)
  
  ## Write new tree in newick
  write.tree(tre, 
             file = paste("/home/spreinalesl/Dropbox/Sandra/Manuscritos/NGS_phyloOchn/analyses/results/phylonetworks/starting_tree/",
                               dnm, ".tre", sep = ""))
  
  ## Reroot tree
  tre <- root.phylo(tre, outgroup = "SRR10951426", resolve.root = T)
  
  for(k in tre$tip.label){
    tre$tip.label[which(tre$tip.label == k)] <- smpdata$spp_name_voucher[which(smpdata$seqID == k)]
  }
  
  pdf(paste("/home/spreinalesl/Dropbox/Sandra/Manuscritos/NGS_phyloOchn/analyses/results/phylonetworks/starting_tree/",
            dnm, ".pdf", sep = ""))
   plot.phylo(ladderize(tre), cex = 0.6)
  dev.off()
}