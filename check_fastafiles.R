library(ape)

flist <- list.files("~/Documents/subsampling/ffed", pattern = "FNA")
length(flist)

summtable <- list()
for(i in 1:length(flist)){
  
  ff <- read.FASTA(paste("~/Documents/subsampling/ffed", flist[i], sep = "/"))
  ltips <- gsub(labels(ff), pattern = " .+", replacement = "")
  
  #if(length(which(grepl("_R_SPRL_145", ltips))) > 0){ 
  #  ltips <- gsub("_R_SPRL_145", ltips, replacement = "SPRL_145")
  #}
  
  #summtable[[i]] <- data.frame(tipnames = ltips)
  #gnm <- data.frame(gname = gsub(pattern = "taper_", flist[i], replacement = ""))
  #gnm <- data.frame(gname = gsub(pattern = "_align.fas", gnm, replacement = ""))
  gnm <- data.frame(gname = gsub(pattern = "\\.FNA", flist[i], replacement = ""))
  summtable[[i]] <- gnm
  
  summtable[[i]]$ntips <- length(ff)
  
  seql <- vector() 
  for(j in 1:length(ff)){
    seql[j] <- length(ff[[j]])
  }
  
  summtable[[i]]$shortest_seq <- min(seql)
  summtable[[i]]$longest_seq <- max(seql)
  summtable[[i]]$median_seq <- median(seql)
  #summtable[[i]]$gname <- rep(gnm, length(ltips))
}

summtable <- do.call(rbind, summtable)

View(summtable)

tiptable <- as.data.frame(table(summtable$tipnames))

tiptable$sppnames <- NA
sppnames <- read.csv("~/Dropbox/Sandra/Manuscritos/NGS_phyloOchn/analyses/sampling/sppidTOsppname.csv")
for(i in unique(tiptable$Var1)){
  tiptable$sppnames[which(tiptable$Var1 == i)] <- sppnames$spp_name_voucher[sppnames$seqID == i]
}

write.table(tiptable, "~/Documents/rmout_infprlgs_tipstable.csv", row.names = F)

View(table(summtable$ntips))

write.csv(summtable, "~/Documents/subsampling/subsampling_summtable.csv", row.names = F)

