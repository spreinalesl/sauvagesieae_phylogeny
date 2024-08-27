map <- read.table("~/Documents/multiseqID_all.map")
map$V1 <- gsub(map$V1, pattern = "_R_", replacement = "")
map$V2 <- gsub(map$V2, pattern = "_R_", replacement = "")

map$V0 <- unlist(lapply(strsplit(map$V1, split = "\\."), FUN = function(x) x[2]))

dnames <- read.csv("~/Dropbox/Sandra/Manuscritos/NGS_phyloOchn/analyses/sampling/sppidTOsppname.csv")

for(i in unique(map$V2)){
  map$V2[which(map$V2 == i)] <- dnames$abbrev[which(dnames$seqID == i)]
}

map$V3 <- paste(map$V2, map$V0, sep = ".")
map$V3 <- gsub(pattern = "\\.NA", map$V3, replacement = "")

write.table(map[,c(4,2)], row.names = F, "~/Documents/multiseqID_all_nm_ed.map")