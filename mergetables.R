## Data path
dpath <- "~/Documents/plot_trees/hybpiper_all_summary.csv" 

## Read table with information about gene recovery and sequence length
sdata <- read.table(dpath, header = T, sep = ";", dec = ",")

hybp <- read.csv("~/Documents/plot_trees/hybphaser.csv",
                 sep = ";")

nms <- read.csv("~/Documents/plot_trees/sppidTOsppname.csv")

for(i in sdata$sample){
  sdata$sppName[which(sdata$sample == i)] <- nms$abbrev[which(nms$seqID == i)]
  sdata$name[which(sdata$sample == i)] <- nms$spp_name_voucher[which(nms$seqID == i)]
}

sdata$AD <- NA
sdata$LH <- NA

for(i in sdata$name){
  sdata$AD[which(sdata$name == i)] <- hybp$AD[which(hybp$sample == i)]
  sdata$LH[which(sdata$name == i)] <- hybp$LH[which(hybp$sample == i)]
}

write.csv(sdata, "~/Documents/plot_trees/hybpiper_all_summary_ed.csv", row.names = F)
