fpath <- "/home/spreinalesl/Documents/taper"
ffiles <- list.files(fpath, ".rfDistances", full.names = F)
summ <- data.frame()

for(i in ffiles){
  fname <- gsub(i, pattern = "_align.fas.raxml.mlTrees.raxml.rfDistances", replacement = "")
  ttable <- read.csv(paste(fpath, "/", i, sep = ""), sep = "\t", header = F)
  ttable$V0 <- fname
  summ <- rbind(summ, ttable)
}

summ2 <- cbind(as.data.frame(tapply(summ$V4, INDEX = as.factor(summ$V0), FUN = min)),
               as.data.frame(tapply(summ$V4, INDEX = as.factor(summ$V0), FUN = max)),
               as.data.frame(tapply(summ$V4, INDEX = as.factor(summ$V0), FUN = mean)))
summ2 <- rbind(summ2, mean = apply(summ2, MARGIN = 2, FUN = mean))
colnames(summ2) <- c("min", "max", "mean")
