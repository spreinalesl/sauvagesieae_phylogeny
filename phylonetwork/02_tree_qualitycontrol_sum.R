###################################
## This script summarize sdsf values
## across trees using the table generated
## by 02_tree_qualitycontrol_gen.sh 
###################################

## Input table path
tpath <- "~/scratch/private/phylonetworks/mbresults"

## Read table
sdsf <- read.table(paste(tpath, "topostability.txt", sep = "/"),
                   header = FALSE, sep = " ", stringsAsFactors = FALSE)

head(sdsf)

## how many files with sdsf above 0.1?
goodsdsf <- sum(sdsf$V3 < 0.1, na.rm = TRUE)

## Files with sdsf>0.1?
badsdsf <- paste(sdsf$V1[which(sdsf$V3 > 0.1)], sdsf$V3[which(sdsf$V3 > 0.1)], sep = "\t") 

## Write log
cat(goodsdsf, " trees of ", nrow(sdsf), " have sdsf<0.1", "\n",
    "Mean sdsf ", mean(sdsf$V3), "\n", "\n",
    "Following trees do not achieved convergence (sdsf>0.1):", "\n", paste0(badsdsf, collapse = "\n"),
    sep = "", file = paste(tpath, "sdsf_warning.out", sep = "/"))

## visualize the dist of sdsf
pdf(paste(tpath, "sdsf_alltrees.pdf", sep = "/"))
boxplot(V3 ~ V1, data = sdsf, xaxt="n",
        xlab = "Trees", ylab = "Average standard deviation of split frequencies")
dev.off()

cat("Run end")
