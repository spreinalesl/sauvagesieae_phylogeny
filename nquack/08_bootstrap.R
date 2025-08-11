## Load packages
library(nQuack)

## Paths
inpath <- "/home/sreinale/scratch/private/nquack/example/sauv/processed/"
modpath <- "/home/sreinale/scratch/private/nquack/example/sauv/models/interpret/"
outpath <- "/home/sreinale/scratch/private/nquack/example/sauv/bootstrap/"

## Iterate files using an array job
sample <- commandArgs(trailingOnly = T)

## Sample name
sample <- gsub("_dd.csv", "", sample)
cat("Begin bootstrap for sample", sample, "\n")

## Read summary table and extract info for the best model
mod <- read.csv(paste0(modpath, "bestModel.csv"))

## To run bootstrap replicates under the best or second best model
## first define the following options:
## Define the number of bootstrap replicates
nboots <- 1000

## Define best (TRUE) or second best (FALSE) model to run
bestm <- TRUE

## Read summary table and extract info for the best or the second best model
if(bestm){
 bestnm <- "best model"
 bestDist <- mod$BestDist[which(mod$sample == sample)]
 bestType <- mod$BestType[which(mod$sample == sample)]
} else{
 bestnm <- "second best model"
 bestDist <- mod$SecondDist[which(mod$sample == sample)]
 bestType <- mod$SecondType[which(mod$sample == sample)]
}

## Define if a uniform distribution for noise will be included in the model
if(grepl("uniform", bestDist)){
 bestDist <- gsub("-uniform", "", bestDist)
 bestUnif <- 1
}else{
 bestUnif <- 0
}

cat(paste0("Running under the ", bestnm, ":"), bestDist, bestType, bestUnif, "\n")

## Read table with allele frequencies      
temp <- as.matrix(read.csv(paste0(inpath, sample, "_dd.csv")))
bout <- quackNboots(temp,
                    nboots = nboots,
                    distribution = bestDist,
                    type = bestType,
                    uniform = bestUnif,
                    mixtures = c("diploid", "triploid", "tetraploid"),
                    samplename = sample)

cat("End bootstrap for sample", sample, "\n")

## Write table
write.csv(bout, file = paste0(outpath, sample, "_boot.csv"), row.names = F)
cat("Writing bootbestModel table", "\n")


############################## Without Bootstrap #########################

## Use the best model to predict their ploidal level with our bestquack() function
#for(i in 1:length(filelist)){
 ## Sample name
# sample <- gsub("_dd.csv", "", filelist[i])

 ## Read summary table and extract info for the best model
# bestDist <- mod$BestDist[which(mod$sample == sample)]
# bestType <- mod$BestType[which(mod$sample == sample)]
# bestUnif <- (grepl("uniform", mod$BestDist[which(mod$sample == sample)])*1)

 ## Read table with allele frequencies
# temp <- as.matrix(read.csv(paste0(inpath, filelist[i])))
# out[[i]] <- bestquack(temp, 
#                       distribution = bestDist,
#                       type = bestType,
#                       uniform = bestUnif,
#                       mixtures = c("diploid", "triploid", "tetraploid"), 
#                       samplename = sample)
#}

## Write table
#out <- do.call(rbind, out)
#write.csv(out, file = paste0(outpath, "runUnderBestModel.csv"), row.names = F)
#cat("Writing bestModel table", "\n")
