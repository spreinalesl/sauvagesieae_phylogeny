## Load packages
library(nQuack)
library(dplyr)

## Paths
inpath <- "/home/sreinale/scratch/private/nquack/example/sauv/models/"
outpath <- "/home/sreinale/scratch/private/nquack/example/sauv/models/interpret/"

## Iterate files
samples <- list.files(pattern = ".csv", inpath)

## Select models based on the BIC score and only considering
## diploid, triploid, and tetraploid mixtures.
summary <- c()
bestModel <- c()
for(i in 1:length(samples)){
  temp <- read.csv(paste0(inpath, samples[i]))

  ## Extract info for all models  
  summary[[i]] <- quackit(model_out =  temp,
                          summary_statistic = "BIC",
                          mixtures = c("diploid", "triploid", "tetraploid"))

  ## Extract info for the most common mixture among models
  mcmx <- table(summary[[i]]$winnerBIC)
  max_mcmx <- max(mcmx)
  mcmx <- names(mcmx[which(mcmx == max_mcmx)])

  ## Extract info for the best model
  pminBIC <- which(temp$BIC == min(temp$BIC))
  bestModel[[i]] <- data.frame(sample = gsub(".csv", "", samples[i]),
                               BIC = min(temp$BIC),
                               BestMixture = temp$mixture[pminBIC],
                               BestDist = temp$distribution[pminBIC],
                               BestType = temp$type[pminBIC],
                               MostCommonMixture = mcmx,
                               MostCommonMixtureTimes = max_mcmx)

  cat("Extracting information for sample", samples[i], "\n")
}

## Write tables
summary <- do.call(rbind, summary)
bestModel <- do.call(rbind, bestModel)

write.csv(summary, file = paste0(outpath, "summary.csv"), row.names = FALSE)
cat("Writing summary table", "\n")

write.csv(bestModel, file = paste0(outpath, "bestModel.csv"), row.names = FALSE)
cat("Writing best model table", "\n")
