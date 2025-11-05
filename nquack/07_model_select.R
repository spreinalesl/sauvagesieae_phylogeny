## Load packages
library(nQuack)

## Paths
inpath <- "~/projects/nquack/models/"
outpath <- "~/projects/nquack/models/interpret/"

## Iterate files
samples <- list.files(pattern = ".csv", inpath)

## Select models based on the BIC score and only considering
## diploid, triploid, and tetraploid mixtures.
summary <- c()
bestModel <- c()
mixtures <- c("diploid", "triploid", "tetraploid")

for(i in 1:length(samples)){
  temp <- read.csv(paste0(inpath, samples[i]))
  
  ## Extract info for all models  
  summary[[i]] <- quackit(model_out =  temp,
                          summary_statistic = "BIC",
                          mixtures = mixtures)
  
  ## Extract info for the most common mixture among models
  mcmx <- table(summary[[i]]$winnerBIC)
  max_mcmx <- max(mcmx)
  mcmx <- names(mcmx[which(mcmx == max_mcmx)])

  # in case of ties
  if(length(mcmx) > 1){
    mcmx <- paste0(mcmx, collapse = "-")
  }
  
  ## Extract info for the best models
  temp <- subset(temp, temp$mixture %in% mixtures) # subset table to include just selected mixtures
  pminBIC <- which(temp$BIC == min(temp$BIC))
  temp2 <- temp[-pminBIC,] # Remove the best model to identify the second best
  psecondBIC <- which(temp2$BIC == min(temp2$BIC))
  bestModel[[i]] <- data.frame(sample = gsub(".csv", "", samples[i]),
                               BIC = min(temp$BIC),
                               BestMixture = temp$mixture[pminBIC],
                               BestDist = temp$distribution[pminBIC],
                               BestType = temp$type[pminBIC],
                               SecondMixture = temp2$mixture[psecondBIC],
                               SecondDist = temp2$distribution[psecondBIC],
                               SecondType = temp2$type[psecondBIC],
                               BICdif = (temp2$BIC[psecondBIC] - temp$BIC[pminBIC]),
                               MostCommonMixture = mcmx,
                               MostCommonMixtureTimes = max_mcmx)
  
  cat("Extracting best model for sample", gsub(".csv", "", samples[i]), "\n")
}

## Write tables
summary <- do.call(rbind, summary)
bestModel <- do.call(rbind, bestModel)

write.csv(summary, file = paste0(outpath, "summary.csv"), row.names = FALSE)
cat("Writing summary table", "\n")

write.csv(bestModel, file = paste0(outpath, "bestModel.csv"), row.names = FALSE)
cat("Writing best model table", "\n")