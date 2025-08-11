## Paths
inpath <- "/home/sreinale/scratch/private/nquack/example/sauv/bootstrap/"
IDtosppnm <- "/home/sreinale/scratch/private/nquack/example/sauv/sppidTOsppname.csv"
modfile <- "/home/sreinale/scratch/private/nquack/example/sauv/models/interpret/bestModel.csv"

## Read species names
spp <- read.csv(IDtosppnm)

## Read file with info about models
models <- read.csv(modfile)

## Iterate files
boots <- list.files(pattern = "boot.csv", inpath)

## Read files
## This table includes two rows, the first show the best model
## for the original data set, and the second tallies the
## bootstrap replicates best model
summtable <- data.frame(sample = NA,
                        BestModel = NA,
                        BestModelPloidialLevel = NA,
                        BSsupport = NA,
                        AlternativePloidalLevel = NA)

for(i in 1:length(boots)){
  booti <- read.csv(paste0(inpath, boots[i]), header = T)
  
  ## Sample name
  summtable[i, "sample"] <- unique(booti$sample)
  
  ## Extract best model info
  summtable[i, "BestModel"] <- paste(models$BestDist[which(models$sample == unique(booti$sample))],
                                     models$BestType[which(models$sample == unique(booti$sample))],
                                     sep = "_")

  ## Extract inferred ploidal and bootstap support
  # Bootstrap replicates
  bsrep <- sum(booti[2,-4], na.rm = T)
  
  # Ploidal level inferred by the best model
  BestModelPloidialLevel <- colnames(booti)[which(!(apply(booti[1,-4], 2, is.na)))]
  summtable[i, "BestModelPloidialLevel"] <- BestModelPloidialLevel
  
  # Check the BS support for the best model ploidal level
  BSsupport <- round((booti[2, BestModelPloidialLevel]/bsrep)*100, digits = 1)
  summtable[i, "BSsupport"] <- BSsupport 
  
  if(BSsupport < 100){
    AlternativePloidalLevel <- colnames(booti)[which(!(apply(booti[2,-4], 2, is.na)))]
    AlternativePloidalLevel <- AlternativePloidalLevel[!(AlternativePloidalLevel %in% BestModelPloidialLevel)]
    if(length(AlternativePloidalLevel) > 1){
      second <- paste(AlternativePloidalLevel[1], "[", 
                      (booti[2, AlternativePloidalLevel[1]]/bsrep)*100, "]", sep = "")
      third <- paste(AlternativePloidalLevel[2], "[", 
                     (booti[2, AlternativePloidalLevel[2]]/bsrep)*100, "]", sep = "")
      summtable[i, "AlternativePloidalLevel"] <- paste0(second, third, collapse = "-")
    }else{
      summtable[i, "AlternativePloidalLevel"] <- paste(AlternativePloidalLevel, "[", 
                                                       (booti[2, AlternativePloidalLevel]/bsrep)*100, "]", sep = "")
    }
  }else{
    summtable[i, "AlternativePloidalLevel"] <- "NA"
  }
  
  cat("Extracting summary BS for", unique(booti$sample), "\n")
}

## Change sample names
for(i in 1:nrow(summtable)){
  summtable$spp[i] <- spp$spp_name_voucher[which(spp$seqID == summtable$sample[i])]
}

## Write summary table
write.csv(summtable, file = paste0(inpath, "BSsummary.csv"), row.names = F)
cat("Writing BSsummary table", "\n")
