## Load packages
library(nQuack)

## Paths
outpath <- "/home/sreinale/scratch/private/nquack/example/sauv/processed/"

## Process data following the filtering approach found in nQuire:
## a minimum depth of 10, allele trunctation minimum of 0.15, and
## allele truncation maximum of 0.85.
## and then, remove noise using one of the following alternatives:
## i) denoise_data() which filters allele frequencies with a normal + uniform mixture model
## ii) Bclean() which utilizes the U-shaped beta distribution 

# choose method
denoiseD <- T
filterm <- "nquire"

if(denoiseD){
 tmin <- 0.15 
 tmax <- 0.85
 mnm <- "dd"  
}else{
 tmin <- 0
 tmax <- 0
 mnm <- "bc"
}

## Iterate files using an array job
samp <- commandArgs(trailingOnly = T)

## Filter data
temp <- process_data(paste0(outpath, samp),
                     min.depth = 10, 
                     max.depth.quantile.prob = 1, 
                     error = 0.01, 
                     trunc = c(tmin, tmax))
  
## write a table
samp <- gsub(".txt", "", samp)
write.csv(temp, file = paste0(outpath, samp, "_", mnm, ".csv"),
          row.names = FALSE)

cat("Writing filtered data for", samp, "\n")
 
## Remove noise 
if(denoiseD){
  ## Using denoise_data
  pdf(paste(outpath, "denoiseD_", samp, ".pdf", sep = ""))
    denoise_data(temp, plot = T, filter = filterm)
  dev.off()
  cat("END denoise", " \n")
}else{
  ## Using Bclean
  pdf(paste(outpath, "bclean_", samp, ".pdf", sep = ""))
   Bclean(temp, plot = TRUE, quick = F)
  dev.off()
  cat("END Bclean", "\n")
}
