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

# define truncation limits
tmin <- 0.15
tmax <- 0.85

# choose method
denoising <- TRUE
dmethod <- "denoisedata" # "bclean" 
filterm <- "both" # "nquire" "nquack"
plotfreq <- FALSE

## Iterate files using an array job
samp <- commandArgs(trailingOnly = T)

## Filter data
temp <- process_data(paste0(outpath, samp),
                     min.depth = 10, 
                     max.depth.quantile.prob = 1, 
                     error = 0.01, 
                     trunc = c(tmin, tmax))

## Change sample name
samp <- gsub(".txt", "", samp)

## If no denosing step will be run, write .csv table after truncation
## conversely, choose one method to denoise the raw data and write the
## table after runnig denoising
if(!denoising){
 write.csv(temp, file = paste0(outpath, samp, ".csv"),
           row.names = FALSE)

 cat("Writing filtered data for", samp, "\n")

}else{

 ## Open a plot window
 if(plotfreq){
  pdf(paste(outpath, "denoising_", samp, ".pdf", sep = ""))
 }

 ## Remove noise 
  if(dmethod == "denoisedata"){
    ## Using denoise_data
    mnm <- "dd"
    denoise_data(temp, plot = plotfreq, filter = filterm)
    cat("END denoise using method", filterm, "\n")
  }

  if(dmethod == "bclean"){
    ## Using Bclean
    mnm <- "bc"
    Bclean(temp, plot = plotfreq, quick = FALSE)
    cat("END Bclean", "\n")
   }

  ## Write table after denoising
  write.csv(temp, file = paste0(outpath, samp, "_", mnm, ".csv"),
            row.names = FALSE)

  cat("Writing filtered and denoised data for", samp, "\n")
  
  ## Close plot window
  if(plotfreq){dev.off()}
}
