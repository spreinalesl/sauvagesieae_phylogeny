###########################################
## This script detects and remove extreme 
## outliers from the alignments
###########################################

fpath <- "~/Dropbox/Sandra/Manuscritos/NGS_phyloOchn/analyses/results/hybpiper/fasta_files/a353/exons/"
outpath <- "~/Documents/"

## Packages
library(ape)

## Input paths
probename <- "a353"
type <- "exons"
fpath <- paste0("~/projects/rbgk/users_area/sreinales/results/sauvagesieae/hybpiper/fasta_files",
                "/", probename, "/", type, "/", collapse = "")

# output files
outpath <- "~/projects/rbgk/users_area/sreinales/results/sauvagesieae/hybpiper/get_seqs_summ/"
out1 <- paste0("lociNOseqs_", probename, ".txt", collapse = "")
out2 <- paste0("samplesTOremove_", probename, ".csv", collapse = "")

## Read data and iterate fasta files containing all sequences per loci, to
# remove samples with shortest sequences
ffiles  <- list.files(fpath, pattern = "\\.FNA")

# To save summary
lstats <- data.frame(loci = vector(),
                     sample_num = vector(),
                     rm_samples = vector(),
                     out_num = vector(),
                     ambp_num = vector(),
                     max_length = vector(),
                     min_allw_length = vector(),
                     mean_length = vector(),
                     rm_samples_names = vector())

# Iterate files
for(i in ffiles){
  
  ## Extract locus name
  locus_name <- gsub(i, pattern = "\\.FNA", replacement = "")
  
  ## Read fasta file
  fasta <- read.FASTA(paste0(fpath, i, collapse = ""))
  
  ## Skip the current iteration because the alignment is empty
  if(length(fasta) == 0){
    write(locus_name, file = paste0(outpath, out1, collapse = ""), sep = "\n", append = T)
    next
  }
  
  ## Calculate extreme outliers using the quantile approach
  # Values above Q3+1.5xIQR or below Q1-1.5xIQR are considered as outliers.
  # Values above Q3+3xIQR or below Q1-3xIQR are considered as extreme outliers.
  # Q1 and Q3 are the first and third quartile, IQR is the interquartile range (IQR = Q3 - Q1).
  # Outliers in a boxplot are not considered as troublesome as those considered extreme points
  # and might even be ignored.
  seq_length <- sapply(fasta, FUN = length)
  q1 <- quantile(seq_length, probs = 0.25)
  q3 <- quantile(seq_length, probs = 0.75)
  diqr <- q3 - q1
  minv <- floor(q1 - 3*diqr)
  
  ## Sequences shorter than minv or shorter than 150 bp, whichever was longer, were removed
  ## For shorter loci, sequences with at least 80bp will be maintained
  if(minv < 150){
    if(median(seq_length) < 650){
      minv <- 80
    }else
      minv <- 150
  }
  
  ## List of samples with sequences shorter than the min allowed length (extreme outliers)
  samples_extremv <- which(seq_length < minv)
  
  ## Count N in each sequence and save which have more than 80% of its length 
  # of ambiguous positions (N)
  samples_ambp <- vector()
  for(j in 1:length(fasta)){
    ambp <- length(which(fasta[[j]] == "f0"))
    if(ambp/length(fasta[[j]]) > 0.8){
      samples_ambp <- c(samples_ambp, names(fasta[j]))
    }
  }
  
  ## If there are samples with extreme values or samples
  ## with many ambiguous positions, remove them and write
  ## edited fasta files
  
  if(length(samples_extremv) > 0 | length(samples_ambp) > 0){
    
    ## List of samples to remove
    samplesToRemove <- unique(c(names(samples_extremv), names(samples_ambp)))
    to_remove <- which(names(fasta) %in% samplesToRemove)
    
    ## Write edited fasta file,
    # loci without outliers are written as well, just to an easy handling forward
    write.FASTA(fasta[-to_remove],
                file = paste0(fpath, gsub(i, pattern = "\\.FNA", replacement = "_ed.FNA"), collapse = ""))
    
    ## lstats contains final number of samples per loci, samples removed (i.e. extreme values),
    ## samples removed for N, and max and min length allowed for this alignment.
    snames <- unlist(lapply(X = samplesToRemove,
                            FUN = function(x) strsplit(x, split = " ")[[1]][1]))
    snames <- gsub(snames, pattern = ">", replacement = "")
    
    lstats <- rbind(lstats, data.frame(loci = locus_name,
                              sample_num = length(fasta) - length(samplesToRemove),
                              rm_samples = length(samples_extremv) + length(samples_ambp),
                              out_num = length(samples_extremv),
                              ambp_num = length(samples_ambp),
                              max_length = max(seq_length),
                              min_allw_length = minv[[1]],
                              mean_length = ceiling(mean(seq_length)),
                              rm_samples_names = paste0(snames, collapse = ",")))
  } else{
    lstats <- rbind(lstats, data.frame(loci = locus_name,
                                       sample_num = length(fasta),
                                       rm_samples = 0, out_num = 0, ambp_num = 0,
                                       max_length = max(seq_length),
                                       min_allw_length = minv[[1]],
                                       mean_length = ceiling(mean(seq_length)),
                                       rm_samples_names = NA))
    
     write.FASTA(fasta, file = paste0(fpath, gsub(i, pattern = "\\.FNA", replacement = "_ed.FNA"), collapse = ""))
  }
  
  cat("Editing", i, "\n")
}

## Write summary table
write.table(lstats, file = paste0(outpath, out2, collapse = ""), row.names = F)
