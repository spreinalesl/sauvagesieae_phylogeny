###########################################
## This script detects and remove extreme 
## outliers from the alignments
###########################################

fpath <- "~/Dropbox/Sandra/Manuscritos/NGS_phyloOchn/analyses/results/hybpiper/fasta_files/a353/exons/"
outpath <- "~/Documents/"

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

for(i in list.files(fpath, pattern = "\\.FNA")){
  
  ## Extract locus name
  locus_name <- gsub(i, pattern = "\\.FNA", replacement = "")
  
  ## Read i fasta file
  fastaf <- readLines(paste0(fpath, i, collapse = ""))
  
  ## Skip the current iteration because the alignment is empty
  if(length(fastaf) == 0){
    write(locus_name, file = paste0(outpath, out1, collapse = ""), sep = "\n", append = T)
    next
  }
  
  ## Positions where seqs are
  nameLines <- grep(">", fastaf)
  up <- nameLines+1
  bottom <- (nameLines[-1])-1
  bottom <- c(bottom, length(fastaf))
  
  ## Calculate extreme outliers using the quantile approach
  # Values above Q3+1.5xIQR or below Q1-1.5xIQR are considered as outliers.
  # Values above Q3+3xIQR or below Q1-3xIQR are considered as extreme outliers.
  # Q1 and Q3 are the first and third quartile, IQR is the interquartile range (IQR = Q3 - Q1).
  # Outliers in a boxplot are not considered as troublesome as those considered extreme points
  # and might even be ignored.
  
  ## Calculate sequence length and count N
  seq_length <- vector()
  samples_ambp <- vector()
  for(j in 1:length(nameLines)){
    nseq <- paste0(fastaf[up[j]:bottom[j]], collapse = "")
    seq_length[j] <- nchar(nseq)
    
    ## Sequences with more than 80% of ambiguous positions (N)
    ambp <- length(gregexpr("N", nseq)[[1]])
    if(ambp/nchar(nseq) > 0.8){
      samples_ambp <- c(samples_ambp, nameLines[j])
    }
  }
  
  ## Samples with ambiguous positions
  ambp <- which(nameLines %in% samples_ambp)
  
  ## Calculate Q1, Q3, IQR and extreme outliers (just shortest the sequences)
  # all sequences shorter than 150bp are removed regardless the alignment length
  diqr <- quantile(seq_length, probs = 0.75) - quantile(seq_length, probs = 0.25)
  minv <- floor(quantile(seq_length, probs = 0.25) - 3*diqr)
  
  if(minv < 150){
    if(median(seq_length) < 650){
      minv <- 80
    }else
      minv <- 150
  }
  
  ## Extreme outliers
  extremv <- which(seq_length < minv)
  
  ## If there are samples with extreme values or samples
  ## with many ambiguous positions, remove them and write
  ## edited fasta files
  
  if(length(extremv) > 0 | length(samples_ambp) > 0){
  
  ## List of samples to remove
  samplesToRemove <- c(fastaf[nameLines[extremv]], fastaf[nameLines[ambp]])
  samplesToRemove <- unlist(lapply(X = samplesToRemove,
                                   FUN = function(x) strsplit(x, split = " ")[[1]][1]))
  samplesToRemove <- gsub(samplesToRemove, pattern = ">", replacement = "")
  
  ## Positions to remove
  to_remove <- c(unlist(sapply(extremv, FUN = function(x) c(nameLines[x]:bottom[x]))),
                 unlist(sapply(ambp, FUN = function(x) c(nameLines[x]:bottom[x]))))
  
  ## Extract unique positions because one sample could be fit both conditions:
  ## has extreme outliers and many N positions at the same time
  to_remove <- unique(to_remove)

  ## Write edited fasta file,
  # loci without outliers are written as well, just to an easy handling forward
  write(fastaf[-to_remove],
        file = paste0(fpath, gsub(i, pattern = "\\.FNA", replacement = "_ed.FNA"), collapse = ""))
  
  ## lstats contains final number of samples per loci, samples removed (i.e. extreme values),
  ## samples removed for N, and max and min length allowed for this alignment. 
  lstats <- paste(length(nameLines) - length(samplesToRemove),
                  length(fastaf[nameLines[extremv]]),
                  length(fastaf[samples_ambp]),
                  max(seq_length), minv, sep = "-")
  } else{
    lstats <- paste(length(nameLines), 0, 0, max(seq_length), minv, sep = "-")
    samplesToRemove <- NA
    write(fastaf, file = paste0(fpath, gsub(i, pattern = "\\.FNA", replacement = "_ed.FNA"), collapse = ""))
    }
  
  cat("Editing", i, "\n")
  
  ## Write summary table
  write(paste(gsub(i, pattern = "\\.FNA", replacement = ""), lstats,
               paste0(samplesToRemove, collapse = ","), sep = "\t"),
        file = paste0(outpath, out2, collapse = ""), sep = "\n", append = T)
}
