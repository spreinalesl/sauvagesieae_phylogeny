###########################################
## This script detects and remove extreme 
## outliers from the alignments
###########################################

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
  
  ## Calculate sequence length
  seq_lenght <- vector()
  for(j in 1:length(nameLines)){
    seq_lenght[j] <- nchar(paste0(fastaf[up[j]:bottom[j]], collapse = ""))  
  }
  
  ## Calculate Q1, Q3, IQR and extreme outliers (just shortest the sequences)
  # just if the sequence is shorter than 250 bp
  diqr <- quantile(seq_lenght, probs = 0.75) - quantile(seq_lenght, probs = 0.25)
  minv <- floor(quantile(seq_lenght, probs = 0.25) - 3*diqr)
  #if(minv > 300){ minv <- 300 }
  
  ## Extreme outliers
  extremv <- which(seq_lenght < minv)
  
  ## Generate list with locus name and sample name which has extreme data
  samplesToRemove <- unlist(lapply(X = fastaf[nameLines[extremv]],
                                   FUN = function(x) strsplit(x, split = " ")[[1]][1]))
  samplesToRemove <- gsub(samplesToRemove, pattern = ">", replacement = "")
  
  ## Positions to remove
  to_remove <- unlist(sapply(extremv, FUN = function(x) c(nameLines[x]:bottom[x])))
  
  ## Write edited fasta file,
  # loci without outliers are written as well, just to an easy handling forward
  if(length(to_remove) > 0){
    write(fastaf[-to_remove],
          file = paste0(fpath, gsub(i, pattern = "\\.FNA", replacement = "_ed.FNA"), collapse = ""))  
  }else{
    write(fastaf, file = paste0(fpath, gsub(i, pattern = "\\.FNA", replacement = "_ed.FNA"), collapse = ""))
  }
  
  cat("Editing", probename, i, "\n")
  
  ## lstats contains final number of samples per loci, samples removed (i.e. extreme values),
  # and max and min length allowed for this alignment. 
  lstats <- paste(length(nameLines) - length(samplesToRemove),
                  length(samplesToRemove),
                  max(seq_lenght), minv, sep = "-")
  
  ## Write summary table
  write(paste(gsub(i, pattern = "\\.FNA", replacement = ""), lstats,
               paste0(samplesToRemove, collapse = ","), sep = "\t"),
        file = paste0(outpath, out2, collapse = ""), sep = "\n", append = T)

 cat("Writing summary table for ", i, "\n")
}

cat("END of the run")
