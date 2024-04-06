#################################################
## This script detects and removes extreme 
## outliers from the alignments, i.e. extremely
## short sequences and sequences with high amount
## of unknown positions 'N'
##################################################

## Packages
library(ape)

## Input paths
namesfile <- "sppidTOsppname.csv"
probename <- "a353"
type <- "exons"
namespath <- "~/projects/rbgk/users_area/sreinales/results/sauvagesieae/name_list/"
fpath <- paste0("~/projects/rbgk/users_area/sreinales/results/sauvagesieae/hybpiper/fasta_files",
                "/", probename, "/", type, "/", collapse = "")

# output files
outpath1 <- "~/projects/rbgk/users_area/sreinales/results/sauvagesieae/hybpiper/get_seqs_summ/"
outpath2 <- paste(fpath, "rm_out/", sep = "")
out1 <- paste0("lociNOseqs_", probename, ".txt", collapse = "")
out2 <- paste0("samplesTOremove_", probename, ".csv", collapse = "")
out3 <- paste0("rm_samples_", probename, ".csv", collapse = "")

## Read data and iterate fasta files containing all sequences per loci, to
# remove samples with shortest sequences
ffiles  <- list.files(fpath, pattern = "\\.FNA$")

# To save summary
lstats <- data.frame(loci = vector(),
                     smp_num_org = vector(),
                     smp_num = vector(),
                     rm_samples = vector(),
                     out_num = vector(),
                     ambp_num = vector(),
                     max_length = vector(),
                     min_allw_length = vector(),
                     mean_length = vector(),
                     rm_out_smp_names = vector(),
                     rm_out_smp_meanlength = vector(),
                     rm_ambp_smp_names = vector(),
                     rm_ambp_smp_meanlength = vector())

# Iterate files
for(i in ffiles){
  
  ## Extract locus name
  locus_name <- gsub(i, pattern = "\\.FNA", replacement = "")
  
  ## Read fasta file
  fasta <- read.FASTA(paste0(fpath, i, collapse = ""))
  
  ## Extract number of samples included in the fasta file
  smp_num_org <- length(fasta)
  
  ## Skip the current iteration because the fasta file is empty (i.e. no sequences were retriueved for this locus)
  if(length(fasta) == 0){
    write(locus_name, file = paste0(outpath1, out1, collapse = ""), sep = "\n", append = T)
    next
  }
  
  ## Check if more than 50% of the sequence length is unknown positions (N), if TRUE, count the number of
  ## known positions (A,C,T,G) between Ns. Tag sequences having less than 150 known positions
  ## between Ns. 
  ambp <- vector()
  length_ambp <- vector()
  
  for(j in 1:length(fasta)){
    sq <- fasta[[j]]
    psum <- 0
    
    if(length(which(fasta[[j]] == "f0"))/length(fasta[[j]]) > 0.5){
      for(k in sq){
        if(!(k == "f0")){
          psum <- psum + 1
          
          if(psum > 150){
            break
          }
          
        }else{psum <- 0}
      }
      
      if(psum < 150){
        ambp <- c(ambp, j)
        length_ambp <- length(fasta[[j]])
      }
    }
  }
    
  ## Remove sequences having less than 150 known positions between Ns from the fasta
  # before the estimation of median length of the sequences
  if(length(ambp) > 0){
    smp_ambp <- names(fasta[ambp])
    fasta <- fasta[-ambp]
  }else{
    smp_ambp <- vector()
  }

  ## Calculate sequences length
  seq_length <- sapply(fasta, FUN = length)
  
  ## Calculate extreme outliers using the quantile approach
  # Values above Q3+1.5xIQR or below Q1-1.5xIQR are considered as outliers.
  # Values above Q3+3xIQR or below Q1-3xIQR are considered as extreme outliers.
  # Q1 and Q3 are the first and third quartile, IQR is the interquartile range (IQR = Q3 - Q1).
  # Outliers in a boxplot are not considered as troublesome as those considered extreme points
  # and might even be ignored.
  #q1 <- quantile(seq_length, probs = 0.25)
  #q3 <- quantile(seq_length, probs = 0.75)
  #diqr <- q3 - q1
  #minv <- floor(q1 - 3*diqr)
  
  ## Establish 1/4 of the median sequence length for that locus as the minimum allowed seq length,
  # just if that minimum is shortest than 150 bp, but for shortest alignments (i.e. 200 bp),
  # the min allowed length could be 70bp. The idea is to avoid removing sequences that can be aligned.
  minv <- ceiling(median(seq_length)*0.25)
  
  if(minv > 150){
      minv <- 150
  }else{
    minv <- 70
  }
  
  ## List of samples with sequences shorter than the min allowed length (extreme outliers)
  smp_extremv <- names(which(seq_length < minv))
  extremvp <- which(names(fasta) %in% smp_extremv) 
  
  ## If there are samples with extreme values remove them and write
  ## the edited fasta files
  
  if(length(smp_extremv) > 0){
    fasta <- fasta[-extremvp]
  }
   
  ## Change sample names from ID to spp names file
  smp_names <- read.csv(paste0(namespath, namesfile, collapse = ""), as.is = T)
  
  if(length(smp_extremv) > 0 | length(smp_ambp) > 0){
    
    nextremv <- vector()
    nambp <- vector()
    
    for(fnx in smp_extremv){
      fname <- gsub(fnx, pattern = " .+", replacement = "")
      nextremv <- c(nextremv, smp_names$abbrev[which(smp_names$seqID == fname)])
    }
    
    for(fnb in smp_ambp){
      fname <- gsub(fnb, pattern = " .+", replacement = "")
      nambp <- c(nambp, smp_names$abbrev[which(smp_names$seqID == fname)])
    }
    
    ## Extract a list of samples to remove
    smp_extremv <- nextremv
    smp_ambp <- nambp 
    samplesToRemove <- unique(c(smp_extremv, smp_ambp))  
  }
  
  ## lstats contains final number of samples per loci, samples removed (i.e. extreme values),
  # samples removed for N, and max and min length allowed for this alignment.  
  lstats <- rbind(lstats, data.frame(loci = locus_name,
                                     smp_num_org = smp_num_org,
                                     smp_num = length(fasta),
                                     rm_samples = length(smp_extremv) + length(smp_ambp),
                                     out_num = length(smp_extremv),
                                     ambp_num = length(smp_ambp),
                                     max_length = max(seq_length),
                                     min_allw_length = minv[[1]],
                                     mean_length = ceiling(mean(seq_length)),
                                     rm_out_smp_names = paste0(smp_extremv, collapse = ","),
                                     rm_out_smp_length = paste0(seq_length[extremvp], collapse = ","),
                                     rm_ambp_smp_names = paste0(smp_ambp, collapse = ","),
                                     rm_ambp_smp_length = paste0(length_ambp, collapse = ",")))
    
  ## Write edited fasta. Loci which no samples were removed, are written as well just to 
  # ease handling of the fasta files in the forward analyses
  write.dna(fasta, format = "fasta", nbcol = 1, colsep = "", colw = 60, indent = 0, blocksep = 0,
            file = paste0(outpath2, gsub(i, pattern = "\\.FNA", replacement = "_ed.FNA"), collapse = ""))

  cat("Editing", i, "\n")
}

## Extract list of removed samples
rm_samples_names <- write.table(as.data.frame(table(c(unlist(strsplit(lstats$rm_out_smp_names, split = ",")),
                                                      unlist(strsplit(lstats$rm_ambp_smp_names, split = ","))))),
                                file = paste0(outpath1, out3, collapse = ""), row.names = F)

## Write summary table
write.table(lstats, file = paste0(outpath1, out2, collapse = ""), row.names = F)
