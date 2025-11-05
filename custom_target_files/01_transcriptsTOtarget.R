#############################################################
## This script reads transcript sequences for a taxon in the
## OneKP projet to establish the correspondence among
## Angiosperms353 and Ochnaceae-specific probes,
## and define a custom target file suitable to assemble reads
## produced by dual-hybridisation
############################################################

## Download and unzip the following folders for selected species from:
# https://datacommons.cyverse.org/browse/iplant/home/shared/commons_repo/curated/oneKP_capstone_2019/transcript_assemblies
# (i) XXXX-SOAPdenovo-Trans-translated.tar.bz2: A folder with transcript assemblies (dna and aa) containing parts of
# the scaffolds/contigs that matched orthogroups generated from annotated plant genomes.
# (ii) XXXX.orthoId_to_geneId.txt: A file with correspondence among assembly IDs and the names of the
# single/low-copy genes (orthologs) including those used to design the Angiospems353 kit, 
# (iii) CKDK-SOAPdenovo-Trans-assembly.fa.bz2: A file containing all scaffolds for Ochna serrulata
# Note: Not all scaffolds matched with orthologous seqs, so, some of the scaffolds could not blast against the
# assemblies (*Trans-assembly.dnas.out).

## Species:
# BNDE - Hypericum perforatum 
# CKDK - Ochna serrulata
# FWCQ - Garcinia oblongifolia
# KKDQ - Salix viminalis

## Define paths
dpath <- "~/projects/custom_target/onekp_data/"

# For edited files
nspath <- "~/projects/custom_target/edited_fasta/" 

## Create a list of samples to analyze
samples <- unique(gsub(pattern = "-.+", replacement = "",
                       x = list.files(dpath, pattern = "SOAPdenovo")))

###################################
## Get the correspondent gene name
## for each assembly
###################################

## Iterate folders to read *-SOAPdenovo-Trans-assembly.dnas.out files
# and match assembly names and gene names using CKDK.orthoId_to_geneId.txt  
folders <- dir(pattern = "Trans-translated", dpath)

for(i in folders){

  ## Create a table to compare assembly and gene names
  sname <- gsub(pattern = "-.+", replacement = "", x = i)
  
  ## Read .orthoId_to_geneId.txt file
  # assembly names are in the second column of names_file and 
  # the correspondent gene name is in the first column of the table
  names_file <- read.table(paste(dpath, sname, ".orthoId_to_geneId.txt", sep = ""),
                           sep = "\t", header = F)
  
  ## Read .dnas.out file containing assembly sequences
  files <- list.files(paste(dpath, i, sep = ""), full.names = T)
  seqs <- readLines(grep(".dnas.out", files, value = T))
  
  ## Get assembly names
  anames <- grep(">", seqs, value = T)
  length(anames)
  
  ## Search every assembly name in the name_file list, find the correspondent
  # gene name and replace it in the seq file
  afail <- NA
    
  for(j in anames){
    
   pt <- which(grepl(pattern = paste(gsub(">", " ", j), " ", sep = ""),
                     x = paste(" ", names_file$V2, " ", sep = "")))
   
   if(length(pt) == 1){
     
     ## If the assembly name is found, paste the correspondent gene name
     # in the seqs file
     seqs[which(seqs == j)] <- paste(">", sname, "-g", names_file$V1[pt], sep = "")
     
   }else{
     
     ## If the assembly name is missing, save the assembly name
     # in a vector
     afail <- c(afail, which(seqs == j))
   }

   cat("Searching", sname, j, "\n")
  }
  
  ## Write a table with missing assemblies
  afail <- afail[-1]
  writeLines(gsub(">", "", seqs[afail]), paste(nspath, sname, "_missing_assemblies.txt", sep = ""))
  
  ## Remove missing assemblies from seq file
  seqs <- seqs[-c(afail, afail+1)]
  
  ## Write a list of genes in the new seqs file
  genes <- unique(gsub(pattern = paste(">", sname, "-g", sep = ""), "", grep(">", seqs, value = T)))
  writeLines(genes, paste(nspath, sname, "_gene_list.txt", sep = ""), sep = "\n")
  
  ## Write new seq file
  writeLines(seqs, paste(nspath, i, "_ed.fas", sep = ""), sep = "\n")
}

###################################
## Read scaffold file for CKDK 
## and edit it for blast
###################################
scaffolds <- list.files(pattern = "SOAPdenovo-Trans-assembly.fa",
                        dpath)

for(i in scaffolds){
  seqs <- readLines(paste0(dpath, i, collapse = ""))
  
  ## Change scaffold names
  namesLines <- grep(">", seqs)
  scffnames <- seqs[namesLines]
  scffnames <- strsplit(scffnames, split = "-")
  scffnames <- unlist(lapply(scffnames,
                             FUN = function(x) paste0(x[c(2,3)], collapse = "-")))
  
  ## Write scaffold names
  writeLines(unique(scffnames), paste0(nspath, gsub(".SOAPdenovo-Trans-assembly.fa", "", i),
                                       "_scaffold_list.txt", collapse = ""))
  
  ## Positions where seqs are
  up <- namesLines+1
  bottom <- (namesLines[-1])-1
  bottom <- c(bottom, length(seqs))
  
  ## Convert scaffold name in fasta format 
  scffnames <- paste(">", scffnames, sep = "")
  
  ## Write edited fasta file
  for(j in 1:length(scffnames)){
 
   ## Paste modified name and seq
   write(paste(scffnames[j], paste0(seqs[up[j]:bottom[j]], collapse = ""), sep = "\n"),
         file = paste0(nspath, gsub(".fa", "", i), "_ed.fas", collapse = ""), append = T)
  
   cat("Editing", scffnames[j] , "\n")
  }
}