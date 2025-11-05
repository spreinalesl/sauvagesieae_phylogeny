######################################################
## This script generates new target files for a353
## and Ochnaceae, removing all shared genes between
## Angiosperms353 and Ochnaceae-specific set of probes
######################################################

#################################
## Ochnaceae-specific target file
#################################

## Read fasta file containing CKDK scaffolds (from OneKP)
filepath <- "~/projects/custom_target/edited_fasta/CKDK-SOAPdenovo-Trans-assembly_ed.fas"
scff <- readLines(filepath)
scff <- gsub("CKDK-","CKDK-g", scff)

## Read Ochnaceae gene IDs previously extracted from
## Schneider et al. (2021) bait sequences file
## and filter the 275 scaffolds from the fasta with
## all CKDK scaffolds (scff)
filepath <- "~/projects/ochnaceaeKit_geneIDs.txt"
ochnIDs <- readLines(filepath)

## First, remove the shared gene 2094005 for ochnaceaeKit_geneIDs
# see details in ./notes
ochnIDs <- ochnIDs[-grep("2094005", ochnIDs)]

## Filter scaffold seqs
filepath <- "~/projects/custom_target/newTarget_files/"

for(i in ochnIDs){
  namep <- grep(i, scff)
  seqp <- namep+1
  
  if(length(namep) == 1){
    p <- sort(c(namep, seqp))
    write(scff[p], file = paste(filepath, "OCHNA-input-seq_Targets_DNA_BaitKit_edit_crop.fas", sep = ""),
          append = T)
  }
  cat("Filtering gene", i, "\n")
}

#################################
## Angiospems353
#################################
## Read mega353 target file download from hybpiper
# after filtering by Malpighiales samples
# for detais see the hybpiper script 'filter_mega353'
filepath <- "~/projects/mega353_filtered.fasta"
mega353 <- readLines(filepath)

## Remove shared genes except g5943
filepath <- "~/projects/custom_target/blast/gene_match.txt"
shared <- read.table(filepath, sep = ",", header = F)
shared <- gsub("g", "", shared$V1[-grep("g5943", shared$V1)])

## List positions to be removed
positions <- list()
for(i in 1:length(shared)){
  p <- grep(paste("-", shared[i], sep = ""), mega353, fixed = T)
  p <- sort(c(p, p+1)) ## Add sequence position
  positions[[i]] <- p
}

positions <- sort(unlist(positions))

## Remove positions
filepath <- "~/projects/custom_target/newTarget_files/"
write(mega353[-positions], file = paste(filepath, "mega353_filtered_crop.fas", sep = ""))

#################################
## Check the number of genes
#################################
## Read each file
filepath <- "~/projects/custom_target/newTarget_files/"
filepath <- list.files(pattern = ".fas", path = filepath, full.names = T)

for(i in 1:length(filepath)){
  
  fas <- readLines(filepath[i])
  
  ## Select just the lines containing the symbol > (i.e. taxon names)
  lines <- grep(">", x = fas, value = T)
  listl <- strsplit(lines, split = "-")
  dataf <- data.frame(sample = rep(NA, length(listl)), gene = NA)
  
  for(j in 1:length(listl)){
    dataf$sample[j] <- listl[[j]][1]
    dataf$gene[j] <- listl[[j]][2] 
  }
  
  ## Write table
  write.table(table(dataf$sample), file = paste(filepath[i], "_genesPeRsample.csv", sep = ""), row.names = F)
  writeLines(unique(dataf$gene), sep = "\n", con = paste(filepath[i], "_geneList.txt", sep = ""))
}