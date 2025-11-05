##########################
## Extract best hits
## from blast results
##########################

## Define paths
filename <- "./blast/blastn_CKDK-ALL.csv"
taxname <- "CKDK-ALL"

## To check which of the query genes are in the a353 kit
# Read angiosperms353 gene IDs
a353IDs <- readLines("~/projects/target_files/a353Kit_geneIDs.txt")

## To check which of the matched scaffolds are in the Schn kit
# Read Schn gene IDs
schnIDs <- readLines("~/projects/target_files/ochnaceaeKit_geneIDs.txt")

## Read blast output
blasto <- read.csv(filename, header = F, as.is = T, sep = ",")
colnames(blasto) <- c("query_acc.ver","subject_acc.ver","identity",
                      "alignment_length","mismatches","gap_opens","q._start",
                      "q._end","s._start","s._end","evalue","bit_score") 

## Mark which genes are in 353 kit and which scaffolds are in schn kit
blasto$a353 <- NA
blasto$schn <- NA

## Transform columns
blasto$query_acc.ver_ed <- gsub(".+-", "", blasto$query_acc.ver)
blasto$subject_acc.ver_ed <- gsub(".+-", "", blasto$subject_acc.ver)

for(i in unique(blasto$query_acc.ver_ed)){
  
  if(i %in% a353IDs){
    blasto$a353[which(blasto$query_acc.ver_ed ==i)] <- "YES"
  }
}

for(i in unique(blasto$subject_acc.ver_ed)){
  
  if(i %in% schnIDs){
    blasto$schn[which(blasto$subject_acc.ver_ed == i)] <- "YES"
  }
}

## Write blast output in a table format
write.csv(blasto, file = paste0(gsub(".csv","", filename), "_ed.csv"), row.names = F)

## Check how many of the 353 genes are represented in this blast
table(a353IDs %in% unique(blasto$query_acc.ver_ed))
missing_genes <- a353IDs[which(!a353IDs %in% unique(blasto$query_acc.ver_ed))]
writeLines(missing_genes, paste(taxname, "missing_a353genes.txt", sep = "_"), sep = "\n")

## Check how many of the 275 scaffolds of the Schn kit are represented in this blast
table(schnIDs %in% unique(blasto$subject_acc.ver_ed))
missing_genes <- schnIDs[which(!schnIDs %in% unique(blasto$subject_acc.ver_ed))]
writeLines(missing_genes, paste(taxname, "missing_schnscaff.txt", sep = "_"), sep = "\n")

################################
## Extract best hits
################################

for(i in unique(blasto$query_acc.ver_ed[which(blasto$a353 == "YES")])){

 scf <- unique(blasto$subject_acc.ver_ed[which(blasto$query_acc.ver_ed == i & blasto$schn == "YES")])
 
 if(length(scf) == 1){
  write(paste(i, unique(blasto$subject_acc.ver[which(blasto$subject_acc.ver_ed == scf)]), sep = ","),
        sep = "\n", file = "gene_match.txt", append = T)
 }
}