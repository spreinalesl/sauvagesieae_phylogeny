#############################################################
## This script checks the congruence
## among Schneider et al. (2021) baits and CKDK scaffolds 
############################################################

## Path
dpath <- "~/projects/target_files/Bait_set_Ochnaceae-120-32-noRM.fas"
rpath <- "~/projects/target_files/custom_target/schn_baits/"

#############################
## Extract bait information
## Sequences names contains the following information
## Taxon ID = CKDK
## scaffold number = ex. 202626833
## Taxon name = Ochna-serrulata
## final position of the seq with respect to the scaffold = ex. 33
#############################
## Read file and extract scaffolds names
# Downloaded from https://zenodo.org/record/4278991
data <- readLines(dpath)
scf_info <- grep(">", data, value = T)

## Summarise bait inf per scaffold
scaffolds <- unlist(lapply(X = scf_info, FUN = function(x) unlist(strsplit(x, split = "-"))[3]))
dtable <- data.frame(scaffold = unique(scaffolds)) 

for(i in 1:nrow(dtable)){
  flength <- unlist(lapply(X = scf_info[grep(dtable$scaffold[i], scf_info)],
                           FUN = function(x) unlist(strsplit(x, split = "_"))[2]))
  dtable$start[i] <- min(as.numeric(flength))
  dtable$end[i] <- max(as.numeric(flength))
  
  dtable$first_base[i] <- min(as.numeric(flength))+1
  dtable$last_base[i] <- max(as.numeric(flength))+120
  dtable$targeted_length[i] <- dtable$last_base[i] - dtable$first_base[i]
}

## Check which scaffolds have match with a353 genes
mgenes <- read.table("~/projects/custom_target/blast/gene_match.txt",
                     header = F, sep = ",")
dtable$share <- NA
for(i in mgenes$V2){
  dtable$share[which(dtable$scaffold == gsub(".+-","", i))] <- "YES"
}

## Write table
write.csv(dtable, file = paste(rpath, "ochnaceae_baits_summ.csv", sep = ""), row.names = F)