#########################################
## This script read bioproject metadata
## from NCBI searching by bioproject code
## and create a table with species names,
## vouches and SRR codes
#########################################

## Define paths
fpath <- "./biosample_result_schneider.txt"

## Read metadata file
data <- readLines(fpath)

## Extract info about Organism
names <- grep("Organism", data, value = T)
names <- lapply(X = strsplit(names, split = ": "), FUN = function(x) x[[2]])
names <- gsub(" ", "_", names)
length(names)

## Extract biosample and SRA
identifier <- grep("Identifiers", data, value = T)

biosample <- lapply(X = strsplit(identifier, split = " "), FUN = function(x) x[[3]])
biosample <-  gsub(";", "", unlist(biosample))
length(biosample)

sra <- unlist(lapply(X = strsplit(identifier, split = " "), FUN = function(x) x[[8]]))
length(sra)

## Extract specimen voucher
svoucher <- grep("specimen voucher", data, value = T)  
svoucher <- strsplit(svoucher, split = "\"")
svoucher <- unlist(lapply(X = svoucher, FUN = function(x) x[[2]]))
length(svoucher)

## Create and write a table
sample_schn <- data.frame(spp = names,
                          voucher = svoucher,
                          biosample = biosample,
                          sra_code = sra,
                          probe = "Schn",
                          library_name = NA)

###########################
## For Shah 2021 Ochna
###########################

## Define paths
fpath <- "./biosample_result_shah_ochna.txt"

## Read metadata file
data <- readLines(fpath)

## Extract info about Organism
names <- grep("Organism", data, value = T)
names <- lapply(X = strsplit(names, split = "\t"), FUN = function(x) x[[2]])
names <- gsub(" ", "_", names)
length(names)

## Extract biosample and SRA
identifier <- grep("Identifiers", data, value = T)

biosample <- lapply(X = strsplit(identifier, split = " "), FUN = function(x) x[[2]])
biosample <-  gsub(";", "", unlist(biosample))
length(biosample)

sra <- unlist(lapply(X = strsplit(identifier, split = "SRA: "), FUN = function(x) x[[2]]))
length(sra)

library_name <- unlist(lapply(X = strsplit(identifier, split = " "), FUN = function(x) x[[5]]))
length(library_name)

## Create and write a table
sample_shah_ochn <- data.frame(spp = names,
                               voucher = NA,
                               biosample = biosample,
                               sra_code = sra,
                               probe = "Schn",
                               library_name = library_name)

###########################
## Filter PAFTOL table
## using a list of names
###########################

## Define paths
fpath <- "./biosample_result_PAFTOL.txt"

## Read metadata file
# For Shah 2021 joining forces
data <- readLines(fpath)

## Extract info about Organism
names <- grep("Organism", data, value = T)
names <- lapply(X = strsplit(names, split = "\t"), FUN = function(x) x[[2]])
names <- gsub(" ", "_", names)
names <- gsub("Sauvagesia_serrata", replacement = "Neckia_serrata", names)
length(names)

## Extract biosample and SRA
identifier <- grep("Identifiers", data, value = T)

biosample <- lapply(X = strsplit(identifier, split = " "), FUN = function(x) x[[2]])
biosample <-  gsub(";", "", unlist(biosample))
length(biosample)

sra <- unlist(lapply(X = strsplit(identifier, split = " "), FUN = function(x) x[[4]]))
length(sra)

## Create and write a table
sample_PAFTOL <- data.frame(spp = names,
                          voucher = NA,
                          biosample = biosample,
                          sra_code = sra,
                          probe = NA,
                          library_name = NA)

## Read voucher table Shah (2021)
vouchers <- read.csv("apendix1_and_tableS1_shah_vouchers.csv")
vouchers$name <- trimws(vouchers$name, "both")
vouchers$name <- gsub("Sauvagesia_serrata", replacement = "Neckia_serrata", vouchers$name)
vouchers$author <- trimws(vouchers$author, "both")
vouchers$country <- trimws(vouchers$country, "both")
vouchers$voucher <- trimws(vouchers$voucher, "both")

## Create a name vector with all names including in the different metadata to filter PAFTOL table
names <- unique(c(sample_schn$spp, sample_shah_ochn$spp, vouchers$name))

for(i in names){
  p <- which(i == sample_PAFTOL$spp)
  
  if(!(length(p) == 0)){
    sample_PAFTOL$probe[p] <- "a353"
  }
}

#########################################
## Merge tables and generate links to
#########################################
tsample <- rbind(sample_PAFTOL[which(sample_PAFTOL$probe == "a353"),],
                 sample_schn,
                 sample_shah_ochn)

tsample$spp <- gsub("sp\\.", replacement = "sp", x = tsample$spp)
tsample$spp <- gsub("cf\\.", replacement = "cf", x = tsample$spp)
tsample$spp <- gsub("var\\.", replacement = "var", x = tsample$spp)

#########################################
## Download metadata from https://www.ebi.ac.uk/
## for the following projects
## PRJEB42299: Shah - PAFTOL
## PRJNA602196: Schn
## PRJNA735954: Shah - Ochna
#########################################
## Merge datasets
system("cat filereport_read_run* > all_run.txt")

## Read file
alldata <- read.csv("all_run.txt", sep = "\t")

## Standardize names
alldata$scientific_name <- gsub("sp\\. ", replacement = "sp_", x = alldata$scientific_name)
alldata$scientific_name <- gsub("cf\\. ", replacement = "cf_", x = alldata$scientific_name)
alldata$scientific_name <- gsub("var\\. ", replacement = "var_", x = alldata$scientific_name)
alldata$scientific_name <- gsub(" ", replacement = "_", x = alldata$scientific_name)
alldata$scientific_name <- gsub("Sauvagesia_serrata", replacement = "Neckia_serrata", x = alldata$scientific_name)
alldata <- alldata[-which(alldata$scientific_name == "scientific_name"),]
alldata$genus <- unlist(lapply(strsplit(alldata$scientific_name, split = "_"), FUN = function(x) x[[1]])) 
  
## Filter just Ochnaceae records
genera <- unique(unlist(lapply(strsplit(names, split = "_"), FUN = function(x) x[[1]])))
alldata <- alldata[which(alldata$genus %in% genera), ]

#########################################
## Write edited table
## Another manual editions will be doing
## in this document
#########################################
write.csv(alldata, "all_run_ed.txt", row.names = F)

#########################################
## Read edited table alldata to complete
## voucher names using the other tables
#########################################
alldata_ed <- read.csv("all_run_edit.txt", sep = ",")
alldata_ed$verif <- NA
alldata_ed$probes <- NA

### Search for voucher info in Schneider dataset
for(i in sample_schn$sra_code){
  l <- apply(alldata_ed, MARGIN = 2, FUN = function(x) grep(i, x))
  p <- as.vector(unlist(l)) 
  
  
  if(unique(p > 0)){

    alldata_ed$verif[p] <- "OK"
    alldata_ed$probes[p] <- "Schn"
      
    if(unique(is.na(alldata_ed$voucher[p]))){
      biosample <- unique(alldata_ed$sample_accession[p])
      p3 <- which(alldata_ed$sample_accession == biosample)
      alldata_ed$voucher[p3] <- sample_schn$voucher[which(sample_schn$sra_code == i)]
    }
  }
}

### Shah tables
for(i in vouchers$library_name){
  l <- apply(alldata_ed, MARGIN = 2, FUN = function(x) grep(i, x))
  p <- as.vector(unlist(l)) 
  
  
  if(length(p) > 0){
    
    alldata_ed$verif[p] <- "OK"
    voucher <- alldata_ed$voucher[p]
    
    if(length(unique(voucher)) == 1){
      if(is.na(unique(voucher))){
        p2 <- alldata_ed$sample_accession[which(alldata_ed$library_name == i)]
        p2 <- which(alldata_ed$sample_accession == p2)
        alldata_ed$voucher[p2] <- vouchers$voucher[which(vouchers$library_name == i)]
        alldata_ed$probes[p2] <- "Schn"
      }
    }
  }
}

#########################################
## write table
#########################################
write.csv(alldata_ed, "all_run_edit.txt", row.names = F)
