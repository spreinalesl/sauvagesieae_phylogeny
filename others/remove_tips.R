library(ape)

probename <- "a353"
type <- "exons"
addf <- "rmout_infprlgs"

fpath <- paste("~/projects/rbgk/users_area/sreinales/results/sauvagesieae/hybpiper/fasta_files",
                probename, type, addf, sep = "/")

# output files
outpath <- paste(fpath, "edfast", sep = "/")

flist <- list.files(fpath, pattern = "FNA")

for(i in 1:length(flist)){
  
  gname <- gsub(pattern = "_ed.FNA", flist[i], replacement = "")
  
  ff <- read.FASTA(paste(fpath, flist[i], sep = "/"))
  
  tips <- gsub(labels(ff), pattern = " .+", replacement = "")
  
  pst <- c(which(tips == "SPRL_291"), which(tips == "SPRL_287"), which(tips == "SPRL_196"))

  if(length(pst) > 0){
    ff <- ff[-pst]
  }
  
  write.dna(ff, format = "fasta", nbcol = 1, colsep = "", colw = 60, indent = 0, blocksep = 0,
            file = paste(outpath, "/", gname, ".FNA", sep = ""))
  
  cat("Writing fasta file", flist[i], "\n")
}
