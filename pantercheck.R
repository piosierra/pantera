#!/usr/bin/env Rscript

options(warn = 0)


## Libraries

if (suppressPackageStartupMessages(!require("this.path", quietly = TRUE))) {
  install.packages("this.path")
}

if (suppressPackageStartupMessages(!require("data.table", quietly = TRUE))) {
  install.packages("data.table")
}


### Function to read a fasta file

ffasta <- function(f) {
  check <- file.size(f)
  if (check == 0 | is.na(check)) {
    return(data.table(name = as.character(), seq = as.character()))
  } else {
  fa_raw <- fread(f,
                  header = F,
                  fill = T,
                  sep = "\n")
  fa_raw[, h := grepl(">", V1)]
  fa_oneline <- fa_raw[, .(paste0(V1, collapse = "")), by = rleid(h)]
  return(data.table(name = fa_oneline[rep(c(TRUE, FALSE), length = .N)]$V1, seq =
                      fa_oneline[rep(c(FALSE, TRUE), length = .N)]$V1))
  }
}

### Function to write a fasta file

wfasta <- function(fasta_data, f) {
  fileconn <- file(f)
  writeLines(t(fasta_data), fileconn)
  close(fileconn)
}

f <- commandArgs(trailingOnly=TRUE)[1]
tes <- ffasta(f)

tes[,len:= nchar(seq),name]
tes[,TIR:= 0]
tes[,LTR:= 0]
tes[,polyA:= 0]
tes[,repp:= 0]


temp <- paste0("temp",make.names(gsub(".*/","",f)))
dir.create(temp)
setwd(temp)

for (i in seq(1:nrow(tes))) {
  wfasta(tes[i, c("name","seq")], paste0("tmp-",i))
  system(paste0("makeblastdb -in tmp-", i," -dbtype nucl 1> /dev/null"))
  te_data <- fread(cmd= paste0("blastn -query tmp-",i," -db tmp-",i," -outfmt 6 -word_size 11 -gapopen 5 -gapextend 2 -reward 2 -penalty -3 | cut -f 7-10"), header = F)
  colnames(te_data) <-c("qstart","qend","sstart","send")
  lente <- tes[i]$len 
  tes[i]$repp <- nrow(te_data)
  te_data[,len:=abs(qend-qstart)]
  te_data <- te_data[order(qstart)][order(-len)]
  te_data[,lgap:=.(min(qstart,qend,sstart,send)-1), by=1:nrow(te_data)]
  te_data[,rgap:=.(lente-max(qstart,qend,sstart,send)), by=1:nrow(te_data)]
  
  te_data[,check:=((qend+qstart-lente)/2)*((send+sstart-lente)/2)<0, by=1:nrow(te_data)]  ### Is each match at one side of the center?
  te_data[,check_type:=((qend-qstart)*(send-sstart))<0, by=1:nrow(te_data)] ### Is it an LTR or a TIR?

  
  if (nrow(te_data[check==T]) > 0) {
    best <- te_data[check==T][1]
    if (best$check_type) {
       tes[i]$TIR <- best$len
    } else {
      tes[i]$LTR <- best$len
    }
  } 
}
setwd("..")
unlink(temp, recursive = TRUE)
tes[,name:=gsub(">","",name)]
fwrite(tes[,c("name", "len", "LTR","TIR","polyA", "repp")], paste0(f, ".benchmark2"), col.names = F,quote =  F, row.names = F, sep =" ")


