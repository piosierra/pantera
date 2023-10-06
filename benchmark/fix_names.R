#!/usr/bin/env Rscript

.libPaths(c("~/R/4.2.2", .libPaths()))

library(data.table)
args <- commandArgs(trailingOnly = TRUE)

ypolib_names <- fread(cmd = eval(paste0("grep '>' ",args[1])), header = F, sep = "\t")
ypolib_seqs <- fread(cmd = eval(paste0("grep -v '>' ",args[1])), header = F, sep = "\t")
ypolib_names[, name:=paste0(substr(gsub(">###",">",V1),1,9),"_",.I)]

fileConn <- file(paste0(args[1],".fixed"))
writeLines(t(data.frame(ypolib_names$name,ypolib_seqs)), fileConn)
close(fileConn)