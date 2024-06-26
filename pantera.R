#!/usr/bin/env Rscript

version <- "0.2.2"

# Uncomment and change the path to a folder to install needed libraries if you 
# don't have access to the default one (e.g. in a cluster)
# .libPaths(c("pathtowritablefolder", .libPaths())) 

options(warn=0)


## Libraries
if (!require("kmer", quietly = TRUE)) {
  install.packages("kmer")
}
if (!require("getopt", quietly = TRUE)) {
  install.packages("getopt")
}
if (!require("parallel", quietly = TRUE)) {
  install.packages("parallel")
}
if (!require("ips", quietly = TRUE)) {
  install.packages("ips")
}
if (!require("stringr", quietly = TRUE)) {
  install.packages("stringr")
}
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!require("Biostrings", quietly = TRUE)) {
  BiocManager::install("Biostrings")
}
if (!require("bioseq", quietly = TRUE)) {
  BiocManager::install("bioseq")
}
if (!require("CellaRepertorium", quietly = TRUE)) {
  BiocManager::install("CellaRepertorium")
}
if (!require("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!require("seqinr", quietly = TRUE)) {
  install.packages("seqinr")
}
if (!require("purrr", quietly = TRUE)) {
  install.packages("purrr")
}
if (!require("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
if (!require("data.table", quietly = TRUE)) {
  install.packages("data.table")
}

## Read args
spec <- matrix(
  c(
    "gfas_list",          "g", 1, "character",
    "output_folder",      "o", 1, "character",
    "cores",              "c", 1, "integer",
    "min_size",           "s", 1, "integer",
    "max_size",           "l", 1, "integer",
    "max_pairs",          "p", 1, "integer",
    "identity",           "i", 1, "double",
    "identity2",          "y", 1, "double",
    "min_cl",             "m", 1, "integer",
    "segments_shared",    "e", 1, "integer",
    "Ns",                 "n", 1, "integer",
    "k",                  "k", 1, "integer",
    "cl_size",            "u", 1, "integer",
    "chrs",               "r", 1, "integer",
    "paths_quantile",     "q", 1, "integer",
    "min_path_len",       "a", 1, "integer",
    "min_edge",           "z", 1, "integer",
    "rep_factor",         "f", 1, "integer",
    "dual",               "d", 0, "logical",
    "keep_temp",          "t", 0, "loginal",
    "help",               "h", 0, "logical"
  ),
  byrow = TRUE,
  ncol = 4
)
opt <- getopt(spec)

if (!is.null(opt$help)) {
  cat(getopt(spec, usage = TRUE))
  q(status = 1)
}

if (is.null(opt$dual)) {
  dual <- FALSE
} else
{ 
  dual <- TRUE
}

if (is.null(opt$keep_temp)) {
  keep_temp <- FALSE
} else
{ 
  keep_temp <- TRUE
}

if (is.null(opt$gfas_list)) {
  print("-gfas_list missing")
  q(status <- 1)
}

if (is.null(opt$output_folder)) {
  opt$output_folder <- "pantera_output"
}


if (is.null(opt$cores)) { # Cores
  opt$cores <- detectCores()
} else {
  opt$cores <- min(opt$cores,detectCores())
}

if (is.null(opt$min_size)) { # Minimum size of TE
  opt$min_size <- 250
}

if (is.null(opt$max_size)) { # Maximum size of TE
  opt$max_size <- 20000
}

if (is.null(opt$max_pairs)) { # Maximum number of paths pairs to examine
  opt$max_pairs <- Inf
}

if (is.null(opt$identity)) { # Minimum identity for clustering
  opt$identity <- 0.95
}

if (is.null(opt$identity2)) { # Minimum identity for clustering
  opt$identity2 <- 0.8
}

if (is.null(opt$min_cl)) { # Minimum number of segments to create a consensus
  opt$min_cl <- 2 
}

if (is.null(opt$segments_shared)) { # Minimium number of segments shared by paths
  opt$segments_shared <- 30
}

if (is.null(opt$Ns)) { # Maximum percentage of Ns in segment.
  opt$Ns <- 0
}

if (is.null(opt$k)) { # kmer size for clustering.
  opt$k <- 7 
}

if (is.null(opt$cl_size)) { # Maximum size of cluster to process
  opt$cl_size <- 1000
}

if (is.null(opt$chrs)) { # Number of primary paths to use (these are compared against all)
  opt$chrs <- Inf
}

if (is.null(opt$paths_quantile)) { # Percentage of paths to use
  opt$paths_quantile <- 100
}

if (is.null(opt$min_path_len)) { # Min path len (in nodes)
  opt$min_path_len <- 10
}

if (is.null(opt$min_edge)) { # Min leading edge (bases)
  opt$min_edge <- 3
}

if (is.null(opt$rep_factor)) { # Max factor of repeated kmers
  opt$rep_factor <- 44
}


## Auxiliary function to convert from dna to binDNA format
reformatDNA <- function(dna) {
  temp <- matrix(as.character(dna),
    nrow = (length(row.names(dna))),
    dimnames = dimnames(dna)
  )
  temp <- apply(temp, 1, function(x) {
    paste0(x, collapse = "")
  })
  return(temp)
}

## Create DNAbin from data.frame
makeDNAbin <- function(df) {
  y <- t(sapply(strsplit(df$seq,""), tolower))
  rownames(y)<- df$name
  return(as.DNAbin(y))
}

## Auxiliary function for log exits
lx <- function(x) {
  cat(paste0(format(Sys.time(), "%y-%m-%d:%H:%M:%S"), " \033[95;1;1m[pantera ",
             version, "]\033[0m ", x, "\n"))
}

## Read gfa files
get_segments <- function(segments_unique) {
  ## Extract unique segments
  pair <- 0 # counter for opt$max_pairs
  if (grepl(".gfa$", opt$gfas_list)) {
    g_list <- opt$gfas_list
  } else {
    g_list <- read.table(opt$gfas_list)$V1
  }
  for (g in g_list) {
    gc()
    lx(paste("Procesing file =", g))
    segments <- fread(cmd = paste0("grep \"^S\" ", g), header = FALSE, fill = T)
    segments <- segments[, c(1:3)]
    colnames(segments) <- c("tag","seg","seq")
    lx(paste("Number of segments =", nrow(segments)))
    segments[, len := nchar(seq)]
    links <- fread(cmd = paste0("grep \"^L\" ", g), header = FALSE, fill = T)
    lstart <-links[,.N,by=V2]
    lend <-links[,.N,by=V4]
    l11 <- merge(lstart[N==1],lend[N==1],by.x="V2",by.y="V4")
    segments <- segments[len>=opt$min_size & len<=opt$max_size]
    segments <- segments[seg %in% l11$V2]
    lx(paste("Number of valid size segments =", nrow(segments)))
    segments[, Ns := str_count(segments$seq, "N")]
    segments <- segments[Ns <= (len * opt$Ns)] 
    lx(paste("Number of valid Ns segments =", nrow(segments)))
    segments[, k5:= 1024-unlist(mclapply(str_split(segments$seq,""), function(x) {apply(kcount(x, residues = "DNA"),1,function(y){sum((y==0))})}))]
    segments[, check:=k5/min(len,1024), by=rownames(segments)]
    segments <- segments[check>(opt$rep_factor /100)]
    lx(paste("Number of valid non repetitive segments =", nrow(segments)))
    seg_candidates <- data.table(name = paste0(">", segments$seg), seq=segments$seq)
      segments_unique <- rbindlist(list(segments_unique,seg_candidates))
      rm(seg_candidates)
      rm(segments)
      rm(l11)
      rm(lend)
      rm(lstart)
      rm(links) 
  }

  
        

return(segments_unique)
}


## Generate consensus from the unique segments for each "zone"
process_zone <- function(zone) {
  start <- as.numeric(zone[1]) - 50
  end <- as.numeric(zone[2]) + 50
  #  identity <- max(0,(opt$max_size-end)%/%2000)
  if (start != end) { # If they have the same value they were processed in another zone
    lx(paste("Procesing segments:", start, "-", end))
    segment_u <- segments_unique[len >= start][len <= end]
    if (nrow(segment_u)>0) {
    ## Divide each zone in subsets if segments larger than
    # segment_sets <- split(segment_u, (1:nrow(segment_u)) %/% opt$cl_size)
    segment_sets <- split(segment_u, 
                          rep(1:(nrow(segment_u) %/% opt$cl_size +1), 
                              length.out = nrow(segment_u), 
                              each = ceiling(nrow(segment_u)/(nrow(segment_u) %/% opt$cl_size +1))))
    lx(paste("Segment sets", start, "-", end, ":", length(segment_sets)))
    for (ss in 1:length(segment_sets)) {
      sg <- segment_sets[[ss]]
      #   sg_seq <- lapply(sg$seq, function(x) {paste0(x, reverseComplement(DNAString(x)))})
      sg_seq <- unlist(lapply(sg$seq, function(x) {gsub("[^A|^T|^G|^C]","N",x)}))
      seqs_splitted <- strsplit(sg_seq, "")
      names(seqs_splitted) <- sg$name
      #    a1 <- makeDNAbin(sg[,c("name","seq")])
      #    seqs_splitted <- cbind(a1,ape::complement(a1))
      # sol <- otu(seqs_splitted, 
      #            k = opt$k, 
      #            threshold = opt$identity / 100, 
      #            nstart = 20, 
      #            method = "farthest"
      #            )
      if (round == 1) {
        clusters <- cdhit(DNAStringSet(sg$seq), identity = opt$identity, G = 0, g =1, b =500, aS = 0.95) # ,min_length = min(sg$len*0.9))
        } else {
        clusters <- cdhit(DNAStringSet(sg$seq), identity = opt$identity, g =1) #, min_length = min(sg$len*0.9))
      }
      
      
      
      rm(seqs_splitted)
      rm(sg_seq)
      # sol_list <- sol[sol %in% names(table(sol)[table(sol) >= opt$min_cl])]
      # lx(paste("# of clusters ", start, "-", end, "-", ss, ":", 
      #          length(unique(sol_list))))
      # names(sol_list) <- gsub("\\*", "", names(sol_list))
      # data_sol <- data.table(
      #   name = names(sol_list),
      #   cluster = sol_list
      # )
      consensi <- data.table(name = as.character(), seq = as.character())
      if (round > 1) {
      clusters <- clusters[clusters$n_cluster>=opt$min_cl,]  
      for (u in unique(clusters$cluster_idx)) {
        seqs_clust <- clusters[clusters$cluster_idx == u,]$seq
        minlen <- min(nchar(seqs_clust))
        clust_temp <- strsplit(seqs_clust, "")
        names(clust_temp) <- rownames(clusters[clusters$cluster_idx == u,])
        seqs <- as.DNAbin(clust_temp)
        if (length(seqs)>1) {
          ali <- mafft(seqs, 
                       #     method = "globalpair", 
                       #     maxiterate = 2, 
                       options = c("--adjustdirection"),
                       ep = 0.123, 
                       #  thread = 1, 
                       exec = mafft_exec)
          cons <- toupper(consensus(as.matrix(reformatDNA(ali)), threshold = cons_threshold))
        } else {
          cons <- toupper(paste0(unlist(as.character(seqs)),collapse=""))
        }
        cons <- gsub("-", "", cons)
        consensi <- rbindlist(list(consensi,data.table(
          name = paste0(">CONS-",start, "-", end,"-",round,"-",
                        length(clust_temp), "-",nchar(cons),"-",u),      
          seq = cons
        )))
        rm(clust_temp)
        rm(cons)
  #      rm(ali)
        rm(seqs)
      } 
      } else {
        consensi <- rbindlist(list(consensi,sg[clusters$n_cluster>=opt$min_cl,c("name","seq")]))
        
      }
      rm(clusters) 
      if (nrow(consensi)>0) {
      fileconn <- file(paste0("consensi_", start, "_", end, "-", ss, ".fa"))
      writeLines(t(consensi), fileconn)
      close(fileconn)
      }
    }
    }
  }
  lx(paste("Zone:", start, "completed."))
  rm(segment_u)
  rm(segment_sets)
  return(0)
}


### cd-hit-est fix
# globalVariables('cluster_idx')
cdhitestC <- function(opts, name, showProgress) {
  .Call('_CellaRepertorium_cdhitestC', PACKAGE = 'CellaRepertorium', opts, name, showProgress)
}
cdhit = function(seqs, identity = NULL, kmerSize = NULL, min_length = 200, 
            #     s = 1, 
                 G = 1,
                 only_index = FALSE, showProgress = interactive(), ...) {
  if(any(width(seqs) < min_length)) stop("Some sequences shorter than `min_length`;
                                           remove these or decrease min_length")
  name = 'CD-Hit'
  uopts = list(...)
  options = list()
  options$i <- tempfile()
#  options$s = s
  writeXStringSet(seqs, options$i)
  on.exit(unlink(options$i))
  type = switch(
    class(seqs),
    AAStringSet = 'cdhitC',
    DNAStringSet = 'cdhitestC',
    stop('seqs must be either AAStringSet or DNAStringSet')
  )
  if(type == 'cdhitestC'){ #DNA
    kmerSize = case_when(identity < .8 ~ 4, identity < .85 ~ 5,
                         identity < .88 ~ 6, identity < .9 ~ 7,
                         identity < .95 ~ 9, identity < 1 ~ 10, TRUE ~ 11)
    #    options = c(options, list(ap = 1, r = 0))
  } else{
    kmerSize = 5
  }
  options$n = kmerSize
  options$G = G
  options$c = identity
  options$l = min_length - 1
  options = c(uopts, options)
  options = options[!duplicated(names(options))]
  options = lapply(options, as.character)
#  print(options)
  if(type == 'cdhitC'){
    seq_cluster_index = cdhitC(options, name, showProgress) + 1
  } else{
    seq_cluster_index = cdhitestC(options, name, showProgress) + 1
  }
  if(only_index) return(seq_cluster_index)
  tibble::tibble(query_name = names(seqs), seq = as.character(seqs),
                 cluster_idx = seq_cluster_index) %>%
    dplyr::group_by(cluster_idx) %>%
    dplyr::mutate(n_cluster = dplyr::n()) %>% ungroup()
}



### MAIN

dir.create(opt$output_folder, showWarnings = FALSE)
sink(paste0(opt$output_folder,"/pantera.log"))

## DF to store the unique segments
segments_unique <- data.table(
  name = as.character(),
  seq = as.character()
)

## Main loop
lx(paste("pantera", version,"\n"))
## Confirm mafft is available
mafft_exec <- system("which mafft", intern = TRUE)
if (length(grep("mafft", mafft_exec)) > 0) {
  lx(paste("mafft exec:", mafft_exec))
} else {
  lx("mafft not found")
  stop()
}
lx(paste("Gfas list:", opt$gfas_list))
lx(paste("Output:", opt$output_folder))
lx(paste("Cores to use:", opt$cores))
lx(paste("Cores detected:", detectCores()))
lx(paste("Min. size:", opt$min_size))
lx(paste("Max. size:", opt$max_size))
lx(paste("Identity for clustering:", opt$identity))
lx(paste("Min. sequences to create a consensus:", opt$min_cl))
lx(paste("Max percentage of Ns:", opt$Ns))
lx(paste("kmer size for clustering:", opt$k))
lx(paste("Max. size of cluster:", opt$cl_size))
lx(paste("Number of paths to use as primary:", opt$chrs))
lx(paste("Minimum number of nodes in path:", opt$min_path_len))
lx(paste("Minimum number of bases at edge of path:", opt$edge))
lx(paste("Fast mode:", dual))
if (!dual) {
lx(paste("Paths pairs to examine:", opt$max_pairs))
lx(paste("Minimum percentaje of shared segments:", opt$segments_shared))
lx(paste("Percentage of paths to use:", opt$paths_quantile))
}

if (!file.exists(paste0(opt$output_folder,"/all_segments.fa"))) {
  segments_unique <- get_segments(segments_unique)
  setwd(opt$output_folder)
  fileConn <- file("all_segments.fa")
  writeLines(t(segments_unique[, c("name", "seq")]), fileConn)
  close(fileConn)
} else {
  setwd(opt$output_folder)
}

round <- 0
prev <- Inf
cons_threshold <- 0.6
while (TRUE) {
  round <- round + 1
  lx(paste("Reading data for round:", round))
  segu_names <- fread(cmd = "grep '>' all_segments.fa", header = F, sep = "\t")
  segu_seqs <- fread(cmd = "grep -v '>' all_segments.fa", header = F, sep = "\t")
  segments_unique <- data.table(name = segu_names$V1, seq = segu_seqs$V1)
  lx(paste("TOTAL Unique segments:", nrow(segments_unique)))
  segments_unique <- segments_unique[!duplicated(segments_unique)]
  segments_unique[, len := nchar(seq)]
  segments_unique <- segments_unique[order(len)]
  lx(paste("TOTAL Unique segments no dups:", nrow(segments_unique)))
  segments_unique[,len :=  nchar(seq)]
  segments_unique <- segments_unique[order(-len)]
  if (nrow(segments_unique)/prev > 0.98) { 
    system("mv all_segments.fa pantera_lib.fa")
    lx(paste("End of process"))
    setwd("..")
    break 
  } else {
    prev <- nrow(segments_unique)
  }
  if (round > 1) {
 #   segments_unique[,times:=as.numeric(gsub("-.*","",gsub(">CONS-[0-9]+-[0-9]+-[0-9]+-","",name)))]
 #   segments_unique <- segments_unique[rep(1:.N,times)][,Indx:=1:.N,by=name]
    opt$identity <- opt$identity2
    opt$min_cl <- 0
    opt$cl_size <- 3000
    cons_threshold <- 0.6
    opt$Ns <- 0
    segments_unique[,name:= paste0(name,"-",.I)]
#    lx(paste("Expanded sequences:", nrow(segments_unique)))
  } 
  lx(paste("Starting loop:", round))
  lx(paste("Processing:", nrow(segments_unique), "segments"))
  lx(paste("Largest segment:", max(segments_unique$len)))
  lx(paste("Smallest segment:", min(segments_unique$len)))
  dir.create(paste0("loop_",round), showWarnings = FALSE)
  setwd(paste0("loop_",round))
  # z <- segments_unique[seq(1,nrow(segments_unique),min(nrow(segments_unique)%/%1,opt$cl_size))]$len
  # if (length(z) < 2) {
  #   zones <- list(data.frame(start = segments_unique[nrow(segments_unique)]$len, end = segments_unique[1]$len))
  # } else {
  #   zones <- asplit(data.frame(start = z[-1], end = z[-length(z)]),1)
  # }
  if (round ==1) {
    z <- seq(min(segments_unique$len), max(segments_unique$len)+100, 100)
    zones <- asplit(data.frame(start = z[-length(z)], end = z[-1]),1)
  } else {
    zones <- asplit(data.frame(start = min(segments_unique$len), end = max(segments_unique$len)),1)
  }
  
  lx(paste("Processing:", length(zones), "windows"))
  gc()
  loop_exit <- mclapply(rev(zones),
                        process_zone,
                        mc.preschedule = FALSE, 
                        mc.cores = opt$cores)
system(paste0("mv ../all_segments.fa ../pantera_lib_", round - 1, ".fa"))
system("cat consen*.fa > ../all_segments.fa")
setwd("..")
gc()
}




