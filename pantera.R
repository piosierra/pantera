#!/usr/bin/env Rscript

version <- "0.0.3"

## For slurm
# local({r <- getOption("repos")
# r["CRAN"] <- "https://cran.r-project.org"
# options(repos = r)
# })
# .libPaths(c("~/R/4.2.2", .libPaths()))

## Libraries
if (!require("data.table", quietly = TRUE)) {
  install.packages("data.table")
}
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

## Reading args
spec <- matrix(
  c(
    "gfas_list",          "g", 1, "character",
    "identity",           "i", 1, "integer",
    "min_cl",             "m", 1, "integer",
    "cl_size",            "c", 1, "integer",
    "output_folder",      "o", 1, "character",
    "min_size",           "s", 1, "integer",
    "max_size",           "l", 1, "integer",
    "max_pairs",          "p", 1, "integer",
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

if (is.null(opt$gfas_list)) {
  print("-gfas_list missing")
  q(status <- 1)
}

if (is.null(opt$output_folder)) {
  print("-output_folder missing")
  q(status <- 1)
}

if (is.null(opt$min_size)) { # Minimun size of TE
  opt$min_size <- 200
}

if (is.null(opt$max_size)) { # Maximum size of TE
  opt$min_size <- 20000
}

if (is.null(opt$max_pairs)) { # Maximum number of paths pairs to examine
  opt$max_pairs <- Inf
}

if (is.null(opt$identity)) { # Minimum identity for clustering
  opt$identity <- 97
}

if (is.null(opt$min_cl)) { # Minimium number of segments to create a consensus
  opt$min_cl <- 2
}

if (is.null(opt$cl_size)) { # Maximum size of cluster to process
  opt$cl_size <- 200
}

## Auxiliary function to convert from dna to binDNA format
reformatDNA <- function(dna) {
  temp <- matrix(as.character(dna),
    nrow = (length(row.names(dna))),
    dimnames = dimnames(dna)
  )
  temp <- dna(apply(temp, 1, function(x) {
    paste0(x, collapse = "")
  }))
  return(temp)
}

# @ Auxiliary function for log exits
lx <- function(x) {
  print(paste0(format(Sys.time(), "%H:%M:%S"), " [pantera ", version, "] ", x))
}

## Confirm mafft is available
mafft_exec <- system("which mafft", intern = TRUE)
if (length(grep("mafft", mafft_exec)) > 0) {
  print(paste("mafft exec:", mafft_exec))
} else {
  lx("mafft not found")
  stop()
}

segments_unique <- data.table(
  path = as.character(),
  pos = as.character(),
  seq = as.character()
)

zones <- list(
  c(200, 300), c(300, 500), c(500, 1000), c(1000, 2000), c(2000, 4000),
  c(4000, 6000), c(6000, 8000), c(8000, 13000), c(13000, 20000)
)


## Main loop
lx(paste("pantera", version))
lx(paste("Cores available:", detectCores()))
lx(paste("Gfas list:", opt$gfas_list))
lx(paste("Output:", opt$output_folder))
lx(paste("Min. size:", opt$min_size))
lx(paste("Max. size:", opt$max_size))
lx(paste("Paths pairs to examine:", opt$max_size))
lx(paste("Identity for clustering:", opt$identity))
lx(paste("Min. sequences to create a consensus:", opt$min_cl))
lx(paste("Sequences to cluster by batch:", opt$cl_size))


## Extract unique segements
path <- 0
for (g in read.table(opt$gfas_list)$V1) {
  segments <- fread(cmd = paste0("grep \"^S\" ", g), header = FALSE)
  segments <- segments[, c(1:3)]
  lx(paste("Number of segments =", nrow(segments)))
  paths <- fread(cmd = paste0("grep \"^P\" ", g), header = FALSE)
  lx(paste("Number of paths =", nrow(paths)))
  paths_s <- mclapply(
    paths$V3,
    function(x) {
      as.numeric(gsub(
        "-", "",
        gsub(
          "\\+", "",
          strsplit(x, ",")[[1]]
        )
      ))
    }
  )
  for (i in 1:(length(paths_s) - 1)) {
    for (j in (i + 1):length(paths_s)) {
      path <- path + 1
      if (path > opt$max_pairs) break
      lx(paste("Paths:", paths$V2[i], "<->", paths$V2[j]))
      if ((sum(unique(paths_s[i][[1]]) %in% unique(paths_s[j][[1]])) /
        min(length(unique(paths_s[i][[1]])), length(unique(paths_s[j][[1]]))))
      > 0.3) {
        p1 <- unique(paths_s[i][[1]])
        p2 <- unique(paths_s[j][[1]])
        pt <- unique(c(p1, p2))
        data <- data.table(seg = pt)
        data[, b1 := (seg %in% p1)]
        data[, b2 := (seg %in% p2)]
        data[, bt := b1 + b2]
        data_unique <- data[bt == 1]$seg
        segments_path <- segments[V2 %in% data_unique]
        colnames(segments_path) <- c("path", "pos", "seq")
        segments_path$path <- paths[i, 2] # Fix!!! This gives the name of the first path, it should be the one that the segment belongs to.
        lx(paste("Unique segments:", nrow(segments_path)))
        segments_unique <- rbindlist(list(segments_unique, segments_path))
      }
    }
  }
}
lx(paste("TOTAL Unique segments:", nrow(segments_unique)))
segments_unique <- segments_unique[!duplicated(segments_unique$pos)]
segments_unique[, Ns := str_count(segments_unique$seq, "N")]
segments_unique[, len := nchar(seq)]
segments_unique[, name := paste0(">", segments_unique$path, "_", segments_unique$pos)]
segments_unique <- segments_unique[len >= opt$min_size &
  len <= opt$max_size &
  Ns < len * 0.05]
segments_unique <- segments_unique[order(len)]
lx(paste("TOTAL Unique segments purged:", nrow(segments_unique)))

dir.create(opt$output_folder, showWarnings = FALSE)
setwd(opt$output_folder)

## Generate consensus from the unique segments for each "zone"
process_zone <- function(zone) {
  start <- zone[1]
  end <- zone[2]
  lx(paste("Procesing segments:", start, "-", end))
  segment_u <- segments_unique[len > start & len < end]

  ## Divide each zone in subsets if segments larger than
  segment_sets <- split(segment_u, (1:nrow(segment_u)) %/% opt$cl_size)
  lx(paste("Segment sets:", length(segment_sets)))
  for (ss in 1:length(segment_sets)) {
    sg <- segment_sets[[ss]]
    der <- strsplit(sg$seq, "")
    names(der) <- sg$name
    sol <- otu(der, k = 7, threshold = opt$identity / 100, nstart = 10)
    sol_rep <- sol[sol %in% sol[grep("\\*", names(sol), invert = T)]][grep("\\*", names(sol[sol %in% sol[grep("\\*", names(sol), invert = T)]]))]
    table(sol[sol %in% sol_rep])[table(sol[sol %in% sol_rep]) >= opt$min_cl]
    lx(paste("# of clusters:", length(sol_rep)))
    sol_list <- sol[sol %in% sol_rep]
    names(sol_list) <- gsub("\\*", "", names(sol_list))
    data_sol <- data.table(
      name = gsub(".*:", "", names(sol_list)),
      cluster = sol_list
    )
    consensi <- data.table(name = as.character(), seq = as.character())
    for (u in unique(data_sol$cluster)) {
      y <- strsplit(sg[name %in% data_sol[cluster == u]$name]$seq, "")
      names(y) <- sg[name %in% data_sol[cluster == u]$name]$name
      seqs <- as.DNAbin(y)
      ali <- mafft(seqs, thread = -1, exec = mafft_exec)
      cons <- seq_consensus(dna(reformatDNA(ali)), gaps = F)
      consensi <- rbindlist(list(
        consensi,
        data.table(
          name = data_sol[cluster == u][1]$name,
          seq = as.character(cons)
        )
      ))
    }
    fileconn <- file(paste0("consensi_", start, "_", end, "-", ss, ".fa"))
    writeLines(t(consensi), fileconn)
    close(fileconn)
  }
  lx(paste("Zone:", zone, "completed."))
  return(0)
}

mclapply(zones, process_zone)
