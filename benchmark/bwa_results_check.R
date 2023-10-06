library(data.table)
library(ggplot2)
library(stringr)
setwd("~/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/pio/projects/20230313_Dmel_pangenome/")
sam <- fread(cmd= "cut -f1-6 drosophila-transposons/releases/mapped.sam")
sam$M <- lapply(sam$V6, function(x) {sum(na.omit(as.numeric(gsub(".*[A-Z]","",str_split(x,"[M]")[[1]]))))})


all_names <- fread(cmd = "grep '>' pantera_output7/all_segments.fa", header = F, sep = "\t")
all_seqs <- fread(cmd = "grep -v '>' pantera_output7/all_segments.fa", header = F, sep = "\t")
all_names[,len:=nchar(all_seqs$V1)]
all_names[,V1:=gsub(">","",V1)]
sam <- merge(sam,all_names)
sam$M <- unlist(sam$M)
sam[,per:=M/len]
all_names[,hit:=(V1 %in% sam[per>=0.98]$V1)]

points <- seq(0,18500,1000)

data <- unlist(lapply(points, function(x) {nrow(all_names[len>=x & hit==TRUE]) / nrow(all_names[len>=x]) }))
plot(points,data)

all_names[,bin:= len %/% 1000]
data2 <- all_names[,.(.N,sum(hit)),bin]
data2[,per:=V2/N]
data2$points <- points
plot(points,data2$per)

svg("segments_dmel.svg")
p <- ggplot(data = data2, aes(x = points, y =per)) +
  geom_col() +
  theme_classic() +
  xlab("Segment size (1Kb bins)") +
  ggtitle("Percentage of unique segments mapping to Dmel TE Ref library (identity >= 98)")
p
dev.off()


# segments mapping with MAPQ60 
nrow(all_names[hit==TRUE])

# segments mapping with MAPQ60 and len > 800

nrow(all_names[len>1000])
nrow(all_names[len>1000 & hit==TRUE])

length(unique(sam[V5==60]$V1))
