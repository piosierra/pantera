library(data.table)
library(ggplot2)
library(ggpubr)
library(gridExtra)

setwd("~/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/pio/projects/20230313_Dmel_pangenome/")
lib1_names <- fread(cmd = "grep '>' pantera21_output3/pantera_lib.fa.classified.sl", header = F, sep = "\t")
lib1_seqs <- fread(cmd = "grep -v '>' pantera21_output3/pantera_lib.fa.classified.sl", header = F, sep = "\t")

lib1 <- data.table(name= gsub(">","",lib1_names$V1), len = nchar(lib1_seqs$V1))
# lib1[,fam:=gsub("_.*$","",name)]
lib1[,sf:=gsub("-[0-9]*_.*$","",gsub("^.*#","",name))]
# lib1[grep("\\[LTR\\]",name),sf:="LTR"]
lib1[,order:=gsub("/.*","",sf)]
lib1[,lib:="pant"]

data2 <- "drosophila-transposons/releases/D_mel_transposon_sequence_set_v10.2.fa.sl"
data2_names <- fread(cmd=paste("cat", data2,"| grep '>'"), sep = "", header = FALSE)
data2_seqs <- fread(cmd=paste("cat", data2,"| grep -v '>'"), sep = "", header = FALSE)
lib2 <- data.table(name= gsub(">","",data2_names$V1), len = nchar(data2_seqs$V1))
# lib1[,fam:=gsub("_.*$","",name)]
lib2[,sf:=gsub("-[0-9]*_.*$","",gsub("^.*#","",name))]
# lib1[grep("\\[LTR\\]",name),sf:="LTR"]
lib2[,order:=gsub("/.*","",sf)]
lib2[,lib:="Reference"]

data3 <- "RM_93922.ThuJun221105402023/consensi.fa.classified.sl"
data3_names <- fread(cmd=paste("cat", data3,"| grep '>'"), sep = "", header = FALSE)
data3_seqs <- fread(cmd=paste("cat", data3,"| grep -v '>'"), sep = "", header = FALSE)
lib3 <- data.table(name= gsub(">","",data3_names$V1), len = nchar(data2_seqs$V1))
lib3[,name:=gsub(" .*","",name)]
# lib1[,fam:=gsub("_.*$","",name)]
lib3[,sf:=gsub("-[0-9]*_.*$","",gsub("^.*#","",name))]
# lib1[grep("\\[LTR\\]",name),sf:="LTR"]
lib3[,order:=gsub("/.*","",sf)]
lib3[,lib:="RepeatModeler"]


data_all <- rbind(lib1, lib2, lib3)

data_all <- data_all[order(-sf)][order(order)]
data_all[, fam:=factor(sf, levels = unique(data_all$sf))]
data_all[sf=="DNA/hAT-hobo", sf:="DNA/hAT"]
data_all[sf=="DNA/CMC-Transib", sf:="DNA/Transib"]
data_all[sf=="DNA/TcMar-Tc1", sf:="DNA/Tc1-Mariner"]
data_all[sf=="LTR/Pao", sf:="LTR/Bel-Pao"]

#svg("Dmelref_vs_pantera_vs_RM.svg", height = 10, width = 20)
p1 <- ggboxplot(data_all, x = "order", y = "len",
                add.params = list(size = 1.5),
                color = "lib", palette = c("#D63384FF","#33AAFFFF","#1122BB"), #,"#5B84B1FF"),
                add = c("jitter"), size = 1, shape="lib") +
#  theme(legend.position = "none") +
  font("y.text", size = 10) +
  font("title", size = 20) +
  font("x.text", size = 10) +
  font("xlab", size = 10) +
#  coord_flip() +
  xlab("") +
  ylab("basepairs") +
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  ggtitle(paste0("Differences in fragmentation of TE libraries:\nPantera library (N = ",nrow(lib1),
                  ") vs \nDmel Ref library (N = ",nrow(lib2), 
                  ") vs \nRepeatModeler generated library (N = ",nrow(lib3),")")) 
# theme_transparent(base_size = 30)
p1
#dev.off()

data1 <- "pantera16d_output_A1A2/pantera_lib.fa.classified.sl"
data1_names <- fread(cmd=paste("cat", data1,"| grep '>'"), sep = "", header = FALSE)
data1_seqs <- fread(cmd=paste("cat", data1,"| grep -v '>'"), sep = "", header = FALSE)

lib1 <- data.table(name= gsub(">","",data1_names$V1), len = nchar(data1_seqs$V1))
# lib1[,fam:=gsub("_.*$","",name)]
lib1[,sf:=gsub("-[0-9]*_.*$","",gsub("^.*#","",name))]
# lib1[grep("\\[LTR\\]",name),sf:="LTR"]
lib1[,order:=gsub("/.*","",sf)]
lib1[,lib:="pant"]
