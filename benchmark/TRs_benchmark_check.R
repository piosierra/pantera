library(data.table)
setwd("~/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/pio/projects/20230313_Dmel_pangenome/")

data <- fread("pant112_all/pantera_lib.fa.classified.sl.benchmark2")
data[,sf:=gsub(".*#","",V1)]
data[,fam:=gsub(".*\\/","",sf)]
data[,order:=gsub("\\/.*","",sf)]
data$check <- FALSE
data[order=="DNA" & V4>0, check :=TRUE]
data[order=="LTR" & V3>0, check :=TRUE]
data[order=="LINE" & V5>0, check :=TRUE]
table(data[order=="LTR"]$check)
table(data[order=="DNA"]$check)
table(data[order=="LINE"]$check)

data <- fread("RM_93922.ThuJun221105402023/consensi.fa.classified.sl.benchmark2")
data[,sf:=gsub(".*#","",V1)]
data[,fam:=gsub(".*\\/","",sf)]
data[,order:=gsub("\\/.*","",sf)]
data$check <- FALSE
data[order=="DNA" & V4>0, check :=TRUE]
data[order=="LTR" & V3>0, check :=TRUE]
data[order=="LINE" & V5>0, check :=TRUE]
data[sf=="DIRS" & V4==0, check := FALSE]
data[sf=="DIRS" & V4>0, check := TRUE]
table(data[order=="LTR"]$check)
table(data[order=="DNA"]$check)
table(data[order=="LINE"]$check)

data <- fread("pantera20_output/pantera_lib.fa.classified.sl.benchmark2")
data[,sf:=gsub(".*#","",V1)]
data[,fam:=gsub(".*\\/","",sf)]
data[,order:=gsub("\\/.*","",sf)]
data$check <- FALSE
data[order=="DNA" & V4>0, check :=TRUE]
data[order=="LTR" & V3>0, check :=TRUE]
data[order=="LINE" & V5>0, check :=TRUE]
data[sf=="DIRS" & V4==0, check := FALSE]
data[sf=="DIRS" & V4>0, check := TRUE]
table(data[order=="LTR"]$check)
table(data[order=="DNA"]$check)
table(data[order=="LINE"]$check)


data <- fread("drosophila-transposons/releases/D_mel_transposon_sequence_set_v10.2.fa.benchmark2")
data[,sf:=gsub(".*#","",V1)]
data[,fam:=gsub(".*\\/","",sf)]
data[,order:=gsub("\\/.*","",sf)]
data$check <- FALSE
data[order=="DNA" & V4>0, check :=TRUE]
data[order=="LTR" & V3>0, check :=TRUE]
data[order=="LINE" & V5>0, check :=TRUE]
data[sf=="DIRS" & V4==0, check := FALSE]
data[sf=="DIRS" & V4>0, check := TRUE]
table(data[order=="LTR"]$check)
table(data[order=="DNA"]$check)
table(data[order=="LINE"]$check)

### ALL SEGMENTS
data <- fread("pantera21_output3/pantera_lib_0.fa.benchmark2")
### LTR
nrow(data[V3>0])
### TIR
nrow(data[V4>0])
### polyA
nrow(data[V5>0])

## NA DNA (palindromes)
print("Palindromes")
nrow(data[V4 > (V2*0.9)])

data[,range:= 1000 *(V2 %/% 1000)]
data_plot <- data[,.(LTR=sum(V3>0),TIR=sum(V4>0),pA=sum(V5>0), none =sum(V3==0 & V4==0 & V5==0) ), by=range]
data_plot <- melt(data_plot, id.vars = c("range"), measure.vars = c("LTR","TIR","pA", "none"))
ggplot(data_plot, aes(range,value, fill = variable)) +
  geom_col()

