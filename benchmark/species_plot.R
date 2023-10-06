library("data.table")
library("ggplot2")
library("ggpubr")

species <- c("Drosophila melanogaster", "Oryza sativa", "Danio rerio")
species_dir <- c("Dmel", "Osav", "Drer")
libs <- c("pantera","Reference","Other")
lib_dir <- c("pantera","ref","alt")

data_summ <- data.table(LTR = as.numeric(),
                        LTR_OK = as.numeric(),
                        TIR = as.numeric(),
                        TIR_OK = as.numeric(),
                        LINE = as.numeric(),
                        LINE_OK = as.numeric(),
                        lib = as.character(),
                        species = as.character())

for (j in 1:3) {
for (i in 1:3) {
  data <- fread(cmd = paste0("cat data/",species_dir[j],"/",lib_dir[i],"/*.benchmark2"))
  data[,sf:=gsub(".*#","",V1)]
  data[,fam:=gsub(".*\\/","",sf)]
  data[,order:=gsub("\\/.*","",sf)]
  data$check <- FALSE
  data[order=="DNA" & sf!= "Crypton" & V4>0, check :=TRUE]
  data[order=="LTR" & V3>0, check :=TRUE]
  data[sf=="DIRS" & V4==0, check := FALSE]
  data[sf=="DIRS" & V4>0, check := TRUE]
  line <- fread(cmd = paste0("cat data/",species_dir,"/",lib_dir[i],"/complete_lines"))
  line[is.na(line)] <- 0
  data_summ <- rbindlist(list(data_summ, 
                              data.table(LTR = sum(data$order=="LTR"),
                                         LTR_OK = sum(data$order=="LTR" & data$sf!="DIRS" & data$V3>0) + sum(data$sf=="DIRS" & data$V4>0),
                                         TIR = sum(data$order=="DNA" & data$sf!= "Crypton" ),
                                         TIR_OK = sum(data$order=="DNA" & data$sf!= "Crypton" & data$V4>0),
                                         LINE = sum(data$order=="LINE"),
                                         LINE_OK = unlist(line[2]),
                                         lib = libs[i],
                                         species = species[j])))
  
}
}

data_summ[,LTR:=LTR-LTR_OK]
data_summ[,TIR:=TIR-TIR_OK]
data_summ[,LINE:=LINE-LINE_OK]

data_plot <- melt(data_summ, id.vars = c("lib","species"),measure.vars = c("LTR","LTR_OK", "TIR", "TIR_OK", "LINE", "LINE_OK") )
data_plot <- data_plot[,order:=gsub("_.*","",variable)]
data_plot <- data_plot[,variable:=gsub(".*_","",variable)]
data_plot <- data_plot[variable!="OK",variable:="Not OK"]

data_plot[order=="TIR",order:="DNA"]

data_plot$order <- factor(data_plot$order, levels = c("DNA","LINE","LTR"))
data_plot$species<- factor(data_plot$species, levels = c("Drosophila melanogaster", "Oryza sativa", "Danio rerio"))
data_plot$lib<- factor(data_plot$lib, levels = c("pantera", "Reference", "Other"))  

                    
p1 <- ggbarplot(data_plot, x = "lib", y = "value", fill = "variable",
                add.params = list(size = 1.5), color = "#FFFFFFFF", scales="free",
                palette = c("#BFBFBFFF","#00CC00FF")) +
  #  theme(legend.position = "none") +
  font("y.text", size = 20) +
  font("title", size = 20) +
  font("x.text", size = 20) +
  font("xlab", size = 10) +
#  coord_flip() +
  xlab("") +
  ylab("basepairs") +
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  ggtitle(paste0("Benchmark"))

p <- facet(p1,facet.by = c("species","order"), scales = "free")
ggsave("plots/benchmark1.svg",p, device = "svg", height =10, width = 6)

data_summ[, c("LTR", "LTR_OK"):= list(LTR/(LTR+LTR_OK), LTR_OK/(LTR+LTR_OK))]
data_summ[, c("TIR", "TIR_OK"):= list(TIR/(TIR+TIR_OK), TIR_OK/(TIR+TIR_OK))]
data_summ[, c("LINE", "LINE_OK"):= list(LINE/(LINE+LINE_OK), LINE_OK/(LINE+LINE_OK))]

data_plot <- melt(data_summ, id.vars = c("lib","species"),measure.vars = c("LTR","LTR_OK", "TIR", "TIR_OK", "LINE", "LINE_OK") )
data_plot <- data_plot[,order:=gsub("_.*","",variable)]
data_plot <- data_plot[,variable:=gsub(".*_","",variable)]
data_plot <- data_plot[variable!="OK",variable:="Not OK"]

data_plot[order=="TIR",order:="DNA"]

data_plot$order <- factor(data_plot$order, levels = c("DNA","LINE","LTR"))
data_plot$species<- factor(data_plot$species, levels = c("Drosophila melanogaster", "Oryza sativa", "Danio rerio"))
data_plot$lib<- factor(data_plot$lib, levels = c("pantera", "Reference", "Other"))  

p1 <- ggbarplot(data_plot, x = "lib", y = "value", fill = "variable",
                add.params = list(size = 1.5), color = "#FFFFFFFF",
                palette = c("#FFFFFFFF","#00CC00FF")) +
  #  theme(legend.position = "none") +
  font("y.text", size = 20) +
  font("title", size = 20) +
  font("x.text", size = 20) +
  font("xlab", size = 10) +
 # coord_flip() +
  xlab("") +
  scale_y_continuous(labels = scales::percent) +
  ylab("basepairs") +
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  ggtitle(paste0("Benchmark"))

p <- facet(p1,facet.by = c("species","order"), scales = "free")
ggsave("plots/benchmark2.svg",p, device = "svg", height =10, width = 6)

data_cover <- data.table(LTR = as.numeric(),
                        DNA = as.numeric(),
                        LINE = as.numeric(),
                        RC = as.numeric(),
                        Unclassified = as.numeric(),
                        lib = as.character(),
                        species = as.character())

for (j in 1:3) {
  for (i in 1:3) {
    data <- fread(cmd = paste0("cat data/",species_dir[j],"/",lib_dir[i],"/*.tbl | tail -n+11 | head -n -11 | awk -F' ' '{print $5 $6}'"), header = F)
    data[,V1:=as.numeric(gsub("bp","",gsub("%","",V1)))]
    data_cover <- rbindlist(list(data_cover, data.table(LTR = data[11]$V1,
                                                       DNA = data[17]$V1,
                                                       LINE = data[1]$V1-data[11]$V1,
                                                       RC = data[27]$V1,
                                                       Unclassified = data[29]$V1,
                                                       lib = libs[i],
                                                       species = species[j]
                                                       )))
  }
}

data_plot <- melt(data_cover, id.vars = c("lib","species"),measure.vars = c("DNA", "LINE", "LTR","RC", "Unclassified") )

data_plot$species<- factor(data_plot$species, levels = c("Drosophila melanogaster", "Oryza sativa", "Danio rerio"))
data_plot$lib<- factor(data_plot$lib, levels = c("pantera", "Reference", "Other"))  
data_plot$value <- data_plot$value / 100
data_plot$variable <- factor(data_plot$variable, levels = c("DNA","LINE","LTR", "RC", "Unclassified"))
p1 <- ggbarplot(data_plot, x = "variable", y = "value", fill = "variable", color ="#FFFFFFFF",
                add.params = list(size = 1.5), palette = c("#FF6600FF","#3366ccFF","#669900FF","#FFCC00FF","#bfbfbfff")) +
  #  theme(legend.position = "none") +
  font("y.text", size = 20) +
  font("title", size = 20) +
  font("x.text", size = 20) +
  font("xlab", size = 10) +
  # coord_flip() +
  xlab("") +
  scale_y_continuous(labels = scales::percent, limits = c(0, .45)) +
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  ggtitle(paste0("Benchmark"))

p <- facet(p1,facet.by = c("species","lib"), scales = "free")
ggsave("plots/cover1.svg",p, device = "svg", height =10, width = 6)


p1 <- ggbarplot(data_plot, x = "lib", y = "value", fill = "lib", color ="#FFFFFFFF",
                add.params = list(size = 1.5), palette = c("#d63384FF","#689624FF","#404040FF")) +
  #  theme(legend.position = "none") +
  font("y.text", size = 20) +
  font("title", size = 20) +
  font("x.text", size = 20) +
  font("xlab", size = 10) +
  # coord_flip() +
  xlab("") +
  scale_y_continuous(labels = scales::percent, limits = c(0, .45)) +
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  ggtitle(paste0("Benchmark"))

p <- facet(p1,facet.by = c("species","variable"), scales = "free")
ggsave("plots/cover2.svg",p, device = "svg", height =10, width = 6)
