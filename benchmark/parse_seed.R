library("data.table")
library("ggplot2")
setwd("~/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/pio/projects/20230313_Dmel_pangenome/")
dist <- fread(cmd = ("grep 'GF CC' pant20_anno/seeds.stk | awk '{print $19}'"), header = F)
type <- fread(cmd = ("grep 'GF TP' pant20_anno/seeds.stk | awk '{print $3}'"), sep = "", header = F)
data_pant <- data.table( type = type[, gsub(".*;","",gsub(";[^;]*$","",V1))], dist = dist$V1, lib = "pant")

dist <- fread(cmd = ("grep 'GF CC' official_lib_anno/seeds.stk | awk '{print $19}'"), header = F)
type <- fread(cmd = ("grep 'GF TP' official_lib_anno/seeds.stk | awk '{print $3}'"), sep = "", header = F)
data_ref <- data.table( type = type[, gsub(".*;","",gsub(";[^;]*$","",V1))], dist = dist$V1, lib = "ref")

dist <- fread(cmd = ("grep 'GF CC' RM_lib_anno/seeds.stk | awk '{print $19}'"), header = F)
type <- fread(cmd = ("grep 'GF TP' RM_lib_anno/seeds.stk | awk '{print $3}'"), sep = "", header = F)
data_RM <- data.table( type = type[, gsub(".*;","",gsub(";[^;]*$","",V1))], dist = dist$V1, lib = "RM")



data_plot <- rbindlist(list(data_pant,data_ref, data_RM))

ggplot(data_plot, aes(y= dist, fill = lib)) +
  geom_boxplot()
