library(data.table)
setwd("~/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/pio/projects/20230313_Dmel_pangenome/")

br <- fread("pantera_output7/blast_results", header = F)
colnames(br) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                 "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen")

br[,slen:=max(c(send,sstart)),by = sseqid]
br[,pl:=length/slen]
br[,score:= (length*pident/100)/qlen]
br_s <- br[, max(score),by= sseqid]


print("Number of Ref matches:")
print(nrow(br_s))

print("Number of matches > 99%")
print(sum(br_s$V1>=0.99))
print("Number of matches > 95%")
print(sum(br_s$V1>=0.95))
print("Number of matches > 80%")
print(sum(br_s$V1>=0.80))


### Now comparing blast to A1

bA1_p7 <- fread("pantera_output7/pantera7_blast_to_A1", header = F)
bA1_dmel <- fread("pantera_output7/Dmelref_blast_to_A1", header = F)
colnames(bA1_p7) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                  "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen")
colnames(bA1_dmel) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                  "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen")

bA1_p7[,pl:=length/qlen]
bA1_p7[,score:= (length*pident/100)/qlen]
bA1_p7_s <- bA1_p7[, max(score),by= qseqid]

print("Number of Ref matches:")
print(nrow(bA1_p7_s))

print("Number of matches > 99%")
print(sum(bA1_p7_s$V1>=0.99))
print("Number of matches > 95%")
print(sum(bA1_p7_s$V1>=0.95))
print("Number of matches > 80%")
print(sum(bA1_p7_s$V1>=0.80))

bA1_dmel[,pl:=length/qlen]
bA1_dmel[,score:= (length*pident/100)/qlen]
bA1_dmel_s <- bA1_dmel[, max(score),by= qseqid]

print("Number of Ref matches:")
print(nrow(bA1_dmel_s))

print("Number of matches > 99%")
print(sum(bA1_dmel_s$V1>=0.99))
print("Number of matches > 95%")
print(sum(bA1_dmel_s$V1>=0.95))
print("Number of matches > 80%")
print(sum(bA1_dmel_s$V1>=0.80))

### List of good Dmel not in pantera
missing <- bA1_dmel_s[V1>=0.98][!(bA1_dmel_s[V1>=0.98]$qseqid %in% br_s[V1>0.8]$sseqid)]
found <- bA1_dmel_s[V1>=0.98][(bA1_dmel_s[V1>=0.98]$qseqid %in% br_s[V1>0.8]$sseqid)]

# With A1A2
b12 <- fread("pantera_output_A1A2/blast_Dmel_A1A2", header = F)
colnames(b12) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                  "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen")

b12[,slen:=max(c(send,sstart)),by = sseqid]
b12[,score:= (length*pident/100)/qlen]
b12[,pl:=length/slen]
b12_s <- b12[, max(score),by= sseqid]


print("Number of Ref matches:")
print(nrow(b12_s))

print("Number of matches > 99%")
print(sum(b12_s$V1>=0.99))
print("Number of matches > 95%")
print(sum(b12_s$V1>=0.95))
print("Number of matches > 80%")
print(sum(b12_s$V1>=0.80))

# With A1Dsim
b1dsim <- fread("pantera_output_A1Dsim/blast_Dmel_A1Dsim", header = F)
colnames(b1dsim) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                   "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen")

b1dsim[,slen:=max(c(send,sstart)),by = sseqid]
b1dsim[,score:= (length*pident/100)/qlen]
b1dsim[,pl:=length/slen]
b1dsim_s <- b1dsim[, max(score),by= sseqid]


print("Number of Ref matches:")
print(nrow(b1dsim_s))

print("Number of matches > 99%")
print(sum(b1dsim_s$V1>=0.99))
print("Number of matches > 95%")
print(sum(b1dsim_s$V1>=0.95))
print("Number of matches > 80%")
print(sum(b1dsim_s$V1>=0.80))

# With RM2
bRM2 <- fread("RM_93922.ThuJun221105402023/blast_Dmel_RM2", header = F)
colnames(bRM2) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                      "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen")

bRM2[,slen:=max(c(send,sstart)),by = sseqid]
bRM2[,score:= (length*pident/100)/qlen]
bRM2[,pl:=length/slen]
bRM2_s <- bRM2[, max(score),by= sseqid]


print("Number of Ref matches:")
print(nrow(bRM2_s))

print("Number of matches > 99%")
print(sum(bRM2_s$V1>=0.99))
print("Number of matches > 95%")
print(sum(bRM2_s$V1>=0.95))
print("Number of matches > 80%")
print(sum(bRM2_s$V1>=0.80))

# With looking at all segs of run7
ball <- fread("pantera_output7/blast_Dmel_allsegs", header = F)
colnames(ball) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                    "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen")

ball[,slen:=max(c(send,sstart)),by = sseqid]
ball[,score:= (length*pident/100)/qlen]
ball[,pl:=length/slen]
ball_s <- ball[, max(score),by= sseqid]


print("Number of Ref matches:")
print(nrow(ball_s))

print("Number of matches > 99%")
print(sum(ball_s$V1>=0.99))
print("Number of matches > 95%")
print(sum(ball_s$V1>=0.95))
print("Number of matches > 80%")
print(sum(ball_s$V1>=0.80))

# With looking at all segs of Dsim
balldsim <- fread("pantera_output_A1Dsim/blast_Dmel_allsegsDsim", header = F)
colnames(balldsim) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                    "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen")

balldsim[,slen:=max(c(send,sstart)),by = sseqid]
balldsim[,score:= (length*pident/100)/qlen]
balldsim[,pl:=length/slen]
balldsim_s <- balldsim[, max(score),by= sseqid]


print("Number of Ref matches:")
print(nrow(balldsim_s))

print("Number of matches > 99%")
print(sum(balldsim_s$V1>=0.99))
print("Number of matches > 95%")
print(sum(balldsim_s$V1>=0.95))
print("Number of matches > 80%")
print(sum(balldsim_s$V1>=0.80))

# Comparing the Ref mappings to Pantera mappings
compe_dmel_pant_in_A1 <- merge(bA1_dmel_s, br_s, by.x="qseqid", by.y= "sseqid", all = T)
colnames(compe_dmel_pant_in_A1) <- c("qseqid","dmel_score", "pant_score")
A1dmel_counts <- bA1_dmel[score>0.80,.N,by=qseqid]
compe_dmel_pant_in_A1 <- merge(compe_dmel_pant_in_A1, A1dmel_counts, all = T)
compe_dmel_pant_in_A1[,d:=(dmel_score-pant_score)^2]

### Pantera in Dmel (that have at least 1 hit better than 0.8 in the A1)
table(!is.na(compe_dmel_pant_in_A1[dmel_score>0.8 & N>0]$pant_score))

### Pantera in Dmel (that have at least 2 hits better than 0.8 in the A1)
table(!is.na(compe_dmel_pant_in_A1[dmel_score>0.8 & N>1]$pant_score))

### Pantera in Dmel (that have at least 10 hits better than 0.8 in the A1)
table(!is.na(compe_dmel_pant_in_A1[dmel_score>0.8 & N>9]$pant_score))

a <- c(0:100)
b <- unlist(lapply(a,function(x){ table(!is.na(compe_dmel_pant_in_A1[dmel_score>0.8 & N>x]$pant_score))[2]/nrow(compe_dmel_pant_in_A1[dmel_score>0.8 & N>x])}))
plot(a,b)

### Now, do we find anything in A1 with pantera that is "real" and not in Dmel?
hits_not_inDmel <-bA1_p7[!(qseqid %in% br[pl>0.8]$qseqid)]



# Inverted Pantera mapping NOT in ref mappings, is there any?
br_s2 <- br[, .(max(pl),max(qlen),max(score)),by= qseqid]
compe_pant_dmel_in_A1 <- merge(bA1_p7_s, br_s2, by.x="qseqid", by.y= "qseqid", all = T)
colnames(compe_pant_dmel_in_A1) <- c("qseqid","p_score", "pant_score", "len", "score")
A1pant_counts <- bA1_p7[score>0.80,.N,by=qseqid]
compe_pant_dmel_in_A1 <- merge(compe_pant_dmel_in_A1, A1pant_counts, all = T)
compe_pant_dmel_in_A1[,d:=(p_score-pant_score)^2]

### Pantera in Dmel (that have at least 1 hit better than 0.8 in the A1)
table(!is.na(compe_pant_dmel_in_A1[p_score>0.8 & N>0]$pant_score))

### Pantera in Dmel (that have at least 2 hits better than 0.8 in the A1)
table(!is.na(compe_pant_dmel_in_A1[p_score>0.8 & N>1]$pant_score))

### Pantera in Dmel (that have at least 10 hits better than 0.8 in the A1)
table(!is.na(compe_pant_dmel_in_A1[p_score>0.8 & N>9]$pant_score))

a <- c(0:100)
b <- unlist(lapply(a,function(x){ table(!is.na(compe_pant_dmel_in_A1[p_score>0.8 & N>x]$pant_score))[2]/nrow(compe_pant_dmel_in_A1[p_score>0.8 & N>x])}))
plot(a,b)

### Now, do we find anything in A1 with pantera that is "real" and not in Dmel?
hits_not_inDmel <-bA1_p7[!(qseqid %in% br[pl>0.8]$qseqid)]


# From start with pansoft
br <- fread("pantera_output7/blast_Dmel_pan7soft", header = F)
colnames(br) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                  "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen")

br[,slen:=max(c(send,sstart)),by = sseqid]
br[,pl:=length/slen]
br[,score:= (length*pident/100)/qlen]
br_s <- br[, max(score),by= sseqid]


print("Number of Ref matches:")
print(nrow(br_s))

print("Number of matches > 99%")
print(sum(br_s$V1>=0.99))
print("Number of matches > 95%")
print(sum(br_s$V1>=0.95))
print("Number of matches > 80%")
print(sum(br_s$V1>=0.80))

### How many sequences in pantlib have a good match in the Ref:
br_pant <- br[, max(score),by= qseqid]
table(br_pant$V1>0.8)

### Now comparing blast to A1

bA1_p7 <- fread("pantera_output7/blast_pan7soft_A1", header = F)
bA1_dmel <- fread("pantera_output7/Dmelref_blast_to_A1", header = F)
colnames(bA1_p7) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                      "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen")
colnames(bA1_dmel) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                        "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen")

bA1_p7[,pl:=length/qlen]
bA1_p7[,score:= (length*pident/100)/qlen]
bA1_p7_s <- bA1_p7[, max(score),by= qseqid]

print("Number of Ref matches:")
print(nrow(bA1_p7_s))

print("Number of matches > 99%")
print(sum(bA1_p7_s$V1>=0.99))
print("Number of matches > 95%")
print(sum(bA1_p7_s$V1>=0.95))
print("Number of matches > 80%")
print(sum(bA1_p7_s$V1>=0.80))

bA1_dmel[,pl:=length/qlen]
bA1_dmel[,score:= (length*pident/100)/qlen]
bA1_dmel_s <- bA1_dmel[, max(score),by= qseqid]

print("Number of Ref matches:")
print(nrow(bA1_dmel_s))

print("Number of matches > 99%")
print(sum(bA1_dmel_s$V1>=0.99))
print("Number of matches > 95%")
print(sum(bA1_dmel_s$V1>=0.95))
print("Number of matches > 80%")
print(sum(bA1_dmel_s$V1>=0.80))

### List of good Dmel not in pantera
missing <- bA1_dmel_s[V1>=0.98][!(bA1_dmel_s[V1>=0.98]$qseqid %in% br_s[V1>0.8]$sseqid)]
found <- bA1_dmel_s[V1>=0.98][(bA1_dmel_s[V1>=0.98]$qseqid %in% br_s[V1>0.8]$sseqid)]

# Comparing the Ref mappings to Pantera mappings
compe_dmel_pant_in_A1 <- merge(bA1_dmel_s, br_s, by.x="qseqid", by.y= "sseqid", all = T)
colnames(compe_dmel_pant_in_A1) <- c("qseqid","dmel_score", "pant_score")
A1dmel_counts <- bA1_dmel[score>0.80,.N,by=qseqid]
compe_dmel_pant_in_A1 <- merge(compe_dmel_pant_in_A1, A1dmel_counts, all = T)
compe_dmel_pant_in_A1[,d:=(dmel_score-pant_score)^2]

### Pantera in Dmel (that have at least 1 hit better than 0.8 in the A1)
table(!is.na(compe_dmel_pant_in_A1[dmel_score>0.8 & N>0]$pant_score))

### Pantera in Dmel (that have at least 2 hits better than 0.8 in the A1)
table(!is.na(compe_dmel_pant_in_A1[dmel_score>0.8 & N>1]$pant_score))

### Pantera in Dmel (that have at least 10 hits better than 0.8 in the A1)
table(!is.na(compe_dmel_pant_in_A1[dmel_score>0.8 & N>9]$pant_score))

### Same, requering that the match of pantera to dmel be good.
# table(compe_dmel_pant_in_A1[!is.na(compe_dmel_pant_in_A1[dmel_score>0.8 & N>0]$pant_score)][pant_score>0.8])



a <- c(0:100)
b <- unlist(lapply(a,function(x){ table(!is.na(compe_dmel_pant_in_A1[dmel_score>0.8 & N>x]$pant_score))[2]/nrow(compe_dmel_pant_in_A1[dmel_score>0.8 & N>x])}))
plot(a,b)

### Now, do we find anything in A1 with pantera that is "real" and not in Dmel?
hits_not_inDmel <-bA1_p7[!(qseqid %in% br[pl>0.8]$qseqid)]

# Inverted Pantera mapping NOT in ref mappings, is there any?
br_s2 <- br[, .(max(pl),max(qlen),max(score)),by= qseqid]
bA1_p7_s2 <- bA1_p7[, max(score),by= qseqid]
compe_pant_dmel_in_A1 <- merge(bA1_p7_s2, br_s2, by.x="qseqid", by.y= "qseqid", all = T)
colnames(compe_pant_dmel_in_A1) <- c("qseqid","p_score", "pant_score", "len", "score")
A1pant_counts <- bA1_p7[score>0.8,.N,by=qseqid]
compe_pant_dmel_in_A1 <- merge(compe_pant_dmel_in_A1, A1pant_counts, all = T)
compe_pant_dmel_in_A1[,d:=(p_score-pant_score)^2]

### Pantera in Dmel (that have at least 1 hit better than 0.8 in the A1)
table(!is.na(compe_pant_dmel_in_A1[p_score>0.8 & N>0]$pant_score))

### Pantera in Dmel (that have at least 2 hits better than 0.8 in the A1)
table(!is.na(compe_pant_dmel_in_A1[p_score>0.8 & N>1]$pant_score))

### Pantera in Dmel (that have at least 10 hits better than 0.8 in the A1)
table(!is.na(compe_pant_dmel_in_A1[p_score>0.8 & N>9]$pant_score))

a <- c(0:100)
b <- unlist(lapply(a,function(x){ table(!is.na(compe_pant_dmel_in_A1[p_score>0.8 & N>x]$pant_score))[2]/nrow(compe_pant_dmel_in_A1[p_score>0.8 & N>x])}))
plot(a,b)

### Now, do we find anything in A1 with pantera that is "real" and not in Dmel?
hits_not_inDmel <-bA1_p7[!(qseqid %in% br[pl>0.8]$qseqid)]


### With pantera11

br <- fread("pant112_all/blast_Dmel_pan11", header = F)
colnames(br) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                  "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen")

br[,slen:=max(c(send,sstart)),by = sseqid]
br[,pl:=length/slen]
br[,score:= (length*pident/100)/qlen]
br_s <- br[, max(score),by= sseqid]


print("Number of Ref matches:")
print(nrow(br_s))

print("Number of matches > 99%")
print(sum(br_s$V1>=0.99))
print("Number of matches > 95%")
print(sum(br_s$V1>=0.95))
print("Number of matches > 80%")
print(sum(br_s$V1>=0.80))

# Now comparing matches to A1
bA1_p7 <- fread("pant112_all/blast_pan11_A1", header = F)
bA1_dmel <- fread("pantera_output7/Dmelref_blast_to_A1", header = F)
colnames(bA1_p7) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                      "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen")
colnames(bA1_dmel) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                        "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen")

bA1_p7[,pl:=length/qlen]
bA1_p7[,score:= (length*pident/100)/qlen]
bA1_p7_s <- bA1_p7[, max(score),by= qseqid]

print("Number of Ref matches:")
print(nrow(bA1_p7_s))

print("Number of matches > 99%")
print(sum(bA1_p7_s$V1>=0.99))
print("Number of matches > 95%")
print(sum(bA1_p7_s$V1>=0.95))
print("Number of matches > 80%")
print(sum(bA1_p7_s$V1>=0.80))

bA1_dmel[,pl:=length/qlen]
bA1_dmel[,score:= (length*pident/100)/qlen]
bA1_dmel_s <- bA1_dmel[, max(score),by= qseqid]

print("Number of Ref matches:")
print(nrow(bA1_dmel_s))

print("Number of matches > 99%")
print(sum(bA1_dmel_s$V1>=0.99))
print("Number of matches > 95%")
print(sum(bA1_dmel_s$V1>=0.95))
print("Number of matches > 80%")
print(sum(bA1_dmel_s$V1>=0.80))

### List of good Dmel not in pantera
missing <- bA1_dmel_s[V1>=0.98][!(bA1_dmel_s[V1>=0.98]$qseqid %in% br_s[V1>0.8]$sseqid)]
found <- bA1_dmel_s[V1>=0.98][(bA1_dmel_s[V1>=0.98]$qseqid %in% br_s[V1>0.8]$sseqid)]


# Comparing the Ref mappings to Pantera mappings
compe_dmel_pant_in_A1 <- merge(bA1_dmel_s, br_s, by.x="qseqid", by.y= "sseqid", all = T)
colnames(compe_dmel_pant_in_A1) <- c("qseqid","dmel_score", "pant_score")
A1dmel_counts <- bA1_dmel[score>0.80,.N,by=qseqid]
compe_dmel_pant_in_A1 <- merge(compe_dmel_pant_in_A1, A1dmel_counts, all = T)
compe_dmel_pant_in_A1[,d:=(dmel_score-pant_score)^2]

### Pantera in Dmel (that have at least 1 hit better than 0.8 in the A1)
table(!is.na(compe_dmel_pant_in_A1[dmel_score>0.8 & N>0]$pant_score))

### Pantera in Dmel (that have at least 2 hits better than 0.8 in the A1)
table(!is.na(compe_dmel_pant_in_A1[dmel_score>0.8 & N>1]$pant_score))

### Pantera in Dmel (that have at least 10 hits better than 0.8 in the A1)
table(!is.na(compe_dmel_pant_in_A1[dmel_score>0.8 & N>9]$pant_score))

### New binned version
### Pantera in Dmel (that have 1 hit better than 0.8 in the A1)
table(!is.na(compe_dmel_pant_in_A1[dmel_score>0.8 & N==1]$pant_score))

### Pantera in Dmel (that have between 2 and 9 hits better than 0.8 in the A1)
table(!is.na(compe_dmel_pant_in_A1[dmel_score>0.8 & N>1 & N<10]$pant_score))

### Pantera in Dmel (that have at least 10 hits better than 0.8 in the A1)
table(!is.na(compe_dmel_pant_in_A1[dmel_score>0.8 & N>9]$pant_score))


a <- c(0:100)
b <- unlist(lapply(a,function(x){ table(!is.na(compe_dmel_pant_in_A1[dmel_score>0.8 & N>x]$pant_score))[2]/nrow(compe_dmel_pant_in_A1[dmel_score>0.8 & N>x])}))
plot(a,b)

### Now, do we find anything in A1 with pantera that is "real" and not in Dmel?
hits_not_inDmel <-bA1_p7[!(qseqid %in% br[pl>0.8]$qseqid)]

# With looking at all segs of p17cA1A2
ball <- fread("pantera17c2_outputA1A2/blast_DmelTEs_p17A1A2final", header = F)
colnames(ball) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                    "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen")

ball[,slen:=max(c(send,sstart)),by = sseqid]
ball[,score:= (length*pident/100)/qlen]
ball[,pl:=length/slen]
ball_s <- ball[, max(score),by= sseqid]


print("Number of Ref matches:")
print(nrow(ball_s))

print("Number of matches > 99%")
print(sum(ball_s$V1>=0.99))
print("Number of matches > 95%")
print(sum(ball_s$V1>=0.95))
print("Number of matches > 80%")
print(sum(ball_s$V1>=0.80))

### With alternate score
ball[,score2:= (length*pident/100)/slen]
ball_s <- ball[, max(score2),by= sseqid]


print("Number of Ref matches:")
print(nrow(ball_s))

print("Number of matches > 99%")
print(sum(ball_s$V1>=0.99))
print("Number of matches > 95%")
print(sum(ball_s$V1>=0.95))
print("Number of matches > 80%")
print(sum(ball_s$V1>=0.80))



### With pantera20

br <- fread("pantera20_output/blast_Dmel_p20", header = F)
colnames(br) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                  "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen")

br[,slen:=max(c(send,sstart)),by = sseqid]
br[,pl:=length/slen]
br[,score:= (length*pident/100)/qlen]
br_s <- br[, max(score),by= sseqid]


print("Number of Ref matches:")
print(nrow(br_s))

print("Number of matches > 99%")
print(sum(br_s$V1>=0.99))
print("Number of matches > 95%")
print(sum(br_s$V1>=0.95))
print("Number of matches > 80%")
print(sum(br_s$V1>=0.80))

# Now comparing matches to A1
bA1_dmel <- fread("pantera_output7/Dmelref_blast_to_A1", header = F)
colnames(bA1_dmel) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                        "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen")


bA1_dmel[,pl:=length/qlen]
bA1_dmel[,score:= (length*pident/100)/qlen]
bA1_dmel_s <- bA1_dmel[, max(score),by= qseqid]

print("Number of Ref matches:")
print(nrow(bA1_dmel_s))

print("Number of matches > 99%")
print(sum(bA1_dmel_s$V1>=0.99))
print("Number of matches > 95%")
print(sum(bA1_dmel_s$V1>=0.95))
print("Number of matches > 80%")
print(sum(bA1_dmel_s$V1>=0.80))

### List of good Dmel not in pantera
missing <- bA1_dmel_s[V1>=0.98][!(bA1_dmel_s[V1>=0.98]$qseqid %in% br_s[V1>0.8]$sseqid)]
found <- bA1_dmel_s[V1>=0.98][(bA1_dmel_s[V1>=0.98]$qseqid %in% br_s[V1>0.8]$sseqid)]


# Comparing the Ref mappings to Pantera mappings
compe_dmel_pant_in_A1 <- merge(bA1_dmel_s, br_s, by.x="qseqid", by.y= "sseqid", all = T)
colnames(compe_dmel_pant_in_A1) <- c("qseqid","dmel_score", "pant_score")
A1dmel_counts <- bA1_dmel[score>0.80,.N,by=qseqid]
compe_dmel_pant_in_A1 <- merge(compe_dmel_pant_in_A1, A1dmel_counts, all = T)
compe_dmel_pant_in_A1[,d:=(dmel_score-pant_score)^2]

### Pantera in Dmel (that have at least 1 hit better than 0.8 in the A1)
table(!is.na(compe_dmel_pant_in_A1[dmel_score>0.8 & N>0]$pant_score))

### Pantera in Dmel (that have at least 2 hits better than 0.8 in the A1)
table(!is.na(compe_dmel_pant_in_A1[dmel_score>0.8 & N>1]$pant_score))

### Pantera in Dmel (that have at least 10 hits better than 0.8 in the A1)
table(!is.na(compe_dmel_pant_in_A1[dmel_score>0.8 & N>9]$pant_score))

### New binned version
### Pantera in Dmel (that have 1 hit better than 0.8 in the A1)
table(!is.na(compe_dmel_pant_in_A1[dmel_score>0.8 & N==1]$pant_score))

### Pantera in Dmel (that have between 2 and 9 hits better than 0.8 in the A1)
table(!is.na(compe_dmel_pant_in_A1[dmel_score>0.8 & N>1 & N<10]$pant_score))

### Pantera in Dmel (that have at least 10 hits better than 0.8 in the A1)
table(!is.na(compe_dmel_pant_in_A1[dmel_score>0.8 & N>9]$pant_score))

### Official lib anno check
dmel_of <- fread(cmd= "sed 's/\\*/ /g' official_lib_anno/A1.sl.fa.out", sep = " ", skip = 3)
dmel_of <- dmel_of[V11!="Simple_repeat"][V11!="Low_complexity"]
unique(dmel_of$V10)
