setwd("/home/jm2294/projects/RNAseq/globus/")
library(data.table)
library(dplyr)

batch1 <- fread("batch1/matrices/counts.FC.unstranded.5281-tic109.txt", data.table=F)
b1_names <- data.frame("RNA_id" = names(batch1), "batch" = 1)
batch2 <- fread("batch2/results-study5591-tic109b/combined/study5591-tic109b-star-genecounts.txt")
b2_names <- data.frame("RNA_id" = names(batch2), "batch" = 2)
batch3 <- fread("batch3/results-study5591-tic109d/combined/study5591-tic109d-star-genecounts.txt", data.table=F)
b3_names <- data.frame("RNA_id" = names(batch3), "batch" = 3)
batch4 <- fread("batch4/results-study5591-tic276/combined/study5591-tic276-star-genecounts.txt")
b4_names <- data.frame("RNA_id" = names(batch4), "batch" = 4)

names <- rbind(b1_names, b2_names, b3_names, b4_names)
names <- names[which(grepl("INT", names$RNA_id)),] # remove rows for columns other than RNA ids
write.csv(names, "RNA_seq_ids.csv", quote = F, row.names = F)