# Get RNA seq identifiers from gene count data for all 8 batches and merge into a single file with ID and batch
setwd("/home/jm2294/rds/rds-jmmh2-pre_qc_data/interval/rna_seq/raw_data/globus")
library(data.table)
library(dplyr)

batch1 <- fread("batch1/matrices/counts.FC.unstranded.5281-tic109.txt", data.table = F)
b1_names <- data.frame("RNA_id" = names(batch1), "batch" = 1)
batch2 <- fread("batch2/results-study5591-tic109b/combined/study5591-tic109b-star-genecounts.txt", data.table = F)
b2_names <- data.frame("RNA_id" = names(batch2), "batch" = 2)
batch3 <- fread("batch3/results-study5591-tic109d/combined/study5591-tic109d-star-genecounts.txt", data.table = F)
b3_names <- data.frame("RNA_id" = names(batch3), "batch" = 3)
batch4 <- fread("batch4/results-study5591-tic276/combined/study5591-tic276-star-genecounts.txt", data.table = F)
b4_names <- data.frame("RNA_id" = names(batch4), "batch" = 4)

batch5 <- fread("batch5/results-study5591-tic328/combined/study5591-tic328-star-genecounts.txt", data.table = F)
b5_names <- data.frame("RNA_id" = names(batch5), "batch" = 5)

batch6 <- fread("batch6/results-study5591-tic297/combined/study5591-tic297-star-genecounts.txt", data.table = F)
b6_names <- data.frame("RNA_id" = names(batch6), "batch" = 6)
batch7 <- fread("batch7/results-study5591-tic329/combined/study5591-tic329-star-genecounts.txt", data.table = F)
b7_names <- data.frame("RNA_id" = names(batch7), "batch" = 7)
batch8 <- fread("batch8/results-study5591-tic364/combined/study5591-tic364-star-genecounts.txt", data.table = F)
b8_names <- data.frame("RNA_id" = names(batch8), "batch" = 8)

names <- rbind(b1_names, b2_names, b3_names, b4_names, b5_names, b6_names, b7_names, b8_names)
names <- names[which(grepl("INT", names$RNA_id)),] # remove rows for columns other than RNA ids
write.csv(names, "RNA_seq_ids.csv", quote = F, row.names = F)