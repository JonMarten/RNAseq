# Get Affy IDs for RNA seq individuals
setwd("U:/projects/RNAseq/")
library(data.table)
library(dplyr)
dat <- fread("RNA_sample_selection/INTERVAL_omics_table_15AUG2019.csv", data.table = F)
out <- dat %>% 
  filter(!is.na(picklistRNA)) %>%
  select(affymetrix_ID)
out$IID <- out$affymetrix_ID
names(out)[1] <- "FID"
write.table(out, row.names = F, col.names = F, quote = F, file = "genetic_PCs/rna_seq_5k_affy_ids.txt")

# Get IDs for only phase 1 RNA seq samples
rna_ids <- fread("U:/Projects/RNAseq/globus/RNA_seq_ids.csv", data.table = F)
out2k <- dat[which(dat$RNA_ID %in% rna_ids$RNA_id),]
out2k$IID <- out2k$affymetrix_ID
out2k$FID <- out2k$affymetrix_ID
out2k <- out2k %>%
  select(IID, FID)
write.table(out2k, row.names = F, col.names = F, quote = F, file = "genetic_PCs/rna_seq_2.7k_affy_ids.txt")
