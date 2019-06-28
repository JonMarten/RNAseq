# Get Affy IDs for RNA seq individuals
setwd("U:/projects/RNAseq/")
library(data.table)
library(dplyr)
dat <- fread("RNA_sample_selection/INTERVAL_omics_table_11JUN2019.csv", data.table=F)
out <- dat %>% 
  filter(!is.na(picklistRNA)) %>%
  select(affymetrix_ID)
out$IID <- out$affymetrix_ID
names(out)[1] <- "FID"
write.table(out, row.names = F, col.names = F, quote = F, file = "genetic_PCs/rna_seq_5k_affy_ids.txt")
