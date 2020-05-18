# Get Affy IDs for RNA seq individuals
setwd("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/04_phase2_full_analysis/")
library(data.table)
library(dplyr)
dat <- fread("covariates/processed/INTERVAL_omics_table_14MAY2020.csv", data.table = F)
out <- dat %>% 
  filter(!is.na(RNA_ID)) %>%
  select(affymetrix_ID)
out$IID <- out$affymetrix_ID
names(out)[1] <- "FID"
write.table(out, row.names = F, col.names = F, quote = F, file = "genotypes/processing/rna_seq_phase1-2_affy_ids.txt")

