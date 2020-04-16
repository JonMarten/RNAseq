setwd("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/04_phase2_full_analysis/covariates")

library(data.table)
library(dplyr)

techCovsPhase1 <- fread("raw/INTERVAL_RNA_Covariates_09-08-2019.csv", data.table = F)
techCovsPhase2 <- fread("raw/INTERVAL_RNA_Covariates_Ph2_15-04-2020.csv", data.table = F)

techCovsPhase2 <- techCovsPhase2 %>% 
  mutate_at(c("Sample_ID", "Lane", "Freezer_Shelf", "Agilent_RINe"), as.character)

techCovsAll <- bind_rows(techCovsPhase1, techCovsPhase2)

fwrite(techCovsAll, "processed/INTERVAL_RNA_technical_covariates_batch1-12_20200416.csv")
