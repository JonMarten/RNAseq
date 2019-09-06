library(dplyr)
library(data.table)
setwd("U:/Projects/RNAseq/covariates")

covs <- fread(data.table = F, "INTERVAL_RNA_batch1-8_covariates_release_2019_08_15.csv")

covout <- covs %>%
  select(sample_id, Agilent_RINe, age_RNA, sequencingBatch, sex) %>%
  mutate(Agilent_RINe = as.numeric(Agilent_RINe)) %>%
  filter(!is.na(Agilent_RINe))

fwrite(covout, file = "INTERVAL_RNAseq_phase1_covariates.txt", sep = "\t")

