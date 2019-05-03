# Reformat covariates for LIMIX pipeline
library(dplyr)
library(data.table)
setwd("U:/Projects/RNAseq")

dat <- fread("covariates/INTERVALdata_02APR2019.csv", data.table=F)
rna_id_mapper <- fread("rna_id_mapper.csv", data.table=F)

dat <- left_join(rna_id_mapper, dat)
dat <- dat %>%
  select(sample_id = RNA_id, sexPulse, agePulse)
fwrite(sep = "\t", dat, file = "test_run/INTERVAL_RNA_batch1_2_covariates_sex_age.txt")
