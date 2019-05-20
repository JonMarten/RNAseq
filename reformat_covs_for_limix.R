# Reformat covariates for LIMIX pipeline
library(dplyr)
library(data.table)
setwd("U:/Projects/RNAseq")

dat <- fread("covariates/sarah_data_2/INTERVALdata_13MAY2019.csv", data.table=F)
datp3 <- fread("covariates/sarah_data_2/INTERVALdata_P3_13MAY2019.csv", data.table=F)
rna_id_mapper <- fread("rna_id_mapper.csv", data.table=F)

dat <- full_join(dat, datp3, by="identifier")
dat <- left_join(rna_id_mapper, dat)
dat <- dat %>%
  select(sample_id = RNA_id, sexPulse, agePulse)
fwrite(sep = "\t", dat, file = "test_run/INTERVAL_RNA_batch1_2_covariates_sex_age.txt")
