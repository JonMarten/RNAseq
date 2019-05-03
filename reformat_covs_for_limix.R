# Reformat covariates for LIMIX pipeline
library(dplyr)
library(data.table)
setwd("U:/Projects/RNAseq")

dat <- fread("covariates/INTERVALdata_02APR2019.csv", data.table=F)
rna_id_mapper <- fread("rna_id_mapper.csv", data.table=F)

