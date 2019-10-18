setwd("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/coviariates")
library(data.table)
library(dplyr)

cov <- fread("INTERVAL_RNA_batch1-8_covariates_release_2019_08_15.csv", data.table = F)
peer <- fread("PEER_factors_100Factors_plusCovariates_1000Iterations_AllBatches_TransformGenes_18373Genes.csv", data.table = F)

merge <- inner_join(cov, peer)