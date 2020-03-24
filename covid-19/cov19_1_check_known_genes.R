# 

setwd("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/phenotype")
library(data.table)
library(dplyr)

phe <- fread(data.table = F, "gene_expr_filtered_TMM_rpkms_genesInvRankTransf_2747Samp_18373genes.csv")

cov19 <- fread("/rds/user/jm2294/rds-jmmh2-projects/interval_rna_seq/scripts/RNAseq/covid-19/covid_genes_b37.csv", data.table = F)

phecov19 <- phe %>% filter(gene_id %in% cov19$ensembl_id)

fullphe <- fread("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/globus_phase2/results/combined/interval_basic-star-fc-genecounts.txt", data.table = F)

phecov19 <- fullphe %>% filter(ENSEMBL_ID %in% cov19$ensembl_id)
