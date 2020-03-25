# 

setwd("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/covid19")
library(data.table)
library(dplyr)

phe <- fread(data.table = F, "UNfilteredSamplesGenes_TMMNormalised_FPKM_Counts.csv")

cov19 <- fread("/rds/user/jm2294/rds-jmmh2-projects/interval_rna_seq/scripts/RNAseq/covid-19/covid_genes_b37.csv", data.table = F)

phecov19 <- phe %>% filter(feature_id %in% cov19$ensembl_id)

fullphe <- fread("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/globus_phase2/results/combined/interval_basic-star-fc-genecounts.txt", data.table = F)

phecov19 <- fullphe %>% filter(ENSEMBL_ID %in% cov19$ensembl_id)



covgenes1 <- phe %>% filter(feature_id %in% cov19$ensembl_id)
idmap <- fread("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/phenotype/sample_mapping_file_gt_to_phe_phase1.txt", data.table = F)

rnaids <- names(covgenes1)[-c(1:6)]
gtids <- idmap$genotype_individual_id[match(rnaids, idmap$phenotype_individual_id)]
names(covgenes1)[-c(1:6)] <- gtids

fwrite(covgenes1, sep = ",", file = "Covid-19_Genes_TMMNormalised_FPKM_Counts.csv")