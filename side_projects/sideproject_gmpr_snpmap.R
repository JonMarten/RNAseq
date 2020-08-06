# Make SNP map file for GMPR side project
library(data.table)
library(dplyr)
setwd("/rds/user/jm2294/rds-jmmh2-projects/interval_rna_seq/analysis/04_phase2_full_analysis/side_projects/gmpr")

snplist <- fread("GMPR_variants_for_lookup_20200717_LS_sorted.txt", data.table = F)

map <- fread("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/04_phase2_full_analysis/genotypes/INTERVAL_RNAseq_Phase1-2_imputed_b38_biallelic_chr6_nofilters.bim", data.table = F)
names(map) <- c("chr","rsID","morg","pos_b38","effect_allele","other_allele")

map2 <- map %>%
  filter(rsID %in% snplist$rsID)

fwrite(map2, file = "GMPR_snp_map.csv")