# Filter normalised gene counts to make input file for limix parameter test

setwd("/rds/user/jm2294/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/phenotype")

library(data.table)
library(dplyr)

phe <- fread("gene_expr_filtered_TMM_rpkms_genesInvRankTransf_2747Samp_18373genes.csv", data.table = F)

pheout <- phe %>%
  filter(chr %in% 1:22) %>%
  select(-(gene_symbol:gene_length), -INT_RNA7427299, -INT_RNA7711096)

fwrite(pheout, sep = "\t", file = "INTERVAL_RNAseq_phase1_filteredSamplesGenes_TMMNormalised_FPKM_Counts_foranalysis.txt")