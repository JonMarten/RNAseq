# Filter normalised gene counts to make input file for limix parameter test

setwd("/home/jm2294/rds/rds-jmmh2-projects/interval_rna_seq/analysis/00_testing/phenotype")
library(data.table)
library(dplyr)

phe <- fread("filteredSamplesGenes_TMMNormalised_FPKM_Counts.csv", data.table = F)

pheout <- phe %>%
  filter(class == "protein_coding") %>%
  select(-(gene_symbol:gene_length), -INT_RNA7427299, -INT_RNA7711096)

fwrite(pheout, sep = "\t", file = "INTERVAL_RNAseq_phase1_filteredSamplesGenes_TMMNormalised_FPKM_Counts_foranalysis.txt")