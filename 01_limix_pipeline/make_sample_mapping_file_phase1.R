library(data.table)
library(dplyr)
library(tidyr)
library(cowplot)

setwd("U:/Projects/RNAseq/analysis/00_testing/phenotype")

a <- fread("U:/Projects/RNAseq/RNA_sample_selection/INTERVAL_omics_table_15AUG2019.csv", data.table = F)

b <- a %>% 
  select(genotype_individual_id = affymetrix_ID, phenotype_individual_id = RNA_ID) %>%
  filter(!is.na(phenotype_individual_id))

write.table(b, sep = "\t", quote = F, col.names = T, row.names = F, file = "sample_mapping_file_gt_to_phe_phase1.txt")

# Write out sample inclusion file for QC tool to make small bgens
samples <- b$genotype_individual_id
write.table(samples, quote = F, col.names = F, row.names = F, file = "../GENETIC_DATA/b37_b38_liftover/rnaseq_affy_ids.txt")