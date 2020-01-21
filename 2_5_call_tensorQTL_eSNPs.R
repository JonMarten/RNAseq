setwd("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/03_tensorqtl/results")

library(dplyr)
library(data.table)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

cis_nom <- fread("tensorqtl_cis_test.cis_qtl_pairs.22.csv", data.table = F)

cis_pheno <- fread("tensorqtl_cis_MAF0.005_cisnominal.tensorQTL.cis_nominal_manualPython.csv", data.table = F)

eGenes <- cis_pheno %>%
  mutate(p.bh = p.adjust(pval_perm, method = "BH")) 

ggplot(eGenes, aes(x = -log10(qval), y = -log10(p.bh))) +
  geom_point()