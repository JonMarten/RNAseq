setwd("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/03_tensorqtl/results")

library(dplyr)
library(data.table)
library(ggplot2)
library(cowplot)
library(readr)
theme_set(theme_cowplot())

#cis_nom <- fread("tensorqtl_cis_MAF0.005_cisnominal.cis_qtl_pairs.22.csv", data.table = F)

cis_nom  <- list.files(pattern = "tensorqtl_cis_MAF0.005_cisnominal.cis_qtl_pairs.*.csv", full.names = TRUE) %>% 
  lapply(read_csv) %>% 
  bind_rows 

fwrite(cis_nom, file = "tensorqtl_cis_MAF0.005_cisnominal.cis_qtl_pairs_merged.csv")

cis_pheno <- fread("tensorqtl_cis_MAF0.005_cisnominal.tensorQTL.cis_nominal_manualPython.csv", data.table = F)

eGenes <- cis_pheno %>%
  mutate(p.bh = p.adjust(pval_perm, method = "BH")) 

cis_nom2 <- cis_nom %>%
  group_by(phenotype_id) %>%
  mutate(p_local_bonf = pval_nominal * n()) %>%
  mutate(p_local_bonf = ifelse(p_local_bonf >= 1, 1, p_local_bonf)) %>%
  data.frame()

eGenes.BH <- cis_nom2 %>%
  group_by(phenotype_id) %>%
  summarise(nSNPs = n(),
            min_p_local_bonf = min(p_local_bonf, na.rm = T)) %>%
  data.frame()

eGenes.BH <- eGenes.BH %>%
  mutate(p_global_BH = p.adjust(min_p_local_bonf, method = "BH"),
         significant = ifelse(p_global_BH < 0.05, 1, 0))




ggplot(eGenes, aes(x = -log10(qval), y = -log10(p.bh))) +
  geom_point()