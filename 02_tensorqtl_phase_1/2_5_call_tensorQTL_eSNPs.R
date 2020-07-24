library(dplyr)
library(data.table)
library(ggplot2)
library(cowplot)
library(readr)
theme_set(theme_cowplot())

setwd("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/03_tensorqtl/results")

#cis_nom  <- list.files(pattern = "tensorqtl_cis_MAF0.005_cisnominal.cis_qtl_pairs.*.csv", full.names = TRUE) %>% 
#  lapply(read_csv) %>% 
#  bind_rows 
#
#fwrite(cis_nom, file = "tensorqtl_cis_MAF0.005_cisnominal.cis_qtl_pairs_merged.csv")
cis_nom <- fread("tensorqtl_cis_MAF0.005_cisnominal.cis_qtl_pairs_merged.csv", data.table = F)

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
eGenes <- eGenes %>%
  mutate(significant = ifelse(p.bh < 0.05, 1, 0))

m <- left_join(eGenes, eGenes.BH, by = "phenotype_id", suffix = c(".perm",".bonf"))

get_bh_threshold <- function(pvals, alpha, mtests = length(pvals)){
  m <- length(pvals)
  pvals <- sort(pvals)
  prejected <- which(pvals <= (1:m)/mtests*alpha)
  ifelse(length(prejected) == 0, 0, pvals[prejected[which.max(prejected)]])
}

bonfBH_threshold <- get_bh_threshold(eGenes.BH$p_global_BH, 0.05)

cis_nom2 <- cis_nom2 %>%
  mutate(eSNP_sig = ifelse(p_local_bonf < bonfBH_threshold, 1, 0))

cis_nom2 <- left_join(cis_nom2, select(cis_pheno, phenotype_id, pval_nominal_threshold))
cis_nom2 <- cis_nom2 %>%
  mutate(eSNP_sig_permthresh = ifelse(pval_nominal < pval_nominal_threshold, 1, 0)) 


table("bonf-BH" = cis_nom2$eSNP_sig, "tensorQTL" = cis_nom2$eSNP_sig_permthresh)

ggplot(eGenes, aes(x = -log10(qval), y = -log10(p.bh))) +
  geom_point()