library(dplyr)
library(data.table)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
setwd("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/03_tensorqtl/results")
a <- fread("tensorqtl_cis_MAF0.005_cisnominal.cis_qtl_pairs.22.csv", data.table = F)

limix22 <- fread("../../01_cis_eqtl_mapping/results/cis_eqtls_18373genes_age_sex_rin_batch_PC10_PEER20_5GenesPerChunk/processed/cis_eqtls_18373genes_age_sex_rin_batch_PC10_PEER20_eSNPs_eQTLgenReplication.txt", data.table = F)

a <- a %>% 
  mutate(matchID = paste0(phenotype_id, ":", variant_id))
plink_map <- fread("/rds/user/jm2294/rds-jmmh2-projects/interval_rna_seq/analysis/03_tensorqtl/genotypes/INTERVAL_b38_autosomes_RNAseqPhase1_biallelic_all_MAF0.005.bim" , data.table = F)
names(plink_map) <- c("chr", "variant_id", "morg", "pos", "A1","A2")
a2 <- left_join(a, plink_map)

m <- inner_join(limix22, a2, by = "matchID")

m <- m %>% 
  mutate(tqFlip = ifelse(A1 == assessed_allele, 0, ifelse(A2 == assessed_allele, 1, NA))) %>%
  mutate(slopeFlip = ifelse(tqFlip == 0, slope, -slope))


eql_lim <- ggplot(m, aes(y =  Zscore.eqtlgen.flipped, x = zscore)) +
  geom_point(pch = 20, cex = 0.8, alpha = 0.5) +
  labs(x = "Limix Zscore", y = "eQTLgenQTL Zscore")

tens_lim <- ggplot(m, aes(y =  slopeFlip/slope_se, x = zscore)) +
  geom_point(pch = 20, cex = 0.8, alpha = 0.5) +
  labs(x = "Limix Zscore", y = "TensorQTL Zscore")

tens_eql <- ggplot(m, aes(y =  slopeFlip/slope_se, x = Zscore.eqtlgen.flipped)) +
  geom_point(pch = 20, cex = 0.8, alpha = 0.5) +
  labs(x = "eQTLgen Zscore", y = "TensorQTL Zscore")

plot_grid(eql_lim, tens_lim, tens_eql)


eQTLgen <- fread("zcat /rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/eQTLgen/cis-eQTLs_full_20180905.txt.gz", data.table = F)
limix <- fread("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/results/cis_eqtls_18373genes_age_sex_rin_batch_PC10_PEER20_5GenesPerChunk/processed/output_results_merged.txt", data.table = F)
eQTLgen2 <- eQTLgen %>% 
  filter(SNPChr == 22) %>%
  mutate(matchID = paste0(Gene, ":", SNP))
m3 <- left_join(a2, eQTLgen2, by = "matchID")

limix2 <- limix %>% 
  filter(feature_chromosome == 22) %>%
  mutate(matchID = paste0(feature_id, ":", rsid))

m4 <- left_join(m3, limix2, by = "matchID")
m5 <- m4 %>% 
  filter(Pvalue < 5e-8) %>%
  filter(!is.na(SNPChr) | !is.na(snp_chromosome)) %>%
  mutate(flipE = ifelse(A1 == AssessedAllele, 0, 1),
         flipL = ifelse(A1 == assessed_allele, 0, 1)) %>%
  mutate(limixBetaFlip = ifelse(flipL == 1, -beta, beta),
         eQTLgenZFlip = ifelse(flipE == 1, -Zscore, Zscore)) %>%
  mutate(tensorZ = slope / slope_se,
         limixZ = limixBetaFlip/beta_se,
         eQTLgenZ = eQTLgenZFlip) %>%
  mutate(tensorZ = (tensorZ - mean(tensorZ, na.rm=T))/sd(tensorZ, na.rm = T),
         limixZ = (limixZ - mean(limixZ, na.rm=T))/sd(limixZ, na.rm = T),
         eQTLgenZ = (eQTLgenZ - mean(eQTLgenZ, na.rm=T))/sd(eQTLgenZ, na.rm = T)) %>%
  mutate(limix_limit = ifelse(abs(limixZ) > 3.439, 1, 0)) %>%
  mutate(limix_limit = as.factor(limix_limit))

LT <- ggplot(filter(m5, limix_limit == 1), aes(y =  tensorZ, x = limixZ, colour = limix_limit)) +
  geom_point(pch = 20, cex = 0.8, alpha = 0.5) +
  labs(x = "Limix Zscore", y = "TensorQTL Zscore")

LE <- ggplot(filter(m5, limix_limit == 1), aes(y =  eQTLgenZ, x = limixZ, colour = limix_limit)) +
  geom_point(pch = 20, cex = 0.8, alpha = 0.5) +
  labs(x = "Limix Zscore", y = "eQTLgen Zscore")

ET <- ggplot(filter(m5, limix_limit == 1), aes(y =  tensorZ, x = eQTLgenZ, colour = limix_limit)) +
  geom_point(pch = 20, cex = 0.8, alpha = 0.5) +
  labs(x = "eQTLgen Zscore", y = "TensorQTL Zscore")

plot_grid(LE, LT, ET)


m5 <- m5 %>% 
  mutate(log10pdiff = log10(pval_nominal) - log10(p_value))

ggplot(filter(m5, abs(log10pdiff) > 15), aes(x = -log10(pval_nominal), y = -log10(p_value), color = log10pdiff)) + 
  labs(x = "TensorQTL P", y = "Limix P") + 
  geom_point(pch = 20, cex = 0.8, alpha = 0.5)

tP <- ggplot(m5, aes(x = pos, y = -log10(pval_nominal))) +
  geom_point()

eP <- ggplot(m5, aes(x = pos, y = -log10(Pvalue))) +
  geom_point()

lP <- ggplot(m5, aes(x = pos, y = -log10(p_value))) +
  geom_point()

plot_grid(ncol = 1, tP, eP, lP)