# Compare to matrix eqtl results
library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)
#library(hexbin)
theme_set(theme_cowplot())

setwd("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/results/")

eSNP_replication <- fread("cis_eqtls_18373genes_age_sex_rin_batch_PC10_PEER20/processed/cis_eqtls_18373genes_age_sex_rin_batch_PC10_PEER20_eSNPs_eqtlgenReplication.txt" ,data.table = F)

# From 05_eQTLgen_comparison.R --------------------------------------------
mismatchGenes <- eSNP_replication %>%
  filter(signMismatch == 1) %>%
  pull(limix_feature_id) %>%
  unique()
# Classify replication type
eSNP_replication <- eSNP_replication %>%
  mutate(Significant = ifelse(sig_bonf_sign == 1, "Bonf. & direction match", 
                              ifelse(sigrep_bonf == 1 & sig_bonf_sign == 0, "Bonf. & direction mismatch",
                                     "NS"))) %>%
  mutate(Significant = factor(Significant, levels = c("NS","Bonf. & direction mismatch", "Bonf. & direction match")))

# Classify SNPs as strand-ambiguous
eSNP_replication$limix_feature_length <- eSNP_replication$limix_feature_end - eSNP_replication$limix_feature_start
eSNP_replication <- eSNP_replication %>%
  mutate(ambig = ifelse(eQTLgen_AssessedAllele == "A" & eQTLgen_OtherAllele == "T" |
                          eQTLgen_AssessedAllele == "T" & eQTLgen_OtherAllele == "A" |
                          eQTLgen_AssessedAllele == "C" & eQTLgen_OtherAllele == "G" |
                          eQTLgen_AssessedAllele == "G" & eQTLgen_OtherAllele == "C", 1, 0)) %>%
  mutate(veryambig = ifelse(ambig == 1 & limix_maf > 0.45, 1, 0))

# get list of features with more than one significant SNPs mismatched  
flexmismatchGenes <- eSNP_replication %>%
  filter(Significant  == "Bonf. & direction mismatch") %>%
  group_by(limix_feature_id) %>%
  summarise(n_mismatched = n()) %>%
  data.frame %>%
  filter(n_mismatched > 1) %>%
  pull(limix_feature_id)


# New code ----------------------------------------------------------------
c22_mateqtl <- fread("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/matrix_eQTL/cis_eQTLs_Chr22.csv", data.table = F)

eSNP_replication <- eSNP_replication %>%
  mutate(matchid2 = paste0(limix_gene_name, ":", limix_rsid))

merge <- inner_join(eSNP_replication, c22_mateqtl) %>%
  mutate(se_matrixeqtl = beta / statistic)

g1 <- ggplot(merge, aes(x = -statistic, y = limix_zscore, colour = Significant)) +
  geom_point(alpha = 0.1) + 
  geom_abline(slope = 1, intercept = 0)+
  geom_vline(xintercept = 0) + 
  geom_hline(yintercept = 0) + 

g2 <- ggplot(merge, aes(x = -statistic, y = eQTLgen_Zscore_flipped, colour = Significant)) +
  geom_point(alpha = 0.1) + 
  geom_abline(slope = 1, intercept = 0)+
  geom_vline(xintercept = 0) + 
  geom_hline(yintercept = 0) + 
g3 <- ggplot(merge, aes(x = limix_zscore, y = eQTLgen_Zscore_flipped, colour = Significant)) +
  geom_point(alpha = 0.1) + 
  geom_abline(slope = 1, intercept = 0)+
  geom_vline(xintercept = 0) + 
  geom_hline(yintercept = 0) + 

plot_grid(g1,g2,g3, nrow = 1)

mm1 <- ggplot(filter(merge, limix_feature_id %in% flexmismatchGenes), 
       aes(x = -statistic, y = eQTLgen_Zscore_flipped, color = eQTLgen_NrCohorts)) + 
  geom_vline(xintercept = 0) + 
  geom_hline(yintercept = 0) + 
  geom_point(cex = 0.8, alpha = 1, pch = 20) +
  theme(aspect.ratio = 1, legend.position = "right") +
  labs(x = "- Matrix eQTL Z-score", y = "eQTLgen Z-score") +
  theme_minimal() +
  scale_colour_viridis_c() +
  facet_wrap(. ~ limix_gene_name)  
save_plot(file = "cis_eqtls_18373genes_age_sex_rin_batch_PC10_PEER20/processed/cis_eqtls_18373genes_age_sex_rin_batch_PC10_PEER20_eSNPs_matrixeQTL_mismatch.png", mm1, base_height = 7, base_width = 12)

merge <- merge %>%
  mutate(weird = ifelse(abs(limix_zscore) > 37.53, 1, 0) %>% as.factor,
         minP = ifelse(limix_p_value == min(merge$limix_p_value), 1, 0)  %>% as.factor)

m1 <- ggplot(merge, aes(x = limix_beta, y = -beta, colour = minP)) + 
  geom_point(alpha = ifelse(merge$minP == 1, 1, 0.1)) + 
  geom_abline(slope = 1, intercept = 0) +
  labs(x = "Limix Beta", y = "Matrix eQTL Beta")

m2 <- ggplot(merge, aes(x = limix_beta_se, y = se_matrixeqtl, colour = minP)) + 
  geom_point(alpha = ifelse(merge$minP == 1, 1, 0.1)) + 
  geom_abline(slope = 1, intercept = 0) +
  labs(x = "Limix SE", y = "Matrix eQTL SE")

m3 <- ggplot(merge, aes(x = limix_zscore, y = -statistic, colour = minP)) + 
  geom_point(alpha = ifelse(merge$minP == 1, 1, 0.1)) + 
  geom_abline(slope = 1, intercept = 0) +
  labs(x = "Limix Z-score", y = "Matrix eQTL Z-score")

m4 <- ggplot(merge, aes(x = -log10(limix_p_value), y = -log10(FDR), colour = minP)) + 
  geom_point(alpha = ifelse(merge$minP == 1, 1, 0.1)) + 
  geom_abline(slope = 1, intercept = 0) +
  labs(x = "-log10 Limix Locally-corrected P-value", y = "-log10 Matrix eQTL FDR-adjusted P-value")

save_plot(plot_grid(m1,m2,m3, m4), file = "cis_eqtls_18373genes_age_sex_rin_batch_PC10_PEER20/processed/cis_eqtls_18373genes_age_sex_rin_batch_PC10_PEER20_eSNPs_matrixeQTL_comparison.png", base_height = 7, base_width = 12)












eSNP_mismatch <- filter(eSNP_replication, limix_feature_id %in% flexmismatchGenes) %>%
  mutate(matchid2 = paste0(limix_gene_name, ":", limix_rsid))
c22_mateqtl <- c22_mateqtl %>% 
  mutate(matchid2 = paste0(gene, ":", snps))
eSNP_mismatch <- inner_join(eSNP_mismatch, c22_mateqtl) 
eSNP_mismatch <- eSNP_mismatch %>%
  mutate(se_matrixeqtl = beta / statistic)

m1 <- ggplot(eSNP_mismatch, aes(x = limix_beta, y = beta, colour = gene)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) + 
  facet_wrap(. ~ gene)

m2 <- ggplot(eSNP_mismatch, aes(x = limix_beta_se, y = limix_beta_se, colour = gene)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) + 
  facet_wrap(. ~ gene)


m3 <- ggplot(eSNP_mismatch, aes(x = limix_zscore, y = statistic, colour = gene)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) + 
  facet_wrap(. ~ gene)
plot_grid(m1,m2,m3)