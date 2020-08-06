library(dplyr)
library(data.table)
library(ggplot2)
library(viridis)
library(cowplot)
theme_set(theme_cowplot())
setwd("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/04_phase2_full_analysis/results/cis_eQTLs/")

tqtl <- fread("tensorqtl_cis_MAF0.005_cisNominal_chr22.csv", data.table = F)

tqtl <- tqtl %>% 
  mutate(matchID = paste0(phenotype_id, ":", variant_id))

plink_map <- fread("../../genotypes/INTERVAL_RNAseq_Phase1-2_imputed_b38_biallelic_MAF0.005_chr22.bim" , data.table = F)
names(plink_map) <- c("chr", "variant_id", "morg", "pos", "A1","A2")
tqtl2 <- left_join(tqtl, plink_map)

eQTLgen <- fread("zcat /rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/eQTLgen/cis-eQTLs_full_20180905.txt.gz", data.table = F)
eQTLgen2 <- eQTLgen %>% 
  filter(SNPChr == 22) %>%
  mutate(matchID = paste0(Gene, ":", SNP))
m <- left_join(tqtl2, eQTLgen2, by = "matchID")

m2 <- m %>% 
  filter(Pvalue < 5e-8) %>%
  filter(!is.na(SNPChr) | !is.na(SNPChr)) %>%
  mutate(flipE = ifelse(A1 == AssessedAllele, 0, 1)) %>%
  mutate(eQTLgenZFlip = ifelse(flipE == 1, -Zscore, Zscore)) %>%
  mutate(tensorZ = slope / slope_se,
         eQTLgenZ = eQTLgenZFlip) %>%
  mutate(tensorZ.norm = (tensorZ - mean(tensorZ, na.rm=T))/sd(tensorZ, na.rm = T),
         eQTLgenZ.norm = (eQTLgenZ - mean(eQTLgenZ, na.rm=T))/sd(eQTLgenZ, na.rm = T))

out <- m2 %>%
  select(phenotype_id,
         gene = GeneSymbol,
         variant_id,
         chr,
         pos_b38 = pos,
         tss_distance,
         assessed_allele = A1,
         other_allele = A2,
         tqtl_maf = maf,
         tqtl_pval_nominal = pval_nominal,
         tqtl_beta = slope,
         tqtl_beta_se = slope_se,
         tqtl_zscore = tensorZ,
         tqtl_zscore_norm = tensorZ.norm,
         eqtlgen_zscore = eQTLgenZ,
         eqtlgen_zscore_norm = eQTLgenZ.norm,
         eqtlgen_pvalue = Pvalue,
         eqtlgen_fdr_pvalue = FDR,
         eqtlgen_num_cohorts = NrCohorts,
         eqtlgen_num_samples = NrSamples)
         
fwrite(out, file = "tensorqtl_cis_MAF0.005_cis_chr22_eQTLgen_comparison.csv")

ET <- ggplot(out, aes(y =  tqtl_zscore_norm, x = eqtlgen_zscore_norm)) +
  geom_point(pch = 20, cex = 0.8, alpha = 0.01) +
  labs(x = "eQTLgen Zscore", y = "TensorQTL Zscore") +
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0)         
         
ETdens <- ggplot(out, aes(y =  tqtl_zscore_norm, x = eqtlgen_zscore_norm)) +
  geom_point(pch = 20, cex = 0.8, alpha = 0.1) +
  stat_density_2d(geom = "polygon", colour="white") +
  labs(x = "eQTLgen Zscore", y = "TensorQTL Zscore") +
  scale_fill_continuous(type = "viridis") +
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0)         

ggsave(ETdens, file = "chr22_eqtlgen_zscore_comparison.png")
