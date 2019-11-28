library(data.table)
library(dplyr)

setwd("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/results/cis_eqtls_18373genes_age_sex_rin_batch_PC10_PEER20_5GenesPerChunk/processed")
prefix <- "cis_eqtls_18373genes_age_sex_rin_batch_PC10_PEER20"

dat <- fread("output_results_merged.txt", data.table = F)
dat <- distinct(dat)
dat <- dat %>% 
  filter(!is.na(snp_chromosome))

#### Multple testing correction
# Get SNP with minimum corrected (at gene level) p-value for each feature
top_snps_per_gene <- dat %>% 
  group_by(feature_id) %>%
  summarise(minP = min(empirical_feature_p_value)) %>%
  data.frame() %>%
  mutate(minp_adjusted_BH = p.adjust(minP, method = "BH")) %>%
  mutate(sig = ifelse(minp_adjusted_BH < 0.05, 1, 0))

eGenes.BH <- top_snps_per_gene %>%
  filter(sig == 1) %>%
  pull(feature_id) 

# Identify eSNPs within significant eGenes by selecting all SNPs with locally-adjusted P lower than the p-value threshold that corresponds to a global BH-adjusted P-value of 0.05
  # from https://rdrr.io/bioc/IHW/src/R/helpers.R Function to work out what the theshold is for rejecting a p-value
get_bh_threshold <- function(pvals, alpha, mtests = length(pvals)){
  m <- length(pvals)
  pvals <- sort(pvals)
  prejected <- which(pvals <= (1:m)/mtests*alpha)
  ifelse(length(prejected) == 0, 0, pvals[prejected[which.max(prejected)]])
}

eSNPThresh <- get_bh_threshold(top_snps_per_gene$minp_adjusted_BH, 0.05)
cat(paste0("\neSNP locally-corrected P-value threshold for significance: ", eSNPThresh))

eSNPs <- dat %>% 
  filter(feature_id %in% eGenes.BH & empirical_feature_p_value < eSNPThresh)

outname.eSNP <- paste0(prefix, "_eSNPs.txt")
fwrite(eSNPs, file = outname.eSNP)

# Plot results
library(qqman)
manhattan(eSNPs, chr = "snp_chromosome", bp = "snp_position", p = "empirical_feature_p_value")
