# Merge results chunks and post-process

setwd( "U:/Projects/RNAseq/analysis/00_testing/results/eqtl_test_13633258/processed")
library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)

files <- list.files()
files <- files[grepl("processed_qtl_results", files)]

for (file in files) {
  
  # if the merged dataset doesn't exist, create it
  if (!exists("dat")) {
    dat <- fread(file, data.table = F)
  }
  
  # if the merged dataset does exist, append to it
  if (exists("dat")) {
    tempDat <- fread(file, data.table = F)
    dat <- rbind(dat, tempDat)
    rm(tempDat)
  }
  
}

# From Marc: sometimes the local gene-based corretion fails, usually because the genes have low variation in the dataset. 
# These have large or small alpha values 
# He proposes setting the empirical p-values back to the original but I prefer just to remove these as they're unlikely to be informative.
datSplit <- dat %>%
  mutate(fail = ifelse((alpha_param < 0.1 | alpha_param > 10), 1, 0)) %>%
  group_by(fail) %>%
  group_split()

# Save list of dropped features
droppedFeatures <- datSplit[[2]] %>% 
  pull(gene_name) %>% 
  unique
write.table(droppedFeatures, row.names = F, col.names = F, quote = F, file = "dropped_features.txt")


# Remove dropped features and reformat results
dat2 <- datSplit[[1]] %>%
  select(-fail) %>%
  data.frame

rsid <- str_split_fixed(dat2$snp_id, ":", 2)
dat2$rsid <- rsid[,2]

dat3 <- dat2 %>% 
  arrange(feature_id) %>%
  select(feature_id, feature_chromosome:beta_param, snp_id, rsid, snp_chromosome:hwe_p, beta, beta_se, p_value, empirical_feature_p_value)

fwrite(dat3, quote = F, file = "results_merged_chr22.txt")

rm(files, dat, dat2, rsid, datSplit)

#dat3 <- fread("results_merged_chr22.txt", data.table = F)

#### Multple testing correction
# Get SNP with minimum corrected (at gene level) p-value for each feature
tops <- dat3 %>% 
  group_by(feature_id) %>%
  filter(empirical_feature_p_value == min(empirical_feature_p_value)) %>%
  data.frame %>%
  mutate(p_adjusted_BH = p.adjust(empirical_feature_p_value, method = "BH"))

bonfThresh <- 0.05 / nrow(tops)

sigTops.bonf <- tops %>%
  filter(empirical_feature_p_value < bonfThresh)
eGenes.bonf <- sigTops.bonf %>% 
  pull(feature_id)

sigTops.BH <- tops %>%
  filter(p_adjusted_BH < 0.05)
eGenes.BH <- sigTops.BH %>%
  pull(feature_id)

# Identify eSNPs within significant eGenes by selecting all SNPs with locally-adjusted P lower than the p-value threshold that corresponds to a global BH-adjusted P-value of 0.05

# from https://rdrr.io/bioc/IHW/src/R/helpers.R Function to work out what the theshold is for rejecting a p-value
get_bh_threshold <- function(pvals, alpha, mtests = length(pvals)){
  m <- length(pvals)
  pvals <- sort(pvals)
  prejected <- which(pvals <= (1:m)/mtests*alpha)
  ifelse(length(prejected) == 0, 0, pvals[prejected[which.max(prejected)]])
}

eSNPThresh <- get_bh_threshold(tops$empirical_feature_p_value, 0.05)

eSNPs <- dat3 %>% 
  filter(feature_id %in% eGenes.BH & empirical_feature_p_value < eSNPThresh)

fwrite(eSNPs, file = "significant_eSNPs_BH.csv")
###############

# Compare to eQTLgen
eQTL <- fread("U:/Projects/RNAseq/eQTLgen/cis-eQTLs_full_20180905.txt", data.table = F)
names(eQTL) <- paste0("eQTLgen_", names(eQTL))
e22 <- eQTL %>% 
  filter(eQTLgen_SNPChr == 22) %>%
  mutate(matchID = paste0(eQTLgen_Gene, "_", eQTLgen_SNP))
dat4 <- dat3
names(dat4) <- paste0("limix_", names(dat4))
dat4 <- dat4 %>%
  mutate(matchID = paste0(limix_feature_id, "_", limix_rsid))
merge <- inner_join(dat4, e22)

mplot <- runif(3000, 1, nrow(merge))

ggplot(merge[mplot,], aes(x = -log10(limix_p_value), y = -log10(eQTLgen_Pvalue))) + geom_point()

cor(merge$limix_p_value, merge$eQTLgen_Pvalue,method = "spearman")

# Check our hits in eqtlgen
# eGenes
eGene_replication <- e22 %>% 
  group_by(eQTLgen_Gene) %>%
  filter(eQTLgen_Gene %in% eGenes.BH) %>%
  summarise(minP = min(eQTLgen_Pvalue)) %>% 
  arrange(minP) %>%
  data.frame() 

# eSNPs
eSNPs <- eSNPs %>%
  mutate(matchID = paste0(feature_id, "_", rsid))
eSNP_replication <- merge %>%
  filter(matchID %in% eSNPs$matchID)

# Count significant SNPs after bonf
eSNP_replication %>%
  filter(eQTLgen_Pvalue < 0.05/nrow(eSNP_replication)) %>%
  nrow
# Count significant SNPs after BH
eSNP_replication <- eSNP_replication %>%
  mutate(p_adjusted_BH = p.adjust(eQTLgen_Pvalue, method = "BH")) 
 
eSNP_replication <- eSNP_replication %>%
  mutate(limix_zscore = limix_beta / limix_beta_se) %>%
  mutate(sigRep = ifelse(eQTLgen_Pvalue < 0.05/nrow(eSNP_replication), 1, 0))

# flip z-score signs for SNPs with different effect alleles
flips <- eSNP_replication$limix_assessed_allele != eSNP_replication$eQTLgen_AssessedAllele
eSNP_replication <- eSNP_replication %>%
  mutate(eQTLgen_Zscore_flipped = ifelse(flips, -eQTLgen_Zscore, eQTLgen_Zscore))

cor(eSNP_replication$limix_zscore, eSNP_replication$eQTLgen_Zscore_flipped)

# Pull genes with mismatched SNP Zscores
eSNP_replication <- eSNP_replication %>%
  mutate(signMismatch = ifelse(sign(limix_zscore) != sign(eQTLgen_Zscore_flipped), 1, 0))

mismatchGenes <- eSNP_replication %>%
  filter(signMismatch == 1) %>%
  pull(limix_feature_id) %>%
  unique()

# Plot just mismatched genes
ggplot(filter(eSNP_replication, limix_feature_id %in% mismatchGenes), aes(x = limix_zscore, y = eQTLgen_Zscore_flipped, colour = limix_feature_id)) + 
  geom_point() +
  geom_vline(xintercept = 0) + 
  geom_hline(yintercept = 0) + 
  theme(aspect.ratio = 1)

# Plot all, highlight mismatches
ggplot(eSNP_replication, aes(x = limix_zscore, y = eQTLgen_Zscore_flipped)) + 
  geom_point(colour = "gray50") +
  geom_vline(xintercept = 0) + 
  geom_hline(yintercept = 0) + 
  geom_point(data = filter(eSNP_replication, limix_feature_id %in% mismatchGenes), aes(x = limix_zscore, y = eQTLgen_Zscore_flipped, colour = limix_feature_id))  +
  theme(aspect.ratio = 1)

