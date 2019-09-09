library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)

setwd("/rds/user/jm2294/rds-jmmh2-projects/interval_rna_seq/analysis/00_testing/test_parameters")

e22 <- fread("eQTLgen/cis-eQTLs_chr22.txt", data.table = F)

jobdf <- fread(data.table = F, "job_ids.txt")
jobdf$taskID <- as.character(jobdf$taskID)
jobdf <- jobdf %>% filter(!is.na(taskID) & Features == "Coding")
jobids <- jobdf$taskID

jobdf$fileprefix <- paste0("results_merged_chr22_window", jobdf$Window,"_Perm", jobdf$Permutations,"_MAF",jobdf$MAF)

jobdf <- jobdf %>%
  mutate(num_SNPs_tested = NA,
         num_features_tested = NA,
         num_SNPs_overlap = NA,
         num_features_overlap = NA,
         eGenes_num = NA,
         eGenes_overlap = NA,
         eGenes_replicate = NA, 
         eSNPs_num = NA,
         eSNPs_overlap = NA,
         eSNPs_replicate = NA, 
         Zscore_correlation = NA,
         eSNP_BH_thresh = NA)

files <- paste0(jobdf$fileprefix, ".txt")

for(i in c(1:5,8)) {

  c22 <- fread(files[i], data.table = F)
  
  #### Multple testing correction
  # Get SNP with minimum corrected (at gene level) p-value for each feature
  top_snps_per_gene <- c22 %>% 
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
  jobdf$eSNP_BH_thresh[i] <- eSNPThresh
  
  eSNPs <- c22 %>% 
    filter(feature_id %in% eGenes.BH & empirical_feature_p_value < eSNPThresh)
  
  outname <- paste0(jobdf$fileprefix[i], "_eSNPs.txt")
  fwrite(eSNPs, file = outname)
  
  # Get table metrics
  jobdf$num_SNPs_tested[i] <- nrow(c22)
  jobdf$num_features_tested[i] <- unique(c22$feature_id) %>% length
  jobdf$eSNPs_num[i] <- nrow(eSNPs)
  jobdf$eGenes_num[i] <- length(eGenes.BH)
  
  # merge with eQTLgen
  c22m <- c22
  names(c22m) <- paste0("limix_", names(c22m))
  c22m <- c22m %>%
    mutate(matchID = paste0(limix_feature_id, "_", limix_rsid))
  merge <- inner_join(c22m, e22)
  
  jobdf$num_SNPs_overlap[i] <- nrow(merge)
  jobdf$num_features_overlap[i] <- unique(merge$limix_feature_id) %>% length
  
  jobdf$eGenes_overlap[i] <- which(eGenes.BH %in% merge$limix_feature_id) %>% length
  
  # check which of our eGenes are significant in eqtlgen
  eGene_replication <- e22 %>% 
    group_by(eQTLgen_Gene) %>%
    filter(eQTLgen_Gene %in% eGenes.BH) %>%
    summarise(minP = min(eQTLgen_FDR)) %>% 
    arrange(minP) %>%
    data.frame() %>%
    mutate(significant = ifelse(minP < 0.05 / jobdf$eGenes_overlap[i], 1, 0))
  
  jobdf$eGenes_replicate[i] <- which(eGene_replication$significant == 1) %>% length
  
  # Check which of our eSNPs are signficant in eqtlgen
  eSNPs <- eSNPs %>%
    mutate(matchID = paste0(feature_id, "_", rsid))
  eSNP_replication <- merge %>%
    filter(matchID %in% eSNPs$matchID)
  jobdf$eSNPs_overlap[i] <- nrow(eSNP_replication)
  
  eSNP_replication <- eSNP_replication %>%
    mutate(eQTLgen_padj_BH = p.adjust(eQTLgen_Pvalue, method = "BH")) %>%
    mutate(sigrep_bonf = ifelse(eQTLgen_Pvalue < 0.05/nrow(eSNP_replication), 1, 0),
           sigrep_BH = ifelse(eQTLgen_padj_BH < 0.05, 1, 0))
  jobdf$eSNPs_replicate[i] <- which (eSNP_replication$sigrep_BH == 1) %>% length
  
  # correlate Z-scores, flip z-score signs for SNPs with different effect alleles
  flips <- eSNP_replication$limix_assessed_allele != eSNP_replication$eQTLgen_AssessedAllele
  eSNP_replication <- eSNP_replication %>%
    mutate(limix_zscore = limix_beta / limix_beta_se) %>%
    mutate(eQTLgen_Zscore_flipped = ifelse(flips, -eQTLgen_Zscore, eQTLgen_Zscore))
  jobdf$Zscore_correlation[i] <- cor(eSNP_replication$limix_zscore, eSNP_replication$eQTLgen_Zscore_flipped)

  # Write out merged eQTLgen replication
  outname_eqtgen <- paste0(jobdf$fileprefix[i], "_eSNPs_eqtlgenReplication.txt")
  fwrite(eSNPs, file = outname_eqtgen)
  
  # Pull genes with mismatched SNP Zscores
  eSNP_replication <- eSNP_replication %>%
    mutate(signMismatch = ifelse(sign(limix_zscore) != sign(eQTLgen_Zscore_flipped), 1, 0))
  
  mismatchGenes <- eSNP_replication %>%
    filter(signMismatch == 1) %>%
    pull(limix_feature_id) %>%
    unique()
  
  g <- ggplot(eSNP_replication, aes(x = limix_zscore, y = eQTLgen_Zscore_flipped)) + 
    geom_point(colour = "gray50", cex = 0.5, alpha = 0.8, pch = 20) +
    geom_vline(xintercept = 0) + 
    geom_hline(yintercept = 0) + 
    theme(aspect.ratio = 1) +
    labs(x = "Limix Z-score", y = "eQTLgen Z-score")
  outname_plot <- paste0(jobdf$fileprefix[i], "_eSNPs_eqtlgenReplication_zscoreplot.png")
  save_plot(file = outname_plot, g)
}

fwrite(jobdf, file = "parameter_comparison_summary.csv")
####



