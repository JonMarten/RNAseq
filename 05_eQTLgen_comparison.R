library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)
library(hexbin)
theme_set(theme_cowplot())

setwd("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/results/")

e22 <- fread("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/eQTLgen/cis-eQTLs_chr22.txt", data.table = F)

#jobdf <- fread(data.table = F, "job_ids.txt")
#jobdf$taskID <- as.character(jobdf$taskID)
#jobdf <- jobdf %>% filter(!is.na(taskID) & Features == "Coding")
#jobids <- jobdf$taskID
#jobdf$fileprefix <- paste0("results_merged_chr22_window", jobdf$Window,"_Perm", jobdf$Permutations,"_MAF",jobdf$MAF)

jobdf <- data.frame("job_id" = c("cis_eqtls_18373genes_age_sex_rin_batch_PC10", "cis_eqtls_18373genes_age_sex_rin_batch_PC10_PEER20", "cis_eqtls_18373genes_age_sex_rin_batch_PC10_NeutPCT_MonoPCT_EoPCT_BasoPCT"))

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

files <- paste0(jobdf$job_id,"/processed/",jobdf$job_id, "_results_merged.txt")

for(i in c(1:length(files))) {

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
  
  outname <- paste0(jobdf$job_id[i],"/processed/",jobdf$job_id[i], "_eSNPs.txt")
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
    mutate(eQTLgen_Zscore_flipped = ifelse(flips, -eQTLgen_Zscore, eQTLgen_Zscore),
           flipped = ifelse(flips, 1, 0))
  jobdf$Zscore_correlation[i] <- cor(eSNP_replication$limix_zscore, eSNP_replication$eQTLgen_Zscore_flipped)

  # Write out merged eQTLgen replication
  outname_eqtgen <- paste0(jobdf$job_id[i],"/processed/",jobdf$job_id[i], "_eSNPs_eqtlgenReplication.txt")
    fwrite(eSNP_replication, file = outname_eqtgen)
  
  # Pull genes with mismatched SNP Zscores
  eSNP_replication <- eSNP_replication %>%
    mutate(signMismatch = ifelse(sign(limix_zscore) != sign(eQTLgen_Zscore_flipped), 1, 0))
  
  mismatchGenes <- eSNP_replication %>%
    filter(signMismatch == 1) %>%
    pull(limix_feature_id) %>%
    unique()
  
  eSNP_replication <- eSNP_replication %>%
    mutate(Significant = ifelse(sigrep_bonf == 1 & sigrep_BH == 0, "Bonf. only", 
                                ifelse(sigrep_bonf == 0 & sigrep_BH == 1, "BH only", 
                                       ifelse(sigrep_bonf == 1 & sigrep_BH == 1, "BH & Bonf.", "NS")))) %>%
    mutate(Significant = factor(Significant, levels = c("NS","BH only", "BH & Bonf.")))
  
  
  eSNP_replication$limix_feature_length <- eSNP_replication$limix_feature_end - eSNP_replication$limix_feature_start
  eSNP_replication <- eSNP_replication %>%
    mutate(ambig = ifelse(eQTLgen_AssessedAllele == "A" & eQTLgen_OtherAllele == "T" |
                            eQTLgen_AssessedAllele == "T" & eQTLgen_OtherAllele == "A" |
                            eQTLgen_AssessedAllele == "C" & eQTLgen_OtherAllele == "G" |
                            eQTLgen_AssessedAllele == "G" & eQTLgen_OtherAllele == "C", 1, 0)) %>%
    mutate(veryambig = ifelse(ambig == 1 & limix_maf > 0.45, 1, 0))
  
# get list of features with more than two SNPs mismatched  
  flexmismatchGenes <- eSNP_replication %>%
    filter(signMismatch == 1) %>%
    group_by(limix_feature_id) %>%
    summarise(n_mismatched = n()) %>%
    data.frame %>%
    filter(n_mismatched > 2) %>%
    pull(limix_feature_id)
  
    g <- ggplot(eSNP_replication, aes(x = limix_zscore, y = eQTLgen_Zscore_flipped, color = eQTLgen_NrSamples)) + 
      geom_vline(xintercept = 0) + 
      geom_hline(yintercept = 0) + 
      geom_point(cex = 0.8, alpha = 0.01, pch = 20) +
      theme(aspect.ratio = 1, legend.position = "right") +
      labs(x = "Limix Z-score", y = "eQTLgen Z-score") 
    
    
  outname_plot <- paste0(jobdf$job_id[i],"/processed/",jobdf$job_id[i], "_eSNPs_eqtlgenReplication_zscoreplot.png")
  save_plot(file = outname_plot, g, base_height = 7, base_width = 6)
  

    # plot mismatches
    g2 <- ggplot(filter(eSNP_replication, limix_feature_id %in% flexmismatchGenes), 
                aes(x = limix_zscore, y = eQTLgen_Zscore_flipped, color = eQTLgen_NrSamples)) + 
      geom_vline(xintercept = 0) + 
      geom_hline(yintercept = 0) + 
      geom_point(cex = 0.8, alpha = 1, pch = 20) +
      theme(aspect.ratio = 1, legend.position = "right") +
      labs(x = "Limix Z-score", y = "eQTLgen Z-score") +
      theme_minimal() +
      scale_color_viridis_c(option = "C") +
      facet_wrap(. ~ limix_gene_name)  
    
    outname_plot2 <- paste0(jobdf$job_id[i],"/processed/",jobdf$job_id[i], "_eSNPs_eqtlgenReplication_zscoreplot_mismatch.png")
    save_plot(file = outname_plot2, g2, base_height = 7, base_width = 12)
}

#mem <- fread("parameter_comparison_joblogs.csv", data.table = F)
#mem <- mem %>% select(-c(Cohort:Permutations)) %>%
#  mutate(taskID = as.character(taskID))
#jobdf2 <- full_join(mem, jobdf)

fwrite(jobdf, file = "parameter_comparison_summary.csv")

# Look at TSS distance

w250 <- fread("results_merged_chr22_window250000_Perm100_MAF0.01_eSNPs_eqtlgenReplication.txt", data.table = F)
w250$window <- "250kb"
w500 <- fread("results_merged_chr22_window500000_Perm100_MAF0.01_eSNPs_eqtlgenReplication.txt", data.table = F)
w500$window <- "500kb"
w1000 <- fread("results_merged_chr22_window1000000_Perm100_MAF0.01_eSNPs_eqtlgenReplication.txt", data.table = F)
w1000$window <- "1Mb"


wAll <- rbind(w250, w500, w1000)

wAll <- wAll %>%
  mutate(eSNP_distance_from_tss = ifelse(limix_feature_strand == 1, 
                                         limix_snp_position - limix_feature_start,
                                         limix_feature_end - limix_snp_position),
         window = factor(window, levels = c("250kb", "500kb", "1Mb")))
ggplot(wAll, aes(x = eSNP_distance_from_tss, fill = window)) +
  geom_density(alpha = 0.5) +
  facet_wrap(. ~ window, ncol = 1)






####

# Test multiple testing scenarios

# mt <- fread("results_merged_chr22_window500000_Perm100_MAF0.01.txt", data.table = F)
# 
# mt_topsnps <- mt %>% 
#   group_by(feature_id) %>%
#   summarise(minP = min(empirical_feature_p_value)) %>%
#   data.frame() %>%
#   mutate(minp_adjusted_BH = p.adjust(minP, method = "BH")) %>%
#   mutate(sig.BH = ifelse(minp_adjusted_BH < 0.05, 1, 0))
# mt_topsnps <- mt_topsnps %>%
#   mutate(sig.bonf = ifelse(minP < 0.05/nrow(mt_topsnps),1,0))
# 
# eGenes.BH <- mt_topsnps %>% filter(sig.BH == 1) %>% pull(feature_id)
# eGenes.bonf <- mt_topsnps %>% filter(sig.bonf == 1) %>% pull(feature_id)
# 
# eSNPs.bonf.bonf <- mt %>% filter(feature_id %in% eGenes.bonf) %>%
#   mutate(sig = ifelse(empirical_feature_p_value < 0.05/length(eGenes.bonf), 1, 0))
# eSNPs.BH.bonf <- mt %>% filter(feature_id %in% eGenes.BH) %>%
#   mutate(sig = ifelse(empirical_feature_p_value < 0.05/length(eGenes.BH), 1, 0))
# 
# 
# thresh.bh
# 
# 
# eSNPThresh <- get_bh_threshold(top_snps_per_gene$minp_adjusted_BH, 0.05)
# jobdf$eSNP_BH_thresh[i] <- eSNPThresh
# 
# eSNPs <- c22 %>% 
#   filter(feature_id %in% eGenes.BH & empirical_feature_p_value < eSNPThresh)



