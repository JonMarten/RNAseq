library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)
#library(hexbin)
theme_set(theme_cowplot())

setwd("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/results/")

e22 <- fread("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/eQTLgen/cis-eQTLs_chr22.txt", data.table = F)

jobdf <- data.frame("job_id" = c("cis_eqtls_18373genes_age_sex_rin_batch_PC10", "cis_eqtls_18373genes_age_sex_rin_batch_PC10_PEER20", "cis_eqtls_18373genes_age_sex_rin_batch_PC10_NeutPCT_LympPCT_MonoPCT_EoPCT_BasoPCT"))

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
         Zscore_correlation_all = NA,
         Zscore_correlation_sig = NA,
         eSNP_BH_thresh = NA)

files <- paste0(jobdf$job_id,"/processed/",jobdf$job_id, "_results_merged.txt")

for(i in c(1:length(files))) {

  c22 <- fread(files[i], data.table = F)
  
  c22 <- distinct(c22)
  
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
  
  # Flip alleles to match eQTLgen and check direction mismatch
  flips <- merge$limix_assessed_allele != merge$eQTLgen_AssessedAllele
  merge <- merge %>%
    mutate(limix_zscore = limix_beta / limix_beta_se) %>%
    mutate(eQTLgen_Zscore_flipped = ifelse(flips, -eQTLgen_Zscore, eQTLgen_Zscore),
           flipped = ifelse(flips, 1, 0)) %>%
    mutate(signMismatch = ifelse(sign(limix_zscore) != sign(eQTLgen_Zscore_flipped), 1, 0))
  
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
           sigrep_BH = ifelse(eQTLgen_padj_BH < 0.05, 1, 0)) %>%
    mutate(sig_bonf_sign = ifelse(signMismatch == 0 & eQTLgen_Pvalue < 0.05/nrow(eSNP_replication), 1, 0))
  jobdf$eSNPs_replicate[i] <- which(eSNP_replication$sig_bonf_sign == 1) %>% length
  
  # correlate Z-scores
  jobdf$Zscore_correlation_all[i] <- cor(eSNP_replication$limix_zscore, eSNP_replication$eQTLgen_Zscore_flipped)
  jobdf$Zscore_correlation_sig[i] <- cor(eSNP_replication$limix_zscore[which(eSNP_replication$sig_bonf_sign == 1)],
                                         eSNP_replication$eQTLgen_Zscore_flipped[which(eSNP_replication$sig_bonf_sign == 1)])
  
  # Write out merged eQTLgen replication
  outname_eqtgen <- paste0(jobdf$job_id[i],"/processed/",jobdf$job_id[i], "_eSNPs_eqtlgenReplication.txt")
    fwrite(eSNP_replication, file = outname_eqtgen)
  
  # Pull genes with mismatched SNP Zscores
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
  
  #Plot Z scores
  g <- ggplot(eSNP_replication, aes(x = limix_zscore, y = eQTLgen_Zscore_flipped, color = Significant)) + 
    geom_point(cex = 0.8, alpha = 0.5, pch = 20) +
    theme(aspect.ratio = 1, legend.position = "right") +
    scale_colour_manual(values = c("gray75", "orangered", "chartreuse3")) +
    geom_vline(xintercept = 0) + 
    geom_hline(yintercept = 0) + 
    labs(x = "Limix Z-score", y = "eQTLgen Z-score") 
  
  outname_plot <- paste0(jobdf$job_id[i],"/processed/",jobdf$job_id[i], "_eSNPs_eqtlgenReplication_zscoreplot.png")
  save_plot(file = outname_plot, g, base_height = 10, base_width = 10)
  

    # plot mismatches coloured by significance 
    g2 <- ggplot(filter(eSNP_replication, limix_feature_id %in% flexmismatchGenes), 
                aes(x = limix_zscore, y = eQTLgen_Zscore_flipped, color = Significant)) + 
      geom_vline(xintercept = 0) + 
      geom_hline(yintercept = 0) + 
      geom_point(cex = 0.8, alpha = 1, pch = 20) +
      theme(aspect.ratio = 1, legend.position = "right") +
      labs(x = "Limix Z-score", y = "eQTLgen Z-score") +
      theme_minimal() +
      scale_colour_manual(values = c("gray75", "orangered", "chartreuse3")) +
      facet_wrap(. ~ limix_gene_name)  
    
    outname_plot2 <- paste0(jobdf$job_id[i],"/processed/",jobdf$job_id[i], "_eSNPs_eqtlgenReplication_zscoreplot_mismatch.png")
    save_plot(file = outname_plot2, g2, base_height = 7, base_width = 12)
    
    # Plot mistmatches coloured by strand-ambiguity
    g3 <- ggplot(filter(eSNP_replication, limix_feature_id %in% flexmismatchGenes), 
                 aes(x = limix_zscore, y = eQTLgen_Zscore_flipped, color = as.factor(ambig))) + 
      geom_vline(xintercept = 0) + 
      geom_hline(yintercept = 0) + 
      geom_point(cex = 0.8, alpha = 1, pch = 20) +
      theme(aspect.ratio = 1, legend.position = "right") +
      labs(x = "Limix Z-score", y = "eQTLgen Z-score") +
      theme_minimal() +
      scale_colour_manual(values = c("gray75", "orangered", "chartreuse3")) +
      facet_wrap(. ~ limix_gene_name)  
    
    outname_plot3 <- paste0(jobdf$job_id[i],"/processed/",jobdf$job_id[i], "_eSNPs_eqtlgenReplication_zscoreplot_mismatch_ambig.png")
    save_plot(file = outname_plot3, g3, base_height = 7, base_width = 12)
    
    # Plot mismatches coloured by eQTLgen number of cohorts size
    g4<- ggplot(filter(eSNP_replication, limix_feature_id %in% flexmismatchGenes), 
                 aes(x = limix_zscore, y = eQTLgen_Zscore_flipped, color = eQTLgen_NrCohorts)) + 
      geom_vline(xintercept = 0) + 
      geom_hline(yintercept = 0) + 
      geom_point(cex = 0.8, alpha = 1, pch = 20) +
      theme(aspect.ratio = 1, legend.position = "right") +
      labs(x = "Limix Z-score", y = "eQTLgen Z-score") +
      theme_minimal() +
      scale_colour_viridis_c() +
      facet_wrap(. ~ limix_gene_name)  
    
    outname_plot4 <- paste0(jobdf$job_id[i],"/processed/",jobdf$job_id[i], "_eSNPs_eqtlgenReplication_zscoreplot_mismatch_eQTLnCohorts.png")
    save_plot(file = outname_plot4, g4, base_height = 7, base_width = 12)
    
}

fwrite(jobdf, file = "peer_comparison_summary.csv")
