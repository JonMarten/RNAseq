library(data.table)
library(dplyr)

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
         Zscore_correlation = NA)


files <- paste0(jobdf$fileprefix, ".txt")

for(i in 1:5) {

  dat3 <- fread(files[i], data.table = F)
  
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
    pull(feature_id) %>%
    unique()
  
  sigTops.BH <- tops %>%
    filter(p_adjusted_BH < 0.05)
  eGenes.BH <- sigTops.BH %>%
    pull(feature_id) %>%
    unique()
  
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
  
  outname <- paste0(jobdf$fileprefix[i], "_eSNPs.txt")
  fwrite(eSNPs, file = outname)
  
  
  # Get table metrics
  jobdf$num_SNPs_tested[i] <- nrow(dat3)
  jobdf$num_features_tested[i] <- nrow(dat3)
  
  jobdf$eSNPs_num[i] <- nrow(eSNPs)
  jobdf$eGenes_num[i] <- length(eGenes.BH)
  
  
  
  
  
  
  jobdf$eSNPs_overlap[i] <- nrow(merge)
  jobdf$eGenes_overlap[i] <- merge %>% pull(limix_feature_id) %>% unique %>% length
  
}

####



