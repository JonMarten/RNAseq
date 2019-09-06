# Merge and filter parameter test results

library(data.table)
library(dplyr)
library(stringr)

setwd("/rds/user/jm2294/rds-jmmh2-projects/interval_rna_seq/analysis/00_testing/test_parameters")

jobdf <- fread(data.table = F, "job_ids.txt")
jobdf$taskID <- as.character(jobdf$taskID)
jobdf <- jobdf %>% filter(!is.na(taskID) & Features == "Coding")
jobids <- jobdf$taskID

for(i in 1:5){

  resdir <- paste0("/home/jm2294/rds/rds-jmmh2-projects/interval_rna_seq/analysis/00_testing/results/eqtl_test_parameters_fulchrt",jobids[i],"/processed")
  setwd(resdir)
  files <- list.files()
  files <- files[grepl("processed_qtl_results", files)]
  
  # Read in limix results
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
  
  if(length(datSplit) < 1) {
    # Save list of dropped features
    print("Filtering variants with alpha < 0.1 or > 10")
    droppedFeatures <- datSplit[[2]] %>% 
      pull(gene_name) %>% 
      unique
    write.table(droppedFeatures, row.names = F, col.names = F, quote = F, file = "dropped_features.txt")
  } else {
    print("No variants with alpha < 0.1 or > 10")
  }
  
  # Remove dropped features and reformat results
  dat2 <- datSplit[[1]] %>%
    select(-fail) %>%
    data.frame
  
  rsid <- str_split_fixed(dat2$snp_id, ":", 2)
  dat2$rsid <- rsid[,2]
  
  dat3 <- dat2 %>% 
    arrange(feature_id) %>%
    select(feature_id, feature_chromosome:beta_param, snp_id, rsid, snp_chromosome:hwe_p, beta, beta_se, p_value, empirical_feature_p_value)
  setwd("/rds/user/jm2294/rds-jmmh2-projects/interval_rna_seq/analysis/00_testing/test_parameters")
  fwrite(dat3, quote = F, file = paste0("results_merged_chr22_window", jobdf$Window[i],"_Perm", jobdf$Permutations[i],"_MAF",jobdf$MAF[i], ".txt"))
  
  dat3$taskID <- jobids[i]
  rm(files, dat, dat2, rsid, datSplit, dat3, resdir)
}