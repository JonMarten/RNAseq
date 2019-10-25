# Merge and filter postprocessed limix output
library(data.table, warn.conflicts = FALSE)
library(dplyr, warn.conflicts = FALSE)
library(stringr, warn.conflicts = FALSE)

setwd("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/results/")

jobids <- c("cis_eqtls_18373genes_age_sex_rin_batch_PC10", "cis_eqtls_18373genes_age_sex_rin_batch_PC10_PEER20", "cis_eqtls_18373genes_age_sex_rin_batch_PC10_NeutPCT_LympPCT_MonoPCT_EoPCT_BasoPCT")

cat("\nThis script merges postprocessed output from LIMIX into a single file.\nEdit the jobids object in the script to change the directories searched for results.\nMerging results from:\n")  

for(i in 1:length(jobids)) {
  resdir <- paste0("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/results/", jobids[i],"/processed")
  cat(paste0(resdir,"\n"))
  setwd(resdir)
  files <- list.files()
  files <- files[grepl("processed_qtl_results", files)]
  
  # Read in limix resul   ts
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
  
  if( length(datSplit) > 1) {
    # Save list of dropped features
    cat("Filtering variants with alpha < 0.1 or > 10\n")
    droppedFeatures <- datSplit[[2]] %>% 
      pull(feature_id) %>% 
      unique
    cat(paste0("Dropping: ", paste0(droppedFeatures, collapse = ", "), " due to alpha filter\n"))
    write.table(datSplit[[2]], row.names = F, col.names = T, quote = F, file = paste0(jobids[i], "_droppedfeatures.txt"))
  } else {
    cat("No variants with alpha < 0.1 or > 10\n")
  }
  
  # Remove dropped features and reformat results
  dat2 <- datSplit[[1]] %>%
    select(-fail) %>%
    data.frame
  
  # Restore RSID column
  rsid <- str_split_fixed(dat2$snp_id, ":", 2)
  dat2$rsid <- rsid[,2]
  dat2 <- dat2 %>%
    mutate(rsid = ifelse(rsid == "", snp_id, rsid))
  
  # Get column for reference allele
  dat2$alleles <- str_extract(dat2$snp_id, "[A-Z]_[A-Z]")
  dat2 <- dat2 %>%
    mutate(alleles = str_extract(snp_id, "[A-Z]+_[A-Z]+"))
  almat <-  str_split_fixed(dat2$alleles, "_", 2)  
  dat2$ref_allele <- ifelse(almat[,1] == dat2$assessed_allele, almat[,2], almat[,1])
  
  # Rearrange for saving
  dat3 <- dat2 %>% 
    arrange(feature_id) %>%
    select(feature_id, feature_chromosome:beta_param, snp_id, rsid, snp_chromosome:assessed_allele, ref_allele, call_rate:hwe_p, beta, beta_se, p_value, empirical_feature_p_value)
  
  # Write out
  outname <- paste0(jobids[i], "_results_merged.txt")
  cat(paste0("Writing out results to ",outname,"\n"))
  fwrite(dat3, quote = F, file = outname)

  rm(files, dat, dat2, rsid, datSplit, dat3, resdir)
}