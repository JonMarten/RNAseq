# Merge and filter postprocessed limix output

library(data.table)
library(dplyr)
library(stringr)

setwd("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/results/")

jobids <- c("cis_eqtls_18373genes_age_sex_rin_batch_PC10", "cis_eqtls_18373genes_age_sex_rin_batch_PC10_PEER20")

for( i in 1:length(jobids)) {

  resdir <- paste0("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/results/", jobids[i],"/processed")
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
    print("Filtering variants with alpha < 0.1 or > 10")
    droppedFeatures <- datSplit[[2]] %>% 
      pull(feature_id) %>% 
      unique
    print(paste0("Dropping: ", paste0(droppedFeatures, collapse = ", "), " due to alpha filter"))
    write.table(datSplit[[2]], row.names = F, col.names = T, quote = F, file = paste0("/rds/user/jm2294/rds-jmmh2-projects/interval_rna_seq/analysis/00_testing/test_parameters/results_merged_chr22_window", jobdf$Window[i],"_Perm", jobdf$Permutations[i],"_MAF",jobdf$MAF[i], "_droppedfeatures.txt"))
  } else {
    print("No variants with alpha < 0.1 or > 10")
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
  fwrite(dat3, quote = F, file = outname)

  rm(files, dat, dat2, rsid, datSplit, dat3, resdir)
}