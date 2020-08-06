library(dplyr)
library(data.table)
#library(ggplot2)
#library(cowplot)
#library(readr)
#theme_set(theme_cowplot())

setwd("/home/jm2294/rds/rds-jmmh2-projects/interval_rna_seq/analysis/04_phase2_full_analysis/results/cis_eQTLs")

for(i in 1:23) {
  cat(paste0("\nReading in Chr", i))
  cis_nom <- fread(data.table = F, paste0("tensorqtl_cis_MAF0.005_cisNominal_chr", i, ".csv"))
  cis <- fread(data.table = F, paste0("tensorqtl_cis_MAF0.005_cisPerGene_chr", i, ".csv"))
  map <- fread(data.table = F, paste0("../../genotypes/INTERVAL_RNAseq_Phase1-2_imputed_b38_biallelic_MAF0.005_chr", i, ".bim"))
  names(map) <- c("chr","variant_id","morg","pos_b38", "effect_allele", "other_allele")
  map <- select(map, -morg)
  
  # remove qval column that may have been added by tensorQTL script (depending if it was run with the version where that line hadn't been commented out)
  if(length(which(names(cis) == "qval")) > 0){
    cis <- select(cis, -qval)
  }
  # Apply BH correction to lowest locally-corrected p-value per gene
  cis <- cis %>%
    mutate(pval_BH = p.adjust(pval_beta, method = "BH"))
  eGenes <- cis %>% filter(pval_BH < 0.05)
  # Add SNP info back into cis results
  eGenes <- left_join(cis, map, by = "variant_id")
  
  # loop over all significant eGenes and get significant eSNPs
  eGenes$num_sig_eSNPs <- NA
  eSNPs <- data.frame()
  cat(paste0("\n\tGetting eSNPs for ", nrow(eGenes), " significant eGenes"))
  for(j in 1:nrow(eGenes)){
    cat(".")
    feature <- eGenes$Phenotype[j]
    cisSNPs <- cis_nom %>% 
      filter(phenotype_id == feature) %>%
      filter(!is.na(pval_nominal)) %>%
      filter(pval_nominal <= eGenes$pval_nominal_threshold[j])
    eGenes$num_sig_eSNPs[j] <- nrow(cisSNPs)
    eSNPs <- rbind(eSNPs, cisSNPs)
    rm(feature, cisSNPs)
  }
  
  fwrite(eGenes, file = paste0("tensorqtl_cis_MAF0.005_cis_chr", i, "_significant_eGenes.csv"))
  fwrite(eSNPs, file = paste0("tensorqtl_cis_MAF0.005_cis_chr", i, "_significant_eSNPs.csv"))
  rm(cis_nom, cis, map, eGenes, eSNPs)
}
              
              
