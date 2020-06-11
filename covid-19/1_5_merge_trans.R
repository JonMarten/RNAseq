# Merge trans results
setwd("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/covid19/results")

library(dplyr)
library(ggplot2)
library(cowplot)
library(data.table)
library(readr)
theme_set(theme_cowplot())

m <- data.frame()
for( i in c(1:22,"X")){
  fin <- fread(data.table = F, file = paste0("tensorqtl_trans_MAF0.005_all_age_sex_rin_batch_readDepth_PC10_PEER20_COVID19_CHR",i,".csv"))
  m <- rbind(m,fin)
  rm(fin)
}

automap <- fread(data.table = F, "/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/04_phase2_full_analysis/genotypes/INTERVAL_RNAseq_Phase1-2_imputed_b38_biallelic_MAF0.005_AllAutosomes.bim")
xmap <- fread("../genotypes/INTERVAL_chrX_merged_cleaned_RNAseq_phase1-2_deduplicated_MAF0.005.bim", data.table = F)

allmap <- rbind(automap, xmap)
names(allmap) <- c("CHR","variant_id","morg","BP", "A1", "A2")
 
m2 <- right_join(allmap, m) %>%
m2 <- m2 %>%
  select(-V1, -morg)

fwrite(m2, file = "tensorqtl_trans_MAF0.005_all_age_sex_rin_batch_readDepth_PC10_PEER20_COVID19_MERGED.csv")

library(qqman)
15,579,156-15,620,271 
ace2SNPs <- m2 %>% filter(CHR == "23", BP > 15579156 - 1000000, BP < 15620271 + 1000000) %>% pull(variant_id)
manhattan(m2, p = "pval", snp = "variant_id", annotatePval = 1e-7, highlight = ace2SNPs, ylim = c(4,8), annotateTop = F)
manhattan(filter(m2, maf > 0.01), p = "pval", snp = "variant_id", annotatePval = 1e-7, highlight = ace2SNPs, ylim = c(4,9), annotateTop = F)

tophits <- m2 %>% filter(pval < 1e-7)
write.csv(tophits, file = "tensorqtl_trans_MAF0.005_all_age_sex_rin_batch_readDepth_PC10_PEER20_COVID19_MERGED_tophits1e-7.csv")
