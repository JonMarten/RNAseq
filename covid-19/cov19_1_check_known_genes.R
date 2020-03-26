# 

setwd("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/covid19")
library(data.table)
library(dplyr)

phe <- fread(data.table = F, "UNfilteredSamplesGenes_TMMNormalised_FPKM_Counts.csv")

cov19 <- fread("/rds/user/jm2294/rds-jmmh2-projects/interval_rna_seq/scripts/RNAseq/covid-19/covid_genes_b37.csv", data.table = F)

phecov19 <- phe %>% filter(feature_id %in% cov19$ensembl_id)

covgenes1 <- phe %>% filter(feature_id %in% cov19$ensembl_id)
idmap <- fread("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/phenotype/sample_mapping_file_gt_to_phe_phase1.txt", data.table = F)

rnaids <- names(covgenes1)[-c(1:6)]
gtids <- idmap$genotype_individual_id[match(rnaids, idmap$phenotype_individual_id)] %>% as.character()

names(covgenes1)[-c(1:6)] <- gtids 
write.csv(col.names = T, row.names = F, covgenes1, file = "Covid-19_Genes_TMMNormalised_FPKM_Counts.csv", quote = F)

# eQTLs
all <- fread("/rds/user/jm2294/rds-jmmh2-projects/interval_rna_seq/analysis/03_tensorqtl/results/python_module_method/tensorqtl_trans_MAF0.005_age_sex_rin_batch_readDepth_PC10_PEER20_merged_annotated.csv", data.table = F)

covAll <- all %>% 
  filter(phenotype_id %in% cov19$ensembl_id) %>%
  mutate(sig = ifelse(pval < (5e-8 / 2), 1, 0))

fwrite(covAll, sep = ",", file = "Covid-19_genes_eQTLs.csv")

library(qqman)

png(filename="covid-19_eqtl_manhattan.png", 
    type="cairo",
    units="in", 
    width=30, 
    height=10, 
    pointsize=18, 
    res=120)
  par(mfrow = c(2,1))
  manhattan(filter(covAll, gene_name == "CTSB") , chr = "snp_chr", bp = "snp_bp", p = "pval", main = "CTSB", ylim = c(0,20)) 
  manhattan(filter(covAll, gene_name == "CTSL") , chr = "snp_chr", bp = "snp_bp", p = "pval", main = "CTSL", ylim = c(0,20)) 
dev.off()
png(filename="covid-19_eqtl_manhattan_CTSB_chr8.png", 
    type="cairo",
    units="in", 
    width=30, 
    height=10, 
    pointsize=18, 
    res=120)
manhattan(filter(covAll, gene_name == "CTSB" & snp_chr == 8) , chr = "snp_chr", bp = "snp_bp", p = "pval", main = "CTSB", ylim = c(0,120))
dev.off()

