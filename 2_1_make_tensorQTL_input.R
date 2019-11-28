# make bed format phenotypes for tensorQTL

#Chr start end ID UNR1 UNR2 UNR3 UNR4 
#chr1 173863 173864 ENSG123 -0.50 0.82 -0.71 0.83
#chr1 685395 685396 ENSG456 -1.13 1.18 -0.03 0.11
#chr1 700304 700305 ENSG789 -1.18 1.32 -0.36 1.26

library(data.table)
library(dplyr)

setwd("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/03_tensorqtl/phenotypes")

phe <- fread("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/phenotype/INTERVAL_RNAseq_phase1_filteredSamplesGenes_TMMNormalised_FPKM_Counts_foranalysis.txt", data.table = F)
phe <- phe %>%
  rename(feature_id = gene_id)

anno <- fread("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/annotation_file/Feature_Annotation_Ensembl_gene_ids_autosomes_b38.txt", data.table = F)

bed <- left_join(phe, anno[,1:4]) %>%
  select(Chr = chromosome, start, end, ID = feature_id, INT_RNA7427205:INT_RNA7959012) %>%
  arrange(Chr, start) %>%
  rename("#Chr" = Chr)

fwrite(bed, sep = "\t", file = "INTERVAL_RNAseq_phase1_filteredSamplesGenes_TMMNormalised_FPKM_Counts_foranalysis.bed")
  
# Make Covariate file
#id UNR1 UNR2 UNR3 UNR4
#PC1 -0.02 0.14 0.16 -0.02
#PC2 0.01 0.11 0.10 0.01
#PC3 0.03 0.05 0.08 0.07
#BIN 1 0 0 1

cov <- fread("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/covariates/INTERVAL_RNAseq_phase1_age_sex_rin_batch_PC10_PEER20.txt", data.table = F)
tcov <- transpose(cov)
names(tcov) <- tcov[1,] 
tcov <- tcov[-1,]
tcov$id <- names(cov)[-1]
tcov <- select(tcov, id, INT_RNA7711053:INT_RNA7878943)

fwrite(tcov, file = "/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/03_tensorqtl/covariates/INTERVAL_RNAseq_phase1_age_sex_rin_batch_PC10_PEER20.txt", sep = "\t")
