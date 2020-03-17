# make bed format phenotypes for tensorQTL

#Chr start end ID UNR1 UNR2 UNR3 UNR4 
#chr1 173863 173864 ENSG123 -0.50 0.82 -0.71 0.83
#chr1 685395 685396 ENSG456 -1.13 1.18 -0.03 0.11
#chr1 700304 700305 ENSG789 -1.18 1.32 -0.36 1.26

library(data.table)
library(dplyr)

setwd("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/03_tensorqtl/phenotypes")

idmap <- fread(data.table = F, "/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/phenotype/sample_mapping_file_gt_to_phe_phase1.txt")

phe <- fread("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/phenotype/INTERVAL_RNAseq_phase1_filteredSamplesGenes_TMMNormalised_FPKM_Counts_foranalysis.txt", data.table = F)
phe <- phe %>%
  rename(feature_id = gene_id)

anno <- fread("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/annotation_file/Feature_Annotation_Ensembl_gene_ids_autosomes_b38.txt", data.table = F)

bed <- left_join(phe, anno[,1:4]) %>%
  select(Chr = chromosome, start, end, ID = feature_id, INT_RNA7427205:INT_RNA7959012) %>%
  arrange(Chr, start) %>%
  filter(!is.na(Chr)) %>%
  rename("#Chr" = Chr)

# Rename IDs to match genotype file
namevec <- base::match(names(bed)[5:ncol(bed)], idmap$phenotype_individual_id)
names(bed)[5:ncol(bed)] <-  as.character(idmap$genotype_individual_id[namevec])

# sort pheno file IDs
sortedids <- names(bed)[-c(1:4)] %>% sort
bed <- bed %>%
  select("#Chr", start, end, ID, sortedids)

fwrite(bed, sep = "\t", file = "INTERVAL_RNAseq_phase1_filteredSamplesGenes_TMMNormalised_FPKM_Counts_foranalysis.bed")
###########################################
## NOTE: The bed file must be compressed and indexed with the commands below:
# module load ceuadmin/tabix/0.2.6
# bgzip INTERVAL_RNAseq_phase1_filteredSamplesGenes_TMMNormalised_FPKM_Counts_foranalysis.bed && tabix -p bed INTERVAL_RNAseq_phase1_filteredSamplesGenes_TMMNormalised_FPKM_Counts_foranalysis.bed.gz
############################################

# Make Covariate file
#id UNR1 UNR2 UNR3 UNR4
#PC1 -0.02 0.14 0.16 -0.02
#PC2 0.01 0.11 0.10 0.01
#PC3 0.03 0.05 0.08 0.07
#BIN 1 0 0 1

# Function to reformat covariate files for TensorQTL
makeCov <- function(file, idmap, bed) {
  cov <- fread(file, data.table = F)
  cov <- cov %>%
    mutate(sample_id =  as.character(idmap$genotype_individual_id[base::match(sample_id, idmap$phenotype_individual_id)]))
  tcov <- transpose(cov)
  names(tcov) <- tcov[1,] 
  tcov <- tcov[-1,]
  tcov$id <- names(cov)[-1]
  # remove IDs missing in phenotype file
  remIDs <- names(tcov)[which(!names(tcov) %in% names(bed))][-1]
  tcov <- select(tcov, -remIDs)
  tcov <- select(tcov, id, sortedids) 
}


covMinimal <- makeCov("../../01_cis_eqtl_mapping/covariates/INTERVAL_RNAseq_phase1_age_sex_rin_batch_readDepth_PC10.txt", idmap, bed)
covPEER <- makeCov("../../01_cis_eqtl_mapping/covariates/INTERVAL_RNAseq_phase1_age_sex_rin_batch_readDepth_PC10_PEER20.txt", idmap, bed)
covCellPct <- makeCov("../../01_cis_eqtl_mapping/covariates/INTERVAL_RNAseq_phase1_age_sex_rin_batch_readDepth_PC10_NeutPCT_LympPCT_MonoPCT_EoPCT_BasoPCT.txt", idmap, bed)

fwrite(covMinimal, file = "../covariates/INTERVAL_RNAseq_phase1_age_sex_rin_batch_readDepth_PC10.txt", sep = "\t")
fwrite(covPEER, file = "../covariates/INTERVAL_RNAseq_phase1_age_sex_rin_batch_readDepth_PC10_PEER20.txt", sep = "\t")
fwrite(covCellPct, file = "../covariates/INTERVAL_RNAseq_phase1_age_sex_rin_batch_readDepth_PC10_NeutPCT_LympPCT_MonoPCT_EoPCT_BasoPCT.txt", sep = "\t")

# GxE file
# currently replcaed NAs with the median, just to see if these are screwing up the GxE
intr <- fread("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/covariates/INTERVAL_RNAseq_phase1_GxE_neutPCT.txt", data.table = F)
intr <- intr %>%
  mutate(sample_id =  as.character(idmap$genotype_individual_id[base::match(sample_id, idmap$phenotype_individual_id)]),
         NEUT_PCT___RNA = ifelse(is.na(NEUT_PCT___RNA), median(NEUT_PCT___RNA, na.rm=T), NEUT_PCT___RNA))
tintr <- transpose(intr)
names(tintr) <- tintr[1,] 
tintr <- tintr[-1,]
tintr$id <- names(intr)[-1]
tintr <- select(tintr, id, sortedids) 
# remove IDs missing in phenotype file
remIDs <- names(tintr)[which(!names(tintr) %in% names(bed))][-1]
tintr <- select(tintr, -remIDs)
fwrite(tintr, file = "/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/03_tensorqtl/covariates/INTERVAL_RNAseq_phase1_GxE_neutPCT.txt", sep = "\t")

