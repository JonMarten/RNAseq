# make bed format phenotypes for tensorQTL for COVID-19

#Chr start end ID UNR1 UNR2 UNR3 UNR4 
#chr1 173863 173864 ENSG123 -0.50 0.82 -0.71 0.83
#chr1 685395 685396 ENSG456 -1.13 1.18 -0.03 0.11
#chr1 700304 700305 ENSG789 -1.18 1.32 -0.36 1.26

library(data.table)
library(dplyr)

setwd("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/03_tensorqtl/phenotypes")

idmap <- fread(data.table = F, "/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/phenotype/sample_mapping_file_gt_to_phe_phase1.txt")

#phe <- fread("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/phenotype/INTERVAL_RNAseq_phase1_filteredSamplesGenes_TMMNormalised_FPKM_#ounts_foranalysis.txt", data.table = F)
phe <- fread("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/covid19/UNfilteredSamplesGenes_TMMNormalised_FPKM_Counts.csv", data.table = F)
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

# Remove IDs not in covariate file
covariates <- fread("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/03_tensorqtl/covariates/INTERVAL_RNAseq_phase1_age_sex_rin_batch_readDepth_PC10_PEER20.txt", data.table = F)
covarids <- covariates[1,-1] %>% as.character()
bedids <- names(bed[5:ncol(bed)])

keepids <- intersect(covarids, bedids)

bed2 <- bed %>%
  select("#Chr", start, end, ID, keepids)
covkeep <- c(1,which(covariates[1,] %in% keepids))
cov2 <- covariates[,covkeep]

#oldbed <- fread("zcat INTERVAL_RNAseq_phase1_filteredSamplesGenes_TMMNormalised_FPKM_Counts_foranalysis.bed.gz", data.table = F)

cov19 <- fread("/rds/user/jm2294/rds-jmmh2-projects/interval_rna_seq/scripts/RNAseq/covid-19/covid_genes_b37.csv", data.table = F)
bedcov19 <- bed2 %>% filter(ID %in% cov19$ensembl_id)


fwrite(bedcov19, sep = "\t", file = "/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/covid19/INTERVAL_RNAseq_phase1_filteredSamplesGenes_TMMNormalised_FPKM_Counts_foranalysis_COVID19.bed")
write.table(quote = F, row.names = F, col.names = F, cov2, sep = "\t", file = "/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/covid19/INTERVAL_RNAseq_COVID19_covariates.txt")

###########################################
## NOTE: The bed file must be compressed and indexed with the commands below:
# module load ceuadmin/tabix/0.2.6
# bgzip INTERVAL_RNAseq_phase1_filteredSamplesGenes_TMMNormalised_FPKM_Counts_foranalysis_COVID19.bed && tabix -p bed INTERVAL_RNAseq_phase1_filteredSamplesGenes_TMMNormalised_FPKM_Counts_foranalysis_COVID19.bed.gz
############################################


# GxE

gxe <- fread("../covariates/INTERVAL_RNAseq_phase1_GxE_neutPCT.txt", data.table = F)
gxe2 <- gxe[,covkeep]
write.table(quote = F, row.names = F, col.names = F, gxe2, sep = "\t", file = "/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/covid19/INTERVAL_RNAseq_COVID19_neutPCT_GxE.txt")


