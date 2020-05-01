# make bed format phenotypes for tensorQTL for COVID-19

#Chr start end ID UNR1 UNR2 UNR3 UNR4 
#chr1 173863 173864 ENSG123 -0.50 0.82 -0.71 0.83
#chr1 685395 685396 ENSG456 -1.13 1.18 -0.03 0.11
#chr1 700304 700305 ENSG789 -1.18 1.32 -0.36 1.26

library(data.table)
library(dplyr)

#setwd("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/covid19/phenotypes")
setwd("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/covid19/phenotypes")

# Create mapping file to match phenotype to genotype
omictable <- fread("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/04_phase2_full_analysis/covariates/processed/INTERVAL_omics_table_02APR2020.csv", data.table = F)
idmap <- omictable %>%
  select(genotype_individual_id = affymetrix_ID, phenotype_individual_id = RNA_ID) %>%
  filter(!is.na(phenotype_individual_id))

#phe <- fread("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/phenotype/INTERVAL_RNAseq_phase1_filteredSamplesGenes_TMMNormalised_FPKM_#ounts_foranalysis.txt", data.table = F)
phe <- fread("UNfilteredSamplesGenes_TMMNormalised_FPKM_Counts_Phase1-2_initialcalling.csv", data.table = F)

anno <- fread("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/04_phase2_full_analysis/annotation_file/Feature_Annotation_Ensembl_gene_ids_autosomesPlusChrX_b38.txt", data.table = F)

bed <- left_join(phe, anno[,1:4]) %>%
  select(Chr = chromosome, start, end, ID = feature_id, INT_RNA7427205:INT_RNA7960548) %>%
  arrange(Chr, start) %>%
  filter(!is.na(Chr)) %>%
  rename("#Chr" = Chr)

# Rename IDs to match genotype file
namevec <- base::match(names(bed)[5:ncol(bed)], idmap$phenotype_individual_id)
names(bed)[5:ncol(bed)] <-  as.character(idmap$genotype_individual_id[namevec])
missvec <- which(is.na(names(bed)))
bed <- bed[,-missvec]

# Remove IDs not in covariate file/genotype file
covariates <- fread("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/04_phase2_full_analysis/covariates/INTERVAL_RNAseq_phase1-2_age_sex_rin_batch_PC10.txt", data.table = F)
covarids <- covariates[1,-1] %>% as.character()
bedids <- names(bed[5:ncol(bed)])

genoids <- fread("../genotypes/INTERVAL_chrX_merged_cleaned_RNAseq_phase1-2.fam", data.table = F)
genoids <- genoids$V1 %>% as.character

keepids <- intersect(covarids, bedids)
keepids <- intersect(keepids, genoids)

bed2 <- bed %>%
  select("#Chr", start, end, ID, keepids)
covkeep <- c(1,which(covariates[1,] %in% keepids))
cov2 <- covariates[,covkeep]

cov19 <- fread("/rds/user/jm2294/rds-jmmh2-projects/interval_rna_seq/scripts/RNAseq/covid-19/covid_genes_b37.csv", data.table = F)
bedcov19 <- bed2 %>% filter(ID %in% cov19$ensembl_id[-5])


fwrite(bedcov19, sep = "\t", file = "/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/covid19/phenotypes/INTERVAL_RNAseq_phase1-2_UNfilteredSamplesGenes_TMMNormalised_FPKM_Counts_foranalysis_COVID19.bed")
write.table(quote = F, row.names = F, col.names = F, cov2, sep = "\t", file = "/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/covid19/covariates/INTERVAL_RNAseq_COVID19_covariates.txt")

###########################################
## NOTE: The bed file must be compressed and indexed with the commands below:
# module load ceuadmin/tabix/0.2.6
# bgzip INTERVAL_RNAseq_phase1-2_UNfilteredSamplesGenes_TMMNormalised_FPKM_Counts_foranalysis_COVID19.bed && tabix -p bed INTERVAL_RNAseq_phase1-2_UNfilteredSamplesGenes_TMMNormalised_FPKM_Counts_foranalysis_COVID19.bed.gz
############################################

# Remove all individuals with no gene counts for ACE2
ace2 <- bedcov19 %>%
  filter(ID == "ENSG00000130234")
ace2.t <- t(ace2[,-(1:4)])
zerovec <- which(ace2.t == min(ace2.t)) + 4
ace2.no0 <- ace2[,-zerovec]
ace2ids <- names(ace2.no0)[-(1:4)]
covkeep.ace2 <- c(1,which(covariates[1,] %in% ace2ids))
cov.ace2 <- covariates[,covkeep.ace2]

fwrite(ace2.no0 , sep = "\t", file = "/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/covid19/phenotypes/INTERVAL_RNAseq_phase1-2_UNfilteredSamplesGenes_TMMNormalised_FPKM_Counts_foranalysis_ACE2_no_zeros.bed")
write.table(cov.ace2, quote = F, row.names = F, col.names = F, sep = "\t", file = "/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/covid19/covariates/INTERVAL_RNAseq_COVID19_covariates_ACE2_no_zeros.txt")
# module load ceuadmin/tabix/0.2.6
# bgzip INTERVAL_RNAseq_phase1-2_UNfilteredSamplesGenes_TMMNormalised_FPKM_Counts_foranalysis_ACE2_no_zeros.bed && tabix -p bed INTERVAL_RNAseq_phase1-2_UNfilteredSamplesGenes_TMMNormalised_FPKM_Counts_foranalysis_ACE2_no_zeros.bed.gz

# Output list of individuals with non-zero ACE2 to filter plink genotypes
plinkout <- data.frame("FID" = ace2ids, "IID" = ace2ids)
write.table(plinkout, sep = "\t", quote = F, col.names = F, row.names = F, file = "/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/covid19/genotypes/ace2_nonzero_ids.txt")


# GxE
gxe <- fread("../covariates/INTERVAL_RNAseq_phase1_GxE_neutPCT.txt", data.table = F)
gxe2 <- gxe[,covkeep]
write.table(quote = F, row.names = F, col.names = F, gxe2, sep = "\t", file = "/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/covid19/INTERVAL_RNAseq_COVID19_neutPCT_GxE.txt")


