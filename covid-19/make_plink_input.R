# make plink files
library(data.table)
library(dplyr)

#setwd("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/covid19/phenotypes")
setwd("/home/jm2294/covid")

# Create mapping file to match phenotype to genotype
omictable <- fread("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/04_phase2_full_analysis/covariates/processed/INTERVAL_omics_table_02APR2020.csv", data.table = F)
idmap <- omictable %>%
  select(genotype_individual_id = affymetrix_ID, phenotype_individual_id = RNA_ID) %>%
  filter(!is.na(phenotype_individual_id))

#phe <- fread("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/phenotype/INTERVAL_RNAseq_phase1_filteredSamplesGenes_TMMNormalised_FPKM_#ounts_foranalysis.txt", data.table = F)
phe <- fread("phenotypes/UNfilteredSamplesGenes_TMMNormalised_FPKM_Counts_Phase1-2_initialcalling.csv", data.table = F)

cov19 <- fread("/rds/user/jm2294/rds-jmmh2-projects/interval_rna_seq/scripts/RNAseq/covid-19/covid_genes_b37.csv", data.table = F)
cov19 <- cov19[-5,]

phe2 <- phe %>% filter(gene_name %in% cov19$gene_name)
phe.t <- phe2[,7:ncol(phe2)] %>% t %>% data.frame
names(phe.t) <- phe2$gene_name
namin <- function(x){
  y = ifelse(x == min(x), NA, x)
  return(y)
}
phe.t.nomin <- apply(X = phe.t, MARGIN = 2, FUN = namin)
phe.t.nomin <- data.frame(phe.t.nomin)
phe.t$ACE2.nomin <- phe.t.nomin$ACE2
phe.t$TMPRSS2.nomin <- phe.t.nomin$TMPRSS2

phe.t$IID <- ownames(phe.t)
phe.t$IID <- idmap$genotype_individual_id[match(phe.t$IID, idmap$phenotype_individual_id)]
phe.t$FID <- phe.t$IID
phe.plink <- phe.t %>%
  select(FID, IID, ACE2, ACE2.nomin, TMPRSS2, TMPRSS2.nomin, CTSL, CTSB)

write.table(phe.plink, row.names = F, col.names = T, quote = F, file = "/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/covid19/plink/covid_RNAseq_pheno.txt", sep = "\t")

# 
covariates <- fread("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/04_phase2_full_analysis/covariates/processed/INTERVAL_RNA_batch1-12_master_covariates_release_2020_04_02.csv", data.table = F)
PCs <- fread("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/04_phase2_full_analysis/genotypes/PCs/rnaseq_phase1.typed.maf1pc.geno3pc.hwe1minus5.complex-region-excl.indep-50-5-1.5.pruning.genome.pca.eigenvec", data.table = F) %>%
  select(-"#FID") %>%
  rename(affymetrix_ID = IID)
cov <- left_join(covariates, PCs)
cov$FID <- cov$affymetrix_ID
cov2 <- cov %>%
  select(FID, IID = affymetrix_ID, age_RNA, sex, Agilent_RINe, sequencingBatch, PC1:PC10)

library(varhandle)
batch <- to.dummy(cov2$sequencingBatch, prefix = "batch") %>%
  data.frame()
names(batch) <- gsub("\\.", "", names(batch))
cov2 <- cbind(cov2, batch)
cov2 <- cov2 %>%select(-sequencingBatch)
cov2[which(cov2 == "", arr.ind = T)] <- NA

write.table(cov2, row.names = F, col.names = T, quote = F, file = "/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/covid19/plink/covid_RNAseq_cov.txt", sep = "\t")
