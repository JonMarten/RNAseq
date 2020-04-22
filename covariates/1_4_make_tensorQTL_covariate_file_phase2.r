# Make Covariate file for TensorQTL
#id UNR1 UNR2 UNR3 UNR4
#PC1 -0.02 0.14 0.16 -0.02
#PC2 0.01 0.11 0.10 0.01
#PC3 0.03 0.05 0.08 0.07
#BIN 1 0 0 1
library(data.table)
library(dplyr)
setwd("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/04_phase2_full_analysis/covariates")

cov <- fread(data.table = F, "processed/INTERVAL_RNA_batch1-12_master_covariates_release_2020_04_02.csv")
phe <- fread(data.table = F, "zcat ../phenotypes/INTERVAL_RNAseq_phase1-2_filteredSamplesGenes_TMMNormalised_FPKM_Counts_initialcalling_foranalysis_chr22.bed.gz")

PCs <- fread("../genotypes/PCs/rnaseq_phase1.typed.maf1pc.geno3pc.hwe1minus5.complex-region-excl.indep-50-5-1.5.pruning.genome.pca.eigenvec", data.table = F) %>%
  select(-"#FID") %>%
  rename(affymetrix_ID = IID)
cov <- left_join(cov, PCs)

# NEED TO ADD IN READ DEPTH
cov2 <- cov %>%
  mutate(affymetrix_ID = as.character(affymetrix_ID)) %>%
  filter(affymetrix_ID %in% names(phe)) %>%
  select(id = affymetrix_ID, age_RNA, sex, Agilent_RINe, sequencingBatch, PC1:PC10) %>%
  mutate(sex = as.numeric(gsub(2, 0, sex)),
         RIN = as.numeric(ifelse(Agilent_RINe == "", NA, Agilent_RINe))) %>%
  filter(!is.na(RIN))

# convert batch to dummy variables
library(varhandle)
batch <- to.dummy(cov2$sequencingBatch, prefix = "batch") %>%
  data.frame()
names(batch) <- gsub("\\.", "", names(batch))
cov2 <- cbind(cov2, batch)

covOut <- cov2 %>%
  select(id:sex, RIN, paste0("batch", 1:12), PC1:PC10)

covOut.t <- covOut %>%
  t %>%
  data.frame(stringsAsFactors = F)


write.table(covOut.t, sep = "\t", quote = F, col.names = F, row.names = T, file = "INTERVAL_RNAseq_phase1-2_age_sex_rin_batch_PC10.txt")