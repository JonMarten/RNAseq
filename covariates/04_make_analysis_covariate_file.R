library(dplyr)
library(data.table)
setwd("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/covariates")

covs <- fread("INTERVAL_RNA_batch1-8_covariates_release_2019_08_15.csv", data.table = F)
peer <- fread("PEER_factors_100Factors_plusCovariates_1000Iterations_AllBatches_TransformGenes_18373Genes.csv", data.table = F)
namemap <- fread("../phenotype/sample_mapping_file_gt_to_phe_phase1.txt", data.table = F)

names(peer)[1] <- "sample_id"
names(namemap)[2] <- "sample_id"
covs <- inner_join(namemap, covs)

allcovs <- inner_join(covs, select(peer, -age_RNA), by = "sample_id")

allcovs$sex <- gsub(2, 0, allcovs$sex)

covout <- allcovs %>%
  select(sample_id, RIN, age_RNA, sex, ReadDepth, batch1:batch8, PC1:PC10) %>%
  filter(!is.na(RIN))

covout.peer <- allcovs %>%
  select(sample_id, RIN, age_RNA, sex, ReadDepth, batch1:batch8, PC1:PC10, PEER1:PEER20) %>%
  filter(!is.na(RIN))

covout.bloodcells <- allcovs %>%
  select(sample_id, RIN, age_RNA, sex, ReadDepth, batch1:batch8, PC1:PC10, NEUT_PCT___RNA, LYMPH_PCT___RNA, MONO_PCT___RNA, EO_PCT___RNA, BASO_PCT___RNA) %>%
  filter(!is.na(RIN) & !is.na(BASO_PCT___RNA))


fwrite(covout, file = "INTERVAL_RNAseq_phase1_age_sex_rin_batch_readDepth_PC10.txt", sep = "\t")
fwrite(covout.peer, file = "INTERVAL_RNAseq_phase1_age_sex_rin_batch_readDepth_PC10_PEER20.txt", sep = "\t")
fwrite(covout.bloodcells, file = "INTERVAL_RNAseq_phase1_age_sex_rin_batch_readDepth_PC10_NeutPCT_LympPCT_MonoPCT_EoPCT_BasoPCT.txt", sep = "\t")

# Make GxE interaction covariate file
gxeOut <- allcovs %>%
  select(sample_id, NEUT_PCT___RNA)
fwrite(gxeOut, file = "INTERVAL_RNAseq_phase1_GxE_neutPCT.txt", sep = "\t")


