library(data.table)
library(dplyr)
library(Hmisc)

setwd("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/covariates")

covs <- fread("INTERVAL_RNA_batch1-8_covariates_release_2019_08_15.csv", data.table = F)
peer <- fread("PEER_factors_100Factors_plusCovariates_1000Iterations_AllBatches_TransformGenes_18373Genes.csv", data.table = F)

names(peer)[1] <- "sample_id"

namemap <- fread("../phenotype/sample_mapping_file_gt_to_phe_phase1.txt", data.table = F)
names(namemap)[2] <- "sample_id"

covs <- inner_join(namemap, covs)

allcovs <- inner_join(covs, select(peer, -age_RNA), by = "sample_id")

allcovs.m <- allcovs %>%
  select(-(genotype_individual_id:Asset_Name), -(Instrument:Well), -(Extr_Sample_ID:Box_Position), -(Extraction_Date:Extraction_Time), -(FluidX_Plate_Position:Freezer_Shelf), -Agilent_RINe, -(Normalization_Plate_ID:Normalization_Plate_Position), -(Rack_ID:Plate_Position), -intervalPhase, -ethnicity, -ABORH, -(attendanceDate___RNA:processTime___RNA)) %>%
  as.matrix() %>%
  na.exclude()

allcovs.m.cor <- allcovs.m %>% rcorr

allcovs[sapply(allcovs, is.character)] <- lapply(allcovs[sapply(allcovs, is.character)], as.factor)

allcovs.df <- allcovs %>%
  select(-(genotype_individual_id:Asset_ID))
  mutate(

lm(PEER1 ~ ., allcovs) %>% summary()




