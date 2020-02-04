library(data.table)
library(dplyr)
library(Hmisc)
library(pheatmap)
library(reshape2)

setwd("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/covariates")

covs <- fread("INTERVAL_RNA_batch1-8_covariates_release_2019_08_15.csv", data.table = F)
peer <- fread("PEER_factors_100Factors_plusCovariates_1000Iterations_AllBatches_TransformGenes_18373Genes.csv", data.table = F)

names(peer)[1] <- "sample_id"

namemap <- fread("../phenotype/sample_mapping_file_gt_to_phe_phase1.txt", data.table = F)
names(namemap)[2] <- "sample_id"

covs <- inner_join(namemap, covs)

allcovs <- inner_join(covs, select(peer, -age_RNA), by = "sample_id")
factors <- c("sequencingBatch", "sex","centre","batch1","batch2","batch3","batch4","batch5","batch6","batch7","batch8","sex2")
allcovs <- allcovs %>%
  select(-factors)

allcovs.m <- allcovs %>%
  select(-(genotype_individual_id:Asset_Name), -(Instrument:Well), -(Extr_Sample_ID:Box_Position), -(Extraction_Date:Extraction_Time), -(FluidX_Plate_Position:Freezer_Shelf), -Agilent_RINe, -(Normalization_Plate_ID:Normalization_Plate_Position), -(Rack_ID:Plate_Position), -intervalPhase, -ethnicity, -ABORH, -(attendanceDate___RNA:processTime___RNA),-intercept) %>%
  as.matrix() %>%
  na.exclude()


allcovs.m.cor <- allcovs.m %>% rcorr(type = "spearman")

  r <- data.frame(allcovs.m.cor$r) 
  row.names(r) <- names(r)
  r$pheno = names(r)
  rtall <- melt(r, id.vars = "pheno")
  rtall <- rtall %>%
    filter(grepl("PEER", pheno) & !grepl("PEER", variable) & !is.na(value) & !is.infinite(value)) %>%
    filter(abs(value) > 0.2)
  bloodNames <- fread("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/covariates/blood_trait_names.csv", data.table = F)
  names(bloodNames)[1] <- "variable"
  rtall2 <- left_join(rtall, bloodNames, by = "variable")
  rtall2 <- rtall2 %>% 
    mutate(variable = ifelse(!is.na(astleName), astleName, ifelse(!is.na(Description), Description, variable))) %>%
    filter(is.na(Use) | Use == 1)
  
  r2 <- dcast(rtall2, pheno ~ variable)
  
  #r2 <- r2 %>% 
  #  filter(pheno %in% paste0("PEER",1:20))
  row.names(r2) <- r2$pheno
  r2 <- select(r2, -pheno)
  r2[which(is.na(r2), arr.ind = T)] <- 0
  pheatmap(t(r2), display_numbers = T)


r <- r%>% 
  filter(!is.na(Asset_ID)) %>%
  select(-Sample_Volume_ul, -pheno)
pheatmap(r)



allcovs[sapply(allcovs, is.character)] <- lapply(allcovs[sapply(allcovs, is.character)], as.factor)

allcovs.df <- allcovs %>%
  select(-(genotype_individual_id:Asset_ID))
  mutate(

lm(PEER1 ~ ., allcovs) %>% summary()




