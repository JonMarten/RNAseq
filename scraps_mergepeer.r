library(data.table)
library(dplyr)

setwd("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/coviariates")

covs <- fread("INTERVAL_RNA_batch1-8_covariates_release_2019_08_15.csv", data.table = F)
peer <- fread("PEER_factors_100Factors_plusCovariates_1000Iterations_AllBatches_TransformGenes_18373Genes.csv", data.table = F)

names(peer)[1] <- "sample_id"

pcs <- fread("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/genetic_PCs/rnaseq_phase1.typed.maf1pc.geno3pc.hwe1minus5.complex-region-excl.indep-50-5-1.5.pruning.genome.pca.eigenvec", data.table = F)
namemap <- fread("../phenotype/sample_mapping_file_gt_to_phe_phase1.txt", data.table = F)
names(namemap)[1] <- "IID"
pcs2 <- inner_join(namemap, pcs)
names(pcs2)[2] <- "sample_id"
covs2 <- inner_join(pcs2, covs)

covs2 <- covs2 %>% 
  mutate(merge_id = paste0(age_RNA, "_", Agilent_RINe,"_", PC1))

peer <- peer %>%
  mutate(merge_id = paste0(age_RNA, "_", RIN,"_", PC1))

covs3 <- inner_join(covs2, peer)