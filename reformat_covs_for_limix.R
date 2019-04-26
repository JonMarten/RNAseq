# Reformat covariates for LIMIX pipeline
library(dplyr)
library(data.table)
setwd("U:/Projects/RNAseq")

dat <- fread("covariates/INTERVALdata_02APR2019.csv", data.table=F)
ids <- fread("covariates/omicsMap.csv", data.table=F)
ids3 <- fread("covariates/omicsMap_P3.csv", data.table=F)
ids[which(ids == "", arr.ind=T)] <- NA
allids <- full_join(ids, ids3)
dat2 <- full_join(allids, dat)

# Filter to only individuals with RNA at either 24 and/or 48 months
ids_RNA <- ids %>% 
  filter(!is.na(RNAseq_QC_24m) | !is.na(RNAseq_QC_48m)) %>%
  mutate(id_RNA = ifelse(is.na(RNAseq_QC_24m), RNAseq_RAW_24m, RNAseq_QC_48m)) %>%
  select(identifier, Affymetrix_QC_bl, id_RNA)
 
dat_RNA <- left_join(ids_RNA, dat)

#map <- fread("v01_188id_preview/sample_mapping_file_gt_to_phe.txt", data.table=F)

rna_seq_ids <- fread("globus/RNA_seq_ids.csv", data.table = F)

a <- apply(dat2, MARGIN = 2, FUN = function(x){names(dat2)[which(map$phenotype_individual_id %in% x)]}) 

RNAcols <- grep("RNA", names(allids))
for(i in RNAcols){
  which(rna_seq_ids$RNA_id %in% allids[,i]) %>% length %>% print
}
