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

# Read in ids from RNA seq data
rna_seq_ids <- fread("globus/RNA_seq_ids.csv", data.table = F)

# filter to just RNA id columns and rows with an RNA id
RNAcols <- grep("RNA", names(allids))
RNArows <- rowSums(is.na(rnaids[,-1])) != ncol(rnaids[,-1]) # list rows with at least one non-NA RNA id
rnaids <- allids %>% 
  select(identifier, RNAcols) %>%
  filter(RNArows) %>%
  mutate(RNAseq_gwasQC_24m = as.character(RNAseq_gwasQC_24m),
         RNAseq_gwasQC_48m = as.character(RNAseq_gwasQC_48m),
         RNAseq_gwasQC_p3 = as.character(RNAseq_gwasQC_p3))

# Make single column for RNA identifier
RNA_any <- NA
for(i in 1:nrow(rnaids)){
  ids <- as.character(rnaids[i,-1])
  RNA_any[i] <- paste(unique(na.exclude(ids)), collapse = ",")
}
rnaids$RNA_any <- RNA_any

# Check overlap with RNA seq data identifiers
rnamiss <- which(!rna_seq_ids$RNA_id %in% rnaids$RNA_any)
phemiss <-which(!rnaids$RNA_any %in% rna_seq_ids$RNA_id )

if(length(rnamiss) == 0){
  cat("\n No RNA seq ids missing in phenotype data.")
} else {
  cat("\n",
      length(rnamiss),
      " RNA seq ids missing in phenotype data:", 
      paste0(rna_seq_ids$RNA_id[rnamiss]))
}

if(length(phemiss) == 0){
  cat("\n No phenotype data ids missing in RNA seq.")
} else {
  cat("\n",
      length(phemiss),
      " phenotype data ids missing in RNA seq:", 
      paste0(rnaids$RNA_any[phemiss], collapse = ", "))
}

rna_id_mapper <- rnaids %>%
  select(identifier, RNA_any)
