# make bed format phenotypes for tensorQTL

#Chr start end ID UNR1 UNR2 UNR3 UNR4 
#chr1 173863 173864 ENSG123 -0.50 0.82 -0.71 0.83
#chr1 685395 685396 ENSG456 -1.13 1.18 -0.03 0.11
#chr1 700304 700305 ENSG789 -1.18 1.32 -0.36 1.26

###########################################
## NOTE: Once created, the bed files must be compressed and indexed with 3_2_index_bed.sh
############################################


library(data.table)
library(dplyr)

setwd("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/04_phase2_full_analysis/phenotypes")

# Create mapping file to match phenotype to genotype
omictable <- fread("../covariates/processed/INTERVAL_omics_table_14MAY2020.csv", data.table = F)
idmap <- omictable %>%
  select(genotype_individual_id = affymetrix_ID, phenotype_individual_id = RNA_ID) %>%
  filter(!is.na(phenotype_individual_id))

# Read in TMM normalised FPKM feature counts
phe <- fread("/home/jm2294/rds/rds-jmmh2-projects/interval_rna_seq/analysis/04_phase2_full_analysis/peer_factors/peer_InputFiles/GeneExpr_PEER_TmmInvRankNormalised_swapsSwapped_mismatchRemoved.csv", data.table = F)
pheT <- phe %>% select(-V1) %>% t %>% data.frame
names(pheT) <- phe$V1
pheT$feature_id <- rownames(pheT)

anno <- fread("/home/jm2294/rds/rds-jmmh2-projects/interval_rna_seq/analysis/04_phase2_full_analysis/annotation_file/Feature_Annotation_Ensembl_gene_ids_autosomesPlusChrX_b38.txt", data.table = F)

bed <- left_join(pheT, anno[,1:4]) %>%
  select(Chr = chromosome, start, end, ID = feature_id, everything()) %>%
  arrange(Chr, start) %>%
  filter(!is.na(Chr)) %>%
  rename("#Chr" = Chr)

# Rename IDs to match genotype file
namevec <- base::match(names(bed)[5:ncol(bed)], idmap$phenotype_individual_id)
names(bed)[5:ncol(bed)] <-  as.character(idmap$genotype_individual_id[namevec])

# sort pheno file IDs
sortedids <- names(bed)[-c(1:4)] %>% sort
bed2 <- bed %>%
  select("#Chr", start, end, ID, sortedids)

fwrite(bed, sep = "\t", file = "INTERVAL_RNAseq_phase1-2_filteredSamplesGenes_TMMNormalised_FPKM_Counts_foranalysis.bed")


for(i in 1:22) {
  bedChr <- bed %>%
    rename(Chr = "#Chr") %>%
    filter(Chr == i)%>%
    rename("#Chr" = Chr)
  
  fwrite(bedChr, sep = "\t", file = paste0("INTERVAL_RNAseq_phase1-2_filteredSamplesGenes_TMMNormalised_FPKM_Counts_foranalysis_chr",i,".bed"))
  rm(bedChr)
}

###########################################

# Make Covariate file
#id UNR1 UNR2 UNR3 UNR4
#PC1 -0.02 0.14 0.16 -0.02
#PC2 0.01 0.11 0.10 0.01
#PC3 0.03 0.05 0.08 0.07
#BIN 1 0 0 1

# Function to reformat covariate files for TensorQTL
makeCov <- function(file, idmap, bed) {
  cov <- fread(file, data.table = F)
  cov <- cov %>%
    mutate(sample_id =  as.character(idmap$genotype_individual_id[base::match(sample_id, idmap$phenotype_individual_id)]))
  tcov <- transpose(cov)
  names(tcov) <- tcov[1,] 
  tcov <- tcov[-1,]
  tcov$id <- names(cov)[-1]
  # remove IDs missing in phenotype file
  remIDs <- names(tcov)[which(!names(tcov) %in% names(bed))][-1]
  tcov <- select(tcov, -remIDs)
  tcov <- select(tcov, id, sortedids) 
}


covMinimal <- makeCov("../../01_cis_eqtl_mapping/covariates/INTERVAL_RNAseq_phase1_age_sex_rin_batch_readDepth_PC10.txt", idmap, bed)
covPEER <- makeCov("../../01_cis_eqtl_mapping/covariates/INTERVAL_RNAseq_phase1_age_sex_rin_batch_readDepth_PC10_PEER20.txt", idmap, bed)
covCellPct <- makeCov("../../01_cis_eqtl_mapping/covariates/INTERVAL_RNAseq_phase1_age_sex_rin_batch_readDepth_PC10_NeutPCT_LympPCT_MonoPCT_EoPCT_BasoPCT.txt", idmap, bed)

fwrite(covMinimal, file = "../covariates/INTERVAL_RNAseq_phase1_age_sex_rin_batch_readDepth_PC10.txt", sep = "\t")
fwrite(covPEER, file = "../covariates/INTERVAL_RNAseq_phase1_age_sex_rin_batch_readDepth_PC10_PEER20.txt", sep = "\t")
fwrite(covCellPct, file = "../covariates/INTERVAL_RNAseq_phase1_age_sex_rin_batch_readDepth_PC10_NeutPCT_LympPCT_MonoPCT_EoPCT_BasoPCT.txt", sep = "\t")

# GxE file
# currently replcaed NAs with the median, just to see if these are screwing up the GxE
intr <- fread("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/covariates/INTERVAL_RNAseq_phase1_GxE_neutPCT.txt", data.table = F)
intr <- intr %>%
  mutate(sample_id =  as.character(idmap$genotype_individual_id[base::match(sample_id, idmap$phenotype_individual_id)]),
         NEUT_PCT___RNA = ifelse(is.na(NEUT_PCT___RNA), median(NEUT_PCT___RNA, na.rm=T), NEUT_PCT___RNA))
tintr <- transpose(intr)
names(tintr) <- tintr[1,] 
tintr <- tintr[-1,]
tintr$id <- names(intr)[-1]
tintr <- select(tintr, id, sortedids) 
# remove IDs missing in phenotype file
remIDs <- names(tintr)[which(!names(tintr) %in% names(bed))][-1]
tintr <- select(tintr, -remIDs)
fwrite(tintr, file = "/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/03_tensorqtl/covariates/INTERVAL_RNAseq_phase1_GxE_neutPCT.txt", sep = "\t")

