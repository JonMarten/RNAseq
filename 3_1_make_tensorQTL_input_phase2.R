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

# Output per-chromosome phenotype files for trans analysis
for(i in 1:23) {
  bedChr <- bed %>%
    rename(Chr = "#Chr") %>%
    filter(Chr == i)%>%
    rename("#Chr" = Chr)
  fwrite(bedChr, sep = "\t", file = paste0("INTERVAL_RNAseq_phase1-2_filteredSamplesGenes_TMMNormalised_FPKM_Counts_foranalysis_chr",i,".bed"))
  rm(bedChr)
}