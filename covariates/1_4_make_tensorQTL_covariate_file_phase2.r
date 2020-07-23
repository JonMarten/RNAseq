# Make Covariate file for TensorQTL
#id UNR1 UNR2 UNR3 UNR4
#PC1 -0.02 0.14 0.16 -0.02
#PC2 0.01 0.11 0.10 0.01
#PC3 0.03 0.05 0.08 0.07
#BIN 1 0 0 1
library(data.table)
library(dplyr)
library(fastDummies)
setwd("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/04_phase2_full_analysis/covariates")

cov <- fread(data.table = F, "processed/INTERVAL_RNA_batch1-12_master_covariates_release_2020_07_01.csv")
phe <- fread(data.table = F, "zcat ../phenotypes/INTERVAL_RNAseq_phase1-2_filteredSamplesGenes_TMMNormalised_FPKM_Counts_foranalysis_chr22.bed.gz")

PCs <- fread("../genotypes/PCs/rnaseq_phase1.typed.maf1pc.geno3pc.hwe1minus5.complex-region-excl.indep-50-5-1.5.pruning.genome.pca.eigenvec", data.table = F) %>%
  select(-"#FID") %>%
  rename(affymetrix_ID = IID)
cov <- left_join(cov, PCs)

# Needs updating to use latest files in peer_30Fact_100Iter_sampleSwapsFixed when this has finished running, this is for prototyping purposes 
peer <- fread("../peer_factors/peer_30Fact_100Iter_sampleSwapsFixed/PEER_factors.txt", data.table = F)
peer <- peer %>%
  rename(sample_id = V1) %>%
  select(sample_id, PEER1:PEER20)

cov <- left_join(cov, peer)

# make dummy cols
batchCols <- cov %>% 
  select(sequencingBatch) %>%
  mutate(sequencingBatch = as.factor(sequencingBatch)) %>%
  dummy_cols

# Create a variable for season that is numerically sequential and equally spaced - an interger corresponding to sample colleciton month and year. This is to stop 201512 having a huge value jump to 201601
datelevels <- paste0(c(rep(2015, 12), rep(2016, 12)), sprintf("%02d", 1:12))
cov <- cov %>%
  mutate(RNA_sample_year = substr(attendanceDate___RNA, 1,4)) %>%
  mutate(yearmonth = paste0(RNA_sample_year, sprintf("%02d", RNA_sample_month))) %>%
  mutate(season = as.numeric(factor(yearmonth, levels = datelevels)))

cov <- cbind(cov, batchCols[,-1])

# Select analysis covariates. Set missing RIN values to the median
cov2 <- cov %>%
  mutate(affymetrix_ID = as.character(affymetrix_ID)) %>%
  filter(affymetrix_ID %in% names(phe)) %>%
  select(id = affymetrix_ID, 
         age_RNA, 
         sex,
         BMI, 
         season,
         Agilent_RINe,
         sequencingBatch_1:sequencingBatch_12,
         Conc_ng_ul,
         BASO_PCT___RNA, 
         EO_PCT___RNA, 
         HCT_PCT___RNA, 
         HGB_g_dL___RNA, 
         IRF_PCT___RNA, 
         LYMPH_PCT___RNA, 
         MCH_pg___RNA, 
         MCHC_g_dL___RNA, 
         MCV_fL___RNA, 
         MONO_PCT___RNA, 
         MPV_fL___RNA, 
         NEUT_PCT___RNA, 
         PCT_PCT___RNA, 
         PDW_fL___RNA, 
         PLT_F_10_9_L___RNA, 
         RBC_10_12_L___RNA, 
         RDW_SD_fL___RNA, 
         RET_PCT___RNA, 
         WBC_10_9_L___RNA,
         PC1:PC10,
         PEER1:PEER20
         ) %>%
  mutate(sex = as.numeric(gsub(2, 0, sex)),
         RIN = as.numeric(ifelse(Agilent_RINe == "", NA, Agilent_RINe))) %>%
  mutate(RIN = ifelse(is.na(RIN), median(RIN, na.rm = T), RIN)) %>%
  select(-Agilent_RINe)

# Remove NAs. Oh yuck, I'm going to do it with a loop, but it's quick
naToMedian <- function(x) {
  x <- ifelse(is.na(x), median(x, na.rm = T), x)
  return(x)
}
naCols <- which(apply(cov2, 2, anyNA))
for(i in naCols){
  cov2[, i] <- naToMedian(cov2[,i])
}

# Reorder cov to match phe 
pheids <- names(phe)[-(1:4)]
ord <- match(pheids, cov2$id) 
cov2 <- cov2[ord,]

covOut.t <- cov2 %>%
  t %>%
  data.frame(stringsAsFactors = F)

write.table(covOut.t, sep = "\t", quote = F, col.names = F, row.names = T, file = "INTERVAL_RNAseq_phase1-2_fullcovariates_foranalysis.txt")