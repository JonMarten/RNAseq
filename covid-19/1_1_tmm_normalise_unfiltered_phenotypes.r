#------------------------------------------------------------------------

# Filtering the count data
 
#------------------------------------------------------------------------

# -------------------------------
# Load libraries
#--------------------------------
library(edgeR)
library(limma)
library(data.table)
library(dplyr)


# ----------------------
# set working directory
# ----------------------
setwd("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/04_phase2_full_analysis/phenotypes/")

# -----------------------------------------------------------------------------
# Read in the data from all the batches, annotation file, metadata file
# -----------------------------------------------------------------------------

# read the meta-data
metadata <- read.csv("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/04_phase2_full_analysis/covariates/processed/INTERVAL_RNA_batch1-12_master_covariates_release_2020_04_02.csv")

# Read in the feature counts
dat <- fread("raw/interval_basic-star-fc-genecounts.txt", data.table = F)

# read the annotation file
#annotation <- read.csv("~/Dropbox/Artika/Projects/INTERVAL/Analysis/Batches1_8/GeneAnnotationFile_GTF/GeneAnnotationFile_EnsembltoGeneIds.csv")
annotation <- fread("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/annotation_file/Feature_Annotation_Ensembl_gene_ids_autosomes_b38.txt", data.table = F)
annotation <- annotation %>%
  mutate(gene_length = end - start)
ACE2 <- data.frame("feature_id" = "ENSG00000130234",
                   "chromosome" = 23,
                   "start" = 15579156,
                   "end" = 15620271,
                   "feature_strand" = -1,
                   "gene_name" = "ACE2",
                   "gene_length" = 41115)
annotation <- rbind(annotation, ACE2)

# Get the gene symbols corresponding to ENSEMBLE IDs and gene length   
counts <- merge(annotation, dat, by.x = "feature_id", by.y = "ENSEMBL_ID")

# -----------------------------------------------------------------------------
# Remove samples with 
# -----------------------------------------------------------------------------
# Remove samples in batch2 with "DO NOT PROCESS" and "EMPTY" information in the metadata                                
# These are N=9 samples (note the empty individual is not in the count data)
remove_batch2 <- c("INT_RNA7709000", "INT_RNA7709004", "INT_RNA7709017", "INT_RNA7709022", "INT_RNA7709151", "INT_RNA7709209", "INT_RNA7709213", "INT_RNA7709279", "INT_RNA7709286")

#Remove samples with RIN < 4
# There are N=15 individuals
metadata$RIN <- as.numeric(as.character(metadata$Agilent_RINe))
remove_RIN_Less4 <- as.character(metadata$sample_id[which(metadata$RIN < 4)])


#Remove samples with ReadDepth Less than 10
# There are N=5 individuals
ReadDepth <- colSums(counts[,12:ncol(counts)])
remove_ReadDepth_lessthan_10Million <- names(ReadDepth)[which(ReadDepth < 10000000)]

# Total remove N=28 individuals
remove_ind <- c(remove_batch2, remove_RIN_Less4)

# -----------------------------------------------------------------------------
# Create DGEList object y and remove outlier genes, using the FC data
# -----------------------------------------------------------------------------
counts1 <- counts[, -which(colnames(counts) %in% remove_ind)]

y <- edgeR::DGEList(counts=counts1[,8:ncol(counts1)], genes=counts1[,1:7])

# Filter lowly expressed genes
# -------------
# Select genes with > 0.5 CPM in at least 10% of the samples
keep <- rowSums(cpm(y) > 0.5) >= round(ncol(y$counts)*0.1) #10% is N=276 individuals

summary(keep)
#y1 <- y[keep, , keep.lib.sizes = FALSE]
y1 <- y
# TMM Normalisation
y2 <- calcNormFactors(y1)


# get RPKM values 
rpkms <- rpkm(y2, normalized.lib.sizes=TRUE, log=TRUE, prior.count=0.25, gene.length = y$genes$Length)

genes <- y2$genes[, c("feature_id", "gene_name", "chromosome", "start", "end","gene_length")]
rpkms1 <- cbind(genes, rpkms)
#fwrite(rpkms1, "/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/covid19/phenotypes/UNfilteredSamplesGenes_TMMNormalised_FPKM_Counts_Phase1-2_initialcalling.csv", sep = ",")
fwrite(rpkms1, "/home/jm2294/covid/phenotypes/UNfilteredSamplesGenes_TMMNormalised_FPKM_Counts_Phase1-2_initialcalling.csv", sep = ",")