#------------------------------------------------------------------------

# Filtering the count data 

#------------------------------------------------------------------------

# -------------------------------
# Load libraries
#--------------------------------
library(edgeR)
library(limma)
library(data.table)


# ----------------------
# set working directory
# ----------------------
setwd("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/covid19")

# -----------------------------------------------------------------------------
# Read in the data from all the batches, annotation file, metadata file
# -----------------------------------------------------------------------------

# read the meta-data
metadata <- read.csv("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/covariates/INTERVAL_RNA_batch1-8_covariates_release_2019_08_15.csv")

path <- "/rds/project/jmmh2/rds-jmmh2-pre_qc_data/interval/rna_seq/raw_data/globus_phase1/"

# read the count data from each batch and store in a list 
batch1 <- fread(paste0(path,  "batch1/matrices/counts.FC.unstranded.5281-tic109.txt"),  data.table = F)
batch2 <- fread(paste0(path, "batch2/results-study5591-tic109b/combined/study5591-tic109b-star-fc-genecounts.txt"), data.table = F)
colnames(batch2) <- gsub(".star.gene.fc.txt", "", colnames(batch2))
batch3 <- fread(paste0(path, "batch3/results-study5591-tic109d/combined/study5591-tic109d-star-fc-genecounts.txt"), data.table = F)
batch4 <- fread(paste0(path, "batch4/results-study5591-tic276/combined/study5591-tic276-star-fc-genecounts.txt"), data.table = F)
batch5 <- fread(paste0(path, "batch5/results-study5591-tic328/combined/study5591-tic328-star-fc-genecounts.txt"), data.table = F)
batch6 <- fread(paste0(path, "batch6/results-study5591-tic297/combined/study5591-tic297-star-fc-genecounts.txt"), data.table = F)
batch7 <- fread(paste0(path, "batch7/results-study5591-tic329/combined/study5591-tic329-star-fc-genecounts.txt"), data.table = F)
batch8 <- fread(paste0(path, "batch8/results-study5591-tic364/combined/study5591-tic364-star-fc-genecounts.txt"), data.table = F)

# read the annotation file
#annotation <- read.csv("~/Dropbox/Artika/Projects/INTERVAL/Analysis/Batches1_8/GeneAnnotationFile_GTF/GeneAnnotationFile_EnsembltoGeneIds.csv")
annotation <- fread("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/annotation_file/Feature_Annotation_Ensembl_gene_ids_autosomes_b38.txt", data.table = F)
annotation <- annotation %>%
  mutate(gene_length = end - start)

#----------------------------------------------------------------
# Combine the batches together 
#----------------------------------------------------------------
exp <- list(batch1, batch2[, -1], batch3[, -1], batch4[,-1], batch5[,-1], batch6[, -1], batch7[, -1], batch8[,-1])
#names(exp) <- c("batch1", "batch2", "batch3", "batch4", "batch5", "batch6", "batch7", "batch8")
dat <- do.call(cbind, exp)
fwrite(dat, "counts_batches_1_8_raw.csv", sep = ",")

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
remove_ind <- c(remove_batch2, remove_RIN_Less4, remove_ReadDepth_lessthan_10Million)

# -----------------------------------------------------------------------------
# Create DGEList object y and remove outlier genes, using the FC data
# -----------------------------------------------------------------------------
counts1 <- counts[, -which(colnames(counts) %in% remove_ind)]

y <- edgeR::DGEList(counts=counts1[,8:ncol(counts1)], genes=counts1[,1:7])

# Filter lowly expressed genes
# -------------
# Select genes with > 0.5 CPM in at least 10% of the samples
keep <- rowSums(cpm(y) > 0.5) >= round(ncol(y$counts)*0.1) #10% is N=276 individuals
#keep3 <- rowSums(cpm(y) > 1) >= round(ncol(y$counts)*0.1)

summary(keep)
#y1 <- y[keep, , keep.lib.sizes = FALSE]
y1 <- y
# TMM Normalisation
y2 <- calcNormFactors(y1)


# get RPKM values 
rpkms <- rpkm(y2, normalized.lib.sizes=TRUE, log=TRUE, prior.count=0.25, gene.length = y$genes$Length)

genes <- y2$genes[, c("feature_id", "gene_name", "chromosome", "start", "end","gene_length")]
rpkms1 <- cbind(genes, rpkms)
fwrite(rpkms1, "UNfilteredSamplesGenes_TMMNormalised_FPKM_Counts.csv", sep = ",")





