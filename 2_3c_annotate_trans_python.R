# Script to aggregate trans results and call significant eGenes and eSNPs. Call with Rscript and trailing arguments for file prefix, with chr replaced with "XX" because I was an idiot when I named my files and this bodge is the best we can do for now.
#Rscript 2_3c_annotate_trans_python.R results/python_module_method/tensorqtl_trans_MAF0.005_chrXX_age_sex_rin_batch_readDepth_PC10
#Rscript 2_3c_annotate_trans_python.R results/python_module_method/tensorqtl_trans_MAF0.005_chrXX_age_sex_rin_batch_readDepth_PC10_PEER20
#Rscript 2_3c_annotate_trans_python.R results/python_module_method/tensorqtl_trans_MAF0.005_chrXX_age_sex_rin_batch_readDepth_PC10_NeutPCT_LympPCT_MonoPCT_EoPCT_BasoPCT
library(data.table)
library(dplyr)
mainpath <- "/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/03_tensorqtl/"
setwd(mainpath)

#args = c("results/python_module_method/tensorqtl_trans_MAF0.005_chrXX_age_sex_rin_batch_readDepth_PC10_NeutPCT_LympPCT_MonoPCT_EoPCT_BasoPCT")
args = commandArgs(trailingOnly=TRUE)
prefix = args[1]
outprefix <- gsub("_chrXX_", "_", prefix)

byChr <- list()
for (i in 1:22){
  file <- paste0(gsub("XX",i, prefix), ".csv")
  byChr[[i]] <- fread(file, data.table = F)
}

all <- bind_rows(byChr)
all <- all %>%
  select(-V1)

# annotate with SNP and gene information
bim <- fread("genotypes/INTERVAL_b38_autosomes_RNAseqPhase1_biallelic_all_MAF0.005.bim", data.table = F)
names(bim) <- c("snp_chr","variant_id", "snp_morgan","snp_bp","snp_A1","snp_A2")
all <- left_join(all, bim, by = "variant_id")

genes <- fread("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/annotation_file/Feature_Annotation_Ensembl_gene_ids_autosomes_b38.txt", data.table = F)
names(genes) <- c("phenotype_id", "feat_chr","feat_start","feat_end","feat_strand","gene_name")
all <- left_join(all, genes, by = "phenotype_id")

all <- all %>% select(-snp_morgan)

all <- all %>% arrange(feat_chr, feat_start, snp_chr, snp_bp)

fwrite(all, file = paste0(outprefix, "_merged_annotated.csv"), sep = ",")

# check snp position is actually trans
all <- all %>%
  mutate(snpDistFromStart = ifelse(snp_chr == feat_chr, snp_bp - feat_start, NA), 
         snpDistFromEnd =  ifelse(snp_chr == feat_chr, snp_bp - feat_end, NA))

# Remove cis SNPs
cis <- all %>%
  filter(abs(snpDistFromStart) <= 500000 & abs(snpDistFromEnd) <= 500000)

trans <- all %>%
  filter(abs(snpDistFromStart) > 500000 & abs(snpDistFromEnd) > 500000)
fwrite(trans, file = paste0(outprefix, "_merged_annotated_transSNPs.csv"), sep = ",")
#trans <- fread("results/python_module_method/tensorqtl_transSNPs_MAF0.005_merged_annotated.csv", data.table = F)

# gene level FDR
topSNPs <- trans %>%
  group_by(phenotype_id) %>%
  summarise(nSNPs = n(), minP = min(pval, na.rm = T)) %>%
  data.frame() %>%
  mutate(minP_GWadj = minP * 10^6)

nphenos <- 17674 # There are this many phenotypes being tested
nmiss <- nphenos - nrow(topSNPs)

# Test dummmy p-value assumption - VALIDATED
##dumTest <- topSNPs %>%
##  select(-minP_GWadj) %>%
##  mutate(p.unif = minP,
##         p.unifBig = minP,
##         p.1 = minP) %>%
##  select(-minP)
##
##dumP.unif <- runif(nmiss, min = 1e-5, max = 1)
##dumP.unifBig <- runif(nmiss, min = 0.5, max = 1)
##dumP.1 <- rep(1, nmiss)
##
##dummy <- data.frame(phenotype_id = paste0("dummy_",1:nmiss),
##                    nSNPs = 1,
##                    p.unif = dumP.unif,
##                    p.unifBig = dumP.unifBig,
##                    p.1 = dumP.1)
##
##dumTest <- rbind(dumTest, dummy) %>%
##  mutate(p.unif.BH = p.adjust(p.unif, method = "BH"),
##         p.unifBig.BH = p.adjust(p.unifBig, method = "BH"),
##         p.1.BH = p.adjust(p.1, method = "BH"))
##
##dumTest <- dumTest %>% 
##  mutate(match = ifelse(p.unif.BH == p.unifBig.BH & p.unif.BH == p.1.BH, 1, 0)) %>%
##  mutate(dummy = ifelse(grepl("dummy", phenotype_id), 1, 0))
##table("match" = dumTest$match, "dummy" = dumTest$dummy)

# Add in dummy SNPs for all genes that had p-vals > 1e-5 and so aren't in the results file
dummy <- data.frame(phenotype_id = paste0("dummy_",1:nmiss),
                    nSNPs = 1,
                    minP = 1, 
                    minP_GWadj = 1) 
topSNPs <- rbind(topSNPs,dummy) %>%
  mutate(minP_GWadj = ifelse(minP_GWadj > 1, 1, minP_GWadj)) %>%
  mutate(minP_GWadj_BH = p.adjust(minP_GWadj, method = "BH")) %>%
  mutate(sig = ifelse(minP_GWadj_BH < 0.05, 1, 0))
#remove dummies
topSNPs <- topSNPs %>%
  filter(!grepl("dummy", phenotype_id))
  
# get signifiant trans eGenes and eSNPs
eGenes.BH <- topSNPs %>%
  filter(sig == 1) %>%
  pull(phenotype_id)

#sigThresh <- 5e-8 / length(eGenes.BH)
# Use largest uncorrected p-value for a significant eGene as threshold 
sigThresh <- topSNPs %>%
  filter(sig == 1) %>%
  pull(minP) %>%
  max

eSNPs.BH <- trans %>%
  filter(phenotype_id %in% eGenes.BH) %>%
  filter(pval < sigThresh)

fwrite(eSNPs.BH, file = paste0(outprefix, "_merged_annotated_transSNPs_eSNPs.csv"), sep = ",")
     
  
eGenes.BH.summary <- eSNPs.BH %>%
  group_by(gene_name) %>%
  summarise(nSNPs = n(),
            nChr = length(unique(snp_chr))) %>% 
  data.frame

# Group SNPs by loci
sentinels <- data.frame()
for(i in seq_along(eGenes.BH)){
  geneSentinels <- data.frame()
  cat(paste0("\nFeature ",i, ": ", eGenes.BH[i]))
  loop <- eSNPs.BH %>% 
    filter(phenotype_id == eGenes.BH[i]) %>%
    arrange(pval)
  j = 1
  # select one top SNP per 1Mb window
  while(nrow(loop) > 0){
    cat(paste0("\n\tLocus ",j))
    row <- loop[1,]
    snpchr <- loop$snp_chr[1]
    maxpos <- loop$snp_bp[1] + 500000
    minpos <- loop$snp_bp[1] - 500000
    locus <- loop %>%
      filter((snp_chr == snpchr & snp_bp > minpos & snp_bp < maxpos))
    row$locus_size = nrow(locus)
    row$locus_start = minpos
    row$locus_end = maxpos
    geneSentinels <- rbind(geneSentinels, row)
    loop <- loop %>% filter(!variant_id %in% locus$variant_id)
    rm(snpchr, maxpos, minpos, row)
    j = j + 1
  }
  # Merge sentinels with overlapping windows
  if(nrow(geneSentinels) > 1){
    superSentinels <- data.frame()
    while(nrow(geneSentinels) > 1){
      row <- geneSentinels[1,]
      locus <- geneSentinels %>%
        filter((snp_chr == row$snp_chr[1] & snp_bp > row$snp_bp - 1000000 & snp_bp < row$snp_bp + 1000000))
      if(nrow(locus) > 1) {
        row$locus_size <- sum(locus$locus_size)
        row$locus_start <- min(locus$locus_start)
        row$locus_end <- max(locus$locus_end)
        superSentinels <- rbind(superSentinels, row)
      } else {
        superSentinels <- rbind(superSentinels, locus)
      }
      geneSentinels <- geneSentinels %>% filter(!variant_id %in% locus$variant_id)
      
      cat("\n\t\t Merging ", nrow(locus), " sentinels with overlapping windows: ", paste0(locus$variant_id, collapse = ", "))
      
      rm(row, locus)
    }
    sentinels <- rbind(sentinels, superSentinels)
  } 
  sentinels <- rbind(sentinels, geneSentinels)
  
}
fwrite(sentinels, file = paste0(outprefix, "_merged_annotated_transSNPs_eSNPs_sentinels.csv"), sep = ",")
