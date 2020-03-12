library(data.table)
library(dplyr)

setwd("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/03_tensorqtl/")

byChr <- list()
for (i in 1:22){
  byChr[[i]] <- fread(paste0("results/python_module_method/tensorqtl_trans_MAF0.005_chr",i, ".csv"), data.table = F)
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

fwrite(all, file = "results/python_module_method/tensorqtl_allSNPs_MAF0.005_merged_annotated.csv", sep = ",")

# check snp position is actually trans
all <- all %>%
  mutate(snpDistFromStart = ifelse(snp_chr == feat_chr, snp_bp - feat_start, NA), 
         snpDistFromEnd =  ifelse(snp_chr == feat_chr, snp_bp - feat_end, NA))

# Remove cis SNPs
cis <- all %>%
  filter(abs(snpDistFromStart) <= 500000 & abs(snpDistFromEnd) <= 500000)

trans <- all %>%
  filter(abs(snpDistFromStart) > 500000 & abs(snpDistFromEnd) > 500000)
fwrite(trans, file = "results/python_module_method/tensorqtl_transSNPs_MAF0.005_merged_annotated.csv", sep = ",")

# gene level FDR
topSNPs <- trans %>%
  group_by(phenotype_id) %>%
  summarise(nSNPs = n(), minP = min(pval, na.rm = T)) %>%
  data.frame() %>%
  mutate(minP_GWadj = minP * 10^6)

# Test dummmy p-value assumption - VALIDATED
##dumTest <- topSNPs %>%
##  select(-minP_GWadj) %>%
##  mutate(p.unif = minP,
##         p.unifBig = minP,
##         p.1 = minP) %>%
##  select(-minP)
##
##nphenos <- 17674 # There are this many phenotypes being tested
##nmiss <- nphenos - nrow(topSNPs)
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

sigThresh <- 5e-8 / length(eGenes.BH)

eSNPs.BH <- trans %>%
  filter(phenotype_id %in% eGenes.BH) %>%
  filter(pval < sigThresh)

eGenes.BH.summary <- eSNPs.BH %>%
  group_by(gene_name) %>%
  summarise(nSNPs = n(),
            nChr = length(unique(snp_chr))) %>% 
  data.frame

# Group SNPs by loci
sentinels <- list()
for(i in seq_along(eGenes.BH)){
  loop <- eSNPs.BH %>% filter(phenotype_id == eGenes.BH[i])
  

loop <- eSNPs.BH()
while(nrow(loop) > 0){
  



trans <- trans %>%
  mutate(sig_bonf = ifelse(pval < 5e-8/nrow(trans), 1, 0))

eGenes_bonf <- trans %>%
  filter(sig_bonf == 1) %>%
  group_by(gene_name) %>%
  summarise(n = n()) %>%
  data.frame()
              
              
