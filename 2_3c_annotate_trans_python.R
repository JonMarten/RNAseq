library(data.table)
library(dplyr)

setwd("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/03_tensorqtl/")

byChr <- list()
for (i in 1:22){
  byChr[[i]] <- fread(paste0("results/python_module_method/tensorqtl_trans_MAF0.005_chr",i, ".csv"), data.table = F)
}

trans <- bind_rows(byChr)
trans <- trans %>%
  select(-V1)

# annotate with SNP and gene information
bim <- fread("genotypes/INTERVAL_b38_autosomes_RNAseqPhase1_biallelic_all_MAF0.005.bim", data.table = F)
names(bim) <- c("snp_chr","variant_id", "snp_morgan","snp_bp","snp_A1","snp_A2")
trans <- left_join(trans, bim, by = "variant_id")

genes <- fread("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/annotation_file/Feature_Annotation_Ensembl_gene_ids_autosomes_b38.txt", data.table = F)
names(genes) <- c("phenotype_id", "feat_chr","feat_start","feat_end","feat_strand","gene_name")
trans <- left_join(trans, genes, by = "phenotype_id")

trans <- trans %>% select(-snp_morgan)

trans <- trans %>% arrange(feat_chr, feat_start, snp_chr, snp_bp)

fwrite(trans, file = "results/python_module_method/tensorqtl_trans_MAF0.005_merged_annotated.csv", sep = ",")

# check snp position is actually trans
trans <- trans %>%
  mutate(snpDistFromStart = ifelse(snp_chr == feat_chr, snp_bp - feat_start, NA), 
         snpDistFromEnd =  ifelse(snp_chr == feat_chr, snp_bp - feat_end, NA))



# gene level FDR
topSNPs <- trans %>%
  group_by(phenotype_id) %>%
  summarise(nSNPs = n(), minP = min(pval, na.rm = T)) %>%
  data.frame() %>%
  mutate(minP_GWadj = minP * 10^6) %>%
  mutate(minP_GWadj_BH = p.adjust(minP_GWadj, method = "BH")) %>%
  mutate(sig = ifelse(minP_GWadj_BH < 0.05, 1, 0))

trans <- trans %>%
  mutate(sig_bonf = ifelse(pval < 5e-8/nrow(trans), 1, 0))

eGenes_bonf <- trans %>%
  filter(sig_bonf == 1) %>%
  group_by(gene_name) %>%
  summarise(n = n()) %>%
  data.frame()
              
              