setwd("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/03_tensorqtl/results")
library(data.table)
library(dplyr)

trans <- fread("tensorqtl_trans_test.trans_qtl_pairs.csv", data.table = F)

bim <- fread("../genotypes/INTERVAL_b38_autosomes_RNAseqPhase1_biallelic_all_MAF0.05.bim", data.table = F)
names(bim) <- c("snp_chr","variant_id", "snp_morgan","snp_bp","snp_A1","snp_A2")

trans2 <- left_join(trans, bim, by = "variant_id")

genes <- fread("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/annotation_file/Feature_Annotation_Ensembl_gene_ids_autosomes_b38.txt", data.table = F)
names(genes) <- c("phenotype_id", "feat_chr","feat_start","feat_end","feat_strand","gene_name")
trans3 <- left_join(trans2, genes, by = "phenotype_id")

fwrite(trans3, file = "tensorqtl_trans_test.trans_qtl_pairs_annotated.csv", sep = ",")