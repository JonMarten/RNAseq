# Convert biomart output into annotation list for making TensorQTL input
setwd("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/04_phase2_full_analysis/annotation_file")

library(data.table)
library(dplyr)

a <- fread("19_9_5_mart_export.txt", data.table = F)
names(a) <- gsub(" ", "_", names(a))
b <- a %>%
  select(feature_id = Gene_stable_ID,
         chromosome = "Chromosome/scaffold_name",
         start = "Gene_start_(bp)",
         end = "Gene_end_(bp)",
         feature_strand = Strand,
         gene_name = HGNC_symbol) %>%
  distinct() %>%
  mutate(chromosome = gsub("X","23", chromosome)) %>%
  filter(chromosome %in% c(1:23))  %>% 
  mutate(chromosome = as.numeric(chromosome)) %>%
  mutate(gene_name = ifelse(gene_name == "", feature_id, gene_name)) %>% # use ensembl gene id for genes with no HGNC symbol
  filter(!duplicated(feature_id)) %>% # remove 3 annotations where one gene has multiple transcripts. Retain one from each pair. 
  arrange(chromosome, start)

fwrite(b, sep = "\t", "Feature_Annotation_Ensembl_gene_ids_autosomesPlusChrX_b38.txt", quote = F)