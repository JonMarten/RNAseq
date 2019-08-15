# Convert biomart output into annotation list for LIMIX use
setwd("/home/jm2294/rds/rds-jmmh2-projects/interval_rna_seq/analysis/00_testing/annotation_file")

library(data.table)
library(dplyr)

a <- fread("zcat mart_export.txt.gz", data.table=F)
names(a) <- gsub(" ", "_", names(a))
b <- a %>%
  select(feature_id = Gene_stable_ID,
         chromosome = "Chromosome/scaffold_name",
         start = "Gene_start_(bp)",
         end = "Gene_end_(bp)",
         feature_strand = Strand,
         gene_name = HGNC_symbol) %>%
  distinct() %>%
  filter(chromosome %in% c(1:22))  %>% 
  mutate(chromosome = as.numeric(chromosome)) %>%
  mutate(gene_name = ifelse(gene_name == "", feature_id, gene_name)) %>% # use ensembl gene id for genes with no HGNC symbol
  filter(!duplicated(feature_id)) %>% # remove 3 annotations where one gene has multiple transcripts. Retain one from each pair. 
  arrange(chromosome, start)

b$chromosome <- as.character(b$chromosome)
vec <- which(nchar(b$chromosome)==1)
b$chromosome[vec] <- paste0("0", b$chromosome[vec])

fwrite(b, sep = "\t", "Feature_Annotation_Ensembl_gene_ids_autosomes_b38.txt", quote = F)