# Convert biomart output into annotation list for LIMIX use
setwd("/home/jm2294/projects/RNAseq/annotation_file/")

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
  filter(chromosome %in% c(1:22, "X", "Y", "MT"))  %>% 
  mutate(gene_name = ifelse(gene_name == "", feature_id, gene_name)) %>% # use ensembl gene id for genes with no HGNC symbol
  filter(!duplicated(feature_id)) # remove 3 annotations where one gene has multiple transcripts. Retain one from each pair.
fwrite(b, sep = "\t", "Feature_Annotation_Ensembl_gene_ids_autosomes_b38.txt")