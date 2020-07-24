# Convert biomart output into annotation list for LIMIX use
setwd("/home/jm2294/projects/RNAseq/test_run")
system("wget http://grch37.ensembl.org/biomart/martresults/179?file=martquery_0411142355_324.txt.gz") # Link will likely time out, see biomart_query_for_annotation.pl for query details

library(data.table)
library(dplyr)

a <- fread("zcat 179?file=martquery_0411142355_324.txt.gz", data.table=F)
names(a) <- gsub(" ", "_", names(a))
b <- a %>%
  select(feature_id = Gene_stable_ID,
         chromosome = "Chromosome/scaffold_name",
         start = "Gene_start_(bp)",
         end = "Gene_end_(bp)",
         ensembl_gene_id = Gene_stable_ID,
         feature_strand = Strand,
         gene_name = HGNC_symbol) %>%
  filter(chromosome %in% 1:22) %>% 
  filter(!duplicated(ensembl_gene_id))

fwrite(b, sep = "\t", "Feature_Annotation_Ensembl_gene_ids_autosomes.txt")