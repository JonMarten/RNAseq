setwd("U:/Projects/RNAseq/annotation_file/")
library(data.table)
library(dplyr)

batch1 <- fread("counts.FC.unstranded.5281-tic109.txt", data.table=F)
features <- batch1$ENSEMBL_ID

a <- fread("mart_export.txt", data.table=F)
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
  mutate(gene_name = ifelse(gene_name == "", feature_id, gene_name)) # use ensembl gene id for genes with no HGNC symbol

c <- b %>% 
  mutate(geneLength = end - start)
c %>% filter(feature_id %in% c("ENSG00000269220", "ENSG00000183785", "ENSG00000093009", "ENSG00000215012", "ENSG00000273164", "ENSG00000093100"))

library(ggplot2)
ggplot(filter(c, chromosome == "22"), aes(x = geneLength)) + geom_density()


fwrite(b, sep = "\t", "Feature_Annotation_Ensembl_gene_ids_autosomes_b38.txt")