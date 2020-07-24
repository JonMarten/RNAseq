setwd("/home/jm2294/rds/rds-jmmh2-projects/interval_rna_seq/analysis/00_testing/annotation_file")

library(data.table)
library(dplyr)
library(ggplot2)

a <- fread("19_9_5_mart_export.txt", data.table = F)
names(a) <- gsub(" ", "_", names(a))
anno <- a %>%
  select(feature_id = Gene_stable_ID,
         chromosome = "Chromosome/scaffold_name",
         start = "Gene_start_(bp)",
         end = "Gene_end_(bp)",
         feature_strand = Strand,
         gene_name = HGNC_symbol,
         gene_type = Gene_type) %>%
  distinct() %>%
  filter(chromosome %in% c(1:22))  %>% 
  mutate(chromosome = as.numeric(chromosome)) %>%
  mutate(gene_name = ifelse(gene_name == "", feature_id, gene_name)) %>% # use ensembl gene id for genes with no HGNC symbol
  filter(!duplicated(feature_id)) %>% # remove 3 annotations where one gene has multiple transcripts. Retain one from each pair. 
  arrange(chromosome, start)

anno <- anno %>% 
  arrange(chromosome, start) %>%
  filter(gene_type == "protein_coding")

anno %>% 
  group_by(chromosome) %>%
  summarise(start_min = min(start),
            start_max = max(start),
            nFeatures = n()) %>%
  data.frame

# Assign each feature to a bin of about 50
anno2 <- data.frame()
for(i in 1:22){
  a <- anno %>% filter(chromosome == i) 
  numbins <- (length(a$start)/50) %>% floor
  b <- cut_number(a$start, n = numbins, labels = F)
  a$bin <- paste0("chr",i,"_",b)
  anno2 <- rbind(anno2,a)
  rm(a, numbins, b)
}

# Get start end end positions for chunks
chunkList <- anno2 %>%
  group_by(bin) %>%
  summarise(chr = unique(chromosome),
            start = min(start),
            end = max(end)) %>%
  mutate(chr = as.numeric(chr)) %>%
  arrange(chr, start) %>%
  select(-bin) %>%
  data.frame()

write.table(chunkList, file = "chunklist_b38_50genes_proteincodingonly.txt", sep = "\t", row.names = F, col.names =F, quote = F)

# Get chunk ranges for chromosome
chrTable <- chunkList %>% 
  mutate(chunkNo = 1:nrow(chunkList)) %>%
  group_by(chr) %>%
  summarise(startChunk = min(chunkNo), endChunk = max(chunkNo)) %>% 
  data.frame %>% 
  mutate(chr = as.numeric(chr)) %>%
  arrange(chr)
write.table(chrTable, row.names = F, col.names = T, file = "chunks_by_chr_50genes_proteincodingonly.txt", quote = F)