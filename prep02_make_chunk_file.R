# Make chunk file for features that pass Artika's expression level filter only. 5 features per chunk.

setwd("/rds/user/jm2294/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/annotation_file")

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

# Read in filtered protein matrix from Artika.
# Note that some proteins are in the gene counts but not in ensembl as they have been retired.
system("awk -F \"\\\"*,\\\"*\" \'{ print $1,$2,$3,$4,$5 }\' ../phenotype/gene_expr_filtered_TMM_rpkms_genesInvRankTransf_2747Samp_18373genes.csv  > ../phenotype/protlist.txt")

prots <- fread("../phenotype/protlist.txt", data.table = F)
anno <- anno %>% 
  arrange(chromosome, start) %>%
  filter(feature_id %in% prots$gene_id)

anno %>% 
  group_by(chromosome) %>%
  summarise(start_min = min(start),
            start_max = max(start),
            nFeatures = n()) %>%
  data.frame

# Assign each feature to a bin of about 5
anno2 <- data.frame()
for(i in 1:22){
  a <- anno %>% filter(chromosome == i) 
  numbins <- (length(a$start)/5) %>% floor
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

write.table(chunkList, file = "chunklist_b38_5genes_filtered.txt", sep = "\t", row.names = F, col.names =F, quote = F)

# Get chunk ranges for chromosome
chrTable <- chunkList %>% 
  mutate(chunkNo = 1:nrow(chunkList)) %>%
  group_by(chr) %>%
  summarise(startChunk = min(chunkNo), endChunk = max(chunkNo)) %>% 
  data.frame %>% 
  mutate(chr = as.numeric(chr)) %>%
  arrange(chr)
write.table(chrTable, row.names = F, col.names = T, file = "chunks_by_chr_5genes_filtered.txt", quote = F)