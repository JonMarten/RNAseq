# Script to generate chunking file to run limix on 100-200 genes at a time
setwd("/rds/user/jm2294/rds-jmmh2-projects/interval_rna_seq/analysis/00_testing/annotation_file")
Packages <- c("dplyr", "data.table","ggplot2")
lapply(Packages, library, character.only = TRUE)

anno <- fread("../annotation_file/Feature_Annotation_Ensembl_gene_ids_autosomes_b38.txt", data.table = F)
anno <- anno %>% 
  arrange(chromosome, start)

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

write.table(chunkList, file = "chunklist_b38_50genes.txt", sep = "\t", row.names = F, col.names =F, quote = F)

# Get chunk ranges for chromosome
chrTable <- chunkList %>% 
  mutate(chunkNo = 1:nrow(chunkList)) %>%
  group_by(chr) %>%
  summarise(startChunk = min(chunkNo), endChunk = max(chunkNo)) %>% 
  data.frame %>% 
  mutate(chr = as.numeric(chr)) %>%
  arrange(chr)
write.table(chrTable, row.names = F, col.names = T, file = "/home/jm2294/rds/hpc-work/projects/RNAseq/scripts/RNAseq/reformat_bgen/chunks_by_chr_50genes.txt", quote = F)