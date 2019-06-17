# Script to generate chunking file to run limix on 100-200 genes at a time
setwd("/home/jm2294/projects/RNAseq/test_run_chunks")
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

# Assign each feature to a bin of about 100
anno2 <- data.frame()
for(i in 1:22){
  a <- anno %>% filter(chromosome == i) 
  numbins <- (length(a$start)/100) %>% floor
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
  arrange(chr, start) %>%
  select(-bin) %>%
  data.frame()

write.table(chunkList, file = "chunklist_b38.txt", sep = "\t", row.names = F, col.names =F)