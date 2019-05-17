# Script to generate chunking file to run limix on 100-200 genes at a time
setwd("/home/jm2294/projects/RNAseq/test_run")
Packages <- c("dplyr", "data.table","ggplot2")
lapply(Packages, library, character.only = TRUE)

anno <- fread("Feature_Annotation_Ensembl_gene_ids_autosomes.txt", data.table = F)
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



breaks <- a[c(seq(1,length(a), 100), length(a))]
b <- cut(a, breaks, labels = F, include.lowest = T)
cut_number(1:length(a), numbins) %>% levels



