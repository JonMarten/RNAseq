# Script to generate chunking file to run limix on 100-200 genes at a time
setwd("/home/jm2294/projects/RNAseq/test_run")
Packages <- c("dplyr", "data.table")
lapply(Packages, library, character.only = TRUE)

anno <- fread("Feature_Annotation_Ensembl_gene_ids_autosomes.txt", data.table = F)
anno <- anno %>% 
  arrange(chromosome, start)

anno %>% 
  group_by(chromosome) %>%
  summarise(chr_start = min(start),
            chr_end = max(end)) %>%
  data.frame