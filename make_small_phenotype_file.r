library(data.table)
library(dplyr)

ann <- fread("Feature_Annotation_Ensembl_gene_ids_autosomes.txt", data.table=F)
ann2 <- ann %>% 
  arrange(chromosome, start) %>%
  filter(chromosome == 22, start > 23500000, end < 24500000)

phe <- fread("phenotype_5281-fc-genecounts.txt", data.table=T)