library(data.table)
library(dplyr)

ann <- fread("Feature_Annotation_Ensembl_gene_ids_autosomes.txt", data.table=F)
ann2 <- ann %>% 
  arrange(chromosome, start) %>%
  filter(chromosome == 22, start > 23500000, end < 24500000)

phe <- fread("phenotype_5281-fc-genecounts.txt", data.table=F)
phe2 <- phe %>%
  filter(ENSEMBL_ID %in% ann2$ensembl_gene_id)

fwrite(phe2, sep = "\t", file = "phenotype_5281-fc-genecounts_chr22_23500000-24500000.txt")