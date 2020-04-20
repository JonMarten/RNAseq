library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)

setwd("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/covid19")

a <- fread("Covid-19_Genes_TMMNormalised_FPKM_Counts.csv", data.table = F)

genes <- a[,1:6]

reads <- a %>%
  select(-c(feature_id, chromosome:gene_length)) %>%
  gather(-gene_name, key = "individual", value = "count")

ggplot(reads, aes(x = gene_name, y = count, colour = gene_name)) +
  geom_boxplot() +
  geom_jitter(width = 0.5, size = 0.1)

ggplot(reads, aes(x = count, fill = gene_name)) +
  geom_density()

raw <- fread("../analysis/04_phase2_full_analysis/phenotypes/raw/interval_basic-star-fc-genecounts.txt", data.table = F)