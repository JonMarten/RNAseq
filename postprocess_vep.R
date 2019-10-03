library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)

setwd("/rds/user/jm2294/rds-jmmh2-projects/interval_rna_seq/analysis/00_testing/test_parameters")

c22 <- fread("results_merged_chr22_window500000_Perm100_MAF0.01_eSNPs.txt", data.table=F)
vep <- fread("vep.txt", data.table = F)
names(vep)[1] <- "Uploaded_variation"         

vep %>% 
  group_by(Uploaded_variation) %>%
  summarise(count = n()) %>%
  filter(count > 1) %>%
  arrange(desc(count)) %>%
  data.frame() %>%
  head(30)
  
22_50621055_A_G:rs76733141

c22$alleles <- str_extract(c22$snp_id, "[A-Z]_[A-Z]")
c22 <- c22 %>%
  mutate(alleles = str_extract(snp_id, "[A-Z]+_[A-Z]+"))
almat <-  str_split_fixed(c22$alleles, "_", 2)  
c22$ref_allele <- ifelse(almat[,1] == c22$assessed_allele, almat[,2], almat[,1])

str_extract("22_50621055_A_G:rs76733141", "_*_*:")
