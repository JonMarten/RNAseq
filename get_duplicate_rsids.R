# get list of duplicated SNPs to filter out
setwd("/home/jm2294/GENETIC_DATA/INTERVAL/RNAseq")
library(dplyr)
library(data.table)
a <- fread("c22_filtered_snp_stats.txt", skip = 8, data.table=F) # skip header of file
#a <- read.table("c22_test_bgen_snp_stats.txt", h=T, stringsAsFactors = F)

a2 <- a %>% 
  group_by(rsid) %>% 
  filter(n() > 1) %>% 
  data.frame

dupes <- data.frame(unique(a2$alternate_ids))
write.table(dupes, row.names = F, col.names = F, quote = F, file = "c22_filter_snps.txt")

#check for duplicates after making new file
b <- read.table("test_bgen_snp_stats_unique.txt", h=T, stringsAsFactors = F)
b2 <- b %>% 
  group_by(rsid) %>% 
  filter(n() > 1) %>% 
  data.frame