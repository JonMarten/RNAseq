# get list of duplicated SNPs to filter out

library(dplyr)
a <- read.table("test_bgen_snp_stats.txt", h=T, stringsAsFactors = F)
b <- read.table("test_bgen_snp_stats_unique.txt", h=T, stringsAsFactors = F)

a2 <- a %>% 
  group_by(rsid) %>% 
  filter(n() > 1) %>% 
  data.frame

b2 <- b %>% 
  group_by(rsid) %>% 
  filter(n() > 1) %>% 
  data.frame

dupes <- data.frame(unique(a2$alternate_ids))
write.table(dupes, row.names = F, col.names = F, quote = F, file = "filter_snps.txt")