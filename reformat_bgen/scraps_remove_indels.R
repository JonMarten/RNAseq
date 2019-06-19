library(data.table)
library(dplyr)
a <- fread("interval_b38_snpstats_filtered_by_rsid.txt", data.table = F, skip = 8)

b <- a %>% mutate(chrpos = paste0(chromosome, "_", position))

c <- b %>%
  group_by(chrpos) %>%
  filter(n() > 1) %>%
  data.frame()

snpsOnly <- b %>% group_by(chrpos) %>% filter(n()==1) %>% data.frame

out <- snpsOnly %>% select("snp_id" = alternate_ids) 

write.table(out, sep = " ", row.names = F, col.names = T, quote = F, file = "b38_biallelic_snps_only_no_indels.txt")

outIndels <- c  %>% select("snp_id" = alternate_ids) 
write.table(outIndels, sep = " ", row.names = F, col.names = T, quote = F, file = "b38_indels_only.txt")
