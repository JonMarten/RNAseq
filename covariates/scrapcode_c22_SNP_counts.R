library(data.table)
library(dplyr)

c22 <- fread("../snp_stats/impute_22_interval_snp_stats_unfiltered.txt", skip = 8, data.table = F)

c22_filt <- c22 %>%
  filter(minor_allele_frequency > 0.002 & info > 0.3)

filter(c22, minor_allele_frequency > 0.001) %>% nrow %>% (function(x){x/nrow(c22)})
filter(c22, minor_allele_frequency > 0.002) %>% nrow %>% (function(x){x/nrow(c22)})
filter(c22, minor_allele_frequency > 0.002) %>% nrow %>% (function(x){x/nrow(c22)})
filter(c22, info > 0.3) %>% nrow %>% (function(x){x/nrow(c22)})
filter(c22, info > 0.4) %>% nrow %>% (function(x){x/nrow(c22)})

filter(c22, info > 0.3 & minor_allele_frequency > 0.002) %>% nrow %>% (function(x){x/nrow(c22)})
filter(c22, info > 0.3 & minor_allele_frequency > 0.005) %>% nrow %>% (function(x){x/nrow(c22)})
