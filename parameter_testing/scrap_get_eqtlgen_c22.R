library(data.table)
library(dplyr)

setwd("/rds/user/jm2294/rds-jmmh2-projects/interval_rna_seq/analysis/00_testing/test_parameters")

# Get chr22 eqtls from eqtlgen
Run below code once
eQTL <- fread("zcat eQTLgen/cis-eQTLs_full_20180905.txt.gz", data.table = F)
names(eQTL) <- paste0("eQTLgen_", names(eQTL))
e22 <- eQTL %>%
  filter(eQTLgen_SNPChr == 22) %>%
  mutate(matchID = paste0(eQTLgen_Gene, "_", eQTLgen_SNP))
fwrite(e22, file = "eQTLgen/cis-eQTLs_chr22.txt")
