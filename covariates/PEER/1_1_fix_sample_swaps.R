setwd("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/04_phase2_full_analysis/peer_factors/peer_InputFiles")

library(dplyr)
library(data.table)

phe <- fread(data.table = F, "GeneExpr_PEER_TmmInvRankNormalised.csv")
swaps <- fread(data.table = F, "SamplesSwaps_to_Change.csv")

removes <- swaps %>% filter(Status == "remove") %>% pull(BamFile_Old)
swap <- swaps %>% filter(Status == "swap") %>% pull(BamFile_Old)


phe2 <- phe %>% 
  filter(!(V1 %in% removes)) %>%
  mutate(V1 = ifelse(V1 %in% swap, 
                     swaps$BamFile_New[match(V1, swaps$BamFile_Old)],
                     V1)) %>% 
  arrange(V1)

fwrite(phe2, file = "GeneExpr_PEER_TmmInvRankNormalised_swapsSwapped_mismatchRemoved.csv")

# Covariates
cov <- fread("Covariates_for_PEER.csv", data.table = F)
cov2 <- cov %>%
  filter(!V1 %in% removes)
fwrite(cov2, file = "Covariates_for_PEER_mismatchRemoved.csv")


# Check order of phe matches order of cov
table(phe2$V1 == cov2$V1)