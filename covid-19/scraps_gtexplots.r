init()
library(tidyr)
cisNom <- fread("INTERVAL_ACE2_MAF0.005_cisNominal_mergedResults.csv", data.table = F)
gtex <- fread("/home/jm2294/rds/rds-jmmh2-projects/covid/ace2/GTEx/gtex_v8_conditionally_independent_ACE2_eQTL_variants_withRSID.csv", data.table = F)
gtexM <- left_join(gtex, cisNom, by = "rsid", suffix = c(".gtex", ".interval"))

slopeplot <- ggplot(gtexM, aes(x = slope.gtex, y = slope.interval, colour = tissue)) + 
  geom_point() + 
  facet_wrap(~tissue_group) +
  geom_hline(yintercept = 0, colour = "grey70") + 
  geom_vline(xintercept = 0, colour = "grey70")


dat <- gtexM %>% 
  select(tss_distance.gtex, pval_nominal.gtex, pval_nominal.interval, tissue_group) %>%
  gather(-tss_distance.gtex, -tissue_group, key = "data", value = "pval_nominal")

thresh <- 0.05/length(unique(gtexM$rsid))

manplot <- ggplot(dat, aes(x = tss_distance.gtex, y = -log10(pval_nominal), colour = tissue_group)) +
  geom_point() +
  facet_wrap(~data, ncol = 1) +
  geom_hline(yintercept = thresh)

bigM <- full_join(gtex, cisNom, by = "rsid", suffix = c(".gtex", ".interval")) %>%
  mutate(gtexHit = ifelse(is.na(pval_nominal.gtex), 0, 1))

bigDat <- bigM %>%
  filter(maf.interval > 0.1) %>%
  select(bp, pval_nominal.gtex, pval_nominal.interval, tissue_group, gtexHit) %>%
  gather(-bp, -tissue_group, -gtexHit, key = "data", value = "pval_nominal")

fullmanplot <- ggplot(bigDat, aes(x = bp, y = -log10(pval_nominal), colour = as.factor(gtexHit))) +
  geom_point(alpha = 0.5) +
  facet_wrap(~data, ncol = 1) +
  geom_hline(yintercept = thresh)

ggsave(slopeplot, file = "slopeplot.png")
ggsave(manplot, file = "manplot.png")
ggsave(fullmanplot, file = "fullmanplot.png")