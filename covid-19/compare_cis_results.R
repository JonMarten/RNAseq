library(data.table)
library(dplyr)
library(ggplot2)
library(qqman)
library(cowplot)
library(tidyr)
theme_set(theme_cowplot())

setwd("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/covid19/results")

allNom <- fread("INTERVAL_ACE2_MAF0.005_cisNominal.cis_qtl_pairs.23.csv", data.table = F)
maleNom <- fread("INTERVAL_ACE2_MAF0.005_cisNominal_male.cis_qtl_pairs.23.csv", data.table = F)
femaleNom <- fread("INTERVAL_ACE2_MAF0.005_cisNominal_female.cis_qtl_pairs.23.csv", data.table = F)
nonNom <- fread("INTERVAL_ACE2_nonzeros_MAF0.005_cisNominal.cis_qtl_pairs.23.csv", data.table = F)

allCis <- fread("INTERVAL_ACE2_MAF0.005_cis.csv", data.table = F)
maleCis  <- fread("INTERVAL_ACE2_MAF0.005_cis_male.csv", data.table = F)
femaleCis  <- fread("INTERVAL_ACE2_MAF0.005_cis_female.csv", data.table = F)
nonCis <- fread("INTERVAL_ACE2_nonzeros_MAF0.005_cis.csv", data.table = F)

map <- fread("../genotypes/INTERVAL_chrX_merged_cleaned_RNAseq_phase1-2_deduplicated_MAF0.005.bim", data.table = F)
names(map) <- c("chr", "variant_id", "morg","bp","A1","A2")

allNom <- right_join(map, allNom, by = "variant_id")

merge <- full_join(maleNom, femaleNom, suffix = c(".male",".female"), by = "variant_id")
merge <- full_join(allNom, merge, by = "variant_id")
merge <- full_join(merge, nonNom, by = "variant_id", suffix = c("", ".nonzeros"))

merge2 <- merge %>%
  select(phenotype_id, chr, variant_id, rsid = rsid.x, bp, A1, A2, tss_distance, maf, maf.male, maf.female, maf.nonzeros, ma_samples, ma_samples.male, ma_samples.female, ma_samples.nonzeros, ma_count, ma_count.male, ma_count.female, ma_count.nonzeros, pval_nominal, pval_nominal.male, pval_nominal.female, pval_nominal.nonzeros, slope, slope.male, slope.female, slope.nonzeros, slope_se, slope_se.male, slope_se.female, slope_se.nonzeros)

# comparison plots
maf <- ggplot(merge2, aes(x = maf.male, y = maf.female)) +
  geom_point()
maf_non <- ggplot(merge2, aes(x = maf, y = maf.nonzeros)) +
  geom_point()

ma_samples <- ggplot(merge2, aes(x = ma_samples.male, y = ma_samples.female)) +
  geom_point()
ma_samples_non <- ggplot(merge2, aes(x = ma_samples, y = ma_samples.nonzeros)) +
  geom_point()

log10p <- ggplot(merge2, aes(x = -log10(pval_nominal.male), y = -log10(pval_nominal.female))) +
  geom_point()
log10p_non <- ggplot(merge2, aes(x = -log10(pval_nominal), y = -log10(pval_nominal.nonzeros))) +
  geom_point()

zscore <- ggplot(merge2, aes(x = slope.male/slope_se.male, y = slope.female/slope_se.female)) +
  geom_point()
zscore_non <- ggplot(merge2, aes(x = slope/slope_se, y = slope.nonzeros/slope_se.nonzeros)) +
  geom_point()


merge2 <- merge2 %>%
  mutate(sigCombined = ifelse(pval_nominal < 0.05/nrow(merge2), 1, 0),
         sig.male = ifelse(pval_nominal.male < 0.05/nrow(merge2), 1, 0),
         sig.female = ifelse(pval_nominal.female < 0.05/nrow(merge2), 1, 0),
         sig.non = ifelse(pval_nominal.nonzeros < 0.05/nrow(merge2), 1, 0))

manplot <- merge2 %>%
  filter(pval_nominal < 0.05  | pval_nominal.male < 0.05 | pval_nominal.female < 0.05) %>%
  select(chr, rsid, bp, A1, A2, pval_nominal, pval_nominal.male, pval_nominal.female, pval_nominal.nonzeros) %>%
  gather(-(chr:A2), key = "model", value = "pval") %>%
  mutate(kb = bp / 10^6)

ace2Start <- 15579156 / 10^6
ace2End <- 15620271 / 10^6
sigThresh = 0.05/nrow(merge2)

manhat <- ggplot(data = manplot, aes(x = kb, y = -log10(pval), colour = model)) +
  geom_point() +
  facet_wrap(~model, ncol = 1) +
  geom_vline(xintercept = ace2Start) + 
  geom_vline(xintercept = ace2End) +
  geom_hline(yintercept = -log10(sigThresh), colour = "grey70")

cov <- fread(data.table = F, "/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/covid19/covariates/INTERVAL_RNAseq_COVID19_covariates.txt") %>%
  t() %>%
  data.frame(stringsAsFactors = F)
names(cov) <- cov[1,] %>% as.character()
cov <- cov[-1,]

cov$sex <- as.numeric(cov$sex)

cowplot::save_plot("ACE2_cis_comparison_manhattan.png", manhat, base_height = 7, base_width = 15)
fwrite(merge2, "INTERVAL_ACE2_MAF0.005_cisNominal_mergedResults.csv")