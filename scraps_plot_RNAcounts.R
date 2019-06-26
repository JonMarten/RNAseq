setwd("U:/Projects/RNAseq/test_run_chunks/output_new_pipeline")
rm(list=ls())
library(dplyr)
library(data.table)
library(qqman)

sumstats <- fread("testqtl_results_22_10736171_17730855.txt", data.table = F)

sumtab <- sumstats %>%
  group_by(feature_id) %>%
  summarise(nsnps = n(),
            maxp = max(p_value),
            minp = min(p_value)) %>%
  data.frame
 
manhattan(filter(sumstats, p_value < 0.001), chr = "snp_chromosome", bp = "snp_position", p = "p_value")

rawdat <- fread("../../test_run/phenotype_5281-fc-genecounts.txt", data.table = F)
dat <- dsumstats %>%
  select(-ENSEMBL_ID) %>%
  transpose
names(dat) <- rawdat[,1]
rownames(dat) <- names(rawdat)[-1]
library(ggplot2)
ggplot(dat, aes(x = ENSG00000000003)) + geom_density()
library(tidyr)
library(tibble)
datTall <- dat %>%
  rownames_to_column("id") %>%
  gather(dat, key = "Feature", value = "Expression", -id)
ggplot(datTall, aes(x = ENSG00000000003)) + geom_density(colour = )
