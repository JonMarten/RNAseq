# Review VEP annotated results
library(data.table)
library(ggplot2)
library(stringr)
library(cowplot)
library(dplyr)

setwd("U:/Projects/RNAseq/analysis/00_testing/test_parameters")

dat <- fread("results_merged_chr22_window500000_Perm100_MAF0.01_eSNPs_vepannotated.csv", data.table = F)

# list un-annoted SNPs
dat.dropped <- dat %>% filter(is.na(IMPACT) | IMPACT == "")

tab <- dat %>%
  group_by(Consequence) %>%
  summarise(count = n()) %>%
  data.frame() %>%
  mutate(Consequence = as.factor(Consequence))

ggplot(dat, aes(x  = Consequence)) +
  geom_bar()+ 
  coord_flip()

ggplot(dat, aes(x = reorder(tab, -count))) +
  geom_bar(stat = "identity") + 
  coord_flip()
