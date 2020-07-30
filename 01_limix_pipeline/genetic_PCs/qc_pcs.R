library(data.table)
library(ggplot2)
library(dplyr)
library(GGally)

setwd("U:/Projects/RNAseq/genetic_PCs")

pcs <- read.table(h = T, stringsAsFactors = F, "rnaseq_phase1.typed.maf1pc.geno3pc.hwe1minus5.complex-region-excl.indep-50-5-1.5.pruning.genome.pca.eigenvec") 

map <- read.table(h = T, file="U:/Projects/RNAseq/analysis/00_testing/phenotype/sample_mapping_file_gt_to_phe_phase1.txt") %>%
  rename(IID = genotype_individual_id, sample_id = phenotype_individual_id)

pcs2 <- left_join(pcs, map)

covs <- fread(data.table = F, "../covariates/INTERVAL_RNA_batch1-8_covariates_release_2019_08_15.csv")

dat <- full_join(covs, pcs2)

datplot <- dat %>%
  select(PC1:PC4, sex, sequencingBatch, age_RNA, BMI, ethnicity) %>%
  mutate(sex = as.factor(sex), sequencingBatch = as.factor(sequencingBatch), ethnicity = as.factor(ethnicity))

ggpairs(datplot)

ggplot(datplot, aes(x = sequencingBatch, y = age_RNA)) + geom_boxplot()




## pc plots
pc3_outliers <- dat %>%
  arrange(-PC3) %>%
  head(2) %>%
  pull(sample_id)

dat <- dat %>%
  mutate(PC3out = ifelse(sample_id %in% pc3_outliers, "red", "black"))

ggplot(dat, aes(x = PC1)) + geom_density(fill = "red")

pclong <- gather(pcs, key = "PC", value = "value", -IID, -FID) %>%
  mutate(PC = factor(PC, levels = paste0("PC",1:20)))
ggplot(pclong, aes(x = value, fill = PC)) + 
  geom_density(alpha= 0.8) +
  facet_wrap(.~PC,scales = "free")

p1 <- ggplot(dat, aes(x = PC1, y = PC2, colour = PC3out)) + geom_point()
p2 <- ggplot(dat, aes(x = PC3, y = PC4, colour = PC3out)) + geom_point()
p3 <- ggplot(dat, aes(x = PC5, y = PC6, colour = PC3out)) + geom_point()
p4 <- ggplot(dat, aes(x = PC7, y = PC8, colour = PC3out)) + geom_point()
p5 <- ggplot(dat, aes(x = PC9, y = PC10, colour = PC3out)) + geom_point()
p6 <- ggplot(dat, aes(x = PC11, y = PC12, colour = PC3out)) + geom_point()
p7 <- ggplot(dat, aes(x = PC13, y = PC14, colour = PC3out)) + geom_point()
p8 <- ggplot(dat, aes(x = PC15, y = PC16, colour = PC3out)) + geom_point()
p9 <- ggplot(dat, aes(x = PC17, y = PC18, colour = PC3out)) + geom_point()
p10 <- ggplot(dat, aes(x = PC19, y = PC20, colour = PC3out)) + geom_point()



plot_grid(plotlist = list(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10))
                