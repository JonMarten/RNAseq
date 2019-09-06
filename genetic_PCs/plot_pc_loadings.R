# Plot PC loadings to look for unfiltered complex regions of the genome

setwd("U:/Projects/RNAseq/genetic_PCs/")

library(dplyr)
library(data.table)
library(ggplot2)
library(cowplot)
library(tidyr)
library(reshape2)

loadings <- fread(data.table = F, "rnaseq_phase1.typed.maf1pc.geno3pc.hwe1minus5.complex-region-excl.indep-50-5-1.5.pruning.genome.pca.eigenvec.var")
names(loadings)[1] <- "CHROM"
loadings <- loadings %>%
  arrange(CHROM, POS) %>%
  mutate(CHROM = as.factor(CHROM))

ggplot(loadings, aes(x = POS, y = PC3, colour = CHROM)) +
  geom_point() +
  facet_wrap(. ~ CHROM) 

chrstarts <- loadings %>% 
  group_by(CHROM) %>%
  summarise(CHRstart = min(POS),
            CHRend = max(POS)) %>% 
  data.frame() %>%
  mutate(cumend = cumsum(CHRend))

loadings$plotpos <- NA
for(i in 1:22){
  rows <- which(loadings$CHROM == i)
  if(i == 1){
    loadings$plotpos[rows] <- loadings$POS[rows]
  } else {
    chrend <- chrstarts$cumend[i - 1]
    loadings$plotpos[rows] <- loadings$POS[rows] + chrend
  }
}

ggplot(loadings, aes(x = plotpos, y = PC3, colour = CHROM)) +
  geom_point() +
  facet_wrap(. ~ CHROM) 

pcs <- fread("rnaseq_phase1.typed.maf1pc.geno3pc.hwe1minus5.complex-region-excl.indep-50-5-1.5.pruning.genome.pca.eigenvec", data.table = F)

load_tall <- gather(loadings, key = "PC", value = "loading", -CHROM, -POS, -ID, -MAJ, -NONMAJ, -plotpos)
load_tall <- load_tall %>%
  mutate(PC = factor(PC, levels = paste0("PC", 1:20)))

lplot <- load_tall %>%
  filter(as.character(PC) %in% paste0("PC", 1:5)) %>%
  mutate(abs_loading = abs(loading))
  filter(abs_loading > 1)

g <- ggplot(lplot, aes(x = plotpos, y = abs_loading, colour = CHROM)) +
  geom_point(cex = 0.8) +
  facet_wrap(. ~ PC) 

save_plot(g, file = "pc_loadings.png", base_width = 6, base_height = 4)

ggplot(pcs, aes(x = PC1, y = PC2)) + geom_point()
ggplot(pcs, aes(x = PC3, y = PC4)) + geom_point()
ggplot(pcs, aes(x = PC5, y = PC6)) + geom_point()

kinmat <- fread("rnaseq_phase1.typed.maf1pc.geno3pc.hwe1minus5.complex-region-excl.indep-50-5-1.5.pruning.genome.pca.rel", data.table = F)
kinmatIDs <- fread("rnaseq_phase1.typed.maf1pc.geno3pc.hwe1minus5.complex-region-excl.indep-50-5-1.5.pruning.genome.pca.rel.id", data.table = F)

kinmatIDs$IID <- as.character(kinmatIDs$IID)

map <- read.table(h = T, file="U:/Projects/RNAseq/analysis/00_testing/phenotype/sample_mapping_file_gt_to_phe_phase1.txt") %>%
  rename(IID = genotype_individual_id, sample_id = phenotype_individual_id) %>%
  mutate(IID = as.character(IID))

map2 <- left_join(kinmatIDs, map)

rownames(kinmat)  <- map2$sample_id
names(kinmat) <- map2$sample_id
kinmat$id <- map2$sample_id

kinlong <- melt(kinmat)

kinlong %>% 
  mutate(id = as.character(id), variable = as.character(variable)) %>% 
  filter(id != variable) %>%st?kes1347

  arrange(desc(value)) %>% head

kinlong %>% filter(!id %in% c("INT_RNA7711096","INT_RNA7427299")) %>% pull(value) %>% summary() 
