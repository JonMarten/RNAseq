# Compare xcell enrichment to actual sysmex cell data
#https://github.com/dviraran/xCell

library(xCell)
library(data.table)
library(dplyr)
library(corrplot)
library(Hmisc)

setwd("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/phenotype/xCell")

phe <- fread("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/phenotype/gene_expr_filtered_TMM_rpkms_genesInvRankTransf_2747Samp_18373genes.csv", data.table = F)


dupegenes <- phe %>% 
  group_by(gene_symbol) %>%
  filter(n() > 1) %>%
  data.frame()

expr <- phe %>% 
  filter(!(gene_id %in% dupegenes$gene_id))
rownames(expr) <- expr$gene_symbol

expr <- expr %>%
  select(-(gene_id:gene_length))

predCells <- xCellAnalysis(expr)
save(predCells, file = "predCells.Rdata")

# Read in sysmex data
phenos <- fread("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/covariates/INTERVAL_RNA_batch1-8_covariates_release_2019_08_15.csv", data.table = F)
phenos <- phenos %>% 
  select(-(Tag_1:processTime___RNA))

predCells.df <- predCells %>%
  t() %>%
  data.frame(stringsAsFactors = F) %>% 
  mutate(sample_id = rownames(.))

merge <- full_join(phenos, predCells.df)
merge.m <- merge %>%
  select(-sample_id) %>%
  as.matrix()
library(Hmisc)
cormat <- rcorr(merge.m, type = "spearman")
sysmex.vec <- grep("___RNA", colnames(cormat$r))
halfcormat <- cormat$r[sysmex.vec, -sysmex.vec] 
halfpmat <- cormat$p[sysmex.vec, -sysmex.vec] 
corrplot(halfcormat, tl.col = "black", tl.cex = 0.55, tl.srt = 40, p.mat = halfpmat, sig.level = 0.05, insig = "blank")


#
library(tidyr)
library(cowplot)
theme_set(theme_cowplot())
predCells.tall <- predCells %>%
  data.frame %>%
  mutate(cellType = rownames(.)) %>% 
  gather(-cellType, key = "sample_id", value = "cellEnrichment")

ggplot(predCells.tall, aes(x = cellType, y = cellEnrichment)) +
  geom_boxplot()
ggplot(predCells.tall, aes(x = cellEnrichment, fill = cellType))+
  geom_density() +
  facet_wrap(. ~ cellType, scales = "free") + 
  theme(legend.position = "blank")
  
predCells.tall %>% 
  group_by(cellType) %>%
  summarise(mean = mean(cell_enrichment))

bloodNames <- fread("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/covariates/blood_trait_names.csv", data.table = F)

a <- rownames(halfcormat)
match(a, bloodNames$"Variable_name")
rownames(halfcormat) <- bloodNames$astleName[match(a, bloodNames$"Variable_name")]
vec <-  which(rownames(halfcormat) %in% (bloodNames %>% filter(Use == 1) %>% pull(astleName)))

plotcormat <- halfcormat[vec,] %>% as.matrix
colnames(plotcormat) <- rownames(predCells)
plotpmat <- halfpmat[vec,]

col<- colorRampPalette(c("navy", "white", "darkred"))(30)

png("corplot.png", height = 10, width = 10, units = "in")
corrplot(plotcormat, tl.col = "black", tl.cex = 0.7, tl.srt = 40, p.mat = plotpmat, sig.level = 0.05, insig = "blank", col=col)
dev.off()

