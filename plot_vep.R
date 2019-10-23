# Review VEP annotated results
library(data.table)
library(ggplot2)
library(stringr)
library(cowplot)
library(dplyr)
library(UpSetR)

setwd("U:/Projects/RNAseq/analysis/00_testing/test_parameters")

dat <- fread("results_merged_chr22_window500000_Perm100_MAF0.01_eSNPs_vepannotated.csv", data.table = F)

# list un-annoted SNPs
dat.dropped <- dat %>% filter(is.na(Existing_variation) | Existing_variation == "")
dat <- dat %>% filter(!is.na(Existing_variation) & Existing_variation != "")


tab <- dat %>%
  group_by(Consequence) %>%
  summarise(count = n()) %>%
  data.frame() %>%
  mutate(Consequence = as.factor(Consequence)) %>% 
  arrange(desc(count))

tab$Consequence <- gsub(pattern = ",", "&", tab$Consequence)

# Make upset plot for consequences 
a <- c(tab$count)
names(a) <- tab$Consequence
upset(fromExpression(a), 
      nsets = 15,
      order.by = "freq",
      mb.ratio = c(0.60, 0.40),
      text.scale = 1.5)

## Annotate with intron number

introns <- str_split_fixed(dat$INTRON, "/", 2) %>%
  data.frame(stringsAsFactors = F)
names(introns) <- c("variant_intron_no", "total_no_introns")
introns <- introns %>%
  mutate(variant_intron_no = as.numeric(variant_intron_no),
         total_no_introns = as.numeric(total_no_introns)) %>%
  mutate(intron_pct = variant_intron_no / total_no_introns) %>%
  mutate(total_no_introns_trunc = ifelse(total_no_introns < 12, total_no_introns, 12))

dat2 <- cbind(dat, introns)      

hist(dat2$variant_intron_no)

ggplot(filter(unique(dat2, by = "rsid"), !is.na(variant_intron_no)), aes(x = as.factor(variant_intron_no))) +
  geom_bar() +
  facet_wrap(.~total_no_introns_trunc, scales = "free")

eqtl <- fread("results_merged_chr22_window500000_Perm100_MAF0.01_eSNPs_eqtlgenReplication.txt", data.table = F)

dat$merge_id <- paste0(dat$feature_id, ":", dat$rsid)
eqtl$merge_id <- paste0(eqtl$limix_feature_id, ":", eqtl$limix_rsid)

m <- inner_join(dat, eqtl)

ggplot(m, aes(x = limix_zscore, y = eQTLgen_Zscore_flipped, color = IMPACT)) + 
  geom_vline(xintercept = 0) + 
  geom_hline(yintercept = 0) + 
  geom_point(cex = 0.8, alpha = 0.8, pch = 20) +
  theme(aspect.ratio = 1, legend.position = "right") +
  labs(x = "Limix Z-score", y = "eQTLgen Z-score") 

# get list of features with more than two SNPs mismatched  
m <- m %>%
  mutate(signMismatch = ifelse(sign(limix_zscore) != sign(eQTLgen_Zscore_flipped), 1, 0))

flexmismatchGenes <- m %>%
  filter(signMismatch == 1) %>%
  group_by(limix_feature_id) %>%
  summarise(n_mismatched = n()) %>%
  data.frame %>%
  filter(n_mismatched > 2) %>%
  pull(limix_feature_id)

plotm <- filter(m, limix_feature_id %in% flexmismatchGenes)
plotvec <- runif(n = 10000, 1, nrow(plotm))
plotm2 <- plotm[plotvec,]

ggplot(plotm, aes(x = limix_zscore, y = eQTLgen_Zscore_flipped, color = BIOTYPE)) + 
  geom_vline(xintercept = 0) + 
  geom_hline(yintercept = 0) + 
  geom_point(cex = 1, alpha = 0.8, pch = 20) +
  theme(aspect.ratio = 1, legend.position = "right") +
  labs(x = "Limix Z-score", y = "eQTLgen Z-score") +
  theme_minimal() +
  facet_wrap(. ~ limix_gene_name) #+
  #scale_colour_gradient(low = "green", high = "red")
  #scale_color_viridis_c(option = "A") +
  #theme(panel.background = element_rect(fill = "gray50"))
