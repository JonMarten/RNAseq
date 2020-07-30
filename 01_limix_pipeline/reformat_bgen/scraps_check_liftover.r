setwd("/home/jm2294/rds/rds-jmmh2-projects/interval_rna_seq/GENETIC_DATA")
library(dplyr)
library(data.table)
library(stringr)


vcflist <- list()

for(i in 1:22) {
  b38 <- fread(paste0("b37_b38_liftover/from_savita_INTERVAL_imputed_liftover_hg38_cleaned/INTERVAL_imputed_liftover_hg38_cleaned_chr",i,".vcf"), data.table = F, skip = 5)
  names(b38) <- c("CHROM.b38","POS.b38","match_id","REF.b38","ALT.b38")
  vcflist[[i]] <- b38
  rm(b38)
}

vcf <- bind_rows(vcflist, .id = "column_label")

b37 <- str_split_fixed(vcf$match_id, "_", 4)
b37_2 <- str_split_fixed(b37[,4], ":", 2)

b37.df <- data.frame(cbind(b37[,-4],b37_2), stringsAsFactors = F)
names(b37.df) <- c("CHROM.b37","POS.b37", "REF.b37","ALT.b37", "rsid")
b37.df <- b37.df %>%
  mutate(CHROM.b37 = as.numeric(CHROM.b37),
         POS.b37 = as.numeric(POS.b37))

liftover <- cbind(b37.df, vcf)
fwrite(liftover, file = "b37_b38_liftover/liftover_map.csv")

library(ggplot2)
library(cowplot)

ltest <- liftover[runif(1000, 1, nrow(liftover)),] %>%
  arrange(CHROM.b37, POS.b37)

ggplot(ltest, aes(x = POS.b37, y = POS.b38, colour = as.factor(CHROM.b37))) +
  geom_point() + 
  facet_wrap(vars(CHROM.b37)) +
  geom_abline(slope = 1, intercept = 0)


ltest2 <- ltest %>%
  mutate(posdif = abs(POS.b37 - POS.b38))

ggplot(ltest2, aes(x = POS.b37, y = posdif, colour = as.factor(CHROM.b37))) + 
  facet_wrap(vars(CHROM.b37)) + 
  geom_point()

#for(i in 1:22){
#  vcflist[[i]] %>% head(n=3) %>% print
#}

snpinfoList <- list()
for(i in 1:22) {
  snpstats <- fread(paste0("snp_stats/impute_",i,"_interval_snp_stats_unfiltered.txt"), skip = 8, data.table=F)
  snpinfoList[[i]] <- snpstats
  rm(snpstats)
}

for(i in 1:22){
  snpinfoList[[i]] %>% head(n=3) %>% print
}



