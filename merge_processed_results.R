setwd( "U:/Projects/RNAseq/analysis/00_testing/results/eqtl_test_13633258/processed")
library(data.table)
library(dplyr)

files <- list.files()

for (file in files){
  
  # if the merged dataset doesn't exist, create it
  if (!exists("dat")){
    dat <- fread(file, data.table = F)
  }
  
  # if the merged dataset does exist, append to it
  if (exists("dat")){
    tempDat <-fread(file, data.table = F)
    dat <- rbind(dat, tempDat)
    rm(tempDat)
  }
  
}

dat2 <- dat %>%
  arrange(feature_id, snp_chromosome, snp_position) 

bonfthresh <- 5e-8 / nrow(dat2)

length(which(dat2$p_value < bonfthresh))


library(qqman)
dat3 <- dat2 %>% filter(p_value < 5e-5)
manhattan(dat3, chr = "snp_chromosome", bp = "snp_position", p = "p_value", genomewideline = -log10(bonfthresh))
