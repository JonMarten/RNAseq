# Take plink variant list as input and output RSIDs which are duplicated in the bgen file
setwd("/home/jm2294/GENETIC_DATA/INTERVAL/RNAseq")
library(dplyr)
library(data.table)

for(i in 1:22){
  cat("\nChr", i,"\n")
  a <- fread(paste0("plink/impute_",i,"_b37_unfiltered.pvar"), data.table = F)
  names(a)[1] <- "CHROM"
  a <- a %>%
    mutate(altid = paste0(CHROM, "_", POS, "_", REF,"_",ALT,":",ID))
  
  b <- a %>% 
    group_by(altid) %>% filter(n() > 1) %>% data.frame
  
  dupes <- b %>%
    distinct() %>%
    mutate(chrpos = paste0(CHROM,":",POS)) %>%
    pull(chrpos)
  
  nfilt <- a %>% filter(POS %in% b$POS) %>% nrow
  cat(paste0("Identified ",length(dupes), " duplicated variants on ",nrow(b)," rows. Filter will remove ", nfilt, " rows.\n" ))
  
  write.table(dupes, quote = F, col.names = F, row.names = F, file = paste0("plink/impute_",i,"_duplicate_positions.txt"))
  rm(a,b,dupes)
}

