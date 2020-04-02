setwd("C:/Users/Jonathan/Documents/INTERVAL_covars")

library(data.table)
library(dplyr)

techCovsPhase1 <- fread("INTERVAL_RNA_Covariates_09-08-2019.csv", data.table = F)
techCovsPhase2 <- fread("INTERVAL_RNA_Covariates_Ph2_27-03-2020.csv", data.table = F)

techCovsPhase2 <- techCovsPhase2 %>% 
  mutate_at(c("Sample_ID", "Lane", "Freezer_Shelf", "Agilent_RINe"), as.character)
  
  
  
  mutate(Sample_ID = as.character(Sample_ID),
         Lane = as.character(Lane))

techCovsAll <- bind_rows(techCovsPhase1, techCovsPhase2)

fwrite(techCovsAll, "INTERVAL_RNA_technical_covariates_batch1-12_20200402.csv")
