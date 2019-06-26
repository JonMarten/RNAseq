# Script to make a table of all INTERVAL participants detailing omics data available on them
library(dplyr)
dataRelease <- "11JUN2019"
setwd(paste0("U:/Projects/RNAseq/RNA_sample_selection/", dataRelease))
ids <- read.csv("omicsMap.csv", stringsAsFactors = F)
ids3 <- read.csv("omicsMap_P3.csv", stringsAsFactors = F)
RNApick <- read.csv("RNAPick.csv", stringsAsFactors = F)
agesex <- read.csv(paste0("INTERVALdata_", dataRelease, ".csv"), stringsAsFactors = F)

ids[which(ids == "", arr.ind=T)] <- NA
ids3[which(ids3 == "", arr.ind=T)] <- NA
allids <- full_join(ids, ids3)

allids2 <- full_join(allids, RNApick)

dat <- full_join(agesex, allids2)


allids %>% select(grep("Aff",names(allids))) %>% head
mismatch <- allids2[which(!allids2$identifier %in% allids$identifier),]

dat <- dat[,order(colnames(dat))]

#dat[which(dat$Affymetrix_gwasQC_24m != dat$Affymetrix_gwasQC_bl),]
#dat[which(dat$Affymetrix_gwasQC_24m == dat$Affymetrix_QC_24m),]

dat2 <-dat %>% 
  rowwise() %>%
  mutate(affymetrix_ID = c(Affymetrix_QC_bl, Affymetrix_gwasQC_bl, Affymetrix_gwasQC_24m, Affymetrix_QC_24m) %>% unique %>% na.exclude %>% paste(collapse = ";"),
         brainshake_ID = c(Brainshake_QC_bl, Brainshake_gwasQC_bl) %>% unique %>% na.exclude %>% paste(collapse = ";"),
         metabolon_ID = c(Metabolon_met_gwasQC_bl, Metabolon_met_QC_bl) %>% unique %>% na.exclude %>% paste(collapse = ";"),
         olink_ID = c(Olink_cvd2_gwasQC_24m, Olink_cvd2_QC_24m, Olink_cvd3_gwasQC_24m, Olink_cvd3_QC_24m, Olink_inf_gwasQC_24m, Olink_inf_QC_24m, Olink_neu_gwasQC_24m, Olink_neu_QC_24m) %>% 
           unique %>% na.exclude %>% paste(collapse = ";"),
         olink_CVD2 = ifelse(is.na(Olink_cvd2_gwasQC_24m) & is.na(Olink_cvd2_QC_24m), 0, 1),
         olink_CVD3 = ifelse(is.na(Olink_cvd3_gwasQC_24m) & is.na(Olink_cvd3_QC_24m), 0, 1),
         olink_INF = ifelse(is.na(Olink_inf_gwasQC_24m) & is.na(Olink_inf_QC_24m), 0, 1),
         olink_NEU = ifelse(is.na(Olink_neu_gwasQC_24m) & is.na(Olink_neu_QC_24m), 0, 1),
         RNA_ID = c(RNAseq_gwasQC_24m, RNAseq_gwasQC_48m, RNAseq_gwasQC_p3, RNAseq_QC_24m, RNAseq_QC_48m, RNAseq_QC_p3, RNAseq_RAW_24m, RNAseq_RAW_48m, RNAseq_RAW_p3) %>% unique %>% na.exclude %>% paste(collapse = ";"),
         soma_ID = c(soma4000_gwasQC_bl, soma4000_QC_bl) %>% unique %>% na.exclude %>% paste(collapse = ";"),
         WES_ID = Wes_gwasQC_bl,
         WGS_ID = c(Wgs_QC_bl, Wgs_gwasQC_bl, Wgs_QC_24m, Wgs_gwasQC_24m) %>% unique %>% na.exclude %>% paste(collapse = ";")
         ) %>%
  data.frame
dat2[which(dat2 == "", arr.ind = T)] <- NA

#Count non-missing ids
apply(dat2, MARGIN = 2, FUN = function(x){length(which(!is.na(x)))})

dat3 <- dat2 %>% select(data_management_project_id = identifier, agePulse, NIHRConsent, picklistRNA, sexPulse, Tempus2Y, TempusP3, Tempusp4, affymetrix_ID, brainshake_ID, metabolon_ID, olink_ID, olink_CVD2, olink_CVD3, olink_INF, olink_NEU, RNA_ID, soma_ID, WES_ID, WGS_ID) %>%
  rowwise() %>%
  mutate(tempus = ifelse(!is.na(Tempus2Y) | !is.na(TempusP3) | !is.na(Tempusp4),1,0)) %>%
  mutate(RNA_pick_priority_1 = ifelse(picklistRNA == 1, ifelse(is.na(RNA_ID), "not_yet_assigned", RNA_ID), NA),
         RNA_pick_priority_2 = ifelse(picklistRNA == 2, ifelse(is.na(RNA_ID), "not_yet_assigned", RNA_ID), NA),
         RNA_pick_priority_3 = ifelse(picklistRNA == 3, ifelse(is.na(RNA_ID), "not_yet_assigned", RNA_ID), NA),
         RNA_pick_priority_4 = ifelse(picklistRNA == 4, ifelse(is.na(RNA_ID), "not_yet_assigned", RNA_ID), NA)) %>%
  data.frame() 

out <- dat3 %>%
  select(data_management_project_id, 
         sex = sexPulse, 
         age_at_baseline = agePulse,
         NIHRConsent,
         affymetrix_ID,
         WES_ID,
         WGS_ID,
         brainshake_ID,
         metabolon_ID,
         soma_ID,
         olink_ID:olink_NEU,
         tempus,
         picklistRNA,
         RNA_ID,
         RNA_pick_priority_1:RNA_pick_priority_4,
         olink_ID)

write.csv(out, row.names = F, quote = T, file = paste0("../INTERVAL_omics_table_",dataRelease,".csv"))
