# Reformat covariates for LIMIX pipeline
library(dplyr)
library(data.table)
setwd("U:/Projects/RNAseq/covariates")
rna_id_mapper <- fread("rna_id_mapper.csv", data.table=F) 
techCov <- fread("rna_technical_covariates_may19.csv", data.table = F)
techCov <- techCov %>% 
  filter(SANGER_Sample_ID != "") %>%
  select(sample_id = SANGER_Sample_ID,INTERVAL_Box,RNA_Extraction_Date, RIN, Normalization_Plate_ID, seq_run_id, seq_tag_index) %>%
  mutate(RNA_Extraction_Date = as.Date(RNA_Extraction_Date, "%d/%m/%Y"))



dat <- fread("data_release_20190531/INTERVALdata_31MAY2019.csv", data.table=F)

# Filter phase 3 data to get the correct timepoint to match sample used for RNA seq data
p3mapper <- rna_id_mapper %>%
  filter(phase == "p3" & !is.na(attendanceDate_p3)) %>%
  mutate(p3_mapper = paste0(identifier,attendanceDate_p3)) %>%
  pull(p3_mapper)

datp3 <- fread("data_release_20190531/INTERVALdata_P3_31MAY2019.csv", data.table=F)
datp3 <- datp3 %>%
  mutate(p3_mapper = paste0(identifier,attendanceDate_p3)) %>%
  filter(p3_mapper %in% p3mapper) %>%
  select(-p3_mapper)

# Merge together all phases and merge in INT_RNA identifiers
dat <- full_join(dat, datp3, by="identifier")
rna_id_mapper <- rna_id_mapper %>% 
  select(-attendanceDate_p3)
dat <- left_join(rna_id_mapper, dat, by="identifier")

# Find all columns with an attendance date
names(dat)[grep("att", names(dat))]

calcAge <- function(dob, testDate){
  (difftime(testDate, dob, unit = "weeks")/52.25) %>% as.numeric %>% round(digits = 1)
}

dat2 <- dat %>%
  select(sample_id = RNA_id, batch, phase, sexPulse, agePulse, monthPulse, yearPulse, attendanceDate, attendanceDate_24m, attendanceDate_48m, attendanceDate_p3) %>% 
  mutate_at(vars(attendanceDate:attendanceDate_p3), funs(as.Date(x = .,format = "%d%b%Y"))) %>%
  mutate(dob = paste0("15-",monthPulse, "-", yearPulse) %>% as.Date(format = "%d-%m-%Y")) %>%
  mutate(age1 = calcAge(dob, attendanceDate),
         age24m = calcAge(dob, attendanceDate_24m),
         age48m = calcAge(dob, attendanceDate_48m),
         agep3 = calcAge(dob, attendanceDate_p3))

dat3 <- dat2 %>%
  mutate(age_RNA = phase) %>%
  mutate(age_RNA = ifelse(phase =="24m", age24m, 
                          ifelse(phase=="48m", age48m,
                                 ifelse(phase=="p3", agep3, NA)))) %>%
  mutate(ageDif = age_RNA - age1)

## Write out test covariate file 
#covOut <- dat3 %>% 
#  select(sample_id, sexPulse, age_RNA, batch, phase) %>% 
#  filter(!is.na(phase) & !is.na(batch))
#
#fwrite(sep = "\t", covOut, file = "INTERVAL_RNA_batch1_2_covariates_sex_age.txt", quote = F, #na = NA)

## Generate phenotype file that only includes cell count data from the same phase as the RNA sample
pheAll <- dat3 %>%
  select(sample_id, age_RNA)
dat <- rename(dat, sample_id = RNA_id)
pheAll <- left_join(pheAll, dat)
pheAll <- pheAll %>%
  select(-identifier)

m24Cols <- names(pheAll)[grep("_24", names(pheAll))] %>% sort
m48Cols <- names(pheAll)[grep("_48", names(pheAll))] %>% sort
p3Cols <- names(pheAll)[grep("_p3", names(pheAll))] %>% sort
blCols <- names(pheAll)[grep("_bl", names(pheAll))] %>% sort
otherCols <- names(pheAll)[which(!names(pheAll) %in% c(blCols, m24Cols, m48Cols, p3Cols))]

# Check cols are the same other than suffix
identical(
  gsub("_24m", "", m24Cols),
  gsub("_48m", "", m48Cols)
)
identical(
  gsub("_24m", "", m24Cols),
  gsub("_p3", "", p3Cols)
)

phe24 <- pheAll %>% 
  filter(phase == "24m") %>%
  select(otherCols, m24Cols) %>%
  rename_at(vars(m24Cols), funs(gsub("_24m","___RNA",.)))

phe48 <- pheAll %>% 
  filter(phase == "48m") %>%
  select(otherCols, m48Cols) %>%
  rename_at(vars(m48Cols), funs(gsub("_48m","___RNA",.)))

phep3 <- pheAll %>% 
  filter(phase == "p3") %>%
  select(otherCols, p3Cols) %>%
  rename_at(vars(p3Cols), funs(gsub("_p3","___RNA",.)))

pheRNA <- rbind(phe24, phe48, phep3) %>%
  filter(!is.na(batch)) %>%
  arrange(sample_id) %>%
  rename(sex=sexPulse, ethnicity=ethnicPulse, sequencingBatch = batch, intervalPhase = phase) %>%
  select(-monthPulse, -yearPulse, -agePulse, -attendanceDate, -appointmentTime, -outCome) 
pheRNA$appointmentTime___RNA[which(pheRNA$appointmentTime___RNA == "")] <- NA # Set missing appointment times to NA

# merge in technical covariates
out <- full_join(techCov, pheRNA, by = "sample_id") %>%
  mutate(attendanceDate___RNA = as.Date(attendanceDate___RNA, format = "%d%b%Y"),
         processDate___RNA = as.Date(processDate___RNA, format = "%d%b%Y"))

write.csv(out, file = "INTERVAL_RNA_batch1_4_covariates_release31MAY2019.csv", row.names = F)
