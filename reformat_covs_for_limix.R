# Reformat covariates for LIMIX pipeline
library(dplyr)
library(data.table)
setwd("U:/Projects/RNAseq")

dat <- fread("covariates/sarah_data_2/INTERVALdata_13MAY2019.csv", data.table=F)
datp3 <- fread("covariates/sarah_data_2/INTERVALdata_P3_13MAY2019.csv", data.table=F)
rna_id_mapper <- fread("rna_id_mapper.csv", data.table=F)

dat <- full_join(dat, datp3, by="identifier")
dat <- left_join(rna_id_mapper, dat)

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
covOut <- dat3 %>% 
  select(sample_id, sexPulse, age_RNA, batch, phase) %>% 
  filter(!is.na(phase) & !is.na(batch))

fwrite(sep = "\t", covOut, file = "covariates/INTERVAL_RNA_batch1_2_covariates_sex_age.txt", quote = F, na = NA)

## Generate phenotype file that only includes cell count data from the same phase as the RNA sample
pheAll <- dat3 %>%
  select(sample_id, age_RNA)
dat <- rename(dat, sample_id = RNA_id)
pheAll <- left_join(pheAll, dat)
pheAll <- pheAll %>%
  select(-identifier)

m24Cols <- names(pheAll)[grep("_24", names(pheAll))] 
m48Cols <- names(pheAll)[grep("_48", names(pheAll))] 
p3Cols <- names(pheAll)[grep("_p3", names(pheAll))] 
blCols <- names(pheAll)[grep("_bl", names(pheAll))] 
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
  filter(!is.na(batch))

### check for dupes:

dat %>% group_by(identifier) %>%
  filter(n()>1) %>%
  data.frame %>% 
  head


write.csv(allOut, file = "covariates/INTERVAL_RNA")
