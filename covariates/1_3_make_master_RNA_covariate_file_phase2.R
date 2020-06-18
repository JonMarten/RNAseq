# Reformat covariates for LIMIX pipeline
library(dplyr)
library(data.table)
setwd("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/04_phase2_full_analysis/covariates")
rna_id_mapper <- fread("processed/rna_id_mapper.csv", data.table=F) 
techCov <- fread("processed/INTERVAL_RNA_technical_covariates_batch1-12_20200416.csv", data.table = F)

techCov <- techCov %>% 
  rename(sample_id = INT_ID) %>%
  mutate(Extraction_Date = as.Date(Extraction_Date, "%d/%m/%Y"))

dat <- fread("raw/INTERVALdata_14MAY2020.csv", data.table=F)

# Filter phase 3 data to get the correct timepoint to match sample used for RNA seq data
p3mapper <- rna_id_mapper %>%
  filter(phase == "p3" & !is.na(attendanceDate_p3)) %>%
  mutate(p3_mapper = paste0(identifier,attendanceDate_p3)) %>%
  pull(p3_mapper)

datp3 <- fread("raw/INTERVALdata_P3_14MAY2020.csv", data.table=F)
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
  select(sample_id = RNA_id, inFeatureCounts, Batch, phase, sexPulse, agePulse, monthPulse, yearPulse, attendanceDate, attendanceDate_24m, attendanceDate_48m, attendanceDate_p3) %>% 
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

## Generate phenotype file that only includes cell count data from the same phase as the RNA sample
pheAll <- dat3 %>%
  select(sample_id, age_RNA, inFeatureCounts)
dat <- rename(dat, 
              sample_id = RNA_id,
              height = ht_bl,
              weight = wt_bl)
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
  #filter(!is.na(batch)) %>%
  arrange(sample_id) %>%
  rename(sex=sexPulse, ethnicity=ethnicPulse, sequencingBatch = Batch, intervalPhase = phase) %>%
  select(-monthPulse, -yearPulse, -agePulse, -attendanceDate, -appointmentTime, -outCome) 
pheRNA$appointmentTime___RNA[which(pheRNA$appointmentTime___RNA == "")] <- NA # Set missing appointment times to NA

# merge in technical covariates
out <- full_join(techCov, pheRNA, by = "sample_id") %>%
  mutate(attendanceDate___RNA = as.Date(attendanceDate___RNA, format = "%d%b%Y"),
         processDate___RNA = as.Date(processDate___RNA, format = "%d%b%Y")) %>%
  select(-sequencingBatch) %>%
  rename(sequencingBatch = Batch)

# Bodgy hack to re-add sequencing batch for samples that have RNA seq data but haven't been added to the phenotype database yet. 
fullbatch <- rna_id_mapper %>%
  select(sample_id=RNA_id, fullbatch = Batch)
out2 <- full_join(out, fullbatch) %>%
  select(-sequencingBatch) %>%
  rename(sequencingBatch = fullbatch)

# Process height asd weight data to remove unlikely values and significant figures and calculate BMI. Add in column for sample collection month, then reorder columns for saving.

# function to set anyone more than 3 SD away from the mean to the median.
fixExtremes <- function(x) {
  minx = mean(x, na.rm = T) - 3 * sd(x, na.rm = T)
  maxx = mean(x, na.rm = T) + 3 * sd(x, na.rm = T)
  medx = median(x[which(x > minx & x < maxx)])
  y <- x
  y[which(y < minx)] <- medx
  y[which(y > maxx)] <- medx
  return(y)
}

# function to report numbers of people outside these extremes
diagExtremes <- function(x, n) {
  print(paste0("Threshold: ", n, " SDs away from the mean"))
  minx = mean(x, na.rm = T) - n * sd(x, na.rm = T)
  maxx = mean(x, na.rm = T) + n * sd(x, na.rm = T)
  medx = median(x[which(x > minx & x < maxx)])
  print(paste0("Trimmed median: ", round(medx,2)))
  y <- x
  toolow <- which(x < minx)
  toohigh <- which(x > maxx)
  print(paste0(length(toolow)," people with value < ",round(minx,2)))
  print(paste0(length(toohigh)," people with value > ",round(maxx,2)))
}
diagExtremes(out$height, 2)
diagExtremes(out$height, 3)
diagExtremes(out$weight, 2)
diagExtremes(out$weight, 3)

medheight <- median(out$height[which(out$height > 1.5 & out$height < 2)])
out3 <- out %>%
  mutate(weight = fixExtremes(weight), 
         height = fixExtremes(height)) %>%
  mutate(weight = round(weight,1), 
         height = round(height,2)) %>%
  mutate(BMI = round(weight / height ^2,2)) %>%
  mutate(BMI = round(BMI, 2)) %>%
  mutate(RNA_sample_month = substr(attendanceDate___RNA,6,7)) %>%
  select(sample_id:age_RNA, inFeatureCounts, sequencingBatch, RNA_sample_month, intervalPhase:weight, BMI, attendanceDate___RNA, appointmentTime___RNA, processDate___RNA, processTime___RNA, BA_D_10_9_L___RNA:PLT_O_10_9_L___RNA, RBC_10_12_L___RNA:WBC_N_10_9_L___RNA)

# Filter out samples which aren't in the sequencing data
out4 <- out3 %>%
  filter(inFeatureCounts == 1) %>%
  select(-Sex, -inFeatureCounts)

# Add in Affy ID
omictable <- fread("processed/INTERVAL_omics_table_02APR2020.csv", data.table = F) %>%
  select(sample_id = RNA_ID, affymetrix_ID) %>% 
  filter(!is.na(sample_id))

out5 <- right_join(omictable, out4)
  
write.csv(out5, file = "processed/INTERVAL_RNA_batch1-12_master_covariates_release_2020_06_11.csv", row.names = F)