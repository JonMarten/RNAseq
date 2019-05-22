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

# plots to check ages
library(tidyr)
library(ggplot2)
library(cowplot)
library(yarrr)

pal <- piratepal(palette = "xmen") %>% unname

ageplot <- ggplot(filter(dat3, !is.na(phase) & !is.na(batch)), 
                  aes(y = ageDif, x = as.factor(phase), fill = as.factor(batch))) +   
  geom_jitter(pch = 21, width = 0.3, height = 0.01, size = 2) +
  scale_fill_manual(values = pal, name = "Batch") +
  labs(y = "Years after initial appointment", x = "INTERVAL 'phase' used for RNA")
ggsave(plot = ageplot, filename = "scripts/RNAseq/ageplot.png")

datTall <- dat2 %>% 
  select(sample_id, agePulse, age1, age24m, age48m, agep3) %>%
  gather(age1:agep3, key = "phase", value = "age", -sample_id)

ggplot(datTall, aes(x = agePulse, y = age, colour = phase)) + 
  geom_point() +
  geom_abline(slope = 1, intercept = 0) 

## Write out covariate file 
covOut <- dat3 %>% 
  select(sample_id, sexPulse, age_RNA, batch, phase) %>% 
  filter(!is.na(phase) & !is.na(batch))

fwrite(sep = "\t", covOut, file = "covariates/INTERVAL_RNA_batch1_2_covariates_sex_age.txt", quote = F, na = NA)
