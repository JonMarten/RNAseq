library(data.table)
library(dplyr)
library(chron) # works with time
library(stringr)
setwd("/rds/user/jm2294/rds-jmmh2-projects/interval_rna_seq/analysis/00_testing/test_parameters")

jobdf <- fread(data.table = F, "job_ids.txt")
jobdf$taskID <- as.character(jobdf$taskID)
jobdf <- jobdf %>% filter(!is.na(taskID)& Features == "Coding")
jobids <- jobdf$taskID

logs <- data.frame()
for(i in 1:nrow(jobdf)) {
  system(paste0("sacct -o JobID%24,AllocCPUS,State,Elapsed,MaxRSS --units=G -j ",jobids[i]," > temp.txt"))
  a <- fread("temp.txt", data.table = F, fill = T)
  a <- a[-1,] %>%
    filter(grepl(".batch", JobID)) %>%
    mutate(AllocCPUS = as.numeric(AllocCPUS),
           Elapsed = chron(times = Elapsed),
           MaxRSS = as.numeric(str_sub(MaxRSS, 1, 5)),
           taskID = str_sub(JobID, 1, 8),
           arrayID = str_sub(JobID, 10, 13))
  a <- right_join(jobdf, a)
  logs <- rbind(logs, a)
  rm(a)
}

# add in re-run chunks 1107 and 1099 for 15171908
system(paste0("sacct -o JobID%24,AllocCPUS,State,Elapsed,MaxRSS --units=G -j 15233831 > temp.txt"))
a <- fread("temp.txt", data.table = F, fill = T)
a <- a[-1,] %>%
  filter(grepl(".batch", JobID)) %>%
  mutate(AllocCPUS = as.numeric(AllocCPUS),
         Elapsed = chron(times = Elapsed),
         MaxRSS = as.numeric(str_sub(MaxRSS, 1, 5)),
         taskID = "15171908",
         arrayID = str_sub(JobID, 10, 13))
a <- right_join(jobdf, a)
logs <- rbind(logs, a) %>%
  filter(State != "CANCELLED")

# summarise
sumtab <- logs %>%
  group_by(taskID) %>%
  summarise(Cohort = unique(Cohort),
            Features = unique(Features),
            Bgen_size = unique(Bgen_size),
            Window = unique(Window), 
            Covariates = unique(Covariates),
            MAF = unique(MAF),
            Permutations = unique(Permutations),
            num_run = n(),
            num_failed = length(which(State=="CANCELLED")),
            mean_time = mean(Elapsed, na.rm = T),
            med_time = median(Elapsed, na.rm = T),
            sd_time = sd(Elapsed, na.rm = T),
            min_time = min(Elapsed, na.rm = T),
            max_time = max(Elapsed, na.rm = T),
            mean_mem = mean(MaxRSS, na.rm = T),
            sd_mem = sd(MaxRSS, na.rm = T),
            med_mem = median(MaxRSS, na.rm = T),
            min_mem = min(MaxRSS, na.rm = T),
            max_mem = max(MaxRSS, na.rm = T)) %>%
  data.frame()

sumtab <- sumtab %>%
  mutate(mean_time = as.character(mean_time),
         med_time = as.character(med_time),
         sd_time = as.character(sd_time),
         min_time = as.character(min_time),
         max_time = as.character(max_time))
fwrite(sumtab, file = "parameter_comparison_joblogs.csv")

l2 <- logs  %>% mutate(MAF = as.factor(MAF), Permutations = as.factor(Permutations), Window = as.factor(Window))
lm(data = l2 , Elapsed ~ Window + MAF + Permutations) %>% summary
lm(data = l2 , MaxRSS ~ Window + MAF + Permutations) %>% summary