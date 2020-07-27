library(data.table)
library(dplyr)
library(chron) # works with time
library(stringr)
#setwd("/rds/user/jm2294/rds-jmmh2-projects/interval_rna_seq/analysis/00_testing/test_parameters")

#jobdf <- fread(data.table = F, "job_ids.txt")

jobdf <- data.frame("taskID" = c("16526165","17014037","17092263","17113553")) 
jobids <- jobdf$taskID

logs <- data.frame()
for(i in 1:nrow(jobdf)) {
  system(paste0("sacct -o JobID%24,AllocCPUS,State,Elapsed,MaxRSS,ReqMem --units=G -j ",jobids[i]," > temp.txt"))
  a <- fread("temp.txt", data.table = F, fill = T)
  a <- a[-1,] %>%
    filter(grepl(".batch", JobID)) %>%
    mutate(AllocCPUS = as.numeric(AllocCPUS),
           Elapsed = chron(times = Elapsed),
           ReqMem = as.numeric(gsub("Gn","",ReqMem)),
           MaxRSS = as.numeric(gsub("G","",MaxRSS)),
           taskID = str_sub(JobID, 1, 8),
           arrayID = str_sub(JobID, str_locate(JobID, "_")[,1] + 1 , str_locate(JobID, "\\.")[,1] - 1))
  a <- right_join(jobdf, a)
  logs <- rbind(logs, a)
  rm(a)
}

logs$arrayID <- as.integer(logs$arrayID)

logs <- logs %>% filter(State == "COMPLETED")
logs2 <- logs %>% 
  group_by(arrayID) %>% 
  arrange(arrayID, rev(ReqMem)) %>%
  filter(row_number() == 1) %>%
  data.frame()

usage <- logs2 %>%
  group_by(ReqMem) %>% 
  summarise(chunks = n(),
            sumTimeDays = chron(times = sum(Elapsed)),
            CPUs = unique(AllocCPUS),
            memUsed = sum(MaxRSS))

