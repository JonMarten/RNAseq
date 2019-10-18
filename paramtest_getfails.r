## Get list of failed jobs to rerun
library(data.table)
library(dplyr)
library(chron) # works with time
library(stringr)
setwd("/rds/user/jm2294/rds-jmmh2-projects/interval_rna_seq/analysis/00_testing/test_parameters")
#job <- "15309885"
#job <- "15760828"
job <- "15831871"
system(paste0("sacct -o JobID%24,AllocCPUS,State,Elapsed,MaxRSS --units=G -j ",job," > temp.txt"))
a <- fread("temp.txt", data.table = F, fill = T)
a <- a[-1,] %>%
  filter(grepl(".batch", JobID)) %>%
  mutate(AllocCPUS = as.numeric(AllocCPUS),
         Elapsed = chron(times = Elapsed),
         MaxRSS = as.numeric(str_sub(MaxRSS, 1, 5)),
         taskID = "15171908",
         arrayID = str_sub(JobID, 10, str_locate(JobID, "\\.")[,1]-1))
#a <- right_join(jobdf, a)

fails <- a %>% filter(State != "COMPLETED") %>% pull(arrayID)
failstring <- paste(fails, collapse = ",")


paste0("sbatch -a ",failstring, " submit_limix_paramtests.sh 500000 100 0.005")  