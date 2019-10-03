## Get list of failed jobs to rerun

setwd("/rds/user/jm2294/rds-jmmh2-projects/interval_rna_seq/analysis/00_testing/test_parameters")
job <- "15309885"
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

fails <- a %>% filter(State != "COMPLETED") %>% pull(arrayID)
failstring <- paste(fails, collapse = ",")