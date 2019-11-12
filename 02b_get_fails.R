## Get list of failed jobs to rerun
library(data.table)
library(dplyr)
library(chron) # works with time
library(stringr)
setwd("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/slurm")
job <- "16526165"

# Get failed jobs from Slurm
system(paste0("sacct -o JobID%24,AllocCPUS,State,Elapsed,MaxRSS --units=G -j ",job," > sacct_temp.txt"))
a <- fread("sacct_temp.txt", data.table = F, fill = T)
a <- a[-1,] %>%
  filter(grepl(".batch", JobID)) %>%
  mutate(AllocCPUS = as.numeric(AllocCPUS),
         Elapsed = chron(times = Elapsed),
         MaxRSS = as.numeric(str_sub(MaxRSS, 1, 5)),
         taskID = "15171908",
         arrayID = str_sub(JobID, 10, str_locate(JobID, "\\.")[,1]-1))

# grep for errors
logpath <- paste0("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/logs")
grepCmd <- paste0("grep -i error ", logpath,"/*",job,"* > grep_temp.txt")
system(grepCmd)
errorLog <- fread("grep_temp.txt", data.table = F, header = F, sep = "\t")

errorLog2 <- str_split_fixed(errorLog$V1, ".log:", 2) %>% 
  data.frame %>%
  rename(arrayID = X1, errorMsg = X2)
errorLog2$arrayID <- gsub(paste0(logpath,"/eqtl_phase1_cis_18373genes_",job,"_"), "", errorLog2$arrayID)
  
jobStatus <- full_join(a, errorLog2, by = "arrayID")

jobFailed <- jobStatus %>%
  filter(State != "COMPLETED" | !is.na(errorMsg))

write.csv(jobFailed, file = paste0("slurm_job_",job, "_failed.csv"))

failstring <- paste(jobFailed$arrayID, collapse = ",")
paste0("sbatch -a ",failstring, " 02_submit_limix.sh")  