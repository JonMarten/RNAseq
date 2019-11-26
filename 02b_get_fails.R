## Get list of failed jobs to rerun
library(data.table)
library(dplyr)
library(chron) # works with time
library(stringr)
library(ggplot2)

setwd("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/slurm")
#job <- "16526165"
#job <- "17014037"
#job <- "17092259"
job <- "17092263"

# Get failed jobs from Slurm
system(paste0("sacct -o JobID%24,AllocCPUS,State,Elapsed,MaxRSS,Start,End --units=G -j ",job," > sacct_temp.txt"))
a <- fread("sacct_temp.txt", data.table = F, fill = T)
a$End <- ifelse(a$End == "", a$Start, a$End)
a <- a[-1,] %>%
  filter(grepl(".batch", JobID)) %>%
  mutate(AllocCPUS = as.numeric(AllocCPUS),
         Elapsed = chron(times = Elapsed),
         MaxRSS = as.numeric(str_sub(MaxRSS, 1, 5)),
         taskID = "15171908",
         arrayID = str_sub(JobID, 10, str_locate(JobID, "\\.")[,1]-1),
         Start = as.POSIXct(gsub("T", " ", Start)),
         End = as.POSIXct(gsub("T", " ", End))) 

# grep for errors
logpath <- paste0("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/logs")
grepCmd <- paste0("grep -i -e error -e \"can't open\" -e fault ", logpath,"/*",job,"* > grep_temp.txt")
system(grepCmd)
errorLog <- fread("grep_temp.txt", data.table = F, header = F, sep = "\t")

errorLog2 <- str_split_fixed(errorLog$V1, ".log:", 2) %>% 
  data.frame %>%
  rename(arrayID = X1, errorMsg = X2)
errorLog2$arrayID <- gsub(paste0(logpath,"/eqtl_phase1_cis_18373genes_",job,"_"), "", errorLog2$arrayID)
  
jobStatus <- full_join(a, errorLog2, by = "arrayID")

jobStatus <- jobStatus %>% 
  mutate(failed = ifelse(State != "COMPLETED" | !is.na(errorMsg), 1, 0) %>% as.factor)


# Classify Errors ---------------------------------------------------------
jobStatus$errorType <- as.character(jobStatus$errorMsg)
jobStatus$errorType[grepl("Segmentation", jobStatus$errorMsg)] <- "SegFault"
jobStatus$errorType[grepl("oom-kill", jobStatus$errorMsg)] <- "OOM-kill"
jobStatus$errorType[grepl("TIME LIMIT", jobStatus$errorMsg)] <- "Timeout"
jobStatus$errorType[grepl("unable to open file", jobStatus$errorMsg)|grepl("Could not open", jobStatus$errorMsg)] <- "File Error"
jobStatus$errorType[grepl("resource busy", jobStatus$errorMsg)] <- "Resource Busy"
jobStatus$errorType[grepl("FileNotFoundError", jobStatus$errorMsg)|grepl("does not exist", jobStatus$errorMsg)] <- "File Not Found"
jobStatus$errorType[grepl("Could not open bgen file", jobStatus$errorMsg)] <- "bgen file error"
jobStatus$errorType[grepl("metadata", jobStatus$errorMsg)] <- "bgen metadata error"
jobStatus$errorType[grepl("negative dimensions", jobStatus$errorMsg)] <- "Negative dimensions"
jobStatus$errorType[grepl("zlib", jobStatus$errorMsg)] <- "Zlib failure"
jobStatus$errorType[grepl("OSError", jobStatus$errorMsg)] <- "OSError"
jobStatus$errorType[grepl("mkdir", jobStatus$errorMsg)] <- "mkdir"
jobStatus$errorType[grepl("run_QTL_analysis.py", jobStatus$errorMsg)] <- "run_QTL_analysis.py not found"
jobStatus$errorType[grepl("slurmstepd", jobStatus$errorMsg)] <- "slurmstepd"
jobStatus$errorType[grepl("BrokenPipeError", jobStatus$errorMsg)] <- "BrokenPipeError"
jobStatus$errorType[grepl("ModuleNotFoundError", jobStatus$errorMsg)] <- "ModuleNotFoundError"
jobStatus$errorType[grepl("HDF5 error", jobStatus$errorMsg)] <- "HDF5 error"
jobStatus$errorType[grepl("is not a directory", jobStatus$errorMsg)] <- "is not a directory"

# Plot errors -------------------------------------------------------------
ggplot(filter(jobStatus), aes(colour = errorType, x = End, y = errorType, pch = as.factor(State))) +
  #geom_boxplot() +
  geom_jitter(alpha = 0.5, width = 0.2) 

ggplot(filter(jobStatus), aes(colour = errorType, x = End, y = failed)) +
  geom_jitter(alpha = 0.5, pch = 20, width = 0.2) 

ggplot(filter(jobStatus), aes(colour = errorType, x = Start, y = End)) +
  geom_point(alpha = 0.5, pch = 20) 


jobFailed <- jobStatus %>%
  filter(failed == 1)

write.csv(jobStatus, file = paste0("slurm_job_",job, "_status.csv"))


# Generate resubmission command
findIntRuns <- function(run){
  rundiff <- c(1, diff(run))
  difflist <- split(run, cumsum(rundiff!=1))
  unlist(lapply(difflist, function(x){
    if(length(x) %in% 1:2) as.character(x) else paste0(x[1], "-", x[length(x)])
  }), use.names=FALSE)
}

failstring <- paste(unique(jobFailed$arrayID), collapse = ",")
s <- "1,2,3,4,8,9,14,15,16,19"
s2 <- as.numeric(unlist(strsplit(failstring, ",")))

failruns <- paste0(findIntRuns(s2), collapse=",")

paste0("sbatch -a ",failruns, " 02_submit_limix.sh")  

# Look at successful jobs
jobSucceeded <- jobStatus %>%
  filter(failed == 0)






