setwd( "/home/jm2294/projects/RNAseq/test_run_chunks/output_newpipeline_nofilter/comparison_test")
library(data.table)
library(dplyr)

ass <- fread("assqtl_results_22_10736171_17730855.txt", data.table=F)

david <- fread("/home/jm2294/projects/RNAseq/v01_188id_preview/davids_calls/20PCs_Removed_eQTLsFDR0.05-SNPLevel.txt", data.table=F)

remap <- fread("/home/jm2294/GENETIC_DATA/INTERVAL/RNAseq/b37_b38_liftover/INTERVAL_chr22_b37_to_b38_map.txt", data.table=F)
names(remap) <- c("altid_b37","rsid_b37","chr_b37","pos_b37","ref_b37", "alt_b37","altid_b38","rsid_b38","chr_b38","pos_b38","ref_b38", "alt_b38")

d2 <- david %>%
  filter(ProbeName %in% ass$feature_id) %>%
  rename(rsid_b37 = SNPName)

d3 <- left_join(d2, remap)  %>%
  mutate(mergeid = paste0(ProbeName, "_", altid_b38))

library(stringr)
tem <- str_split_fixed(d3$"Beta (SE)", "\\(",2)
tem[,1] <- substr(tem[,1],1,nchar(tem[,1])-1)
tem[,2] <- substr(tem[,2],1,nchar(tem[,2])-1)
tem2 <- data.frame("beta" = as.numeric(tem[,1]), "SE" = as.numeric(tem[,2]))

d4 <- cbind(d3, tem2)


ass2 <- ass %>% 
  rename(altid_b38 = snp_id) %>%
  mutate(mergeid = paste0(feature_id, "_", altid_b38))
d5 <- left_join(d4, ass2, by = "mergeid")

d5 <- d5 %>% filter(!is.na(beta.y))

d5$p_david <- pchisq()
d6 <- filter(d5, p_value < 0.05)



