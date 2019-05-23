# scrap code to check ages are the same for 24, 48 and p3 samples where an id is duplicated.
# Uses rnaids object from make_id_mapper_file.R and dat3 from reformat_covs_for_limix.R

rnaids <- rnaids %>% 
  mutate(has_id_24 = ifelse(is.na(RNAseq_RAW_24m),0,1),
         has_id_48 = ifelse(is.na(RNAseq_RAW_48m),0,1),
         has_id_p3 = ifelse(is.na(RNAseq_RAW_p3),0,1)) %>%
  mutate(has_id_sum = has_id_24 + has_id_48 + has_id_p3)

dupes <- rnaids %>% 
  filter(has_id_sum > 1) %>%
  rename(sample_id = RNA_any)

dat4 <- left_join(dupes, dat3)
dat4 %>% 
  select(sample_id, has_id_24, has_id_48, has_id_p3, age24m, age48m, agep3)
