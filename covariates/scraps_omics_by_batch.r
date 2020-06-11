init()
setwd("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/04_phase2_full_analysis")
dataRelease <- "14MAY2020"
otab <- fread(data.table = F, paste0("covariates/processed/INTERVAL_omics_table_",dataRelease,".csv"))

cov <- fread("covariates/processed/INTERVAL_RNA_batch1-12_master_covariates_release_2020_05_18.csv", data.table = F)

cov2 <- cov %>%
  select(sample_id, sequencingBatch) %>%
  rename(RNA_ID = sample_id)

m <- right_join(otab, cov2, by = "RNA_ID")


sumTab <- m %>% group_by(sequencingBatch) %>%
  summarise(RNA = n(),
            Affymetrix = length(which(!is.na(affymetrix_ID))),
            WES = length(which(!is.na(WES_ID))),
            WGS = length(which(!is.na(WGS_ID))),
            Nightingale = length(which(!is.na(brainshake_ID))),
            Metabolon = length(which(!is.na(metabolon_ID))),
            Somalogic = length(which(!is.na(soma_ID))),
            Olink = length(which(!is.na(olink_ID))),
  )
           
library(tidyr)  
library(viridis)
sumTall <- sumTab %>%
  #select(-Affy, -Nightingale) %>% 
  gather(-sequencingBatch, key = "Ome", value ="n") %>%
  mutate(Ome = factor(Ome, levels = c("Nightingale", "Metabolon", "WGS","WES","Affymetrix", "Olink","Somalogic","RNA")))

ggplot(sumTall, aes(x = as.factor(sequencingBatch), y = n, fill = Ome)) +
  geom_bar(stat = "identity", color="black", position=position_dodge()) +
  scale_fill_viridis(option = "plasma", discrete = TRUE) +
  scale_y_continuous(expand = c(0, 0))

omicplot <- ggplot(sumTall, aes(x = Ome, y = n, fill = Ome)) +
  geom_bar(stat = "identity", color="black", position=position_dodge()) +
  scale_fill_viridis(option = "plasma", discrete = TRUE) +
  facet_wrap(~as.factor(sequencingBatch)) +
  coord_flip() +
  scale_y_continuous(expand = c(0, 0))

ggsave(omicplot, file = "covariates/omic_plot.png")