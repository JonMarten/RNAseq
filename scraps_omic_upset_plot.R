init()
setwd("/rds/user/jm2294/rds-jmmh2-projects/interval_rna_seq/analysis/04_phase2_full_analysis/covariates/processed")
a <- fread("INTERVAL_omics_table_02APR2020.csv", data.table = F)

a <- a %>% filter(!is.na(RNA_ID))
 
library(UpSetR)

affymetrix_ID <- a %>% 
  filter(!is.na(affymetrix_ID)) %>%
  pull(data_management_project_id)

WES_ID <- a %>% 
  filter(!is.na(WES_ID)) %>%
  pull(data_management_project_id)

WGS_ID <- a %>% 
  filter(!is.na(WGS_ID)) %>%
  pull(data_management_project_id)

brainshake_ID <- a %>% 
  filter(!is.na(brainshake_ID)) %>%
  pull(data_management_project_id)

metabolon_ID <- a %>% 
  filter(!is.na(metabolon_ID)) %>%
  pull(data_management_project_id)

soma_ID <- a %>% 
  filter(!is.na(soma_ID)) %>%
  pull(data_management_project_id)

olink_ID <- a %>% 
  filter(!is.na(olink_ID)) %>%
  pull(data_management_project_id)

RNA_ID <- a %>% 
  filter(!is.na(RNA_ID)) %>%
  pull(data_management_project_id)

idlist <- list("WES" = WES_ID, 
               "WGS" = WGS_ID, 
               "Brainshake" = brainshake_ID, 
               "Metabolon" = metabolon_ID, 
               "SomaLogic" = soma_ID, 
               "Olink" = olink_ID, 
               "RNA-seq" = RNA_ID)

upset(fromList(idlist), order.by = "freq", nsets = 7, text.scale = 2.5, nintersects = 7, point.size = 5, mb.ratio = c(0.5, 0.5), shade.color = "gray50")

library(VennDiagram)
venn.diagram(x = list("SomaLogic" = soma_ID, "Olink" = olink_ID, "RNA-seq" = RNA_ID), filename = "omic_venn.png")