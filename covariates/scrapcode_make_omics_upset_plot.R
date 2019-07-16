# Uses out object from make_omics_table.R
#install.packages("UpSetR")
library(UpSetR)

binary <- function(x){ifelse(is.na(x),0,1)}

upset.df <- out %>%
  select("NIHR Consent" = NIHRConsent, 
         "WES" = WES_ID, 
         "WGS" = WGS_ID, 
         "Brainshake" = brainshake_ID, 
         "Metabolon" = metabolon_ID, 
         "Somalogic" = soma_ID, 
         "Olink" = olink_ID, 
         "RNA seq" = picklistRNA)
upset.df.bin <-apply(upset.df, MAR = 2, FUN = binary) %>% data.frame

upset(upset.df.bin, sets = c("WES","WGS","Metabolon", "Somalogic","Olink","RNA.seq"))


