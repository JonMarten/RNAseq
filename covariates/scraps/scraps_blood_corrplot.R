init()
allcov <- fread("processed/INTERVAL_RNA_batch1-12_master_covariates_release_2020_05_18.csv", data.table = F)
bloodnames <- fread("blood_trait_names.csv", data.table = F)
bloodcol <- which(colnames(allcov) %in% bloodnames$Variable_name[bloodnames$Use == 1])
blood <- allcov[,bloodcol] %>%
  na.exclude

names(blood) <- bloodnames$humanName[match(names(blood), bloodnames$Variable_name)]

library(corrplot)

corrmat <- cor(blood)
corrplot(corrmat)