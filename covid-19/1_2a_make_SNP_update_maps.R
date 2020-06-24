# Create file to update ChrX to b38 and rename SNPs
library(dplyr)
library(data.table)
library(stringr)

setwd("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/covid19/genotypes")

# make update file for b38 positions
b38x <- fread("~/rds/rds-jmmh2-projects/interval_imputation/chrX_imputation/28_04/liftover/INTERVAL_chrX_hg38_out.bed")
b38x$snpname <- gsub(pattern = "_b37","", b38x$V4)
outmap <- b38x %>%
  select(snpname, V2)
write.table(outmap, col.names = F, row.names = F, quote = F, sep = "\t", file = "INTERVAL_chrX_update_b38pos.txt")

# Make update file to add RSIDs
alleles <- str_split_fixed(b38x$snpname, "_", 3)
b38x$A1 <- alleles[,2]
b38x$A2 <- alleles[,3]
b38x <- b38x %>%
  mutate(b38CPTID = paste0("X:",V2,"_", A1, "_", A2))
xids <- fread(data.table = F, "~/rds/rds-jmmh2-projects/covid/ace2/interval_genetic_data/interval_imputed_data/INTERVAL_imp_rsIDs.txt")
names(xids)[1] <- "snpname"

b38x2 <- left_join(b38x, xids, by = c("snpname"))
b38x2 <- b38x2 %>% 
  mutate(newSNPname = ifelse(dbSNP144_rsID != ".", dbSNP144_rsID, b38CPTID))

dupes <- b38x2$newSNPname[duplicated(b38x2$newSNPname)] %>% unique
b38x2 <- b38x2 %>%
  mutate(newnewSNPname = ifelse(newSNPname %in% dupes , paste0(newSNPname, "_", A1, "_", A2), newSNPname))
stilldupes <- b38x2$newnewSNPname[duplicated(b38x2$newnewSNPname)] %>% unique
# remove one duplicated SNP with different b37 names
b38x2 <- b38x2 %>% filter(newnewSNPname != stilldupes)

outname2 <- b38x2 %>%
  select(snpname, newnewSNPname)
write.table(outname2, col.names = F, row.names = F, quote = F, sep = "\t", file = "INTERVAL_chrX_update_rsid.txt")

# filter file to retain only SNPs with b38 positions and new RSids.
include <- outname2$newnewSNPname
write.table(include, col.names = F, row.names = F, quote = F, sep = "\t", file = "INTERVAL_chrX_b38snpfilter.txt")


