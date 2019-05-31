# Make new alternate ID for bgen 1.2 files
setwd("/home/jm2294/GENETIC_DATA/INTERVAL/RNAseq")
library(dplyr)
library(data.table)

chr <- 22

# Read in SNP stats output from qctool and add column to check if variant is an indel
snpstats <- fread(paste0("/home/jm2294/GENETIC_DATA/INTERVAL/RNAseq/snp_stats/impute_",chr,"_interval_snp_stats_unfiltered.txt"), skip = 8, data.table=F)
snpstats <- snpstats %>%
  mutate(indel = ifelse(nchar(alleleA)>1 | nchar(alleleB)>1, 1, 0)) %>%
  filter(rsid != ".") %>%
  mutate(alleleA_sort = ifelse(alleleA <= alleleB, alleleA, alleleB),
         alleleB_sort = ifelse(alleleA <= alleleB, alleleB, alleleA)) %>%
  mutate(match_id = paste0(rsid,"_",alleleA_sort,"_",alleleB_sort))

# read in b38 map
b38 <- fread(paste0("/home/jm2294/GENETIC_DATA/INTERVAL/RNAseq/b37_b38_liftover/INTERVAL_24_8_18_imputed_chr",chr,"_hg38.vcf"), data.table = F, skip = 5)
names(b38) <- c("CHROM.b38","POS.b38","rsid","REF.b38","ALT.b38")
b38 <- b38 %>%
  filter(rsid != ".") %>%
  mutate(REF.b38_sort = ifelse(REF.b38 < ALT.b38, REF.b38, ALT.b38),
         ALT.b38_sort = ifelse(REF.b38 < ALT.b38, ALT.b38, REF.b38))  %>%
  mutate(match_id = paste0(rsid,"_",REF.b38_sort,"_",ALT.b38_sort))

snpsmerge <- inner_join(snpstats, b38, by = "match_id")
no_b38 <- snpstats$match_id[which(!snpstats$match_id %in% snpsmerge$match_id)]

# Create CPTIDs for SNPs (sort alleles in alphabetical order)
snps <- snpsmerge %>%
  filter(indel == 0) %>%
  filter(!is.na(CHROM.b38)) %>%
  mutate(REF.b38_sort = ifelse(REF.b38 < ALT.b38, REF.b38, ALT.b38),
         ALT.b38_sort = ifelse(REF.b38 < ALT.b38, ALT.b38, REF.b38)) %>%
  mutate(cptid.b38 = paste0(CHROM.b38, ":", POS.b38, "_", REF.b38_sort, "_", ALT.b38_sort))

# Create CPTIDs for indels (sort alleles by smallest first)
indels <- snpsmerge %>%
  filter(indel == 1) %>%
  filter(!is.na(CHROM.b38)) %>%
  mutate(REF.b38_sort = ifelse(nchar(REF.b38) <= nchar(ALT.b38), REF.b38, ALT.b38),
         ALT.b38_sort = ifelse(nchar(REF.b38) > nchar(ALT.b38), REF.b38, ALT.b38)) %>%
  mutate(cptid.b38 = paste0(CHROM.b38, ":", POS.b38, "_", REF.b38_sort, "_", ALT.b38_sort))

# Check for duplicated ids on SNPs
snpdupe <- snps %>% 
  group_by(cptid.b38) %>%
  filter(n() > 1) %>%
  data.frame()%>%
  mutate(alleleA_sort = ifelse(alleleA <= alleleB, alleleA, alleleB),
         alleleB_sort = ifelse(alleleA <= alleleB, alleleB, alleleA)) %>%
  mutate(cptid = paste0(CHROM.b38, ":", POS.b38, "_", alleleA_sort, "_", alleleB_sort))

# Check for duplicated ids on indels
dupes <- indels %>% 
  group_by(cptid) %>%
  filter(n() > 1) %>%
  data.frame()

# Append rsid to indels with duplicated cptids
if(nrow(dupes > 0)){
  dupevec <- which(indels$cptid %in% dupes$cptid)
  indels$cptid[dupevec] <- paste0(indels$cptid[dupevec], "___",indels$rsid[dupevec])
}

# Merge SNPs and indels back together
b <- rbind(snps, indels) %>%
  arrange(chromosome, position)




