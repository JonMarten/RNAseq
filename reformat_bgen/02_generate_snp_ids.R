# Make new alternate ID for bgen 1.2 files
setwd("/home/jm2294/rds/rds-jmmh2-projects/interval_rna_seq/GENETIC_DATA")
library(dplyr)
library(data.table)
library(stringr)

args <- commandArgs(trailingOnly=TRUE)
chr <- as.numeric(args[1])
#for(chr in 1:22){

# Read in SNP stats output from qctool and add column to check if variant is an indel
snpstats <- fread(paste0("snp_stats/impute_",chr,"_interval_snp_stats_unfiltered.txt"), skip = 8, data.table=F)
snpstats <- snpstats %>%
  mutate(indel = ifelse(nchar(alleleA)>1 | nchar(alleleB)>1, 1, 0)) %>%
  #filter(rsid != ".") %>%
  mutate(alleleA_sort = ifelse(alleleA <= alleleB, alleleA, alleleB),
         alleleB_sort = ifelse(alleleA <= alleleB, alleleB, alleleA)) %>%
  mutate(match_id = paste0(str_pad(chr,2,pad = "0"),"_", position, "_", alleleA, "_", alleleB, ":", rsid))

# read in b38 map
b38 <- fread(paste0("b37_b38_liftover/from_savita_INTERVAL_imputed_liftover_hg38_cleaned/INTERVAL_imputed_liftover_hg38_cleaned_chr",chr,".vcf"), data.table = F, skip = 5)
names(b38) <- c("CHROM.b38","POS.b38","match_id","REF.b38","ALT.b38")

#b38 <- b38 %>%
#  distinct %>%
#  mutate(REF.b38_sort = ifelse(REF.b38 < ALT.b38, REF.b38, ALT.b38),
#         ALT.b38_sort = ifelse(REF.b38 < ALT.b38, ALT.b38, REF.b38))  %>%
#  mutate(match_id = paste0(rsid,"_",REF.b38_sort,"_",ALT.b38_sort))

snpsmerge <- inner_join(snpstats, b38, by = "match_id")

# Make list of SNPs that don't have b38 positions
no_b38 <- snpstats$match_id[which(!snpstats$match_id %in% snpsmerge$match_id)]
no_snpstats <- b38$match_id[which(!b38$match_id %in% snpsmerge$match_id)]

cat("\n",length(no_b38), "variants do not match to b38.")
cat("\n",length(no_snpstats), "variants are missing from qctool snpstats.")

#drops <- snpstats %>%
#  filter(match_id %in% no_b38) %>%
#  mutate(chrpos = paste0(chromosome,":",position)) %>%
#  pull(chrpos)

# Create CPTIDs for SNPs (sort alleles in alphabetical order) USING ALLELES FROM b37 HERE.
snps <- snpsmerge %>%
  filter(indel == 0) %>%
  mutate(cptid.b38 = ifelse(rsid == ".",
                            paste0(CHROM.b38, "_", POS.b38, "_", alleleA_sort, "_", alleleB_sort),
                            paste0(CHROM.b38, "_", POS.b38, "_", alleleA_sort, "_", alleleB_sort,":",rsid)
                            )
         )

# Create CPTIDs for indels (sort alleles by smallest first)
indels <- snpsmerge %>%
  filter(indel == 1) %>%
  filter(!is.na(CHROM.b38)) %>%
  mutate(alleleA_sort = ifelse(nchar(alleleA) <= nchar(alleleB), alleleA, alleleB),
         alleleB_sort = ifelse(nchar(alleleA) > nchar(alleleB), alleleA, alleleB)) %>%
  mutate(cptid.b38 = ifelse(rsid == ".",
                            paste0(CHROM.b38, "_", POS.b38, "_", alleleA_sort, "_", alleleB_sort),
                            paste0(CHROM.b38, "_", POS.b38, "_", alleleA_sort, "_", alleleB_sort,":",rsid)
                            )
         )

# Merge SNPs and indels back together
b <- rbind(snps, indels) %>%
  arrange(chromosome, position)

# remove duplicate rows. This appears to be more or less equivalent to identifying duplicates in plink, and since the final run includes only those SNPs included in this list, there should be no risk of missing any due to duplicates not mapping to b38 and so not being in the b38 vcf.
dupes <- b %>% 
  group_by(cptid.b38) %>%
  filter(n() > 1) %>%
  data.frame()

if(nrow(dupes > 0)){
  b <- b %>%
    filter(!match_id %in% dupes$match_id)
}

# For qctool -map-id-data	
# Update the chromosome, position, IDs and/or alleles of a set of SNPs with new values. 
newMap <- b %>%
  mutate(rsid=rsid, alternate_id.b38=cptid.b38,rsid.38=rsid, chromosome.38=CHROM.b38, position.b38=POS.b38,alleleA.b38=alleleA, alleleB.b38=alleleB) %>%
  select(alternate_ids, rsid, chromosome, position, alleleA, alleleB, alternate_id.b38, rsid.38, chromosome.38, position.b38, alleleA.b38, alleleB.b38) %>%
  mutate(rsid.38 = ifelse(rsid.38 == ".", alternate_id.b38, rsid.38))

if(chr < 10){
  newMap$chromosome <- paste0("0", as.character(newMap$chromosome))
  newMap$chromosome.38 <- paste0("0", as.character(newMap$chromosome.38))
}

write.table(newMap, 
            quote = F, 
            sep = " ", 
            col.names = T, 
            row.names = F, 
            file = paste0("b37_b38_liftover/INTERVAL_chr",chr,"_b37_to_b38_map.txt"))

# As not all SNPs are currently mappable in the b38 data due to ambiguous rsids (eg snps coded as "."), SNPs have to be filtered to retain only those which have unambiguous b38 positions.
#--incl-snps	Exclude all SNPs not in the given file(s) from the analysis. The format of this file is the same as that output by the -write-snp-excl-list option. It must have six columns interpreted as SNPID, rsid, chromosome, position, first and second alleles. 
snplist <- newMap %>%
  select(alternate_id.b38)
write.table(snplist, row.names = F, col.names = F, quote = F, file = paste0("b37_b38_liftover/c",chr,"_b38_filter_snps.txt"))


# Alternative for faster filtering: use bgenix. Problem - filters on RSID, not alt_id. Stupid me should have made those the same but I didn't because I thought I was being clever, and I have no compute budget left to do it. So instead of an inclusion filter, I'm going to use an exclusion filter. This might mean some SNPs are lost that would otherwise be retained, but this is better than the reverse. This will remove all remaining SNPs with "." as rsid but that's okay since the SNPs to be retained have had this changed to cptid instead.
no_b38_df <- snpstats[which(!snpstats$match_id %in% snpsmerge$match_id),]
filter_rsids <- no_b38_df$rsid.b38 %>% unique
lostsnps <- which(filter_rsids %in% newMap$rsid.38)
cat("\n",length(lostsnps), "variants will be lost if filtering by rsid in bgenix.")
write.table(filter_rsids, row.names = F, col.names = F, quote = F, file = paste0("b37_b38_liftover/c",chr,"_b38_rsids_to_remove_bgenix.txt"))

                  